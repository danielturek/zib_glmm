##---
##title: 'Adam Zeilinger ZIB model (v.3) Month Fixed'
##author: Daniel Turek
##publish: true
##---





#### Updates (v.3)

##- Changed the month covariate into a factor, that is, a fixed value 1, 2, 3, ..., 12 corresponding to each month, instead of a standardized continuous variable as before.

##- Part I: Used this month factor as a fixed effect in the model, and removed the old month variables: month, quadratic month, and interactions between month and year

##- Part II: Added a second model, *without* the monthly fixed effets, since none of them had any significant effect.

##- Also included 90% credible intervals for the second model (Part II), since several of the `beta[i]` become significant at the 90% level.







#### Updates (v.2)

##- Using new data set with N = 958 observations (lists >= 4 species)

##- Removed quadratic year^2 effect, and added year:month and year:month^2 interactions to occpancy process

##- Transformed new siteID (raster cell numbers) into sequential indexing of random effects

##- Changed the prior for sigma_alpha form uniform(0, 500) to uniform(0, 1000), to capture more of the long right-tail of the posterior.








#### Preliminaries

##This time, we'll jump straight to the most efficient MCMCs.

##There are two parts.  Part I considers the use of fixed effects for each month, and Part II drops these fixed effects. 


##```{r, message = FALSE}
require(nimble)
require(rjags)
require(coda)

#### this is the new data set, with month as a factor: 1, 2, 3, ..., 12
inputData <- readRDS('data/data_nimble_zib.rds')

source('zib_glmm/code/definitions.R')

fromR <- TRUE
##```

##```{r, echo=FALSE}
fromR <- FALSE
##```






####Part I: Month as fixed effects


##Same as before, except month will be a fixed effect.

##```{r define-dOcc-model-season}
code <- nimbleCode({
    mu_alpha ~ dnorm(0, 0.001)
    sigma_alpha ~ dunif(0, 1000)
    for(j in 1:nsite) { 
        alpha[j] ~ dnorm(mu_alpha, sd = sigma_alpha)  ## site random effect
    }
    for(i in 1:7) {
        beta[i] ~ dnorm(0, 0.001)
    }
    for(i in 1:4) {
        betaseason[i] ~ dnorm(0, 0.001)    ## new fixed effects for each season
    }
    for(i in 1:N) {
        logit(p_occ[i]) <- alpha[siteID[i]] + beta[4]*aet[i] + beta[5]*tmn[i] + beta[6]*tmx[i] + beta[7]*year[i] + betaseason[season[i]]
        logit(p_obs[i]) <- beta[1] + beta[2]*list_length[i] + beta[3]*year_list_length[i]
        y[i] ~ dOccupancy(p_occ[i], p_obs[i])
    }
})

constants <- with(inputData,
                  list(N=N, nsite=nsite, aet=aet, tmn=tmn, tmx=tmx, year=year, season=season, list_length=list_length, year_list_length=year_list_length, siteID=siteID))

data <- with(inputData, list(y=y))

inits <- list(mu_alpha=0, sigma_alpha=1, alpha=rep(0,inputData$nsite), beta=rep(0,7), betaseason=rep(0,12))

modelInfo_season <- list(code=code, constants=constants, data=data, inits=inits, name='season_model')
##```





#### Part I: MCMC runs with month fixed effects


##```{r, eval = FALSE}
Rmodel <- nimbleModel(modelInfo_season$code,
                      modelInfo_season$constants,
                      modelInfo_season$data,
                      modelInfo_season$inits)

Cmodel <- compileNimble(Rmodel)

spec <- configureMCMC(Rmodel)

#### configure this specification the same as log_shft_blk as before
spec$removeSamplers('beta[1:7]')
spec$addSampler('beta[1:3]', 'RW_block')
spec$addSampler('beta[4:7]', 'RW_block')
spec$removeSamplers('sigma_alpha')
spec$addSampler('sigma_alpha', 'RW_log_shift', list(shiftNodes='alpha'))

Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
##```

##Now we'll run two chains of this MCMC algorithm, each with 500,000 MCMC iterations.

##We'll manually remove the first 100,000 samples as burnin.

##```{r, eval = FALSE}
niter <- 500000

set.seed(1)
Cmcmc$run(niter)
samples1 <- as.matrix(Cmcmc$mvSamples)[100001:niter, ]

set.seed(2)
Cmcmc$run(niter)
samples2 <- as.matrix(Cmcmc$mvSamples)[100001:niter, ]

mcmc1 <- coda::as.mcmc(samples1)
mcmc2 <- coda::as.mcmc(samples2)
mcmcs <- coda::mcmc.list(mcmc1, mcmc2)

save(samples1, samples2, mcmcs, file = 'zib_glmm/final_MCMC_season.RData')
##```

##```{r, echo = FALSE}
load(file = '../cached/final_MCMC_fixedmonth.RData')
##```


#### Part I: Assessing convergence

##Now we assess convergence using these two MCMC chains.

##We'll use the mainstream Gelman & Rubin convergence diagnostic from the `coda` package.  We want to see all values near to one.

##```{r gelman-diag}
coda::gelman.diag(mcmcs, autoburnin = FALSE)
##```

##Convergence looks excellent.

##```{r final-mcmc-ess}
apply(samples1, 2, coda::effectiveSize)
##```





#### Part I: Posterior Inferences

##### Part I: Posterior Mean

##```{r}
cbind(apply(samples1, 2, mean),
      apply(samples2, 2, mean))
##```

##### Part I: Posterior Median

##```{r}
cbind(apply(samples1, 2, median),
      apply(samples2, 2, median))
##```

##### Part I: 95% Credible Intervals

##```{r}
cbind(apply(samples1, 2, function(x) quantile(x, 0.025)),
      apply(samples1, 2, function(x) quantile(x, 0.975)))
##```



##Note that *none* of the monthly fixed effects `betamonth` are significant.  In fact, they're fairly normally distributed around zero.

##So we'll simplify the model in Part II, and drop these monthly fixed effects.







#### Part II: Without monthly fixed effects

##Since none of the monthly fixed effects had any significant effect, we'll remove them from the model and re-run.  This saves us from fitting these 12 parameters, and should make all the other inferences more accurate.

##Same as before, except without the month fixed effects.

##```{r define-dOcc-model-no-month}
code <- nimbleCode({
    mu_alpha ~ dnorm(0, 0.001)
    sigma_alpha ~ dunif(0, 1000)
    for(j in 1:nsite) { 
        alpha[j] ~ dnorm(mu_alpha, sd = sigma_alpha)  ## site random effect
    }
    for(i in 1:7) {
        beta[i] ~ dnorm(0, 0.001)
    }
    for(i in 1:N) {
        logit(p_occ[i]) <- alpha[siteID[i]] + beta[4]*aet[i] + beta[5]*tmn[i] + beta[6]*tmx[i] + beta[7]*year[i]
        logit(p_obs[i]) <- beta[1] + beta[2]*list_length[i] + beta[3]*year_list_length[i]
        y[i] ~ dOccupancy(p_occ[i], p_obs[i])
    }
})

constants <- list(N=N, nsite=nsite, aet=aet, tmn=tmn, tmx=tmx, year=year, list_length=list_length, year_list_length=year_list_length, siteID=siteID)

data <- list(y=y)

inits <- list(mu_alpha=0, sigma_alpha=1, alpha=rep(0,nsite), beta=rep(0,7))

modelInfo_nomonth <- list(code=code, constants=constants, data=data, inits=inits, name='nomonth')
##```




#### Part II: MCMC runs without month fixed effects


##```{r, eval = FALSE}
Rmodel <- nimbleModel(modelInfo_nomonth$code,
                      modelInfo_nomonth$constants,
                      modelInfo_nomonth$data,
                      modelInfo_nomonth$inits)

Cmodel <- compileNimble(Rmodel)

spec <- configureMCMC(Rmodel)

#### configure this specification the same as log_shft_blk as before
spec$removeSamplers('beta[1:7]')
spec$addSampler('beta[1:3]', 'RW_block')
spec$addSampler('beta[4:7]', 'RW_block')
spec$removeSamplers('sigma_alpha')
spec$addSampler('sigma_alpha', 'RW_log_shift', list(shiftNodes='alpha'))

Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
##```

##Now we'll run two chains of this MCMC algorithm, each with 500,000 MCMC iterations.

##We'll manually remove the first 100,000 samples as burnin.

##```{r, eval = FALSE}
niter <- 500000

set.seed(1)
Cmcmc$run(niter)
samples1 <- as.matrix(Cmcmc$mvSamples)[100001:niter, ]

set.seed(2)
Cmcmc$run(niter)
samples2 <- as.matrix(Cmcmc$mvSamples)[100001:niter, ]

mcmc1 <- coda::as.mcmc(samples1)
mcmc2 <- coda::as.mcmc(samples2)
mcmcs <- coda::mcmc.list(mcmc1, mcmc2)

save(samples1, samples2, mcmcs, file = '../cached/final_MCMC_nomonth.RData')
##```

##```{r, echo = FALSE}
load(file = '../cached/final_MCMC_nomonth.RData')
##```


#### Part II: Assessing convergence

##Now we assess convergence using these two MCMC chains.

##We'll use the mainstream Gelman & Rubin convergence diagnostic from the `coda` package.  We want to see all values near to one.

##```{r gelman-diag-no-month}
coda::gelman.diag(mcmcs, autoburnin = FALSE)
##```

##Convergence looks excellent.

##```{r final-mcmc-ess-no-month}
apply(samples1, 2, coda::effectiveSize)
##```


##Notice the *improved* effective sample size for `sigma_alpha`.  The mixing improved noticeably when we removed the extra 12 parameters that weren't contributing anything to the model fit.







#### Part II: Posterior Inferences

##### Part II: Posterior Mean

##```{r}
cbind(apply(samples1, 2, mean),
      apply(samples2, 2, mean))
##```

##### Part II: Posterior Median

##```{r}
cbind(apply(samples1, 2, median),
      apply(samples2, 2, median))
##```

##### Part II: 95% Credible Intervals

##```{r}
cbind(apply(samples1, 2, function(x) quantile(x, 0.025)),
      apply(samples1, 2, function(x) quantile(x, 0.975)))
##```


##### Part II: 90% Credible Intervals

##I'm also including 90% credible intervals here, since both `beta[2]` and `beta[3]` are significant at the 90% level.

##```{r}
cbind(apply(samples1, 2, function(x) quantile(x, 0.05)),
      apply(samples1, 2, function(x) quantile(x, 0.95)))
##```





##### Part II: Posterior Density Plots

##```{r, eval = FALSE}
plot(mcmcs[[1]], ask = FALSE)
##```

##```{r posterior-plots-no-month, echo = FALSE}
for(i in 1:3)
    plot(coda::as.mcmc(samples1[, (3*i-2):(3*i)]), ask = FALSE)
##```












#### Wrap Up

##I think this last version of the model (Part II) is our best yet.  The handling of the month variable truly needed to be changed.

##Also, now we see significant effects from a number of the variables, with some nice temporal trends as you had hoped.

##Let me know when you've had a chance to look over this.

##Cheers!  - Daniel

