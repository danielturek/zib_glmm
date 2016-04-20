##---
##title: 'Adam Zeilinger ZIB model'
##author: Daniel Turek
##publish: true
##---


##```{r, echo = FALSE, eval = FALSE}
####```{r, echo = FALSE, fig.with=10, fig.height=3.5}
####if(fromR) dev.new()
##```




##### Introduction

##We'll look at two different ways of parametrizing your model.  The first uses latent states, and is the traditional form for occupancy models (the only to fit these models in most software).  The second parametrization will be NIMBLE-only, and we'll see is much faster.


##I've did some processing of your data, to remove the huge numble of NA's from the raw observation matrix.  Instead of large 2-dimensional arrays, with only roughly about 600 actual observations, the data is now in 'flat' structures: 1-dimensional vectors containing *only* the 600 actual obsevations, with a new site-membership indicator variable `siteID` for correctly specifying the random effects. To see the data processing that takes place, take a look at `code/create_data.R`.

##```{r, message = FALSE}
#### load the nimble library
library(nimble)

#### load the processed data
load('../data/zib_data.RData')

#### some other function definitions for custom samplers, ploting, etc.
source('definitions.R')

#### this has to do with how the page is being made
fromR <- TRUE
##```







##### Latent state model

##The latent state version of the model includes the latent variables `z`, indicating true presence / absense.  Including these introduces about 600 additional unknown variables into the model which must be sampled, which slows down the MCMC a lot.

##In the second formulation of the model, we'll use a custom distribution in NIMBLE to remove these latent states.  But to use JAGS (for the first comparison), we can't use custom distributions, and therefore use the latent state model representation.

##```{r}
#### latent state model code
code <- nimbleCode({
    mu_alpha ~ dnorm(0, 0.001)
    sigma_alpha ~ dunif(0, 500)
    tau_alpha <- 1 / (sigma_alpha * sigma_alpha)
    for(j in 1:nsite) { 
        alpha[j] ~ dnorm(mu_alpha, tau_alpha)  ## site random effect
    }
    for(i in 1:10) {
        beta[i] ~ dnorm(0, 0.001)
    }
    ## notice the single index i, using the new 'flat' data structures
    ## also the siteID group membership, used to assign the correct random effect
    for(i in 1:N) {
        logit(p_occ[i]) <- alpha[siteID[i]] + beta[4]*aet[i] + beta[5]*tmn[i] + beta[6]*tmx[i] + beta[7]*year[i] + beta[8]*month[i] + beta[9]*year2[i] + beta[10]*month2[i]
        logit(p_obs[i]) <- beta[1] + beta[2]*list_length[i] + beta[3]*year_list_length[i]
        z[i] ~ dbern(p_occ[i])
        mu_y[i] <- z[i] * p_obs[i]
        y[i] ~ dbern(mu_y[i])
    }
})

#### notice everything is the 'flat' structure
constants <- list(N=N, nsite=nsite, aet=aetflat, tmn=tmnflat, tmx=tmxflat, year=yearflat, year2=year2flat, month=monthflat, month2=month2flat, list_length=list_lengthflat, year_list_length=year_list_lengthflat, siteID=siteIDflat)

#### notice everything is the 'flat' structure
data <- list(y=yflat)

inits <- list(mu_alpha=0, sigma_alpha=1, alpha=rep(0,nsite), beta=rep(0,10), z=rep(1,N))

#### LS indicates model information for the 'latent state' version
modelInfo_LS <- list(code=code, constants=constants, data=data, inits=inits, name='LS')
##```



##We'll use both NIMBLE and JAGS to fit the latent state model

##We'll use 22,000 iterations, after 2000 burnin will leave 20,000 samples

##```{r, eval = FALSE}
niter <- 22000

set.seed(0)

comp_LS <- compareMCMCs(modelInfo_LS,
                        niter = niter,
                        summary = TRUE,
                        MCMCs = c('jags','nimble'))

save(comp_LS, file = '../cached/comp_LS.RData')
##```

##```{r, echo = FALSE}
load(file = '../cached/comp_LS.RData')
##```

## Take a look at how quickly that ran, and the resulting effective sample size (ESS).  ESS takes into account the autocorrelation in the posterior chains, and is a measure of the actual number of *independent* posterior samples.

##```{r}
#### runtime (in minutes) of each algorithm
comp_LS$LS$timing / 60
##```

##JAGS takes about 5 minutes, NIMBLE takes about 30 seconds, to produce 20,000 samples.

## Next is the *minimum efficiency* for all sampled parameters, which is the number of effective samples produces per second of runtime.

##```{r}
comp_LS$LS$efficiency$min
##```

##JAGS has minimum efficiency about 0.01, meaning it produces about 1 sample every 100 seconds. NIMBLE has minimum efficiency of 0.14, meaning it produces about 1 sample every 7 seconds.

##This is too slow, since we'll want upwards of 100,000 posterior samples, at least.  We'll stop looking at the latent state model here, and move on to the next, faster, model representation.






#####Custom distribution representation

##Now we'll try writing your model in a slightly different way.  I've written a custom distribution for use in the model specification, called `dOccupancy()`.  You can see the definition of this in `code/definitions.R`.  As parameters it takes the occupancy probability and the observation probability, and it returns the probability density.

##You can see the new model specification below, without the latent `z` variables, and using this new custom distribution for the observations `y`.  This model can only be used with NIMBLE.

##```{r}
#### latent state model code
code <- nimbleCode({
    mu_alpha ~ dnorm(0, 0.001)
    sigma_alpha ~ dunif(0, 500)
    for(j in 1:nsite) { 
        alpha[j] ~ dnorm(mu_alpha, sd = sigma_alpha)  ## site random effect
    }
    for(i in 1:10) {
        beta[i] ~ dnorm(0, 0.001)
    }
    for(i in 1:N) {
        logit(p_occ[i]) <- alpha[siteID[i]] + beta[4]*aet[i] + beta[5]*tmn[i] + beta[6]*tmx[i] + beta[7]*year[i] + beta[8]*month[i] + beta[9]*year2[i] + beta[10]*month2[i]
        logit(p_obs[i]) <- beta[1] + beta[2]*list_length[i] + beta[3]*year_list_length[i]
        y[i] ~ dOccupancy(p_occ[i], p_obs[i])
    }
})

#### still using the 'flat' data structures
constants <- list(N=N, nsite=nsite, aet=aetflat, tmn=tmnflat, tmx=tmxflat, year=yearflat, year2=year2flat, month=monthflat, month2=month2flat, list_length=list_lengthflat, year_list_length=year_list_lengthflat, siteID=siteIDflat)

#### still using the 'flat' data structures
data <- list(y=yflat)

#### Note: no longer have initial values for 'z' variables,
#### since they've been removed from the model.
inits <- list(mu_alpha=0, sigma_alpha=1, alpha=rep(0,nsite), beta=rep(0,10))

modelInfo_dOccupancy <- list(code=code, constants=constants, data=data, inits=inits, name='dOccupancy')
##```



##Now we'll do some runs with this model.  Since it will run much faster, we'll do 502,000 iterations.  I'm only choosing that number such that after the burnin of 2,000 iterations, we'll retain 500,000 posterior samples.

##This time, we'll consider a large number of different MCMC algorithms.  This serves two purposes:

##- Finding the most efficient algorithm

##- Comparing posteriors to help us assess convergence

##I did a fair bit of experimenting (not shown here) to arrive at this list of well-performing candidate MCMCs, but it's still intersting to look at the performance of all of them.  Just to briefly describe each MCMC algorithm:

##- a

##- b

##```{r, eval = FALSE}
niter <- 502000

set.seed(0)

comp_dOcc <- compareMCMCs(modelInfo_dOccupancy,
                          MCMCs = c('nimble', 'block', 'log_shft', 'log_shft_blk', 'log_shft2_blk'),
                          niter = niter,
                          summary = FALSE,
                          MCMCdefs = list(
                              block = quote({
                                  spec <- configureMCMC(Rmodel)
                                  spec$removeSamplers('beta[1:10]')
                                  spec$addSampler('beta[1:3]', 'RW_block')
                                  spec$addSampler('beta[4:10]', 'RW_block')
                                  spec }),
                              log_shft = quote({
                                  spec <- configureMCMC(Rmodel)
                                  spec$removeSamplers('sigma_alpha')
                                  spec$addSampler('sigma_alpha', 'RW_log_shift', list(shiftNodes='alpha'))
                                  spec }),
                              log_shft_blk = quote({
                                  spec <- configureMCMC(Rmodel)
                                  spec$removeSamplers('beta[1:10]')
                                  spec$addSampler('beta[1:3]', 'RW_block')
                                  spec$addSampler('beta[4:10]', 'RW_block')
                                  spec$removeSamplers('sigma_alpha')
                                  spec$addSampler('sigma_alpha', 'RW_log_shift', list(shiftNodes='alpha'))
                                  spec }),
                              log_shft2_blk = quote({
                                  spec <- configureMCMC(Rmodel)
                                  spec$removeSamplers('beta[1:10]')
                                  spec$addSampler('beta[1:3]', 'RW_block')
                                  spec$addSampler('beta[4:10]', 'RW_block')
                                  spec$removeSamplers('sigma_alpha')
                                  spec$addSampler('sigma_alpha', 'RW_log_shift', list(shiftNodes='alpha'))
                                  spec$removeSamplers('mu_alpha')
                                  spec$addSampler('mu_alpha', 'RW_shift', list(shiftNodes='alpha', skipDependencies=FALSE))
                                  spec })
                          ))
                                  
save(comp_dOcc, file = '../cached/comp_dOcc.RData')
##```

##```{r, echo = FALSE}
load(file = '../cached/comp_dOcc.RData')
##```



##```{r, eval = FALSE}
comp_dOcc$dOccupancy$timing
comp_dOcc$dOccupancy$efficiency

make_MCMC_comparison_pages(comp_dOcc, dir = '../html')
system('open ../html/dOccupancy.html')


samples <- comp_dOcc$dOccupancy$samples
dim(samples)

iPlot <- 400000:498000
tsplot(comp_dOcc$dOccupancy$samples['nimble',       'sigma_alpha', iPlot])
tsplot(comp_dOcc$dOccupancy$samples['block',        'sigma_alpha', ])
tsplot(comp_dOcc$dOccupancy$samples['log',          'sigma_alpha', ])
tsplot(comp_dOcc$dOccupancy$samples['log_shft',     'sigma_alpha', ])
tsplot(comp_dOcc$dOccupancy$samples['log_shft_blk', 'sigma_alpha', ])

samples <- comp_dOcc$dOccupancy$samples
iPlot <- 400000:500000


MCMCs <- c('nimble', 'block', 'log_shft', 'log_shft_blk', 'log_shft2_blk')
nMCMCs <- length(MCMCs)
dev.new(height=8, width=8)
par(mfrow = c(nMCMCs, 1))
for(i in 1:nMCMCs) {
    tsplot(samples[MCMCs[i], 'sigma_alpha', iPlot], main=MCMCs[i])
}
##```





##make_MCMC_comparison_pages(comp_LS, dir = '../html')








