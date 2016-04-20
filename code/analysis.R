##---
##title: 'Adam Zeilinger ZIB model'
##author: Daniel Turek
##publish: true
##---






#### Introduction

##We'll look at two different ways of parametrizing your model.  The first uses latent states, and is the traditional form for occupancy models (the only to fit these models in most software).  The second parametrization will be NIMBLE-only, and we'll see is much faster.


##I've did some processing of your data, to remove the huge numble of NA's from the raw observation matrix.  Instead of large 2-dimensional arrays, with only roughly about 600 actual observations, the data is now in 'flat' structures: 1-dimensional vectors containing *only* the 600 actual obsevations, with a new site-membership indicator variable `siteID` for correctly specifying the random effects. To see the data processing that takes place, take a look at `code/create_data.R`.

##```{r, message = FALSE}
#### load libraries
require(nimble)
require(rjags)
require(coda)

#### load the processed data
load('../data/zib_data.RData')

#### some other function definitions for custom samplers, ploting, etc.
source('definitions.R')

#### this has to do with how the page is being made
fromR <- TRUE
##```

##```{r, echo=FALSE}
fromR <- FALSE
##```






#### Latent state model

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



##We'll use both NIMBLE and JAGS to fit the latent state model.

##We'll use 22,000 iterations, after 2000 burnin will leave 20,000 samples.

##```{r, eval = FALSE}
niter <- 22000

set.seed(0)

comp_LS <- compareMCMCs(modelInfo_LS,
                        niter = niter,
                        summary = TRUE,
                        MCMCs = c('jags','nimble'))

save(comp_LS, file = '../cached/comp_LS.RData')
##```

##```{r load-comp_LS, echo = FALSE}
load(file = '../cached/comp_LS.RData')
##```

## We'll take a look at how quickly that ran, and the resulting effective sample size (ESS).

##```{r}
#### runtime (in minutes) of each algorithm
comp_LS$LS$timing / 60
##```

##JAGS takes about 5 minutes, NIMBLE takes about 30 seconds, to produce 20,000 samples.

##ESS takes into account the autocorrelation in the posterior chains, and is a measure of the actual number of *independent* posterior samples.  At this point, it's worth you taking a few minutes to read about this metric of MCMC performance.  There's a good description of it [available here](https://nature.berkeley.edu/~pdevalpine/MCMC_comparisons/nimble_MCMC_comparisons.html).

## Next is the *minimum efficiency* for all sampled parameters, which is the number of effective samples produces per second of runtime.

##```{r}
comp_LS$LS$efficiency$min
##```

##JAGS has minimum efficiency about 0.01, meaning it produces about 1 sample every 100 seconds. NIMBLE has minimum efficiency of 0.14, meaning it produces about 1 sample every 7 seconds.

##This is too slow, since we'll want upwards of 100,000 posterior samples, at least.  We'll stop looking at the latent state model here, and move on to the next, faster, model representation.






####Custom distribution representation

##Now we'll try writing your model in a slightly different way.  I've written a custom distribution for use in the model specification, called `dOccupancy()`.  You can see the definition of this in `code/definitions.R`.  As parameters it takes the occupancy probability and the observation probability, and it returns the probability density.

##You can see the new model specification below, without the latent `z` variables, and using this new custom distribution for the observations `y`.  This model can only be used with NIMBLE.

##```{r define-dOcc-model}
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

##- nimble: NIMBLE's default MCMC.

##- block: Putting block samplers on the coefficients in each linear predictor term, which are very likely to be correlated in the posterior.

##- log_shft: We will see that `sigma_alpha` is the slowest mixing parameter, so we'll try to speed that up (since the mixing of this parameter will therefore limit our inferential power for the rest of the model). This MCMC samples `sigma_alpha` on a log scale, and also helps it mix by shifting the `alpha` random effects proportionally for every change in `sigma_alpha`.

##- log_shft_blk: Combination of the above two ideas. Using block sampling for coefficients, and also the log-shift sampling for `sigma_alpha`.

##- log_shft2_blk: Same as above, but also does shift sampling for `mu_alpha`, that is it shifts the `alpha` random effect values by the same amount for every change in `mu_alpha`, which will help the mixing of `mu_alpha`.

##If you take a minute to look at the `MCMCdefs` argument list below, you'll get some insight about how MCMC algorithms can be customized in NIMBLE.

##```{r run-dOcc-model, eval = FALSE}
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

##Here's a really nice feature of NIMBLE's MCMC comparisons feature: the MCMC comparisons pages.  This next command will generate nice visual summaries of the performance of all the MCMC algorithms we just ran.

##```{r make-dOcc-comparison-page, eval = FALSE}
make_MCMC_comparison_pages(comp_dOcc, dir = '../html')
##```

##It's worth spending at least a few minutes looking over [these comparison pages](http://danielturek.github.io/zib_glmm/html/dOccupancy.html).  At this point we can start making some conclusions.

##Looking at the Posterior Summaries section at the bottom of the comparisons page, we see that all MCMC algorithms appear to have converged nicely to basically the same posterior distribution.  The only one with any noticeable differences is the posterior for `sigma_alpha`, but the differences are minor.  Also, for the purposes of your biological anlysis, the standard deviation of the site random effects is probably not the most critical model parameter.  We'll check convergence otherwise, but this gives us good confidence in the results.

##Looking at the minimum Efficiency for each MCMC, the log_shft_blk is the best overall, where the mixing is still limited by `sigma_alpha`, producing only 6 effectively independent samples per second of runtime.

##We'll accept this MCMC algorithm, log_shft_blk, as the "best" hereafter, although most all of these algorithms would be just fine to use.

##Before we go ahead an run it longer, let's examine a few more things.

##See how long it took to actually run these, for 500,000 samples.

##```{r print-dOcc-timing}
comp_dOcc$dOccupancy$timing / 60
##```

##log_shft_blk algorithm took about 3 minutes, not too bad at all.

##What was the actual Effective Sample Size (ESS) it produced for all parameters?  Let's make sure it's somewhat reasonable.

##```{r print-dOcc-ess}
comp_dOcc$dOccupancy$summary['log_shft_blk', 'ess',]
##```

##The smallest ESS is for `sigma_alpha` (as we knew), with about 1,200 independent posterior samples.  The others are all more.  This is ok, but we'll run the algorithm for longer the final time.

##Let's also just take a look at the mixing for `sigma_alpha`, to make sure it looks ok.

##```{r make-sigma-alpha-traceplot, fig.width=8, fig.height=3.5}
if(fromR) dev.new(width=8, height=3.5)
tsplot(comp_dOcc$dOccupancy$samples['log_shft_blk', 'sigma_alpha', 300000:500000], main='sigma_alpha', xlab='', ylab='')
##```

##The mixing is good. Not the best we can imagine, and it's definitely upper-truncated at 500, but this is pretty good overall.







#### Final MCMC runs

##Finally, we'll go ahead and run two chains of this MCMC, double-check convergence, and look at some posterior statistics.

##The code below uses NIMBLE's core MCMC functionality.  We'll create and compile this one MCMC algorithm, and run two chains each of 1,000,000 iterations.  This will give us reasonably good inferentialy power for all parameters.

##```{r, eval = FALSE}
#### create the R model
Rmodel <- nimbleModel(modelInfo_dOccupancy$code,
                      modelInfo_dOccupancy$constants,
                      modelInfo_dOccupancy$data,
                      modelInfo_dOccupancy$inits)

#### compile the model to C++
Cmodel <- compileNimble(Rmodel)

#### create a default MCMC specification
spec <- configureMCMC(Rmodel)

#### configure this specification the same as log_shft_blk from above
spec$removeSamplers('beta[1:10]')
spec$addSampler('beta[1:3]', 'RW_block')
spec$addSampler('beta[4:10]', 'RW_block')
spec$removeSamplers('sigma_alpha')
spec$addSampler('sigma_alpha', 'RW_log_shift', list(shiftNodes='alpha'))

#### build the R MCMC algorithm, and compile it to C++
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
##```

##Now we'll run two chains of this MCMC algorithm, each with one million MCMC iterations.

##We'll also manually remove the first 100,000 samples as burnin.

##```{r, eval = FALSE}
niter <- 1000000

#### run the first chain and collect samples
set.seed(1)
Cmcmc$run(niter)
samples1 <- as.matrix(Cmcmc$mvSamples)[100001:niter, ]

#### run the second chain and collect samples
set.seed(2)
Cmcmc$run(niter)
samples2 <- as.matrix(Cmcmc$mvSamples)[100001:niter, ]

## create mcmc and mcmc.list objects for plotting and convergence diagnostics
mcmc1 <- coda::as.mcmc(samples1)
mcmc2 <- coda::as.mcmc(samples2)
mcmcs <- coda::mcmc.list(mcmc1, mcmc2)

save(samples1, samples2, mcmcs, file = '../cached/final_MCMC.RData')
##```

##```{r, echo = FALSE}
load(file = '../cached/final_MCMC.RData')
##```


#### Assessing convergence

##Here we'll assess convergence from these two MCMC chains.

##We'll use the mainstream Gelman & Rubin convergence diagnostic from the `coda` package.  We want to see values near to 1 here.

##```{r gelman-diag}
coda::gelman.diag(mcmcs, autoburnin = FALSE)
##```

##Convergence looks excellent, on all counts.

##Let's also just take a quick look at the Effective Sample Size (ESS) that we've actually collected for each parameter.  We'll just look at the first chain.

##```{r final-mcmc-ess}
apply(samples1, 2, coda::effectiveSize)
##```

##These vary a lot, but we have over 2,000 independent samples for `sigma_alpha`, and more for all the other parameters.  Perhaps we could desire more, but this will be enough for our inferences.  Of course if you want more samples, you can use this code to generate more.


#### Posterior Inferences

##Finally!  We have lots of samples that we're happy with.  Now we can now make posterior inferences.  For extra assurance, we'll make the inferences from both chain1 and chain2, and expect to see (essentially) the same results.

## Let me know if you want any help on how to interpet all these numbers.  As far as interpreting the actual `beta[i]` parameters, you'll have to just look back at the model specifications above, but that's easy enough.


##### Posterior Mean

##```{r}
cbind(apply(samples1, 2, mean),
      apply(samples2, 2, mean))
##```


##### Posterior Median

##```{r}
cbind(apply(samples1, 2, median),
      apply(samples2, 2, median))
##```

##We see the posterior distributions are only mildly skewed, the one exception being `sigma_alpha`.


##### Standard Error

##```{r}
nsamples <- dim(samples1)[1]

cbind(apply(samples1, 2, sd),
      apply(samples2, 2, sd)) / sqrt(nsamples)
##```

##Very pleasingly small standard error for each parameter.  That's a direct result of having been able to collect so many samples from a highly-efficient MCMC algorithm  =)


##### Posterior Density Plots

##We'll only plot posterior densities from the first chain.

##```{r posterior-plots}
plot(mcmcs[[1]], ask = FALSE)
##```






#### Wrap Up

##I suspect there might still be some changes or modifications to your model or dataset.  If that's the case, it's ok.  Just keep in touch.

##Also, I'll be happy to contribute a short bit on fitting the model for any publications that come out of this work.

##Cheers!  - Daniel

