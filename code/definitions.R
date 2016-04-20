
##################################################
### Latent state version of zip model ###########
##################################################

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
    ##for(i in 1:nlist) {        # i: events (year-months)
    ##    for(j in 1:nsite) {    # j: sites
    for(i in 1:N) {
        logit(p_occ[i]) <- alpha[siteID[i]] + beta[4]*aet[i] + beta[5]*tmn[i] + beta[6]*tmx[i] + beta[7]*year[i] + beta[8]*month[i] + beta[9]*year2[i] + beta[10]*month2[i]
        logit(p_obs[i]) <- beta[1] + beta[2]*list_length[i] + beta[3]*year_list_length[i]
        z[i] ~ dbern(p_occ[i])
        mu_y[i] <- z[i] * p_obs[i]
        y[i] ~ dbern(mu_y[i])
    }
})

constants <- list(N=N, nsite=nsite, aet=aetflat, tmn=tmnflat, tmx=tmxflat, year=yearflat, year2=year2flat, month=monthflat, month2=month2flat, list_length=list_lengthflat, year_list_length=year_list_lengthflat, siteID=siteIDflat)

data <- list(y=yflat)

inits <- list(mu_alpha=0, sigma_alpha=1, alpha=rep(0,nsite), beta=rep(0,10), z=rep(1,N))

modelInfo_LS <- list(code=code, constants=constants, data=data, inits=inits, name='LS')

rm(list = c('code', 'constants', 'data', 'inits'))




##################################################
### dOccupancy version of zip model ##############
##################################################

dOccupancy <- nimbleFunction(
    run = function(x = double(), p_occ = double(), p_obs = double(), log.p = double()) {
        if(x == 0) L <- 1 - p_occ * p_obs else L <- p_occ * p_obs
        returnType(double())
        if(log.p) return(log(L)) else return(L)
    }
)

rOccupancy <- nimbleFunction(
    run = function(n = integer(), p_occ = double(), p_obs = double()) {
        print('not implemented')
        returnType(double())
        return(0)
    }
)

registerDistributions(list(
    dOccupancy = list(
        BUGSdist = 'dOccupancy(p_occ, p_obs)',
        discrete = TRUE
    )
))

code <- nimbleCode({
    mu_alpha ~ dnorm(0, 0.001)
    sigma_alpha ~ dunif(0, 500)
    for(j in 1:nsite) { 
        alpha[j] ~ dnorm(mu_alpha, sd = sigma_alpha)  ## site random effect
    }
    for(i in 1:10) {
        beta[i] ~ dnorm(0, 0.001)
    }
    ##for(i in 1:nlist) {        # i: events (year-months)
    ##    for(j in 1:nsite) {    # j: sites
    for(i in 1:N) {
        logit(p_occ[i]) <- alpha[siteID[i]] + beta[4]*aet[i] + beta[5]*tmn[i] + beta[6]*tmx[i] + beta[7]*year[i] + beta[8]*month[i] + beta[9]*year2[i] + beta[10]*month2[i]
        logit(p_obs[i]) <- beta[1] + beta[2]*list_length[i] + beta[3]*year_list_length[i]
        y[i] ~ dOccupancy(p_occ[i], p_obs[i])
    }
})

constants <- list(N=N, nsite=nsite, aet=aetflat, tmn=tmnflat, tmx=tmxflat, year=yearflat, year2=year2flat, month=monthflat, month2=month2flat, list_length=list_lengthflat, year_list_length=year_list_lengthflat, siteID=siteIDflat)

data <- list(y=yflat)

inits <- list(mu_alpha=0, sigma_alpha=1, alpha=rep(0,nsite), beta=rep(0,10))

modelInfo_dOccupancy <- list(code=code, constants=constants, data=data, inits=inits, name='dOccupancy')

rm(list = c('code', 'constants', 'data', 'inits'))



##################################################
### time-series plot #############################
##################################################


tsplot <- function(y, ...) plot(x = seq_along(y), y = y, type = 'l', ...)




##################################################
### RW_shift and RW_shift_log samplers ##########
##################################################

## Sampler for a target node that offsets a set of other nodes by an opposite amount as the proposed random walk step for the target node
## Control element skipDepenendencies can be TRUE if you are SURE that the resulting likelihood will be completely unmodified
## Otherwise skipDependencies should be FALSE.

sampler_RW_shift <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ###  control list extraction  ###
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        scale         <- control$scale
        shiftNodes    <- model$expandNodeNames(control$shiftNodes)
        skipDependencies <- control[['skipDependencies']]
        if(is.null(skipDependencies)) skipDependencies <- FALSE
        ###  node list generation  ###
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        if(length(targetAsScalar) > 1)     stop('more than one target; cannot use RW sampler, try RW_block sampler')

        if(!skipDependencies) {
            calcNodes  <- model$getDependencies(c(target, shiftNodes))
            determCalcNodes <- character()
        }
        else {
            calcNodes <- c(target, shiftNodes)
            determCalcNodes <- model$getDependencies(c(target, shiftNodes), determOnly = TRUE, self = FALSE)
        }
        ###  numeric value generation  ###
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory          <- c(0, 0)
        acceptanceRateHistory <- c(0, 0)
        ## variables previously inside of nested functions:
        optimalAR <- 0.44
        gamma1    <- 0
    },
    
    run = function() {
        propShift <-  rnorm(1, mean = 0, sd = scale)
        model[[target]] <<- model[[target]] + propShift
        values(model, shiftNodes) <<- values(model, shiftNodes) - propShift
        logMHR <- calculateDiff(model, calcNodes)
        jump <- decide(logMHR)
        if(jump) {
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
            if(skipDependencies) { ## This ensures that deterministic dependencies are updated.
                calculate(model, determCalcNodes)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = determCalcNodes, logProb = TRUE)
            }
        }
        else
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        if(adaptive)     adaptiveProcedure(jump)
    },
    
    methods = list(
        
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                setSize(scaleHistory,          timesAdapted)
                setSize(acceptanceRateHistory, timesAdapted)
                scaleHistory[timesAdapted] <<- scale
                acceptanceRateHistory[timesAdapted] <<- acceptanceRate
                gamma1 <<- 1/((timesAdapted + 3)^0.8)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        
        reset = function() {
            scale <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            scaleHistory          <<- scaleHistory          * 0
            acceptanceRateHistory <<- acceptanceRateHistory * 0
            gamma1 <<- 0
        }
    ), where = getLoadingNamespace()
)


## Sampler for a standard deviation (typically) on a log scale that shifts the magnitude of a set of other nodes by the same value, also on a log scale.
## This helps mixing for a set of random effects and their shared standard deviation

sampler_RW_log_shift <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ###  control list extraction  ###
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        scale         <- control$scale
        shiftNodes    <- model$expandNodeNames(control$shiftNodes)
        ###  node list generation  ###
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        if(length(targetAsScalar) > 1)     stop('more than one target; cannot use RW sampler, try RW_block sampler')

        numTargets <- 1 + length(shiftNodes)
        calcNodes  <- model$getDependencies(c(target, shiftNodes))

        ###  numeric value generation  ###
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory          <- c(0, 0)
        acceptanceRateHistory <- c(0, 0)
        ## variables previously inside of nested functions:
        optimalAR <- 0.44
        gamma1    <- 0
    },
    
    run = function() {
        propLogShift <-  rnorm(1, mean = 0, sd = scale)
        propMult <- exp(propLogShift)
        model[[target]] <<- model[[target]] * propMult
        newVals <- values(model, shiftNodes) * propMult
        values(model, shiftNodes) <<- newVals 
        logMHR <- calculateDiff(model, calcNodes) + numTargets * propLogShift
        jump <- decide(logMHR)
        if(jump)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        else
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        if(adaptive)     adaptiveProcedure(jump)
    },
    
    methods = list(
        
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                setSize(scaleHistory,          timesAdapted)
                setSize(acceptanceRateHistory, timesAdapted)
                scaleHistory[timesAdapted] <<- scale
                acceptanceRateHistory[timesAdapted] <<- acceptanceRate
                gamma1 <<- 1/((timesAdapted + 3)^0.8)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        
        reset = function() {
            scale <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            scaleHistory          <<- scaleHistory          * 0
            acceptanceRateHistory <<- acceptanceRateHistory * 0
            gamma1 <<- 0
        }
    ), where = getLoadingNamespace()
)




