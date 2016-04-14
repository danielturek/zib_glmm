rrr
library(nimble)
load('../data/zib_data.RData')
source('definitions.R')
niter <- 20000


## Latent State model (suitable for JAGS)
if(FALSE) {
    set.seed(0)
    comp_LS <- compareMCMCs(modelInfo_LS, MCMCs=c('jags','nimble'), niter=niter, summary=TRUE)
    comp_LS$LS$timing
    comp_LS$LS$efficiency
    make_MCMC_comparison_pages(comp_LS, dir = '../html')
}


## dOccupancy model
set.seed(0)
comp_dOcc <- compareMCMCs(modelInfo_dOccupancy,
                          MCMCs = c('nimble', 'block', 'log', 'log_shft', 'log_shft_blk'),
                          niter = niter,
                          summary = FALSE,
                          MCMCdefs = list(
                              block = quote({
                                  spec <- configureMCMC(Rmodel)
                                  spec$removeSamplers('beta[1:10]')
                                  spec$addSampler('beta[1:2]', 'RW_block')
                                  spec$addSampler('beta[3:10]', 'RW_block')
                                  spec }),
                              log = quote({
                                  spec <- configureMCMC(Rmodel)
                                  spec$removeSamplers('sigma_alpha')
                                  spec$addSampler('sigma_alpha', 'RW', list(log=TRUE))
                                  spec }),
                              log_shft = quote({
                                  spec <- configureMCMC(Rmodel)
                                  spec$removeSamplers('sigma_alpha')
                                  spec$addSampler('sigma_alpha', 'RW_log_shift', list(shiftNodes='alpha'))
                                  spec }),
                              log_shft_blk = quote({
                                  spec <- configureMCMC(Rmodel)
                                  spec$removeSamplers('beta[1:10]')
                                  spec$addSampler('beta[1:2]', 'RW_block')
                                  spec$addSampler('beta[3:10]', 'RW_block')
                                  spec$removeSamplers('sigma_alpha')
                                  spec$addSampler('sigma_alpha', 'RW_log_shift', list(shiftNodes='alpha'))
                                  spec })
                              ))
                                  

comp_dOcc$dOccupancy$timing
comp_dOcc$dOccupancy$efficiency

tsplot(comp_dOcc$dOccupancy$samples['nimble',       'sigma_alpha', ])
tsplot(comp_dOcc$dOccupancy$samples['block',        'sigma_alpha', ])
tsplot(comp_dOcc$dOccupancy$samples['log',          'sigma_alpha', ])
tsplot(comp_dOcc$dOccupancy$samples['log_shft',     'sigma_alpha', ])
tsplot(comp_dOcc$dOccupancy$samples['log_shft_blk', 'sigma_alpha', ])


make_MCMC_comparison_pages(comp_dOcc, dir = '../html')















##combined <- combine_MCMC_comparison_results(comp[[1]], comp2[[1]])
## 
##make_MCMC_comparison_pages(combined, dir = '../html')
## 
##comp3 <- compareMCMCs(modelInfo,
##                     MCMCs = c('custom2'),
##                     MCMCdefs = list(custom2 = quote({
##                         spec <- configureMCMC(Rmodel)
##                         spec$removeSamplers('sigma_alpha')
##                         spec$addSampler('sigma_alpha', 'slice')
##                         spec$removeSamplers('z')
##                         for(node in Rmodel$expandNodeNames('z'))   spec$addSampler(node, 'slice', print=FALSE)
##                         spec
##                     })),
##                     niter = niter,
##                     summary = TRUE)

##my_packages<-c('data.table', 'snow', 'dclone', 'rjags', 'R2jags')
##lapply(my_packages, require, character.only=TRUE)
##jagsGLMMdata <- readRDS("Data_JAGS_GLMM.rds")

##codeJagsModelFile <- function(code, file) {
##    sink(file)
##    cat(paste0('123model ', paste0(deparse(code, width.cutoff=500L), collapse='\n'), '\n'))
##    sink()
##}
##dir <- tempdir()
##dir <- '.'
##jagsModelFile <- paste0(dir, '/', 'zib_JAGS.model')
##codeToJagsModelFile(code, jagsModelFile)

#### Non-parallel jags()
##ni=81000; nt=10; nc=3
##n.adapt <- 500; n.update <- 500
##nb <- n.adapt + n.update
## glmmOutput <- jags(data = jagsGLMMdata, 
##                 inits = inits, 
##                 parameters.to.save = params, 
##                 model.file = "climate_zib_glmm_model.jags", 
##                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
##                 working.directory = getwd())    
## print(glmmOutput)

##system.time(Rmodel <- nimbleModel(code, constants, data, inits))
##system.time(spec <- configureMCMC(Rmodel))
##spec$getSamplers()
##system.time(Rmcmc <- buildMCMC(spec))
##system.time(Cmodel <- compileNimble(Rmodel))
##system.time(Cmcmc <- compileNimble(Rmcmc, project = Rmodel))
##set.seed(0)
##Cmcmc$run(10000)
##samples <- as.matrix(Cmcmc$mvSamples)
##apply(samples, 2, mean)
## 
##out <- MCMCsuite(code, constants, data, inits,
##                 MCMCs = c('nimble', 'jags'),
##                 calculateEfficiency = TRUE,
##                 makePlot = FALSE)
##out$timing
##out$efficiency

