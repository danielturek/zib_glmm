rrr

createNamedObjectsFromList <- function(lst, writeToFile = NULL, envir = parent.frame()) {
    for(i in seq_along(lst)) {
        objName <- names(lst)[i]
        obj <- eval(lst[[i]])
        assign(objName, obj, envir = envir)
    }
    if(!is.null(writeToFile)) {
        write('', file = writeToFile)
        for(i in seq_along(lst)) {
            expr <- substitute(VAR <- VALUE, list(VAR = as.name(names(lst)[i]), VALUE = lst[[i]]))
            deparseExpr <- deparse(expr, control=c())
            deparseExpr <- gsub('\"', '\'', deparseExpr)
            write(deparseExpr, file = writeToFile, append = TRUE)
            write('\n\n\n', file = writeToFile, append = TRUE)
        }
    }
}

##jagsGLMMdata <- readRDS('../orig_files/Data_JAGS_GLMM_v1.rds')
jagsGLMMdata <- readRDS('../orig_files/Data_JAGS_GLMM_v2.rds')

createNamedObjectsFromList(jagsGLMMdata)

y <- detectionMatrix
year2 <- year^2     ## quadratic year effect
month2 <- month^2   ## quadratic month effect
year_list_length <- year * list_length   ## year:list_length interaction term

ind <- which(!is.na(y))
N <- length(ind)

aetflat <- aet[ind]
##cwdflat <- cwd[ind]   ## dropped 'cwd' covarite in v2 of model
tmxflat <- tmx[ind]
tmnflat <- tmn[ind]
yearflat <- year[ind]
year2flat <- year2[ind]
monthflat <- month[ind]
month2flat <- month2[ind]
list_lengthflat <- list_length[ind]
year_list_lengthflat <- year_list_length[ind]  ## year:list_length interaction term
yflat <- y[ind]

siteMatrix <- array(rep(1:nsite, each=nlist), c(nlist, nsite))
siteIDflat <- siteMatrix[ind]

save(list = c(
         'N', 'nsite',
         ## 'cwdflat',     ## dropped 'cwdflat' covariate
         'aetflat', 'tmxflat', 'tmnflat',
         'yearflat', 'year2flat', 'monthflat', 'month2flat',
         'list_lengthflat',
         'year_list_lengthflat',  ## year:list_length interaction term
         'yflat',
         'siteIDflat'
),
     file = '../data/zib_data.RData')

qqq



