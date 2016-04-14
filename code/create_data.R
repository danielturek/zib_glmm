
rm(list=ls())

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

jagsGLMMdata <- readRDS('../orig_files/Data_JAGS_GLMM.rds')
createNamedObjectsFromList(jagsGLMMdata)

y <- detectionMatrix
year2 <- year^2
month2 <- month^2

##save(list = c('nsite', 'nlist', 'aet', 'cwd', 'tmx', 'tmn', 'year', 'year2', 'month', 'month2', 'list_length', 'y'), file = '../data/zib_data.RData')

ind <- which(!is.na(y))
N <- length(ind)

aetflat <- aet[ind]
cwdflat <- cwd[ind]
tmxflat <- tmx[ind]
tmnflat <- tmn[ind]
yearflat <- year[ind]
year2flat <- year2[ind]
monthflat <- month[ind]
month2flat <- month2[ind]
list_lengthflat <- list_length[ind]
yflat <- y[ind]

siteMatrix <- array(rep(1:nsite, each=nlist), c(nlist, nsite))
siteIDflat <- siteMatrix[ind]

save(list = c(
         'N', 'nsite',
         'aetflat', 'cwdflat', 'tmxflat', 'tmnflat',
         'yearflat', 'year2flat', 'monthflat', 'month2flat',
         'list_lengthflat',
         'yflat',
         'siteIDflat'
),
     file = '../data/zib_data.RData')

q('no')



