
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
##jagsGLMMdata <- readRDS('../orig_files/Data_JAGS_GLMM_v2.rds')
jagsGLMMdata <- readRDS('../orig_files/data_nimble_zib_v3.rds')
##potato <- readRDS('../orig_files/potato_psyllid_detection_dataset.rds')

createNamedObjectsFromList(jagsGLMMdata)

##y <- detectionMatrix
##year2 <- year^2     ## quadratic year effect
##month2 <- month^2   ## quadratic month effect
##year_list_length <- year * list_length   ## year:list_length interaction term

##ind <- which(!is.na(y))
##N <- length(ind)

##aetflat <- aet[ind]
##cwdflat <- cwd[ind]   ## dropped 'cwd' covarite in v2 of model
##tmxflat <- tmx[ind]
##tmnflat <- tmn[ind]
##yearflat <- year[ind]
##year2flat <- year2[ind]   ## dropped quadratic year effect
##monthflat <- month[ind]
##month2flat <- month2[ind]
##list_lengthflat <- list_length[ind]
##year_list_lengthflat <- year_list_length[ind]  ## year:list_length interaction term
##yflat <- y[ind]

##siteMatrix <- array(rep(1:nsite, each=nlist), c(nlist, nsite))
##siteIDflat <- siteMatrix[ind]

## change new siteID into indices from 1:nsite
siteID <- as.numeric(factor(siteID, levels = unique(siteID)))


save(list = c(
         'N', 'nsite',
         ## 'cwdflat',     ## dropped 'cwdflat' covariate
         'aet', 'tmx', 'tmn',
         'year', 'month', 'month2',
         'list_length',
         'year_list_length',  ## year:list_length interaction term
         'year_month', 'year_month2',
         'y',
         'siteID'
),
     file = '../data/zib_data_v3.RData')


## v4, with fixed effects for months

## changed the standardized 'month' variable into factors 1,2,3,...,12
monthfixed <- as.numeric(factor(month, levels = sort(unique(month))))

save(list = c(
         'N', 'nsite',
         'aet', 'tmx', 'tmn',
         'year',
         ## 'month', 'month2',   ## dropped standardized month variables
         'monthfixed',           ## new fixed effects for months
         'list_length',
         'year_list_length',
         ## 'year_month', 'year_month2',   ## dropped standardized month variables
         'y',
         'siteID'
),
     file = '../data/zib_data_v4_monthfixed.RData')


qqq


