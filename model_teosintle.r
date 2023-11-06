library(tidyverse)

Dat <- read_csv("data/teosintle_maxent_input.csv")
bioclim <- terra::rast("data/teosintle_bioclim_raster.tif")
Dat <-Dat %>%
  bind_cols(terra::extract(x = bioclim,
                           y = Dat %>% select(lon, lat),
                           ID = FALSE))

me_g1 <- dismo::maxent(x = as(bioclim, "Raster") ,
                       p = Dat %>%
                         filter(ID == "g1") %>%
                         select(lon, lat) %>%
                         as.data.frame())
me_g2 <- dismo::maxent(x = as(bioclim, "Raster") ,
                       p = Dat %>%
                         filter(ID == "g2") %>%
                         select(lon, lat) %>%
                         as.data.frame())
me_all <- dismo::maxent(x = as(bioclim, "Raster") ,
                       p = Dat %>%
                         filter(ID == "all") %>%
                         select(lon, lat) %>%
                         as.data.frame())
dismo::plot(me_g1)
dismo::plot(me_g2)
dismo::plot(me_g3)


myfun <- function(x, p, args=NULL, path, ...) {
  
  stopifnot(dismo::maxent(silent=TRUE))
  
  cat("Hello\n")
  x <- cbind(p, x)
  cat("bye\n")
  x <- stats::na.omit(x)
  x[is.na(x)] <- -9999  # maxent flag for NA, unless changed with args(nodata= ), so we should check for that rather than use this fixed value.
  
  p <- x[,1]
  x <- x[, -1 ,drop=FALSE]
  
  factors <- NULL
  for (i in 1:ncol(x)) {
    if (inherits(x[,i], 'factor')) {
      factors <- c(factors, colnames(x)[i])
    }
  }
  
  if (!missing(path)) {
    path <- trim(path)
    dir.create(path, recursive=TRUE, showWarnings=FALSE)
    if (!file.exists(path)) {
      stop('cannot create output directory: ', path)
    }
    dirout <- path			
  } else {
    dirout <- .meTmpDir()
    f <- paste(round(runif(10)*10), collapse="")
    dirout <- paste(dirout, '/', f, sep='')
    dir.create(dirout, recursive=TRUE, showWarnings=FALSE)
    if (! file.exists(dirout)) {
      stop('cannot create output directory: ', f)
    }
  }
  
  pv <- x[p==1, ,drop=FALSE]
  av <- x[p==0, ,drop=FALSE]
  me <- new('MaxEnt')
  me@presence <- pv
  me@absence <- av
  me@hasabsence <- TRUE
  me@path <- dirout
  
  pv <- cbind(data.frame(species='species'), x=1:nrow(pv), y=1:nrow(pv), pv)
  av <- cbind(data.frame(species='background'), x=1:nrow(av), y=1:nrow(av), av)
  
  pfn <- paste(dirout, '/presence', sep="")
  afn <- paste(dirout, '/absence', sep="")
  write.table(pv, file=pfn, sep=',', row.names=FALSE)
  write.table(av, file=afn, sep=',', row.names=FALSE)
  
  mxe <- rJava::.jnew("mebridge")
  
  names(args) = NULL
  replicates <- .getreps(args) 
  args <- c("-z", args)
  
  if (is.null(factors)) {
    str <- rJava::.jcall(mxe, "S", "fit", c("autorun", "-e", afn, "-o", dirout, "-s", pfn, args)) 
  } else {
    str <- rJava::.jcall(mxe, "S", "fit", c("autorun", "-e", afn, "-o", dirout, "-s", pfn, args), rJava::.jarray(factors))
  }
  if (!is.null(str)) {
    stop("args not understood:\n", str)
  }
  
  
  if (replicates > 1) {
    
    mer <- new('MaxEntReplicates')
    d <- t(read.csv(paste(dirout, '/maxentResults.csv', sep='') ))
    d1 <- d[1,]
    d <- d[-1, ,drop=FALSE]
    dd <- matrix(as.numeric(d), ncol=ncol(d))
    rownames(dd) <- rownames(d)
    colnames(dd) <- d1
    mer@results <- dd
    f <- paste(dirout, "/species.html", sep='')
    html <- readLines(f)
    html[1] <- "<title>Maxent model</title>"
    html[2] <- "<CENTER><H1>Maxent model</H1></CENTER>"
    html[3] <- sub("model for species", "model result", html[3])
    newtext <- paste("using 'dismo' version ", packageDescription('dismo')$Version, "& Maxent version")
    html[3] <- sub("using Maxent version", newtext, html[3])
    f <- paste(dirout, "/maxent.html", sep='')
    writeLines(html, f)	
    mer@html <- f
    
    for (i in 0:(replicates-1)) {	
      mex <- me
      mex@lambdas <- unlist( readLines( paste(dirout, '/species_', i, '.lambdas', sep='') ) )
      
      f <- paste(mex@path, "/species_", i, ".html", sep='')
      html <- readLines(f)
      html[1] <- "<title>Maxent model</title>"
      html[2] <- "<CENTER><H1>Maxent model</H1></CENTER>"
      html[3] <- sub("model for species", "model result", html[3])
      newtext <- paste("using 'dismo' version ", packageDescription('dismo')$Version, "& Maxent version")
      html[3] <- sub("using Maxent version", newtext, html[3])
      f <- paste(mex@path, "/maxent_", i, ".html", sep='')
      writeLines(html, f)
      mex@html <- f
      mer@models[[i+1]] <- mex
      mer@models[[i+1]]@results <- dd[, 1+1, drop=FALSE]				
    }
    
    return(mer)
    
  } else {
    
    me@lambdas <- unlist( readLines( paste(dirout, '/species.lambdas', sep='') ) )
    d <- t(read.csv(paste(dirout, '/maxentResults.csv', sep='') ))
    d <- d[-1, ,drop=FALSE]
    dd <- matrix(as.numeric(d))
    rownames(dd) <- rownames(d)
    me@results <- dd
    
    f <- paste(me@path, "/species.html", sep='')
    html <- readLines(f)
    html[1] <- "<title>Maxent model</title>"
    html[2] <- "<CENTER><H1>Maxent model</H1></CENTER>"
    html[3] <- sub("model for species", "model result", html[3])
    newtext <- paste("using 'dismo' version ", packageDescription('dismo')$Version, "& Maxent version")
    html[3] <- sub("using Maxent version", newtext, html[3])
    f <- paste(me@path, "/maxent.html", sep='')
    writeLines(html, f)	
    me@html <- f
  }
  
  me
}


me_g1 <- myfun(x = Dat %>%
                         filter(ID == "g1") %>%
                         select(-ID, -lon, -lat) %>%
                         as.data.frame(),
                       p = rep(1, sum(Dat$ID == "g1")))

