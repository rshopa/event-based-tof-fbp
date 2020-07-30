############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################
# functions for attenuation correction taken from an article in a 'strange' journal:
# Ruru Li et al. Journal of Information & Computational Science 7:14 (2010) 3125

SetEnvForAttenuation <- function(Att.file){
  # import attenuation map as ASCII or .rds
  # (binary means map is stored in .rds)
  cat("Loading attenuation map... ")
  attenuation.map <- importImage(Att.file)
  cat("Done!\nRedefine map as factors... ")
  # modify the data into unrolled integers (it is faster!) and levels
  # can be expensive for a large map!
  att.map.unrolled      <- as.factor(c(attenuation.map[[4]]))
  attenuation.map[[4]]  <- NULL
  levels.of.attenuation <- as.numeric(levels(att.map.unrolled))
  att.map.unrolled      <- as.integer(att.map.unrolled)
  cat("Done!\n")
  # length of an unrolled attenuation map
  .size.unrolled <- length(att.map.unrolled)
  # dimension of FOV subvolume where phantom is drawn
  .phantom.box.dims <- as.integer(sapply(seq_len(3), function(i) 
                                  length(attenuation.map[[i]])))
  # voxel size -- must be symmetric!
  .d.xyz <- sapply(seq_len(3), function(i) round(mean(diff(attenuation.map[[i]])),3))
  # Errors if found
  .errors <- character()
  if(sum(.phantom.box.dims %% 2) != 0)
    .errors <- c(.errors, paste0("ERROR: all dimensions (",
                                 paste(.phantom.box.dims, collapse=" "),
                                 ") must be even. Please redefine attenuation map!"))
  if(diff(diff(.d.xyz)) > 1e-3) # prevent small difference (weird R behaviour)
    .errors <- c(.errors, paste0("ERROR: asymmetric voxel (",
                                 paste(as.numeric(.d.xyz),collapse=" "),
                                 ") is not supported. Please redefine attenuation map!"))
  else .d.xyz <- as.numeric(.d.xyz[1]) # leave only one element and round
  if(length(.errors)!=0) stop(paste(.errors, collapse = "\n  "))
  rm(attenuation.map) # remove redundant

  # --- main function ---
  # estimates attenuation path in cm and estimated attenuation:
  # returns c(AttPath,exp(-AttCoeff))
  estimateAttPathAndCoeff <- function(LOR, Calculate.path=FALSE){
    special.return <- if(Calculate.path) c(0,1) else 1 # for exceptions - c(path, exp(-AttCoeff))
    px.1 <- findPixelID(LOR[1:3], .d.xyz)
    px.2 <- findPixelID(LOR[5:7], .d.xyz)
    l.px      <- distanceBetweenPoints(px.1, px.2)
    l.ranges  <- .getLambdaRanges(px.1, px.2, .phantom.box.dims, l.px)
    if(sum(is.finite(l.ranges))==0) return(special.return) # exception
    l.min.max <- .getMinMaxLambdas(l.ranges)
    # exclusion for ranges (lambda_min - lambda_max)
    if(diff(l.min.max)<=0) return(special.return) # exception
    # entering voxel
    e.vox   <- .enteringVoxel(px.1, px.2, l.ranges, l.min.max[1], .phantom.box.dims, l.px)
    m.first <- .estimateFirstVoxel(e.vox, px.1, px.2, l.px, l.min.max[1])
    # return special (no attenuation)
    if(sum(abs(m.first)<=.phantom.box.dims/2)!=3) return(special.return) # exception
    else {
      path <- .combineLambdas(px.1, px.2, l.px, e.vox, l.min.max) 
      att.coeff <- 0 # init
      if(Calculate.path) att.path  <- 0 # optional
      starting.vox <- m.first + as.integer(.phantom.box.dims/2L)+1L
      # Exclusions for short or zero lengths
      if(is.null(n.steps <- length(path[[1]]))) return(special.return)
      if(n.steps==2L) path[[2]] <- t(as.matrix(path[[2]])) # prevent error: incorrect number of dimensions
      # c(path, AttCoeff) - if Calculate.path=T
      for(i in seq_len(n.steps)){
        # don't adjust first voxel
        if(i>1) starting.vox[path[[2]][i-1L,1]] <- 
                starting.vox[path[[2]][i-1L,1]] + path[[2]][i-1L,2]
        map.row <- detect3DRow(starting.vox[1],
                               starting.vox[2],
                               starting.vox[3], .phantom.box.dims[1], .phantom.box.dims[2])
        if(map.row>=1 & map.row <= .size.unrolled){
          # select one of the coefficients from phantom materials (levels) 
          # according to ID in att.map.unrolled
          map.element  <- as.numeric(levels.of.attenuation[att.map.unrolled[map.row]])
          path.element <- path[[1]][i]*.d.xyz
          # attenuation coefficient
          att.coeff <- att.coeff + map.element*path.element
          # optional: exclude air/no attenuation
          if(Calculate.path) att.path <- att.path + (map.element!=0)*path.element
        }
      }
      if(!Calculate.path) return(exp(-att.coeff)) else return(c(att.path,exp(-att.coeff)))
    }
  }
  
  # Functions assumed to be PRIVATE: denoted by (.) dots
  ######################################################
  # ------- specific function from the Ref.:
  # ------- Ruru Li et al.
  .lambdaParam <- function(I,X.1,X.2,L) (I-X.1)/(X.2-X.1)*L
  # Axis ranges for X,Y,Z (by row)
  .getLambdaRanges <- function(Px.1, Px.2, Phantom.dims, L.px){
    n.xyz <- Phantom.dims/2
    out.matrix <- array(0, dim=c(3,2))
    for(i in seq_len(3)){
      rng <- .lambdaParam(-n.xyz[i]:n.xyz[i], Px.1[i], Px.2[i], L.px)
      out.matrix[i,] <- c(min(rng), max(rng)) # faster than range()
    }
    return(out.matrix)
  }
  # three lambdas for max and min (WARNING: voxels are symmetrical!)
  .getMinMaxLambdas <- function(Lambda.rng){
    # c(lmin,lmax)
    return(c(max(Lambda.rng[is.finite(Lambda.rng[,1]),1], na.rm = T),
             min(Lambda.rng[is.finite(Lambda.rng[,2]),2], na.rm = T)))
  }
  # 1D: detect the first pixel where the ray enters the phantom box
  # cannot be vectorised
  .enteringPixel <- function(X.1, X.2, X.lambda.min, Lambda.min, X.phantom.span, L.px){
    if(X.1 == X.2) i.enter <- X.1
    else if(X.1 < X.2){
      if(Lambda.min == X.lambda.min) i.enter <- -X.phantom.span/2+1L
      else i.enter <- ceiling(X.1 + Lambda.min/L.px*(X.2 - X.1))
    } else {
      if(Lambda.min == X.lambda.min) i.enter <- X.phantom.span/2-1L
      else i.enter <- floor(X.1 + Lambda.min/L.px*(X.2 - X.1))
    }
    return(as.integer(i.enter))
  }
  # voxel where the ray enters the phantom box
  .enteringVoxel <- function(Px.1, Px.2, Lambda.rng, Lambda.min, Phantom.dims, L.px){
    e.voxel <- rep(NA,3)
    for(i in seq_len(3))
      e.voxel[i] <- .enteringPixel(Px.1[i], Px.2[i], 
                                   Lambda.rng[i,1], Lambda.min, Phantom.dims[i], L.px)
    return(e.voxel)
  }
  # first voxel to estimate m_first (intersecting length)
  .estimateFirstVoxel <- function(Enter.vox, Px.1, Px.2, L.px, Lambda.min){
    min.enter.lambda <- min(.lambdaParam(Enter.vox, Px.1, Px.2, L.px), na.rm = T)
    return(as.integer(floor(Px.1+(min.enter.lambda+Lambda.min)/(2*L.px)*(Px.2-Px.1))))
  }
  # see Ruru Li et al. Journal of Information & Computational Science 7:14 (2010) 3125
  # fo the next two functions
  .buildSetOfLambdas <- function(Enter.vox, Px.1, Px.2, L.px, Lambda.rng){
    # estimate first lambda(0) - on entering phantom box
    entering.lambdas <- .lambdaParam(Enter.vox, Px.1, Px.2, L.px)
    # eq(26-27)
    delta.lambdas <- abs(L.px/(Px.2-Px.1))
    # max indices (denoted as l,m,n)
    i.max <- as.integer(floor((Lambda.rng[2] - entering.lambdas)/delta.lambdas))
    i.max[is.na(i.max) | is.infinite(i.max) | i.max<0] <- 0L # avoid NaN/NA
    out <- list()
    for(i in seq_len(3)){
      lambda.i <- entering.lambdas[i] + (0:i.max[i])*delta.lambdas[i]
      out[[i]] <- lambda.i[!is.na(lambda.i) & inInterval(lambda.i, Lambda.rng)]
    }
    return(out)
  }
  
  .combineLambdas <- function(Px.1, Px.2, L.px, Enter.vox, Lambda.rng){
    l.set <- .buildSetOfLambdas(Enter.vox, Px.1, Px.2, L.px, Lambda.rng)
    i.delta <- (Px.1 < Px.2)*1L + (Px.1 >= Px.2)*(-1L) # eq.30
    step.matrix <- matrix(c(rep(c(1L,i.delta[1]), length(l.set[[1]])),
                            rep(c(2L,i.delta[2]), length(l.set[[2]])),
                            rep(c(3L,i.delta[3]), length(l.set[[3]]))), ncol=2, byrow=T)
    # return ordered 
    l.set <- sort.int(unlist(l.set),index.return = T,method = "quick")
    # unnamed list
    return(list(c(diff(c(Lambda.rng[1], l.set[[1]], Lambda.rng[2]))),
                  step.matrix[l.set[[2]],]))
  }
  # --- remove redundant ---
  rm(Att.file)
  environment() # return environment
}
