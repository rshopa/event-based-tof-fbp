############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

# main functions that actually make single-emission reconstruction
LORbyLORCalculator <- function(Params.JSON){
  # set initial environments and compile cpp functions
  source("source/SetEnvWithIniParams.R", local = TRUE)
  cat("Loading cpp functions... ")
  sourceCpp("cpp/IntensityEstimator.cpp", env = environment())
  cat("Done!\n")
  # import parameters from JSON list
  Ini.parameters <- SetEnvWithIniParams(Params.JSON) # constants
  # if non-null, attenuate
  .attenuate.or.not <- !is.null(Ini.parameters[["att.environment"]])
  # create parameters for each kernel component
  .kernel.params <- list(FBP = list(Ini.parameters[[".delta.s"]],
                                    Ini.parameters[[".dummy.FBP"]],
                                    as.integer((length(Ini.parameters[[".dummy.FBP"]])-1)/2)))
  # types that define the shape of kernel component
  .kernel.types  <- list(TOF="gauss", z="gauss") # default kernel types
  # verify for other types of kernels and update if needed
  invisible(lapply(c("TOF", "z"), function(n){
    .kernel.params[[n]] <<- list(Ini.parameters[[paste0(".sigma.",n)]])
    if(!is.null(dummy <- Ini.parameters[[paste0(".dummy.",n)]])){
      .kernel.types[[n]]       <<- "inverse-gauss"
      .kernel.params[[n]][[2]] <<- dummy
      .kernel.params[[n]][[3]] <<- as.integer((length(dummy)-1)/2)
    }
    if(!is.null(half.bin.size <- Ini.parameters[[paste0(".half.bin.",n,".size")]])){
      .kernel.types[[n]]       <<- "cdf"
      .kernel.params[[n]][[2]] <<- half.bin.size
    }
  }))
  # list made of unrolled Z-slices
  map.big.list <- lapply(seq_len(Ini.parameters[["size.z"]]), 
                         function(i) rep(0L, Ini.parameters[["size.xy"]]^2))
  # --- main function ---
  addLORToMap <- function(LOR){ # if FBP, normalize also by path
    # attenuation correction
    norm.att.constant <- (if(.attenuate.or.not)
      1/Ini.parameters[["att.environment"]][["estimateAttPathAndCoeff"]](LOR) else 1)
    if(!is.null(LOR.modified <- .getNewLOR(LOR))){
      # if LOR fits inside one Z-slice
      if(as.logical(LOR.modified[8]))
        # set all variables as separate parameters to boost the performance 
        # and avoid multiple subsetting from LOR.modified vector
        .addOneSlice(LOR.modified[1:3],
                     LOR.modified[4],LOR.modified[5],
                     LOR.modified[6],LOR.modified[7],
                     norm.att.constant)
      # else apply general function
      else .addMultipleSlices(LOR.modified[1:3],
                              LOR.modified[4],LOR.modified[5],
                              LOR.modified[6],LOR.modified[7],
                              norm.att.constant)
    }
  }
  # returns map as array (can be cut out from the whole FOV),
  # e.g. Axes.cuts[2] = 15.0 means Y will be restricted to [-15,15] range
  extract3DImage <- function(Axes.cuts = rep(NA,3), 
                             Verbose = FALSE, Purge.original = FALSE){
    n.slices <- Ini.parameters[["size.z"]]
    # default - complete FOV
    ranges.xyz <- rep(list(seq_len(Ini.parameters[["size.xy"]])), 2)
    ranges.xyz[[3]] <- seq_len(n.slices)
    # update ranges of indices if needed
    for(i in 1:2){
      if(!is.na(Axes.cuts[i])) 
        ranges.xyz[[i]] <- which(abs(Ini.parameters[["axis.xy"]]) <= Axes.cuts[i])
    }
    if(!is.na(Axes.cuts[3])){
      ranges.xyz[[3]] <- which(abs(Ini.parameters[["axis.z"]]) <= Axes.cuts[3])
      n.slices <- length(ranges.xyz[[3]])
    }
    # output image
    map.out <- array(0, dim=c(length(ranges.xyz[[1]]), length(ranges.xyz[[2]]), n.slices))
    # iterate by slices (list elements in map.big.list)
    for(k in seq_len(n.slices)){
      slice.z <- array(map.big.list[[ranges.xyz[[3]][k]]], 
                       dim = rep(Ini.parameters[["size.xy"]], 2)) # full transverse FOV
      map.out[,,k] <- slice.z[ranges.xyz[[1]], ranges.xyz[[2]]]
      if(Verbose) 
        cat(paste("Slice:",k,",\tz =", Ini.parameters[["axis.z"]][ranges.xyz[[3]][k]],"\r"))
    }
    # clear out initial map list (risky!)
    if(Purge.original) .purgeMapBigList(Reinitialize = FALSE)
    return(list(x         = Ini.parameters[["axis.xy"]][ranges.xyz[[1]]],
                y         = Ini.parameters[["axis.xy"]][ranges.xyz[[2]]],
                z         = Ini.parameters[["axis.z"]][ranges.xyz[[3]]],
                intensity = map.out))
  }
  # -----------------------------------------------------------------------
  # . (dot) denote private functions
  #----- TOF functions -----
  # annihilation point
  .detectAnnihilationPosition <- function(Hit.1, Hit.2, Time.1, Time.2){
    # speed of light here
    return(299792458e-10*(Time.1-Time.2)/sqrt(sum((Hit.1-Hit.2)^2))*
                         (Hit.2 - Hit.1)/2 + (Hit.1 + Hit.2)/2)
  }
  # transforms LORs into annihilation position and spherical angles:
  # c(X_src, Y_src, Z_src, Theta, Phi)
  .getNewLOR <- function(LOR){
    hit.1 <- LOR[1:3]
    hit.2 <- LOR[5:7]
    ann.pt <- .detectAnnihilationPosition(hit.1, hit.2, LOR[4], LOR[8])
    # outside FOV
    if((distanceBetweenPoints(hit.1, ann.pt) < Ini.parameters[[".semi.axis.TOF"]] & 
        distanceBetweenPoints(hit.2, ann.pt) < Ini.parameters[[".semi.axis.TOF"]]) | 
       abs(ann.pt[3]) > Ini.parameters[["l.pet"]]/2 | 
       distanceBetweenPoints(ann.pt[1:2], c(0,0)) > Ini.parameters[["r.pet"]]) return(NULL) 
    else {
      deltas <- hit.1 - hit.2 # c(Dx,Dy,Dz)
      # spherical coordinates
      theta.polar  <- acos(deltas[3]/sqrt(deltas[1]^2+deltas[2]^2+deltas[3]^2))
      phi.azimutal <- atan(deltas[2]/deltas[1])%%(2*pi)
      return(c(ann.pt, 
               cos(theta.polar), 
               sin(theta.polar),
               cos(phi.azimutal),
               sin(1e-15*(phi.azimutal==0) + phi.azimutal), # small positive non-zero (for sign())
               abs(deltas[3]) < Ini.parameters[["delta.z"]])) # whether LOR fits in one slice
    }
  }
  ########### MAPPING FUNCTIONS (use outer variables) ##################
  # conventions:
  # a - semi-axis along LOR (r)
  # b - semi-axis along X_r
  # c - semi-axis along Y_r
  # detects axis b
  .detectXrSemiAxis <- function(Cos.Q, Sin.Q){
    # abs to conserve positive! If ellipse impossible - get from TOF only
    if(abs(Ini.parameters[[".semi.axis.z"]]*Cos.Q) >= Ini.parameters[[".semi.axis.TOF"]]) 
      Ini.parameters[[".semi.axis.TOF"]]/abs(Sin.Q)
    else Ini.parameters[[".semi.axis.TOF"]]*Ini.parameters[[".semi.axis.z"]]*abs(Sin.Q) /
         sqrt(Ini.parameters[[".semi.axis.z"]]^2*(Sin.Q^4-Cos.Q^4) +
              Ini.parameters[[".semi.axis.TOF"]]^2*Cos.Q^2)
  }
  # detects span for X_r (rotated transverse axis)
  .detectXrSpan <- function(Z.r, a, b, Cos.Q, Sin.Q){
    A <- a^2*Cos.Q^2 + b^2*Sin.Q^2
    B <- 2*Z.r*Sin.Q*Cos.Q*(b^2-a^2)
    C <- Z.r^2*(b^2*Cos.Q^2 + a^2*Sin.Q^2)-a^2*b^2
    solution <- solveQuadraticEq(A,B,C)
    return(solution - Z.r*Sin.Q/Cos.Q) # subtract Delta Xr (distance from Z axis to LOR)
  }
  # detects span for Y_r (rotated transverse axis)
  .detectYrSpan <- function(Z.r, a, c, Cos.Q){
    # a-axis for the ellipse Y_r==Y_e,Z_e
    a <- a*Cos.Q
    wing.Y <- c*sqrt(1-Z.r^2/a^2) # produces NaNs if Z.r>=a
    if(is.nan(wing.Y)) return(NA) else return(wing.Y*c(-1,1))
  }
  # detects span for Z axis (no boundaries)
  .detectZSpan <- function(Z.src, Cos.Q) Z.src + Ini.parameters[[".semi.axis.TOF"]]*Cos.Q*c(-1,1)
  # finds span of a rectangle, biased along x and rotated by Phi
  .detectRectangleNewSpan <- function(X.span, Y.span, Cos.F, Sin.F){
    # - see .getNewLOR()
    Sign.factor <- sign(Sin.F*Cos.F)   
    return(list(X.span*Cos.F + Y.span*Sin.F * Sign.factor,     # x_range
                X.span*Sin.F + Y.span*Cos.F * Sign.factor))    # y_range
  }
  # assigns list of ranges for a particular Z-slice
  .findKernRangeXY <- function(Z, a, b, c, Cos.Q, Sin.Q, Cos.F, Sin.F, X.src, Y.src){
    span.Xr <- .detectXrSpan(Z, a, b, Cos.Q, Sin.Q) 
    span.Yr <- .detectYrSpan(Z, a, c, Cos.Q) 
    if(length(span.Xr) < 2 | is.na(span.Yr[1])) return(NA) # if outside ellipsoid or at the edge
    else {
      Output <- .detectRectangleNewSpan(span.Xr, span.Yr, Cos.F, Sin.F)
      Output[[1]] <- Output[[1]] + X.src # x_range
      Output[[2]] <- Output[[2]] + Y.src # y_range
      return(Output)
    }
  }
  # vectorised functions for rotation of axis XY by angle Phi with adjustment of the centre
  .adjustAndRotateX <- function(X, Y, Centre.XY, Cos.F, Sin.F) 
    (X-Centre.XY[1])*Cos.F + (Y-Centre.XY[2])*Sin.F
  .adjustAndRotateY <- function(X, Y, Centre.XY, Cos.F, Sin.F) 
    -(X-Centre.XY[1])*Sin.F + (Y-Centre.XY[2])*Cos.F
  # compose a factor for axis (IDs instead of coordinates)
  .matchedAxisFactor <- function(Axis, Span, Min.half.span){
    if(Span[1] > Span[2]) Span <- c(Span[2], Span[1]) # if descending
    if((Span[2] - Span[1]) < Min.half.span) 
      return(which(abs(Axis-Span[1]-Min.half.span) < (Min.half.span))) # case of thin span
    else return(which(inInterval(Axis, Span)))
  }
  # compose vector of Z's from Axis
  # adjusts centre to new Z coordinate along LOR
  .moveCentreOfAnnihilation <- function(Centre, Z, Cos.Q, Sin.Q, Cos.F, Sin.F){
    return(Centre+(Z-Centre[3])*c((c(Sin.Q*Cos.F,
                                     Sin.Q*Sin.F)/Cos.Q),
                                  1))
  }
  # creates XY-grid for the selected Z-slice, restricted by ROI ellipsoid 
  .createTransverseGrid <- function(XY.span, 
                                    Slice.centre, 
                                    A.elps, B.elps, C.elps,
                                    Cos.Q, Sin.Q, 
                                    Cos.F, Sin.F){
    # for multi-slice: modify b-axis for the ellipsoid:
    # a semi-axis in XY plane (not X_e,Y_e==Y_r!!!)
    if(!is.null(A.elps)) B.elps <- A.elps*B.elps/sqrt(A.elps^2*Cos.Q^2 + B.elps^2*Sin.Q^2)
    # create ID-grid for XY where truncated ROR is
    xy.IDs <- lapply(seq_len(2), function(i)
      .matchedAxisFactor(Ini.parameters[["axis.xy"]], XY.span[[i]], Ini.parameters[["delta.xy"]]/2)
    )
    # expand grid (IDs)
    grid.i <- rep(xy.IDs[[1]], length(xy.IDs[[2]]))
    grid.j <- rep(xy.IDs[[2]], each=length(xy.IDs[[1]]))
    if(length(grid.i)==0) return(NULL) # skip if no matches
    # expanded grid (Cartesians)
    grid.x <- Ini.parameters[["axis.xy"]][grid.i]
    grid.y <- Ini.parameters[["axis.xy"]][grid.j]
    # use points only within this ellipse,
    # i.e. inside.xy.ellipse <- x.r^2/(B.elps^2) + y.r^2/(C.elps^2) <= 1
    inside.xy.ellipse <- .adjustAndRotateX(grid.x,grid.y,Slice.centre,Cos.F,Sin.F)^2/(B.elps^2) + 
      .adjustAndRotateY(grid.x,grid.y,Slice.centre,Cos.F,Sin.F)^2/(C.elps^2) <= 1
    if((n.voxels <- sum(inside.xy.ellipse))==0) return(NULL) # no voxels to update
    else {
      row.IDs.for.map.big <- detectRow(grid.i[inside.xy.ellipse],
                                       grid.j[inside.xy.ellipse], Ini.parameters[["size.xy"]])
      # return as unnamed list: 1-2 - grid, 3 - n.voxels, 4 - row.IDs
      return(list(grid.x[inside.xy.ellipse], grid.y[inside.xy.ellipse],
                  n.voxels, row.IDs.for.map.big))
    }
  }
  # special case - LOR fits within one Z-slice 
  .addOneSlice <- function(Ann.pt, Cos.theta, Sin.theta, Cos.phi, Sin.phi, Norm.att.const){
    # find number of slice and avoid exception
    z.ID <- findVectorID(Ann.pt[3] - (Ini.parameters[["axis.z"]][1] - 
                                      Ini.parameters[["delta.z"]]/2), Ini.parameters[["delta.z"]])
    # slice ID should be within FOV
    if(z.ID > 0L & z.ID <= Ini.parameters[["size.z"]]){
      z.r  <- Ini.parameters[["axis.z"]][z.ID]
      # estimate XY span
      xy.span <- .detectRectangleNewSpan(Ini.parameters[[".semi.axis.TOF"]]*c(-1,1), 
                                         Ini.parameters[[".semi.axis.RL"]]*c(-1,1), 
                                         Cos.phi, Sin.phi)
      xy.span[[1]] <- xy.span[[1]] + Ann.pt[1]
      xy.span[[2]] <- xy.span[[2]] + Ann.pt[2]
      # create unnamed list: 1-2 - unrolled vectors for X and Y, 3 - n.voxels, 4 - row.IDs
      xy.slice.list <- .createTransverseGrid(xy.span, Ann.pt, 
                                             A.elps = NULL, # prevent update of B semi-axis
                                             B.elps = Ini.parameters[[".semi.axis.TOF"]],
                                             C.elps = Ini.parameters[[".semi.axis.RL"]],
                                             Cos.theta, Sin.theta, Cos.phi, Sin.phi)
      if(is.null(xy.slice.list)) break() # skip if empty
        intensities <- estimateIntensity_cpp(xy.slice.list[[1]],
                                             xy.slice.list[[2]],
                                             rep(z.r, xy.slice.list[[3]]),
                                             Ann.pt,
                                             Cos.theta, Sin.theta, Cos.phi, Sin.phi,
                                             .kernel.params[["TOF"]],
                                             .kernel.params[["z"]],
                                             .kernel.params[["FBP"]],
                                             .kernel.types[["TOF"]],
                                             .kernel.types[["z"]]
        ) * Norm.att.const # attenuation
        # update intensity in one slice
        map.big.list[[z.ID]][xy.slice.list[[4]]] <<-
          map.big.list[[z.ID]][xy.slice.list[[4]]] + intensities
    }
  }
  # general case
  .addMultipleSlices <- function(Ann.pt, Cos.theta, Sin.theta, Cos.phi, Sin.phi, Norm.att.const){
    # find span along Z
    z.span <- .detectZSpan(Ann.pt[3], Cos.theta)
    z.IDs  <- .matchedAxisFactor(Ini.parameters[["axis.z"]], z.span, Ini.parameters[["delta.z"]])
    # exception (outside FOV, length = 0)
    if(length(z.IDs)==0) break()
    z.grid.set <- Ini.parameters[["axis.z"]][z.IDs]
    # reassign for convenience
    a.elps <- Ini.parameters[[".semi.axis.TOF"]]
    b.elps <- .detectXrSemiAxis(Cos.theta, Sin.theta)
    c.elps <- Ini.parameters[[".semi.axis.RL"]]
    # --- loop over Z-slices ---
    for(slice.ID in seq_len(length(z.IDs))){
      z.r <- z.grid.set[slice.ID]
      # estimate semiaxis along X_r
      slice.centre <- .moveCentreOfAnnihilation(Ann.pt, z.r, Cos.theta, Sin.theta, Cos.phi, Sin.phi)
      # pass if outside truncated ROR
      if(distanceBetweenPoints(slice.centre, Ann.pt) < a.elps){
        # find span in transverse plane
        xy.span <- .findKernRangeXY(z.r-Ann.pt[3],
                                    a.elps, b.elps, c.elps,
                                    Cos.theta, Sin.theta, Cos.phi, Sin.phi,
                                    slice.centre[1], slice.centre[2]) 
        # create unnamed list: 1-2 - unrolled vectors for X and Y, 3 - n.voxels, 4 - row.IDs
        xy.slice.list <- .createTransverseGrid(xy.span, slice.centre, 
                                               a.elps, b.elps, c.elps,
                                               Cos.theta, Sin.theta, Cos.phi, Sin.phi)
        if(is.null(xy.slice.list)) break() # skip if empty
        intensities <- estimateIntensity_cpp(xy.slice.list[[1]],
                                             xy.slice.list[[2]],
                                             rep(z.r, xy.slice.list[[3]]),
                                             Ann.pt, 
                                             Cos.theta, Sin.theta, Cos.phi, Sin.phi,
                                             .kernel.params[["TOF"]],
                                             .kernel.params[["z"]],
                                             .kernel.params[["FBP"]],
                                             .kernel.types[["TOF"]],
                                             .kernel.types[["z"]]
        ) * Norm.att.const # attenuation
        # name for a slice in a list (MapBig is a list of unrolled 2D slices)
        slice.tag <- z.IDs[slice.ID]
        # update intensity in one slice
        map.big.list[[slice.tag]][xy.slice.list[[4]]] <<-
          map.big.list[[slice.tag]][xy.slice.list[[4]]] + intensities
      }
    }
  }
  # save RAM by removing or reinitialising (integers) output map
  .purgeMapBigList <- function(Reinitialize = TRUE){
    if(Reinitialize) map.big.list <<- lapply(seq_len(Ini.parameters[["size.z"]]), 
                                             function(i) rep(0L, Ini.parameters[["size.xy"]]^2))
    else map.big.list <<- NULL
    invisible(lapply(rep(F,5),gc)) # garbage collector
  }
  # remove redundant
  rm(SetEnvWithIniParams, Params.JSON)
  environment()
}