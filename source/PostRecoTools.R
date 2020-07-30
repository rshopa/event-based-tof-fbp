############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

# loads functions to be applied after the reconstruction (if given)
PostRecoTools <- function(Params.JSON){
  # read sensitivity if present and assign functions to update
  if(!is.null(Params.JSON[["corrections"]][["sensitivity"]])){
    sensitivity.map <- importImage(Params.JSON[["corrections"]][["sensitivity"]], Header = TRUE)
    # updates intensity using sensitivity
    applySensitivityMap <- function(Image.3D){
      if(.validateSensMapSize(Image.3D)){
        if(is.array(sensitivity.map)) Image.3D[[4]] <- Image.3D[[4]] / sensitivity.map
        else {
          dim.sens <- dim(sensitivity.map[[4]])
          dim.img  <- dim(Image.3D[[4]])
          diff.dims <- as.integer((dim.sens - dim.img)/2)
          ranges.updated <- lapply(seq_len(3), function(i)
            seq_len(dim.sens[i])[seq_len(dim.img[i]) + diff.dims[i]])
          Image.3D[[4]] <- Image.3D[[4]] / sensitivity.map[[4]][ranges.updated[[1]], 
                                                                ranges.updated[[2]], 
                                                                ranges.updated[[3]]]
          
        }
        Image.3D[[4]][is.nan(Image.3D[[4]])] <- 0 # escape 0/0 case
        return(Image.3D)
      } else {
        cat("ERROR: cannot apply sensitivity map. Please check dimensions. Return unchanged.\n")
        return(Image.3D)
      }
    }
    # whether the map can be applied
    .validateSensMapSize <- function(Image.3D, Tol.voxel = 1e-3){ # tolerance for voxel size
      # same dimensions (default)
      if(is.array(sensitivity.map)) 
        return(Reduce("&",dim(sensitivity.map) == dim(Image.3D[[4]])))
      else {
        # no. of voxels must be odd for each axis and the same or larger than for the image
        return(Reduce("&", dim(sensitivity.map[[4]]) >= dim(Image.3D[[4]])) &
                 Reduce("&", sapply(1:3, function(i) 
                   length(sensitivity.map[[i]]) %% 2 == 1 &
                     abs(mean(diff(sensitivity.map[[i]])) - mean(diff(Image.3D[[i]]))) < Tol.voxel))
        )
      }
    }
  }
  # post-filter if present
  if(!is.null(Params.JSON[["filters"]][["post-filter"]])){
    # load library and compile functions
    cat("Loading cpp functions for post-filtering... ")
    require("RcppArmadillo")
    sourceCpp("cpp/ConvolutionalFilers3D.cpp", env = environment())
    cat("Done!\n")
    # applies convolutional median or mean average filter
    applyConvolutionalFilter <- function(Img,
                                         Radius,
                                         Rounded,
                                         Filter.type,
                                         Edge.method,
                                         Fill.by = NULL,
                                         Verbose = FALSE){
      dims <- dim(Img)
      Img.extended <- .expandImageForKernel(Img, Radius, Edge.method, Fill.by, Verbose)
      Img.extended <- filterImage_cpp(Img.extended, Radius, Rounded, Filter.type, Verbose)
      if(Verbose) cat("\rDone!\n")
      return(Img.extended[(Radius+1L):(Radius+dims[1]),
                          (Radius+1L):(Radius+dims[2]), (Radius+1L):(Radius+dims[3])])
    }
    # convolutional filter needs larger image to handle edges
    .expandImageForKernel <- function(Img, Radius, Method, Fill.by=NULL, Verbose=FALSE){
      if(missing(Method)) Method <- "zeros";
      Method <- match.arg(Method, c("zeros", "neighbour", "mean", "median", "fill"))
      if(Verbose) cat(paste("Method chosen to fill boundaries:", 
                            ifelse(Method=="fill", paste0(Method," by ",Fill.by), Method),"\n"))
      dims <- dim(Img)
      new.dims <- dims + 2*Radius
      out.image <- array(NA, dim = new.dims) # new image
      scales.ijk <- lapply(seq_len(3), function(n) (Radius + 1L):(Radius + dims[n]))
      out.image[scales.ijk[[1]], scales.ijk[[2]], scales.ijk[[3]]] <- Img
      # unroll IDs for axes and find voxels 'outside' the initial image
      ids.tab <- expand.grid(1:new.dims[1], 1:new.dims[2], 1:new.dims[3], KEEP.OUT.ATTRS = FALSE)
      outside.ini.img.factor <- ids.tab[[1]]<=Radius | ids.tab[[2]]<=Radius | ids.tab[[3]]<=Radius |
        ids.tab[[1]]>(Radius+dims[1]) | ids.tab[,2]>(Radius+dims[2]) | ids.tab[,3]>(Radius+dims[3]) 
      # switch between methods 
      if(Method != "neighbour"){
        # clear central volume exept of edging voxels (assign NA's and cut)
        out.image[scales.ijk[[1]], 
                  scales.ijk[[2]], 
                  scales.ijk[[3]]][2:(dims[1]-1),2:(dims[2]-1),2:(dims[3]-1)] <- NA
        # various cases
        if(Method=="mean")   out.image[outside.ini.img.factor] <- mean(out.image[!is.na(out.image)])
        if(Method=="median") out.image[outside.ini.img.factor] <- median(out.image[!is.na(out.image)])
        if(Method=="fill")   out.image[outside.ini.img.factor] <- Fill.by
        if(Method=="zeros")  out.image[outside.ini.img.factor] <- 0
        # reassign from ini image central voxels
        out.image[scales.ijk[[1]], scales.ijk[[2]], scales.ijk[[3]]] <- Img
      } else invisible(apply(ids.tab[outside.ini.img.factor,],1,function(row){
          i <- row[1]; j <- row[2]; k <- row[3]
          out.image[i,j,k] <<- out.image[scales.ijk[[1]][findPixelIDByScale(i,scales.ijk[[1]])],
                                         scales.ijk[[2]][findPixelIDByScale(j,scales.ijk[[2]])],
                                         scales.ijk[[3]][findPixelIDByScale(k,scales.ijk[[3]])]]}))
      return(out.image)
    }
  }
  # remove redundant
  rm(Params.JSON)
  environment()
}