############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

# generates sensitivity map from transverse .png image,
# expanding it along Z using analytical method by A. Strzelecki (2016)

require("methods", quietly = TRUE)    # essential for Rscript
require("jsonlite", quietly = TRUE)
require("EBImage", quietly = TRUE)

# Check if the parameters (.json file) have been given
# reassign 'args' variable
if(length(args <- commandArgs(trailingOnly = TRUE)) != 2){
  cat("\nUsage: ")
  cat("Rscript [--vanilla] GenerateSystemMap.R <json_params> <xy_sensitivity_png> \n\n")
  stop("No input file!")
}
# import params and data
params       <- read_json(normalizePath(args[1]))
path.to.png  <- normalizePath(args[2])
output.name  <- paste(dirname(path.to.png), paste0(params[["output-name"]],".rds"), sep="/")
sens.xy      <- readImage(path.to.png)
# avoid asymmetric transverse map
if(diff(dim(sens.xy@.Data))==0) n.px.png <- dim(sens.xy@.Data)[1] else stop("Asymmetric XY-map!")

# absolute span of XY-sensitivity image (in cm)
span.xy.png <- params[["png-xy-size"]]*n.px.png
# estimate XY-axis vector, its length and its span
delta.xy  <- pi*params[["scanner-radius"]]/params[["number-of-strips"]]/params[["zoom"]]
max.xy.id <- as.integer(floor(span.xy.png/delta.xy/2)) # max index 
axis.xy   <- (-max.xy.id:max.xy.id)*delta.xy
size.xy   <- length(axis.xy)
# same for Z-axis, but not related to transverse sensitivity
delta.z <- params[["scanner-length"]]/params[["number-of-rings"]]/2
size.z  <- as.integer(params[["number-of-rings"]]*2-1)
axis.z  <- seq_len(size.z)*delta.z - params[["scanner-length"]]/2

# resize XY-sensitivity (EBImage) and reassign as unrolled vector
unrolled.sens.xy <- c(resize(sens.xy@.Data, size.xy))
# unrolled vector of radii, i.e. sqrt(x^2+y^2)
unrolled.rs <- sqrt(rep(axis.xy,length(axis.xy))^2 + rep(axis.xy,each=length(axis.xy))^2)
# create sensitivity matrix
cat("Expanding matrix to 3D...\n")
sens.matrix.3D <- array(0, dim=c(rep(size.xy,2),size.z))
invisible(lapply(1:length(axis.z),function(Slice.No){
  cat(paste("Slice No:", Slice.No, "\r"))
  sens.matrix.3D[,,Slice.No] <<- as.numeric(lapply(1:length(unrolled.rs),function(i){
    # analytical approach from A.Strzelecki formulas (6.11-12)    
    denominator.sum  <- params[["scanner-radius"]] + unrolled.rs[i]
    denominator.diff <- params[["scanner-radius"]] - unrolled.rs[i]
    numerator.sum  <- params[["scanner-length"]]/2 + axis.z[Slice.No]
    numerator.diff <- params[["scanner-length"]]/2 - axis.z[Slice.No]
    arg.Fi.min <- max(c(-numerator.sum/denominator.diff,
                        -numerator.diff/denominator.sum))
    arg.Fi.max <- min(c(numerator.diff/denominator.diff,
                        numerator.sum/denominator.sum))
    return(unrolled.sens.xy[i]*(atan(arg.Fi.max)-atan(arg.Fi.min))/pi)
  }))
}))
# export to .rds
cat(paste("Done!\nSaving to file:", output.name,"\n"))
saveRDS(list(x=axis.xy, y=axis.xy, z=axis.z, intensity=sens.matrix.3D), output.name, version = 2)
