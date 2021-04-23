############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

# main file for the reconstruction
require("methods", quietly = TRUE)    # essential for Rscript
require("jsonlite", quietly = TRUE)
require("Rcpp", quietly = TRUE)

# extract directory of a script
args <- commandArgs(trailingOnly = FALSE)
# solution from here: stackoverflow.com/a/32016824/538403
script.dir  <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))

# Check if the parameters (.json file) have been given
# reassign 'args' variable
if(length(args <- commandArgs(trailingOnly = TRUE)) != 1){
  cat("\nUsage: ")
  cat("Rscript [--vanilla] LaunchReconstruction.R <json_params>\n\n")
  stop("Improper input!")
}
# will be used later
current.dir <- getwd()
params.dir  <- normalizePath(args[1])

setwd(paste0(script.dir,"/.."))
source("modules/GeneralMath.R")
source("modules/SystemTools.R")
source("source/ReadFromJSON.R")
source("source/PostRecoTools.R")

# load calculator which uses external Rcpp functions
cat("Loading event-by-event calculator and cpp functions... ")
source("source/LORbyLORCalculator.R")
sourceCpp("cpp/IntensityEstimator.cpp", env = environment())
cat("Done!\n")

cat(paste("Script directory:",script.dir,"\n"))
cat(paste("Current_directory:",getwd(),"\n"))

reconstruction.params <- JSONReader(json_params_file=params.dir)
FBP.TOF.calculator    <- LORbyLORCalculator(reconstruction.params)
post.calculator       <- PostRecoTools(reconstruction.params)

cat("Importing data file... ")
B2B.emissions <- if(grepl("\\.rds$", input.path <- reconstruction.params[["input_path"]],ignore.case = T))
  readRDS(input.path) else importASCII(input.path, FALSE)
cat("Done!\n")

cat("--------------------------\nInitiating reconstruction:\n")
Start.time <- Sys.time()
for(iter in seq_len(nrow(B2B.emissions))){ #
  FBP.TOF.calculator$addLORToMap(B2B.emissions[iter,])
  if(iter%%100==0) cat(paste("Events:",iter,"\r"))
}
cat("Done! Total time:\n")
Sys.time() - Start.time # Time difference
cat(paste("Total events:", iter,"\n"))

# extract image
output.3D.image  <- FBP.TOF.calculator$extract3DImage(reconstruction.params[["axes_ranges"]],
                                                      Verbose = TRUE, Purge.original = TRUE)
# apply sensitivity correction
if(!is.null(reconstruction.params[["corrections"]][["sensitivity"]])){
  output.3D.image <- post.calculator$applySensitivityMap(output.3D.image)
  cat("Sensitivity correction has successfully been applied.\n")
}
# apply post-filter
if(!is.null(post.filter <- reconstruction.params[["filters"]][["post-filter"]])){
  cat("Applying post-filter:\n")
  cat(paste("\tType:", post.filter[["type"]],"\n"))
  cat(paste("\tRadius:", post.filter[["radius"]],"\n"))
  cat(paste("\tBall-shaped:", post.filter[["ball-shaped"]],"\n"))
  output.3D.image[[4]] <- 
    post.calculator$applyConvolutionalFilter(output.3D.image[[4]],
                                             post.filter[["radius"]],
                                             post.filter[["ball-shaped"]],
                                             post.filter[["type"]],
                                             "median", # boundaries method
                                             Verbose=TRUE)
}

# save the output
export.path <- paste(dirname(reconstruction.params[["input_path"]]),
                     reconstruction.params[["output_name"]],sep="/")
cat(paste0("\r\nSaving to file:", export.path,
          if(reconstruction.params[[".export_as_nrrd"]]) ".nrrd" else ".rds", "\n"))
exportImage(output.3D.image, export.path, Rds = !(reconstruction.params[[".export_as_nrrd"]]))
setwd(current.dir)
