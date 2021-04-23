############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

# main file for the reconstruction
require("methods",   quietly = TRUE)    # essential for Rscript
require("jsonlite",  quietly = TRUE)
require("Rcpp",      quietly = TRUE)
require("parallel",  quietly = TRUE)

# extract directory of a script
args <- commandArgs(trailingOnly = FALSE)
# solution from here: stackoverflow.com/a/32016824/538403
script.dir  <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))

# Check if the parameters (.json file) have been given
# reassign 'args' variable
if(length(args <- commandArgs(trailingOnly = TRUE)) != 1){
  cat("\nUsage: ")
  cat("Rscript [--vanilla] LaunchReconstructionParallel.R <json_params>\n\n")
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

cat(paste("Script directory:",script.dir,"\n"))
cat(paste("Current_directory:",getwd(),"\n"))

reconstruction.params <- JSONReader(json_params_file=params.dir)
post.calculator       <- PostRecoTools(reconstruction.params)

cat("Importing data file... ")
B2B.emissions <- if(grepl("\\.rds$", input.path <- reconstruction.params[["input_path"]], ignore.case = T))
  readRDS(input.path) else importASCII(input.path, FALSE)
cat("Done!\n")

cat("--------------------------\nInitiating reconstruction:\n")

########### SPLIT THE DATA INTO CHUNKS AND RUN IN PARALLEL ############
n.cores    <- detectCores()
n.events   <- nrow(B2B.emissions)
chunk.size <- ceiling(n.events / n.cores)
ids.ranges <- lapply(1 : as.integer( ceiling(n.events / chunk.size) ), 
                     function(i){
                       min.i <- (i - 1) * chunk.size + 1
                       max.i <- min(n.events, i * chunk.size)
                       return(min.i:max.i)
                     })

cat("Creating temporary dir... ")
system("mkdir -p __TEMP_DIR/",ignore.stdout = T)
cat("done!\n")
cat(paste0("Launching in parallel (",n.cores," cores)...\n")) # ??? change later ???
cat("Loading calculator and cpp functions... ")
source("source/LORbyLORCalculator.R")
sourceCpp("cpp/IntensityEstimator.cpp", env = environment())
cat("Done!\n")

Start.time <- Sys.time()
invisible(mclapply(ids.ranges, function(rng){
  # cut chunk and create a separate calculator for each one (temp invisible!!!)
  B2B.chunk <- B2B.emissions[rng,]
  invisible(capture.output(FBP.TOF.calculator <- LORbyLORCalculator(reconstruction.params)))
  # run loop
  for(iter in seq_len(nrow(B2B.chunk))) FBP.TOF.calculator$addLORToMap(B2B.chunk[iter,])
  # extract image
  output.3D.image  <- FBP.TOF.calculator$extract3DImage(reconstruction.params[["axes_ranges"]],
                                                        Verbose = FALSE, Purge.original = TRUE)
  # save .rds to temp dir (using the first digit in IDs range)
  name.to.save <- paste0(sprintf("__TEMP_DIR/chunk_%09d", rng[1]),".rds") 
  saveRDS(output.3D.image, file = name.to.save)
  cat(paste(name.to.save, "completed.\n"))
}, mc.cores = n.cores, mc.cleanup = TRUE))
rm(ids.ranges)
invisible(lapply(1:5, function(i) gc())) # slightly awkward!!!

cat("Reconstruction completed! Total time:\n")
Sys.time() - Start.time # Time difference

# extract image
cat("Merging chunks into one image...")
fls <- system("ls __TEMP_DIR/*.rds", intern = TRUE)
output.3D.image <- readRDS(fls[1])
for(f in fls[2:length(fls)]){
  chunk.img <- readRDS(f)
  output.3D.image[["intensity"]] <- output.3D.image[["intensity"]] + chunk.img[["intensity"]]
}
rm(chunk.img)
cat("Done!\n")
# purge temporal directory
system("rm -rf __TEMP_DIR/", ignore.stdout = T)

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
