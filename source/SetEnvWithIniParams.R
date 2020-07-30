############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################
# this function "unrolls" constants/parameters into a dedicated environment
# to boost the performance and avoid too complex structures (named numeics, lists etc)
# the script may not DRY!

SetEnvWithIniParams <- function(Params.JSON){
  # temporary environment (will be deleted at the end)
  Dummies <- new.env()
  source("modules/Dummies.R", local = Dummies)
  # image space: dimensions for virtual scanner 
  r.pet     <- Params.JSON[["scanner"]][["radius"]]
  l.pet     <- Params.JSON[["scanner"]][["length"]]
  # kernel components profile parameters (sigmas etc)
  # --- FBP ---
  .delta.s  <- Params.JSON[["filters"]][["FBP"]][["delta-s"]]
  .dummy.FBP <- Dummies$createFBPdummy(Params.JSON[["filters"]][["FBP"]][["dummy-FFT-span"]],
                                       Params.JSON[["filters"]][["FBP"]][["omega-cut"]],
                                       Params.JSON[["filters"]][["FBP"]][["alpha"]],
                                       TOF=ifelse(Params.JSON[["filters"]][["FBP"]][["tau-regularisation"]]>=15,F,T),
                                       Tau.reg=Params.JSON[["filters"]][["FBP"]][["tau-regularisation"]])
  # --- TOF ---
  .sigma.TOF <- Params.JSON[["filters"]][["TOF"]][["sigma"]]
  if(Params.JSON[["filters"]][["TOF"]][["filter-type"]] == "inverse-gauss")
    .dummy.TOF <- Dummies$createInverseGaussDummy(.sigma.TOF, 
                                                  Params.JSON[["filters"]][["TOF"]][["dummy-FFT-span"]],
                                                  Params.JSON[["filters"]][["TOF"]][["nu-cut-intensity-factor"]], 
                                                  Alpha=Params.JSON[["filters"]][["TOF"]][["alpha"]]) 
  if(tolower(Params.JSON[["filters"]][["TOF"]][["filter-type"]]) == "cdf")
    .half.bin.TOF.size <- Params.JSON[["filters"]][["TOF"]][["half-bin-width"]]
  # --- Z ---
  .sigma.z <- Params.JSON[["filters"]][["Z"]][["sigma"]]
  if(Params.JSON[["filters"]][["Z"]][["filter-type"]] == "inverse-gauss")
    .dummy.z <- Dummies$createInverseGaussDummy(.sigma.TOF, 
                                                Params.JSON[["filters"]][["Z"]][["dummy-FFT-span"]],
                                                Params.JSON[["filters"]][["Z"]][["nu-cut-intensity-factor"]], 
                                                Alpha=Params.JSON[["filters"]][["Z"]][["alpha"]]) 
  if(tolower(Params.JSON[["filters"]][["Z"]][["filter-type"]]) == "cdf")
    .half.bin.z.size <- Params.JSON[["filters"]][["Z"]][["half-bin-width"]]

  # voxel size (transverse voxel size is the same for x and y)
  delta.xy <- .delta.s/Params.JSON[["filters"]][["FBP"]][["zoom"]]
  delta.z  <- 0.5*Params.JSON[["scanner"]][["length"]]/Params.JSON[["scanner"]][["number-of-rings"]]
  # axis and their length
  size.xy <- as.integer(floor(Params.JSON[["filters"]][["FBP"]][["zoom"]]*r.pet/.delta.s))*2L+1L 
  size.z  <- as.integer(Params.JSON[["scanner"]][["number-of-rings"]])*2L-1L
  axis.xy <- (-(size.xy-1)/2):((size.xy-1)/2)*delta.xy
  axis.z  <- seq_len(size.z)*delta.z-l.pet/2
  # set span for truncated ROR
  .semi.axis.TOF <- .sigma.TOF*Params.JSON[["filters"]][["TOF"]][["semi-axis-span-sigma-factor"]] #HalfBinSize*3.3 (11 for HP filters)
  .semi.axis.RL  <- .delta.s*Params.JSON[["filters"]][["FBP"]][["semi-axis-span-sigma-factor"]]  #9, 15 and more for apodised RL !
  .semi.axis.z   <- .sigma.z*Params.JSON[["filters"]][["Z"]][["semi-axis-span-sigma-factor"]] 
  # load attenuation if present
  if(!is.null(att.path <- Params.JSON[["corrections"]][["attenuation"]])){
    source("source/SetEnvForAttenuation.R", local = TRUE)
    att.environment <- SetEnvForAttenuation(att.path)
    rm(att.path, SetEnvForAttenuation)
  }
  rm(Params.JSON, Dummies)
  environment() # return environment
}