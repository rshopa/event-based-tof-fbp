############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

# A reader for .json files containing parameters
JSONReader <- setRefClass(
  "JSONReader",
  fields = list(
    input_path     = "character",
    output_name    = "character",
    scanner        = "list",
    filters        = "list",
    corrections    = "list",
    axes_ranges    = "numeric",
    # intended to be private:
    .initial_dir    = "character",
    .export_as_nrrd = "logical",
    .is_imported    = "logical"
  ),
  methods = list(
    # assign parameters from JSON file
    initialize = function(..., json_params_file = NULL, show_info = TRUE){
      .self[[".is_imported"]] <- FALSE
      # initial directory
      .self[[".initial_dir"]] <- getwd()
      if(!is.null(json_params_file)){
        callSuper(...)
        cat("JSON directory: ")
        cat(paste0(.self$.splitFileAndDir(json_params_file)[1],"\n"))
        # read list
        json_list <- read_json(json_params_file, simplifyVector = TRUE)
        # assign input/output file names and scanner
        .self$.assignIO(json_list)
        # change axes ranges if needed
        .self$.updateFOVspan(json_list)
        # import scanner parameters
        .self[["scanner"]] <- json_list[["virtual-scanner"]]
        # set as successfully imported
        .self[["filters"]] <- list(FBP = json_list[["FBP-filter"]],
                                   TOF = json_list[["TOF-filter"]],
                                   Z   = json_list[["Z-filter"]])
        if(!is.null(json_list[["post-filter"]]))
          .self[["filters"]][["post-filter"]] <- json_list[["post-filter"]] # post-filter
        validated_fields <- .self$.validateFBPfilter() &
                            .self$.validateVirtualScanner() &
                            .self$.validateTOFfilter(json_list) &
                            .self$.validateZfilter(json_list)   &
          .self$.nonNULLnonZeroLength(is.logical(.self[[".export_as_nrrd"]]))
        if(validated_fields){
          .self[[".is_imported"]] <- TRUE
          # filter out redundant fields in filters
          .self[["scanner"]]["CRT"]     <- NULL
          .self[["scanner"]]["sigma-z"] <- NULL
          for(n in names(.self[["filters"]])) 
            .self[["filters"]][[n]] <- .self$.filteredList(.self[["filters"]][[n]])
        }
      }
      # show parameters or the usage if needed, with possible memory issues
      if(show_info) .self$.usageInfo()
      # get back to initial directory
      setwd(.self[[".initial_dir"]])
    },
    
    # ----- dot (.) denotes 'private' methods -----
    # ----- IO -----
    # assigns IO parameters
    .assignIO = function(json_list){
      IO_list <- json_list[["input-output"]]
      .self[["input_path"]]  <- .self$.fullPath(IO_list[["input-data-path"]]) # full path
      .self[["output_name"]] <- IO_list[["output-file-name"]]                 # only a name
      # check if any attenuations
      .self[["corrections"]] <- list(attenuation=.self$.fullPath(IO_list[["attenuation-map-path"]]),
                                     sensitivity=.self$.fullPath(IO_list[["sensitivity-map-path"]]))
      # how to export the reconstructed image
      .self[[".export_as_nrrd"]] <- IO_list[["save-as-nrrd"]]
    },
    # shrinks FOV to the volume limited by [-axes_ranges, axes_ranges]
    .updateFOVspan = function(json_list){
      .self[["axes_ranges"]] <- as.numeric(rep(NA,3))
      if(!is.null(axes_ranges <- json_list[["restrict-FOV-axes-range"]]))
        invisible(lapply(seq_len(3), function(i){
          if(!is.null(axes_ranges[[i]])) .self[["axes_ranges"]][i] <<- axes_ranges[[i]]
        }))
    },
    .validateVirtualScanner = function(){
      return(Reduce("&",lapply(names(.self[["scanner"]]), function(n){
        if(n!="sigma-z" & n!="CRT") .self$.nonNULLnonZeroLength(.self[["scanner"]][[n]])
        else TRUE
      })))
    },
    .validateFBPfilter = function(){
      # assign Delta s = R*pi/N if null
      if(!is.numeric(.self[["filters"]][["FBP"]][["delta-s"]]))
        .self[["filters"]][["FBP"]][["delta-s"]] <- 
          pi*.self[["scanner"]][["radius"]]/.self[["scanner"]][["number-of-strips"]]
      # default zoom is set to 2
      if(!is.numeric(.self[["filters"]][["FBP"]][["zoom"]])) 
        .self[["filters"]][["FBP"]][["zoom"]] <- 2L
      # default dummy span is set to 35L
      if(!is.numeric(.self[["filters"]][["FBP"]][["dummy-FFT-span"]])) 
        .self[["filters"]][["FBP"]][["dummy-FFT-span"]] <- 35L
      # validate only these fields (others can be NULL)
      names_to_validate <- c("delta-s","alpha","omega-cut",
                             "tau-regularisation","zoom","semi-axis-span-sigma-factor")
      # combine and concatenate logically
      return(Reduce("&",lapply(names_to_validate, function(n){
        .self$.nonNULLnonZeroLength(.self[["filters"]][["FBP"]][[n]])
      })))
    },
    # validate TOF - from filter type, CRT?,  hallf bin...
    .validateGeneralFilter = function(FilterName){
      # assign first from ["gauss","CDF","inverse-gauss"] if unassigned
      .self[["filters"]][[FilterName]][["filter-type"]] <- 
        .self[["filters"]][[FilterName]][["filter-type"]][1]
      # CDF means you need to validate bin size (different from sigma)
      if(tolower(.self[["filters"]][[FilterName]][["filter-type"]]) == "cdf" &
         is.null(.self[["filters"]][[FilterName]][["half-bin-width"]]))
        .self[["filters"]][[FilterName]][["half-bin-width"]] <-
          .self[["filters"]][[FilterName]][["sigma"]]*sqrt(2*log(2))
      # validate only these fields (others can be NULL)
      names_to_validate <- c("filter-type","sigma","semi-axis-span-sigma-factor")
      # alpha is required for inverse filter only
      if(tolower(.self[["filters"]][[FilterName]][["filter-type"]]) == "inverse-gauss")
        names_to_validate <- c(names_to_validate,"alpha","dummy-FFT-span")
      # combine and concatenate logically
      return(Reduce("&",lapply(names_to_validate, function(n){
        .self$.nonNULLnonZeroLength(.self[["filters"]][[FilterName]][[n]])
      })))
    },
    # TOF-specific validate
    .validateTOFfilter = function(JSONlist){
      # try to estimate sigma from CRT (must be in ps)
      if(is.null(.self[["filters"]][["TOF"]][["sigma"]]))
        .self[["filters"]][["TOF"]][["sigma"]] <-
          JSONlist[["virtual-scanner"]][["CRT"]]/(4*sqrt(2*log(2)))*1e-10*299792458 # speed of light
      return(.self$.validateGeneralFilter("TOF"))
    },
    # validate Z (if present sigma-z)
    # Z-specific validate
    .validateZfilter = function(JSONlist){
      # try to estimate sigma from sigma_WLS, i.e. sigma_Z = sigma_WLS/sqrt(2)
      if(is.null(.self[["filters"]][["Z"]][["sigma"]])){
        if(is.numeric(JSONlist[["virtual-scanner"]][["sigma-z"]]))
          .self[["filters"]][["Z"]][["sigma"]] <-
            JSONlist[["virtual-scanner"]][["sigma-z"]]/sqrt(2)
        # try to estimate sigma from CRT (must be in ps)
        # effective speed of optical signal propagating through the scintillator
        # is 12.6 cm/ns = 0.0126 cm/ps
        else .self[["filters"]][["Z"]][["sigma"]] <-
            JSONlist[["virtual-scanner"]][["CRT"]]*0.0126/(4*sqrt(2*log(2)))
      }
      return(.self$.validateGeneralFilter("Z"))
    },
    # helpful function
    .nonNULLnonZeroLength = function(X) !is.null(X) & length(X)>0,
    # filtered list
    .filteredList = function(List){
      List[sapply(List, is.null)] <- NULL
      return(List)
    },
    # ---- INFO ----
    # returns usage if incorrectly imported
    .usageInfo = function(){
      if(.self[[".is_imported"]]){
        filters <- .self[["filters"]]
        scanner <- .self[["scanner"]]
        cat("----- Input parameters -----\n")
        cat(paste("Input file path:",          .self[["input_path"]],"\n"))
        cat(paste("Output file name:",         .self[["output_name"]],"\n"))
        cat("Effective (1-layer cylindrical) scanner:\n")
        cat(paste("\tradius:",              scanner[["radius"]],"cm\n"))
        cat(paste("\tlength:",              scanner[["length"]],"cm\n"))
        cat(paste("\tno. of strips:",       scanner[["number-of-strips"]],"\n"))
        cat("FBP filter:\n")
        cat(paste("\talpha:",              filters[["FBP"]][["alpha"]],"\n"))
        cat(paste("\tomega_cut:",          filters[["FBP"]][["omega-cut"]],"(cycles)\n"))
        cat(paste("\tTOF regulatisation:", filters[["FBP"]][["tau-regularisation"]],"\n"))
        for(Filt in c("TOF","Z")){
          cat(paste(Filt,"filter:\n"))
          cat(paste("\tfilter type:",        filters[[Filt]][["filter-type"]],"\n"))
          cat(paste("\tsigma:",              filters[[Filt]][["sigma"]],"cm\n"))
          if(tolower(filters[[Filt]][["filter-type"]]) == "cdf")
            cat(paste("\tbin width:", 2*filters[[Filt]][["half-bin-width"]],"cm\n"))
          if(tolower(filters[[Filt]][["filter-type"]]) == "inverse-gauss"){
            cat(paste("\talpha:", filters[[Filt]][["alpha"]],"\n"))
            cat(paste("\tcut-off intensity factor:", 
                      ifelse(is.null(nu_cut <- filters[[Filt]][["nu-cut-intensity-factor"]]),
                             "NULL", nu_cut),"\n"))
          }
          cat(paste("\t(+/-) kernel span / sigma:",filters[[Filt]][["semi-axis-span-sigma-factor"]],"\n"))
        }
        cat("Corrections:\n")
        for(n in names(.self[["corrections"]]))
          cat(paste0("\t",n,": ",ifelse(is.null(correction <- .self[["corrections"]][[n]]),
               "NULL", correction),"\n"))
        cat("\nAll parameters have successfully been imported.\n")
      } else cat("Something's gone wrong. Please check the input parameters.\n")
      cat("-----------------------------------------------\n\n")
    },
    # --- functions related to handling full path to file ---
    # returns full file path
    .fullPath = function(incomplete_path){
      if(is.null(incomplete_path)) return(NULL)
      else return(system(paste0("readlink -f ", incomplete_path), intern = TRUE))
    },
    # splits incomplete file path to full directory and filename
    .splitFileAndDir = function(file_path){
      full_path <- .self$.fullPath(file_path)
      directory <- dirname(full_path)
      # cut out directory path
      file_name <- unlist(strsplit(full_path, paste0(directory,"/")))[2]
      return(c(directory, file_name))
    }
  )
)
