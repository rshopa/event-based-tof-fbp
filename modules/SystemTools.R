############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################
# system tools for importing and exporting data

# imports image as .rds, .nrrd or ASCII
importImage <- function(Path.to.image, Header=TRUE){ # can contain header variable
  if(grepl("\\.rds$", ext <- tolower(Path.to.image))) return(readRDS(Path.to.image))
  else return(if(grepl("\\.nrrd$", tolower(Path.to.image)))
                nat::read.nrrd(Path.to.image) # try reading as .nrrd binary
              else import3DImageASCII(Path.to.image, Header)) # header is required
}
# allows using data.table library
importASCII <- function(Path.to.ASCII, Header=FALSE){
  rows <- if("data.table" %in% rownames(installed.packages()))
    as.matrix(data.table::fread(Path.to.ASCII, header=Header))
  else as.matrix(read.table(Path.to.ASCII, header=Header))
  attributes(rows)["dimnames"] <- NULL # will be processed faster
  return(rows)
}
# for corrections
import3DImageASCII <- function(Path.to.ASCII, Header=TRUE){
  img <- importASCII(Path.to.ASCII, Header)
  img <- list(x = sort(unique(img[,1])), 
              y = sort(unique(img[,2])), 
              z = sort(unique(img[,3])), intensity=img[,4])
  img[["intensity"]] <- array(img[["intensity"]], dim=c(length(img[["x"]]), 
                                                        length(img[["y"]]), length(img[["z"]])))
  return(img)
}
# export image in .rds or .nrrd format
exportImage <- function(Image.3D, File.name, Rds=TRUE){
  if(Rds) saveRDS(Image.3D, file = paste0(File.name,".rds"), version = 2) # for older R before 3.5.0
  else {
    out.path <- paste0(File.name,".nrrd")
    exportNRRD(Image.3D[[4]], out.path, "little")     # only map, no axes
  }
}
# exports image as .nrrd binary
exportNRRD <- function(Map3D, Output.path, Endian=.Platform$endian){
  data.file <- tools::file_path_sans_ext(Output.path)
  header <- list(type = "float", encoding = "raw", endian = Endian,
                 dimension = length(dim(Map3D)), sizes = dim(Map3D),
                 datafile = basename(data.file))
  # writes header
  cat("NRRD0004\n", file = Output.path)
  # write.nrrd.header(header, Output.path) 
  for (n in names(header)){
    f <- header[[n]]
    cat(paste0(n, ": ", 
               if(length(f) > 1) paste(f, collapse = " ") else f, 
               "\n"), file = Output.path, append = TRUE)
  }
  cat("\n", file = Output.path, append = TRUE)
  # writes data
  fc <- file(data.file, open = "wb")
  writeBin(as.vector(Map3D, mode = "numeric"), fc, size = 4L, # 4 bytes, float
           endian = Endian)
  close(fc)
}
