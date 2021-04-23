# Image reconstruction by event-based 3D TOF FBP
This application, written in [R language](https://cran.r-project.org/), is devoted to image reconstruction for Jagiellonian PET (J-PET), using a combination of three independent kernel functions in image domain, similar to [multivariate kernel density estimation (KDE)](https://en.wikipedia.org/wiki/Multivariate_kernel_density_estimation "Wikipedia"), but with asymmetric definition.

## Prerequisites
Tested on Ubuntu 16.04 LTE with R version 3.6.3 installed, as well as on Scientific Linux CERN SLC release 6.10, with R version 3.4.1. The script operates with some shell commands
(for example, ```readlink -f <incomplete_path_to_file>```), but. basically, all of them are standard.

Additional R packages are required:
* [jsonlite](https://cran.r-project.org/web/packages/jsonlite/index.html)
* [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)
* [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html)
* [data.table](https://github.com/Rdatatable/data.table/wiki) (optional)
* [nat](https://cran.r-project.org/web/packages/nat/index.html) (optional - for reading ```.nrrd``` sensitivity map)
<!-- ** [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html) (optional - for plotting images) -->

There is no need to install anything. Just put directories ```modules/```, ```source/``` and ```cpp/``` into a single location, in order to match the paths properly.

## Input data and parameters
The application proceeds with [GOJA output format](https://github.com/JPETTomography/j-pet-gate-tools/tree/master/goja#goja-output) (16-column ASCII) as the input, but uses only first 8 columns: ```(hit1_x, hit1_y, hit1_z, time_of_hit1, hit2_x, hit2_y, hit2_z, time_of_hit2)```. Cartesian coordinates ```hit*_*``` must be in centimetres, time tags -- in picoseconds.

The parameters for the reconstruction, sensitivity map generation, as well as the geometry of the scanner, are stored in JSON format: it is self-describing and easy to undersand. Please refer to the ```examples/``` directory for the details.

## Map formats for corrections
Attenuation and sensitivity corrections need external files with dedicated maps. By default, they are stored as [```.rds``` R format](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/readRDS)  which contains [R list](https://www.r-tutor.com/r-introduction/list) of the following structure:

```
List of 4
 $ x        : num [1:223] -19.9 -19.7 -19.5 -19.3 -19.1 ...
 $ y        : num [1:223] -19.9 -19.7 -19.5 -19.3 -19.1 ...
 $ z        : num [1:191] -24.7 -24.5 -24.2 -24 -23.7 ...
 $ intensity: num [1:223, 1:223, 1:191] 0 0 0 0 0 0 0 0 0 0 ...
```

Alternatively, [```.nrrd```](http://teem.sourceforge.net/nrrd/format.html) can be used for sensitivity map (must be of the same dimensions as the reconstructed image), while plain ASCII is allowed for attenuation map, stored as follows:
```
"x"	"y"	"z"	"att.coeff"
-15.125	-11.725	-11.3	0
-14.925	-11.725	-11.3	0
...
15.075	-11.725	-11.3	0.096
-15.125	-11.525	-11.3	0.096
...
15.075	11.675	-11.3	0.096
-15.125	-11.725	-11.1	0.097
...
14.875	11.675	11.3	0
15.075	11.675	11.3	0
```
Here, a proper order should be preserved for al coordinates.

## Output formats
Available output file formats are ```.rds``` and ```.nrrd```. The former is native to R (has the same list structure as [described above](https://github.com/rshopa/event-based-tof-fbp#map-formats-for-corrections)), while the latter is also used in [J-PET MLEM](https://github.com/JPETTomography/j-pet-mlem) application. You need external app to view ```.nrrd``` images or dedicated libraries, e.g. [MRIcroGL](https://www.nitrc.org/projects/mricrogl), [Aliza](https://www.aliza-dicom-viewer.com/), [pynrrd](https://pypi.org/project/pynrrd/) etc. Store the outcome as ```.rds``` only if you are familiar with basic R: you need to run it interactively or write a dedicated script to view the images.

## Script architecture
The script is designed to work as quick as possible, hence [Rcpp](http://www.rcpp.org/) package is used and R environments are utilised for the [encapsulation](https://r6.r-lib.org/articles/Performance.html). The only exception has been made for JSON reader ```source/ReadFromJSON.R```, which uses Reference class.

## Usage
The implementation can be run in a single-thread mode or in parallel on multi-core nodes ([Open MPI](https://www.open-mpi.org/) is required). Assuming all script files are in the same place, execute the following to run the reconstruction:
```
$ cd <TOF_FBP_app_dir>
$ Rscript [--vanilla] source/LaunchReconstruction.R <path_to_json_params_file> [> log_file.log]
```
or (in multi-thread):
```
$ Rscript [--vanilla] source/LaunchReconstructionParallel.R <path_to_json_params_file> [> log_file.log]
```

Here, ```<path_to_json_params_file>``` denotes the path to JSON file with parameters (see ```examples/```). The option ```--vanilla``` [prevents Rscript from reading R history, profile, or environment files, as well as reloading data or objects from previous sessions](https://stat.ethz.ch/R-manual/R-devel/library/base/html/Startup.html).
