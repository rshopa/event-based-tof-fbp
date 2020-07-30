# Setting the parameters for the event-based 3D TOF FBP
The parameters for the reconstruction, sensitivity map generation, as well as the geometry of the scanner, are set in [JSON format](https://www.json.org/json-en.html). In order to define them properly, you need to understand the key principles of the event-based TOF FBP.

## Files needed for sensitivity correction
Sensitivity correction is used to take into account the imperfections of the scanner. It is generally done in a straightforward way:

<img src="https://latex.codecogs.com/svg.latex?
\mathrm{Img}_{corrected}=\mathrm{Img}_{ini}/\mathrm{Sens.map}."/>

You can create such map e.g. using Monte Carlo simulation in [GATE](http://www.opengatecollaboration.org/) or other software. The app allows the following extensions for sensitivity map <img src="https://latex.codecogs.com/svg.latex?\mathrm{Sens.map}"/> : ```.rds``` (R native), ```.nrrd``` or any ASCII format (not tested). The dimensions of ```.nrrd``` map must be the same as for the output image <img src="https://latex.codecogs.com/svg.latex?\mathrm{Img}_{ini}"/>. For ```.rds``` and ASCII, it is assumed axis vectors are explicitly given (refer to the main README for the example), while <img src="https://latex.codecogs.com/svg.latex?\mathrm{Sens.map}"/> covers the same or larger field-of-viev (FOV) as <img src="https://latex.codecogs.com/svg.latex?\mathrm{Img}_{ini}"/>.

There is also a tool to create full 3D sensitivity map using 2D transverse cross-section only, simulated beforehand. An exemplary ```.png``` files generated by Monte Carlo are given in ```examples/sensitivity_map/```: for ideal cylindric scanner, big barrel J-PET and modular 2-layer J-PET. You need an additional JSON file which constitutes parameters that define pixel size ```png-xy-size```, as well as the size of the voxel (given indirectly, as the geometry set for the effective ideal scanner, which TOF FBP parameters are estimated from). For instance, in ```ideal_jpet.json```:

```
"output-name": "sensitivity_map",
"scanner-radius": 43.73,
"scanner-length": 50.0,
"number-of-strips": 384,
"number-of-rings": 96,
"zoom": 2,
"png-xy-size": 0.2
```
Two features here are used in the same way as in [STIR](http://stir.sourceforge.net/) framework: ```number-of-rings``` artificially splits J-PET scintillator strips into axial "rings", meaning that the number of voxels will amount to ```2*number-of-rings – 1``` in axial direction. ```zoom``` is used to increase the sampling in transverse direction compared to projection space, so that XY-size of the voxel will be <img src="https://latex.codecogs.com/svg.latex?\Delta xy = \pi R/(N_{strips}\cdot \mathrm{zoom})"/>.

To run a generation of 3D sensitivity map using XY ```.png```-dummy, execute the following:
```
$ Rscript [--vanilla] source/GenerateSystemMap.R <json_params> <xy_sensitivity_png>
```
This will save the map as ```.rds``` in the same directory as 2D image.

**Caution: the examples presented here don't cover the whole transverse FOV!** You need to manually restrict XY-span of the image using parameter ```restrict-FOV-axes-range``` in the main JSON file (see below).

<img src="https://github.com/rshopa/event-based-tof-fbp/blob/master/examples/images/FOV_sizes.png" width="500">

## JSON parameters for the reconstruction