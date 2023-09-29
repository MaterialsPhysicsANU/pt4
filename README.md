# pt4
A program to describe 4D-CT phantoms. 4D-phantoms ('.pt4' files) can be rendered to volumes or projection data. The scanning tracjectory and noise model can be customised when simulating projections.

This project requires:
* libTIFF   ( https://gitlab.com/libtiff/libtiff            , http://www.libtiff.org                              )
* netCDF    ( https://github.com/Unidata/netcdf-c           , https://www.unidata.ucar.edu/software/netcdf/       )
* exprtk    (                                               , http://www.partow.net/programming/exprtk/index.html )
* simplecpp ( https://github.com/danmar/simplecpp           ,                                                     )
* GNU Bison (                                               , https://www.gnu.org/software/bison/ )


## Usage

The program reads in a `.pt4` and then generates a time-series of volumes `-v` or projections `-p` of the volume in the X,Y,Z axes in the ncf format. The program can be called with:
```
pt4.exe [-p] [-v] (path to .pt4 without extension)
```
and will output to the `(path to .pt4 without extension)/` directory.

The projection data will output floating point TIFF files, as the floating points will generally be outside the range [0,1] it will not be viewable in standard image viewers. This is okay.

Please check out the phantoms folder for a list of examples.

## Dependencies

When cloning with `git clone --recurse-submodules` the following dependancies should clone automatically into the `./lib/` directory:
* libTIFF
* simplecpp
Otherwise they can be manually downloaded to `./lib/`

### netCDF
 netCDF's C interface should be installed in a location known to cmake. The build process does not automatically build netCDF, so prepackaged binaries are sufficient.

### exprtk
 Please download exprtk and extract to `./lib/exprtk`

## Building

The program can be built using CMake. It is recommended to build the release version due to the significant speedup. The following instructions show how to build on Windows with CMake and Ninja. The build process assumes that GNU Bison is installed in a WSL instance. 
```
 mkdir build
 cd build
 cmake .. -DCMAKE_BUILD_TYPE=Release -G"Ninja"
 ninja
```

For questions relating to pt4 please contact: [todo]@anu.edu.au
