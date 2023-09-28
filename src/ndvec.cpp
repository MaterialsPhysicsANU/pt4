
#include "ndvec.h"
#include "quaternion.h"
#include <cstdint>
#include <algorithm>
#include <tuple>

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include "tiffio.h"

#define PACK(c0, c1, c2, c3) \
    (((uint32_t)(uint8_t)(c0) << 24) | \
    ((uint32_t)(uint8_t)(c1) << 16) | \
    ((uint32_t)(uint8_t)(c2) << 8) | \
    ((uint32_t)(uint8_t)(c3)))

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

#pragma warning(push)
#pragma warning(disable : 4838)//if you are writing an array where size_t becomes an issue you may have other issues
//this function is not thread safe as netCDF isn't
int write_ncf(const char* fp, const void* data, const size_t3 size, int NC_TYPE) {
    #pragma omp critical (ncf)
    {
    int retval;
    const int ndims = 3;
    /* When we create netCDF variables and dimensions, we get back an
        * ID for each one. */
    int ncid, x_dimid, y_dimid, z_dimid, varid;
    int dimids[ndims];
    
    auto [nx,ny,nz] = size;

    static int grid_size[] = {nx, ny, nz};

    /* Create the file. NC_CLOBBER overwrites file*/
    if ((retval = nc_create(fp, NC_NETCDF4|NC_CLOBBER, &ncid))){
        ERR(retval);
    }
    /* Define the dimensions. NetCDF will hand back an ID for each. */
    if ((retval = nc_def_dim(ncid, "z", nz, &x_dimid))){
        ERR(retval);
    }
    if ((retval = nc_def_dim(ncid, "y", ny, &y_dimid))){
        ERR(retval);
    }
    if ((retval = nc_def_dim(ncid, "x", nx, &z_dimid))){
        ERR(retval);
    }
    
    /* The dimids array is used to pass the IDs of the dimensions of
        * the variable. */
    dimids[0] = x_dimid;
    dimids[1] = y_dimid;
    dimids[2] = z_dimid;
    
    if ((retval = nc_def_var(ncid, "data", NC_TYPE, ndims,
                    dimids, &varid))){
        ERR(retval);
    }
        if ((retval = nc_put_att_int(ncid,  NC_GLOBAL, "total_grid_size_xyz",
                     NC_INT, ndims, grid_size))){
        ERR(retval);
    }
    /* End define mode no more metadata. */
    if ((retval = nc_enddef(ncid))){
        ERR(retval);
    }
    /* Write data to the file.*/

    switch (NC_TYPE)
    {
    case NC_FLOAT:
        retval = nc_put_var_float(ncid, varid, (float*)data);
        break;
    case NC_UINT:
        retval = nc_put_var_uint(ncid, varid, (uint32_t*)data);
        break;
    }
    if (retval){
        ERR(retval);
    }

    if ((retval = nc_close(ncid))){
        ERR(retval);
    }
    return 0;
    }
}
#pragma warning(pop)

int vol2ncf(const char* fp, const volume<float> vol) {
    return write_ncf(fp, (void*)vol.data.data(), vol.size, NC_FLOAT);
}

#pragma warning(push)
#pragma warning(disable : 4838)//if you are writing an array where size_t becomes an issue you may have other issues
//this function is not thread safe as netCDF isn't
volume<float> ncf2vol(const char* fp) {
    #pragma omp critical (ncf)
    {
    int retval;
    int ncid;
   /* Open the file. */
   if ((retval = nc_open(fp, NC_NETCDF4|NC_NOWRITE, &ncid))){
      ERR(retval)
   }

    int ndims_in, nvars_in, ngatts_in, unlimdimid_in;
    if ((retval = nc_inq(ncid, &ndims_in, &nvars_in, &ngatts_in,
            &unlimdimid_in))){
      ERR(retval);
    }
    if (ndims_in != 3 || nvars_in != 1){
        printf("nfc file wrong shape\n");
        exit(2);
    }

    int varid;
    if ((retval = nc_inq_varid(ncid, "data", &varid))){
      ERR(retval);
    }

    //get size
    int x_dimid, y_dimid, z_dimid;
    char dimname[NC_MAX_NAME];
    size_t nx,ny,nz;
    if ((retval = nc_inq_dimid(ncid, "x", &x_dimid))){      ERR(retval);    }
    if ((retval = nc_inq_dim(ncid, x_dimid, dimname, &nx))){      ERR(retval);    }

    if ((retval = nc_inq_dimid(ncid, "y", &y_dimid))){      ERR(retval);    }
    if ((retval = nc_inq_dim(ncid, y_dimid, dimname, &ny))){      ERR(retval);    }

    if ((retval = nc_inq_dimid(ncid, "z", &z_dimid))){      ERR(retval);    }
    if ((retval = nc_inq_dim(ncid, z_dimid, dimname, &nz))){      ERR(retval);    }
    

    //read in data
    auto vol = volume<float>({nx,ny,nz});
    if ((retval = nc_get_var_float(ncid, varid, vol.data.data()))){
      ERR(retval);
    }

    if ((retval = nc_close(ncid))){
      ERR(retval);
    }

    return vol;
    }
}
#pragma warning(pop)

//reads tiff and gives slice
template<typename T>
slice<T> tiff2slice(const char* fp){
    TIFF* tiff = TIFFOpen(fp, "r");

    slice<float> slice0;
    if (tiff) {
        uint32_t width, height;
        size_t npixels;
        float* raster;

        TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &height);
        slice0 = slice<float>({width,height});

        // Write the float data to the TIFF file
        for (size_t row = 0; row < slice0.size[1]; ++row) {
            if (TIFFReadScanline(tiff, &slice0[{0,row}], row, 0) < 0) {
            }
        }

        TIFFClose(tiff);
    } else {
        out.log(ERR) << "Error opening TIFF file for reading:" << fp << std::endl;
    }

    return slice0;
}

template slice<float> tiff2slice<float>(const char* fp);

//saves slice to tiff
template<typename T>
int slice2tiff(const char* fp, slice<T>& slice0 ){
    TIFF* tiff = TIFFOpen(fp, "w");

    if (tiff) {
        // Set TIFF tags to indicate that this is a float image
        TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
        TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, sizeof(T) * 8);
        TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
        TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, slice0.size[0]);
        TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, slice0.size[1]);
        TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

        // Write the float data to the TIFF file
        for (size_t row = 0; row < slice0.size[1]; ++row) {
            if (TIFFWriteScanline(tiff, &slice0[{0,row}], row, 0) < 0) {
                return 1;
            }
        }
        TIFFClose(tiff);
    } else {
        out.log(ERR) << "Error opening TIFF file for writing:" << fp << std::endl;
        return 1;
    }
    return 0;
}

//saves slice to tiff
template<typename T>
int slice2tiff(const char* fp, slice<std::array<T,3>>& slice0){
    TIFF* tiff = TIFFOpen(fp, "w");

    T absmax = T();
    for(std::array<T,3>& v : slice0.data){
        absmax = std::max(
                            absmax, std::max(
                            std::abs(v[0]), std::max(
                            std::abs(v[1]), 
                            std::abs(v[2])
                        )));
    }
    T absmaxi = T(UINT8_MAX)/(T(2)*absmax);

    auto rgba_slice = slice<std::array<uint8_t,4>>(slice0.size);
    for(size_t i = 0; i < slice0.data.size(); ++i){
        const std::array<T,3> v = slice0.data[i];
        rgba_slice.data[i] = {
            static_cast<uint8_t>(v[0]*absmaxi + T((UINT8_MAX+1)/2)),
            static_cast<uint8_t>(v[1]*absmaxi + T((UINT8_MAX+1)/2)),
            static_cast<uint8_t>(v[2]*absmaxi + T((UINT8_MAX+1)/2)),
            255u
        };
    }


    if (tiff) {
        // Set TIFF tags to indicate that this is a float image
        TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
        TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, sizeof(uint8_t) * 8);
        TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 4);
        TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, rgba_slice.size[0]);
        TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, rgba_slice.size[1]);
        TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
        // TIFFSetField(tiff, TIFFTAG_EXTRASAMPLES, EXTRASAMPLE_ASSOCALPHA);//this should be in here but it is bugged in libtiff

        // Write the float data to the TIFF file
        for (size_t row = 0; row < rgba_slice.size[1]; ++row) {
            if (TIFFWriteScanline(tiff, &rgba_slice[{0,row}], row, 0) < 0) {
                return 1;
            }
        }
        TIFFClose(tiff);
    } else {
        out.log(ERR) << "Error opening TIFF file for writing:" << fp << std::endl;
        return 1;
    }
    return 0;
}


template int slice2tiff<float>(const char* fp, slice<float>& slice0);
template int slice2tiff<float>(const char* fp, slice<float3>& slice0);
