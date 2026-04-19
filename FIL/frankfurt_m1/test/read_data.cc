#include <hdf5.h>

#include "test_implicit_rootfinder.h"

static double readvar(std::string const& vname, std::string const& fname){
  hid_t file_id, dspace_id, dset_id, group_id ;
  herr_t h5err ;
  file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT ) ;

  dset_id = H5Dopen( file_id, vname.c_str(), H5F_ACC_RDONLY ) ;
  double value ;
  h5err = H5Dread( dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value ) ;
  return value ;
}
