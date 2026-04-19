#pragma once

#include <H5Cpp.h>
#include <vector>
#include <string>
#include <memory>
#include <stdlib.h>
#include <array>


namespace M1_eas_tables{


  enum SPECIES { NUE=0, NUA, NUX, NUM_SPECIES };
  enum VARS { Q =0 , R, K_A, K_S, K_N, NUM_VARS };
  enum DIMS { SPECIES=0,Y_E,TEMP,RHO, NUM_DIMS };
  
  
  template <size_t dim,typename var_t>
  class eas_table{
  private: 
    using dim_t = double*;
    using tab_t = std::vector<var_t> ;
    std::array<size_t,dim> npoints;
    std::array<dim_t,dim> x;
    tab_t alltables;
    
    void write_1D_dataset(H5::H5File* file, const H5std_string & DSET_NAME, hsize_t const *dimf, size_t const var_index) ;
    void write_ND_dataset(H5::H5File* file, const H5std_string & DSET_NAME, hsize_t const *dimf, size_t const var_index) ;
    
  public:
    
    void write_table(std::string const & fname) ;
    template<typename tab_t, typename ... x_t>
    eas_table(tab_t && alltables_, std::array<size_t,dim>&& num_points_, x_t && ... x_ ): npoints(std::move(num_points_)), alltables(std::move(alltables_)) {
      int n = 0;
      int tmp[] = {(x[n] = std::forward<x_t>(x_), ++n)...};
      (void)tmp;
    }
  };
}
