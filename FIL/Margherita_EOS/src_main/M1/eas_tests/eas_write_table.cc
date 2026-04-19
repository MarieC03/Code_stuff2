#pragma once
#include<iostream>
#include<string>

#include<H5Cpp.h>

#include"eas_write_table.hh"

namespace M1_eas_tables{

template<size_t dim,typename T>
void eas_table<dim,T>::write_1D_dataset(H5::H5File* file, const H5std_string & DSET_NAME, hsize_t const *dimf, size_t const var_index) {
  using namespace H5;
  FloatType datatype(PredType::NATIVE_DOUBLE);
  DataSpace dataspace( 1, dimf) ;
  DataSet dataset = file->createDataSet( DSET_NAME, datatype, dataspace);
  dataset.write(x[var_index], PredType::NATIVE_DOUBLE);
}

  template<size_t dim,typename T>
  void eas_table<dim,T>::write_ND_dataset(H5::H5File* file, const H5std_string & DSET_NAME, hsize_t const * dimf, size_t const var_index)
{
  using namespace H5;  
  const size_t RANK = dim;
  FloatType datatype(PredType::NATIVE_DOUBLE);
  DataSpace dataspace(RANK,dimf);
  DataSet dataset = file->createDataSet(DSET_NAME, datatype, dataspace);
  dataset.write(alltables[var_index], PredType::NATIVE_DOUBLE);
}


template<size_t dim,typename T>
void eas_table<dim,T>::write_table(std::string const& fname){
  using namespace H5;

  FloatType datatype(PredType::NATIVE_DOUBLE);
  
  const H5std_string FILE_NAME(fname.c_str()) ;
  const H5std_string DSET_K_A_NAME("absorption_opacity");
  const H5std_string DSET_K_S_NAME("scattering_opacity");
  const H5std_string DSET_K_N_NAME("number_opacity");

  const H5std_string DSET_Q_NAME("emissivity");
  const H5std_string DSET_R_NAME("number_emissivity");

  const H5std_string DSET_RHO_NAME("pointsrho");
  const H5std_string DSET_TEMP_NAME("pointstemp");
  const H5std_string DSET_Y_E_NAME("pointsye");

  H5File file(FILE_NAME, H5F_ACC_TRUNC );

  hsize_t dimf = npoints[RHO];
  write_1D_dataset(&file,DSET_RHO_NAME,&dimf,RHO);
  std::cout << "Done with rho" << std::endl;

  dimf = npoints[TEMP];
  write_1D_dataset(&file,DSET_TEMP_NAME,&dimf,TEMP);
  std::cout << "Done with temp" << std::endl;
  
  dimf = npoints[Y_E];
  write_1D_dataset(&file,DSET_Y_E_NAME,&dimf,Y_E);
  std::cout << "Done with ye" << std::endl;
  
  hsize_t dimf_tab[dim];
  dimf_tab[0] = 3;
  for(int i=1; i<dim; i++) {
    std::cerr << npoints[i] << std::endl; 
    dimf_tab[i] = npoints[i];
  }
  write_ND_dataset(&file,DSET_K_A_NAME,dimf_tab,K_A);
  write_ND_dataset(&file,DSET_K_S_NAME,dimf_tab,K_S);
  write_ND_dataset(&file,DSET_K_N_NAME,dimf_tab,K_N);

  write_ND_dataset(&file,DSET_Q_NAME,dimf_tab,Q);
  write_ND_dataset(&file,DSET_R_NAME,dimf_tab,R);
  
}
}
