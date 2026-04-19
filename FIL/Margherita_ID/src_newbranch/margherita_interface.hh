
template<class eos> 
class Pressure_Interface{
  public:
  static inline double press_eps_yle_ymu__beta_eq__rho_temp(double &eps, double &ye, double& ymu,
                                                double &rho, double &temp,
                                                typename eos::error_type &error){
    return eos::press_eps_yle_ymu__beta_eq__rho_temp(eps,ye,ymu, rho,temp,error);
  }
};

template<>
class Pressure_Interface<EOS_Tabulated>{
  public:
  static inline double press_eps_yle_ymu__beta_eq__rho_temp(double &eps, double &ye, double& ymu,
                                                double &rho, double &temp,
                                                typename EOS_Tabulated::error_type &error){

    typename Hot_Slice::error_t error_slice;
    const auto interp = Hot_Slice::get_extra_quantities(rho,error_slice);   
    eps = exp(interp[Hot_Slice::v_index::EPS]) + Hot_Slice::energy_shift;
    ye = interp[Hot_Slice::v_index::YE];
    temp = exp(interp[Hot_Slice::v_index::TEMP]);
    
    return exp(interp[Hot_Slice::v_index::PRESS]);
  }
};
