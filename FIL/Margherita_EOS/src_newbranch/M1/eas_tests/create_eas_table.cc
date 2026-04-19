#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

// Let's have some debug output
#define STANDALONE
#define COLDTABLE_SETUP

//#define DEBUG
#include "../FIL_M1_headers.h"
#include "../../3D_Table/readtable.cc"
#include "../M1.hh"

#include "eas_write_table.cc"


std::string const table_name = std::string("/Users/carlomusolino/numrel/EOS_Tables/EOSs/DD2/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
std::string const fname = "DD2_eas.h5";

constexpr int num_rho = 256;
constexpr int num_ye = 60;
constexpr int num_temp = 180;

bool const do_BB_correction = true ;

double const rhomin{1e-10}, rhomax{5e-03};
double const tempmin{1e-02}, tempmax{150};
double const yemin{0.01}, yemax{0.6};

template<size_t N>
std::array<double,N> generate_linspace(double const &min,double const& max);

template<size_t N>
std::array<double,N> generate_logspace(double const &min,double const& max);



int main(){

  //std::string("/zhome/academic/HLRS/xfp/xfpmusol/EOS_Tables/DD2/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");

  constexpr bool use_energy_shift = true;
  constexpr bool recompute_mu_nu = false;
  EOS_Tabulated::readtable_scollapse(table_name.c_str(), use_energy_shift, recompute_mu_nu);

  std::array<size_t,4> npoints{3,num_ye,num_temp,num_rho}; // order matters here !
  
  std::array<double,num_rho> lrho = generate_logspace<num_rho>(rhomin,rhomax);
  std::array<double,num_ye> ye = generate_linspace<num_ye>(yemin,yemax);
  //std::array<double,num_temp> ltemp = generate_logspace<num_temp>(tempmin,tempmax);
  std::array<double,num_temp> ltemp = generate_linspace<num_temp>(tempmin,tempmax);

  //for(auto const& r: lrho) std::cout << pow(10,r) << "\t";
  //std::cout << std::endl ;
  
  auto kappa_a = new double[M1_eas_tables::NUM_SPECIES][num_ye][num_temp][num_rho];
  auto kappa_s = new double[M1_eas_tables::NUM_SPECIES][num_ye][num_temp][num_rho];
  auto kappa_n = new double[M1_eas_tables::NUM_SPECIES][num_ye][num_temp][num_rho];
  auto Q = new double[M1_eas_tables::NUM_SPECIES][num_ye][num_temp][num_rho];
  auto R = new double[M1_eas_tables::NUM_SPECIES][num_ye][num_temp][num_rho];

  
  for (int nn_rho = 0; nn_rho < num_rho; nn_rho++)
    for (int nn_ye = 0; nn_ye < num_ye; nn_ye++)
      for (int nn_temp = 0; nn_temp < num_temp; nn_temp++) {
	
	double const rhoL = pow(10,lrho[nn_rho]);
	double const tempL = ltemp[nn_temp];
	double const yeL = ye[nn_ye];

	// calc eas with Margherita at yeL tempL rhoL 
	const auto tau_init= compute_analytic_opacity(rhoL);
	std::array<double,3> tau_n {{tau_init,tau_init,tau_init}};
	auto F= Fugacities<double>(rhoL,tempL,yeL, std::move(tau_n));
	auto eas = EAS<double>(F,do_BB_correction);
	eas.calc_eas();

	//save data
	Q[NUE][nn_ye][nn_temp][nn_rho] = eas.Q_cgs(NUE);
	Q[NUA][nn_ye][nn_temp][nn_rho] = eas.Q_cgs(NUA);
	Q[NUX][nn_ye][nn_temp][nn_rho] = eas.Q_cgs(NUX);

	R[NUE][nn_ye][nn_temp][nn_rho] = eas.R_orig(NUE);
	R[NUA][nn_ye][nn_temp][nn_rho] = eas.R_orig(NUA);
	R[NUX][nn_ye][nn_temp][nn_rho] = eas.R_orig(NUX);

	kappa_a[NUE][nn_ye][nn_temp][nn_rho] = eas.kappa_a_cgs(NUE);
	kappa_a[NUA][nn_ye][nn_temp][nn_rho] = eas.kappa_a_cgs(NUA);
	kappa_a[NUX][nn_ye][nn_temp][nn_rho] = eas.kappa_a_cgs(NUX);

	kappa_s[NUE][nn_ye][nn_temp][nn_rho] = eas.kappa_s_cgs(NUE);
	kappa_s[NUA][nn_ye][nn_temp][nn_rho] = eas.kappa_s_cgs(NUA);
	kappa_s[NUX][nn_ye][nn_temp][nn_rho] = eas.kappa_s_cgs(NUX);

	kappa_n[NUE][nn_ye][nn_temp][nn_rho] = eas.kappa_n_cgs(NUE);
	kappa_n[NUA][nn_ye][nn_temp][nn_rho] = eas.kappa_n_cgs(NUA);
	kappa_n[NUX][nn_ye][nn_temp][nn_rho] = eas.kappa_n_cgs(NUX);
	  
	}

  auto logrho_ptr = (new double[num_rho]);
  auto logtemp_ptr =(new double[num_temp]);
  auto ye_ptr = (new double[num_ye]);
  auto dummy = (new double[3]);

  for (int i=0; i<4; i++) dummy[i]=i;
  for (int i = 0; i < num_rho; ++i) logrho_ptr[i] = lrho[i];
  for (int i = 0; i < num_temp; ++i) logtemp_ptr[i] = ltemp[i];
  for (int i = 0; i < num_ye; ++i) ye_ptr[i] = ye[i];


  M1_eas_tables::eas_table<4,decltype(Q)> tab(std::vector<decltype(Q)>{Q,R,kappa_a,kappa_s,kappa_n},
					      std::move(npoints),
					      std::move(dummy),
					      std::move(ye_ptr),
					      std::move(logtemp_ptr),
					      std::move(logrho_ptr)
					      );

  tab.write_table(fname);

  delete[] Q; delete[] R; delete[] kappa_a; delete[] kappa_s; delete[] kappa_n;
  delete[] logrho_ptr; delete[] logtemp_ptr; delete[] ye_ptr;
}

template<size_t N>
std::array<double,N> generate_linspace(double const &min,double const& max){
  std::array<double,N> out;
  for(int i=0;i<N;i++){
    out[i] = min + i * (max-min)/N;
  }
  return out;
}

template<size_t N>
std::array<double,N> generate_logspace(double const &min,double const& max){
  std::array<double,N> out;
  double lmin{log10(min)}, lmax{log10(max)};
  for(int i=0;i<N;i++){
    out[i] = lmin + i * (lmax-lmin)/N;
  }
  return out;
}
