//
//
// Quick and dirty test to benchmark the inversion process
// and check it for consistency. This file is anything but clean!
//

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

// Let's have some debug output
#define STANDALONE
#define COLDTABLE_SETUP
#define PWPOLY_SETUP
#define DEBUG

typedef double CCTK_REAL;

#define C2P_TYPE C2P_MHD_KASTAUN

#include "../src/utils/int_pow.hh"

#include "../src/margherita.hh"
#include "../src/3D_Table/readtable.cc"

#include "../src/margherita_c2p.hh"
#include "../src/c2p/register_c2p.C"

constexpr int MAXIT =
    1; // Need to sample many inversions for performance measurement
// Don't forget to disable STANDALONE when benchmarking over many iterations..

// This is from
// https://stackoverflow.com/questions/1861294/how-to-calculate-execution-time-of-a-code-snippet-in-c
#include <sys/time.h>
typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp() {
  struct timeval now;
  gettimeofday(&now, NULL);
  return now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

using namespace Margherita_C2P_conventions;
int main() {

  std::cout << std::setprecision(15);

  // choose whether you want to test 3D or 4D
  EOS_Leptonic::use_muonic_eos = false;

  // EOS
  std::string table_name =
//      std::string("/Users/emost/Downloads/"
//                  "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5");
    //std::string("/Users/emost/Downloads/FTYSS_EOS_rho220_temp180_ye65_version2.0_20171101.h5");
    //std::string("/home/relastro-shared/EOS/scollapse/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5") ;
    //std::string("/mnt/rafast/hng/BHAC-FIL_handover/EOS/HSDD2_new.h5") ;
    std::string("/mnt/rafast/hng/BHAC-FIL_handover/EOS/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5") ;

  constexpr bool use_energy_shift = true;
  constexpr bool recompute_mu_nu = false;

  if (EOS_Leptonic::use_muonic_eos) {
  	EOS_Leptonic::readtable_scollapse(table_name.c_str(), use_energy_shift, false);
  } else {
  	EOS_Tabulated::readtable_scollapse(table_name.c_str(), use_energy_shift, false);
  }
  std::cout << "Read scollapse  Table" << std::endl;

  //EOS_Leptonic::readtable_compose(table_name.c_str());
  //std::cout << "Read compose  Table" << std::endl;

  if (EOS_Leptonic::use_muonic_eos) {
     std::string leptonic_table_name = 
       //std::string("/mnt/raarchive/hng/BHAC-FIL_handover/HSDD2_eos_ymumin_5d-4_ylemaxmin_equal_yemaxmin_Muon_electron_tables.h5") ;
       std::string("/mnt/raarchive/hng/BHAC-FIL_handover/LS220_eos_ymumin_5d-4_ylemaxmin_equal_yemaxmin_Muon_electron_tables.h5") ;
     EOS_Leptonic::readtable_leptonic(leptonic_table_name.c_str());
     std::cout << "Read leptonic  Table" << std::endl;
  }



  // Try flat space

  std::array<double, 6> metric{{1.0, 0.0, 0.0, 1.0, 0.0, 1.0}};
  std::array<double, 6> invmetric{{1.0, 0.0, 0.0, 1.0, 0.0, 1.0}};

  metric_c METRIC(metric, invmetric, 1.0);


  std::array<double, 3> Bvec = {{
      0.0, 0.0, 0.0
  }};

     std::array<double, Margherita_C2P_Registration::MAX_NUM_LIMITS> limitsL{{
         1.0e-14, // rho_atm,
         0.1,     // temp_atm,
         -1.,     // ye_atm,
         -1.,     // ymu_atm,
         1.e20,    // eps_max,
         1.e-14,  // tauenergy_min,
         1.e5,     // gamma_max,
         100      /// psi6_threshold
  }};


  // Making arrays for log Pm/P and log W parameter space
  double Pm_to_P_init_low = 1e-8;
  double Pm_to_P_init_upper = 1e10;

  double W_low = 1.0e-4 + 1.0;
  double W_upper = 1.0e4 + 1.0;

  int n_W = 200;
  int n_pm_to_p = 200;
 
  double *Pm_to_P_logtable = new double[n_pm_to_p];
  double *W_logtable = new double[n_W];

  Pm_to_P_logtable[0] = log10(Pm_to_P_init_low);
  Pm_to_P_logtable[n_pm_to_p-1] = log10(Pm_to_P_init_upper);

  W_logtable[0] = log10(W_low - 1.0);
  W_logtable[n_W-1] = log10(W_upper - 1.0);

  std::cout << W_logtable[0] << " " << W_logtable[n_W-1] << "\n";

  for (int i = 1; i < n_W-1; ++i) {
     W_logtable[i] = (W_logtable[n_W-1] - W_logtable[0])/(double)n_W + W_logtable[i-1];
  }
  for (int i = 1; i < n_pm_to_p-1; ++i) {
     Pm_to_P_logtable[i] = (Pm_to_P_logtable[n_pm_to_p-1] - Pm_to_P_logtable[0])/(double)n_pm_to_p + Pm_to_P_logtable[i-1];
  }


  std::ofstream myfile2 ("iW_of_table_merged_branch.log");
  for (int iW = 0; iW < n_W; ++iW) {
    double W;
    W = pow(10, W_logtable[iW]) + 1.0;
    myfile2 <<std::scientific  << W << "\n";
  }
  myfile2.close();

  std::ofstream myfile ("c2p_iteration_test2_fixedrhoT_3D_merged_branch.log");
  myfile << ", num of iteration , W, Pmag_to_P  : n_W, n_pm_to_p = " <<" " << n_W <<" "<< n_pm_to_p <<" "<< "\n";

  std::cout << "fk" << "\n";


  for (int iW = 0; iW < n_W; ++iW) {
     for (int ipm = 0; ipm < n_pm_to_p; ++ipm) {


        EOS_Leptonic::c2p_roofinding_call_iteration_masfunconly = 0;

        if (EOS_Leptonic::use_muonic_eos) {

           double ZVECX_local, Pm_to_P_init, W_init, b2;
           Pm_to_P_init = :: pow(10,Pm_to_P_logtable[ipm]);
           W_init = pow(10,W_logtable[iW]) + 1.0; 
           // In FIL, ZVECX is Wv^x

           ZVECX_local = W_init * sqrt(1.0 - (1.0/ ::int_pow<2>(W_init)  )   );



           // Set primitives
           std::array<double, C2P_MHD<EOS_Leptonic>::numprims> PRIMS{{
               1e11*Margherita_constants::RHOGF, // RHOB
               0.0,       // PRESSURE
               ZVECX_local,     // ZVECX
               0.0,     // ZVECY
               0.0,     // ZVECZ
               1.e-0,      // EPS
               0.0,       // CS2
               //0.0,       // ENTROPY
               5.0,       // TEMP
               0.1,        // Ye
               0.05        // Ymu
           }};
      
      
      
           EOS_Leptonic::atm_muonic_beta_eq = true;
           Margherita_register_limits<EOS_Leptonic>(limitsL);
           std::array<double, C2P_MHD<EOS_Leptonic>::numprims> PRIMSN;
           EOS_Leptonic::error_type error; // TODO How do we handle this?
           PRIMS[EPS] = EOS_Leptonic::eps__temp_rho_yle_ymu(PRIMS[TEMP],PRIMS[RHOB],PRIMS[YE],PRIMS[YMU],error);
           double h,s;
           auto epsrange =
               EOS_Leptonic::eps_range__rho_yle_ymu(PRIMS[RHOB], PRIMS[YE], PRIMS[YMU], error);
           PRIMS[EPS] = std::max(epsrange[0], std::min(epsrange[1], PRIMS[EPS]));
           // 4. Compute a by making a pressure call
           // Update all vars here, cs2, temp, etc..
           PRIMS[PRESSURE] = EOS_Leptonic::press_h_csnd2_temp_entropy__eps_rho_yle_ymu(
               h, PRIMS[CS2], PRIMS[TEMP], s, // out all but temp (inout)
               PRIMS[EPS], PRIMS[RHOB], PRIMS[YE], PRIMS[YMU], error); // in


           b2 = Pm_to_P_init * 2.0 * PRIMS[PRESSURE];
           Bvec[0] = 0.0;
           Bvec[1] = sqrt(b2) * W_init;
           Bvec[2] = 0.0;

           // Declare conservatives
           std::array<double, C2P_MHD<EOS_Leptonic>::numcons> CONS;
      
           // Need to compute Lorentz factor and zvec_low
           std::array<double, 3> zvec_lowL = METRIC.lower_index<ZVECX, NUM_PRIMS>(PRIMS);
           const double z2L = PRIMS[ZVECX] * zvec_lowL[0] + PRIMS[ZVECY] * zvec_lowL[1] +
                              PRIMS[ZVECZ] * zvec_lowL[2];
           const double WL = sqrt(1. + z2L);
      
           // Compute conservatives
           CONS = C2P_MHD<EOS_Leptonic>::compute_conservatives(PRIMS, METRIC, Bvec);
           PRIMSN[YMU] = PRIMS[YMU];
           PRIMSN[YE] = PRIMS[YE];
           PRIMSN[TEMP] = PRIMS[TEMP];
      
           ////////    Invert conservatives
             double secs = 0.;
      
           // Start timing
           timestamp_t t0 = get_timestamp();
           // Invert conservatives for MAXIT times to benchmark the inversion
           for (int i = 0; i < MAXIT; ++i) {
      
      
             // Compute conservatives
             CONS = C2P_MHD<EOS_Leptonic>::compute_conservatives(PRIMS, METRIC, Bvec);
             PRIMSN[YMU] = PRIMS[YMU];
             PRIMSN[YE] = PRIMS[YE];
             PRIMSN[TEMP] = PRIMS[TEMP];
             error_bits_t error_bits;
                Margherita_Con2Prim<C2P_TYPE, EOS_Leptonic>(CONS, PRIMSN, METRIC,
            						 error_bits);
           }
           
           timestamp_t t1 = get_timestamp();
           // Add to global timer
           secs += (t1 - t0) / 1000000.0L;
           secs /= MAXIT;
           
           const std::array<std::string, NUM_PRIMS> names{
             {"RHOB", "PRESSURE", "ZVECX", "ZVECY", "ZVECZ", "EPS", "CS2", 
             //{"RHOB", "PRESSURE", "ZVECX", "ZVECY", "ZVECZ", "EPS", "CS2", "ENTROPY",
               "TEMP", "YE", "YMU"}};
           
           std::cout << std::endl;
           std::cout << "Specific enthalpy * lorentz: " << h * WL << std::endl;
           std::cout << std::endl;
           // Output old and new primitives
           double *relative_diff_prim = new double[NUM_PRIMS];
           for (int i = 0; i < NUM_PRIMS; ++i) {
             std::cout << names[i] << ": " << PRIMS[i] << " , " << PRIMSN[i]
                       << std::endl;
             relative_diff_prim[i] = std::abs( (PRIMS[i] - PRIMSN[i])/PRIMS[i] );
           }
           double rela_diff_average;

           // one of them is B^i: rela diff = 0
           rela_diff_average = 1./5. * ( relative_diff_prim[RHOB] + 
                                         relative_diff_prim[ZVECX] +
                                         relative_diff_prim[YE] +
                                         0.0 + 
                                         relative_diff_prim[EPS]
                                       );
           // this identified to be unsuccessful, tol_recov =  5e-9 * 1e1
           if ( rela_diff_average > 5e-8) {
              EOS_Leptonic::c2p_roofinding_call_iteration_masfunconly = 1000000;
           }
           // Output average timing
           std::cout << "\n Average time[s]: " << secs << std::endl;
      
           myfile << EOS_Leptonic::c2p_roofinding_call_iteration_masfunconly <<" "<< W_init <<" "<< Pm_to_P_init <<" "<< "\n";
      
      
        } else {

           double ZVECX_local, Pm_to_P_init, W_init, b2;
           Pm_to_P_init = :: pow(10,Pm_to_P_logtable[ipm]);
           W_init = pow(10,W_logtable[iW]) + 1.0;
           // In FIL, ZVECX is Wv^x

           ZVECX_local = W_init * sqrt(1.0 - (1.0/ ::int_pow<2>(W_init)  )   );

           // Set primitives
           std::array<double, C2P_MHD<EOS_Tabulated>::numprims> PRIMS{{
               1e11*Margherita_constants::RHOGF, // RHOB
               0.0,       // PRESSURE
               ZVECX_local,     // ZVECX
               0.0,     // ZVECY
               0.0,     // ZVECZ
               1.e-0,      // EPS
               0.0,       // CS2
               //0.0,       // ENTROPY
               5.0,       // TEMP
               0.1,        // Ye
               0.0        // Ymu
           }};

           EOS_Tabulated::atm_beta_eq = true;
           Margherita_register_limits<EOS_Tabulated>(limitsL);
           std::array<double, C2P_MHD<EOS_Tabulated>::numprims> PRIMSN;
           EOS_Tabulated::error_type error; // TODO How do we handle this?
           double dummy;
           PRIMS[EPS] = EOS_Tabulated::eps__temp_rho_yle_ymu(PRIMS[TEMP],PRIMS[RHOB],PRIMS[YE],PRIMS[YMU],error);
           double h,s;
           auto epsrange =
               EOS_Tabulated::eps_range__rho_yle_ymu(PRIMS[RHOB], PRIMS[YE], PRIMS[YMU], error);
           PRIMS[EPS] = std::max(epsrange[0], std::min(epsrange[1], PRIMS[EPS]));
           // 4. Compute a by making a pressure call
           // Update all vars here, cs2, temp, etc..
           PRIMS[PRESSURE] = EOS_Tabulated::press_h_csnd2_temp_entropy__eps_rho_yle_ymu(
               h, PRIMS[CS2], PRIMS[TEMP], s,  // out all but temp (inout)
               PRIMS[EPS], PRIMS[RHOB], PRIMS[YE], PRIMS[YMU], error); // in

           b2 = Pm_to_P_init * 2.0 * PRIMS[PRESSURE];
           Bvec[0] = 0.0;
           Bvec[1] = sqrt(b2) * W_init;
           Bvec[2] = 0.0;

           // Declare conservatives
           std::array<double, C2P_MHD<EOS_Tabulated>::numcons> CONS;
      
           // Need to compute Lorentz factor and zvec_low
           std::array<double, 3> zvec_lowL = METRIC.lower_index<ZVECX, NUM_PRIMS>(PRIMS);
           const double z2L = PRIMS[ZVECX] * zvec_lowL[0] + PRIMS[ZVECY] * zvec_lowL[1] +
                              PRIMS[ZVECZ] * zvec_lowL[2];
           const double WL = sqrt(1. + z2L);
      
           // Compute conservatives
           CONS = C2P_MHD<EOS_Tabulated>::compute_conservatives(PRIMS, METRIC, Bvec);
           PRIMSN[YMU] = 0.0;
           PRIMSN[YE] = PRIMS[YE];
           PRIMSN[TEMP] = PRIMS[TEMP];
      
           //- ------------------not added to muonic eos  yet
      
           // Fun test with hydro conservatives
           /*  auto CONS2 = C2P_Hydro<EOS_Polytropic>::compute_conservatives(PRIMS,
               METRIC);
               for(int j=0; j<C2P_Hydro<EOS_Polytropic>::numcons; ++j) {
               std::cout << CONS[j] << "  ,  " << CONS2[j] <<std::endl;
               CONS[j]=CONS2[j];
               }
           */
           ////////    Invert conservatives
             double secs = 0.;
      
           // Start timing
           timestamp_t t0 = get_timestamp();
           // Invert conservatives for MAXIT times to benchmark the inversion
           for (int i = 0; i < MAXIT; ++i) {
      
      
                    // Compute conservatives
                CONS = C2P_MHD<EOS_Tabulated>::compute_conservatives(PRIMS, METRIC, Bvec);
                PRIMSN[YMU] = 0.0;
             // Fun test with hydro conservatives
             /*  auto CONS2 = C2P_Hydro<EOS_Polytropic>::compute_conservatives(PRIMS,
               METRIC);
               for(int j=0; j<C2P_Hydro<EOS_Polytropic>::numcons; ++j) {
                 std::cout << CONS[j] << "  ,  " << CONS2[j] <<std::endl;
                 CONS[j]=CONS2[j];
               }
           */
             PRIMSN[YE] = PRIMS[YE];
             PRIMSN[TEMP] = PRIMS[TEMP];
             error_bits_t error_bits;
                Margherita_Con2Prim<C2P_TYPE, EOS_Tabulated>(CONS, PRIMSN, METRIC,
            						 error_bits);
           }
           
           timestamp_t t1 = get_timestamp();
           // Add to global timer
           secs += (t1 - t0) / 1000000.0L;
           secs /= MAXIT;
           
           const std::array<std::string, NUM_PRIMS> names{
             {"RHOB", "PRESSURE", "ZVECX", "ZVECY", "ZVECZ", "EPS", "CS2", 
             //{"RHOB", "PRESSURE", "ZVECX", "ZVECY", "ZVECZ", "EPS", "CS2", "ENTROPY",
               "TEMP", "YE", "YMU"}};
           
           std::cout << std::endl;
           std::cout << "Specific enthalpy * lorentz: " << h * WL << std::endl;
           std::cout << std::endl;
           // Output old and new primitives
           double *relative_diff_prim = new double[NUM_PRIMS];
           for (int i = 0; i < NUM_PRIMS; ++i) {
             std::cout << names[i] << ": " << PRIMS[i] << " , " << PRIMSN[i]
                       << std::endl;
             relative_diff_prim[i] = std::abs( (PRIMS[i] - PRIMSN[i])/PRIMS[i] );
           }
           double rela_diff_average;

           // one of them is B^i: rela diff = 0
           rela_diff_average = 1./5. * ( relative_diff_prim[RHOB] +
                                         relative_diff_prim[ZVECX] +
                                         relative_diff_prim[YE] +
                                         0.0 +
                                         relative_diff_prim[EPS] 
                                       );
           // this identified to be unsuccessful, tol_recov =  5e-9 * 1e1
           if ( rela_diff_average > 5e-8) {
              EOS_Leptonic::c2p_roofinding_call_iteration_masfunconly = 1000000;
           }

           // Output average timing
           std::cout << "\n Average time[s]: " << secs << std::endl;

           myfile << EOS_Leptonic::c2p_roofinding_call_iteration_masfunconly <<" "<< W_init <<" "<< Pm_to_P_init <<" "<< "\n";
        } // endif of use_muonic_eos
     }
  } // for loop of iW

  myfile.close();
  
  return 0;
}
