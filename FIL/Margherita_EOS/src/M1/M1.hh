#include "../margherita.hh"
#include "../Weakhub/Weakhub.hh"
#include "fermi.hh"
#include "fugacities.hh"
#include "helpers.hh"
#include "constexpr_utils.hh"
#include <cassert> 

#ifdef DEBUG
#include <string>
#include <iostream>
#endif

#ifndef __H_M1
#define __H_M1

#ifdef DEBUG
enum m1_species_t { NUE=0, NUE_BAR, NUMU, NUMU_BAR, NUX, NUMSPECIES} ;
enum m1_fluid_index_t { RHO=0, YE, YMU, TEMP, REQ_FLUID_VARS} ;
#endif 

namespace Margherita_M1_EAS { 
using namespace m1_helpers; 
using namespace M1_Constants;
using namespace Margherita_helpers;
using namespace Margherita_constants;

#ifdef DEBUG
const std::array<std::string, 5> names { 
    "NUE",
    "NUE_BAR",
    "NUMU",
    "NUMU_BAR",
    "NUX"
} ; 
#endif

template <typename T> 
class EAS
{
private:
    // template shorthands 
    template <int N>
    using FD = Fermi_Dirac<N> ; 
    template <int N, int M> 
    using FDR = Fermi_Dirac_Ratio<N,M> ; 
    using FG  = Fugacities<T> ;

    using arr_t = std::array<T, NUMSPECIES> ;

    // number of DoFs 
//Harry: becareful , 3 spec: g_nu = {1, 1, 0, 0, 4}
		 //  5 spec: g_nu = {1, 1, 1, 1, 2}
    std::array<T,5> const g_nu_3spec {1, 1, 0, 0, 4} ;
    std::array<T,5> const g_nu_5spec {1, 1, 1, 1, 2} ;
    // which interactions are active
    std::array<bool, NINTERACTIONS> which_interactions   ;
    // have the rates been computed yet? 
    bool computed = false ; 

    FG F  ;

    // Rates

    // Emission: MeV / s / cm^3
    arr_t Q {0., 0., 0., 0., 0.} ;
    arr_t Q_diff_inv {0., 0., 0., 0., 0.} ;
    arr_t Q_eff {0., 0., 0., 0., 0.} ;

    // Number emission: number / s / cm^3
    arr_t R {0., 0., 0., 0., 0.} ;
    arr_t R_diff_inv {0., 0., 0., 0., 0.} ;
    arr_t R_eff {0., 0., 0., 0., 0.} ;


    arr_t tau_n {0., 0., 0., 0., 0.} ;
    arr_t tau_e {0., 0., 0., 0., 0.} ;

    // Mean energy (squared): MeV^2
    arr_t eps2_leak {0., 0., 0., 0., 0.} ;
    arr_t eps2_fluid {0., 0., 0., 0., 0.} ;
    arr_t eps {0., 0., 0., 0., 0.} ; // MeV 

    // Opacities: 1/cm
    arr_t kappa_s {0., 0., 0., 0., 0.} ;
    arr_t kappa_a {0., 0., 0., 0., 0.} ;
    arr_t kappa_n {0., 0., 0., 0., 0.} ;

    arr_t Y_nu {0., 0., 0., 0., 0.} ;
    arr_t Z_nu {0., 0., 0., 0., 0.} ;

public:
    // Ctors
    // ========================================================
    EAS(FG &_F, 
        std::array<bool,NINTERACTIONS> const& wi) : 
        F(_F), which_interactions(wi) 
        {  
            tau_n = _F.tau_n; 
        } ;
    
    EAS(const T &rho, 
        const T &ye, const T &ymu, 
        const T &temp, 
        std::array<bool,NINTERACTIONS> const& wi)
    : F(FG(rho, temp, ye, ymu)),  which_interactions(wi) {};
    // ========================================================
    // Dtor
    virtual ~EAS() = default;
    // ========================================================
    // Setter for fugacities
    void set_F(FG const& __F){
        F = __F ; 
        computed=false ;
    }
    // ========================================================
    // Access operators 
    // ========================================================
    template<int i>
    inline decltype(auto) get_tau() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        return tau_n[i] ;
    }
    // Emissivities 
    template< int i >
    inline decltype(auto) Q_cgs() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        assert(computed) ;
        return Q[i] * mev_to_erg ; // erg / s / cm^3 
    }

    template< int i >
    inline decltype(auto) Q_cactus() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        assert(computed) ;
        return Q[i] * mev_to_erg * RHOGF * EPSGF / TIMEGF ; 
    }
    // Number Emissivities 
    template< int i >
    inline decltype(auto) R_cactus() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        assert(computed) ;
        return R[i] * mnuc_cgs * RHOGF / TIMEGF ; 
    }
    // Opacities
    template< int i >
    inline decltype(auto) kappa_s_cgs() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        assert(computed) ;
        return kappa_s[i] ; 
    }
    template< int i >
    inline decltype(auto) kappa_s_cactus() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        assert(computed) ;
        return kappa_s[i] / LENGTHGF ; 
    }

    template< int i >
    inline decltype(auto) kappa_a_cgs() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        assert(computed) ;
        return kappa_a[i] ; 
    }
    template< int i >
    inline decltype(auto) kappa_a_cactus() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        assert(computed) ;
        return kappa_a[i] / LENGTHGF ; 
    }

    template< int i >
    inline decltype(auto) kappa_n_cgs() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        assert(computed) ;
        return kappa_n[i] ; 
    }
    template< int i >
    inline decltype(auto) kappa_n_cactus( ) const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        assert(computed) ;
        return kappa_n[i] / LENGTHGF; 
    }
    // BlackBody functions 
    template< int i, int n> 
    inline __attribute__((always_inline)) 
    decltype(auto) black_body_mev() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        static_assert( (n==0) || (n==1) ) ;
        if (EOS_Leptonic::use_muonic_eos) {
        	auto const bb = g_nu_5spec[i] * 4.*pi*clite/::int_pow<3>(hc_mevcm) * ::int_pow<3+n>(F.temp) * FD<2+n>::get(F.eta[i]) ; 
        	return std::max(bb, 1.e-30)  ;
        } else {
        	auto const bb = g_nu_3spec[i] * 4.*pi*clite/::int_pow<3>(hc_mevcm) * ::int_pow<3+n>(F.temp) * FD<2+n>::get(F.eta[i]) ; 
        	return std::max(bb, 1.e-30)  ;
        }
    } 

    template< int i, int n> 
    inline __attribute__((always_inline)) 
    decltype(auto) black_body_cactus() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        static_assert( (n==0) || (n==1) ) ;
        auto const bb = black_body_mev<i,n>() * mev_to_erg * EPSGF * RHOGF / TIMEGF * LENGTHGF ;
        return std::max(bb, 1.e-30)  ;
    } 

    // Mean Energies
  template<int i, int n> 
    inline T eps2_fluid_cgs() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        return eps2_fluid_mev<i,n>() * ::int_pow<2>(mev_to_erg);
    }
  template<int i, int n> 
    inline T eps2_fluid_mev() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        return FDR<4+n,2+n>::get(F.eta[i] ) * ::int_pow<2>(F.temp) ;;
    }
    template<int i> 
    inline T eps2_leak_cgs() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        return eps2_leak[i] * ::int_pow<2>(mev_to_erg);
    }
    template<int i> 
    inline T eps_cgs() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        return eps[i] * mev_to_erg;
    }
    template<int i> 
    inline T eps_mev() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        return eps[i] ;
    }
    template<int i> 
    inline T eps_cactus() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        return eps[i] * mev_to_erg 
               * RHOGF * EPSGF * ::int_pow<3>(LENGTHGF);
    }
    
    // Compute the temperature of neutrinos according to FD distribution
    static T neutrino_temperature(T const& eps, T const& eta){
       return FDR<2,3>::get(eta) * eps ;
    }
    // Compute very approximate optical depth based on density
    static inline double compute_analytic_opacity(const double rho) {
        // This computes an empirical value of the opacity for cold neutron stars
        // as given in Deaton et al. 2013. (10)
        const auto rho_cgs = rho * INVRHOGF;
        const auto log_tau_approx = 0.96 * (log(rho_cgs) / log(10.) - 11.7);
        return exp(log(10.) * log_tau_approx);
    }

    template<int i, int n>
    inline __attribute__((always_inline))
    T get_neutrino_density() {
        auto const rho_lim = 1e11 ; // g / cm^3 
        auto Ynu = 4. * pi  / ::int_pow<3>(hc_mevcm) * 
                ::int_pow<3+n>(F.temp) * FD<2+n>::get(F.eta[i]) * exp(-rho_lim/F.rho) ;
        if ( n==1 ) { // ! This is in mev / cm^3 
            Ynu *= mev_to_erg * EPSGF * RHOGF ;
        } else { //! This is adimensional -- rho is CGS inside of Fugacities 
            Ynu *= mnuc_cgs / F.rho ; 
        }
        return Ynu ;  
    }
    // ========================================================
    // ========================================================
    // ========================================================
    // Computation of the rates
    // ========================================================
    // ========================================================
    // ========================================================
    
    // ! The following functions add the relevant 
    // ! rates to the input 
    // ========================================================
    // Opacities 
    // ========================================================
    inline __attribute__((always_inline))
    void add_charged_current_absorption_opacity(arr_t& ka,
                                                 arr_t& kn ) ; 
    
    inline __attribute__((always_inline))
    void add_Kirchoff_absorption_opacity( arr_t& ka,
                                          arr_t& kn,
                                          const arr_t& QQ,
                                          const arr_t& RR, double& rho_local) ;

    inline __attribute__((always_inline))
    void add_scattering_opacity(arr_t& ks) ; 

    // ========================================================
    // Emissivities
    // ========================================================
    inline __attribute__((always_inline))
    void add_charged_current_emission(arr_t& QQ,
                                       arr_t& RR ) ;

    inline __attribute__((always_inline))
    void add_plasmon_decay_emission(arr_t& QQ,
                                     arr_t& RR ) ;

    inline __attribute__((always_inline))
    void add_pair_process_emission(arr_t& QQ,
                                    arr_t& RR ) ;

    inline __attribute__((always_inline))
    void add_brems_emission(arr_t& QQ,
                             arr_t& RR ) ;


    // ! NOTE: this OVERWRITES the rates of NUE, NUE_BAR
    inline __attribute__((always_inline))
    void calc_Kirchoff_emission(arr_t& QQ,
                                arr_t& RR,
                                const arr_t& ka,
                                const arr_t& kn, bool need_nux_emission ) ;


    // Compute leakage effective rates and 
    // mean energies 
    inline __attribute__((always_inline))
    void calc_diffusion_rates() {
        constexpr auto factor = 3. / a_diff * clite / ::int_pow<3>(hc_mevcm) * normfact ;
        // Rosswog+ 2002, Sekiguchi 2002
        for(int ii=NUE; ii<=NUX; ++ii){
            const auto E_sq = ::int_pow<2>(F.temp) * FDR<4,2>::get(F.eta[ii]) ; // MeV

            if (EOS_Leptonic::use_muonic_eos) {
                 R_diff_inv[ii] =
                     ::int_pow<2>( tau_n[ii] ) / ( factor*g_nu_5spec[ii] * kappa_n[ii] * F.temp * FD<0>::get(F.eta[ii]) * E_sq ) ;

                 const auto E_sq_e = ::int_pow<2>(F.temp) * FDR<5,3>::get(F.eta[ii]) ; // MeV
                 Q_diff_inv[ii] =
                     ::int_pow<2>(tau_n[ii]) / (factor*g_nu_5spec[ii]*kappa_a[ii]*::int_pow<2>(F.temp)*FD<1>::get(F.eta[ii]) * E_sq_e ) ;
            } else {
           	 R_diff_inv[ii] = 
           	     ::int_pow<2>( tau_n[ii] ) / ( factor*g_nu_3spec[ii] * kappa_n[ii] * F.temp * FD<0>::get(F.eta[ii]) * E_sq ) ;

           	 const auto E_sq_e = ::int_pow<2>(F.temp) * FDR<5,3>::get(F.eta[ii]) ; // MeV
           	 Q_diff_inv[ii] = 
           	     ::int_pow<2>(tau_n[ii]) / (factor*g_nu_3spec[ii]*kappa_a[ii]*::int_pow<2>(F.temp)*FD<1>::get(F.eta[ii]) * E_sq_e ) ;
	    }
        } // for loop
     }

    inline __attribute__((always_inline))
    void calc_effective_rates() {
        calc_diffusion_rates() ;
        for(int i=NUE; i<=NUX; ++i){
            if ( std::isfinite(R[i] * R_diff_inv[i]) &&
                 std::isfinite( 10. * R_diff_inv[i])   ) 
                 {
                     R_eff[i] = R[i] / (1. + normfact * R[i] * R_diff_inv[i]) ;
                     if( !std::isnormal( R_diff_inv[i] ) )
                        R_eff[i] = R[i] ; 
                 }  else {
                     R_eff[i] = 0. ; 
                 }

            if ( std::isfinite(Q[i] * Q_diff_inv[i]) &&
                 std::isfinite( 10. * Q_diff_inv[i])   ) 
                 {
                     Q_eff[i] = Q[i] / (1. + normfact * Q[i] * Q_diff_inv[i]) ;
                     if( !std::isnormal( Q_diff_inv[i] ) )
                        Q_eff[i] = Q[i] ; 
                 }  else {
                     Q_eff[i] = 0. ; 
                 }
        }
     }

    inline __attribute__((always_inline))
    void calc_energies() {
        calc_effective_rates() ;
        for(int i=NUE; i<=NUX; ++i){
            if( R_eff[i] > 0. ) {
                eps2_leak[i] = ::int_pow<2>( Q_eff[i] / R_eff[i] ) ;
            } else {
                eps2_leak[i] = 0. ;
            }
            eps2_fluid[i] = FDR<5,3>::get(F.eta[i] ) * ::int_pow<2>(F.temp) ;
            eps[i] = FDR<5,4>::get(F.eta[i]) * F.temp ; 
        }

    }

    // ========================================================
    // Weakhub routines
    // ========================================================
    inline __attribute__((always_inline))
    void kappa_ast_kappa_s__temp_rho_yle_ymu_M1(
        arr_t& kappa_a_st_en, arr_t& kappa_a_st_num, arr_t& kappa_s, int num_spec, double &temp,
        double &rho, double &yle, double &ymu);


    // ========================================================
    // Computation of the rates
    // ========================================================
    inline void calc_eas(bool recompute_emission_from_absorption=true) {
        Q = {0., 0., 0., 0., 0.}; R={0.,0.,0.,0.,0.} ;
        kappa_s = {0., 0., 0., 0., 0.}; kappa_a = {0.,0.,0.,0.,0.}; kappa_n = {0., 0., 0., 0., 0.}; 


	// Harry: note that we should ignore pair processes contribution to nue nue_bar, 
	          // due to the ignored neutrino blocking factor of nue and nue_bar in these processes
                  // Weakhub did this automically

        if (M1_Weakhub::use_Weakhub_eas) {
        // Inside Weakhub, everything related to micrphysics is already decided in the table
        // declare with F.  is it okay? FIXME
           double ye_pt = F.ye;
           double ymu_pt = F.ymu;
           //double ymu_pt = F.ymu;
           double rho_pt = F.rho;
           double temp_pt = F.temp;

           // 1. get kappa_a, kappa_n and kappa_s for all species
           if (EOS_Leptonic::use_muonic_eos) {
             kappa_ast_kappa_s__temp_rho_yle_ymu_M1(
                   kappa_a, kappa_n, kappa_s, 5, temp_pt,
                   rho_pt, ye_pt, ymu_pt);
           } else {
             double dummy = 0.0;
             kappa_ast_kappa_s__temp_rho_yle_ymu_M1(
                   kappa_a, kappa_n, kappa_s, 3, temp_pt,
                   rho_pt, ye_pt, dummy);
           }

           // Weakhub only provides opacity for all species --> before adding numu numubar nux from pair process 
           // apply the Kirchoff first, later PP could use different fugacities
           calc_Kirchoff_emission(Q, R, kappa_a, kappa_n, true) ;

           // 2.  Plasma process in Weakhub does not have approximated form yet. so get Q,R for this interaction
           arr_t Q_tmp {0., 0., 0., 0., 0. } ;
           arr_t R_tmp {0., 0., 0., 0., 0. } ;
           if ( which_interactions[PLASMON_DECAY] ) {
               add_plasmon_decay_emission(Q_tmp, R_tmp) ;
           }
           if ( which_interactions[BREMS] )
               add_brems_emission(Q_tmp,R_tmp) ;
           if ( which_interactions[PAIR])
               add_pair_process_emission(Q_tmp, R_tmp) ;

           arr_t kappa_a_tmp {0., 0., 0., 0., 0. } ;
           arr_t kappa_n_tmp {0., 0., 0., 0., 0. } ;

           // 3. get kappa_a/n for plasma process only for numu numu_bar nux
           add_Kirchoff_absorption_opacity(kappa_a_tmp, kappa_n_tmp, Q_tmp, R_tmp, rho_pt) ;

           // 4. add up all kappa_a/n
           for( int i=NUE; i<=NUX; i++ ) {
               kappa_a[i] += kappa_a_tmp[i];
               kappa_n[i] += kappa_n_tmp[i];
               Q[i] += Q_tmp[i];
               R[i] += R_tmp[i];
           }


           // FIXME: I should not put here, since will suppress pair processes as well, but still, ee PP and NNbrem PP also included inside....
           if (EOS_Leptonic::use_muonic_eos) {
             rho_pt = F.rho;
             temp_pt = F.temp;
             //// Suppression factor from Bollig thesis eq.9.31, and Bollig told only need density suppression for beta processes
             //// but Bollig said density suppression is enough, and 2.5 mev suppression temp is quite high??
             //kappa_a[NUMU]     = kappa_a[NUMU]     / (1.0 + ::int_pow<5>(1.0e11 / (rho_pt) ) ) / (1.0 + ::int_pow<6>(2.5 / (temp_pt) ) );
             //kappa_a[NUMU_BAR] = kappa_a[NUMU_BAR] / (1.0 + ::int_pow<5>(1.0e11 / (rho_pt) ) ) / (1.0 + ::int_pow<6>(2.5 / (temp_pt) ) );
             //kappa_n[NUMU]     = kappa_n[NUMU]     / (1.0 + ::int_pow<5>(1.0e11 / (rho_pt) ) ) / (1.0 + ::int_pow<6>(2.5 / (temp_pt) ) );
             //kappa_n[NUMU_BAR] = kappa_n[NUMU_BAR] / (1.0 + ::int_pow<5>(1.0e11 / (rho_pt) ) ) / (1.0 + ::int_pow<6>(2.5 / (temp_pt) ) );

             //Q[NUMU]     = Q[NUMU] / (1.0 + ::int_pow<5>(1.0e11 / (rho_pt) ) )     / (1.0 + ::int_pow<6>(2.5 / (temp_pt) ) );
             //Q[NUMU_BAR] = Q[NUMU_BAR] / (1.0 + ::int_pow<5>(1.0e11 / (rho_pt) ) ) / (1.0 + ::int_pow<6>(2.5 / (temp_pt) ) );
             //R[NUMU]     = R[NUMU] / (1.0 + ::int_pow<5>(1.0e11 / (rho_pt) ) )     / (1.0 + ::int_pow<6>(2.5 / (temp_pt) ) );
             //R[NUMU_BAR] = R[NUMU_BAR] / (1.0 + ::int_pow<5>(1.0e11 / (rho_pt) ) ) / (1.0 + ::int_pow<6>(2.5 / (temp_pt) ) );

 	    if ( rho_pt < 1.0e10 || temp_pt < 2.5) {
	    	kappa_a[NUMU] = 0.0;
	    	kappa_a[NUMU_BAR] = 0.0;
	    	kappa_n[NUMU] = 0.0;
	    	kappa_n[NUMU_BAR] = 0.0;

	    	Q[NUMU] = 0.0;
	    	Q[NUMU_BAR] = 0.0;
	    	R[NUMU] = 0.0;
	    	R[NUMU_BAR] = 0.0;
            }
           }

           // Ensure opacity >= 0.0 and non NaN
           for (int i = NUE; i <= NUX; i++) {
             kappa_a[i] = std::max(kappa_a[i], 1.e-60);
             kappa_n[i] = std::max(kappa_n[i], 1.e-60);
             kappa_s[i] = std::max(kappa_s[i], 1.e-60);
	     Q[i] = std::max(Q[i], 1.e-60);
	     R[i] = std::max(R[i], 1.e-60);
             if (kappa_a[i] != kappa_a[i]) {kappa_a[i] = 1.e-60;}
             if (kappa_n[i] != kappa_n[i]) {kappa_n[i] = 1.e-60;}
             if (kappa_s[i] != kappa_s[i]) {kappa_s[i] = 1.e-60;}
           }

//harry:no need this, 
           //calc_energies() ;

           computed = true ;
     
        } else { //other m1_eas methods (analytical)
           double rho_pt = F.rho;

           // free absorption rates 
           if ( which_interactions[BETA] )
               add_charged_current_absorption_opacity(kappa_a, kappa_n) ;
           // scattering 
           add_scattering_opacity(kappa_s) ;

           // Emission rates: 
           arr_t Q_beta{0., 0., 0., 0., 0.}, R_beta{0., 0., 0., 0., 0.} ;
           if( which_interactions[BETA] )
               add_charged_current_emission(Q_beta, R_beta) ;
           
           if ( which_interactions[PLASMON_DECAY] ) 
               add_plasmon_decay_emission(Q, R) ;
           if ( which_interactions[BREMS] )
               add_brems_emission(Q,R) ;
           if ( which_interactions[PAIR])
               add_pair_process_emission(Q, R) ;
           
           arr_t kappa_a_tmp {0., 0., 0., 0., 0.} ;
           arr_t kappa_n_tmp {0., 0., 0., 0., 0.} ;

           // get kappa_a/n of nux by putting Q_nux R_nux
           add_Kirchoff_absorption_opacity(kappa_a_tmp, kappa_n_tmp, Q, R, rho_pt) ;

           // use kappa_a (only nue nue_bar) to get Q, R nue nue_bar
           if ( recompute_emission_from_absorption )
               calc_Kirchoff_emission(Q, R, kappa_a, kappa_n, false) ;

           for( int i=NUE; i<=NUX; i++ ) {
               kappa_a[i] += kappa_a_tmp[i];
               kappa_n[i] += kappa_n_tmp[i];
           }

           calc_energies() ;


           computed = true ;         
        } // endif of m1_which_eas


        #ifdef DEBUG
        std::cout << "Total" << std::endl;
        std::cout << "kappa_a : ";

        std::cout << kappa_a[0] << " , ";
        std::cout << kappa_a[1] << " , ";
        std::cout << kappa_a[2] << " , ";
        std::cout << kappa_a[3] << " , ";
        std::cout << kappa_a[4] << " , ";
        std::cout << std::endl;

        std::cout << "Q : ";

        std::cout << Q[0] << " , ";
        std::cout << Q[1] << " , ";
        std::cout << Q[2] << " , ";
        std::cout << Q[3] << " , ";
        std::cout << Q[4] << " , ";
        std::cout << std::endl;

        std::cout << "kappa_n : ";

        std::cout << kappa_n[0] << " , ";
        std::cout << kappa_n[1] << " , ";
        std::cout << kappa_n[2] << " , ";
        std::cout << kappa_n[3] << " , ";
        std::cout << kappa_n[4] << " , ";
        std::cout << std::endl;

        std::cout << "R : ";
        std::cout << R[0] << " , ";
        std::cout << R[1] << " , ";
        std::cout << R[2] << " , ";
        std::cout << R[3] << " , ";
        std::cout << R[4] << " , ";

        std::cout << std::endl;

        std::cout << "Total [code units]" << std::endl;
        std::cout << "kappa_a : ";

        std::cout << kappa_a[0] / LENGTHGF<< " , ";
        std::cout << kappa_a[1] / LENGTHGF<< " , ";
        std::cout << kappa_a[2] / LENGTHGF<< " , ";
        std::cout << kappa_a[3] / LENGTHGF<< " , ";
        std::cout << kappa_a[4] / LENGTHGF<< " , ";
        std::cout << std::endl;

        std::cout << "Q : ";

        std::cout << Q[0] * mev_to_erg * EPSGF * RHOGF / TIMEGF << " , ";
        std::cout << Q[1] * mev_to_erg * EPSGF * RHOGF / TIMEGF << " , ";
        std::cout << Q[2] * mev_to_erg * EPSGF * RHOGF / TIMEGF << " , ";
        std::cout << Q[3] * mev_to_erg * EPSGF * RHOGF / TIMEGF << " , ";
        std::cout << Q[4] * mev_to_erg * EPSGF * RHOGF / TIMEGF << " , ";
        std::cout << std::endl;
        
        std::cout << "kappa_n : ";

        std::cout << kappa_n[0] / LENGTHGF<< " , ";
        std::cout << kappa_n[1] / LENGTHGF<< " , ";
        std::cout << kappa_n[2] / LENGTHGF<< " , ";
        std::cout << kappa_n[3] / LENGTHGF<< " , ";
        std::cout << kappa_n[4] / LENGTHGF<< " , ";
        std::cout << std::endl;

        std::cout << "R : ";
        std::cout << R[0] * EPSGF * RHOGF / TIMEGF << " , ";
        std::cout << R[1] * EPSGF * RHOGF / TIMEGF << " , ";
        std::cout << R[2] * EPSGF * RHOGF / TIMEGF << " , ";
        std::cout << R[3] * EPSGF * RHOGF / TIMEGF << " , ";
        std::cout << R[4] * EPSGF * RHOGF / TIMEGF << " , ";

        std::cout << std::endl;
        #endif 
    }

    #ifdef DEBUG 
    #endif
};

// using analytical rates, no accounting for numu and numubar in charged_current
template<typename T>
inline __attribute__((always_inline))
void EAS<T>::add_charged_current_absorption_opacity(std::array<T, NUMSPECIES>& ka, 
                                                    std::array<T, NUMSPECIES>& kn ) { 

    // ! Note that heavy lepton neutrinos do not undergo these processes 
    // ! ( We don't have data about muons and tau in the EOS tables )                                             
    arr_t zeta_tmp {0., 0., 0., 0., 0.} ;
    
    constexpr auto const abs_const = 0.25 * ( 1. + 3. * ::int_pow<2>(alpha) ) * sigma_0 ; // ? alpha or g_A ? 
    
    // nue + n --> p + e- 
    auto const block_factor_n = 
                         1. + exp( F.eta[FG::ELECTRON] - FDR<5,4>::get(F.eta[NUE]) ) ;

    const auto absorb_e = F.eta_np * abs_const / block_factor_n ;       
    
    if( std::isfinite(absorb_e) && (absorb_e > 0 ) )
        zeta_tmp[NUE] += absorb_e ;

    // nue_bar + p --> n + e+ 
    auto const block_factor_p = 
                         1. + exp( -F.eta[FG::ELECTRON] - FDR<5,4>::get(F.eta[NUE_BAR]) ) ;
    auto const absorb_a = F.eta_pn * abs_const / block_factor_p ;
    
    if( std::isfinite(absorb_a) && (absorb_a > 0 ) )
        zeta_tmp[NUE_BAR] += absorb_a ;

    // Go from zeta to kappa and include factors of T/m_e c^2 
    #pragma unroll(2)
    for( int i=NUE; i<NUX; ++i){
        // TODO: try out ILEAS "improved rates"
        auto nfac =  ::int_pow<2>( F.temp / me_mev ) *
                     FDR<4,2>::get(F.eta[i]) ; 
        kn[i] += zeta_tmp[i] * nfac ;
        auto efac =  ::int_pow<2>( F.temp / me_mev ) *
                     FDR<5,3>::get(F.eta[i]) ; 
        ka[i] += zeta_tmp[i] * efac ;
    }

    #ifdef DEBUG
    std::cout << "Absorption (Charged currents)" << std::endl;
    std::cout << "kappa_a : \n";
    for( int i=NUE; i<=NUX; ++i){
        std::cout << names[i] << '\t' << ka[i] << '\n';
    }
    std::cout << "\n\n" ;
    std::cout << "kappa_n : \n";
    for( int i=NUE; i<=NUX; ++i){
        std::cout << names[i] << '\t' << kn[i] << '\n';
    }
    std::cout << "\n\n" ;
    #endif

}

template<typename T>
inline __attribute__((always_inline))
void EAS<T>::add_scattering_opacity(std::array<T, NUMSPECIES>& ks) { 

    // Ruffert+ A1 (1996)
    constexpr auto const Cs_n = (1. + 5. * ::int_pow<2>(alpha) ) / 24. * sigma_0 ;
    constexpr auto const Cs_p = (4. * ::int_pow<2>(Cv-1.) + 5. * ::int_pow<2>(alpha)) / 24. * sigma_0 ;
    // Rosswog+ 2003 (A18)
    auto C_nucleus = 1. / 16. * sigma_0 * F.Ym[FG::ABAR] *
                     ::int_pow<2>(1. - F.Ym[FG::ZBAR] / F.Ym[FG::ABAR]);
    if (!std::isnormal(C_nucleus)) C_nucleus = 0;

    // Ruffert+ A8 (1996)   
    // Approximate Pauli blocking to account for degenerate matter
    auto ymp = F.Ym[FG::PROTON] / ( 1. + 2. / 3. * std::max(0., F.eta[FG::PROTON] ) ) ;
    auto ymn = F.Ym[FG::NEUTRON] / (1. + 2. / 3. * max(0., F.eta[FG::NEUTRON])); 

    if (!std::isnormal(ymp)) ymp = F.Ym[FG::PROTON];
    if (!std::isnormal(ymn)) ymn = F.Ym[FG::NEUTRON];

//Harry unroll(5)??    
    #pragma unroll(3)
    for( int i=NUE; i<=NUX; ++i) {
        auto efac =  ::int_pow<2>( F.temp / me_mev ) *
                     FDR<5,3>::get(F.eta[i]) ; 
        ks[i] += F.nb * Cs_n * ymn ;
        ks[i] += F.nb * Cs_p * ymp ;
        ks[i] += C_nucleus * F.nb * F.Ym[FG::HEAVY] ;
        ks[i] *= efac ;
    }

    #ifdef DEBUG
    std::cout << F.eta_np << std::endl;
    std::cout << "Scattering" << std::endl;
    std::cout << "kappa_s : \n";
    for( int i=NUE; i<=NUX; ++i){
        std::cout << names[i] << '\t' << ks[i] << '\n';
    }
    std::cout << "\n\n" ;
    #endif

}

// no numu numubar here, should use Weakhub
template<typename T>
inline __attribute__((always_inline))
void EAS<T>::add_charged_current_emission( std::array<T, NUMSPECIES>& QQ,
                                           std::array<T, NUMSPECIES>& RR) 
{
    // ! NOTE: heavy lepton neutrinos are left untouched by this
    // ! routine
    // TODO: might want to try ILEAS's "improved rates"
    arr_t block_factor { 1., 1., 1., 1., 1.} ;

    // Ruffert+ (B3-4)
    block_factor[NUE] = 
         1. + exp( F.eta[NUE] - FDR<5,4>::get(F.eta[FG::ELECTRON]) ) ; 

    block_factor[NUE_BAR] = 
         1. + exp( F.eta[NUE_BAR] - FDR<5,4>::get(-F.eta[FG::ELECTRON]) ) ;  

    // Ruffert+ (B1-2)
    auto const Re = beta * F.eta_pn * ::int_pow<5>( F.temp ) * FD<4>::get(F.eta[FG::ELECTRON]) / block_factor[NUE] ;
    auto const Ra = beta * F.eta_pn * ::int_pow<5>( F.temp ) * FD<4>::get(-F.eta[FG::ELECTRON]) / block_factor[NUE_BAR] ;

    // Ruffert+ (B15)
    auto const Qe = beta*F.eta_pn * ::int_pow<6>(F.temp) * FD<5>::get(F.eta[FG::ELECTRON]) / block_factor[NUE];
    auto const Qa = beta*F.eta_pn * ::int_pow<6>(F.temp) * FD<5>::get(-F.eta[FG::ELECTRON]) / block_factor[NUE_BAR];

    
    if( Re > 0 && std::isfinite(Re) && 
        Ra > 0 && std::isfinite(Ra) ) {    
        RR[NUE] += Re ;
        RR[NUE_BAR] += Ra ;
    }
       
    if( Qe > 0 && std::isfinite(Qe) && 
        Qa > 0 && std::isfinite(Qa) ) {    
        QQ[NUE] += Qe ;
        QQ[NUE_BAR] += Qa ;
    }   
    
    #ifdef DEBUG
    std::cout << "Beta Decays" << std::endl;
    std::cout << "Q : \n";
    std::cout << "NUE: " << Qe << std::endl ;
    std::cout << "NUE_BAR: " << Qa << std::endl ;
    std::cout << "\n\n" ;
    std::cout << "R : \n";
    std::cout << "NUE: " << Re << std::endl ;
    std::cout << "NUE_BAR: " << Ra << std::endl ;
    std::cout << "\n\n" ;
    #endif

}

template<typename T>
inline __attribute__((always_inline))
void EAS<T>::add_pair_process_emission(std::array<T, NUMSPECIES>& QQ,
                                           std::array<T, NUMSPECIES>& RR)
{
    arr_t block_factor{1.,1.,1.,1.,1.} ;
    // Ruffert+ (B9)
    #pragma unroll(3)
    for (int i = NUE; i <= NUX; ++i) {
      block_factor[i] =
          1. + exp(F.eta[i] - 0.5 * (FDR<4, 3>::get(F.eta[FG::ELECTRON]) +
                                     FDR<4, 3>::get(-F.eta[FG::ELECTRON])));
    }
    constexpr auto eps_const = 8. * pi / ::int_pow<3>(hc_mevcm);
    // Ruffert et al. (B5)
    const auto eps_m =
        eps_const * ::int_pow<4>(F.temp) * FD<3>::get(F.eta[FG::ELECTRON]);
    const auto eps_p =
        eps_const * ::int_pow<4>(F.temp) * FD<3>::get(-F.eta[FG::ELECTRON]);

    // Ruffert et al. (B6)
    const auto eps_tilde_m =
        eps_const * ::int_pow<5>(F.temp) * FD<4>::get(F.eta[FG::ELECTRON]);
    const auto eps_tilde_p =
        eps_const * ::int_pow<5>(F.temp) * FD<4>::get(-F.eta[FG::ELECTRON]);

    
    const auto eps_fraction = 0.5 * F.temp *
                              (FDR<4, 3>::get(F.eta[FG::ELECTRON]) +
                               FDR<4, 3>::get(-F.eta[FG::ELECTRON]));


    // Pair constant in Ruffert (B8)
    const auto pair_const =
        sigma_0 * clite / ::int_pow<2>(me_mev) * eps_m * eps_p;
    const auto R_pair = pair_const *
                        (::int_pow<2>(Cv - Ca) + ::int_pow<2>(Cv + Ca)) /
                        (36. * block_factor[NUE] * block_factor[NUE_BAR]);
    
    if (std::isfinite(R_pair)) {
//code test, even blocking factor used here, still overproducing nue nue bar, so dont use it
      //RR[NUE] += R_pair;
      //RR[NUE_BAR] += R_pair;

      //// Ruffert et al. (B16)
      //QQ[NUE] += R_pair * eps_fraction;
      //QQ[NUE_BAR] += R_pair * eps_fraction;
    }

    // numu, numu_bar, nux
    // Ruffert et al. (B10)
    if (EOS_Leptonic::use_muonic_eos) {
    	const auto R_pair_numu_numubar = 1. / 36. * pair_const *
    	                                  (::int_pow<2>(Cv - Ca) + ::int_pow<2>(Cv + Ca - 2.)) /
    	                                  (block_factor[NUMU] * block_factor[NUMU_BAR]);

    	const auto R_pair_x = 1. / 18. * pair_const *
    	                      (::int_pow<2>(Cv - Ca) + ::int_pow<2>(Cv + Ca - 2.)) /
    	                      ::int_pow<2>(block_factor[NUX]);

        // Ruffert et al. (B16)
        if (std::isfinite(R_pair_numu_numubar)) {
          RR[NUMU] += R_pair_numu_numubar;
          RR[NUMU_BAR] += R_pair_numu_numubar;
          QQ[NUMU] += R_pair_numu_numubar * eps_fraction;
          QQ[NUMU_BAR] += R_pair_numu_numubar * eps_fraction;
        }
        if (std::isfinite(R_pair_x)) {
          RR[NUX] += R_pair_x;
          QQ[NUX] += R_pair_x * eps_fraction;
        }
    } else {
    	const auto R_pair_x = 1. / 9. * pair_const *
    	                      (::int_pow<2>(Cv - Ca) + ::int_pow<2>(Cv + Ca - 2.)) /
    	                      ::int_pow<2>(block_factor[NUX]);
    	// Ruffert et al. (B16)
    	RR[NUMU] += 0.0;
    	RR[NUMU_BAR] += 0.0;
    	QQ[NUMU] += 0.0;
    	QQ[NUMU_BAR] += 0.0;
    	if (std::isfinite(R_pair_x)) {
    	  RR[NUX] += R_pair_x;
    	  QQ[NUX] += R_pair_x * eps_fraction;
    	}
    }

    #ifdef DEBUG
    std::cout << "Electrion Positron pair annihilation" << std::endl;
    std::cout << "Q : ";

    std::cout << R_pair * eps_fraction << " , ";
    std::cout << R_pair * eps_fraction << " , ";
    std::cout << R_pair_x * eps_fraction << " , ";
    std::cout << std::endl;

    std::cout << "R : ";
    std::cout << R_pair << " , ";
    std::cout << R_pair << " , ";
    std::cout << R_pair_x << " , ";

    std::cout << std::endl << std::endl ;
    #endif
}

template<typename T>
inline __attribute__((always_inline))
void EAS<T>::add_plasmon_decay_emission(std::array<T, NUMSPECIES>& QQ,
                                        std::array<T, NUMSPECIES>& RR)
{
    arr_t block_factor{1.,1.,1.,1.,1.} ;
    // Ruffert et al. see comment after (B12)
    const auto gamma =
    gamma_0 *
    sqrt(1. / 3. *
            (::int_pow<2>(pi) + 3.0 * ::int_pow<2>(F.eta[FG::ELECTRON])));
    
    // Ruffert et al. (B11)
    const auto gamma_const =
    ::int_pow<3>(pi) * sigma_0 * clite /
    (::int_pow<2>(me_mev) * 3.0 * fsc * ::int_pow<6>(hc_mevcm)) *
    ::int_pow<6>(gamma) * exp(-gamma) * (1.0 + gamma) * ::int_pow<8>(F.temp);

    #pragma unroll(3)
    for (int i = NUE; i <= NUX; ++i) {   
        block_factor[i] =
          1. + exp(F.eta[i] - (1. + 0.5 * ::int_pow<2>(gamma) / (1. + gamma)));
    }   

    // Ruffert et al. (B11)
    const auto R_gamma =
    ::int_pow<2>(Cv) * gamma_const / (block_factor[NUE] * block_factor[NUE_BAR]);
    const auto Q_gamma = F.temp * 0.5 * (2. + ::int_pow<2>(gamma) / (1. + gamma));
    
    if (std::isfinite(R_gamma)) {
// even blocking factor here for nue/nue_bar still overproducing, dont use
        //RR[NUE] += R_gamma ;
        //RR[NUE_BAR] += R_gamma ;
        //QQ[NUE] += Q_gamma * R_gamma;
        //QQ[NUE_BAR] += Q_gamma * R_gamma;
    }

    // Ruffert et al. (B12)
    if (EOS_Leptonic::use_muonic_eos) {
        const auto R_gamma_numu =
        ::int_pow<2>(Cv - 1.) * gamma_const / (block_factor[NUMU] * block_factor[NUMU_BAR]);
        if (std::isfinite(R_gamma_numu)) {
                RR[NUMU] += R_gamma_numu;
                QQ[NUMU] += Q_gamma * R_gamma_numu;
        }
        const auto R_gamma_numu_bar =
        ::int_pow<2>(Cv - 1.) * gamma_const / (block_factor[NUMU] * block_factor[NUMU_BAR]);
        if (std::isfinite(R_gamma_numu_bar)) {
                RR[NUMU_BAR] += R_gamma_numu_bar;
                QQ[NUMU_BAR] += Q_gamma * R_gamma_numu_bar;
        }
        // nux only counts nutau and nutau_bar now
        const auto R_gamma_x =
        2. * ::int_pow<2>(Cv - 1.) * gamma_const / ::int_pow<2>(block_factor[NUX]);
        if (std::isfinite(R_gamma_x)) {
                RR[NUX] += R_gamma_x;
                QQ[NUX] += Q_gamma * R_gamma_x;
        }
    } else {
    	const auto R_gamma_x =
    	4. * ::int_pow<2>(Cv - 1.) * gamma_const / ::int_pow<2>(block_factor[NUX]);
    	if (std::isfinite(R_gamma_x)) {
    	        RR[NUX] += R_gamma_x;
    	        // Ruffert et al. (B17)
    	        QQ[NUX] += Q_gamma * R_gamma_x;
    	}
    }
    #ifdef DEBUG
    std::cout << "Plasmon decay" << std::endl;
    std::cout << "Q : ";

    std::cout << R_gamma * Q_gamma << " , ";
    std::cout << R_gamma * Q_gamma << " , ";
    std::cout << R_gamma_x * Q_gamma << " , ";
    std::cout << std::endl;

    std::cout << "R : ";
    std::cout << R_gamma << " , ";
    std::cout << R_gamma << " , ";
    std::cout << R_gamma_x << " , ";

    std::cout << std::endl << std::endl ;
    #endif
}

template<typename T>
inline __attribute__((always_inline))
void EAS<T>::add_brems_emission(std::array<T, NUMSPECIES>& QQ,
                                std::array<T, NUMSPECIES>& RR)
{

    // ! Note: this only affects heavy lepton neutrinos ! 
    const auto R_brems =
          0.231 * (2.0778e2 * erg_to_mev) * 0.5 *
          (::int_pow<2>(F.Ym[FG::NEUTRON]) + ::int_pow<2>(F.Ym[FG::PROTON]) +
           28.0 / 3.0 * F.Ym[FG::NEUTRON] * F.Ym[FG::PROTON]) *
          ::int_pow<2>(F.rho) * ::int_pow<4>(F.temp) * sqrt(F.temp);

    const auto Q_brems = R_brems * F.temp / 0.231 * 0.504;
    // will it over producing numu and numubar in npemu matter?
    if ( std::isfinite( Q_brems ) && ( Q_brems > 0 ) ) {
        if (EOS_Leptonic::use_muonic_eos) {
		QQ[NUMU] += g_nu_5spec[NUMU] * Q_brems;
		RR[NUMU] += g_nu_5spec[NUMU] * R_brems;
		QQ[NUMU_BAR] += g_nu_5spec[NUMU_BAR] * Q_brems;
		RR[NUMU_BAR] += g_nu_5spec[NUMU_BAR] * R_brems;
        	QQ[NUX] += g_nu_5spec[NUX] * Q_brems ;
        	RR[NUX] += g_nu_5spec[NUX] * R_brems  ;
        } else {
		QQ[NUMU] += 0.0;
		RR[NUMU] += 0.0;
		QQ[NUMU_BAR] += 0.0;
		RR[NUMU_BAR] += 0.0;
        	QQ[NUX] += g_nu_3spec[NUX] * Q_brems ;
        	RR[NUX] += g_nu_3spec[NUX] * R_brems  ;
        }
    }
    #ifdef DEBUG
      std::cout << "Nucleon-Nuclean Bremsstrahlung" << std::endl;
      std::cout << "Q : ";
      std::cout << Q_brems << " , ";
      std::cout << std::endl;

      std::cout << "R : ";
      std::cout << R_brems  << " , ";

      std::cout << std::endl;

    #endif
}
// only for heavy leptons
template<typename T>
inline __attribute__((always_inline))
void EAS<T>::add_Kirchoff_absorption_opacity(std::array<T, NUMSPECIES>& ka,
                                             std::array<T, NUMSPECIES>& kn,
                                             const std::array<T, NUMSPECIES>& QQ,
                                             const std::array<T, NUMSPECIES>& RR, double& rho_cgs)
{
    auto const k_a_thermo = QQ[NUX] / black_body_mev<NUX,1>() ;
    auto const k_n_thermo = RR[NUX] / black_body_mev<NUX,0>() ;

    if (EOS_Leptonic::use_muonic_eos) {
        // code test:  if rho < 3e11, set BB of eta_numu and eta_numu = eta_nux
        // To prevent over de-muonization with eta(-100) 
        double black_body_mev_1 = 0.0;
        double black_body_mev_0 = 0.0;
	if (rho_cgs < 3.e0) {
	  black_body_mev_1 = black_body_mev<NUX,1>();
	  black_body_mev_0 = black_body_mev<NUX,0>();
        } else {
	  black_body_mev_1 = black_body_mev<NUMU,1>();
	  black_body_mev_0 = black_body_mev<NUMU,0>();
        }
    	auto const k_a_thermo_numu = QQ[NUMU] / black_body_mev_1 ;
    	auto const k_n_thermo_numu = RR[NUMU] / black_body_mev_0 ;
        if (rho_cgs < 3.e0) {
          black_body_mev_1 = black_body_mev<NUX,1>();
          black_body_mev_0 = black_body_mev<NUX,0>();
        } else { 
          black_body_mev_1 = black_body_mev<NUMU_BAR,1>();
          black_body_mev_0 = black_body_mev<NUMU_BAR,0>();
        }
    	auto const k_a_thermo_numu_bar = QQ[NUMU_BAR] / black_body_mev_1 ;
    	auto const k_n_thermo_numu_bar = RR[NUMU_BAR] / black_body_mev_0 ;

	if( std::isfinite(k_a_thermo_numu) && k_a_thermo_numu > 0 &&
	    std::isfinite(k_n_thermo_numu) && k_n_thermo_numu > 0 )
	{
	    ka[NUMU] = k_a_thermo_numu ;
	    kn[NUMU] = k_n_thermo_numu ;
	}

	if( std::isfinite(k_a_thermo_numu_bar) && k_a_thermo_numu_bar > 0 &&
	    std::isfinite(k_n_thermo_numu_bar) && k_n_thermo_numu_bar > 0 )
	{
	    ka[NUMU_BAR] = k_a_thermo_numu_bar ;
	    kn[NUMU_BAR] = k_n_thermo_numu_bar ;
	}
    }

    if( std::isfinite(k_a_thermo) && k_a_thermo > 0 &&
        std::isfinite(k_n_thermo) && k_n_thermo > 0 ) 
    {
        ka[NUX] = k_a_thermo ;
        kn[NUX] = k_n_thermo ;
    }
    #ifdef DEBUG
    std::cout << "Absorption (Kirchoff law)" << std::endl;
    std::cout << "kappa_a : \n";
    std::cout << names[NUX] << '\t' << k_a_thermo << '\n';
    std::cout << "\n\n" ;
    std::cout << "kappa_n : \n";
    std::cout << names[NUX] << '\t' << k_n_thermo << '\n';
    std::cout << "\n\n" ;
    #endif
}

template<typename T>
inline __attribute__((always_inline))
void EAS<T>::calc_Kirchoff_emission(std::array<T, NUMSPECIES>& QQ,
                                             std::array<T, NUMSPECIES>& RR,
                                             const std::array<T, NUMSPECIES>& ka,
                                             const std::array<T, NUMSPECIES>& kn, 
                                             bool need_nux_emission)
{
    auto Q_thermo = ka[NUE]  * black_body_mev<NUE,1>() ;
    auto R_thermo = kn[NUE]  * black_body_mev<NUE,0>() ;

    if( std::isfinite(Q_thermo) && Q_thermo > 0 &&
        std::isfinite(R_thermo) && R_thermo > 0 ) 
    {
        QQ[NUE] = Q_thermo ;
        RR[NUE] = R_thermo ;
    }

    Q_thermo = ka[NUE_BAR]  * black_body_mev<NUE_BAR,1>() ;
    R_thermo = kn[NUE_BAR]  * black_body_mev<NUE_BAR,0>() ;

    if( std::isfinite(Q_thermo) && Q_thermo > 0 &&
        std::isfinite(R_thermo) && R_thermo > 0 ) 
    {
        QQ[NUE_BAR] = Q_thermo ;
        RR[NUE_BAR] = R_thermo ;
    }

    // Weakhub needs this by providing kappa_a/n of nux
    if (need_nux_emission) {
        if (EOS_Leptonic::use_muonic_eos) {
	        Q_thermo = ka[NUMU] * black_body_mev<NUMU,1>() ;
	        R_thermo = kn[NUMU] * black_body_mev<NUMU,0>() ;
	        if( std::isfinite(Q_thermo) && Q_thermo > 0 &&
	            std::isfinite(R_thermo) && R_thermo > 0 )
	        {
	            QQ[NUMU] = Q_thermo ;
	            RR[NUMU] = R_thermo ;
	        }

                Q_thermo = ka[NUMU_BAR] * black_body_mev<NUMU_BAR,1>() ;
                R_thermo = kn[NUMU_BAR] * black_body_mev<NUMU_BAR,0>() ;
                if( std::isfinite(Q_thermo) && Q_thermo > 0 &&
                    std::isfinite(R_thermo) && R_thermo > 0 )
                {
                    QQ[NUMU_BAR] = Q_thermo ;
                    RR[NUMU_BAR] = R_thermo ;
                }

        }

      	Q_thermo = ka[NUX] * black_body_mev<NUX,1>() ;
      	R_thermo = kn[NUX] * black_body_mev<NUX,0>() ;
      	if( std::isfinite(Q_thermo) && Q_thermo > 0 &&
      	    std::isfinite(R_thermo) && R_thermo > 0 )
      	{
      	    QQ[NUX] = Q_thermo ;
      	    RR[NUX] = R_thermo ;
      	}
    }

    #ifdef DEBUG
    std::cout << "Emission (Kirchoff law)" << std::endl;
    std::cout << "Q : \n";
    for( int i=NUE; i<NUX; ++i){
        std::cout << names[i] << '\t' << QQ[i] << '\n';
    }
    std::cout << "\n\n" ;
    std::cout << "R : \n";
    for( int i=NUE; i<NUX; ++i){
        std::cout << names[i] << '\t' << RR[i] << '\n';
    }
    std::cout << "\n\n" ;
    #endif
}

// Weakhub related, Weakhub tables are in cgs units, but entries need rho to be code unit
template<typename T>
inline __attribute__((always_inline))
void EAS<T>::kappa_ast_kappa_s__temp_rho_yle_ymu_M1(
      std::array<T, NUMSPECIES>& kappa_a_st_en, std::array<T, NUMSPECIES>& kappa_a_st_num, 
      std::array<T, NUMSPECIES>& kappa_s, int num_spec, double &temp,
      double &rho, double &yle, double &ymu){

  // the input rho is cgs unit

  double lrho = log(rho * RHOGF);
  double ltemp = log(temp);

  //code test it passed, it means, opacity reading has problem! FIXME
  //   for (int i = NUE; i < NUX; i++) {
  //    kappa_a_st_en[i] = 0.0;
  //    kappa_a_st_num[i] = 0.0;
  //    kappa_s[i] = 0.0;
  //   }
  //   return;

  // rho
  if ( rho * RHOGF <= exp(M1_Weakhub::logrho_min_IV) ) {
     for (int i = NUE; i <= NUX; i++) {
        kappa_a_st_en[i] = 1.e-60;
        kappa_a_st_num[i] = 1.e-60;
        kappa_s[i] = 1.e-60;
     }
     return;
  } else if (rho * RHOGF > exp(M1_Weakhub::logrho_max_IV) ) {
   	//std::cout << "lrho ,  M1_Weakhub::logrho_max_IV" << "\n";
   	//std::cout << lrho << ", " << M1_Weakhub::logrho_max_IV << "\n";
   	//CCTK_ERROR("lrho  > maximum in Weakhub table") ;
   	lrho = std::min( M1_Weakhub::logrho_max_IV, lrho);
  }

  // temp
  if (temp > exp(M1_Weakhub::logtemp_max_IV) ) {
     	//std::cout << "ltemp , M1_Weakhub::logtemp_max_IV" << "\n";
     	//std::cout << ltemp << ", " << M1_Weakhub::logtemp_max_IV << "\n";
     	//CCTK_ERROR("ltemp > maximum in Weakhub table") ;
        ltemp = std::min( M1_Weakhub::logtemp_max_IV, ltemp);
  } else if (temp <= exp(M1_Weakhub::logtemp_min_IV) ) {
        ltemp = std::max( M1_Weakhub::logtemp_min_IV, ltemp);
  }

  // yle
  if (yle > M1_Weakhub::ye_max_IV) { 
     yle = M1_Weakhub::ye_max_IV;
  }

  if (yle < M1_Weakhub::ye_min_IV) { 
     yle = M1_Weakhub::ye_min_IV;
  }

  // ymu
  if (EOS_Leptonic::use_muonic_eos) {
 	 if (ymu > exp(M1_Weakhub::logymu_max_IV) ) {
 	    ymu = exp(M1_Weakhub::logymu_max_IV);
 	 } 

 	 if (ymu < exp(M1_Weakhub::logymu_min_IV) ) {
 	    ymu = exp(M1_Weakhub::logymu_min_IV);
 	 }
  }

  // For safety code test:
  lrho = std::max( M1_Weakhub::logrho_min_IV, std::min(lrho, M1_Weakhub::logrho_max_IV) );
  ltemp = std::max( M1_Weakhub::logtemp_min_IV, std::min(ltemp, M1_Weakhub::logtemp_max_IV) );

  if (num_spec == 3) {
     auto const kappa_a_st_en_local =
        M1_Weakhub::kappa_a_st_en_table_3spec.interpolate<M1_Weakhub::NEU_NUM::i_nu_e, M1_Weakhub::NEU_NUM::i_nu_e_bar, 
           M1_Weakhub::NEU_NUM::i_nu_mu>(lrho, ltemp, yle);
     auto const kappa_a_st_num_local =
        M1_Weakhub::kappa_a_st_num_table_3spec.interpolate<M1_Weakhub::NEU_NUM::i_nu_e, M1_Weakhub::NEU_NUM::i_nu_e_bar, 
           M1_Weakhub::NEU_NUM::i_nu_mu>(lrho, ltemp, yle);
     auto const kappa_s_local =
        M1_Weakhub::kappa_s_table_3spec.interpolate<M1_Weakhub::NEU_NUM::i_nu_e, M1_Weakhub::NEU_NUM::i_nu_e_bar, 
           M1_Weakhub::NEU_NUM::i_nu_mu>(lrho, ltemp, yle);
     for (int i = NUE; i <= NUE_BAR; i++) {
       kappa_a_st_en[i]  = kappa_a_st_en_local[i];
       kappa_a_st_num[i] = kappa_a_st_num_local[i];
       kappa_s[i]        = kappa_s_local[i];
     }
     //Marie: ok, i_nu_mu=2 and in 3spec-table nux=2
     kappa_a_st_en[NUMU] = 1.e-60;
     kappa_a_st_en[NUMU_BAR] = 1.e-60;
     kappa_a_st_en[NUX] = kappa_a_st_en_local[2];
     kappa_a_st_num[NUMU] = 1.e-60;
     kappa_a_st_num[NUMU_BAR] = 1.e-60;
     kappa_a_st_num[NUX] = kappa_a_st_num_local[2];
     kappa_s[NUMU] = 1.e-60;
     kappa_s[NUMU_BAR] = 1.e-60;
     kappa_s[NUX] = kappa_s_local[2];
     // kappa_a_nux no needs to times 4 in this code, only emissivty needs
     //kappa_a_st_en[NUX]  *= 4.0;
     //kappa_a_st_num[NUX]  *= 4.0;
     // Harry: should I * 4 for kappa_s_nux?
  }
  else if (num_spec == 5) {
    // this option treats kappa[i_nu_tau] = kappa[i_nu_tau_bar]
    double lymu = log(ymu);
    lymu = std::max( M1_Weakhub::logymu_min_IV, std::min(lymu, M1_Weakhub::logymu_max_IV) );

    auto const kappa_a_st_en_local =
       M1_Weakhub::kappa_a_st_en_table_6spec.interpolate
         <M1_Weakhub::NEU_NUM::i_nu_e, M1_Weakhub::NEU_NUM::i_nu_e_bar, M1_Weakhub::NEU_NUM::i_nu_mu, 
          M1_Weakhub::NEU_NUM::i_nu_mu_bar, M1_Weakhub::NEU_NUM::i_nu_tau>(lrho, ltemp, yle, lymu);
    auto const kappa_a_st_num_local =
       M1_Weakhub::kappa_a_st_num_table_6spec.interpolate
         <M1_Weakhub::NEU_NUM::i_nu_e, M1_Weakhub::NEU_NUM::i_nu_e_bar, M1_Weakhub::NEU_NUM::i_nu_mu, 
          M1_Weakhub::NEU_NUM::i_nu_mu_bar, M1_Weakhub::NEU_NUM::i_nu_tau>(lrho, ltemp, yle, lymu);
    auto const kappa_s_local =
       M1_Weakhub::kappa_s_table_6spec.interpolate
         <M1_Weakhub::NEU_NUM::i_nu_e, M1_Weakhub::NEU_NUM::i_nu_e_bar, M1_Weakhub::NEU_NUM::i_nu_mu, 
          M1_Weakhub::NEU_NUM::i_nu_mu_bar, M1_Weakhub::NEU_NUM::i_nu_tau>(lrho, ltemp, yle, lymu);
     for (int i = NUE; i <= NUX; i++) {
       kappa_a_st_en[i]  = kappa_a_st_en_local[i];
       kappa_a_st_num[i] = kappa_a_st_num_local[i];
       kappa_s[i]        = kappa_s_local[i];
     }
     // kappa_a_nux no needs to times 2 in this code, only emissivity needs!
     //kappa_a_st_en[NUX]  *= 2.0;
     //kappa_a_st_num[NUX]  *= 2.0;
  }
  else {
     CCTK_ERROR("\n Unknown num_spec inside kappa_ast_kappa_s (Weakhub) in side Margherita_M1 \n") ;
  }

  // Ensure opacity >= 0.0 and non NaN
  for (int i = NUE; i <= NUX; i++) {
   kappa_a_st_en[i] = std::max(kappa_a_st_en[i], 1.e-60);
   kappa_a_st_num[i] = std::max(kappa_a_st_num[i], 1.e-60);
   kappa_s[i] = std::max(kappa_s[i], 1.e-60);
   if (kappa_a_st_en[i] != kappa_a_st_en[i]) {kappa_a_st_en[i] = 1.e-60;}
   if (kappa_a_st_num[i] != kappa_a_st_num[i]) {kappa_a_st_num[i] = 1.e-60;}
   if (kappa_s[i] != kappa_s[i]) {kappa_s[i] = 1.e-60;}
   
  }
  return;

} // End of Weakhub table interpolation

} // End of Namespace 
#endif
