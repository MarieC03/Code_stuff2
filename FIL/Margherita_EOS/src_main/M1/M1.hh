#include "../margherita.hh"
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
enum m1_species_t { NUE=0, NUA, NUX, NUMSPECIES} ;
enum m1_fluid_index_t { RHO=0, YE, TEMP, REQ_FLUID_VARS} ;
#endif 

namespace Margherita_M1_EAS { 
using namespace m1_helpers; 
using namespace M1_Constants;
using namespace Margherita_helpers;
using namespace Margherita_constants;

#ifdef DEBUG
const std::array<std::string, 3> names { 
    "NUE",
    "NUA",
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
    std::array<T,3> const g_nu {1, 1, 4} ;
    // which interactions are active
    std::array<bool, NINTERACTIONS> which_interactions   ;
    // have the rates been computed yet? 
    bool computed = false ; 

    FG F  ;

    // Rates

    // Emission: MeV / s / cm^3
    arr_t Q {0., 0., 0.} ;
    arr_t Q_diff_inv {0., 0., 0.} ;
    arr_t Q_eff {0., 0., 0.} ;

    // Number emission: number / s / cm^3
    arr_t R {0., 0., 0.} ;
    arr_t R_diff_inv {0., 0., 0.} ;
    arr_t R_eff {0., 0., 0.} ;


    arr_t tau_n {0., 0., 0.} ;
    arr_t tau_e {0., 0., 0.} ;

    // Mean energy (squared): MeV^2
    arr_t eps2_leak {0., 0., 0.} ;
    arr_t eps2_fluid {0., 0., 0.} ;
    arr_t eps {0., 0., 0.} ; // MeV 

    // Opacities: 1/cm
    arr_t kappa_s {0., 0., 0.} ;
    arr_t kappa_a {0., 0., 0.} ;
    arr_t kappa_n {0., 0., 0.} ;

    arr_t Y_nu {0., 0., 0.} ;
    arr_t Z_nu {0., 0., 0.} ;

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
        const T &ye, 
        const T &temp, 
        std::array<bool,NINTERACTIONS> const& wi)
    : F(FG(rho, temp, ye)),  which_interactions(wi) {};
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
  template< int i >
    inline decltype(auto) Q_mev() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        assert(computed) ;
        return Q[i] ; 
    }
    // Number Emissivities 
    template< int i >
    inline decltype(auto) R_cactus() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        assert(computed) ;
        return R[i] * mnuc_cgs * RHOGF / TIMEGF ; 
    }
  template< int i >
    inline decltype(auto) R_mev() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        assert(computed) ;
        return R[i] ; 
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
    inline decltype(auto) kappa_a_mev() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        assert(computed) ;
        return kappa_a[i]; 
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
      template< int i >
    inline decltype(auto) kappa_n_mev( ) const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        assert(computed) ;
        return kappa_n[i] ; 
    }
    // BlackBody functions 
    template< int i, int n> 
    inline __attribute__((always_inline)) 
    decltype(auto) black_body_mev() const {
        static_assert( !(i<0) && (i<NUMSPECIES) ) ;
        static_assert( (n==0) || (n==1) ) ;
        auto const bb = g_nu[i] * 4.*pi*clite/::int_pow<3>(hc_mevcm) * ::int_pow<3+n>(F.temp) * FD<2+n>::get(F.eta[i]) ; 
        return std::max(bb, 1.e-30)  ;
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
                                          const arr_t& RR ) ;

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


    // ! NOTE: this OVERWRITES the rates of NUE, NUA 
    inline __attribute__((always_inline))
    void calc_Kirchoff_emission(arr_t& QQ,
                                arr_t& RR,
                                const arr_t& ka,
                                const arr_t& kn ) ;


    // Compute leakage effective rates and 
    // mean energies 
    inline __attribute__((always_inline))
    void calc_diffusion_rates() {
        constexpr auto factor = 3. / a_diff * clite / ::int_pow<3>(hc_mevcm) * normfact ;
        // Rosswog+ 2002, Sekiguchi 2002
        for(int ii=NUE; ii<=NUX; ++ii){
            const auto E_sq = ::int_pow<2>(F.temp) * FDR<4,2>::get(F.eta[ii]) ; // MeV
            R_diff_inv[ii] = 
                ::int_pow<2>( tau_n[ii] ) / ( factor*g_nu[ii] * kappa_n[ii] * F.temp * FD<0>::get(F.eta[ii]) * E_sq ) ;

            const auto E_sq_e = ::int_pow<2>(F.temp) * FDR<5,3>::get(F.eta[ii]) ; // MeV
            Q_diff_inv[ii] = 
                ::int_pow<2>(tau_n[ii]) / (factor*g_nu[ii]*kappa_a[ii]*::int_pow<2>(F.temp)*FD<1>::get(F.eta[ii]) * E_sq_e ) ;
        }
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
    // Computation of the rates
    // ========================================================
    inline void calc_eas(bool recompute_emission_from_absorption=true) {
        Q = {0., 0., 0.}; R={0.,0.,0.} ;
        kappa_s = {0., 0., 0.}; kappa_a = {0.,0.,0.}; kappa_n = {0., 0., 0.}; 
        // free absorption rates 
        if ( which_interactions[BETA] )
            add_charged_current_absorption_opacity(kappa_a, kappa_n) ;
        // scattering 
        add_scattering_opacity(kappa_s) ;

        // Emission rates: 
        arr_t Q_beta{0., 0., 0.}, R_beta{0., 0., 0.} ;
        if( which_interactions[BETA] )
            add_charged_current_emission(Q_beta, R_beta) ;
        
        if ( which_interactions[PLASMON_DECAY] ) 
            add_plasmon_decay_emission(Q, R) ;
        if ( which_interactions[BREMS] )
            add_brems_emission(Q,R) ;
        if ( which_interactions[PAIR])
            add_pair_process_emission(Q, R) ;
        
        arr_t kappa_a_tmp {0., 0., 0. } ;
        arr_t kappa_n_tmp {0., 0., 0. } ;

        // thermal absorption for nu_x 
        add_Kirchoff_absorption_opacity(kappa_a_tmp, kappa_n_tmp, Q, R) ;

        if ( recompute_emission_from_absorption )
            calc_Kirchoff_emission(Q, R, kappa_a, kappa_n) ;

        for( int i=NUE; i<=NUX; i++ ) {
            kappa_a[i] += kappa_a_tmp[i];
            kappa_n[i] += kappa_n_tmp[i];
        }

	calc_energies() ;

        #ifdef DEBUG
        std::cout << "Total" << std::endl;
        std::cout << "kappa_a : ";

        std::cout << kappa_a[0] << " , ";
        std::cout << kappa_a[1] << " , ";
        std::cout << kappa_a[2] << " , ";
        std::cout << std::endl;

        std::cout << "Q : ";

        std::cout << Q[0] << " , ";
        std::cout << Q[1] << " , ";
        std::cout << Q[2] << " , ";
        std::cout << std::endl;

        std::cout << "kappa_n : ";

        std::cout << kappa_n[0] << " , ";
        std::cout << kappa_n[1] << " , ";
        std::cout << kappa_n[2] << " , ";
        std::cout << std::endl;

        std::cout << "R : ";
        std::cout << R[0] << " , ";
        std::cout << R[1] << " , ";
        std::cout << R[2] << " , ";

        std::cout << std::endl;

        std::cout << "Total [code units]" << std::endl;
        std::cout << "kappa_a : ";

        std::cout << kappa_a[0] / LENGTHGF<< " , ";
        std::cout << kappa_a[1] / LENGTHGF<< " , ";
        std::cout << kappa_a[2] / LENGTHGF<< " , ";
        std::cout << std::endl;

        std::cout << "Q : ";

        std::cout << Q[0] * mev_to_erg * EPSGF * RHOGF / TIMEGF << " , ";
        std::cout << Q[1] * mev_to_erg * EPSGF * RHOGF / TIMEGF << " , ";
        std::cout << Q[2] * mev_to_erg * EPSGF * RHOGF / TIMEGF << " , ";
        std::cout << std::endl;
        
        std::cout << "kappa_n : ";

        std::cout << kappa_n[0] / LENGTHGF<< " , ";
        std::cout << kappa_n[1] / LENGTHGF<< " , ";
        std::cout << kappa_n[2] / LENGTHGF<< " , ";
        std::cout << std::endl;

        std::cout << "R : ";
        std::cout << R[0] * EPSGF * RHOGF / TIMEGF << " , ";
        std::cout << R[1] * EPSGF * RHOGF / TIMEGF << " , ";
        std::cout << R[2] * EPSGF * RHOGF / TIMEGF << " , ";

        std::cout << std::endl;
        #endif 

        computed = true ;         
    }

    #ifdef DEBUG 
    #endif
};

template<typename T>
inline __attribute__((always_inline))
void EAS<T>::add_charged_current_absorption_opacity(std::array<T, NUMSPECIES>& ka, 
                                                    std::array<T, NUMSPECIES>& kn ) { 

    // ! Note that heavy lepton neutrinos do not undergo these processes 
    // ! ( We don't have data about muons and tau in the EOS tables )                                             
    arr_t zeta_tmp {0., 0., 0.} ;
    
    constexpr auto const abs_const = 0.25 * ( 1. + 3. * ::int_pow<2>(alpha) ) * sigma_0 ; // ? alpha or g_A ? 
    
    // nue + n --> p + e- 
    auto const block_factor_n = 
                         1. + exp( F.eta[FG::ELECTRON] - FDR<5,4>::get(F.eta[NUE]) ) ;

    const auto absorb_e = F.eta_np * abs_const / block_factor_n ;       
    
    if( std::isfinite(absorb_e) && (absorb_e > 0 ) )
        zeta_tmp[NUE] += absorb_e ;

    // nua + p --> n + e+ 
    auto const block_factor_p = 
                         1. + exp( -F.eta[FG::ELECTRON] - FDR<5,4>::get(F.eta[NUA]) ) ;
    auto const absorb_a = F.eta_pn * abs_const / block_factor_p ;
    
    if( std::isfinite(absorb_a) && (absorb_a > 0 ) )
        zeta_tmp[NUA] += absorb_a ;

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

template<typename T>
inline __attribute__((always_inline))
void EAS<T>::add_charged_current_emission( std::array<T, NUMSPECIES>& QQ,
                                           std::array<T, NUMSPECIES>& RR) 
{
    // ! NOTE: heavy lepton neutrinos are left untouched by this
    // ! routine
    // TODO: might want to try ILEAS's "improved rates"
    arr_t block_factor { 1., 1., 1. } ;

    // Ruffert+ (B3-4)
    block_factor[NUE] = 
         1. + exp( F.eta[NUE] - FDR<5,4>::get(F.eta[FG::ELECTRON]) ) ; 

    block_factor[NUA] = 
         1. + exp( F.eta[NUA] - FDR<5,4>::get(-F.eta[FG::ELECTRON]) ) ;  

    // Ruffert+ (B1-2)
    auto const Re = beta * F.eta_pn * ::int_pow<5>( F.temp ) * FD<4>::get(F.eta[FG::ELECTRON]) / block_factor[NUE] ;
    auto const Ra = beta * F.eta_pn * ::int_pow<5>( F.temp ) * FD<4>::get(-F.eta[FG::ELECTRON]) / block_factor[NUA] ;

    // Ruffert+ (B15)
    auto const Qe = beta*F.eta_pn * ::int_pow<6>(F.temp) * FD<5>::get(F.eta[FG::ELECTRON]) / block_factor[NUE];
    auto const Qa = beta*F.eta_pn * ::int_pow<6>(F.temp) * FD<5>::get(-F.eta[FG::ELECTRON]) / block_factor[NUA];

    
    if( Re > 0 && std::isfinite(Re) && 
        Ra > 0 && std::isfinite(Ra) ) {    
        RR[NUE] += Re ;
        RR[NUA] += Ra ;
    }
       
    if( Qe > 0 && std::isfinite(Qe) && 
        Qa > 0 && std::isfinite(Qa) ) {    
        QQ[NUE] += Qe ;
        QQ[NUA] += Qa ;
    }   
    
    #ifdef DEBUG
    std::cout << "Beta Decays" << std::endl;
    std::cout << "Q : \n";
    std::cout << "NUE: " << Qe << std::endl ;
    std::cout << "NUA: " << Qa << std::endl ;
    std::cout << "\n\n" ;
    std::cout << "R : \n";
    std::cout << "NUE: " << Re << std::endl ;
    std::cout << "NUA: " << Ra << std::endl ;
    std::cout << "\n\n" ;
    #endif

}

template<typename T>
inline __attribute__((always_inline))
void EAS<T>::add_pair_process_emission(std::array<T, NUMSPECIES>& QQ,
                                           std::array<T, NUMSPECIES>& RR)
{
    arr_t block_factor{1.,1.,1.} ;
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
                        (36. * block_factor[NUE] * block_factor[NUA]);
    
    if (std::isfinite(R_pair)) {
      RR[NUE] += R_pair;
      RR[NUA] += R_pair;

      // Ruffert et al. (B16)
      QQ[NUE] += R_pair * eps_fraction;
      QQ[NUA] += R_pair * eps_fraction;
    }

    // Ruffert et al. (B10)
    const auto R_pair_x = 1. / 9. * pair_const *
                          (::int_pow<2>(Cv - Ca) + ::int_pow<2>(Cv + Ca - 2.)) /
                          ::int_pow<2>(block_factor[NUX]);

    if (std::isfinite(R_pair_x)) {
      RR[NUX] += R_pair_x;
      // Ruffert et al. (B16)
      QQ[NUX] += R_pair_x * eps_fraction;
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
    arr_t block_factor{1.,1.,1.} ;
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
    ::int_pow<2>(Cv) * gamma_const / (block_factor[NUE] * block_factor[NUA]);
    const auto Q_gamma = F.temp * 0.5 * (2. + ::int_pow<2>(gamma) / (1. + gamma));
    
    if (std::isfinite(R_gamma)) {
        RR[NUE] += R_gamma ;
        RR[NUA] += R_gamma ;
        QQ[NUE] += Q_gamma ;
        QQ[NUA] += Q_gamma ;
    }

    // Ruffert et al. (B12)
    const auto R_gamma_x =
    4. * ::int_pow<2>(Cv - 1.) * gamma_const / ::int_pow<2>(block_factor[NUX]);
    if (std::isfinite(R_gamma_x)) {
	    RR[NUX] += R_gamma_x;
	    // Ruffert et al. (B17)
	    QQ[NUX] += Q_gamma * R_gamma_x;
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
    if ( std::isfinite( Q_brems ) && ( Q_brems > 0 ) ) {
        QQ[NUX] += g_nu[NUX] * Q_brems ;
        RR[NUX] += g_nu[NUX] * R_brems  ;
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

template<typename T>
inline __attribute__((always_inline))
void EAS<T>::add_Kirchoff_absorption_opacity(std::array<T, NUMSPECIES>& ka,
                                             std::array<T, NUMSPECIES>& kn,
                                             const std::array<T, NUMSPECIES>& QQ,
                                             const std::array<T, NUMSPECIES>& RR)
{
    auto const k_a_thermo = QQ[NUX] / black_body_mev<NUX,1>() ;
    auto const k_n_thermo = RR[NUX] / black_body_mev<NUX,0>() ;

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
                                             const std::array<T, NUMSPECIES>& kn)
{
    auto Q_thermo = ka[NUE]  * black_body_mev<NUE,1>() ;
    auto R_thermo = kn[NUE]  * black_body_mev<NUE,0>() ;

    if( std::isfinite(Q_thermo) && Q_thermo > 0 &&
        std::isfinite(R_thermo) && R_thermo > 0 ) 
    {
        QQ[NUE] = Q_thermo ;
        RR[NUE] = R_thermo ;
    }

    Q_thermo = ka[NUA]  * black_body_mev<NUA,1>() ;
    R_thermo = kn[NUA]  * black_body_mev<NUA,0>() ;

    if( std::isfinite(Q_thermo) && Q_thermo > 0 &&
        std::isfinite(R_thermo) && R_thermo > 0 ) 
    {
        QQ[NUA] = Q_thermo ;
        RR[NUA] = R_thermo ;
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
} // End of Namespace 
#endif
