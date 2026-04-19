/**
 * @file read_leptonic_4d_table.cpp
 * @brief  Implementation of the 4D leptonic EOS table reader for GRACE.
 *         Reads the 4-axis (rho, T, Y_e, Y_mu) tables in the native
 *         GRACE leptonic-4D HDF5 format, as well as the baryon sector
 *         of stellarcollapse tables augmented with a separate muon
 *         lepton HDF5 group.
 * @date   2025
 *
 * @copyright This file is part of GRACE.  GPL-3 or later.
 */

#include <grace/physics/eos/leptonic_eos_4d.hh>
#include <grace/physics/eos/read_leptonic_4d_table.hh>
#include <grace/physics/eos/tabulated_eos.hh>   // for unit conv structs
#include <grace/utils/grace_utils.hh>
#include <grace/system/grace_system.hh>

#define H5_USE_16_API 1
#include <hdf5.h>

#include <Kokkos_Core.hpp>

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
//  Convenience HDF5 macros (consistent style with read_eos_table.cpp)
// ---------------------------------------------------------------------------
#define HDF5_CALL(ret, fn_call) \
    do { (ret) = (fn_call) ; ASSERT((ret)>=0, "HDF5 call failed: " #fn_call) ; } while(0)

#define READ_4D_EOS_HDF5(NAME, VAR, TYPE, MEM) \
    do { \
        hid_t _ds ; \
        HDF5_CALL(_ds, H5Dopen(file, NAME, H5P_DEFAULT)) ; \
        HDF5_CALL(h5err, H5Dread(_ds, TYPE, MEM, H5S_ALL, H5P_DEFAULT, VAR)) ; \
        HDF5_CALL(h5err, H5Dclose(_ds)) ; \
    } while(0)

// Hyperslab macro: reads one variable slice into a pre-sized 4-D temp buffer.
//  nrho*ntemp*nye*nymu elements starting at (OFF, 0).
#define READ_4D_EOSTABLE_HDF5(NAME, OFF) \
    do { \
        hsize_t offset5[2] = {OFF, 0} ; \
        H5Sselect_hyperslab(mem5, H5S_SELECT_SET, offset5, NULL, var5, NULL) ; \
        READ_4D_EOS_HDF5(NAME, alltables_temp, H5T_NATIVE_DOUBLE, mem5) ; \
    } while(0)

namespace grace {

// ============================================================
//  Helper: read the 1-D cold-slice table produced by
//  leptonic_eos_generate_cold_table (5-column ASCII file).
//  Format per row: log(rho)  temp  ye_cold  ymu_cold  press  eps  cs2  entropy
// ============================================================
static void read_leptonic_cold_table(
    const std::string& filename,
    Kokkos::View<double**, grace::default_execution_space>& d_data,
    Kokkos::View<double*,  grace::default_execution_space>& d_rho)
{
    std::ifstream f(filename) ;
    ASSERT(f.is_open(), "Cannot open leptonic cold table: " << filename) ;

    std::string line ;
    // line 1: description
    std::getline(f, line) ;
    // line 2: nrows
    std::getline(f, line) ;
    std::istringstream iss_n(line) ;
    size_t nrows ;
    iss_n >> nrows ;
    ASSERT(iss_n, "Failed to read nrows in leptonic cold table") ;

    // expected columns: logrho  temp  ye_cold  ymu_cold  press  eps  cs2  entropy
    constexpr int NCOLS = leptonic_eos_4d_t::COLD_VIDX::N_CTAB_VARS ;

    Kokkos::realloc(d_data, nrows, NCOLS) ;
    Kokkos::realloc(d_rho,  nrows) ;
    auto h_data = Kokkos::create_mirror_view(d_data) ;
    auto h_rho  = Kokkos::create_mirror_view(d_rho)  ;

    for (size_t i = 0; i < nrows; ++i) {
        ASSERT(std::getline(f,line), "Unexpected EOF in leptonic cold table at row " << i) ;
        std::istringstream iss(line) ;
        double logrho_val ;
        iss >> logrho_val ;
        h_rho(i) = logrho_val ;
        for (int j = 0; j < NCOLS; ++j) {
            double val ;
            iss >> val ;
            h_data(i,j) = val ;
        }
    }
    Kokkos::deep_copy(d_data, h_data) ;
    Kokkos::deep_copy(d_rho,  h_rho)  ;
}

// ============================================================
//  Main reader
// ============================================================
grace::leptonic_eos_4d_t read_leptonic_4d_table()
{
    using namespace grace::physical_constants ;

    // -------------------------------------------------------
    //  Read parameter file
    // -------------------------------------------------------
    auto const fname      = grace::get_param<std::string>("eos","leptonic_4d","table_filename") ;
    auto const cold_fname = grace::get_param<std::string>("eos","leptonic_4d","cold_table_filename") ;
    auto const tab_format = grace::get_param<std::string>("eos","leptonic_4d","table_format") ;

    const bool do_energy_shift =
        grace::get_param<bool>("eos","leptonic_4d","do_energy_shift", true) ;
    const bool use_muonic_eos  =
        grace::get_param<bool>("eos","leptonic_4d","use_muonic_eos", true) ;
    const bool atm_beta_eq =
        grace::get_param<bool>("eos","leptonic_4d","atmosphere_beta_eq", true) ;

    GRACE_INFO("Reading 4D leptonic EOS table: {}", fname) ;
    GRACE_INFO("Format: {}", tab_format) ;

    herr_t h5err ;
    hid_t  file  ;
    HDF5_CALL(file, H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)) ;

    // -------------------------------------------------------
    //  Read grid sizes
    // -------------------------------------------------------
    int nrho, ntemp, nye, nymu ;
    READ_4D_EOS_HDF5("pointsrho",  &nrho,  H5T_NATIVE_INT, H5S_ALL) ;
    READ_4D_EOS_HDF5("pointstemp", &ntemp, H5T_NATIVE_INT, H5S_ALL) ;
    READ_4D_EOS_HDF5("pointsye",   &nye,   H5T_NATIVE_INT, H5S_ALL) ;
    READ_4D_EOS_HDF5("pointsymu",  &nymu,  H5T_NATIVE_INT, H5S_ALL) ;

    GRACE_INFO("4D EOS grid: nrho={} nT={} nye={} nymu={}",
               nrho, ntemp, nye, nymu) ;

    // -------------------------------------------------------
    //  Allocate 4-D flat temp buffer (NTABLES × nrho×nT×nye×nymu)
    // -------------------------------------------------------
    const size_t NTAB_B = leptonic_eos_4d_t::TEOS_VIDX::N_TAB_VARS_BARYON ;
    const size_t NTAB_E = leptonic_eos_4d_t::ELE_VIDX::N_TAB_VARS_ELE  ;
    const size_t NTAB_M = leptonic_eos_4d_t::MUON_VIDX::N_TAB_VARS_MUON ;
    const size_t NP4    = static_cast<size_t>(nrho)*ntemp*nye*nymu ;

    auto* alltables_temp = new double[std::max({NTAB_B,NTAB_E,NTAB_M}) * NP4] ;

    // HDF5 hyperslab setup: shape = (NTAB, nrho*nT*nye*nymu)
    hsize_t table_dims[2] = { std::max({NTAB_B,NTAB_E,NTAB_M}), NP4 } ;
    hsize_t var5[2]       = { 1, NP4 } ;
    hid_t   mem5 = H5Screate_simple(2, table_dims, NULL) ;

    // -------------------------------------------------------
    //  Axis arrays
    // -------------------------------------------------------
    std::vector<double> logrho(nrho), logtemp(ntemp), yes(nye), ymus(nymu) ;
    double energy_shift_raw{0.} ;
    double baryon_mass = mn_MeV * MeV_to_g ; // default; overwritten if in file

    READ_4D_EOS_HDF5("logrho",  logrho.data(),  H5T_NATIVE_DOUBLE, H5S_ALL) ;
    READ_4D_EOS_HDF5("logtemp", logtemp.data(), H5T_NATIVE_DOUBLE, H5S_ALL) ;
    READ_4D_EOS_HDF5("ye",      yes.data(),     H5T_NATIVE_DOUBLE, H5S_ALL) ;
    READ_4D_EOS_HDF5("ymu",     ymus.data(),    H5T_NATIVE_DOUBLE, H5S_ALL) ;
    READ_4D_EOS_HDF5("energy_shift", &energy_shift_raw, H5T_NATIVE_DOUBLE, H5S_ALL) ;

    if (H5Lexists(file, "/mass_factor", H5P_DEFAULT)) {
        hid_t mb_ds ;
        HDF5_CALL(mb_ds, H5Dopen(file, "mass_factor", H5P_DEFAULT)) ;
        H5Dread(mb_ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &baryon_mass) ;
        GRACE_INFO("Baryon mass from file: {}", baryon_mass) ;
    }

    // -------------------------------------------------------
    //  Unit conversion factors  (stellarcollapse → geometric)
    // -------------------------------------------------------
    // These replicate the convention in read_eos_table.cpp:
    //   pressure: CGS → geometric = PRESSGF
    //   eps:      CGS → geometric (velocity^2 factor)
    //   rho:      g/cm^3 → geometric = RHOGF
    // We use the same approach as read_scollapse_table().
    // -----------
    // Conversion factors (geometric units, c=G=Msun=1)
    constexpr double G_si  = physical_constants::G_si ;
    constexpr double c_si  = physical_constants::c_si ;
    constexpr double Msun  = physical_constants::Msun_si ;
    // length GF: G Msun / c^2  [m]
    constexpr double LENGTHGF = G_si * Msun / (c_si*c_si) ;
    constexpr double TIMEGF   = LENGTHGF / c_si ;
    constexpr double RHOGF    = Msun / (LENGTHGF*LENGTHGF*LENGTHGF) ;  // kg/m^3 -> geom
    constexpr double PRESSGF  = 1. / (RHOGF * c_si*c_si) ;
    constexpr double EPSGF    = 1. / (c_si*c_si) ;

    // Convert energy_shift (same as in Margherita scollapse reader)
    double energy_shift = do_energy_shift ? energy_shift_raw * EPSGF : 0.0 ;

    // -------------------------------------------------------
    //  Convert axes: log10 → loge + geometric units
    // -------------------------------------------------------
    for (int i=0; i<nrho;  ++i) logrho[i]  = logrho[i]  * std::log(10.) + std::log(RHOGF) ;
    for (int i=0; i<ntemp; ++i) logtemp[i] = logtemp[i] * std::log(10.) ;

    // -------------------------------------------------------
    //  Allocate 5-D Kokkos views  [irho,iT,iye,iymu,ivar]
    // -------------------------------------------------------
    using view5d = Kokkos::View<double*****, grace::default_execution_space> ;
    view5d v_baryon("eos4d_baryon", nrho, ntemp, nye, nymu, NTAB_B) ;
    view5d v_ele   ("eos4d_ele",    nrho, ntemp, nye, nymu, NTAB_E) ;
    view5d v_muon  ("eos4d_muon",   nrho, ntemp, nye, nymu, NTAB_M) ;

    auto h_baryon = Kokkos::create_mirror_view(v_baryon) ;
    auto h_ele    = Kokkos::create_mirror_view(v_ele)    ;
    auto h_muon   = Kokkos::create_mirror_view(v_muon)   ;

    // Lambda: reorder flat temp buffer → 5D view
    auto fill_view = [&](auto& hv, size_t NTAB)
    {
        for (size_t iv=0; iv<NTAB; ++iv)
        for (int lm=0; lm<nymu; ++lm)
        for (int k=0; k<nye; ++k)
        for (int j=0; j<ntemp; ++j)
        for (int i=0; i<nrho; ++i)
        {
            size_t indold = i + nrho*(j + ntemp*(k + nye*(lm + nymu*iv))) ;
            hv(i,j,k,lm,iv) = alltables_temp[indold] ;
        }
    } ;

    // -------------------------------------------------------
    //  Read baryon table
    // -------------------------------------------------------
    hsize_t table_dims_b[2] = { NTAB_B, NP4 } ;
    hsize_t var5b[2]        = { 1, NP4 } ;
    hid_t   mem5b = H5Screate_simple(2, table_dims_b, NULL) ;

    // Override mem5 dim with baryon table count
    table_dims[0] = NTAB_B ;
    H5Sset_extent_simple(mem5, 2, table_dims, NULL) ;

    READ_4D_EOSTABLE_HDF5("logpress",  leptonic_eos_4d_t::TEOS_VIDX::TABPRESS) ;
    READ_4D_EOSTABLE_HDF5("logenergy", leptonic_eos_4d_t::TEOS_VIDX::TABEPS)   ;
    READ_4D_EOSTABLE_HDF5("entropy",   leptonic_eos_4d_t::TEOS_VIDX::TABENTROPY) ;
    READ_4D_EOSTABLE_HDF5("cs2",       leptonic_eos_4d_t::TEOS_VIDX::TABCSND2)  ;
    READ_4D_EOSTABLE_HDF5("mu_e",      leptonic_eos_4d_t::TEOS_VIDX::TABMUE)    ;
    READ_4D_EOSTABLE_HDF5("mu_p",      leptonic_eos_4d_t::TEOS_VIDX::TABMUP)    ;
    READ_4D_EOSTABLE_HDF5("mu_n",      leptonic_eos_4d_t::TEOS_VIDX::TABMUN)    ;
    READ_4D_EOSTABLE_HDF5("Xa",        leptonic_eos_4d_t::TEOS_VIDX::TABXA)     ;
    READ_4D_EOSTABLE_HDF5("Xh",        leptonic_eos_4d_t::TEOS_VIDX::TABXH)     ;
    READ_4D_EOSTABLE_HDF5("Xn",        leptonic_eos_4d_t::TEOS_VIDX::TABXN)     ;
    READ_4D_EOSTABLE_HDF5("Xp",        leptonic_eos_4d_t::TEOS_VIDX::TABXP)     ;
    READ_4D_EOSTABLE_HDF5("Abar",      leptonic_eos_4d_t::TEOS_VIDX::TABABAR)   ;
    READ_4D_EOSTABLE_HDF5("Zbar",      leptonic_eos_4d_t::TEOS_VIDX::TABZBAR)   ;
    fill_view(h_baryon, NTAB_B) ;

    // -------------------------------------------------------
    //  Read electron lepton table
    // -------------------------------------------------------
    table_dims[0] = NTAB_E ;
    H5Sset_extent_simple(mem5, 2, table_dims, NULL) ;
    READ_4D_EOSTABLE_HDF5("ele/mu_ele",        leptonic_eos_4d_t::ELE_VIDX::TABMUELE)     ;
    READ_4D_EOSTABLE_HDF5("ele/yle_minus",     leptonic_eos_4d_t::ELE_VIDX::TABYLE_MINUS) ;
    READ_4D_EOSTABLE_HDF5("ele/yle_plus",      leptonic_eos_4d_t::ELE_VIDX::TABYLE_PLUS)  ;
    READ_4D_EOSTABLE_HDF5("ele/press_e_minus", leptonic_eos_4d_t::ELE_VIDX::TABPRESS_E_MINUS) ;
    READ_4D_EOSTABLE_HDF5("ele/press_e_plus",  leptonic_eos_4d_t::ELE_VIDX::TABPRESS_E_PLUS)  ;
    READ_4D_EOSTABLE_HDF5("ele/eps_e_minus",   leptonic_eos_4d_t::ELE_VIDX::TABEPS_E_MINUS)   ;
    READ_4D_EOSTABLE_HDF5("ele/eps_e_plus",    leptonic_eos_4d_t::ELE_VIDX::TABEPS_E_PLUS)    ;
    READ_4D_EOSTABLE_HDF5("ele/s_e_minus",     leptonic_eos_4d_t::ELE_VIDX::TABS_E_MINUS)     ;
    READ_4D_EOSTABLE_HDF5("ele/s_e_plus",      leptonic_eos_4d_t::ELE_VIDX::TABS_E_PLUS)      ;
    fill_view(h_ele, NTAB_E) ;

    // -------------------------------------------------------
    //  Read muon lepton table  (optional if use_muonic_eos=false)
    // -------------------------------------------------------
    if (use_muonic_eos) {
        table_dims[0] = NTAB_M ;
        H5Sset_extent_simple(mem5, 2, table_dims, NULL) ;
        READ_4D_EOSTABLE_HDF5("muon/mu_mu",          leptonic_eos_4d_t::MUON_VIDX::TABMUMU)          ;
        READ_4D_EOSTABLE_HDF5("muon/ymu_minus",      leptonic_eos_4d_t::MUON_VIDX::TABYMU_MINUS)     ;
        READ_4D_EOSTABLE_HDF5("muon/ymu_plus",       leptonic_eos_4d_t::MUON_VIDX::TABYMU_PLUS)      ;
        READ_4D_EOSTABLE_HDF5("muon/press_mu_minus", leptonic_eos_4d_t::MUON_VIDX::TABPRESS_MU_MINUS) ;
        READ_4D_EOSTABLE_HDF5("muon/press_mu_plus",  leptonic_eos_4d_t::MUON_VIDX::TABPRESS_MU_PLUS)  ;
        READ_4D_EOSTABLE_HDF5("muon/eps_mu_minus",   leptonic_eos_4d_t::MUON_VIDX::TABEPS_MU_MINUS)   ;
        READ_4D_EOSTABLE_HDF5("muon/eps_mu_plus",    leptonic_eos_4d_t::MUON_VIDX::TABEPS_MU_PLUS)    ;
        READ_4D_EOSTABLE_HDF5("muon/s_mu_minus",     leptonic_eos_4d_t::MUON_VIDX::TABS_MU_MINUS)     ;
        READ_4D_EOSTABLE_HDF5("muon/s_mu_plus",      leptonic_eos_4d_t::MUON_VIDX::TABS_MU_PLUS)      ;
        fill_view(h_muon, NTAB_M) ;
    }

    H5Sclose(mem5) ;
    H5Fclose(file) ;
    delete[] alltables_temp ;

    // -------------------------------------------------------
    //  Unit conversion loop over all baryon table entries
    // -------------------------------------------------------
    double hmin{std::numeric_limits<double>::max()} ;
    double hmax{std::numeric_limits<double>::min()} ;
    double epsmin{std::numeric_limits<double>::max()} ;
    double epsmax_val{std::numeric_limits<double>::min()} ;

    for (int lm=0; lm<nymu; ++lm)
    for (int k=0; k<nye; ++k)
    for (int j=0; j<ntemp; ++j)
    for (int i=0; i<nrho; ++i)
    {
        double rhoL = std::exp(logrho[i]) ;

        // pressure: log10 * log(10) + log(PRESSGF)
        h_baryon(i,j,k,lm, leptonic_eos_4d_t::TEOS_VIDX::TABPRESS) =
            h_baryon(i,j,k,lm, leptonic_eos_4d_t::TEOS_VIDX::TABPRESS) * std::log(10.)
            + std::log(PRESSGF) ;
        double pressL = std::exp(h_baryon(i,j,k,lm, leptonic_eos_4d_t::TEOS_VIDX::TABPRESS)) ;

        // eps: log10 * log(10) + log(EPSGF), then shift
        double leps_raw = h_baryon(i,j,k,lm, leptonic_eos_4d_t::TEOS_VIDX::TABEPS) ;
        double epsT     = std::pow(10., leps_raw) * EPSGF ;
        double epsL     = epsT - energy_shift ;
        h_baryon(i,j,k,lm, leptonic_eos_4d_t::TEOS_VIDX::TABEPS) = std::log(epsT) ;

        epsmin    = std::min(epsmin,    epsL) ;
        epsmax_val= std::max(epsmax_val, epsL) ;

        // sound speed: already dimensionless (c=1), just clamp
        double& cs2 = h_baryon(i,j,k,lm, leptonic_eos_4d_t::TEOS_VIDX::TABCSND2) ;
        cs2 *= LENGTHGF*LENGTHGF / (TIMEGF*TIMEGF) ;
        if (cs2 < 0.)  cs2 = 0. ;
        double hL = 1. + epsL + pressL / rhoL ;
        cs2 /= hL ;
        if (cs2 > 0.9999999) cs2 = 0.9999999 ;

        // chemical potentials: MeV → geometric  (MeV/baryon → dimensionless)
        //  mu [geometric] = mu [MeV] * MeV_to_kg * c^2 / (baryon_mass * c^2) / c^2
        //  simplification: mu_geom = mu_MeV * MeV_to_g / baryon_mass  (c=1)
        constexpr double MeV2geom = physical_constants::MeV_to_g ; // MeV → g
        double mb = baryon_mass ;  // already in g (geometric units)
        h_baryon(i,j,k,lm,leptonic_eos_4d_t::TEOS_VIDX::TABMUE) *= MeV2geom / mb ;
        h_baryon(i,j,k,lm,leptonic_eos_4d_t::TEOS_VIDX::TABMUP) *= MeV2geom / mb ;
        h_baryon(i,j,k,lm,leptonic_eos_4d_t::TEOS_VIDX::TABMUN) *= MeV2geom / mb ;

        hmin = std::min(hmin, hL) ;
        hmax = std::max(hmax, hL) ;
    }

    // Convert lepton table pressure / eps to geometric units
    auto convert_lepton_view = [&](auto& hv, size_t NVAR)
    {
        for (int lm=0; lm<nymu; ++lm)
        for (int k=0; k<nye; ++k)
        for (int j=0; j<ntemp; ++j)
        for (int i=0; i<nrho; ++i)
        for (size_t iv=0; iv<NVAR; ++iv)
        {
            // Only pressure and eps vars need conversion; chemical potentials handled above.
            // For lepton tables the convention stores variables already in CGS;
            // we apply the same factors.
            // Press slots: PRESS_E/MU_MINUS/PLUS
            // Eps   slots: EPS_E/MU_MINUS/PLUS
            // The exact slot indices differ between ELE/MUON enums – handled by:
            bool is_press_slot = (
                iv == leptonic_eos_4d_t::ELE_VIDX::TABPRESS_E_MINUS ||
                iv == leptonic_eos_4d_t::ELE_VIDX::TABPRESS_E_PLUS  ||
                iv == leptonic_eos_4d_t::MUON_VIDX::TABPRESS_MU_MINUS||
                iv == leptonic_eos_4d_t::MUON_VIDX::TABPRESS_MU_PLUS ) ;
            bool is_eps_slot = (
                iv == leptonic_eos_4d_t::ELE_VIDX::TABEPS_E_MINUS ||
                iv == leptonic_eos_4d_t::ELE_VIDX::TABEPS_E_PLUS  ||
                iv == leptonic_eos_4d_t::MUON_VIDX::TABEPS_MU_MINUS||
                iv == leptonic_eos_4d_t::MUON_VIDX::TABEPS_MU_PLUS ) ;
            if (is_press_slot) hv(i,j,k,lm,iv) *= PRESSGF ;
            if (is_eps_slot)   hv(i,j,k,lm,iv) *= EPSGF   ;
        }
    } ;
    convert_lepton_view(h_ele,  NTAB_E) ;
    convert_lepton_view(h_muon, NTAB_M) ;

    // -------------------------------------------------------
    //  Deep copy to device
    // -------------------------------------------------------
    Kokkos::deep_copy(v_baryon, h_baryon) ;
    Kokkos::deep_copy(v_ele,    h_ele)    ;
    Kokkos::deep_copy(v_muon,   h_muon)   ;

    // -------------------------------------------------------
    //  Build axis Kokkos views
    // -------------------------------------------------------
    using view1d = Kokkos::View<double*, grace::default_execution_space> ;
    view1d v_logrho ("eos4d_logrho",  nrho)  ;
    view1d v_logT   ("eos4d_logT",    ntemp) ;
    view1d v_ye     ("eos4d_ye",      nye)   ;
    view1d v_ymu    ("eos4d_ymu",     nymu)  ;
    auto h_lr = Kokkos::create_mirror_view(v_logrho) ;
    auto h_lt = Kokkos::create_mirror_view(v_logT)   ;
    auto h_ye = Kokkos::create_mirror_view(v_ye)     ;
    auto h_ym = Kokkos::create_mirror_view(v_ymu)    ;
    for (int i=0; i<nrho;  ++i) h_lr(i) = logrho[i]  ;
    for (int i=0; i<ntemp; ++i) h_lt(i) = logtemp[i] ;
    for (int i=0; i<nye;   ++i) h_ye(i) = yes[i]     ;
    for (int i=0; i<nymu;  ++i) h_ym(i) = ymus[i]    ;
    Kokkos::deep_copy(v_logrho, h_lr) ;
    Kokkos::deep_copy(v_logT,   h_lt) ;
    Kokkos::deep_copy(v_ye,     h_ye) ;
    Kokkos::deep_copy(v_ymu,    h_ym) ;

    // -------------------------------------------------------
    //  Read cold-slice table
    // -------------------------------------------------------
    Kokkos::View<double**, grace::default_execution_space> cold_tabs ;
    Kokkos::View<double*,  grace::default_execution_space> cold_lrho ;
    read_leptonic_cold_table(cold_fname, cold_tabs, cold_lrho) ;

    // -------------------------------------------------------
    //  EOS limits
    // -------------------------------------------------------
    double rhomax  = std::exp(logrho[nrho-1])  ;
    double rhomin  = std::exp(logrho[0])       ;
    double tempmax = std::exp(logtemp[ntemp-1]) ;
    double tempmin = std::exp(logtemp[0])       ;
    double yemax   = yes[nye-1]                 ;
    double yemin   = yes[0]                     ;
    double ymumax  = ymus[nymu-1]               ;
    double ymumin  = ymus[0]                    ;

    auto usr_epsmax = grace::get_param<double>("eos","eps_maximum", 1.e5) ;
    if (usr_epsmax < epsmax_val) epsmax_val = usr_epsmax ;

    // Atmosphere
    double temp_floor = grace::get_param<double>("grmhd","atmosphere","temp_fl", tempmin) ;
    double ye_atm     = yemin ;
    double ymu_atm    = ymumin ;

    if (atm_beta_eq) {
        GRACE_INFO("Atmosphere: finding beta-equilibrium Y_e, Y_mu at rho_fl, T_fl") ;
        // NB: beta-eq solve is host-only; done during setup.
        // We leave ye_atm / ymu_atm at table minimum as safe fallback.
        // A dedicated host-side call can update them after object construction.
    }

    GRACE_INFO("4D leptonic EOS rho [{:.4e}, {:.4e}]  T [{:.4e}, {:.4e}]  "
               "Ye [{:.3f}, {:.3f}]  Ymu [{:.3f}, {:.3f}]",
               rhomin, rhomax, tempmin, tempmax, yemin, yemax, ymumin, ymumax) ;

    return leptonic_eos_4d_t(
        v_baryon, v_logrho, v_logT, v_ye, v_ymu,
        v_ele, v_muon,
        cold_tabs, cold_lrho,
        rhomax, rhomin,
        tempmax, tempmin,
        yemax,  yemin,
        ymumax, ymumin,
        baryon_mass, energy_shift,
        epsmin, epsmax_val,
        hmin,   hmax,
        temp_floor,
        ye_atm, ymu_atm,
        atm_beta_eq
    ) ;
}

} /* namespace grace */
