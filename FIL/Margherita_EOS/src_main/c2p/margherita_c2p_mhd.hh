#include "margherita_c2p_mhd_formulations.hh"
#include "margherita_c2p_kastaun_f.hh"
#include "margherita_c2p_mhd_t.hh"

// Instantiate
template <typename eos>
using C2P_MHD = C2P_MHD_t<eos, C2P_MHD_Palenzuela_f>;

template <typename eos>
using C2P_MHD_ENTROPY = C2P_MHD_t<eos, C2P_MHD_Entropy_Fix_f>;

template <typename eos>
using C2P_MHD_ENTROPY_NEW = C2P_MHD_t<eos, C2P_MHD_Entropy_Fix_New_f>;

template <typename eos>
using C2P_MHD_NEWMAN = C2P_MHD_t<eos, C2P_MHD_Newman_f>;

template <typename eos>
using C2P_MHD_KASTAUN = C2P_MHD_t<eos,C2P_MHD_kastaun_f>;
