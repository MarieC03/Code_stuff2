#ifndef HH_INT_POW
#define HH_INT_POW
template<int N>
inline double int_pow(double const& x) {
  return x * int_pow<N-1>(x) ;
}
template<> inline double int_pow<0>(double const& x) { return 1. ; }
#endif
