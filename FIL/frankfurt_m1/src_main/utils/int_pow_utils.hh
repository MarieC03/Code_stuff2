#ifndef M1_INT_POW__HH_
#define M1_INT_POW__HH_

template<int N>
static inline __attribute__((always_inline)) double int_pow(double const& x) {
  return x * int_pow<N-1>(x) ;
}

template<>
inline __attribute__((always_inline)) double int_pow<1>(double const& x) { return x ; }

static inline constexpr __attribute__((always_inline)) int kronecker_delta(int i, int j) {
  return i==j ;
}

#endif
