////
//// This file is part of Margherita, the light-weight EOS framework
////
////  Copyright (C) 2017, Elias Roland Most
////                      <emost@th.physik.uni-frankfurt.de>
////
////  This program is free software: you can redistribute it and/or modify
////  it under the terms of the GNU General Public License as published by
////  the Free Software Foundation, either version 3 of the License, or
////  (at your option) any later version.
////
////  This program is distributed in the hope that it will be useful,
////  but WITHOUT ANY WARRANTY; without even the implied warranty of
////  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
////  GNU General Public License for more details.
////
////  You should have received a copy of the GNU General Public License
////  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef _HH_METRIC
#define _HH_METRIC

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>


class metric_c {
 public:
  typedef std::array<double, 6> metric_ptr;
  typedef std::array<double, 3> vec_ptr;

 private:
  metric_ptr invmetric, metric ;
  vec_ptr shift ;

  static constexpr int WVX = 0;
  static constexpr int WVY = 1;
  static constexpr int WVZ = 2;

  static constexpr int GXX = 0;
  static constexpr int GXY = 1;
  static constexpr int GXZ = 2;
  static constexpr int GYY = 3;
  static constexpr int GYZ = 4;
  static constexpr int GZZ = 5;

  static inline double SQ(const double &x) { return x * x; };

 public:
  //! Square root of the determinant
  /*!
    Store the square root of the 3-metric \f$\sqrt{\gamma}\f$
  */
  double sqrtdet;
  double lapse = 1.;

  //! Compute determinant
  /*!
    Compute the determinant of the three-metric \f$\gamma_{ij}\f$
  */
  static inline double det(const metric_ptr &__metric) {
    //   using namespace Margherita_C2P_conventions;
    return -SQ(__metric[GXZ]) * __metric[GYY] +
           2.0 * __metric[GXY] * __metric[GXZ] * __metric[GYZ] -
           __metric[GXX] * SQ(__metric[GYZ]) - SQ(__metric[GXY]) * __metric[GZZ] +
           __metric[GXX] * __metric[GYY] * __metric[GZZ];
  }

  static inline metric_ptr compute_inv_metric(const metric_ptr &__metric,
                                              const double &__sqrtdet) {
    metric_ptr inv;
    const double detL = SQ(__sqrtdet);
    inv[GXX] = (-SQ(__metric[GYZ]) + __metric[GYY] * __metric[GZZ]) / detL;
    inv[GXY] = ((__metric[GYZ] * __metric[GXZ]) - __metric[GXY] * __metric[GZZ]) / detL;
    inv[GYY] = (-SQ(__metric[GXZ]) + __metric[GXX] * __metric[GZZ]) / detL;
    inv[GXZ] =
        (-(__metric[GXZ] * __metric[GYY]) + __metric[GXY] * __metric[GYZ]) / detL;
    inv[GYZ] = ((__metric[GXY] * __metric[GXZ]) - __metric[GXX] * __metric[GYZ]) / detL;
    inv[GZZ] = (-SQ(__metric[GXY]) + __metric[GXX] * __metric[GYY]) / detL;

    return inv;
  }

  inline metric_c(metric_ptr &metric, metric_ptr &invmetric,
                  const double sqrtdet, vec_ptr &shift, double const lapse_)
      : metric(std::move(metric)),
        invmetric(std::move(invmetric)), shift(std::move(shift)), lapse(lapse_),
        sqrtdet(fabs(sqrtdet)) {}  //!< a constructor
					       
  inline metric_c(metric_ptr &__metric, metric_ptr &invmetric, vec_ptr &shift, double const lapse_)
    : metric(std::move(__metric)), invmetric(std::move(invmetric)),
      shift(std::move(shift)), lapse(lapse_){
    sqrtdet = sqrt(fabs(det(metric)));
  }

  inline metric_c(metric_ptr &&__metric, vec_ptr && __shift, double const lapse_) : metric(std::move(__metric)), shift(std::move(__shift)), lapse(lapse_) {
    sqrtdet = sqrt(std::abs(det(metric)));
    invmetric = compute_inv_metric(metric, sqrtdet);
  }

  inline metric_c(metric_ptr const &__bssn_metric, double const phi, vec_ptr& shift, double const lapse_) :
    shift(std::move(shift)), lapse(lapse_) {
    double const psi4 = exp(4.*phi) ;
    sqrtdet = exp(6.*phi) ;
    for( int ii=0; ii<6; ii++)
      metric[ii] = __bssn_metric[ii]*psi4 ;

    invmetric = compute_inv_metric(metric,sqrtdet) ;
  }

  //! Raise index function
  /*
    this is a description of the raising function.
  */
  template <int offset, int num, typename T>
  inline auto raise_index(const std::array<T, num> &a) const -> std::array<T,3> {
    //    using namespace Margherita_C2P_conventions;
    constexpr int X = 0;
    constexpr int Y = 1;
    constexpr int Z = 2;

    std::array<T,3> a_up{
        {a[offset + X] * invmetric[GXX] + a[offset + Y] * invmetric[GXY] +
             a[offset + Z] * invmetric[GXZ],
         a[offset + X] * invmetric[GXY] + a[offset + Y] * invmetric[GYY] +
             a[offset + Z] * invmetric[GYZ],
         a[offset + X] * invmetric[GXZ] + a[offset + Y] * invmetric[GYZ] +
             a[offset + Z] * invmetric[GZZ]}};
    return a_up;
  }

  template <int offset, int num, typename T>
  inline auto lower_index(const std::array<T, num> &a) const -> std::array<T,3>{
    //    using namespace Margherita_C2P_conventions;
    constexpr int X = 0;
    constexpr int Y = 1;
    constexpr int Z = 2;
    std::array<T,3> a_up{{a[offset + X] * metric[GXX] + a[offset + Y] * metric[GXY] +
                      a[offset + Z] * metric[GXZ],
                  a[offset + X] * metric[GXY] + a[offset + Y] * metric[GYY] +
                      a[offset + Z] * metric[GYZ],
                  a[offset + X] * metric[GXZ] + a[offset + Y] * metric[GYZ] +
                      a[offset + Z] * metric[GZZ]}};
    return a_up;
  }

  template <int offset, int num,typename T>
  inline auto square_3_upper(const std::array<T, num> &E) const -> T const{
    //    using namespace Margherita_C2P_conventions;
    return metric[GXX] * E[offset + WVX] * E[offset + WVX] +
           metric[GYY] * E[offset + WVY] * E[offset + WVY] +
           metric[GZZ] * E[offset + WVZ] * E[offset + WVZ] +
           2.0 * metric[GXY] * E[offset + WVX] * E[offset + WVY] +
           2.0 * metric[GXZ] * E[offset + WVX] * E[offset + WVZ] +
           2.0 * metric[GYZ] * E[offset + WVY] * E[offset + WVZ];
  }

  template <int offset, int num, typename T>
  inline auto square_3_lower(
      const std::array<T, num> &E) const -> T const {
    // using namespace Margherita_C2P_conventions;
    return invmetric[GXX] * E[offset + WVX] * E[offset + WVX] +
           invmetric[GYY] * E[offset + WVY] * E[offset + WVY] +
           invmetric[GZZ] * E[offset + WVZ] * E[offset + WVZ] +
           2.0 * invmetric[GXY] * E[offset + WVX] * E[offset + WVY] +
           2.0 * invmetric[GXZ] * E[offset + WVX] * E[offset + WVZ] +
           2.0 * invmetric[GYZ] * E[offset + WVY] * E[offset + WVZ];
  }

  //! The scalar product of two vectors
  /*!
    This function computes the scalar product of two vectors \f$a\cdot b\f$
    \param a vector 1
    \param b vector 2
  */
  template<typename Ta, typename Tb>
  inline static auto scalar_product(const Ta &a, const Tb &b) -> decltype(a[0]*b[0]) {
    decltype(a[0]*b[0]) sp = 0;
#pragma unroll
    for (int i = 0; i < 3; ++i) sp += a[i] * b[i];
    return sp;
  }

  template<typename Ta, typename Tb>
  inline auto scalar_product_ll(const Ta &a, const Tb &b) -> decltype(a[0]*b[0]) {
    auto tmp = raise_index<0,3>(a) ;
    return scalar_product(tmp,b) ;
  }

  template<typename Ta, typename Tb>
  inline auto scalar_product_uu(const Ta &a, const Tb &b) -> decltype(a[0]*b[0]) {
    auto tmp = lower_index<0,3>(a);
    return scalar_product(tmp,b);
  }

  inline static auto contract_two_tensor_symm( const std::array<double,6> & A, const std::array<double,6> & B ) {
    constexpr int XX=0, XY=1, XZ=2, YY=3, YZ=4, ZZ=5 ;
    return A[XX]*B[XX] + A[YY]*B[YY] + A[ZZ]*B[ZZ]
      + 2. * ( A[XY] * B[XY] + A[XZ]*B[XZ] + A[YZ]*B[YZ] ) ;
  }

template<bool low, size_t offseta, size_t offsetb, size_t numa,size_t numb, typename Ta, typename Tb>
inline auto vector_Product(
    const std::array<Ta, numa> &a, std::array<Tb, numb> const &b) const -> std::array<decltype(a[0]*b[0]),3> const {

    auto const gammafac = low ? sqrtdet : 1./sqrtdet;

    std::array<decltype(a[0]*b[0]),3> c {
    (a[1+offseta]*b[2+offsetb]-a[2+offseta]*b[1+offsetb])*gammafac,
    (a[2+offseta]*b[0+offsetb]-a[0+offseta]*b[2+offsetb])*gammafac,
    (a[0+offseta]*b[1+offsetb]-a[1+offseta]*b[0+offsetb])*gammafac};

    return c;

  };

  inline __attribute__((always_inline)) auto metric_comp(size_t comp) const {
    assert( comp < 6 ) ;
    return metric[comp] ;
  }

  inline __attribute__((always_inline)) auto invmetric_comp(size_t comp) const {
    assert( comp < 6 ) ;
    return invmetric[comp] ;
  }
  
  inline __attribute__((always_inline)) auto get_shift(size_t dirn) const {
    assert( dirn < 3 ) ;
    return shift[dirn] ;
  }

  inline __attribute__((always_inline)) auto get_lapse() const {
    return lapse ;
  } 

};
#endif
