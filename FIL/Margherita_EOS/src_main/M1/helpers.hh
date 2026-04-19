/*
 * =====================================================================================
 *
 *       Filename:  helpers.hh
 *
 *    Description:  Useful utility functions
 *
 *        Version:  1.0
 *        Created:  13/05/2017 19:40:44
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#ifndef __HH_M1_HELPERS
#define __HH_M1_HELPERS
namespace m1_helpers {

template <typename T>
constexpr inline T horner(const T x, const T A) {
  return A;
}

template <typename T, typename... TArgs>
constexpr inline T horner(const T x, const T A, const T B,
                          const TArgs... Args) {
  return horner(x, B, Args...) * x + A;
}
}

#endif
