#ifndef ROOTFIND_UTILS_HH_
#define ROOTFIND_UTILS_HH_


namespace utils {

  template<typename T>
  inline __attribute__((always_inline)) int sign(T const& x) {
    return static_cast<int>( ( x>static_cast<T>(0.) ) - ( x < static_cast<T>(0.) ) ) ;
  }

  template<typename T>
  inline __attribute__((always_inline)) bool contains( T const&a, T const& b, T const& x) {
    return static_cast<bool>( (x>a)*(x<b) ) ;
  }
}

#endif
