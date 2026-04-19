#include <stdlib.h>
#include <iostream>
#include <array>
#include<cmath>

template<size_t N>
std::array<double,N> generate_linspace(double const &min,double const& max){
  std::array<double,N> out;
  for(int i=0;i<N;i++){
    out[i] = min + i * (max-min)/N;
  }
  return out;
}

template<size_t N>
std::array<double,N> generate_logspace(double const &min,double const& max){
  std::array<double,N> out;
  double lmin{log10(min)}, lmax{log10(max)};
  for(int i=0;i<N;i++){
    out[i] = lmin + i * (lmax-lmin)/N;
  }
  return out;
}

constexpr int num_rho = 10;
constexpr int num_ye = 60;
constexpr int num_temp = 10;

double const rhomin{1e-10}, rhomax{5e-03};
double const tempmin{1e-02}, tempmax{150};
double const yemin{0.01}, yemax{0.6};

int main() {

  std::array<double,num_rho> lrho = generate_logspace<num_rho>(rhomin,rhomax);
  std::array<double,num_ye> ye = generate_linspace<num_ye>(yemin,yemax);
  std::array<double,num_temp> ltemp = generate_logspace<num_temp>(tempmin,tempmax);

  std::cout << "rho: " << std::endl;
  for(auto const& r: lrho) std::cout << r << "\t";
  std::cout << std::endl ;

  std::cout << "T: " << std::endl;
  for(auto const& r: ltemp) std::cout << r << "\t";
  std::cout << std::endl ;

  std::cout << "Y_e: " << std::endl;
  for(auto const& r: ye) std::cout << r << "\t";
  std::cout << std::endl ;
}
