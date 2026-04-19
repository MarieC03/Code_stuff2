#include <benchmark/benchmark.h>
#include <cmath>

#include "../src/utils/brent.hh"
#include "../src/utils/rootfind.hh"

static void BM_brent_orig(benchmark::State& state)
{

  auto const func_to_rootfind = [&] (double const & x)
  {
    return log(x) ;
  } ;
  double root ;
  double const tol = 1e-15 ;
  for(auto _: state )
    {
      root = zero_brent(1e-12, 10, tol, func_to_rootfind) ;
      benchmark::DoNotOptimize(root) ;
    }
  benchmark::ClobberMemory() ;
  
}


static void BM_new_brent(benchmark::State& state)
{

  auto const func_to_rootfind = [&] (double const & x)
  {
    return log(x) ;
  } ;
  std::pair<double,double> root ;
  double const tol = 1e-15 ;
  unsigned long iter {10000} ; 
  for(auto _: state )
    {
      root = rootfind::rootfind_brent(1e-12, 10.0, func_to_rootfind, tol, iter) ;
      benchmark::DoNotOptimize(root) ;
    }
  benchmark::ClobberMemory() ;
  
}


BENCHMARK(BM_brent_orig) ;
BENCHMARK(BM_new_brent)  ;

BENCHMARK_MAIN()         ; 
