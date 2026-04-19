#include "../src/utils/metric.hh"


#include "../src/frankfurt_m1.h"

#include <iostream>
#include <iomanip>


using namespace std ;

int main() {


  std::array<double,10> metric ;

  metric[GXX] = 1.1;
  metric[GXY] = 0. ;
  metric[GXZ] = 0.01 ;
  metric[GYY] = 1.1 ;
  metric[GYZ] = 1.01 ;
  metric[GZZ] = 1.1 ;

  metric[SHIFTX] = metric[SHIFTY] = metric[SHIFTZ] = 0.;
  metric[LAPSE]  = 1. ;

  metric_c gamma( {metric[GXX],
		   metric[GXY],
		   metric[GXZ],
		   metric[GYY],
		   metric[GYZ],
		   metric[GZZ]}, {metric[SHIFTX],metric[SHIFTY],metric[SHIFTZ]}, metric[LAPSE]) ;


  std::array<double,9> g,gup ;
  
  cout << "Metric and inverse: \n";
  int idx = 0;
  for( int ii=0; ii<3; ii++)
    for( int jj=ii; jj<3; jj++){
    cout << gamma.metric_comp(idx) << ", " << gamma.invmetric_comp(idx) << std::endl ;
    gup[jj+3*ii] = gup[ii+3*jj] = gamma.invmetric_comp(idx) ;
    g[jj+3*ii] = g[ii+3*jj] = gamma.metric_comp(idx) ;
    idx++ ;
  }

  
  cout << "Detg: " << gamma.sqrtdet << endl ;

  std::array<double,3> vector_up{0.5,0.7,1.2} ;
  auto const vector_low = gamma.lower_index<0,3>(vector_up) ;
  auto sanity_check = gamma.raise_index<0,3>(vector_low) ;

  cout << "Check inverse metric: \n";

  std::array<double,9> delta{0} ;
  for( int ii=0; ii<3; ii++){
    for( int jj=0; jj<3; jj++){
      for( int kk=0; kk<3; kk++) {
	delta[jj+3*ii] += g[3*kk+jj]*gup[kk+3*ii];
      }
      cout << delta[jj+3*ii] << '\t' ;
    }
    cout << '\n';
  }


  

  for( int ii=0; ii<3; ii++)
    sanity_check[ii] -= vector_up[ii] ;
 
  cout << fixed << setprecision(16) ;

  cout << "Sanity check, everything should be 0\n";
  for( int ii=0; ii<3; ii++)
    cout << sanity_check[ii] << '\t';
  cout << "\n\n";



  cout << "Checking scalar products: \n";

  auto const v2_metric = gamma.scalar_product(vector_low,vector_up) ;
  auto const v2_metric_l = gamma.square_3_lower<0,3>(vector_low) ;
  auto const v2_metric_u = gamma.square_3_upper<0,3>(vector_up) ;

  auto const v2_metric_ll = gamma.scalar_product_ll(vector_low,vector_low) ;
  auto const v2_metric_uu = gamma.scalar_product_uu(vector_up,vector_up) ;


  cout << "Everything should be the same: \n"
       << v2_metric << ", " << v2_metric_l << ", " << v2_metric_u << ", " << v2_metric_ll << ", " << v2_metric_uu << "\n\n" ;

  
		      
  
  
}
