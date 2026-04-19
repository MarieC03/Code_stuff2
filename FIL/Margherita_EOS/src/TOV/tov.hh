/*
 * =====================================================================================
 *
 *       Filename:  tov.hh
 *
 *    Description:  Margherita TOV: Simple TOV solving routines
 *
 *        Version:  1.0
 *        Created:  27/12/2017 17:33:46
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#include<array>
#include<cmath>
#include<memory>
#include<vector>
#include<iostream>

#include "eos_wrapper.hh"

template<typename EOS>
class MargheritaTOV {
  private:
    enum {RHOB=0,PHI,PRESS,MASSR,MASSB,TIDALH,TIDALBETA,RADIUS,RADIUS_ISOTROPIC,NUM_STATE_VARS};
    using state_array_t = std::array<std::vector<double>, NUM_STATE_VARS>;
    using state_t = std::array<double, NUM_STATE_VARS>;


    size_t iteration=0;

    size_t max_iterations = 40000;

    double radiusstep0 = 20./static_cast<double>(max_iterations);
    double radiusstep = radiusstep0;

    bool finished = false;

    bool init_tidal = false;



  public:

    state_array_t state;

    double mass=0;
    double radius=0;
    double baryon_mass=0;
    double tidal_love_k2=0;

    double press_c;
    
    bool compact_output = false;
    bool adaptive_step = true;

    void reset(){
	mass = 0;
	radius = 0;
	baryon_mass =0;
	tidal_love_k2=0;
	iteration=0;
	init_tidal = false;
	finished = false;

        radiusstep = radiusstep0;

	auto error = typename EOS::error_t();
	state[RHOB][0] = EOS::rho__press_cold(press_c,error);
	state[PRESS][0] = press_c;
	state[PHI][0] = 0; //Can really be anything here
	state[MASSR][0] = 0;
	state[MASSB][0] = 0;
	state[TIDALH][0] = 0;
	state[TIDALBETA][0] = 0;
	state[RADIUS][0] = 0;
	state[RADIUS_ISOTROPIC][0] = 0;
    };

    MargheritaTOV(double press_c) : press_c(press_c){
     
        for(int nn=0; nn<NUM_STATE_VARS; ++nn){
	  state[nn].reserve(max_iterations);
	};
    };


    state_t const compute_rhs(state_t const &stateL) {

      auto rhs = state_t {0, 0, 0, 0, 0, 0, 0, 0};

      if(stateL[RADIUS]==0) return rhs;

	auto error = typename EOS::error_t();
	double energy,dedp;
	auto press = stateL[PRESS];
	auto rho = EOS::rho_energy_dedp__press_cold(energy,dedp,press,error);	

	auto const rhoh = energy + press;

	auto const r = stateL[RADIUS];
	auto const r2 = r*r;
	auto const r3 = r2*r;

	auto const mr = stateL[MASSR];

	rhs[PRESS] = - rhoh * (mr + 4.*M_PI * r3 * press ) / (r*(r-2.*mr));
	rhs[PHI] = - rhs[PRESS]/rhoh;
	rhs[MASSR] = 4.*M_PI*r2*energy;

	auto H = stateL[TIDALH];
	auto beta = stateL[TIDALBETA];

	if(iteration==0){
	  H = r2;
	  beta = 2.*r;

	  init_tidal = true;
	}

	rhs[TIDALH] = beta;

	auto schwarzschild_fac = r/(r-2*mr);
	rhs[MASSB]  = rhs[MASSR]*sqrt(schwarzschild_fac);

//	auto const dpdrho = EOS::dpress_cold_drho__rho(rho,error);

	auto const cs2 = 1./dedp; //dpdrho*rho/rhoh;
	
	auto const tmp1 = (4. * M_PI * r *press + mr/r2);

      //https://arxiv.org/pdf/0911.3535.pdf
	rhs[TIDALBETA] = 2.*H*schwarzschild_fac *
	  ( - 2. * M_PI * ( 5.* energy  + 9. * press + rhoh/cs2 )
	    + 3./r2 + 2.*tmp1*tmp1*schwarzschild_fac)
	  +2.*beta/r * schwarzschild_fac *
	  (-1. + mr/r + 2.*M_PI*r2*(energy - press));

	rhs[RADIUS_ISOTROPIC] = stateL[RADIUS_ISOTROPIC]/stateL[RADIUS] * sqrt(schwarzschild_fac);

	return rhs;
    };

    state_t get_state(){
      return state_t {
	state[RHOB][iteration],
	state[PHI][iteration],
	state[PRESS][iteration],
	state[MASSR][iteration],
	state[MASSB][iteration],
	state[TIDALH][iteration],
	state[TIDALBETA][iteration],
	state[RADIUS][iteration],
	state[RADIUS_ISOTROPIC][iteration]
      };
    }

    void insert_state(state_t const &stateL){
	for(int nn = 0; nn < NUM_STATE_VARS; ++nn){
	  state[nn][iteration+1] = stateL[nn];
	}
    };

    inline state_t add_states(state_t const &a, double const &lambda, state_t const &b){
      state_t result;
        auto const radius = a[RADIUS];
	for(int nn = 0; nn < NUM_STATE_VARS; ++nn){
	  result[nn] = a[nn] + lambda*b[nn];
	}
	result[RADIUS] = radius + lambda;
	return result;
    };


    void check_progress(){

      if(state[PRESS][iteration] <= 0 ){
	//Reset iteration and adjust radius step so as to go to zero!
	iteration--;

	auto stateL = get_state();

	auto rhs = compute_rhs(stateL);

	//Use simple Euler step here...
	auto const delta_r = -stateL[PRESS]/rhs[PRESS];

	auto const radius = stateL[RADIUS];

	state_t new_state;

	for(int nn = 0; nn<NUM_STATE_VARS; ++nn){
	  new_state[nn] = stateL[nn] + delta_r*rhs[nn];
	}

	new_state[RADIUS] = radius + delta_r;

	auto error = typename EOS::error_t();
	new_state[RHOB] = EOS::rho__press_cold(new_state[PRESS],error);

	insert_state(new_state);

	finished = true;
	
      }

    };

    state_t const evolve_state_rk45(){
      auto stateL = get_state();
      auto const radius = stateL[RADIUS];

      auto eps = 5.e-8*stateL[PRESS];
 
      state_t new_state;

      while(true){

	      auto k1 = compute_rhs(stateL);
	      state_t s1;
	      for(int nn = 0; nn<NUM_STATE_VARS; ++nn){
		 s1[nn] = stateL[nn] + 2./9.*radiusstep*k1[nn];
	      }
	      s1[RADIUS] = stateL[RADIUS] + 2./9.*radiusstep;

	      auto k2 = compute_rhs(s1);

	      state_t s2;
	      for(int nn = 0; nn<NUM_STATE_VARS; ++nn){
		 s2[nn] = stateL[nn] + radiusstep*(1./12.*k1[nn] + 0.25*k2[nn]);
	      }
	      s2[RADIUS] = stateL[RADIUS] + 1./3.*radiusstep;

	      auto k3 = compute_rhs(s2);

	      state_t s3;
	      for(int nn = 0; nn<NUM_STATE_VARS; ++nn){
		 s3[nn] = stateL[nn] + radiusstep*(69./128.*k1[nn] - 243./128.*k2[nn] + 135./64.*k3[nn]);
	      }
	      s3[RADIUS] = stateL[RADIUS] + 3./4.*radiusstep;

	      auto k4 = compute_rhs(s3);

	      state_t s4;
	      for(int nn = 0; nn<NUM_STATE_VARS; ++nn){
		 s4[nn] = stateL[nn] + radiusstep*(-17./12.*k1[nn] + 27./4.*k2[nn] - 27./5.*k3[nn] + 16./15.*k4[nn]);
	      }
	      s4[RADIUS] = stateL[RADIUS] + radiusstep;

	      auto k5 = compute_rhs(s4);

	      state_t s5;
	      for(int nn = 0; nn<NUM_STATE_VARS; ++nn){
		 s5[nn] = stateL[nn] + radiusstep*(65./432.*k1[nn] - 5./16.*k2[nn] + 13./16.*k3[nn] + 4./27.*k4[nn] + 5./144.*k5[nn]);
	      }
	      s5[RADIUS] = stateL[RADIUS] + 5./6.*radiusstep;

	      auto k6 = compute_rhs(s5);

	      for(int nn = 0; nn<NUM_STATE_VARS; ++nn){
		 new_state[nn] = stateL[nn] + radiusstep*(1./9.*k1[nn] + 9./20.*k3[nn] + 16./45.*k4[nn] + 1./12.*k5[nn]);
	      }
	      new_state[RADIUS] = stateL[RADIUS] + radiusstep;


	      //compute error
	      state_t error;
	      for(int nn = 0; nn<NUM_STATE_VARS; ++nn){
		 error[nn] = radiusstep*1./300.*(-2.*k1[nn] + 9.*k3[nn] - 64.*k4[nn] - 15.*k5[nn] + 72.*k6[nn]);
	      };

	      //Is the error too large?
	      if(std::abs(error[PRESS]) > eps){
		 if(radiusstep > 2.*radiusstep0){
		    radiusstep = radiusstep*0.5; 	     
		 }else{
			 break;
		 }

		 continue;
	      }
	 
	      if(std::abs(error[PRESS]) < 1./20.*eps) radiusstep*=2.;

	      break;
     };

     auto press = stateL[PRESS];

     auto error = typename EOS::error_t();
     new_state[RHOB] = EOS::rho__press_cold(press,error);

     return new_state;

    };


    state_t const evolve_state(){

      auto stateL = get_state();
      auto const radius = stateL[RADIUS];

      auto k1 = compute_rhs(stateL);

      auto k2 = compute_rhs(add_states(stateL,0.5*radiusstep,k1));

      auto k3 = compute_rhs(add_states(stateL,0.5*radiusstep,k2));
      
      auto k4 = compute_rhs(add_states(stateL,radiusstep,k3));

	//Simple RK4 update
	for(int nn = 0; nn<NUM_STATE_VARS; ++nn){
	  stateL[nn] +=  1./6.*radiusstep*(
	      k1[nn] + 2.*k2[nn] + 2.*k3[nn] + k4[nn]);
	}

	stateL[RADIUS] = radius + radiusstep;

	auto press = stateL[PRESS];

	auto error = typename EOS::error_t();
	stateL[RHOB] = EOS::rho__press_cold(press,error);

	return stateL;

    };

    void solve(){
      reset();

      while(!finished){
	auto new_state = state_t{};

	if(adaptive_step){
          new_state = evolve_state_rk45();
	}else{
	  new_state = evolve_state();
	}

	if(iteration==0){
	  auto const r= new_state[RADIUS];
	  new_state[TIDALH] = r*r;
	  new_state[TIDALBETA] = 2.*r;
	}

	insert_state(new_state);
	iteration++;
	check_progress();
/*  
	std::cout << std::endl;
	std::cout << "Iteration " << iteration << std::endl;
	std::cout << std::endl;
	for(int nn = 0; nn < NUM_STATE_VARS; ++nn){
	  std::cout << state[nn][iteration] <<std::endl;
	}
*/	
	if(iteration == max_iterations -1){
	   if(!compact_output)
	   std::cout << "Did not converge, aborting!" << std::endl;
	   return;
	};


      }

      mass = state[MASSR][iteration];
      baryon_mass = state[MASSB][iteration];
      radius = state[RADIUS][iteration];


      auto const C = mass/radius; //compactness
      auto const C2 = C*C;
      auto const C5 = C2*C2*C;
      auto const tmp2 = (1.-2.*C) * (1.-2.*C);

      //The way the equations are written, the metric should be
      //ds^2=...+r^2*dOmega^2. Then, r=r(P=0) is the circumferential
      //radius according to definition. Hence, what you finally call 
      //coord_radius should be the circumferential radius. Also according
      //to RNS this is the circumferential radius (I ran several tests to
      //check this).

      //Adjust PHI by matching to Schwarzschild
      auto const correction = 0.5*log(1.-2.* C) - state[PHI][iteration];

      for(auto & value : state[PHI])
	value += correction;

      //Now compute the tidal deformability k2


      auto const y = radius*state[TIDALBETA][iteration]/state[TIDALH][iteration];

      //https://arxiv.org/pdf/0911.3535.pdf
      tidal_love_k2 = 8./5.*C5*tmp2*(2. + 2. * C * (y-1.)  -y)/
	      (2.*C*(6.-3.*y + 3.*C*(5.*y-8)) 
	       + 4.*C2*C * (13. - 11.*y + C*(3.*y -2.) + 2.*C2*(1.+y))
	       + 3.*tmp2*(2.-y + 2.*C*(y-1.))*log(1.-2.*C));


    };


    inline size_t const get_iterations(){
	return iteration;
    };


    friend std::ostream& operator<< (std::ostream& stream, const MargheritaTOV& tov) {

      if(tov.compact_output){
      stream << std::setprecision(12);
      stream << tov.press_c << "\t"
      	     << tov.mass << "\t"
      	     << tov.baryon_mass << "\t"
      	     << tov.radius << "\t"
      	     << tov.tidal_love_k2 <<std::endl;
      }else{

      stream << " Central pressure: " << tov.press_c << std::endl;
      stream << " Mass: " << tov.mass << std::endl;
      stream << " Baryon mass: " << tov.baryon_mass << std::endl;
      stream << " Radius: " << tov.radius << std::endl;
      stream << " Love number k2: " << tov.tidal_love_k2 << std::endl;
      };

      return stream;
    };

    
};

