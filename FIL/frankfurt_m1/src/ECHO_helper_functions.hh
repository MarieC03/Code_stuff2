#ifndef _HH_HELPER_FUNCTIONS
#define _HH_HELPER_FUNCTIONS

#include "constexpr_utils.hh"

/*
 *
 *   Some helper function to do higher order flux corrections
 *   as in ECHO.
 *
 *   Written by Elias R. Most in Feb. 2017
 *
 */

template<int order>
inline constexpr CCTK_REAL InterpToFaces_LowLevel(const CCTK_REAL *F){
                return (order==0)*0.5*(F[0]+F[1])
                      +(order==1)*(-1./16.*(F[0]+F[3])+9./16.*(F[1]+F[2]))
                      +(order==2)*(3./256.*(F[0]+F[5])-25./256.*(F[1]+F[4])+150./256.*(F[2]+F[3]));
}

template<int dir, int order>
inline static CCTK_REAL InterpToFaces(CCTK_POINTER_TO_CONST _GH, const CCTK_REAL *F, int index[3]){
    cGH const * const cctkGH = (cGH const * const)_GH;
	CCTK_REAL B[2*(order+1)];
	for(int nn=-order-1;nn<order+1;++nn){
	    const int indexL = CCTK_GFINDEX3D(cctkGH,index[0]+nn*(dir==0),index[1]+nn*(dir==1),index[2]+nn*(dir==2));
	    B[1+order+nn]= F[indexL];
	}
	return InterpToFaces_LowLevel<order>(B);
}

template<int dir, int order>
inline static CCTK_REAL InterpToFaces_B(CCTK_POINTER_TO_CONST _GH, const int* lsh, const CCTK_REAL *F, int index[3]){
    if(order==0)
	return InterpToFaces<dir,0>(_GH,F,index);

    if(index[dir]> order && index[dir]< lsh[dir]-order){
	return InterpToFaces<dir,order>(_GH,F,index);
    }else{
    	if(index[dir]> order-1 && index[dir]< lsh[dir]-order+1){
		return InterpToFaces<dir,(order-1)*(order>0)>(_GH,F,index);
	}
	else{
    		if(index[dir]> order-2 && index[dir]< lsh[dir]-order+2 ){
			return InterpToFaces<dir,(order-2)*(order>1)>(_GH,F,index);
		}
	}
     }
	return InterpToFaces<dir,0>(_GH,F,index);
}
/* 
 * Old routine, not stable against round-off errors...
 * 
template<int order>
inline constexpr CCTK_REAL DER(CCTK_REAL* flux, int index[2*order+1]){
	return flux[index[0+order]]
	      - (order>0)* 1./24.*(flux[index[0+order-1]] + flux[index[2+order-1]]-2.*flux[index[1+order-1]])
	      + (order>1)* 3./640.*(flux[index[0]]-4.*flux[index[1]]+6.*flux[index[2]]-4.*flux[index[3]]+flux[index[4]]);
}
*/
template<int order>
inline constexpr CCTK_REAL DER(CCTK_REAL* flux, int index[2*order+1]){
        return (order==0)*flux[index[0+order]]
              + (order==1)* (-1./24.*(flux[index[0+order-1]] + flux[index[2+order-1]])+13./12.*flux[index[1+order-1]])
              + (order==2)* (3./640.*(flux[index[0]]+flux[index[4]])-29./480.*(flux[index[1]]+flux[index[3]]) + 1067./960.* flux[index[2]]);
}

template<int order,class T>
inline constexpr T DER_General(T* flux, int index[2*order+1]){
        return (order==0)*flux[index[0+order]]
              + (order==1)* (-1./24.*(flux[index[0+order-1]] + flux[index[2+order-1]])+13./12.*flux[index[1+order-1]])
              + (order==2)* (3./640.*(flux[index[0]]+flux[index[4]])-29./480.*(flux[index[1]]+flux[index[3]]) + 1067./960.* flux[index[2]]);
}


template<int order>
inline constexpr CCTK_REAL DER_comb(CCTK_REAL* flux, int index[2*(order+1)]){
        return (order==0)*(flux[index[1]]-flux[index[0]])
              + (order==1)* (9./8.*(flux[index[2]] - flux[index[1]])-1./24.*(flux[index[3]]-flux[index[0]]))
              + (order==2)* (75./64.*(flux[index[3]]-flux[index[2]])-25./384.*(flux[index[4]]-flux[index[1]]) + 3./640.* (flux[index[5]]-flux[index[0]]));
}

template<int dir, int order, class T>
inline static CCTK_REAL DER_Fd_General(CCTK_POINTER_TO_CONST _GH, CCTK_REAL* flux, int index[3]){
    cGH const * const cctkGH = (cGH const * const)_GH;

	int index_fd[2*order+1];
	for(int nn=-order;nn<order+1;++nn){
	    index_fd[order+nn] = CCTK_GFINDEX3D(cctkGH,index[0]+nn*(dir==0),index[1]+nn*(dir==1),index[2]+nn*(dir==2));
	}
	T fluxT[2*order+1];
	int indexT[2*order+1];
	for(int i=0;i<2*order+1;++i){fluxT[i]= (T) flux[index_fd[i]]; indexT[i]=i; }
	return DER_General<order,T>(fluxT,indexT);
}


template<int dir, int order>
inline static CCTK_REAL DER_Fd(CCTK_POINTER_TO_CONST _GH, CCTK_REAL* flux, int index[3]){
    cGH const * const cctkGH = (cGH const * const)_GH;

	int index_fd[2*order+1];
	for(int nn=-order;nn<order+1;++nn){
	    index_fd[order+nn] = CCTK_GFINDEX3D(cctkGH,index[0]+nn*(dir==0),index[1]+nn*(dir==1),index[2]+nn*(dir==2));
	}
	return DER<order>(flux,index_fd);
}
template<int dir, int order>
inline static CCTK_REAL DER_Fd_comb(CCTK_POINTER_TO_CONST _GH,CCTK_REAL* flux, int index[3]){
    cGH const * const cctkGH = (cGH const * const)_GH;

      int index_fd[2*(order+1)];
      for(int nn=-order-1;nn<order+1;++nn){
           index_fd[order+nn+1] = CCTK_GFINDEX3D(cctkGH, index[0]+nn*(dir==0),index[1]+nn*(dir==1),index[2]+nn*(dir==2));
      }
      return DER_comb<order>(flux,index_fd);
}

template<int dir, int order, int shiftL, int shiftR>
inline static CCTK_REAL DER_Bnd(CCTK_POINTER_TO_CONST _GH, const int *lsh, CCTK_REAL* flux, int index[3]){

    if(index[dir]-shiftL>= order && index[dir] +shiftR< lsh[dir]-order){
	return DER_Fd<dir, order>(_GH, flux,index);	
    }
    else{ //We are close to a boundary! decrease order
    	if(index[dir]-shiftL>= order-1 && index[dir]+ shiftR< lsh[dir]-order+1){
	   return DER_Fd<dir, (order-1)*(order>0)>(_GH, flux,index);	
        }
        else{ //We are close to a boundary! decrease order
    	   if(index[dir]-shiftL>= order-2 && index[dir] + shiftR< lsh[dir]-order+2 ){
	      return DER_Fd<dir, (order-2)*(order>1)>(_GH, flux,index);	
           }
        }
    }
    cGH const * const cctkGH = (cGH const * const)_GH;
    return flux[CCTK_GFINDEX3D(cctkGH,index[0],index[1],index[2])];
}

template<int dir, int order>
inline static CCTK_REAL DER_BndC(CCTK_POINTER_TO_CONST _GH, const int *lsh, CCTK_REAL* flux, int index[3]){
    if(order==0)
	return DER_Fd_comb<dir, 0>(_GH, flux,index);	

    if(index[dir]-1>= order && index[dir] < lsh[dir]-order){
	return DER_Fd_comb<dir, order>(_GH, flux,index);	
    }
    else{ //We are close to a boundary! decrease order
    	if(index[dir]-1>= order-1 && index[dir] < lsh[dir]-order+1){
	   return DER_Fd_comb<dir, (order-1)*(order>0)>(_GH, flux,index);	
        }
        else{ //We are close to a boundary! decrease order
    	   if(index[dir]-1>= order-2 && index[dir] < lsh[dir]-order+2 ){
	      return DER_Fd_comb<dir, (order-2)*(order>1)>(_GH, flux,index);	
           }
        }
    }
    return 1e300; //Something went really wrong....
}



#endif
