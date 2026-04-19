#ifndef __HH__REC_HELPERS
#define __HH__REC_HELPERS



/*
 *
 *   Some helper function to do higher order reconstructions
 *   as in ECHO. Mostly unused now...
 *
 *   Written by Elias R. Most in Feb. 2017
 *
 */


static inline constexpr int sgn(const CCTK_REAL& A){return (0 < A) - (A < 0);};

static inline constexpr CCTK_REAL max(const CCTK_REAL& A, const CCTK_REAL& B){ return (A>B) ? A : B;};

template<typename... Targs>
inline constexpr double max(const double& A, const double & B, Targs... Fargs) 
{
        return max(A,max(B,Fargs...));
}


inline constexpr CCTK_REAL min(const CCTK_REAL& A, const CCTK_REAL& B){ return (A>B) ? B : A;};

template<typename... Targs>
inline constexpr double min(const double& A, const double & B, Targs... Fargs) 
{
        return min(A,min(B,Fargs...));
}


static inline CCTK_REAL minmod(const CCTK_REAL& A, const CCTK_REAL& B){ return 0.5*(sgn(A)+sgn(B))*min(std::fabs(A),std::fabs(B));};

static inline  CCTK_REAL minmod4(const CCTK_REAL& A, const CCTK_REAL& B,const CCTK_REAL& C, const CCTK_REAL& D){
	return 0.125*(sgn(A)+sgn(B))*fabs((sgn(A)+sgn(C))*(sgn(A)+sgn(D)))*min(std::fabs(A),std::fabs(B),std::fabs(C), std::fabs(D));};

static inline  CCTK_REAL minmod3(const CCTK_REAL& A, const CCTK_REAL& B,const CCTK_REAL& C){
	return 0.25*(sgn(A)+sgn(B))*fabs((sgn(A)+sgn(C)))*min(std::fabs(A),std::fabs(B),std::fabs(C));};

static inline CCTK_REAL median(const CCTK_REAL& A, const CCTK_REAL& B, const CCTK_REAL& C){ return A + minmod(B-A,C-A);};

template<bool left,int offset>
static inline CCTK_REAL minmod_REC(CCTK_REAL U[MAXNUMVARS][MAXNUMINDICES], const int whichvar ){
	  const CCTK_REAL deltap = U[whichvar][PLUS1-2*left+offset]- U[whichvar][PLUS0 + offset];
	  const CCTK_REAL deltam = U[whichvar][PLUS0 + offset]- U[whichvar][MINUS1+2*left + offset];
 	  return U[whichvar][PLUS0 + offset] +0.5*minmod(deltam,deltap);
};

static inline CCTK_REAL mc2(const CCTK_REAL& A, const CCTK_REAL& B){ 
  return 0.5*(sgn(A)+sgn(B))*min(2.*std::fabs(A),2.*std::fabs(B), 0.5*std::fabs(A+B));};

template<bool left,int offset>
static inline CCTK_REAL mc2_REC(CCTK_REAL U[MAXNUMVARS][MAXNUMINDICES], const int whichvar ){
	  const CCTK_REAL deltap = U[whichvar][PLUS1-2*left+offset]- U[whichvar][PLUS0 + offset];
	  const CCTK_REAL deltam = U[whichvar][PLUS0 + offset]- U[whichvar][MINUS1+2*left + offset];
 	  return U[whichvar][PLUS0 + offset] +0.5*mc2(deltam,deltap);
};


inline constexpr CCTK_REAL SQ(const CCTK_REAL &X){ return X*X;};



#endif
