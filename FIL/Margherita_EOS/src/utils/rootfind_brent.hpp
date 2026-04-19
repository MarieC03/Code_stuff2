#ifndef BRENT_ROOTFIND_HH__
#define BRENT_ROOTFIND_HH__

#include <cmath>
#include <limits>

namespace rootfind 
{
enum retcode_t {
    SUCCESS=0,
    INVINPUT,
    NOCONV
} ;

constexpr unsigned long const maxiter = 150 ;

template<typename F>
inline retcode_t 
zero_brent(const double& a, const double& b, const double tol, F&& func, double& z)
{
    double c;
    double d;
    double e;
    double fa;
    double fb;
    double fc;
    double m;
    double p;
    double q;
    double r;
    double s;
    double sa;
    double sb;
    double tol;
    //
    //  Make local copies of A and B.
    //
    sa = a;
    sb = b;
    fa = f(sa);
    fb = f(sb);

    c = sa;
    fc = fa;
    e = sb - sa;
    d = e;

    if( fa == 0 )
    {
        z = a; return SUCCESS ;
    } else if ( fb == 0 )
    {
        z = b; return SUCCESS ; 
    } else if ( fa*fb > 0 )
    {
        z = std::numeric_limits<double>::signaling_NaN ;
        return INVINPUT ;
    }

    constexpr double macheps = std::numeric_limits<double>::epsilon();
    unsigned long niter = 0 ;

    while(niter < maxiter)
    {
        if (std::fabs(fc) < std::fabs(fb)) {
            sa = sb;
            sb = c;
            c = sa;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        tol = 2.0 * macheps * std::fabs(sb) + t;
        m = 0.5 * (c - sb);

        if (std::fabs(m) <= tol || fb == 0.0) {
        break;
        }

        if (std::fabs(e) < tol || std::fabs(fa) <= std::fabs(fb)) {
            e = m;
            d = e;
        } else {
            s = fb / fa;

            if (sa == c) {
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * m * q * (q - r) - (sb - sa) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }

            if (0.0 < p) {
                q = -q;
            } else {
                p = -p;
            }

            s = e;
            e = d;

            if (2.0 * p < 3.0 * m * q - std::fabs(tol * q) &&
                p < std::fabs(0.5 * s * q)) {
                d = p / q;
            } else {
                e = m;
                d = e;
            }
        }
        sa = sb;
        fa = fb;

        if (tol < std::fabs(d)) {
        sb = sb + d;
        } else if (0.0 < m) {
        sb = sb + tol;
        } else {
        sb = sb - tol;
        }

        fb = f(sb);

        if ((0.0 < fb && 0.0 < fc) || (fb <= 0.0 && fc <= 0.0)) {
            c = sa;
            fc = fa;
            e = sb - sa;
            d = e;
        }
        niter ++ ; 
    }

    if( niter < maxiter )
    {
        z = sb; return SUCCESS ;
    }
    return NOCONV;
}

}


#endif