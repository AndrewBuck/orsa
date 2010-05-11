#ifndef _ORSA_LEGENDRE_
#define _ORSA_LEGENDRE_

#include <orsa/cache.h>
#include <orsa/double.h>

#include <vector>

namespace orsa {
  
    //! Legendre polynomials
    //! NOTE: always 0 <= m <= l, i.e. non-negative m
    class Legendre {
    public:
        Legendre(const double x) : _x(x), _sqrt_one_minus_x2(sqrt(1-x*x)) {
            // ORSA_DEBUG("RecursiveLegendre created, x = %g",_x);
        }
    public:
        Legendre(const double x, const unsigned int hint_order) : _x(x), _sqrt_one_minus_x2(sqrt(1-x*x)) { 
            // ORSA_DEBUG("RecursiveLegendre created, x = %g",_x);
            _P.resize(1+hint_order);
            _dP.resize(1+hint_order);
            for (unsigned int j=0; j<(1+hint_order); ++j) {
                _P[j].resize(1+j);
                _dP[j].resize(1+j);
            }
        }
    protected:
        void __check_P__(const unsigned int l, 
                         const unsigned int m) const {
            // ORSA_DEBUG("called check_P(%i,%i)",l,m);
            if (m > l) {
                return;
            }
            if (_P.size() < (1+l)) {
                _P.resize(1+l);
                for (unsigned int j=0; j<(1+l); ++j) {
                    _P[j].resize(1+j);
                }
            }
            if (_P[l][m].isSet()) {
                // ORSA_DEBUG("check: P(%i,%i) is set to: %g",l,m,_P[l][m].get());
                return;
            }
            if (l==m) {
                //
                // first sign convention
                // _P[l][m].set(power_sign(l)*__fact_fact__(2.0*l-1.0)*__int_pow__(_sqrt_one_minus_x2,l));
                // _P[l][m].set(power_sign(l)*bi_factorial(2.0*l-1.0)*__int_pow__(_sqrt_one_minus_x2,l));
                // _P[l][m].set(power_sign(l)*bi_factorial(2*mpz_class(l)-1)*int_pow(_sqrt_one_minus_x2,l));
                // second sign convention (-1)^m
                // _P[l][m].set(power_sign(m)*power_sign(l)*__fact_fact__(2.0*l-1.0)*__int_pow__(_sqrt_one_minus_x2,l));
                // _P[l][m].set(power_sign(m)*power_sign(l)*bi_factorial(2.0*l-1)*int_pow(_sqrt_one_minus_x2,l));
                // _P[l][m].set(power_sign(m)*power_sign(l)*bi_factorial(2*l-1)*int_pow(_sqrt_one_minus_x2,l));
                _P[l][m].set(power_sign(m)*power_sign(l)*bi_factorial(2*mpz_class(l)-1).get_d()*int_pow(_sqrt_one_minus_x2,l));
                //
                // ORSA_DEBUG("check: P(%i,%i) = %g",l,m,_P[l][m].get());
                return;
            }
            if (l==(m+1)) {
                __check_P__(m,m);
                //
                // there is not a sign dependence here
                // first sign convention
                _P[l][m].set(_x*(2.0*m+1)*_P[m][m].getRef());
                // second sign convention (-1)^m
                // _P[l][m].set(_x*(2.0*m+1.0)*_P[m][m].get());
                //
                // ORSA_DEBUG("check: P(%i,%i) = %g",l,m,_P[l][m].get());
                return;
            }
            /* NO negative m
               if (_P[l][-m].isSet()) {
               // _P[l][m].set(power_sign(m)*(__fact__(l-m)/__fact__(l+m))*_P[l][-m].get());
               _P[l][m].set(power_sign(m)*(factorial(l-m)/__fact__(l+m))*_P[l][-m].get());
               return;
               }
            */
            __compute_P__(l,m);
        }
    protected:
        void __check_dP__(const unsigned int l, 
                          const unsigned int m) const {
            // ORSA_DEBUG("called check_dP(%i,%i)",l,m);
            if (m > l) {
                return;
            }
            if (_dP.size() < (1+l)) {
                _dP.resize(1+l);
                for (unsigned int j=0; j<(1+l); ++j) {
                    _dP[j].resize(1+j);
                }
            }
            if (_dP[l][m].isSet()) {
                // ORSA_DEBUG("check_dP: is set...");
                return;
            }
            __compute_dP__(l,m);
        }
    protected:
        void __compute_P__(const unsigned int l, 
                           const unsigned int m) const {
            // ORSA_DEBUG("called compute_P(%i,%i)",l,m);
            if (_P[l][m].isSet()) {
                return;
            }
            /* 
               __check_P__(l-1,m); 
               __check_P__(l-2,m); 
               _P[l][m].set((_x*(2.0*l-1.0)*_P[l-1][m].get()-(l+m-1.0)*_P[l-2][m].get())/(l-m));
            */
            //
            if ((l>=2) && ((l+2)>=m)) {
                __check_P__(l-1,m); 
                __check_P__(l-2,m); 
                _P[l][m].set((_x*(2.0*l-1.0)*_P[l-1][m].getRef()-(l+m-1.0)*_P[l-2][m].getRef())/(l-m));
            } else if ((l>=1) && ((l+1)>=m)) {
                __check_P__(l-1,m); 
                _P[l][m].set((_x*(2.0*l-1.0)*_P[l-1][m].getRef())/(l-m));
            } else {
                ORSA_ERROR("this case should have been handled already somewhere else");
            }	
            // ORSA_DEBUG("computed value: P(%i,%i) = %g",l,m,_P[l][m].get());
        }
    protected:	
        void __compute_dP__(const unsigned int l, 
                            const unsigned int m) const {
            //  ORSA_DEBUG("called compute_dP(%i,%i)",l,m);
            if (_dP[l][m].isSet()) {
                return;
            }
            if ((l>0) && (l>m)) {
                __check_P__(l,m);
                __check_P__(l-1,m);
                // _dP[l][m].set((l*_x*_P[l][m].getRef()-(l+m)*_P[l-1][m].getRef())/_sqrt_one_minus_x2);
                if (fabs(_sqrt_one_minus_x2) > epsilon()) {
                    _dP[l][m].set((l*_x*_P[l][m].getRef()-(l+m)*_P[l-1][m].getRef())/_sqrt_one_minus_x2);
                } else {
                    ORSA_DEBUG("PROBLEM: division by zero... FIX FIX FIX (how?)");
                    //
                    _dP[l][m].set(0);
                }
            } else {
                __check_P__(l,m);
                // _dP[l][m].set((l*_x*_P[l][m].getRef())/_sqrt_one_minus_x2);
                if (fabs(_sqrt_one_minus_x2) > epsilon()) {
                    _dP[l][m].set((l*_x*_P[l][m].getRef())/_sqrt_one_minus_x2);
                } else {
                    ORSA_DEBUG("PROBLEM: division by zero... FIX FIX FIX (how?)");
                    //
                    _dP[l][m].set(0);
                }
            }
            // ORSA_DEBUG("computed value: dP(%i,%i) = %g",l,m,_dP[l][m].get());
        }
    public:
        const double P(const unsigned int l, 
                       const unsigned int m) const {
            // ORSA_DEBUG("called P(%i,%i)",l,m);
            if (m > l) {
                ORSA_ERROR("incorrect values: m > l (%i > %i)",m,l);
                return 0;
            }
            if (_P.size() < (1+l)) {
                _P.resize(1+l);
                for (unsigned int j=0; j<(1+l); ++j) {
                    _P[j].resize(1+j);
                }
            }
            /* 
               if (0) {
               // debug
               for (unsigned int j=0; j<_P.size(); ++j) {
               for (unsigned int k=0; k<_P[j].size(); ++k) {
               if (_P[j][k].isSet()) {
               ORSA_DEBUG("P(%i,%i) = %g",j,k,_P[j][k].get());
               } else {
               ORSA_DEBUG("P(%i,%i) is unset...",j,k);
               }
               }
               }
               }
            */
            //
            if (_P[l][m].isSet()) {
                return _P[l][m].getRef();
            }
            __check_P__(l,m);
            return _P[l][m].getRef();
        }
    public:
        //! dP returns the value of dP/d theta, assuming P=P(cos theta);
        const double dP(const unsigned int l, 
                        const unsigned int m) const {
            if (m > l) {
                ORSA_ERROR("incorrect values: l < m (%i < %i)",l,m);
                return 0;
            }
            if (_dP.size() < (1+l)) {
                _dP.resize(1+l);
                for (unsigned int j=0; j<(1+l); ++j) {
                    _dP[j].resize(1+j);
                }
            }
            if (_dP[l][m].isSet()) {
                return _dP[l][m].getRef();
            }
            __check_dP__(l,m);
            return _dP[l][m].getRef();
        }
    public:
        /* 
           static double norm(const unsigned int l,
           const unsigned int m) {
        */
        //
        static double norm(const mpz_class & l,
                           const mpz_class & m);
        /* 
           static double norm(const mpz_class & l,
           const mpz_class & m) {
           // return sqrt((__fact__(l-m)*(2.0*l+1.0)*(2.0-__kronecker__(0,m)))/(__fact__(l+m)));
           // return sqrt((factorial(l-m)*(2.0*l+1.0)*(2.0-__kronecker__(0,m)))/(__fact__(l+m)));
           //
           // return sqrt(double(factorial(l-m)*(2-kronecker(0,m)))/double(factorial(l+m)));
           //
           // return sqrt(double(factorial(l-m)*(2-kronecker(0,m))*(2*l-1))/double(factorial(l+m)));
           //
           // ORSA_DEBUG("pre-sqrt: %Fg",double(double(factorial(l+m))/double(factorial(l-m)*(2-kronecker(0,m))*(2*l+1)))());
           //
           // almost same as paper
           // return sqrt(double(factorial(l+m))/double(factorial(l-m)*(2-kronecker(0,m))*(2*l+1)));
           //
           // another test
           return sqrt(1/double((2-kronecker(0,m))*(2*l+1)));
           } 
        */
   
    private:
        // Legendre(l,m)
        mutable std::vector<std::vector<orsa::Cache<double> > > _P;
        mutable std::vector<std::vector<orsa::Cache<double> > > _dP;
    private:
        const double _x;
        const double _sqrt_one_minus_x2;
    };
  
} // namespace orsa

#endif // _ORSA_LEGENDRE_
