#ifndef _ORSA_DATETIME_H_
#define _ORSA_DATETIME_H_

#include <orsa/double.h>
#include <orsa/debug.h>
#include <orsa/unit.h>

namespace orsa {

    class Time {
    public:
        Time() { }
    public:
        Time(const Time & t) : _mu_sec(t._mu_sec) { }
    public:
        Time(const mpz_class mu_sec) : _mu_sec(mu_sec) { }
    public:
        Time(const mpz_class d,
             const mpz_class H,
             const mpz_class M,
             const mpz_class S,
             const mpz_class mu_S) {
            _mu_sec = mu_S + mpz_class("1000000") *
                (S + mpz_class("60") *
                 (M + mpz_class("60") *
                  (H + mpz_class("24") * d)));
        }
    public:
        /* inline Time & operator = (const mpz_class & rhs) {
           if (rhs != 0) {
           ORSA_DEBUG("warning: setting time from non-zero integer = %Zi",rhs.get_mpz_t());
           }
           _mu_sec = rhs;
           return * this;
           }
        */
    public:
        ~Time() { }
    public:
        bool set(const mpz_class d,
                 const mpz_class H,
                 const mpz_class M,
                 const mpz_class S,
                 const mpz_class mu_S) {
            _mu_sec = mu_S + mpz_class("1000000") *
                (S + mpz_class("60") *
                 (M + mpz_class("60") *
                  (H + mpz_class("24") * d)));
            return true;
        }
    public:
        const mpz_class & getMuSec() const {
            return _mu_sec;
        }
    public:
        double get_d() const {
            return (FromUnits(_mu_sec.get_d(),Unit::MICROSECOND));
        }
    public:
        inline Time & operator += (const Time & rhs) {
            _mu_sec += rhs._mu_sec;
            return * this;
        }
    public:
        inline Time & operator -= (const Time & rhs) {
            _mu_sec -= rhs._mu_sec;
            return * this;
        }
    public:
        inline Time & operator *= (const mpz_class & rhs) {
            _mu_sec *= rhs;
            return * this;
        }
    
    public:
        inline const Time operator + (const Time & rhs) const {
            Time _t(*this);
            _t += rhs;
            return _t;
        }
    public:
        inline const Time operator - (const Time & rhs) const {
            Time _t(*this);
            _t -= rhs;
            return _t;
        }
    
    public:
        inline const Time operator / (const mpz_class & rhs) const {
            Time _t(*this);
            _t._mu_sec /= rhs;
            return _t;
        }

    public:
        inline const Time operator + () const {
            Time _t(*this);
            return _t;
        }
    public:
        inline const Time operator - () const {
            Time _t(0);
            _t -= (*this);
            return _t;
        }

    public:
        inline bool operator < (const Time & rhs) const {
            return (_mu_sec < rhs._mu_sec);
        }
    public:
        inline bool operator <= (const Time & rhs) const {
            return (_mu_sec <= rhs._mu_sec);
        }
        inline bool operator > (const Time & rhs) const {
            return (_mu_sec > rhs._mu_sec);
        }
    public:
        inline bool operator >= (const Time & rhs) const {
            return (_mu_sec >= rhs._mu_sec);
        }
    public:
        inline bool operator == (const Time & rhs) const {
            return (_mu_sec == rhs._mu_sec);
        }
    public:
        inline bool operator != (const Time & rhs) const {
            return (_mu_sec != rhs._mu_sec);
        }
    
    protected:
        mpz_class _mu_sec;
    };

    inline const Time operator * (const Time      & lhs,
                                  const mpz_class & rhs) {
        return Time(lhs.getMuSec()*rhs);
    }

    inline const Time operator * (const mpz_class & lhs,
                                  const Time      & rhs) {
        return Time(lhs*rhs.getMuSec());
    }

} // namespace orsa

#endif // _ORSA_DATETIME_H_
