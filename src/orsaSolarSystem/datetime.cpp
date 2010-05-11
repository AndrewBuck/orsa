#include <orsaSolarSystem/datetime.h>

#include "sdncal.h"

#include <cmath>

#include <QDateTime>

using namespace orsa;
using namespace orsaSolarSystem;

// TimeScale data

// TODO: improve TimeScale code using Intervals and orsa::Time

struct TAI_minus_UTC_element {
    int day,month,year;
    double TAI_minus_UTC;
};

inline bool operator == (const TAI_minus_UTC_element & x, 
                         const TAI_minus_UTC_element & y) {
    if (x.day   != y.day)   return false;
    if (x.month != y.month) return false;
    if (x.year  != y.year)  return false;
    if (x.TAI_minus_UTC != y.TAI_minus_UTC) return false;
    return true;
}

inline bool operator != (const TAI_minus_UTC_element & x,
                         const TAI_minus_UTC_element & y) {
    return (!(x==y));
}

/*
 *  The TAI_minus_UTC_table should be updated periodically, 
 * for instance by consulting:
 * 
 * http://hpiers.obspm.fr/webiers/general/communic/publi/PUBLI.html
 * 
 * The present version is updated through Bulletin C17 (Jan 28, 1999)
 * Day, month, year, TAI-UTC (s)
 */

const TAI_minus_UTC_element TAI_minus_UTC_table_final_element = {0,0,0,0};
const TAI_minus_UTC_element TAI_minus_UTC_table[] = {
    {1,1,1972,10},
    {1,7,1972,11},
    {1,1,1973,12},
    {1,1,1974,13},
    {1,1,1975,14},
    {1,1,1976,15},
    {1,1,1977,16},
    {1,1,1978,17},
    {1,1,1979,18},
    {1,1,1980,19},
    {1,7,1981,20},
    {1,7,1982,21},
    {1,7,1983,22},
    {1,7,1985,23},
    {1,1,1988,24},
    {1,1,1990,25},
    {1,1,1991,26},
    {1,7,1992,27},
    {1,7,1993,28},
    {1,7,1994,29},
    {1,1,1996,30},
    {1,7,1997,31},
    {1,1,1999,32},
    {1,1,2006,33},
    {1,1,2009,34},
    {0,0,0,0} // TAI_minus_UTC_table_final_element
};

struct ET_minus_UT_element {
    int day,month,year;
    double ET_minus_UT;
};

inline bool operator == (const ET_minus_UT_element & x,
                         const ET_minus_UT_element & y) {
    if (x.day   != y.day)   return false;
    if (x.month != y.month) return false;
    if (x.year  != y.year)  return false;
    if (x.ET_minus_UT != y.ET_minus_UT) return false;
    return true;
}

inline bool operator != (const ET_minus_UT_element & x,
                         const ET_minus_UT_element & y) {
    return (!(x==y));
}

/* 
 * Values ET - UT at 0h UT of the date (specified as day, month, year)
 * from The Astronomical Almanac 1999, pages K8-9
 * or: http://www.stjarnhimlen.se/comp/time.html
 */

const ET_minus_UT_element ET_minus_UT_table_final_element = {0,0,0,0};
const ET_minus_UT_element ET_minus_UT_table[] = {
    {1,1,1800,13.7},
    {1,1,1801,13.4},
    {1,1,1802,13.1},
    {1,1,1803,12.9},
    {1,1,1804,12.7},
    {1,1,1805,12.6},
    {1,1,1806,12.5},
    {1,1,1816,12.5},
    {1,1,1817,12.4},
    {1,1,1818,12.3},
    {1,1,1819,12.2},
    {1,1,1820,12.0},
    {1,1,1821,11.7},
    {1,1,1822,11.4},
    {1,1,1823,11.1},
    {1,1,1824,10.6},
    {1,1,1825,10.2},
    {1,1,1826, 9.6},
    {1,1,1827, 9.1},
    {1,1,1828, 8.6},
    {1,1,1829, 8.0},
    {1,1,1830, 7.5},
    {1,1,1831, 7.0},
    {1,1,1832, 6.6},
    {1,1,1833, 6.3},
    {1,1,1834, 6.0},
    {1,1,1835, 5.8},
    {1,1,1836, 5.7},
    {1,1,1837, 5.6},
    {1,1,1838, 5.6},
    {1,1,1839, 5.6},
    {1,1,1840, 5.7},
    {1,1,1841, 5.8},
    {1,1,1842, 5.9},
    {1,1,1843, 6.1},
    {1,1,1844, 6.2},
    {1,1,1845, 6.3},
    {1,1,1846, 6.5},
    {1,1,1847, 6.6},
    {1,1,1848, 6.8},
    {1,1,1849, 6.9},
    {1,1,1850, 7.1},
    {1,1,1851, 7.2},
    {1,1,1852, 7.3},
    {1,1,1853, 7.4},
    {1,1,1854, 7.5},
    {1,1,1855, 7.6},
    {1,1,1856, 7.7},
    {1,1,1857, 7.7},
    {1,1,1858, 7.8},
    {1,1,1859, 7.8},
    {1,1,1860, 7.88},
    {1,1,1861, 7.82},
    {1,1,1862, 7.54},
    {1,1,1863, 6.97},
    {1,1,1864, 6.40},
    {1,1,1865, 6.02},
    {1,1,1866, 5.41},
    {1,1,1867, 4.10},
    {1,1,1868, 2.92},
    {1,1,1869, 1.82},
    {1,1,1870, 1.61},
    {1,1,1871, 0.10},
    {1,1,1872,-1.02},
    {1,1,1873,-1.28},
    {1,1,1874,-2.69},
    {1,1,1875,-3.24},
    {1,1,1876,-3.64},
    {1,1,1877,-4.54},
    {1,1,1878,-4.71},
    {1,1,1879,-5.11},
    {1,1,1880,-5.40},
    {1,1,1881,-5.42},
    {1,1,1882,-5.20},
    {1,1,1883,-5.46},
    {1,1,1884,-5.46},
    {1,1,1885,-5.79},
    {1,1,1886,-5.63},
    {1,1,1887,-5.64},
    {1,1,1888,-5.80},
    {1,1,1889,-5.66},
    {1,1,1890,-5.87},
    {1,1,1891,-6.01},
    {1,1,1892,-6.19},
    {1,1,1893,-6.64},
    {1,1,1894,-6.44},
    {1,1,1895,-6.47},
    {1,1,1896,-6.09},
    {1,1,1897,-5.76},
    {1,1,1898,-4.66},
    {1,1,1899,-3.74},
    {1,1,1900,-2.72},
    {1,1,1901,-1.54},
    {1,1,1902,-0.02},
    {1,1,1903, 1.24},
    {1,1,1904, 2.64},
    {1,1,1905, 3.86},
    {1,1,1906, 5.37},
    {1,1,1907, 6.14},
    {1,1,1908, 7.75},
    {1,1,1909, 9.13},
    {1,1,1910,10.46},
    {1,1,1911,11.53},
    {1,1,1912,13.36},
    {1,1,1913,14.65},
    {1,1,1914,16.01},
    {1,1,1915,17.20},
    {1,1,1916,18.24},
    {1,1,1917,19.06},
    {1,1,1918,20.25},
    {1,1,1919,20.95},
    {1,1,1920,21.16},
    {1,1,1921,22.25},
    {1,1,1922,22.41},
    {1,1,1923,23.03},
    {1,1,1924,23.49},
    {1,1,1925,23.62},
    {1,1,1926,23.86},
    {1,1,1927,24.49},
    {1,1,1928,24.34},
    {1,1,1929,24.08},
    {1,1,1930,24.02},
    {1,1,1931,24.00},
    {1,1,1932,23.87},
    {1,1,1933,23.95},
    {1,1,1934,23.86},
    {1,1,1935,23.93},
    {1,1,1936,23.73},
    {1,1,1937,23.92},
    {1,1,1938,23.96},
    {1,1,1939,24.02},
    {1,1,1940,24.33},
    {1,1,1941,24.83},
    {1,1,1942,25.30},
    {1,1,1943,25.70},
    {1,1,1944,26.24},
    {1,1,1945,26.77},
    {1,1,1946,27.28},
    {1,1,1947,27.78},
    {1,1,1948,28.25},
    {1,1,1949,28.71},
    {1,1,1950,29.15},
    {1,1,1951,29.57},
    {1,1,1952,29.97},
    {1,1,1953,30.36},
    {1,1,1954,30.72},
    {1,1,1955,31.07},
    {1,1,1956,31.35},
    {1,1,1957,31.68},
    {1,1,1958,32.18},
    {1,1,1959,32.68},
    {1,1,1960,33.15},
    {1,1,1961,33.59},
    {1,1,1962,34.00},
    {1,1,1963,34.47},
    {1,1,1964,35.03},
    {1,1,1965,35.73},
    {1,1,1966,36.54},
    {1,1,1967,37.43},
    {1,1,1968,38.29},
    {1,1,1969,39.20},
    {1,1,1970,40.18},
    {1,1,1971,41.17},
    {1,1,1972,42.23},
    {1,7,1972,42.80},
    {1,1,1973,43.37},
    {1,7,1973,43.93},
    {1,1,1974,44.49},
    {1,7,1974,44.99},
    {1,1,1975,45.48},
    {1,7,1975,45.97},
    {1,1,1976,46.46},
    {1,7,1976,46.99},
    {1,1,1977,47.52},
    {1,7,1977,48.03},
    {1,1,1978,48.53},
    {1,7,1978,49.06},
    {1,1,1979,49.59},
    {1,7,1979,50.07},
    {1,1,1980,50.54},
    {1,7,1980,50.96},
    {1,1,1981,51.38},
    {1,7,1981,51.78},
    {1,1,1982,52.17},
    {1,7,1982,52.57},
    {1,1,1983,52.96},
    {1,7,1983,53.38},
    {1,1,1984,53.79},
    {1,7,1984,54.07},
    {1,1,1985,54.34},
    {1,7,1985,54.61},
    {1,1,1986,54.87},
    {1,7,1986,55.10},
    {1,1,1987,55.32},
    {1,7,1987,55.57},
    {1,1,1988,55.82},
    {1,7,1988,56.06},
    {1,1,1989,56.30},
    {1,7,1989,56.58},
    {1,1,1990,56.86},
    {1,7,1990,57.22},
    {1,1,1991,57.57},
    {1,7,1991,57.94},
    {1,1,1992,58.31},
    {1,7,1992,58.72},
    {1,1,1993,59.12},
    {1,7,1993,59.55},
    {1,1,1994,59.98},
    {1,7,1994,60.38},
    {1,1,1995,60.78},
    {1,7,1995,61.20},
    {1,1,1996,61.63},
    {1,7,1996,61.96},
    {1,1,1997,62.29},
    {1,7,1997,62.63},
    {1,1,1998,62.97},
    {1,7,1998,63.22},
    {1,1,1999,63.47},
    {1,7,1999,63.66},
    {1,1,2000,63.82},
    {1,7,2000,63.98},
    {1,1,2001,64.09},
    {1,7,2001,64.21},
    {1,1,2002,64.30},
    {1,7,2002,64.41},
    {1,1,2003,64.47},
    {1,7,2003,64.55},
    {1,1,2004,64.57},
    {1,7,2004,64.65},
    {1,1,2005,64.69},
    {1,7,2005,64.80},
    {1,1,2006,64.85},
    {1,7,2006,64.99},
    {1,1,2007,65.15},
    {1,7,2007,65.34},
    {1,1,2008,65.45},
    {1,7,2008,65.63},
    {1,1,2009,65.70},
    {0,0,0,0} // ET_minus_UT_table_final_element
};

static double deltaSeconds(int y, int m, int d, 
                           const orsaSolarSystem::TimeScale ts) {
  
    double ds=0;
  
    switch (ts) {
    
        case TS_TDT: 
            ds=0;
            break;
    
        case TS_TAI:
            ds=32.184;
            break;
    
        case TS_UTC:
            ds = 32.184;
            if (y>=TAI_minus_UTC_table[0].year) {
                unsigned int j=0;
                TAI_minus_UTC_element candidate = TAI_minus_UTC_table[0];
                while (TAI_minus_UTC_table[j] != TAI_minus_UTC_table_final_element) {
                    if (TAI_minus_UTC_table[j].year   <=y) {
                        if (TAI_minus_UTC_table[j].month<=m) {
                            if (TAI_minus_UTC_table[j].day<=d) {
                                if ( ( (TAI_minus_UTC_table[j].year >  candidate.year) ) ||
                                     ( (TAI_minus_UTC_table[j].year == candidate.year) && (TAI_minus_UTC_table[j].month >  candidate.month) ) ||
                                     ( (TAI_minus_UTC_table[j].year == candidate.year) && (TAI_minus_UTC_table[j].month == candidate.month) && (TAI_minus_UTC_table[j].day > candidate.day) ) ) { 
                                    candidate = TAI_minus_UTC_table[j];
                                }		
                            }
                        }
                    }
                    j++;
                }
                ds += candidate.TAI_minus_UTC;
            }
            break;
    
        case TS_UT1:
            ds=0;
            if (y>=ET_minus_UT_table[0].year) {
                unsigned int j=0;
                ET_minus_UT_element candidate=ET_minus_UT_table[0];
                while (ET_minus_UT_table[j] != ET_minus_UT_table_final_element) {
                    if (ET_minus_UT_table[j].year   <=y) {
                        if (ET_minus_UT_table[j].month<=m) {
                            if (ET_minus_UT_table[j].day<=d) {
                                if ( ( (ET_minus_UT_table[j].year >  candidate.year) ) ||
                                     ( (ET_minus_UT_table[j].year == candidate.year) && (ET_minus_UT_table[j].month >  candidate.month) ) ||
                                     ( (ET_minus_UT_table[j].year == candidate.year) && (ET_minus_UT_table[j].month == candidate.month) && (ET_minus_UT_table[j].day > candidate.day) ) ) {
                                    candidate = ET_minus_UT_table[j];
                                }
                            }
                        }
                    }
                    j++;
                }
                ds += candidate.ET_minus_UT;
            }
            break;
    
        case TS_GPS:
            ds = 32.184 + 19.0;
            break;
    }
  
    // ORSA_DEBUG("delta-seconds = %g",ds);
  
    return (ds);
}


orsa::Time orsaSolarSystem::FromTimeScale(const orsa::Time & t,
                                          const orsaSolarSystem::TimeScale ts) {
    int y,m,d,H,M,S,ms;
    gregorDay(t,y,m,d,H,M,S,ms);
    const orsa::Time dt = orsa::Time(0,0,0,0,1000000*deltaSeconds(y,m,d,ts));
    return (t+dt);
}

orsa::Time orsaSolarSystem::ToTimeScale(const orsa::Time & t,
                                        const orsaSolarSystem::TimeScale ts) {
    int y,m,d,H,M,S,ms;
    gregorDay(t,y,m,d,H,M,S,ms);
    const orsa::Time dt = orsa::Time(0,0,0,0,1000000*deltaSeconds(y,m,d,ts));
    return (t-dt);
}

orsa::Time orsaSolarSystem::gregorTime(int y, 
                                       int m, 
                                       int d, 
                                       int H, 
                                       int M, 
                                       int S, 
                                       int ms) {
  
    // checks...
    while (ms >= 1000) {
        S  += 1;
        ms -= 1000;
    }
    while (S >= 60) {
        M += 1;
        S -= 60;
    }
    while (M >= 60) {
        H += 1;
        M -= 60;
    }
    while (H >= 24) {
        d += 1;
        H -= 24;
    }
  
    // more...
    while (ms < 0) {
        S  -= 1;
        ms += 1000;
    }
    while (S < 0) {
        M -= 1;
        S += 60;
    }
    while (M < 0) {
        H -= 1;
        M += 60;
    }
    while (H < 0) {
        d -= 1;
        H += 24;
    }
  
    return Time(GregorianToSdn(y,m,d),
                0,
                0,
                0,
                mpz_class(mpz_class(mpz_class(H*60 + M)*60 + S)*1000 + ms)*1000);
}

orsa::Time orsaSolarSystem::gregorTime(int y, 
                                       int m, 
                                       double d) {
    // test with negative years!!
    // const int int_d = (int)d;
    const int int_d = (int)(floor(d));
    return gregorTime(y,
                      m,
                      int_d,
                      0,
                      0,
                      0,
                      (int)((d-int_d)*86400000)); // this las number, the ms, should be always in the valid range for the int type
}

void orsaSolarSystem::gregorDay(const orsa::Time & t,
                                int & y,
                                int & m,
                                int & d,
                                int & H,
                                int & M,
                                int & S,
                                int & ms) {
    const mpz_class muDay("86400000000");
    const mpz_class muHour("3600000000");
    const mpz_class muMinute("60000000");
    const mpz_class muSecond( "1000000");
  
    SdnToGregorian(mpz_class(t.getMuSec() / muDay).get_si(),
                   &y,
                   &m,
                   &d);
    mpz_class dayFractionMuSec = mpz_class(t.getMuSec() % muDay);
    //
    H = mpz_class(dayFractionMuSec / muHour).get_si();
    dayFractionMuSec -= H * muHour;
    M = mpz_class(dayFractionMuSec / muMinute).get_si();
    dayFractionMuSec -= M * muMinute;
    S = mpz_class(dayFractionMuSec / muSecond).get_si();
    dayFractionMuSec -= S * muSecond;
    ms = mpz_class(dayFractionMuSec / 1000).get_si();
    dayFractionMuSec -= ms * 1000;
}

void orsaSolarSystem::gregorDay(const orsa::Time & t,
                                int & y,
                                int & m,
                                int & d,
                                double & fd) {
  
    const mpz_class muDay("86400000000");
  
    SdnToGregorian(mpz_class(t.getMuSec() / muDay).get_si(),
                   &y,
                   &m,
                   &d);
    mpz_class dayFractionMuSec = mpz_class(t.getMuSec() % muDay);
    fd = dayFractionMuSec.get_d() / muDay.get_d();
}

double orsaSolarSystem::timeToJulian(const orsa::Time & t) {
    return orsa::FromUnits(FromUnits(t.getMuSec().get_d(), orsa::Unit::MICROSECOND), orsa::Unit::DAY,-1) - 0.5;
}

orsa::Time orsaSolarSystem::julianToTime(const double & jd) {
    return orsa::Time(orsa::FromUnits(orsa::FromUnits(jd+0.5, orsa::Unit::DAY), orsa::Unit::MICROSECOND,-1));
}

orsa::Time orsaSolarSystem::J2000() {
  
    /*! 
      IAU definition: [2000 Jan 1d 12h TDT]
      http://en.wikipedia.org/wiki/Terrestrial_Time 
      http://en.wikipedia.org/wiki/Month 
      http://aa.usno.navy.mil/software/novas/novas_c/novasc_info.html
      http://nanotitan.com/software/api/suite/diamond/nT/quantity/constant/TIME_INSTANT.html#J2000
    */
  
    return orsaSolarSystem::FromTimeScale(gregorTime(2000,
                                                     1,
                                                     1,
                                                     12,
                                                     0,
                                                     0,
                                                     0),
                                          orsaSolarSystem::TS_TDT);
}

/* 
   orsa::Time orsaSolarSystem::now() {
   const time_t tt_now = time(0);
   struct tm * tm_struct = gmtime(&tt_now);
   return orsaSolarSystem::FromTimeScale(gregorTime(1900+tm_struct->tm_year,
   1+tm_struct->tm_mon,
   tm_struct->tm_mday,
   tm_struct->tm_hour,
   tm_struct->tm_min,
   tm_struct->tm_sec,
   0),
   orsaSolarSystem::TS_UTC);
   }
*/

orsa::Time orsaSolarSystem::now() {
    const QDateTime datetime = QDateTime::currentDateTime().toTimeSpec(Qt::UTC);
    const QDate         date = datetime.date();
    const QTime         time = datetime.time();
    return orsaSolarSystem::FromTimeScale(gregorTime(date.year(),
                                                     date.month(),
                                                     date.day(),
                                                     time.hour(),
                                                     time.minute(),
                                                     time.second(),
                                                     time.msec()),
                                          orsaSolarSystem::TS_UTC);
}

/*****/

void orsaSolarSystem::ApproximatedLunation(const orsa::Time & t,
                                           orsa::Time & begin,
                                           orsa::Time & end,
                                           int & lunationID) {
  
    // data for new moon epoch: http://en.wikipedia.org/wiki/New_moon
  
    const orsa::Time t0 =  
        orsaSolarSystem::FromTimeScale(orsaSolarSystem::gregorTime(2000,1,1,0,0,0,0), 
                                       orsaSolarSystem::TS_TT);

    const double d0 = 5.597661;
    const double d1 = 29.5305888610;
    const double d2 = 102.026e-12;
    // better initial guess?
    int N=((t-t0).get_d()-orsa::FromUnits(d0,orsa::Unit::DAY))/orsa::FromUnits(d1,orsa::Unit::DAY);
    while (1) {
        // ORSA_DEBUG("N: %i",N);
        // original, new moon d:
        // d = d0 + d1*N + d2*N*N; 
        //
        // modified, full moon d: N becomes N-1/2 or N+1/2
        const double d_begin = d0 + d1*(N-0.5) + d2*(N*N-N+0.25); 
        const double d_end   = d0 + d1*(N+0.5) + d2*(N*N+N+0.25); 
        //
        begin = t0 + orsa::Time(orsa::FromUnits(orsa::FromUnits(d_begin,orsa::Unit::DAY),orsa::Unit::MICROSECOND,-1));
        end   = t0 + orsa::Time(orsa::FromUnits(orsa::FromUnits(d_end,  orsa::Unit::DAY),orsa::Unit::MICROSECOND,-1));
        //
        lunationID = N;
        //
        if (begin > t) {
            N--; continue;
        } else if (end < t) {
            N++; continue;
        } else {
            break;
        }
    }
}

/*****/

double orsaSolarSystem::fractionalYear(const orsa::Time & t) {
    int y,m,d;
    double fd;
    orsaSolarSystem::gregorDay(t,y,m,d,fd);
    // NOTE: might not work around year zero....
    orsa::Time t0 = orsaSolarSystem::gregorTime(y,  1,1);
    orsa::Time t1 = orsaSolarSystem::gregorTime(y+1,1,1);
    return (y+((t-t0).get_d()/(t1-t0).get_d()));
}
