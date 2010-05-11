#ifndef _ORSA_PAUL_H_
#define _ORSA_PAUL_H_

#include <orsa/bodygroup.h>

namespace orsa {
  
    class Attitude;
    class PaulMoment;
  
    //! A class implementing the equations in M.K. Paul (1988), "An expansion in power series of mutual potential for gravitating bodies with finite sizes", Celestial Mechanics, 44, 49-59.
  
    class Paul {
    public:
        // Eq. (1')
        static double gravitationalPotential(const orsa::PaulMoment * M1,
                                             const orsa::Matrix     & A1_g2l,
                                             const orsa::PaulMoment * M2,
                                             const orsa::Matrix     & A2_g2l,
                                             const orsa::Vector     & R);
    
    public:
        static orsa::Vector gravitationalForce(const orsa::PaulMoment * M1,
                                               const orsa::Matrix     & A1_g2l,
                                               const orsa::PaulMoment * M2,
                                               const orsa::Matrix     & A2_g2l,
                                               const orsa::Vector     & R,
                                               const bool               thisIsRepeatedCall=false);
    
    public:
        static orsa::Vector gravitationalTorque(const orsa::PaulMoment * M1,
                                                const orsa::Matrix     & A1_g2l,
                                                const orsa::PaulMoment * M2,
                                                const orsa::Matrix     & A2_g2l,
                                                const orsa::Vector     & R);
    
    public:
        // Eq. (21b)
        static double C_lmn(const int l,
                            const int m,
                            const int n);
    
    public:	
        // Eq. (23)
        class t_lmnLMN {
      
        public:
            static t_lmnLMN * instance();
        protected:
            static t_lmnLMN * _instance;
      
        protected:
            t_lmnLMN();
        public:
            virtual ~t_lmnLMN();
      
        public:
            double get(const int l,
                       const int m,
                       const int n,
                       const int L,
                       const int M,
                       const int N) const;
        private:
            double trueGet(const int l,
                           const int m,
                           const int n,
                           const int L,
                           const int M,
                           const int N) const;
      
        private:
            void resize(const size_t order) const;
        private:
            /* void sort(int & l_out,
               int & m_out,
               int & n_out,
               const int l_in,
               const int m_in,
               const int n_in) const;
            */
        private:
            mutable 
            std::vector<
            std::vector<
                std::vector<
                    std::vector<
                        std::vector<
                            std::vector<
                                orsa::Cache<double> 
                                > 
                            > 
                        > 
                    > 
                > 
            > _data;
        };
    
    };
  
}; // namespace orsa

#endif // _ORSA_PAUL_H_
