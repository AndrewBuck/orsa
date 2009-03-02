#ifndef _ORSAQT_EVENT_
#define _ORSAQT_EVENT_

#include <QEvent>

#include <orsaSolarSystem/orbit.h>

namespace orsaQt {
  
  // IMPORTANT: make sure that each event has its own type ID (the argument of QEvent constructor)
  
  // define all user events here
  
#define __ORSAQT_EVENT_DEBUG__          QEvent::Type(QEvent::User + 1)
#define __ORSAQT_EVENT_ORBIT_MULTIFIT__ QEvent::Type(QEvent::User + 2)
  
  class DebugEvent : public QEvent {
  public:
    DebugEvent(const QString & s) : 
      QEvent(__ORSAQT_EVENT_DEBUG__), 
      msg(s) { }
  public:
    const QString & debugMessage() const { 
      return msg; 
    }
  private:
    const QString msg;
  };  
  
  class OrbitMultifitEvent : public QEvent {
  public:
    OrbitMultifitEvent(const orsaSolarSystem::OrbitWithEpoch & o,
		       const double f,        // f = max(1,sqrt(chi/sqrt(dof))
		       const gsl_matrix * covarianceMatrix) :
      QEvent(__ORSAQT_EVENT_ORBIT_MULTIFIT__),
      orbit(o),
      factor(f) {
      
      ORSA_DEBUG("covar size: %i x %i +++++++++++++++++++++",
		 covarianceMatrix->size1,
		 covarianceMatrix->size2);
      
      covar = gsl_matrix_alloc(covarianceMatrix->size1,
			       covarianceMatrix->size2);
      
      gsl_matrix_memcpy(covar,covarianceMatrix);
    }
  public:
    ~OrbitMultifitEvent() {
      gsl_matrix_free(covar);
    }
  public:
    const orsaSolarSystem::OrbitWithEpoch orbit;
    const double factor;
    gsl_matrix * covar;
  };
  
} // namespace orsaQt

#endif // _ORSAQT_EVENT_
