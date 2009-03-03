#ifndef _ORSA_REFERENCE_FRAME_
#define _ORSA_REFERENCE_FRAME_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <list>

#include <orsa/datetime.h>
#include <orsa/matrix.h>

namespace orsa {
  
  class ReferenceFrame {
  public:
    virtual ReferenceFrame * instance() = 0;
  public:
    virtual bool instanciated() = 0; 
  protected:
    ReferenceFrame() { }
  public:
    virtual ~ReferenceFrame() { }
  protected:
    
  public:
    virtual orsa::Matrix relativeToAbsolute(const orsa::Time &) const = 0;  
    virtual orsa::Matrix absoluteToRelative(const orsa::Time &) const = 0;
  };
  
  /* Use this macro every time you subclass from ReferenceFrame,
   * and remember to declare:
   * ReferenceFrame * className::_instance = 0;
   * in the .cpp file
   * see the AbsoluteReferenceFrame example below
   */
#define __ORSA_REFERENCE_FRAME_MACRO__(className)	\
  public:						\
    ReferenceFrame * instance() {			\
    if (_instance == 0) {				\
      _instance = new className;			\
    }							\
    return _instance;					\
  }							\
 public:						\
  bool instanciated() {					\
    return (_instance != 0);				\
  }							\
 protected:						\
  className();						\
 public:						\
  virtual ~className() {				\
    _instance = 0;					\
  }							\
 protected:						\
  static ReferenceFrame * _instance;			
  
  /*****/
  
  class AbsoluteReferenceFrame : public ReferenceFrame {
    
    __ORSA_REFERENCE_FRAME_MACRO__(AbsoluteReferenceFrame);
    
  public:
    orsa::Matrix relativeToAbsolute(const orsa::Time &) const {
      return _identity;
    }
    
  public:
    orsa::Matrix absoluteToRelative(const orsa::Time &) const {
      return _identity;
    }
    
  private:
    const orsa::Matrix _identity;
  };
  
} // namespace orsa

#endif // _ORSA_REFERENCE_FRAME_
