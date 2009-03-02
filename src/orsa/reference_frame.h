#ifndef _ORSA_REFERENCE_FRAME_
#define _ORSA_REFERENCE_FRAME_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <list>

#include <orsa/datetime.h>
#include <orsa/matrix.h>

namespace orsa {
  
  /* 
     class Body;
     
     class ReferenceFrame {
     
     public:
     enum ReferenceFrameType {
     RELATIVE_RF,
     ABSOLUTE_RF
     };
     
     protected:
     osg::ref_ptr<Body> _parentBody;
     };
  */
  
  /* 
   *
   * IDEA per implementare ReferenceFrame: qui si definisce solo il sistema "assoluto",
   * e poi in SolarSystem si inseriscono 2 sistemi: Ecliptic == assoluto, ed Equatorial,
   * con le relative chiamate from() e to() ad una determinata epoca;
   *
   */
  
  /* 
     class ReferenceFrameDefinition : public osg::Referenced {
     public:
     ReferenceFrameDefinition() : Referenced() { }
     protected:
     virtual ~ReferenceFrameDefinition() { }
     
     public:
     virtual bool setName(const std::string & name) {
     return _name.conditionalSet(name);
     }
     public:
     virtual const std::string & getName() const {
     return _name.getRef();
     }
     protected:
     mutable orsa::Cache<std::string> _name;
     
     public:
     virtual orsa::Matrix localToGlobal(const orsa::Time &) const = 0;  
     virtual orsa::Matrix globalToLocal(const orsa::Time &) const = 0;
     };
  */
  
  /* 
     class ReferenceFrame {
     public:
     static ReferenceFrame * instance() {
     if (_instance == 0) {
     _instance = new ReferenceFrame;
     }
     return _instance;
     }
     public:
     static bool instanciated() {
     return (_instance != 0);
     }
     protected:
     ReferenceFrame();
     public:
     virtual ~ReferenceFrame() {
     _instance = 0;
     }
     protected:
     static ReferenceFrame * _instance;
     private:
     void _init();
     
     public:
     virtual bool insertDefinition(const ReferenceFrameDefinition * d) {
     definitionType::iterator it = _definition.begin();
     while (it != _definition.end()) {
     // check here if definition is already present in list...
     ++it;
     }
     _definition.push_back(d);
     return true;
     }
     protected:
     typedef std::list< osg::ref_ptr<const ReferenceFrameDefinition> > definitionType;
     definitionType _definition;
     };
  */
  
  /* 
     class ReferenceFrame : public osg::Referenced {
     public:
     ReferenceFrame() : Referenced() { }
     protected:
     virtual ~ReferenceFrame() { }
     public:
     virtual orsa::Matrix relativeToAbsolute(const orsa::Time &) const = 0;  
     virtual orsa::Matrix absoluteToRelative(const orsa::Time &) const = 0;
     };
  */
  
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
  
  
  /* 
     class AbsoluteReferenceFrame : public ReferenceFrame {
     public:
     ReferenceFrame * instance() {
     if (_instance == 0) {
     _instance = new AbsoluteReferenceFrame;
     }
     return _instance;
     }
     public:
     bool instanciated() {
     return (_instance != 0);
     }
     protected:
     AbsoluteReferenceFrame();
     public:
     virtual ~AbsoluteReferenceFrame() {
     _instance = 0;
     }
     protected:
     static ReferenceFrame * _instance;
     
     public:
     orsa::Matrix relativeToAbsolute(const orsa::Time &) const {
     return Matrix::identity();
     }
     public:
     orsa::Matrix absoluteToRelative(const orsa::Time &) const {
     return Matrix::identity();
     }
     };
  */
  
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
