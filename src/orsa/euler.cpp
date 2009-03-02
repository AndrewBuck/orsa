#include <orsa/euler.h>

#include <orsa/print.h>
#include <orsa/util.h>

using namespace orsa;

bool orsa::Euler(orsa::Vector       & omegaDot,
		 const orsa::Vector & omega,
		 const orsa::Matrix & I,
		 const orsa::Vector & T,
		 const bool usePrincipalAxisRotationTransformation) {
  
  // test
  /* 
     #warning test
     omegaDot = orsa::Vector(0,0,0);
     return true;
  */
  
  if (usePrincipalAxisRotationTransformation) {
    
    orsa::Matrix Ip, genericToPrincipal, principalToGeneric;
    
    orsa::principalAxis(genericToPrincipal,
			Ip,
			I);
    
    if (!orsa::Matrix::invert(genericToPrincipal, principalToGeneric)) {
      ORSA_DEBUG("problems...");
      return false;
    }
    
    /* 
       print(I);
       print(Ip);
       print(genericToPrincipal*Ip*principalToGeneric);
       print(principalToGeneric*I *genericToPrincipal);
       
       print(genericToPrincipal);
       print(principalToGeneric);
       print(principalToGeneric*genericToPrincipal);
       print(genericToPrincipal*principalToGeneric);
    */
    
    orsa::Matrix inv_Ip;
    
    if (orsa::Matrix::invert(Ip, inv_Ip)) {
      
      /* 
	 print(Ip);
	 print(inv_Ip);
	 print(Ip*inv_Ip);
	 print(inv_Ip*Ip);
      */
      
      const orsa::Vector Lp = Ip * genericToPrincipal * omega;
      
      omegaDot = principalToGeneric*inv_Ip*(genericToPrincipal*T - orsa::externalProduct(genericToPrincipal*omega,Lp));
      
      ORSA_DEBUG("omegaDot...");
      print(omegaDot);
      
      return true;
      
    } else {
      
      ORSA_DEBUG("problems...");
      
      return false;
    }
    
  } else {
    
    orsa::Matrix inv_I;
    
    if (orsa::Matrix::invert(I, inv_I)) {
      
      const orsa::Vector L = I * omega;
      
      omegaDot = inv_I*(T-orsa::externalProduct(omega,L));
      
      ORSA_DEBUG("omegaDot...");
      print(omegaDot);
      
      return true;
      
    } else {
      
      ORSA_DEBUG("problems...");
      
      return false;
    }
    
  }
  
}
