#ifndef BxLocalFunctions_h
#define BxLocalFunctions_h 1
#include "G4ThreeVector.hh"

#include "globals.hh"
#include <stdio.h>
#include <iostream>
#include <vector>


using namespace std;
/**This class extracts randomly a point inside a spherical region
 * which is contained between a minimum and maximum radius.
 */

class BxLocalFunctions {
   private:
    G4double fRadiusMin;
    G4double fRadiusMax;
    G4ThreeVector*        fPosition;
    G4ThreeVector*        fCenter;


   public:
    
    BxLocalFunctions() ;
    virtual ~BxLocalFunctions() ;
    
    void SetRadiusMin(G4double a)     { fRadiusMin = a; }
    void SetRadiusMax(G4double a)     { fRadiusMax = a; }
    ///Default value ofthe centre is the origin
    void SetCenter(G4ThreeVector *a)  { fCenter = a; }
   ///It extracts the point randomly in the region between two concentrical spheres with radii fRadiusMin and fRadiusMax
    G4ThreeVector* SDSphere();


};

#endif
