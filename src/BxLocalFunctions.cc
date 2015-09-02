//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#include "globals.hh"
#include "BxLocalFunctions.hh"
#include <stdio.h>
#include <iostream>
#include "Randomize.hh"

using namespace std;

BxLocalFunctions::BxLocalFunctions() {
  fCenter = new G4ThreeVector(0.,0.,0.) ;
}

BxLocalFunctions::~BxLocalFunctions() {
  delete fCenter;
}

//It extracts a random point in the region between two concentrical spheres with radii fRadiusMin and fRadiusMax
G4ThreeVector* BxLocalFunctions::SDSphere() {
  G4double x = 0.;
  G4double y = 0.;
  G4double z = 0.;
  while(((x*x)+(y*y)+(z*z)) < (fRadiusMin*fRadiusMin) 
     || ((x*x)+(y*y)+(z*z)) > (fRadiusMax*fRadiusMax)) {
    x = 2*G4UniformRand()*fRadiusMax - fRadiusMax;
    y = 2*G4UniformRand()*fRadiusMax - fRadiusMax;
    z = 2*G4UniformRand()*fRadiusMax - fRadiusMax;

  }

  fPosition = new G4ThreeVector(x,y,z);
  *fPosition += *fCenter;
  return  fPosition ;
}
