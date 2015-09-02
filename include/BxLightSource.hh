//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef BXLIGHTSOURCE_H
#define BXLIGHTSOURCE_H

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>

class   G4PhysicsOrderedFreeVector;
/** 
 * This class is responsible of the generation of optical-photons during the absorption and reemission of light
*/
 class	BxLightSource  {

public:
  BxLightSource();
  ~BxLightSource() { }

  enum SourceType {
    ScintillatorEmission,
    PCenergyTransfertoPPOEmission,
    PPOEmission,
    PCEmission,
    PCinDMPEmission
  };
  
  static BxLightSource* Get();

  G4double	  DefineEmissionEnergy(SourceType type);
  G4double	  DefineEmissionTime(SourceType type, G4int ParticleType);
  G4ThreeVector   DefineRandomPolarization(G4ThreeVector& Direction);
  G4ThreeVector   DefineRandomDirection3D();
  G4ThreeVector   DefineRandomDirectionSemi3D();

private:
  
  static BxLightSource *me;

  G4PhysicsOrderedFreeVector*	  ScintillationIntegral[4];
  G4double			  CIImax[4];
  
  G4double   IsAlpha ;
  std::vector<G4double>   fDecayTimeConstant ;
  std::vector<G4double>   fDecayWeight ;
  G4double tppo;
  G4double tpc;
  G4double tpc_to_ppo;
};

#endif
