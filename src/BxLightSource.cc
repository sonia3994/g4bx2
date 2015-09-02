//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#include "BxLightSource.hh"
#include "BxLogger.hh"
#include "G4PhysicsTable.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Random/RandExponential.h"
#include "Randomize.hh"
#include "BxPropertyCollection.hh"
#include "BxReadParameters.hh"
#include <vector>

namespace CLHEP {} 
using namespace CLHEP;


BxLightSource* BxLightSource::me = 0;


BxLightSource::BxLightSource() {

  G4PhysicsTable  *thePhysicsTable;
  G4Material      *theMaterial;
  G4int		   theMaterialIndex;

  thePhysicsTable = BxPropertyCollection::Get()->GetPhysicsTable();
  if(!thePhysicsTable)	{
    G4cout << "Couldn't load Physics Table for Light Source" << G4endl;
    exit(-1);
  }

  theMaterial = G4Material::GetMaterial("PPO");
  theMaterialIndex = theMaterial->GetIndex();
  ScintillationIntegral[0] =    (G4PhysicsOrderedFreeVector*)((*thePhysicsTable)(theMaterialIndex));
  CIImax[0] = ScintillationIntegral[0]->GetMaxValue();


  theMaterial = G4Material::GetMaterial("PC");
  theMaterialIndex = theMaterial->GetIndex();
  ScintillationIntegral[1] =     (G4PhysicsOrderedFreeVector*)((*thePhysicsTable)(theMaterialIndex));
  CIImax[1] = ScintillationIntegral[1]->GetMaxValue();

  tppo	=	BxReadParameters::Get()->GetTimePPOEmission();
  tpc	=	BxReadParameters::Get()->GetTimePCEmission();
  tpc_to_ppo =	BxReadParameters::Get()->GetTimePCtoPPOTransfer();
  
  /*OLD CODE
  tppo = 1.1;  // was tppo = 1.6
  tpc  = 28.0;
  tpc_to_ppo = 2.7; // was tpc_to_ppo = 3.574
*/
  }

// Singleton
BxLightSource* BxLightSource::Get() {
  if (!me) 
    me = new BxLightSource();
       
  return me;
}

G4double BxLightSource::DefineEmissionEnergy(SourceType type)
{
  G4double theCIIvalue;
  G4double theEnergy=0.0;

  switch (type) {
    case ScintillatorEmission: 
    case PCenergyTransfertoPPOEmission:
    case PPOEmission:
      	theCIIvalue = G4UniformRand()*CIImax[0];
    	theEnergy = ScintillationIntegral[0]->GetEnergy(theCIIvalue);
	break;
    case PCEmission:
    case PCinDMPEmission:
	theCIIvalue = G4UniformRand()*CIImax[1];
    	theEnergy = ScintillationIntegral[1]->GetEnergy(theCIIvalue);
	break;
  }

  if(theEnergy == 0.0)
    G4cout << "Error in BxLightSource.cc, Energy=0"<<G4endl;
  return theEnergy;
}

G4double BxLightSource::DefineEmissionTime(SourceType type, G4int ParticleType)
{
  IsAlpha = ParticleType;
  G4int NumOfExpo;
  
  if(IsAlpha) {
    fDecayTimeConstant = BxReadParameters::Get()->GetAlphaDecayTimeConstant();
    fDecayWeight       = BxReadParameters::Get()->GetAlphaDecayWeight();
  } else {
    fDecayTimeConstant = BxReadParameters::Get()->GetBetaDecayTimeConstant();
    fDecayWeight       = BxReadParameters::Get()->GetBetaDecayWeight();
  }

  if (fDecayTimeConstant.size()<fDecayWeight.size())
	    	NumOfExpo=fDecayTimeConstant.size();
    else
	    	NumOfExpo=fDecayWeight.size();
  
  G4double normalization=0;
  for (G4int i=0; i<NumOfExpo; i++)
	  normalization+=fDecayWeight[i];

  if (abs(normalization-1.)>4*std::numeric_limits<G4double>::epsilon())
	  BxLog(fatal) << "ERROR! Normalization of scintillation decay exponentials is not 1!" << endlog;

  G4double delay = 0.0;
  G4double rand = G4UniformRand();
  G4double sumq = 0;	
       
  switch (type) {
    
    case ScintillatorEmission:	
      
      sumq = 0;
      for(G4int i=0; i<NumOfExpo; i++) {
	if ( (sumq < rand) && (rand < (sumq+fDecayWeight[i])) ) 
	               delay = RandExponential::shoot(fDecayTimeConstant[i]);
	
	sumq += fDecayWeight[i];
      }

      break;

    case PPOEmission:	
      delay = RandExponential::shoot(tppo);
      break;
    
    case PCEmission:
      delay = RandExponential::shoot(tpc);
      break;

    case PCenergyTransfertoPPOEmission:
      delay = RandExponential::shoot(tpc_to_ppo);
      break;

    case PCinDMPEmission:
      delay = RandExponential::shoot(tpc/10.0);
      break;
  }
  return delay;
}


G4ThreeVector BxLightSource::DefineRandomPolarization(G4ThreeVector&   PhotonDirection) {
  G4ThreeVector PhotonPolarization,
  Perp;
  G4double	Phi,
  SinPhi,
  CosPhi;

  PhotonPolarization = PhotonDirection.orthogonal();

  Perp = PhotonDirection.cross(PhotonPolarization);

  Phi = 2*M_PI*G4UniformRand();
  SinPhi = sin(Phi);
  CosPhi = cos(Phi);

  PhotonPolarization = CosPhi * PhotonPolarization + SinPhi * Perp;

  PhotonPolarization = PhotonPolarization.unit();

  return PhotonPolarization;
}

G4ThreeVector	BxLightSource::DefineRandomDirection3D() {
  G4ThreeVector PhotonDirection;

  G4double	CosTheta,
  SinTheta,
  Phi,
  SinPhi,
  CosPhi;

  CosTheta = 1. - 2. * G4UniformRand();
  SinTheta = sqrt(1.0 - (CosTheta * CosTheta));

  Phi = 2*M_PI*G4UniformRand();
  SinPhi = sin(Phi);
  CosPhi = cos(Phi);

  PhotonDirection.setX(SinTheta * CosPhi);
  PhotonDirection.setY(SinTheta * SinPhi);
  PhotonDirection.setZ(CosTheta);

  return  PhotonDirection;
}

G4ThreeVector	BxLightSource::DefineRandomDirectionSemi3D() {
  G4ThreeVector PhotonDirection;

  G4double	CosTheta,
  SinTheta,
  Phi,
  SinPhi,
  CosPhi;

  CosTheta = 1. - G4UniformRand();
  SinTheta = sqrt(1.0 - (CosTheta * CosTheta));

  Phi = 2*M_PI*G4UniformRand();
  SinPhi = sin(Phi);
  CosPhi = cos(Phi);

  PhotonDirection.setX(SinTheta * CosPhi);
  PhotonDirection.setY(SinTheta * SinPhi);
  PhotonDirection.setZ(CosTheta);

  return  PhotonDirection;
}
