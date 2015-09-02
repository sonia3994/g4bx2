/** ////////////////////////////////////////////////////////////////////////
* // Optical Photon Absorption and Reemission Class Implementation
* ////////////////////////////////////////////////////////////////////////
* File:        BxOpAbsorptionReemission.cc
* Description: Discrete Process -- Absorption and reemission of Optical Photons  
* Version:     1.0
* Created:     1999-07-2
* Author:      MOSCOW "KURCHATOV'S INSTITUTE"
* Revised by A. Caminata and S. Marcocci, Sept. 2014
* ////////////////////////////////////////////////////////////////////////
* //        THE FOLLOWING PHYSICS OF LIGHT PROPAGATION IS EMPLOYED:
* ////////////////////////////////////////////////////////////////////////
*        I.  For PC+1.5g/l PPO (Scintillator):
*        
*        a) PPO Absorption:
*        1. The reemission probability for the case of
*        PPO Absorption equals to 0.82.
*        2. Scintillator decay time equals to 1.6 ns
*        for the process of direct PPO reemission. 
*        3. The energy specterun of reemited photons
*        equals to PPO emission spectrum.
*
*        b) PC Absorption:
*        1. The reemission probability equals to 0.82
*        for the process of direct PC excitement with
*        energy transfer to PPO molecules.
*        2. Scintillator decay time equals to 3.6 ns
*        for the process of direct PC excitement with
*        energy transfer to PPO molecules.
*        3. The energy specterun of reemited photons
*        equals to PPO emission spectrum.
*
*        II. For PC+DMP (DMP-buffer):
*
*        a) DMP Absorption:
*        In the case of DMP Absorption there is no
*        light reemission.
*
*        b) PC Absorption:
*        It is necessary to employ the following
*        model for DMP-buffer:
*        1. The reemission probability for the case of
*           PC Absorption equals to 0.40/10.
*        2. The reemission PC time equals to 28/10 ns
*        3. The energy spectrum of reemited photons
*        equals to PC emission spectrum
*        The reference person for this model is Barbara
*        Caccianiga (who refers to Genova PC+DMP studies)
*
*
*        III.Interaction of light with PC (for both Scintillator
*            and DMP-buffer cases)
*
*        a) The Railey scattering on PC starts after 310 nm
*        b) Absorption in PC exists below 310 nm
*        c) The reemission probability for the case of
*           PC Absorption in DMP-buffer equals to 0.40/10
*	  d) The reemission probability for the case of
*           PC absorption in scintillator equals to 0.82
*
*/

//////////////////////////////////////////////////////////////////////////
#include "G4ios.hh"
#include <iostream>
#include "BxOpAbsorptionReemission.hh"
#include "BxPropertyCollection.hh"
#include "BxDataCollection.hh"
#include "BxReadParameters.hh"
#include "BxLightSource.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicsTable.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "BxIO.hh"
#include "BxDetectorConstruction.hh"

namespace CLHEP {} 
using namespace CLHEP;

using namespace std;


BxOpAbsorptionReemission::BxOpAbsorptionReemission( const G4String& processName):G4VDiscreteProcess(processName) {
  
					  
  if (verboseLevel>0) G4cout << GetProcessName() << " is created " << G4endl;
  
  fMaterialIndexPPO  = G4Material::GetMaterial("PPO")->GetIndex();
  fMaterialIndexPC   = G4Material::GetMaterial("PC")->GetIndex();
  fMaterialIndexDMP  = G4Material::GetMaterial("DMP")->GetIndex();
  
  fScintillatorIndex = G4Material::GetMaterial("Scintillator")->GetIndex();
  fBufferIndex       = G4Material::GetMaterial("DMPbuffer")->GetIndex();
  
  fAbsorber = 0;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->size();

  G4cout << "Number of materials = " << numOfMaterials << G4endl;

  for (G4int i=0 ; i < numOfMaterials; i++)  {

    G4Material* theMaterial = (*theMaterialTable)[i];

    G4MaterialPropertiesTable*   theMaterialPropertyTable = theMaterial->GetMaterialPropertiesTable();
    
    if (theMaterialPropertyTable) {
      AbsorptionLengthVector[i] =  theMaterialPropertyTable->GetProperty("ABSLENGTH");
    } else	{
      G4cout << "No MaterialPropertyTable for " << theMaterial->GetName()<< G4endl;
      exit(1);
    }
  }

  G4PhysicsTable* thePhysicsTable;
  


  thePhysicsTable = BxPropertyCollection::Get()->GetPhysicsTable();
  if(!thePhysicsTable)	{
    G4cout << "Couldn't load Physics Table" << G4endl;
    exit(-1);
  }

  G4MaterialPropertiesTable*   theMaterialPropertyTable =  G4Material::GetMaterial("PPO")->GetMaterialPropertiesTable();

  if (theMaterialPropertyTable)
    ReemissionProbabilityVector1 = theMaterialPropertyTable->GetProperty("REEMISPROB");
  else	{
    G4cout << "No MaterialPropertyTable for PPO" << G4endl;
    exit(1);
  }

  theMaterialPropertyTable  = G4Material::GetMaterial("PC")->GetMaterialPropertiesTable();

  if (theMaterialPropertyTable ) 
    ReemissionProbabilityVector2 = theMaterialPropertyTable ->GetProperty("REEMISPROB");
  else	{
    G4cout << "No MaterialPropertyTable for PC" << G4endl;
    exit(1);
  }

}
BxOpAbsorptionReemission::~BxOpAbsorptionReemission() { }


//-----------------------------------------------------------------------------------
G4VParticleChange* BxOpAbsorptionReemission::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) {
  G4bool AbsIsGot; 
  enum  RainbowProcess {OrdinaryModeCase, DMPbufferCase, ScintillatorCase};

  RainbowProcess AbsReemMode;

  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle* aDynamicParticle = aTrack.GetDynamicParticle();

  G4double ReemissionProbability1 = 0.0;
  //G4double ReemissionProbability2 = 0.0;

  const G4Material* theMaterial = aTrack.GetMaterial();

  // To chose the medium processes
  AbsReemMode=OrdinaryModeCase;

  if(theMaterial->GetIndex () == fScintillatorIndex)  AbsReemMode  = ScintillatorCase;

  if(theMaterial->GetIndex()  == fBufferIndex)  AbsReemMode  = DMPbufferCase;
  switch(AbsReemMode) {
    case OrdinaryModeCase: DoAbsorption();
    break;

    case DMPbufferCase: {
      G4double thePhotonMomentum = aDynamicParticle->GetTotalMomentum();
      G4double thePhotonWave = (h_Planck*c_light)/thePhotonMomentum;

      if(ReemissionProbabilityVector2)	{
	//Check for Photon momentum range
	if (thePhotonMomentum < ReemissionProbabilityVector2-> GetMinLowEdgeEnergy() ||
		    thePhotonMomentum >ReemissionProbabilityVector2->GetMaxLowEdgeEnergy() ) {
	  G4cout << "Photon momentum out of range in reemission probability vector" << G4endl;
	  G4cout << "DMPbufferCase" << G4endl;
	  exit(1);
	}
	//ReemissionProbability2 = ReemissionProbabilityVector2->GetValue(thePhotonMomentum,AbsIsGot);
      } else {
	G4cout << "No reemission probability specified for PC" << G4endl;
	exit(1);
      }
      if (fAbsorber == 3)  DoAbsorption();	// 3 ---> DMP

      if (fAbsorber == 2) {	// 2 ---> PC

	if (thePhotonWave/nanometer >= 310) { 	//  Scattering on PC

	  //if(BxReadParameters::Get()->GetScattering()) DoScattering (aTrack);
          
	  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	} else  {        
	  //HERE I put Reemission by PC mixed with DMP (quenching = 0.40/10, time = 28ns/10,
	  //the Emission spectrum is the same like in pure PC)
	  //
	  if (G4BooleanRand(0.04)){
	  
	    DoReemission(BxLightSource::PCinDMPEmission,aTrack);
	  } else {
	    DoAbsorption();
	  }
	}
      }
    }
    break;

    case ScintillatorCase: {
      G4double thePhotonMomentum = aDynamicParticle->GetTotalMomentum();
      G4double thePhotonWave = (h_Planck*c_light)/thePhotonMomentum;
      if(ReemissionProbabilityVector1)	{
	//Check for Photon momentum range
	if (thePhotonMomentum <
	ReemissionProbabilityVector1-> GetMinLowEdgeEnergy() ||
		    thePhotonMomentum > ReemissionProbabilityVector1->GetMaxLowEdgeEnergy() ) {
	  exit(1);
	}
	ReemissionProbability1 = ReemissionProbabilityVector1->GetValue(thePhotonMomentum,AbsIsGot);
      }  else {
	      G4cout << "No reemission probability specified for PPO" << G4endl;
	      exit(1);
      }
      if(ReemissionProbabilityVector2)	{
	      //Check for Photon momentum range
	      if (thePhotonMomentum < ReemissionProbabilityVector2->GetMinLowEdgeEnergy() ||
			      thePhotonMomentum > ReemissionProbabilityVector2->GetMaxLowEdgeEnergy() ) {
		      G4cout << "Photon momentum \
			      out of range in \
			      reemission \
			      probability  \
			      vector" << G4endl;
		      exit(1);
	      }

	//ReemissionProbability2 = ReemissionProbabilityVector2->GetValue(thePhotonMomentum,AbsIsGot);
      } else {
	G4cout << "No reemission probability specified \
	for PC" << G4endl;
	exit(1);
      }
      
      // ****************************************************************************
      if (fAbsorber == 1) {	// 1 ---> PPO rremission
	// t=1.6 ns, quenching=0.82, PPO spectrum
	if ( G4BooleanRand(ReemissionProbability1) ) {

	  if(BxReadParameters::Get()->GetReemission()) DoReemission(BxLightSource::PPOEmission,aTrack);
	  //HERE I fill the flag for PPO reemissions in PhotonData
	} else {
          DoAbsorption();
	}
      }  else if (fAbsorber == 2) {	// 2 ---> PC
 	 if (thePhotonWave/nanometer >= 310) { 	

	  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	  	  
	} else {	// reemission by energy transfer from PC to PPO
	  // t=3.57 ns, quenching=0.82 , PPO spectrum
	  if (G4BooleanRand(ReemissionProbability1)){
	    DoReemission(BxLightSource::PCenergyTransfertoPPOEmission,aTrack);
	  } else {
	    DoAbsorption();
	  }
	}
      }
    }
    break;

  default :
  G4cout << "No choice for Medium in Absorption-Reemisiion module" << G4endl;
}

  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

}
//
// GetMeanFreePath
// ---------------
//
G4double BxOpAbsorptionReemission::GetMeanFreePath(const G4Track& aTrack,G4double , G4ForceCondition* ) {
 G4bool AbsIsGot;
  const G4DynamicParticle* aDynamicParticle = aTrack.GetDynamicParticle();
  const G4Material* theMaterial = aTrack.GetMaterial();

  G4double AbsorptionLength  = DBL_MAX;
  G4double AbsorptionLength1 = DBL_MAX;
  G4double AbsorptionLength2 = DBL_MAX;

  G4double thePhotonMomentum = aDynamicParticle->GetTotalMomentum();


  if(theMaterial->GetIndex() != fScintillatorIndex && theMaterial->GetIndex() != fBufferIndex )	{
    G4int  theMaterialIndex = theMaterial->GetIndex();
    if (AbsorptionLengthVector[theMaterialIndex])	{
      AbsorptionLength =AbsorptionLengthVector[theMaterialIndex]->GetValue (thePhotonMomentum,AbsIsGot);
    } else {
      if(verboseLevel > 0)	
        G4cout << "No Absorption length specified for "<< theMaterial->GetName() << G4endl;  
    }
  }

  if(theMaterial->GetIndex() == fScintillatorIndex)       {


    if (AbsorptionLengthVector[fMaterialIndexPPO])	{
      AbsorptionLength1 = AbsorptionLengthVector[fMaterialIndexPPO]->GetValue (thePhotonMomentum,AbsIsGot);
    } else {
      G4cout << "No Absorption length specified for PPO" <<  G4endl;
    }

    if (AbsorptionLengthVector[fMaterialIndexPC]) {
      AbsorptionLength2 = AbsorptionLengthVector[fMaterialIndexPC]->GetValue (thePhotonMomentum,AbsIsGot);
    } else {
      G4cout << "No Absorption length specified for PC" << G4endl;
    }

    G4double path1 = RandExponential::shoot(AbsorptionLength1);//PPO
    G4double path2 = RandExponential::shoot(AbsorptionLength2); //PC

    if (path1 < path2) {
      AbsorptionLength = AbsorptionLength1;
      fAbsorber = 1;
    } else {
      AbsorptionLength = AbsorptionLength2;
      fAbsorber = 2;
    }  
    return AbsorptionLength;

  }

  if(theMaterial->GetIndex() == fBufferIndex)       {

    if (AbsorptionLengthVector[fMaterialIndexDMP])	{
      AbsorptionLength1 = AbsorptionLengthVector[fMaterialIndexDMP]->GetValue (thePhotonMomentum,AbsIsGot);
    } else {
      G4cout << "No Absorption length specified for DMP" << G4endl;
    }

    if (AbsorptionLengthVector[fMaterialIndexPC])	{
      AbsorptionLength2 =AbsorptionLengthVector[fMaterialIndexPC]->GetValue (thePhotonMomentum,AbsIsGot);
    } else {
      G4cout << "No Absorption length specified for PC" << G4endl;
    }

    G4double path1 = RandExponential::shoot(AbsorptionLength1);//DMP
    G4double path2 = RandExponential::shoot(AbsorptionLength2);//PC 

    if (path1 < path2) {
      AbsorptionLength = AbsorptionLength1;
      fAbsorber = 3;
    } else {
      AbsorptionLength = AbsorptionLength2;
      fAbsorber = 2;
    }
    return AbsorptionLength;
  }      
  
   return AbsorptionLength;

}

void BxOpAbsorptionReemission::DoAbsorption(){ 
  aParticleChange.ProposeTrackStatus(fStopAndKill);

  if (verboseLevel>0)   G4cout << "Photon Absorbed!" << G4endl;
  
}

void BxOpAbsorptionReemission::DoReemission(BxLightSource::SourceType type,const G4Track& aTrack) {
  G4double	thePhotonEnergy = DBL_MAX*eV;
const G4DynamicParticle* aDynamicParticle=aTrack.GetDynamicParticle();
  if (verboseLevel>0) {
    G4cout << "Photon Absorpted and Reemitted!" << G4endl;
    if (verboseLevel>1) {
      G4cout << "Old Momentum Direction: "
      << aDynamicParticle->GetMomentumDirection() << G4endl;
      G4cout << "Old Polarization: "
      << aDynamicParticle->GetPolarization() << G4endl;
    }
  }


  G4ThreeVector NewMomentumDirection = BxLightSource::Get()->DefineRandomDirection3D();

  G4ThreeVector NewPolarization = BxLightSource::Get()->DefineRandomPolarization(NewMomentumDirection);

  aParticleChange.ProposePolarization(NewPolarization.unit());

  aParticleChange.ProposeMomentumDirection(NewMomentumDirection.unit());

  thePhotonEnergy = BxLightSource::Get()->DefineEmissionEnergy(type);
  thePhotonEnergy *= 1E6;
  thePhotonEnergy *= eV;

  aParticleChange.ProposeEnergy(thePhotonEnergy);

  //_______________HERE I TAKE INTO ACCOUNT THE TIME FOR REEMISSION________________
  //
  G4double delay = BxLightSource::Get()->DefineEmissionTime(type,0);

  G4double aTime = aTrack.GetGlobalTime();

  aTime = aTime + delay;

  aParticleChange.ProposeGlobalTime(aTime);

  if (verboseLevel>1) {
    G4cout << "New Polarization: " << NewPolarization << G4endl;
    G4cout << "Polarization Change: "
    << *(aParticleChange.GetPolarization()) << G4endl;
    G4cout << "New Momentum Direction: "
    << NewMomentumDirection << G4endl;
    G4cout << "Momentum Change: "
    << *(aParticleChange.GetMomentumDirection()) << G4endl;
  }
}

void BxOpAbsorptionReemission::DoScattering(const G4Track& aTrack)  {

  const G4DynamicParticle* aDynamicParticle = aTrack.GetDynamicParticle();

  if (verboseLevel>0) {
    G4cout << "Here is the Scattering Photon!" << G4endl;

    G4cout << "Old Momentum Direction: "
    << aDynamicParticle->GetMomentumDirection() << G4endl;
    G4cout << "Old Polarization: "
    << aDynamicParticle->GetPolarization() << G4endl;
  }



  // find polar angle w.r.t. old polarization vector

  G4double rand = G4UniformRand();

  G4double CosTheta = pow(rand, 1./3.);
  G4double SinTheta = sqrt(1.-CosTheta*CosTheta);

  if(G4UniformRand() < 0.5)CosTheta = -CosTheta;

  // find azimuthal angle w.r.t old polarization vector

  rand = G4UniformRand();

  G4double Phi = twopi*rand;
  G4double SinPhi = sin(Phi);
  G4double CosPhi = cos(Phi);

  G4double unit_x = SinTheta * CosPhi;
  G4double unit_y = SinTheta * SinPhi;
  G4double unit_z = CosTheta;

  G4ThreeVector NewPolarization (unit_x,unit_y,unit_z);

  G4ThreeVector OldPolarization = aDynamicParticle->GetPolarization();
  OldPolarization = OldPolarization.unit();

  NewPolarization.rotateUz(OldPolarization);
  NewPolarization = NewPolarization.unit();

  // -- new momentum direction is normal to the new polarization
  // vector (components below expressed in reference system where
  // new polarization vector is aligned with the z axis)

  G4ThreeVector NewMomentumDirection = OldPolarization - NewPolarization * CosTheta;

  if(G4UniformRand() < 0.5)NewMomentumDirection = -NewMomentumDirection;
  NewMomentumDirection = NewMomentumDirection.unit();

  aParticleChange.ProposePolarization(NewPolarization);

  aParticleChange.ProposeMomentumDirection(NewMomentumDirection);

  if (verboseLevel>0) {
    G4cout << "New Polarization: "
    << NewPolarization << G4endl;
    G4cout << "Polarization Change: "
    << *(aParticleChange.GetPolarization()) << G4endl;
    G4cout << "New Momentum Direction: "
    << NewMomentumDirection << G4endl;
    G4cout << "Momentum Change: "
    << *(aParticleChange.GetMomentumDirection()) << G4endl;
  }
}
