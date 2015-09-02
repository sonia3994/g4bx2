// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4GeneralParticleSource.hh"
#include "G4SPSRandomGenerator.hh"
#include "BxReadParameters.hh"
#include "G4Electron.hh"
#include <fstream>
#include <iostream>
#include "BxLogger.hh"
#include "BxDataCollection.hh"
#include "BxOutputVertex.hh"
#include "BxGeneratorPositronium.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
using namespace std;


BxGeneratorPositronium::BxGeneratorPositronium(): BxVGenerator("BxGeneratorPositronium") {
  
  fTau = 3.12 ;
  fOrthoProb = 0.5 ;
  fEnergy = 100.*keV ;
  fSpectrum=false;
  isRead = false;
  fVolumeFlag = false ;
  G4ThreeVector zero(0., 0., 0.) ;
  fSPSAng = new G4SPSAngDistribution ;
  fSPSPos = new G4SPSPosDistribution ;


  G4SPSRandomGenerator *RndGen = new G4SPSRandomGenerator;
  fSPSPos->SetBiasRndm(RndGen);
  fSPSAng->SetBiasRndm(RndGen);
  SetRadius0(0.);
  SetCentreCoords(zero);    
  
  fParticleTable = G4ParticleTable::GetParticleTable();
  
  fElectron =  fParticleTable->FindParticle(11) ;
  fGamma    =  fParticleTable->FindParticle(22) ;

  fSPSAng->SetAngDistType("iso");    
  fSPSAng->SetPosDistribution(fSPSPos);

  fTheMessenger = new BxGeneratorPositroniumMessenger(this);
 
  BxLog(routine) << "Positronium Generator Constructed." << endlog;
  

}

//---------------------------------------------------------------------------//

BxGeneratorPositronium::~BxGeneratorPositronium()
{
  delete fTheMessenger;
  delete fSPSPos;
  delete fSPSAng;
  
}

//---------------------------------------------------------------------------//


void BxGeneratorPositronium::BxGeneratePrimaries(G4Event *event) {
  if(isFirstTime) {
    isFirstTime = false ;
  }

  if(fVolumeFlag)  {
    fSPSAng->SetVerbosity(0);
    SetParticlePosition(fSPSPos->GenerateOne());  
    fSPSAng->SetAngDistType("iso");    
    fSPSAng->SetPosDistribution(fSPSPos);
    SetParticleMomentumDirection(fSPSAng->GenerateOne());    
  }
  
  if(fSpectrum){ 
    if(!isRead){
      ifstream ff("../data/borex/c11_cumulativa.dat");
      G4double ee=0, cc=0;
      while(!ff.eof()){
	ff >> ee >> cc;
	if(ff.eof()) break;
	fShootEnergy.push_back(ee);
	fShootRandom.push_back(cc);
      }
      isRead=true;
    }  
    fEnergy = ShootEnergy()/MeV ; 
  }
  
  G4double mass = fElectron->GetPDGMass();
  G4double EneTot = mass + fEnergy;
  G4double pmom = std::sqrt(EneTot*EneTot - mass*mass);
  G4double px = pmom*fDirection.x();
  G4double py = pmom*fDirection.y();
  G4double pz = pmom*fDirection.z();
  
  G4PrimaryVertex*   vertexe   = new G4PrimaryVertex(fPosition,0.0*ns);
  G4PrimaryParticle* particlee = new G4PrimaryParticle(fElectron,px,py,pz);

  vertexe->SetPrimary( particlee );
  event->AddPrimaryVertex( vertexe );    
  
  BxOutputVertex::Get()->SetPosition(fPosition);  
  BxOutputVertex::Get()->SetDirection(fDirection);
  BxOutputVertex::Get()->SetTime(0.);
  BxOutputVertex::Get()->SetEnergy(fEnergy) ;

 
  G4double delay=0;
  
  if(CLHEP::RandFlat::shoot() < fOrthoProb) 
    delay = CLHEP::RandExponential::shoot(fTau);
  else 
    delay = 0;


  const G4double energyg = 0.511*MeV;

  SetParticleMomentumDirection(fSPSAng->GenerateOne());    

  G4double pxg = energyg*fDirection.x();
  G4double pyg = energyg*fDirection.y();
  G4double pzg = energyg*fDirection.z();

  G4PrimaryVertex*   vertexg   = new G4PrimaryVertex(fPosition,delay*ns);
  G4PrimaryParticle* gamma1 = new G4PrimaryParticle(fGamma,pxg,pyg,pzg);
  G4PrimaryParticle* gamma2 = new G4PrimaryParticle(fGamma,-pxg,-pyg,-pzg);
  
  vertexg->SetPrimary( gamma1 );
  vertexg->SetPrimary( gamma2 );
  event->AddPrimaryVertex( vertexg );    


  BxOutputVertex::Get()->SetEnergy(fEnergy/MeV);	
  BxOutputVertex::Get()->SetDId(2);
  BxOutputVertex::Get()->SetDPosition(fPosition);
  BxOutputVertex::Get()->SetDDirection(fDirection);
  BxOutputVertex::Get()->SetDTime(delay*ns);
  BxOutputVertex::Get()->SetDEnergy((fEnergy+2*energyg)/MeV);
  BxOutputVertex::Get()->SetDaughters();
	

}


G4double BxGeneratorPositronium::ShootEnergy() {
  G4double val = 0;
  while(val==0) {
    G4double _val = G4UniformRand()  ;
    if(_val > fShootRandom[0]) val = _val;
  }
  for(G4int i=0;i<G4int(fShootRandom.size());i++) {
    if(fShootRandom[i] >= val) {
      if(i == 0) return fShootEnergy[0];
      G4double deltaX = val - fShootRandom[i-1] ;
      G4double dy = fShootEnergy[i] - fShootEnergy[i-1] ;
      G4double dx = fShootRandom[i] - fShootRandom[i-1] ;
      
      return fShootEnergy[i-1] + dy/dx*deltaX;
    }
  }
  return 0;
}

