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
#include "G4SPSRandomGenerator.hh"
#include "G4UImanager.hh"
#include "BxLogger.hh"
#include "BxReadParameters.hh"
#include "BxOutputVertex.hh"
#include "BxDataCollection.hh"
#include "G4RadioactiveDecay.hh"
#include "G4IonTable.hh"
#include "BxGeneratorRDMDecayGun.hh"
#include "BxGeneratorRDMDecayGunMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
//---------------------------------------------------------------------------//


BxGeneratorRDMDecayGun::BxGeneratorRDMDecayGun(): BxVGenerator("BxGeneratorRDMDecayGun") {
  fVolumeFlag = false ;
  G4ThreeVector zero(0., 0., 0.) ;
  fSPSAng = new G4SPSAngDistribution ;
  fSPSPos = new G4SPSPosDistribution ;
 fRndGen = new G4SPSRandomGenerator;
  fSPSPos->SetBiasRndm(fRndGen);
  fSPSAng->SetBiasRndm(fRndGen);
  SetRadius0(0.);
  SetCentreCoords(zero);    
  fParticleGun  = new G4ParticleGun ;
  IsFirst = true;


  BxOutputVertex::Get()->SetPDG(fParticleGun->GetParticleDefinition()->GetPDGEncoding());
  BxOutputVertex::Get()->SetEnergy(0.);
  BxOutputVertex::Get()->SetDirection(G4ThreeVector(0.0, 0.0, 1.0));

  fParticleGun->SetParticleEnergy(0.*eV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0));
  fTheMessenger = new BxGeneratorRDMDecayGunMessenger(this);
  if(!fParticleGun) {
    BxLog(error) << "Could not allocate G4ParticleGun! Out of memory?"<<endlog;
    BxLog(fatal) << endlog;
  }
  BxLog(routine) << "G4ParticleGun Constructed." << endlog;
}

//---------------------------------------------------------------------------//

void BxGeneratorRDMDecayGun::SetNucleus (BxGeneratorRDMNucleus theIon1)
{

  fA = theIon1.GetA();
  fZ = theIon1.GetZ();
  fE = theIon1.GetE();
}


BxGeneratorRDMDecayGun::~BxGeneratorRDMDecayGun()
{
  delete fTheMessenger;
  delete fParticleGun;
  delete fSPSPos;
  delete fSPSAng;
  delete fRndGen; 
}

//---------------------------------------------------------------------------//
 

void BxGeneratorRDMDecayGun::BxGeneratePrimaries(G4Event *event) {

  if(IsFirst) {
 //   G4IonTable *theIonTable = new G4IonTable ;
//    theIonTable = G4ParticleTable::GetParticleTable()->GetIonTable();
//G4ParticleDefinition *aIon=G4IonTable::GetIonTable()->GetIon(9,18,0*keV);
//G4ParticleDefinition *aIon=G4IonTable::GetIonTable()->GetIon(84, 212, 0*keV);   
G4ParticleDefinition *aIon=G4IonTable::GetIonTable()->GetIon(fZ, fA, fE);   
 //G4ParticleDefinition *aIon = theIonTable->GetIon (fZ, fA, fE);
    fParticleGun->SetParticleDefinition(aIon);
    BxLog(routine) << "**************************" << endlog ;
    BxLog(routine) << "  Isotope: " << aIon->GetParticleName() << endlog ;
    BxLog(routine) << "**************************" << endlog ;
    IsFirst = false ;
  }
  
  if(fVolumeFlag)  {
    fSPSAng->SetVerbosity(0);
    SetParticlePosition(fSPSPos->GenerateOne());  
    fSPSAng->SetAngDistType("iso");    
    fSPSAng->SetPosDistribution(fSPSPos);    
  }
  BxOutputVertex::Get()->SetPosition(fParticleGun->GetParticlePosition());  
  BxOutputVertex::Get()->SetDirection(fParticleGun->GetParticleMomentumDirection()); 
  BxOutputVertex::Get()->SetPDG(fParticleGun->GetParticleDefinition()->GetPDGEncoding());
  BxOutputVertex::Get()->SetTime(0.);
  BxOutputVertex::Get()->SetEnergy(fParticleGun->GetParticleEnergy()/MeV);		

  BxDataCollection::Get()->SetEnergy(fParticleGun->GetParticleEnergy());
  BxDataCollection::Get()->SetPosition(fParticleGun->GetParticlePosition());
  BxDataCollection::Get()->SetDirection(fParticleGun->GetParticleMomentumDirection());
  //BxDataCollection::Get()->SetShootDist(sqrt(fParticleGun->GetParticlePosition().x()*fParticleGun->GetParticlePosition().x()+
 //                                      fParticleGun->GetParticlePosition().y()*fParticleGun->GetParticlePosition().y()));
  fParticleGun->SetParticleTime (0.0*ns);
  fParticleGun->GeneratePrimaryVertex(event);

}


