// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4SystemOfUnits.hh"
#include "BxLogger.hh"
#include "BxDataCollection.hh"
#include "BxOutputVertex.hh"
#include "BxGeneratorMultiEvent.hh"
#include "BxGeneratorMultiEventMessenger.hh"
//---------------------------------------------------------------------------//


BxGeneratorMultiEvent::BxGeneratorMultiEvent(): BxVGenerator("BxGeneratorMultiEvent") {

  fRead  = false ;
  fVolumeFlag = false ;
  fSPSAng = new G4SPSAngDistribution ;
  fSPSPos = new G4SPSPosDistribution ;
  fNVertex = 0;
  
  G4SPSRandomGenerator *RndGen = new G4SPSRandomGenerator;

  fSPSPos->SetBiasRndm(RndGen);
  fSPSAng->SetBiasRndm(RndGen);
  fSPSAng->SetAngDistType("iso");    

  SetRadius0(0.);
  SetCentreCoords(G4ThreeVector(0.,0.,0.));    
  fPosition = G4ThreeVector(0.,0.,0.) ;
  
  fParticleTable = G4ParticleTable::GetParticleTable();  ;
  fMessenger = new BxGeneratorMultiEventMessenger(this);
  
  
  
  
  
    BxLog(routine) << "G4MultiEvent Constructed." << endlog;
}
//---------------------------------------------------------------------------//

BxGeneratorMultiEvent::~BxGeneratorMultiEvent()
{
  delete fMessenger;
  delete fSPSPos;
  delete fSPSAng;
}

//---------------------------------------------------------------------------//


void BxGeneratorMultiEvent::BxGeneratePrimaries(G4Event *event) {
  if(!fRead) {
    fRead = true ;
    BxLog(routine)<< "Generated particles in each vertex:" << endlog;
    
    if(!int(fListOfParticles.size())) {
      BxLog(error) << "Number of selected particles equals to 0"<<endlog;
      BxLog(fatal) << endlog;  
    }

    // Branching ratio normalizer
    G4int Id = 0;
    G4double BRNorm = 0;
    for(G4int i=0; i<G4int(fListOfParticles.size()); i++) { 
      if(fListOfParticles[i].Id > Id) {
	Id = fListOfParticles[i].Id;
	BRNorm += fListOfParticles[i].BRTOT;
      }   
    }
        
    // Build the cumulative probability ditribution
    Id = 0;    
    G4double BRCum = 0;
    for(G4int i=0; i<G4int(fListOfParticles.size()); i++) {
      G4int pdg =  fListOfParticles[i].PDG;
      BxLog(routine) << i << " " 
                     << fParticleTable->FindParticle(pdg)->GetParticleName() << " "
                     << fListOfParticles[i].Energy << " MeV"
  		     << endlog;
      if(fListOfParticles[i].Id > Id) {
        fNVertex++; 
	Id = fListOfParticles[i].Id;
	BRCum += fListOfParticles[i].BRTOT;
	fBR.push_back(BRCum/BRNorm);
      }
    }
    fSPSAng->SetPosDistribution(fSPSPos);

  }


  // Choose the primary event
  G4int source = 0;
  G4double rand = G4UniformRand() ;  
  for(G4int i=0;i<fNVertex;i++) {
    if(rand < fBR[i]) {
      source = i ;
      break;
    }
  }
  
  
  // Fish the primary position
  if(fVolumeFlag)  fPosition = fSPSPos->GenerateOne();

  for(G4int i = 0; i <G4int(fListOfParticles.size()); i++) {
    if(fListOfParticles[i].Id != source + 1) continue ;
    if(G4UniformRand() > fListOfParticles[i].BR) continue;  
 
    const G4int pdg = fListOfParticles[i].PDG ;
    fParticle = fParticleTable->FindParticle(pdg);     

    const G4double mass = fParticle->GetPDGMass();
    const G4double energy = fListOfParticles[i].Energy*MeV  + mass;

    fDirection = fSPSAng->GenerateOne();
    
    G4double pmom = std::sqrt(energy*energy - mass*mass);
    G4double px = pmom*fDirection.x();
    G4double py = pmom*fDirection.y();
    G4double pz = pmom*fDirection.z();

    G4PrimaryVertex*   vertex   = new G4PrimaryVertex(fPosition,0);
    G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle,px,py,pz);

    vertex->SetPrimary( particle );
    event->AddPrimaryVertex( vertex );    
    
    BxOutputVertex::Get()->SetEnergy(BxOutputVertex::Get()->GetEnergy() + fListOfParticles[i].Energy*MeV);		
    BxOutputVertex::Get()->SetPosition(fPosition);		
   
    BxOutputVertex::Get()->SetDId(source+1);     
    BxOutputVertex::Get()->SetDPosition(fPosition);  
  //  BxOutputVertex::Get()->SetDDirection(pmom); 
    BxOutputVertex::Get()->SetDPDG(pdg);
    BxOutputVertex::Get()->SetDTime(0.);
    BxOutputVertex::Get()->SetDEnergy(fListOfParticles[i].Energy*MeV);		
    BxOutputVertex::Get()->SetDaughters();  
  }
  if(!event->GetNumberOfPrimaryVertex ()) {
    const G4int pdg =  12;
    fParticle = fParticleTable->FindParticle(pdg);     

    const G4double mass = fParticle->GetPDGMass();
    const G4double energy = mass;

    fDirection = fSPSAng->GenerateOne();
    
    G4double pmom = std::sqrt(energy*energy - mass*mass);
    G4double px = pmom*fDirection.x();
    G4double py = pmom*fDirection.y();
    G4double pz = pmom*fDirection.z();

    G4PrimaryVertex*   vertex   = new G4PrimaryVertex(fPosition,0);
    G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle,px,py,pz);
    vertex->SetPrimary( particle );
    event->AddPrimaryVertex( vertex );    
  }
}

void BxGeneratorMultiEvent::SetfParticles( G4ThreeVector fPart) {
  fPDG.push_back(G4int(fPart.x()));        // PDG Code number 
  fPDGEnergy.push_back(fPart.y()*MeV);   // Energy in MeV
  fPDGBR.push_back(fPart.z());           // Branching Ratio
}   
