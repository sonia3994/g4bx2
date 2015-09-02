// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
* Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#include "HistoManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SPSRandomGenerator.hh"
#include "BxLogger.hh"
#include "G4SystemOfUnits.hh"
#include "BxDataCollection.hh"
#include "BxOutputVertex.hh"
#include "BxGeneratorG4Gun.hh"
#include "BxGeneratorG4GunMessenger.hh"
#include "BxLightSource.hh"
//---------------------------------------------------------------------------//


BxGeneratorG4Gun::BxGeneratorG4Gun(): BxVGenerator("BxGeneratorG4Gun") {
	fVolumeFlag = false ;
	fEnergyDistribution = false ;
	fDebugGenerator=true;
	G4ThreeVector zero(0., 0., 0.) ;
	fSPSAng = new G4SPSAngDistribution ;
	fSPSPos = new G4SPSPosDistribution ;
	fSPSEne = new G4SPSEneDistribution ;
	G4SPSRandomGenerator *RndGen = new G4SPSRandomGenerator;
	fSPSPos->SetBiasRndm(RndGen);
	fSPSAng->SetBiasRndm(RndGen);
	fSPSEne->SetBiasRndm(RndGen);
	SetRadius(0.);
	SetCentreCoords(zero);    
	fParticleGun= new G4ParticleGun;
	fTheMessenger = new BxGeneratorG4GunMessenger(this);

	if(!fParticleGun) {
		BxLog(error) << "Could not allocate G4ParticleGun! Out of memory?"<<endlog;
		BxLog(fatal) << endlog;
	}
	BxLog(routine) << "G4ParticleGun Constructed." << endlog;
}

//---------------------------------------------------------------------------//

BxGeneratorG4Gun::~BxGeneratorG4Gun()
{
	delete fTheMessenger;
	delete fParticleGun;
	delete fSPSPos;
	delete fSPSAng;
	delete fSPSEne;
}

//---------------------------------------------------------------------------//


void BxGeneratorG4Gun::BxGeneratePrimaries(G4Event *event) {
	
	if(fVolumeFlag)  {
		fSPSAng->SetVerbosity(0);
		SetParticlePosition(fSPSPos->GenerateOne()); 
		fSPSAng->SetAngDistType("iso");    
		fSPSAng->SetPosDistribution(fSPSPos);
		SetParticleMomentumDirection(fSPSAng->GenerateOne());    
	}
	if(fEnergyDistribution) {
		fSPSEne->SetVerbosity(3);
		SetParticleEnergy(fSPSEne->GenerateOne(fParticleGun->GetParticleDefinition()));
	}
	if((fParticleGun->GetParticleDefinition()->GetPDGEncoding())==50){
		G4ThreeVector OpDirection=fParticleGun->GetParticleMomentumDirection();
		fParticleGun->SetParticlePolarization(BxLightSource::Get()->DefineRandomPolarization(OpDirection));
	}
	fParticleGun->SetNumberOfParticles(1);
	fParticleGun->GeneratePrimaryVertex(event);

	//to be activated to debug particle generation

	if(fDebugGenerator){
		BxLog(debugging)<<"----------------------------"<<endlog;
		BxLog(debugging)<<"Primary particle has direction "<<fParticleGun->GetParticleMomentumDirection()<<endlog;
		BxLog(debugging)<<"Primary particle has polarization "<<fParticleGun->GetParticlePolarization()<<endlog;
		BxLog(debugging)<<"Primary particle has energy "<<fParticleGun->GetParticleEnergy()<<" MeV"<<endlog;
		BxLog(debugging)<<"Primary particle generated in "<<fParticleGun->GetParticlePosition()<<" mm"<<endlog;
		BxLog(debugging)<<"Particle PDG code = "<<fParticleGun->GetParticleDefinition()->GetPDGEncoding()<<endlog;
	}
	BxDataCollection::Get()->SetEnergy(fParticleGun->GetParticleEnergy());
	BxDataCollection::Get()->SetEID(event->GetEventID ());
	BxDataCollection::Get()->SetPosition(fParticleGun->GetParticlePosition());
	BxDataCollection::Get()->SetDirection(fParticleGun->GetParticleMomentumDirection());

//AC: why z not included??
//  BxDataCollection::Get()->SetShootDist(sqrt(fParticleGun->GetParticlePosition().x()*fParticleGun->GetParticlePosition().x()+fParticleGun->GetParticlePosition().y()*fParticleGun->GetParticlePosition().y()));


  BxOutputVertex::Get()->SetPosition(fParticleGun->GetParticlePosition());  
  BxOutputVertex::Get()->SetDirection(fParticleGun->GetParticleMomentumDirection()); 
  BxOutputVertex::Get()->SetPDG(fParticleGun->GetParticleDefinition()->GetPDGEncoding());
  BxOutputVertex::Get()->SetTime(0.);
  BxOutputVertex::Get()->SetEnergy(fParticleGun->GetParticleEnergy()/MeV);		
}

