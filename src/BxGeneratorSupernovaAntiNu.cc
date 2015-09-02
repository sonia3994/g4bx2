#include "BxEventAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4Neutron.hh"
#include "BxLogger.hh"
#include "BxDataCollection.hh"
#include "BxOutputVertex.hh"
#include "BxGeneratorSupernovaAntiNu.hh"
#include "G4SystemOfUnits.hh"
#include "BxGeneratorSupernovaAntiNuMessenger.hh"
//---------------------------------------------------------------------------//
/**Please use the stacking action command in macro
 * file /bx/stack/select to postpone the neutron capture
 * event. This is necessary for bx-elec simultion to
 * process the neutron capture as the new event.
 */
///Antineutrino + Proton reaction simulation


BxGeneratorSupernovaAntiNu::BxGeneratorSupernovaAntiNu(): BxVGenerator("BxGeneratorSupernovaAntiNu") {


  isFirstTime = true ; 
  //fNeutrinoType = -1;
  


  /*fVolumeFlag = false ;
  G4ThreeVector zero(0., 0., 0.) ;
  fSPSAng = new G4SPSAngDistribution ;
  fSPSPos = new G4SPSPosDistribution ;
  G4SPSRandomGenerator *RndGen = new G4SPSRandomGenerator;
  fSPSPos->SetBiasRndm(RndGen);
  fSPSAng->SetBiasRndm(RndGen);
  SetRadius0(0.);
  SetCentreCoords(zero);    
  fPosition = G4ThreeVector(0.,0.,0.) ;
  fParticleTable = G4ParticleTable::GetParticleTable(); */   

  fTheMessenger = new BxGeneratorSupernovaAntiNuMessenger(this);
 
}
//---------------------------------------------------------------------------//


BxGeneratorSupernovaAntiNu::~BxGeneratorSupernovaAntiNu()
{
  delete fTheMessenger;
  //delete fSPSPos;
  //delete fSPSAng;
  
}

//---------------------------------------------------------------------------//

void BxGeneratorAntiNeutrino::BxGeneratePrimaries(G4Event *event) {


	if(isFirstTime) {
		if(fNeutrinoType < 0) {
      	    BxLog(error) << "Set the antineutrino source type"<<endlog;
      	    BxLog(fatal) << endlog;
    	}

    	  BxLog(routine) << "AntiNeutrino source " << fNeutrinoType << endlog ;

    	  isFirstTime = false ;
	}

	if ((event->GetEventID())%2 == 0) {

		fParticle = fParticleTable->FindParticle(-11); // -11 for the positron
		G4double anglePos = ShootAnglePositron();;
	}

}
