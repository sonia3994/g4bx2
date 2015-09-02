//Created by I. Machulin
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#include "BxSteppingAction.hh"
#include "BxLogger.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "BxReadParameters.hh"
#include "BxDataCollection.hh"
#include "BxOutputVertex.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

/**
*	This part is developed to count for photon reflections from tyvec 
*	or steel surfaces after the consultance with Peter Gruimpheld
*	(Geant4 optic photon expert).
*	Contact: machulin@in2p3.fr
*/

extern   G4int NbOfPhReflections;

BxSteppingAction::BxSteppingAction() {
}

void BxSteppingAction::UserSteppingAction(const G4Step* theStep) {
	if(theStep->GetTrack()->GetDefinition()->GetPDGEncoding() != 50) {
		if(BxLogger::GetSeverity() == BxLogger::development) {
			BxLog(development)
				<< " " <<theStep->GetTrack()->GetDefinition()->GetParticleName() << " "
				<< " " <<theStep->GetTrack()->GetDefinition()->GetPDGEncoding() << " "
				<<  " E: " <<theStep->GetTrack()->GetKineticEnergy()/keV<< " keV; "
				<<  " " <<theStep->GetTrack()->GetVolume()->GetName()<< " "
				<<  " step " <<   theStep->GetStepLength()/um << " "   
				<<  " ID: " <<theStep->GetTrack()->GetTrackID()
				<<  " "    <<theStep->GetTrack()->GetParentID()
				<<  " gtime: " <<theStep->GetTrack()->GetGlobalTime() << " "
				<<  " ltime: " <<theStep->GetTrack()->GetLocalTime()/ns  << " "
				<<  " steps: " <<theStep->GetTrack()->GetCurrentStepNumber() << " "
				<<theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName ()
				<<endlog ;
		}
	}

	if(theStep->GetTrack()->GetDefinition()->GetPDGEncoding() == 50) { 
		// Kill optical photons outside range definition
		if(theStep->GetTrack()->GetKineticEnergy()  < 0.5*eV ||
				theStep->GetTrack()->GetKineticEnergy()  > 17*eV) 
			theStep->GetTrack()->SetTrackStatus(fStopAndKill);         

		// Kill trapped photons
		if(theStep->GetTrack()->GetCurrentStepNumber() > 200)
			theStep->GetTrack()->SetTrackStatus(fStopAndKill);


	}


	if(BxReadParameters::Get()->IsRDMDecay()) {

		if(theStep->GetTrack()->GetParentID() <= 0 && theStep->GetTrack()->GetDefinition()->GetAtomicNumber () > 2 ) 
			theStep->GetTrack()->GetDefinition()->SetPDGLifeTime(BxReadParameters::Get()->GetRealPDGMeanLife());

		if(theStep->GetTrack()->GetDefinition()->GetAtomicNumber () > 2 && !BxReadParameters::Get()->IsRDMChain()) 
theStep->GetTrack()->SetTrackStatus(fStopAndKill);         
	}


	if(  BxOutputVertex::Get()->PostponeFlag() 
			&& theStep->GetPostStepPoint()->GetPosition().r() < BxOutputVertex::Get()->GetKillRadius()) 
		BxOutputVertex::Get()->SetIsPostponed(true);

if (BxOutputVertex::Get()->GetWriteDeposits() && !BxOutputVertex::Get()->GetWriteEB()){
	if(theStep->GetTotalEnergyDeposit () >0) {
		BxOutputVertex::Get()->SetDepPDG(theStep->GetTrack()->GetDefinition()->GetPDGEncoding());
		BxOutputVertex::Get()->SetDepEnergy(theStep->GetTotalEnergyDeposit ()/MeV);
		BxOutputVertex::Get()->SetDepPosition(theStep->GetPostStepPoint()->GetPosition()/cm);
		BxOutputVertex::Get()->SetDeposits();
	}
}
}
