//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#include "BxReadParameters.hh"
#include "BxOutputVertex.hh"
#include "BxDataCollection.hh"
#include "BxStackingRDMMessenger.hh"
#include "BxStackingRDM.hh"
#include "BxLogger.hh"
#include "G4StackManager.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4ParticleTable.hh"
#include "G4NeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4NeutrinoTau.hh"
#include "G4SystemOfUnits.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4OpticalPhoton.hh" 
#include "G4GenericIon.hh" 
#include "globals.hh"
#include "CLHEP/Random/RandExponential.h"
#include "BxReadParameters.hh"

using namespace std;

BxStackingRDM::BxStackingRDM() {

   IsReclassify = false ;
   fRadius = 425.*cm ;
   fKillOuterDirected = false ;
   fMessenger = new BxStackingRDMMessenger(this);
   
   BxLog(routine) << "RDM Stacking Methode Active" << endlog ;  
 
}


BxStackingRDM::~BxStackingRDM(){
delete fMessenger;
}

G4ClassificationOfNewTrack BxStackingRDM::BxClassifyNewTrack (const G4Track* aTrack) {
 if(aTrack->GetDefinition()->GetAtomicNumber () > 2 && aTrack->GetParentID() != 0){
    return fKill ;
}

  BxOutputVertex::Get()->SetIsPostponed(false);
  
   if(aTrack->GetDefinition()->GetParticleName()!="opticalphoton" && aTrack->GetParentID() != 0)
   (const_cast<G4Track *>(aTrack))->SetGlobalTime( 0. ); 
/*
  BxLog(development) 
   <<   aTrack->GetDefinition()->GetParticleName() << " "
   <<  " gtime: " << aTrack->GetGlobalTime()/ns << " "
   <<  " ltime: " << aTrack->GetLocalTime()/ns  << " "
   <<  " steps: " << aTrack->GetCurrentStepNumber() << " "
   <<  " E: "     <<aTrack->GetKineticEnergy()/eV<< " "
   <<  " ID: "    <<aTrack->GetTrackID() << " "
   <<  " PID: "   <<aTrack->GetParentID() << " " 
  << endlog ;

  if(fKillOuterDirected && aTrack->GetDefinition()->GetPDGEncoding() == 22) {
   G4double cosangle = aTrack->GetMomentumDirection ().dot(aTrack->GetPosition());
   if(aTrack->GetMomentumDirection ().mag() &&aTrack->GetPosition().mag() ) {
     cosangle /= aTrack->GetMomentumDirection ().mag()*aTrack->GetPosition().mag() ;
     if(M_PI - acos(cosangle) > fAngle/rad ) return fKill ; 
   }
  }
*/

 if(aTrack->GetParentID() == 0) {
    BxReadParameters::Get()->SetRealPDGMeanLife(aTrack->GetDefinition()->GetPDGLifeTime());
    BxDataCollection::Get()->SetfDTEMPOASS(CLHEP::RandExponential::shoot(BxReadParameters::Get()->GetRealPDGMeanLife()));
    aTrack->GetDefinition()->SetPDGLifeTime(0*s);
    return fUrgent;
  }

  /*
  if(    (aTrack->GetCurrentStepNumber() == 0) 
	 &&  (aTrack->GetDefinition()->GetPDGEncoding() != 50) 
     && (abs(aTrack->GetParentID()) == 1) ) {

    for(G4int i = 0; i < G4int( fPDGToBeKilled.size() ); i++ ){ 
	    if (aTrack->GetDefinition()->GetPDGEncoding() == fPDGToBeKilled[i])
		    return fKill;
    }

    for(G4int i = 0; i < G4int( fLEPDGToBeKilled.size() ); i++ ){ 
	    if (  aTrack->GetDefinition()->GetPDGEncoding() == fLEPDGToBeKilled[i] 
			    && aTrack->GetKineticEnergy() < fLEnergyToBeKilled[i]) 
		    return fKill;
    }

    BxLog(development)  <<   aTrack->GetDefinition()->GetParticleName() << " "
	  <<  " gtime: " << aTrack->GetGlobalTime() << " "
	  <<  " ltime: " << aTrack->GetLocalTime()/ns  << " "
	  <<  " steps: " << aTrack->GetCurrentStepNumber() << " "
	  <<  " E: " <<aTrack->GetKineticEnergy()/eV<< " "
	  <<  " ID: " <<aTrack->GetTrackID() << " "
	  <<  " PID: " <<aTrack->GetParentID() << " "
	  <<  " Proper time: " <<aTrack->GetProperTime ()
	  << endlog ;
     
     BxOutputVertex::Get()->SetDId(fDaughters);		
     BxOutputVertex::Get()->SetDPDG(aTrack->GetDefinition()->GetPDGEncoding());	     
     BxOutputVertex::Get()->SetDTime(0.0);	     
     BxOutputVertex::Get()->SetDEnergy(aTrack->GetKineticEnergy ()/MeV);
     BxOutputVertex::Get()->SetDPosition(aTrack->GetPosition ()/cm );
     BxOutputVertex::Get()->SetDDirection(aTrack->GetMomentumDirection());
     BxOutputVertex::Get()->SetDaughters();
     fDaughters++;

  }
*/
  if (aTrack->GetDefinition() == G4NeutrinoE::NeutrinoE())             return fKill;
  if (aTrack->GetDefinition() == G4NeutrinoMu::NeutrinoMu())           return fKill;
  if (aTrack->GetDefinition() == G4NeutrinoTau::NeutrinoTau())         return fKill;
  if (aTrack->GetDefinition() == G4AntiNeutrinoE::AntiNeutrinoE())     return fKill;
  if (aTrack->GetDefinition() == G4AntiNeutrinoMu::AntiNeutrinoMu())   return fKill;
  if (aTrack->GetDefinition() == G4AntiNeutrinoTau::AntiNeutrinoTau()) return fKill;

/*
  if(BxOutputVertex::Get()->GetWriteEB() && aTrack->GetDefinition()->GetPDGEncoding() == 11) {
     BxOutputVertex::Get()->SetDepPDG(aTrack->GetDefinition()->GetPDGEncoding());
     BxOutputVertex::Get()->SetDepEnergy(aTrack->GetKineticEnergy ()/MeV);
     BxOutputVertex::Get()->SetDepPosition(aTrack->GetPosition ()/cm);
     BxOutputVertex::Get()->SetDeposits();
     if(aTrack->GetPosition ().mag() < BxOutputVertex::Get()->GetKillRadius()) 
                                                 BxOutputVertex::Get()->SetWriteEBevent(true); 
     
          return fKill;
  }
*/

  G4ClassificationOfNewTrack status = fUrgent;
/*
  if(IsReclassify) {
    switch(stage) {
      case 0: // Stage 0 
	// Here, I identified and suspended all photons at their generation point. 
	// If at least a compton does not reach the bulk, all photons are thrown away

	if(aTrack->GetParentID()==0) return fUrgent ;
	 if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhoton()) {
	  if(aTrack->GetGlobalTime()/ns > 200. ) return fKill ;
	  fCounter++;
	  return  fWaiting ;
	}



	if(aTrack->GetPosition ().r() < fRadius) IsValid = true ;


	break;

      case 1: // Stage 1

	// Here, you can set some further cut on the photons belonging to the waiting stack
	break; 
    }
  }    
*/
  return status;

}

void BxStackingRDM::BxNewStage() { 
/*  if(IsReclassify) {
    stage++;
    if(!IsValid)  BxStackClearAll();
    else          BxStackReClassify(); 
}
*/
  return ;
}

void BxStackingRDM::BxPrepareNewEvent() {  
  //fCounter = 0;
  //IsValid = false ;
  //stage = 0;
  //fDaughters = 0;
  //fAlpha = 0;
  //fBeta = 0;
  //fGamma = 0;
  //IsTheFirstDaughter = false ;
}
