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

#include "BxStackingMuon.hh"
#include "BxStackingAction.hh"
#include "BxOutputVertex.hh"
#include "BxLogger.hh"
#include "G4StackManager.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4ParticleTable.hh"
#include "G4NeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4OpticalPhoton.hh" 
#include "G4SystemOfUnits.hh"
#include "globals.hh"
using namespace std;

BxStackingMuon::BxStackingMuon() {
  /* BxLog(routine) << "Muon Stacking Action Active" << endlog ;   
   fManager = new G4StackManager ;
   NSTACK = 1000;*/
}


BxStackingMuon::~BxStackingMuon(){;}

G4ClassificationOfNewTrack BxStackingMuon::BxClassifyNewTrack (const G4Track* aTrack) {
 G4ClassificationOfNewTrack status = fUrgent;
BxLog(development) << "Muon kinetic energy (MeV) " << aTrack->GetKineticEnergy() << endlog;
 /*  //if(aTrack->GetGlobalTime() > 10*ns ) return fWaiting ;

      //cout << aTrack->GetGlobalTime()/ns<< " "  << aTrack->GetDefinition()->GetParticleName()<< endl ;
      //BxStackCheckStatus();
  if (aTrack->GetDefinition()->GetPDGEncoding() != 50 && aTrack->GetKineticEnergy ()/MeV > 2. ) {
    BxOutputVertex::Get()->SetDId(fDaughters);         
    BxOutputVertex::Get()->SetDPDG(aTrack->GetDefinition()->GetPDGEncoding());      
    BxOutputVertex::Get()->SetDTime(aTrack->GetGlobalTime());           
    BxOutputVertex::Get()->SetDEnergy(aTrack->GetKineticEnergy ()/MeV);
    BxOutputVertex::Get()->SetDPosition(aTrack->GetPosition ()/cm );
    BxOutputVertex::Get()->SetDDirection(aTrack->GetMomentumDirection());
    BxOutputVertex::Get()->SetDaughters();
    fDaughters++;
  }

  if (aTrack->GetDefinition() == G4NeutrinoE::NeutrinoE())             return  fKill;
  if (aTrack->GetDefinition() == G4NeutrinoMu::NeutrinoMu())           return  fKill;
  if (aTrack->GetDefinition() == G4NeutrinoTau::NeutrinoTau())         return  fKill;
  if (aTrack->GetDefinition() == G4AntiNeutrinoE::AntiNeutrinoE())     return  fKill;
  if (aTrack->GetDefinition() == G4AntiNeutrinoMu::AntiNeutrinoMu())   return  fKill;
  if (aTrack->GetDefinition() == G4AntiNeutrinoTau::AntiNeutrinoTau()) return  fKill;

  //cout << fCounter << " " << aTrack->GetGlobalTime()/ns<< " "  << aTrack->GetDefinition()->GetParticleName()
  // << " " << aTrack->GetKineticEnergy()/MeV  << endl ;
 */ return status;
}

void BxStackingMuon::BxNewStage() { 
  //cout << "Wait " << endl ;
   //BxStackCheckStatus();
  //BxStackClearAll();
  //BxStackClearUrgent();
  fDaughters = 0;
  return ;
}

void BxStackingMuon::BxPrepareNewEvent() {  
  fCounter = 0;
  IsValid = false ;
  stage = 0;
  MaxTime = 0*ns ;
}


