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
// Created by D. Franco
// Revised by A. Caminata and S. Marcocci, Sept. 2014
//
#include "BxReadParameters.hh"
#include "BxDataCollection.hh"
#include "BxStackingEBMessenger.hh"
#include "BxStackingEB.hh"
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
#include "G4GenericIon.hh" 
#include "globals.hh"
using namespace std;

BxStackingEB::BxStackingEB() {  
/*
   IsReclassify = false ;
   fRadius = 425.*cm ;
   fKillOuterDirected = false ;
   fSupDef = false;
   fBufDef = false;
   ffRhocut = false;
   ffZcut = false;
   fKillVessel = false ;
   isRead = false;
   fVthick = 0.0000125*m;
   VSloaded = false; 
   frhoCut = 10*m;
   fzCut = 10*m;  
   fMessenger = new BxStackingEBMessenger(this);
   BxLog(routine) << "EB Stacking Methode Active" << endlog ;   
*/}


BxStackingEB::~BxStackingEB(){;}


G4ClassificationOfNewTrack BxStackingEB::BxClassifyNewTrack (const G4Track* aTrack) {
//dummy cout to avoid warning before everything is decommented
	BxLog(development) << "Particle tracked with energy (MeV)" << aTrack->GetKineticEnergy() << endlog;
	/* if(ffRhocut){
    zz = aTrack->GetPosition().z();  
    rr = aTrack->GetPosition().mag();
    if(sqrt(rr*rr - zz*zz) > frhoCut) return fKill;
    
  }
  if(ffZcut){
    zz = aTrack->GetPosition().z();  
    if(zz > fzCut) return fKill;
  }
  
  if(fSupDef){
    if(aTrack->GetParentID()==0){
      rr = aTrack->GetPosition().mag();
      zz = aTrack->GetPosition().z(); 
      if (VSloaded == false){ 
         BxReadParameters::Get()->ReadVesselGeometry();
         RZDim = G4int((BxReadParameters::Get()->GetZDVessel()).size()); 
	 VSloaded = true;}
      G4double Rvessel = 0;
      G4double Rhovessel = 0; 
      for(G4int i = 0; i < RZDim; i++) { //il file parte da z massimo e scende
          if(zz < (BxReadParameters::Get()->GetZDVessel())[i] && zz > (BxReadParameters::Get()->GetZDVessel())[i+1]){
          G4double dr = (BxReadParameters::Get()->GetRDVessel())[i+1]-(BxReadParameters::Get()->GetRDVessel())[i];
          G4double dz = (BxReadParameters::Get()->GetZDVessel())[i+1]-(BxReadParameters::Get()->GetZDVessel())[i];
          Rhovessel = (BxReadParameters::Get()->GetRDVessel())[i] + dr/dz*(zz-(BxReadParameters::Get()->GetZDVessel())[i]);
          Rvessel = pow(Rhovessel*Rhovessel + zz*zz,0.5);
	  break;
        } 
      }
      
      if(rr > Rvessel-(1/1000000.)*m || rr < Rvessel-fVthick-(1/1000000.)*m ) { 
        return fKill ;}
      
    } 
  }
  
  
  if(fBufDef){
      if(aTrack->GetParentID() == 0){
      rr = aTrack->GetPosition().mag();
      zz = aTrack->GetPosition().z(); 
      G4double Rmax = 4.7*m;
      if (VSloaded == false){ 
         BxReadParameters::Get()->ReadVesselGeometry();
         RZDim = G4int((BxReadParameters::Get()->GetZDVessel()).size()); 
	 VSloaded = true;}
      if(rr < Rmax){
        G4double Rvessel = 0;
        G4double Rhovessel = 0; 
        for(G4int i = 0; i < RZDim; i++) { //il file parte da z massimo e scende
          if(zz < (BxReadParameters::Get()->GetZDVessel())[i] && zz > (BxReadParameters::Get()->GetZDVessel())[i+1]){
            G4double dr = (BxReadParameters::Get()->GetRDVessel())[i+1]-(BxReadParameters::Get()->GetRDVessel())[i];
            G4double dz = (BxReadParameters::Get()->GetZDVessel())[i+1]-(BxReadParameters::Get()->GetZDVessel())[i];
            Rhovessel = (BxReadParameters::Get()->GetRDVessel())[i] + dr/dz*(zz-(BxReadParameters::Get()->GetZDVessel())[i]);
            Rvessel = pow(Rhovessel*Rhovessel + zz*zz,0.5);
            break;
          }          
        }
        if(rr < Rvessel) { return fKill ;}
      }
    }
  }
  

   if (aTrack->GetDefinition() == G4NeutrinoE::NeutrinoE())             return fKill;
   if (aTrack->GetDefinition() == G4NeutrinoMu::NeutrinoMu())           return fKill;
   if (aTrack->GetDefinition() == G4NeutrinoTau::NeutrinoTau())         return fKill;
   if (aTrack->GetDefinition() == G4AntiNeutrinoE::AntiNeutrinoE())     return fKill;
   if (aTrack->GetDefinition() == G4AntiNeutrinoMu::AntiNeutrinoMu())   return fKill;
   if (aTrack->GetDefinition() == G4AntiNeutrinoTau::AntiNeutrinoTau()) return fKill;


*/
  
  G4ClassificationOfNewTrack status = fUrgent;
/*
  switch(stage) {
    case 0: // Stage 0 
      // Here, I identified and suspended all photons at their generation point. 
      // If at least a compton does not reach the bulk, all photons are thrown away

      if(aTrack->GetParentID()==0) return fUrgent ;
       if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhoton()) {
	if(aTrack->GetGlobalTime()/ns > 200 ) return fKill ;
	return  fWaiting ;
      }
     


      if(aTrack->GetVolume()->GetName() == "ZoneI") IsValid = true ;


      break;

    case 1: // Stage 1

      // Here, you can set some further cut on the photons belonging to the waiting stack
      break; 
  } */ 
  return status;
}


void BxStackingEB::BxNewStage() { 
/*
 stage++;
  if(!IsValid)  BxStackClearAll();
  else          BxStackReClassify(); 
  //ALREADY COMMENTED
   candidates
  BxStackAbort() 	       
  BxStackClearAll()	       
  BxStackClearUrgentAndWaiting()
  BxStackClearWaiting()         
  BxStackClearUrgent()	       
  BxStackClearPostponed()       
  BxStackCheckStatus()	       
  BxStackReClassify()	    
    
  
  return ;*/
}

void BxStackingEB::BxPrepareNewEvent() {  
 //  IsValid = false ;
 // stage = 0;
}


