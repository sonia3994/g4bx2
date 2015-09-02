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

#include "BxSensitiveDetector.hh"
#include "BxOutputVertex.hh"
#include "G4SDManager.hh"
#include "BxLogger.hh"
#include "G4SystemOfUnits.hh"
using namespace std;

BxSensitiveDetector::BxSensitiveDetector(G4String name): G4VSensitiveDetector(name)
{   
  BxLog(routine) << "Sensistive Detector " << name << " Builded" << endlog ;
}

BxSensitiveDetector::~BxSensitiveDetector()
{
}

void BxSensitiveDetector::Initialize(G4HCofThisEvent*)
{ // Create hits collection from Alessio-XENON
//	  fHitsCollection= new XeEndcapHitsCollection(SensitiveDetectorName, collectionName[0]);

	    // Add this collection in hce
//G4int hcID       = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
//hce->AddHitsCollection( hcID, fHitsCollection );
}

void BxSensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
}

void BxSensitiveDetector::clear()
{
  trackID.clear();
}


G4bool BxSensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  if(aStep->GetTrack()->GetKineticEnergy()/MeV < 1) return false ;
  
  int pdg = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  if(pdg == -12 || pdg == 12 || pdg == 14 || pdg == -14) return false ;
  for(G4int i=0;i< G4int(trackID.size()); i++) 
    if(aStep->GetTrack()->GetTrackID() == trackID[i]) return false;
  
  trackID.push_back(aStep->GetTrack()->GetTrackID());
  
  BxOutputVertex::Get()->SetDId(BxOutputVertex::Get()->GetDId());
  BxOutputVertex::Get()->SetDPDG(aStep->GetTrack()->GetDefinition()->GetPDGEncoding());
  BxOutputVertex::Get()->SetDEnergy(aStep->GetTrack()->GetKineticEnergy()/MeV);
  BxOutputVertex::Get()->SetDTime(aStep->GetTrack()->GetLocalTime()/ns); 
  BxOutputVertex::Get()->SetDPosition(aStep->GetPreStepPoint ()->GetPosition());       
  BxOutputVertex::Get()->SetDDirection(aStep->GetPreStepPoint ()->GetMomentumDirection());	
  BxOutputVertex::Get()->SetNDaughters(BxOutputVertex::Get()->GetNDaughters() + 1);	
  BxOutputVertex::Get()->SetDaughters();
  
  return false ;

}
