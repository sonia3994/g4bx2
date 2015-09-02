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
#include "BxStackingNeutron.hh"
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
#include "G4Gamma.hh" 
#include "globals.hh"
#include "CLHEP/Random/RandExponential.h"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include "BxEventAction.hh"

using namespace std;

BxStackingNeutron::BxStackingNeutron() {
   IsReclassify = false ;
   BxLog(routine) << "Neutron Stacking Methode Active" << endlog ;   
}


BxStackingNeutron::~BxStackingNeutron(){;}

G4ClassificationOfNewTrack BxStackingNeutron::BxClassifyNewTrack (const G4Track* aTrack) {
const G4VProcess* creator = aTrack->GetCreatorProcess();

	G4int CaptureFlag=0;


  if(creator)
    {
    if (creator->GetProcessName() == "nCapture")  CaptureFlag=1;   
    }
   
 
	if(aTrack->GetDefinition ()->GetParticleName() == "gamma"  && CaptureFlag==1  && aTrack->GetParentID() != -1) 
{

return fPostpone;

}

if(aTrack->GetDefinition ()->GetParticleName() == "gamma"  && CaptureFlag==1 ) 
{
	BxDataCollection::Get()->SetNeutronTime(aTrack->GetGlobalTime()/ns);
	BxDataCollection::Get()->SetNeutronCaptureFlag(1);
	BxOutputVertex::Get()->SetNSequence(0);
	BxOutputVertex::Get()->SetIsotopeCoinc(1);
	BxOutputVertex::Get()->SetEnergy(aTrack->GetKineticEnergy()/MeV);
	BxOutputVertex::Get()->SetPosition(aTrack->GetPosition());
	BxOutputVertex::Get()->SetDirection(aTrack->GetMomentumDirection());
	BxOutputVertex::Get()->SetTime(aTrack->GetGlobalTime()/ns);
	BxOutputVertex::Get()->SetPDG(22);

}

  return fUrgent;
}

void BxStackingNeutron::BxNewStage() { 
  return ;
}

void BxStackingNeutron::BxPrepareNewEvent() {  
 }


