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
#include "BxStackingAction.hh"
#include "BxVStackingAction.hh"
#include "BxStackingActionMessenger.hh"
#include "G4StackManager.hh"
using namespace std;

BxStackingAction::BxStackingAction() {
   fMessenger = new BxStackingActionMessenger(this);
}


BxStackingAction::~BxStackingAction(){
    delete fMessenger;
}


G4ClassificationOfNewTrack  BxStackingAction::ClassifyNewTrack(const G4Track* aTrack) {
  return fStacking->BxClassifyNewTrack(aTrack);
}

void BxStackingAction::NewStage() {
  fStacking->BxNewStage();
}

void BxStackingAction::PrepareNewEvent() {
  fStacking->BxPrepareNewEvent();

}
