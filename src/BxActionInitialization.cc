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
//Created by A. Caminata and S. Marcocci, Sept. 2014
#include "BxActionInitialization.hh"
#include "BxPrimaryGeneratorAction.hh"
#include "BxRunAction.hh"
#include "BxEventAction.hh"
#include "BxSteppingAction.hh"
#include "BxTrackingAction.hh"
#include "BxStackingAction.hh"
#include "BxLogger.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BxActionInitialization::BxActionInitialization()
 : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BxActionInitialization::~BxActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BxActionInitialization::BuildForMaster() const
{
  SetUserAction(new BxRunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BxActionInitialization::Build() const
{
  // Register event generator.
  BxLog(trace) << "Creating and registering event generator" << endlog;
  SetUserAction(new BxPrimaryGeneratorAction);

  // Register run action. What to do at beginning and end of each run.
  BxLog(trace) << "Registering G4 run action." << endlog;
  SetUserAction(new BxRunAction);
  
  // Register event action, ie. what to save/compute for each event.
  BxLog(trace) << "Registering G4 event action." << endlog;
  SetUserAction(new BxEventAction);
  
  // Register tracking action
  BxLog(trace) << "Registering G4 tracking action." << endlog;
  SetUserAction(new BxTrackingAction);
  
  // Register stepping action, ie. what to save/compute for each step. 
  BxLog(trace) << "Registering G4 stepping action." << endlog;
  SetUserAction(new BxSteppingAction);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
