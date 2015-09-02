//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
// Created by D. Franco
// Revised by A. Caminata and S. Marcocci, Sept. 2014
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "CLHEP/Random/Random.h"
#include "BxOutputVertex.hh"

#include "BxRunActionMessenger.hh"
#include "BxRunAction.hh"
#include "BxLogger.hh"
#include "BxIO.hh"
#include "BxEventAction.hh"

namespace CLHEP {} 
using namespace CLHEP;

BxRunActionMessenger::BxRunActionMessenger
(BxRunAction* runAct) : runAction(runAct){
  fDirectory = new G4UIdirectory("/run/");

  fSetAutoSeedCmd = new G4UIcmdWithABool("/run/autoSeed",this);
  fSetAutoSeedCmd->SetGuidance("Switch on/off time-based random seeds");
  fSetAutoSeedCmd->SetGuidance("true: run seeds determined by system time");
  fSetAutoSeedCmd->SetGuidance("false: use command 'random/resetEngineFrom'");
  fSetAutoSeedCmd->SetGuidance("Default = true");
  //fSetAutoSeedCmd->SetParameterName("autoSeed", false);
  //fSetAutoSeedCmd->AvailableForStates(G4State_Idle);

  fHEPRandomSeedCmd = new G4UIcmdWithAnInteger("/run/heprandomseed",this);
  fHEPRandomSeedCmd->SetGuidance("Sets random number generator seed.");


  // /run/rootfilename
  fSetFileNameCmd = new G4UIcmdWithAString("/run/filename", this);
  fSetFileNameCmd->SetGuidance("Name for output files");
 
  fSetIsBinaryCmd = new G4UIcmdWithABool("/run/setbinary", this);
  fSetIsBinaryCmd->SetGuidance("Set Binary");

  fSetG4BxNameCmd = new G4UIcmdWithAString("/run/g4bx", this);
  fSetG4BxNameCmd->SetGuidance("Name of the  output files");

  //fSetFileNameCmd->SetGuidance("Decay nuclide");

  fRunCmd = new G4UIcmdWithAnInteger("/run/id",this);
  
  fRateCmd = new G4UIcmdWithADoubleAndUnit("/run/rate",this);
  fRateCmd->SetUnitCandidates("Hz kHz MHz");
  
  fWriteDepositsCmd = new G4UIcmdWithABool("/run/writedeposits", this);

  fWriteIsotopesCmd = new G4UIcmdWithABool("/run/writeisotopes", this);

  fWriteEBCmd = new G4UIcmdWithABool("/run/writeeb", this);

  fSetCommentCmd = new G4UIcmdWithAString("/run/comment", this);

}


BxRunActionMessenger::~BxRunActionMessenger() {

  delete fSetAutoSeedCmd;  
  delete fHEPRandomSeedCmd;
  delete fSetFileNameCmd;
  delete fSetIsBinaryCmd;
  delete fDirectory;
  delete fRunCmd;
  delete fRateCmd;
  delete fWriteDepositsCmd;
  delete fSetCommentCmd;
  delete fSetG4BxNameCmd;
  delete fWriteIsotopesCmd;
  delete fWriteEBCmd;
}


void BxRunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue) { 

  if(command == fSetAutoSeedCmd) {
    runAction->SetAutoSeed(fSetAutoSeedCmd->ConvertToBool(newValue));    
  } else if (command == fHEPRandomSeedCmd) {
    runAction->SetAutoSeed(false);
    HepRandom::setTheSeed(fHEPRandomSeedCmd->GetNewIntValue(newValue));
    BxLog(trace) << "HepRandom seed set to: "<< fHEPRandomSeedCmd->GetNewIntValue(newValue) << endlog;
  } else if(command == fSetFileNameCmd) {      
    ; // Already set in g4bx.cc
  } else if(command == fSetIsBinaryCmd) {
    BxIO::Get()->SetIsBinary(fSetIsBinaryCmd->ConvertToBool(newValue));
  } else if(command == fSetG4BxNameCmd) {      
      BxIO::Get()->SetG4BxFileName(newValue);
      BxIO::Get()->SetIsG4Bx(true);  
  } else if(command == fRunCmd) { 
      BxOutputVertex::Get()->SetRun(fRunCmd->ConvertToInt(newValue));
  } else if(command == fSetCommentCmd) { 
      BxOutputVertex::Get()->SetComment(newValue);
  } else if(command == fRateCmd) { 
      BxOutputVertex::Get()->SetRate(fRateCmd->ConvertToDouble(newValue));
  } else if(command == fWriteDepositsCmd) {
      BxLog(routine) << "Write Deposits: "<< newValue << endlog;
      BxOutputVertex::Get()->SetWriteDeposits(fWriteDepositsCmd->ConvertToBool(newValue));    
  } else if(command == fWriteIsotopesCmd) {
      BxLog(routine) << "Write Isotopes: "<< newValue << endlog;
      BxOutputVertex::Get()->SetWriteIsotope(fWriteIsotopesCmd->ConvertToBool(newValue));    
  } else if(command == fWriteEBCmd) {
      BxLog(routine) << "Write External Background: "<< newValue << endlog;
      BxOutputVertex::Get()->SetWriteEB(fWriteEBCmd->ConvertToBool(newValue));    
  }
}
