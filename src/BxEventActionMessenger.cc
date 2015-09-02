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
//
//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "BxEventActionMessenger.hh"
#include "CLHEP/Random/Random.h"
#include <iostream>
#include <time.h>
#include <sys/times.h>
#include "BxEventAction.hh"
#include "BxRunAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcommand.hh"
#include "globals.hh"
#include "BxLogger.hh"
#include "BxOutputVertex.hh"
#include "BxReadParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BxEventActionMessenger::BxEventActionMessenger(BxEventAction* EvAct)
:eventAction(EvAct) { 
  fDrawCmd = new G4UIcmdWithAString("/event/drawTracks",this);
  fDrawCmd->SetGuidance("Draw the tracks in the event");
  fDrawCmd->SetGuidance("  Choice : none, charged, neutral, all");
  fDrawCmd->SetParameterName("choice",true);
  fDrawCmd->SetDefaultValue("all");
  fDrawCmd->SetCandidates("none charged neutral all");
  fDrawCmd->AvailableForStates(G4State_Idle);
  
  fPrintCmd = new G4UIcmdWithAnInteger("/event/printModulo",this);
  fPrintCmd->SetGuidance("Print events modulo n");
  fPrintCmd->SetParameterName("EventNb",false);
  fPrintCmd->SetRange("EventNb>0");
  fPrintCmd->AvailableForStates(G4State_Idle);     

  fEmptyCmd = new G4UIcmdWithAnInteger("/event/npethreshold",this);
  
  fCounterCmd = new G4UIcmdWithAnInteger("/event/eventreport",this);
  
  fVerbosityCmd= new G4UIcmdWithAnInteger("/event/verbosity",this);
  
  fSetDepoRadiusCmd = new G4UIcmdWithADoubleAndUnit("/event/deporadialcut", this);
  fSetDepoRadiusCmd->SetDefaultUnit("m");
  fSetDepoRadiusCmd->SetUnitCandidates("m cm mm");

  fSetVisEneCmd = new G4UIcmdWithADoubleAndUnit("/event/visibleenergycut", this);
  fSetVisEneCmd->SetDefaultUnit("MeV");
  fSetVisEneCmd->SetUnitCandidates("eV keV MeV GeV");

  fSetTimeCutCmd = new G4UIcmdWithADoubleAndUnit("/event/timecut", this);
  fSetTimeCutCmd->SetDefaultUnit("ns");
  fSetTimeCutCmd->SetUnitCandidates("ns mus ms s");
 
 
  fBitCmd = new G4UIcmdWithABool("/event/64bit", this);
  fBitCmd->SetGuidance("Able or disable writing of events on 64bit architecture");
  fBitCmd->SetGuidance("Default: false");
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BxEventActionMessenger::~BxEventActionMessenger() {
  delete fDrawCmd;
  delete fPrintCmd;
  delete fCounterCmd;  
  delete fEmptyCmd;
  delete fVerbosityCmd;
  delete fSetDepoRadiusCmd;
  delete fSetVisEneCmd;
  delete fSetTimeCutCmd;
  delete fBitCmd;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BxEventActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue){ 
  if(command == fDrawCmd){
    eventAction->SetDrawFlag(newValue);
  } else if(command == fPrintCmd) {
    eventAction->SetPrintModulo(fPrintCmd->GetNewIntValue(newValue));
  } else if(command == fCounterCmd) {
    eventAction->SetCounter(fCounterCmd->GetNewIntValue(newValue));
  } else if(command == fEmptyCmd) {
    BxLog(routine) << "Discard events with npe <  "<< newValue  << endlog ;
    eventAction->SetWriteLimit(fEmptyCmd->GetNewIntValue(newValue));
  } else if(command ==fVerbosityCmd ) {
    BxOutputVertex::Get()->SetVerbosity(fVerbosityCmd->GetNewIntValue(newValue));
  } else if(command == fSetDepoRadiusCmd) {
    eventAction->SetDepoRadius(fSetDepoRadiusCmd->ConvertToDimensionedDouble(newValue));
    BxOutputVertex::Get()->SetWriteDeposits(true);
    BxLog(routine) << "Discard events without any deposit with r <  "<< newValue  << endlog ;
  } else if(command == fSetVisEneCmd) {
    eventAction->SetVisEnergyCut(fSetVisEneCmd->ConvertToDimensionedDouble(newValue));
    BxLog(routine) << "Discard events with visible energy < "<< newValue  << endlog ;
  } else if(command == fSetTimeCutCmd) {
    BxReadParameters::Get()->SetIsTimeCut(true);
    BxReadParameters::Get()->SetTimeCut(fSetTimeCutCmd->ConvertToDimensionedDouble(newValue));
    BxLog(routine) << "Kill particles with global time > "<< newValue  << endlog ;
  }else if(command == fBitCmd) {
    eventAction->Set64Bit(fBitCmd->ConvertToBool(newValue));
    BxLog(routine) << "Activation for event writing on 64bit architecture: " << newValue  <<  endlog;   
  }
}
