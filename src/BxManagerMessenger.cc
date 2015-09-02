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
//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014


#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"
#include "BxLogger.hh"
#include "BxManagerMessenger.hh"
#include "BxManager.hh"
#include "BxIO.hh"


BxManagerMessenger::BxManagerMessenger(BxManager *){

  fBxLogCmd = new G4UIcmdWithAString("/bxlog", this);
  fBxLogCmd->SetGuidance("Set severity of logs to report to stdout.");
  fBxLogCmd->SetGuidance("Options, in ascending order of severity, are:");
  fBxLogCmd->SetGuidance("debugging: Displays all logs ");
  fBxLogCmd->SetGuidance("trace: All logs, except debugging(default)");
  fBxLogCmd->SetGuidance("routine: All logs, except debugging and trace");
  fBxLogCmd->SetGuidance("warning: All logs, except trace, debugging and routine:");
  fBxLogCmd->SetGuidance("error: Only error and fatal logs.");
  fBxLogCmd->SetGuidance("fatal: Only fatal logs.");
}


BxManagerMessenger::~BxManagerMessenger() {

  delete fBxLogCmd;
}


void BxManagerMessenger::SetNewValue
(G4UIcommand* command, G4String newValue) { 
	if(command == fBxLogCmd) {
    if(newValue == "debugging")
      BxLogger::SetSeverity(BxLogger::debugging);
    else if(newValue == "trace")
      BxLogger::SetSeverity(BxLogger::trace);
    else if(newValue == "routine")
      BxLogger::SetSeverity(BxLogger::routine);
    else if(newValue == "warning")
      BxLogger::SetSeverity(BxLogger::warning);
    else if(newValue == "error")
      BxLogger::SetSeverity(BxLogger::error);
    else if(newValue == "fatal")
      BxLogger::SetSeverity(BxLogger::fatal);
    else
      BxLog(error) << "Unknown option." << endlog;
    BxIO::Get()->GetStreamLogFile() << "Severity Log output set to " << newValue<< endl;

    G4cout << "     ####################################" <<  endl ;
    G4cout << "     ######## Logger Severity: " << BxLogger::GetSeverity() << " ########" << endl ;
    G4cout << "     ####################################" <<  endl ;
  }      
 }
