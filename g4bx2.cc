//---------------------------------------------------------------------------//
//Bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//                                                                           //
//                                                                           //
//                         G4Bx2 Simulation                                   //
//                                                                           //
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//
// File:        g4bx2.cc
// Description: Test of Continuous Process G4Cerenkov
//              and RestDiscrete Process G4Scintillation
//              -- Generation Cerenkov Photons --
//              -- Generation Scintillation Photons --
//              -- Transport of optical Photons --
// Version:     5.0
// Created:     2002
// Major Revision: A. Caminata and S. Marcocci, Sept. 2014
// --------------------------------------------------------------


#include "BxLogger.hh"
#include "BxManager.hh"
#include "BxPropertyCollection.hh"
#include "BxReadParameters.hh"
#include "G4VUserDetectorConstruction.hh"

#include "BxRunAction.hh"
#include "BxActionInitialization.hh"
#include "BxPhysicsList.hh"
#include "BxDetectorConstruction.hh"
#include "BxPrimaryGeneratorAction.hh"
#include "BxEventAction.hh"
#include "BxTrackingAction.hh"
#include "BxSteppingAction.hh"
#include "BxStackingAction.hh"

#include "BxPropertyCollection.hh"
#include "BxDataCollection.hh"

//#include "BxLightSource.hh"
#include "BxIO.hh"
#include "G4UImanager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include "G4ios.hh"
#include <stdlib.h>
#include "HistoManager.hh"


//Global variable to count the amount of photons for each event
G4int thePhotonID=0;
G4int NbOfPhReflections=0;

//BxLightSource*                  theLightSource;


// Functions called from main().
void PrintHeader();
void PrintUsage(void);
void WriteSeparationOnLogFile();
void WriteHeaderOnLogfile();

using namespace std;


int main(int argc,char** argv) {
	PrintHeader();
	G4bool IsStack = false ;
	//---------------------------------------------------------------------------//
	// Definition of primary parameters 
	BxLogger::SetSeverity(BxLogger::routine);
	BxIO::Get()->CheckFileName("output");
	if (argc==2) {
		ifstream ifs(argv[1]);
		string st ;
		while (getline (ifs, st)) {
			if (!st.find("/bxlog")) {
				BxLogger::SetSeverity(BxLogger::toEnum(st.substr(st.find(" ")+1)));
			} else if (!st.find("/run/filename")) {
				BxIO::Get()->CheckFileName(st.substr(st.find(" ")+1));
			} else if (!st.find("/run/borexproperty")) {
				BxIO::Get()->SetBorexProperty(st.substr(st.find(" ")+1));
				BxLog(routine) << "BorexProperty file: " << BxIO::Get()->GetBorexProperty() << endlog;
			}  else if (!st.find("/bx/stack")) {
				IsStack = true ;
			}  else if (st.find("/bx/detector/run_number")!=std::string::npos) {
				BxReadParameters::Get()->SetRunNumber(atoi(st.substr(st.find(" ")+1).c_str()));
			} else if (st.find("/scintillator/PPOattenuationLengthFactor")!=std::string::npos){
				BxReadParameters::Get()->SetPPOAttenuationLengthFactor(atof(st.substr(st.find(" ")+1).c_str()));
			} else if (st.find("/scintillator/PCattenuationLengthFactor")!=std::string::npos){
				BxReadParameters::Get()->SetPCAttenuationLengthFactor(atof(st.substr(st.find(" ")+1).c_str()));
			} else if (st.find("/scintillator/DMPattenuationLengthFactor")!=std::string::npos){
				BxReadParameters::Get()->SetDMPAttenuationLengthFactor(atof(st.substr(st.find(" ")+1).c_str()));
			} else if (st.find("/scintillator/NylonattenuationLengthFactor")!=std::string::npos){
				BxReadParameters::Get()->SetNylonAttenuationLengthFactor(atof(st.substr(st.find(" ")+1).c_str()));
			} else if (st.find("/scintillator/PPOAbsReemProbXRays")!=std::string::npos){
				BxReadParameters::Get()->SetPPOAbsReemProbXRays(atof(st.substr(st.find(" ")+1).c_str()));
			} else if (st.find("/scintillator/PPOAbsReemProbThXRays")!=std::string::npos){
				BxReadParameters::Get()->SetPPOAbsReemProbXRaysTh(atof(st.substr(st.find(" ")+1).c_str()));
			}
		} ifs.close();
	}	
	
	BxLog(routine) << "Output file name: " << BxIO::Get()->GetFileName() <<   endlog;
	HistoManager::Get()->SetFileName(BxIO::Get()->GetFileName()+G4String("sim.root"));
	HistoManager::Get()->book();
	BxIO::Get()->OpenLogFile();

	WriteHeaderOnLogfile();
	
	if (argc==2) {
		ifstream ifst(argv[1]);
		string st ;
	
		while (getline (ifst, st)) 
			BxIO::Get()->GetStreamLogFile() << st << endl  ;
	
		WriteSeparationOnLogFile();
		ifst.close();
	}  

	//---------------------------------------------------------------------------//

	BxLog(trace) << "Creating G4 Run Manager" << endlog;
	BxManager* runManager   = BxManager::Get();     
        G4UImanager* UImanager = G4UImanager::GetUIpointer();
	
  // Register detector geometry and materials.
  BxLog(trace) << "Creating and registering G4 geometry" << endlog;
 runManager->SetUserInitialization(new BxDetectorConstruction());

  // Register Geant4 physics processes
  BxLog(trace) << "Creating and registering G4  physics processes" << endlog;
  runManager->SetUserInitialization(new BxPhysicsList());
  
 // Initialize UserActions
   runManager->SetUserInitialization(new BxActionInitialization());
 
   // Register stacking action, ie. what to save/compute for each step.
  if(IsStack) {
    BxLog(trace) << "Registering G4 stacking action." << endlog;
    runManager->SetUserAction(new BxStackingAction());
  }

#ifdef G4VIS_USE
  // Initialize visualization
     G4VisManager* visManager = new G4VisExecutive;
           visManager->Initialize();
           #endif
if (argc!=1) {
		// batch mode
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UImanager->ApplyCommand(command+fileName);
		runManager->Initialize();
}
   else { // interactive mode : define UI session
	   runManager->Initialize();
#ifdef G4UI_USE
	   G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
	   UImanager->ApplyCommand("/control/execute init_vis.mac");
#else
	   UImanager->ApplyCommand("/control/execute init.mac");
#endif
	   ui->SessionStart();
	   delete ui;
#endif
   }
  // job termination
#ifdef G4VIS_USE
   delete visManager;
#endif
   delete runManager;
HistoManager::Get()->save();
   //delete theLightSource;  
	}


//---------------------------------------------------------------------------//

void PrintHeader(void)
{
  G4cout << "Borexino Monte Carlo Simulation" << G4endl;
  G4cout << "------------------------------------------" << G4endl;
  G4cout << "Version: 1.05, CVS Tag  G4Bx-01-05-00" << G4endl;
  G4cout << "Last Update: 11-09-2006" << G4endl;
}

//---------------------------------------------------------------------------//

void PrintUsage(void)
{
  G4cout << "Usage:" << G4endl;
  G4cout << "g4bx2 -h : Displays this message" << G4endl;
  G4cout << "g4bx2 <filename> : Executes script <filename>" << G4endl;
  G4cout << "g4bx2 : Executes G4Bx interactively" << G4endl << G4endl;
}

void WriteHeaderOnLogfile(){
	BxIO::Get()->GetStreamLogFile() << endl;
	BxIO::Get()->GetStreamLogFile() <<"------------------------------------------" << endl;
	BxIO::Get()->GetStreamLogFile() <<"                   G4Bx                   " << endl;
	BxIO::Get()->GetStreamLogFile() << endl;
	BxIO::Get()->GetStreamLogFile() <<"     The Borexino Geant4 Simulator        " << endl;
	BxIO::Get()->GetStreamLogFile() <<"------------------------------------------" << endl;
	BxIO::Get()->GetStreamLogFile() <<  endl;
	BxIO::Get()->GetStreamLogFile() <<  endl;
	BxIO::Get()->GetStreamLogFile() <<"------------------------------------------" << endl;
	BxIO::Get()->GetStreamLogFile() <<"                   Input                  " << endl;
	BxIO::Get()->GetStreamLogFile() <<"------------------------------------------" << endl;
	BxIO::Get()->GetStreamLogFile() <<  endl;
}


void WriteSeparationOnLogFile(){
		BxIO::Get()->GetStreamLogFile() <<  endl;
		BxIO::Get()->GetStreamLogFile() <<  endl;
		BxIO::Get()->GetStreamLogFile() <<"------------------------------------------" << endl;
		BxIO::Get()->GetStreamLogFile() <<"                   Output                 " << endl;
		BxIO::Get()->GetStreamLogFile() <<"------------------------------------------" << endl;
		BxIO::Get()->GetStreamLogFile() <<  endl;
}
	//---------------------------------------------------------------------------//
