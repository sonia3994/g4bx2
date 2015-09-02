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

// Make this appear first!
#include "G4Timer.hh"
#include "BxDataCollection.hh"
#include "BxManager.hh"

#include "BxG4BxReader.hh"
#include "BxRunAction.hh"
#include "BxRunActionMessenger.hh"
#include "BxLogger.hh"
#include "BxOutputVertex.hh"
#include "BxIO.hh"
#include "BxReadParameters.hh"
#include "G4ios.hh"
#include "G4Run.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include <time.h>
#include "HistoManager.hh"

class  BxIO;

BxRunAction::BxRunAction() {
  timer = new G4Timer;
  autoSeed = true;
  BxIO::Get()->SetIsBinary(false); //New output format
  BxIO::Get()->SetIsG4Bx(false);
  BxOutputVertex::Get()->ClearHeader();
  runMessenger = new BxRunActionMessenger(this);
}

BxRunAction::~BxRunAction() {
  delete timer;
  delete runMessenger;
}

void BxRunAction::BeginOfRunAction(const G4Run* ) {
 G4UImanager *UI = G4UImanager::GetUIpointer();
 UI->ApplyCommand("/vis/scene/notifyHandlers");
  BxLog(routine) << "### Run " << BxOutputVertex::Get()->GetRun() << " start " << endlog;
  
  timer->Start();

  if(autoSeed) {
    BxLog(routine) << "*******************" << endlog;
    BxLog(routine) << "*** AUTOSEED ON ***" << endlog;
    BxLog(routine) << "*******************" << endlog;
    long seeds[2];
    time_t systime = time(NULL);
    seeds[0] = (long) systime;
    seeds[1] = (long) (systime*G4UniformRand());
    BxLog(routine) << "Seed: " << seeds[1] << endlog;
    CLHEP::HepRandom::setTheSeed(seeds[1]);
  } else {
    BxLog(routine) << "********************" << endlog;
    BxLog(routine) << "*** AUTOSEED OFF ***" << endlog;
    BxLog(routine) << "********************" << endlog; 
  
  }
  
BxIO::Get()->GetStreamLogFile() << "Random seed: " << CLHEP::HepRandom::getTheSeed() << endlog ;
  BxIO::Get()->GetStreamLogFile() << "Random seed: " << CLHEP::HepRandom::getTheSeed() << endlog ;
  
  if(BxIO::Get()->IsG4Bx()) {
    BxIO::Get()->OpenG4BxFile();
    if( BxIO::Get()->GetG4BxFile().fail()) {
      BxLog(error) << "G4Bx file does not exist!" << endlog ;
      BxLog(fatal) << "Fatal " << endlog ;
    }
    BxLog(routine) << "G4Bx file " << BxIO::Get()->GetG4BxFileName() << " opened" << endlog ;
    BxG4BxReader::Get()->ReadHeader();
  }

  BxOutputVertex::Get()->SetEvents(BxManager::Get()->GetCurrentRun()->GetNumberOfEventToBeProcessed());

  BxOutputVertex::Get()->SetPhotonYield(BxReadParameters::Get()->GetLightYield());  
  BxOutputVertex::Get()->SetKB(BxReadParameters::Get()->GetBirksAlpha());     
  BxOutputVertex::Get()->SetKB2(BxReadParameters::Get()->GetBirksSecondOrderAlpha()) ;   

  BxIO::Get()->OpenBinaryFile();
  
  if(BxIO::Get()->GetIsBinary()) BxLog(routine) << "OLD Output Format " << endlog;      
  else BxLog(routine) << "NEW Output Format " << endlog;
   
  BxLog(routine) << "Initialized Binary File: " << BxIO::Get()->GetBinaryFileName() << endlog;      
  if(!BxIO::Get()->GetIsBinary()) { // Header writing
    int SIZE =  sizeof(HeaderStructure) - (10000 - BxOutputVertex::Get()->GetCommentLength())*sizeof(char);
    BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&SIZE)   ,sizeof( int  ));
    BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetHeader()),SIZE);
    BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&SIZE)   ,sizeof( int  ));
  }
}
    
void BxRunAction::EndOfRunAction(const G4Run* aRun) {
  if (G4VVisManager::GetConcreteInstance()) 
        G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  BxIO::Get()->CloseBinaryFile();
  BxIO::Get()->CloseLogFile();  
  BxLog(routine) << "Binary File: " << BxIO::Get()->GetBinaryFileName()<< " closed" << endlog;      
  timer->Stop();
  BxLog(routine) << "Number of event = " << aRun->GetNumberOfEvent() << " " << *timer << endlog;
}
