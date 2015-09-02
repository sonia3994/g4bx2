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

#include "BxEventAction.hh"
#include "BxOutputVertex.hh"
#include "BxEventActionMessenger.hh"
#include "BxDataCollection.hh"
#include "BxPropertyCollection.hh"
#include "BxReadParameters.hh"
#include "BxLogger.hh"
#include "BxIO.hh"
#include "Randomize.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4String.hh"
#include "G4Timer.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <fstream>
#include <string>
#include "HistoManager.hh"
#include <vector>
#include <algorithm>

using namespace std;

extern  G4int thePhotonID;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BxEventAction::BxEventAction():drawFlag("all"),printModulo(1),eventMessenger(0){
  totnpe = 0;
  Counter = 1000;  
  fWriteLimit = 1; //threshold for writing on binary file
  f64bit = true;
  fCheckDepositPosition = false ;
  timer = new G4Timer;
  timer->Start();
  eventMessenger = new BxEventActionMessenger(this);
  fVisibleEnergy = 0. ;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BxEventAction::~BxEventAction(){
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BxEventAction::BeginOfEventAction(const G4Event* evt){
	BxDataCollection::Get()->SetEID(evt->GetEventID()); 
	//thePhotonID=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BxEventAction::EndOfEventAction(const G4Event* evt){
  /*BinaryStructure	      theBinaryStructure;
  BinaryStructureFake	      theBinaryStructureFake;
  */
  G4int		theEID;
  theEID=evt->GetEventID();
  

  vector<PhotonData>& thePhotonData = BxDataCollection::Get()->GetPhotonData();
  G4ThreeVector BarLight(0.,0.,0.);
  //it computes the barycenter of light
  //G4cout<<"HIT PMT = "<<BxDataCollection::Get()->GetNumberOfHitPMT()<<G4endl;
  //G4cout<<"HIT PMTmu = "<<BxDataCollection::Get()->GetNumberOfHitPMTmu()<<G4endl;
  for(G4int i=0;i<BxDataCollection::Get()->GetNumberOfPhotons();i++)	 
  {    BarLight  += BxReadParameters::Get()->GetPMTPosition()[thePhotonData[i].GetPMTNumber()]/G4float(BxDataCollection::Get()->GetNumberOfPhotons());
  } 
  G4int n_trajectories = 0;
  if (evt->GetTrajectoryContainer()) n_trajectories = evt->GetTrajectoryContainer()->entries();

  BxOutputVertex::Get()->SetEventID(evt->GetEventID());
  BxOutputVertex::Get()->SetNDeposits(int(BxOutputVertex::Get()->GetVDeposits().size()));
  BxOutputVertex::Get()->SetNDaughters(int(BxOutputVertex::Get()->GetVDaughters().size()));
  BxOutputVertex::Get()->SetNUsers(int(BxOutputVertex::Get()->GetVUsers().size()));
  BxOutputVertex::Get()->SetBarycenter(BarLight);
  BxOutputVertex::Get()->SetNPE(int(BxOutputVertex::Get()->GetVPhotons().size()));
  BxOutputVertex::Get()->SetMuNPE(int(BxOutputVertex::Get()->GetVMuPhotons().size()));
 
  totnpe += BxOutputVertex::Get()->GetNPE();
    
  if (!theEID)  { timer->Stop();    timer->Start(); }
 
  //It does something every 1000 events 
  if (theEID%Counter == 0 && theEID) {  
    timer->Stop();
    BxLog(routine) << ">>> Event " << evt->GetEventID()  << ";  NPE = " <<  BxOutputVertex::Get()->GetNPE() << ";  NPE/event = " <<  G4float(totnpe)/G4float(Counter) <<endlog;
   BxLog(trace)   << "    Starting Position: " <<  BxOutputVertex::Get()->GetPosition()/cm << " cm" << endlog;
    BxLog(trace)   << "    Light Barycenter : " << BarLight/cm << " cm" << endlog;
    BxLog(trace)   << "    CPUTime/event    : "<<  timer->GetRealElapsed()/G4float(Counter) << " s" << endlog ;
    BxLog(trace)   << endlog;
  
    if(BxLogger::GetSeverity() <= BxLogger::trace) {
     if(BxOutputVertex::Get()->GetVerbosity() > 0) BxOutputVertex::Get()->DumpVertex();
     if(BxOutputVertex::Get()->GetVerbosity() > 1) BxOutputVertex::Get()->DumpDaughter();
     if(BxOutputVertex::Get()->GetVerbosity() > 2) BxOutputVertex::Get()->DumpDeposit();
     if(BxOutputVertex::Get()->GetVerbosity() > 3) BxOutputVertex::Get()->DumpPhoton();
     if(BxOutputVertex::Get()->GetVerbosity() > 3) BxOutputVertex::Get()->DumpMuPhoton();
     if(BxOutputVertex::Get()->GetVerbosity() > 2) BxOutputVertex::Get()->DumpUser();
   }  

   BxLog(debugging) << "    " << n_trajectories 
      << " trajectories stored in this event." << endlog;
    totnpe = 0;
    timer->Start();
  }
 

  G4int TotalHitPMT;  
  TotalHitPMT = BxDataCollection::Get()->GetNumberOfHitPMT()+ BxDataCollection::Get()->GetNumberOfHitPMTmu();

  HistoManager::Get()->FillHisto(1, TotalHitPMT, 1);
  HistoManager::Get()->FillHisto(3,BxOutputVertex::Get()->GetNPhotons(),1);
  
  if(TotalHitPMT < fWriteLimit) { 
    thePhotonID=0;
    BxDataCollection::Get()->Clear();
    BxOutputVertex::Get()->ClearAll();
    return ; 
  }
  
//Draw trajectory if the graphical UI is istantiated
  if (G4VVisManager::GetConcreteInstance()) {
    G4int n_trajectories_new = 0;
    if (evt->GetTrajectoryContainer()) n_trajectories_new = evt->GetTrajectoryContainer()->entries();
    for (G4int i=0; i<n_trajectories_new; i++)  { 
      G4Trajectory* trj = (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
      if (trj->GetCharge() != 0. && trj->GetPDGEncoding()==50) { 
        trj->DrawTrajectory(); 
      } else if ( i%1000 == 0 && trj->GetCharge() == 0. && trj->GetPDGEncoding()==50 ) {   
	trj->DrawTrajectory();
      }
    }
  }

  BxOutputVertex::Get()->SetNPE(int(BxOutputVertex::Get()->GetVPhotons().size()));
  BxOutputVertex::Get()->SetMuNPE(int(BxOutputVertex::Get()->GetVMuPhotons().size()));

  bool fWriteEvent  = false  ;
 
  
  // Write events if number of npe is larger than the threshold
  if ( BxOutputVertex::Get()->GetNPE() >= fWriteLimit 
    || BxOutputVertex::Get()->GetVisEnergy() >= fVisibleEnergy 
    ||  BxOutputVertex::Get()->GetMuNPE() >= fWriteLimit)  fWriteEvent = true ; 
  

  if(BxOutputVertex::Get()->PostponeFlag() && !BxOutputVertex::Get()->IsPostponed()) fWriteEvent = false;
  
  if(BxOutputVertex::Get()->GetWriteEB() && BxOutputVertex::Get()->GetWriteEBevent()) fWriteEvent = true;
 
  //New output format
  if(!BxIO::Get()->GetIsBinary() && fWriteEvent) { 

  int SIZE =                                          sizeof(VertexStructureDiskFormat) 
      + BxOutputVertex::Get()->GetNDaughters()*44//sizeof(DaughterStructure)
               + BxOutputVertex::Get()->GetNDeposits() *sizeof(DepositStructure)
               + BxOutputVertex::Get()->GetNUsers()    *sizeof(UserStructure)
               + BxOutputVertex::Get()->GetNPE()       *sizeof(PhotonStructure)
               + BxOutputVertex::Get()->GetMuNPE()     *sizeof(MuPhotonStructure);
	       
    BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&SIZE)   ,sizeof( int  ));

    BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVertex())           ,sizeof(VertexStructureDiskFormat));


    if(!f64bit){ 
       for(int i = 0; i <BxOutputVertex::Get()->GetNDaughters(); i++ ) 
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDaughters()[i])  ,sizeof(DaughterStructure));
       for(int i = 0; i <BxOutputVertex::Get()->GetNDeposits(); i++ ) 
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDeposits()[i])   ,sizeof(DepositStructure));  
    }
       
    if(f64bit){ 
       for(int i = 0; i <BxOutputVertex::Get()->GetNDaughters(); i++ ){
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDaughters()[i].Id)  ,sizeof(int));		
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDaughters()[i].PDG)  ,sizeof(int));		
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDaughters()[i].Time)  ,sizeof(double));	
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDaughters()[i].Energy)  ,sizeof(float));	
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDaughters()[i].Position[0])  ,sizeof(float));  
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDaughters()[i].Position[1])  ,sizeof(float));  
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDaughters()[i].Position[2])  ,sizeof(float));  
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDaughters()[i].Direction[0])  ,sizeof(float)); 
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDaughters()[i].Direction[1])  ,sizeof(float)); 
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDaughters()[i].Direction[2])  ,sizeof(float));}
       for(int i = 0; i <BxOutputVertex::Get()->GetNDeposits(); i++ ){
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDeposits()[i].PDGParent)  ,sizeof(int));
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDeposits()[i].Energy)  ,sizeof(float));
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDeposits()[i].Position[0])  ,sizeof(float));
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDeposits()[i].Position[1])  ,sizeof(float));
    	  BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVDeposits()[i].Position[2])  ,sizeof(float));}
    }
     
    
    for(int i = 0; i <BxOutputVertex::Get()->GetNUsers(); i++ ) 
      BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVUsers()[i])      ,sizeof(UserStructure));
   
    for(int i = 0; i <BxOutputVertex::Get()->GetNPE(); i++ ) 
      BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVPhotons()[i])    ,sizeof(PhotonStructure));

    for(int i = 0; i <BxOutputVertex::Get()->GetMuNPE(); i++ ) 
      BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&BxOutputVertex::Get()->GetVMuPhotons()[i])  ,sizeof(MuPhotonStructure));
    
    BxIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&SIZE)	,sizeof( int    ));
  }
 BxOutputVertex::Get()->DumpVertex();
  BxOutputVertex::Get()->SetWriteEBevent(false);
  BxOutputVertex::Get()->ClearAll();

//  To add the Photon Number
  thePhotonID=0;
  BxDataCollection::Get()->Clear();
  BxOutputVertex::Get()->ClearAll();


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
