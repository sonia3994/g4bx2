//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#include "BxTrackingAction.hh"
#include "BxDataCollection.hh"
#include "BxPropertyCollection.hh"
#include "BxLogger.hh"
#include "BxReadParameters.hh"
#include "BxOutputVertex.hh"
#include "G4SystemOfUnits.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TrackingManager.hh"         
#include "G4Track.hh"
#include "Randomize.hh"
#include "G4Timer.hh"
#include <stdlib.h>
#include "HistoManager.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


extern  G4int thePMTconfig;
extern  G4int thePhotonID;
extern  G4int NbOfPhReflections;


BxTrackingAction::BxTrackingAction() {
  // PMT Cathode Quantum Efficiency
  ReadPMTQE();
  PMTQEVector = PMTQETable->GetProperty("PMTQE");
}


void BxTrackingAction::PostUserTrackingAction(const G4Track* aTrack) {
  G4LogicalVolume*	theLogicalVolume;
  G4VPhysicalVolume*	thePhysicalVolume;
  G4String		theName;
  G4double    		thePhotonEnergy;
  G4int                 BxPMTNumber;
  G4double              PMTQE;
  G4double              RelPMTQE;

  PhotonData		thePhotonData;
  if(BxOutputVertex::Get()->GetWriteIsotope() && aTrack->GetCreatorProcess ()  ) {
      double ene = aTrack->GetVertexKineticEnergy ()/MeV;
      if( 
          (
	    (abs(aTrack->GetDefinition()->GetPDGEncoding()) > 10000 
                 && !aTrack->GetDefinition()->GetPDGStable ()
           ) || (
	     aTrack->GetDefinition()->GetPDGEncoding() == 22 
	     &&  ene > 2.2 && ene < 2.25 &&  aTrack->GetCreatorProcess ()->GetProcessName () == "nCapture")
	   )	
	   
	   && aTrack->GetVolume()->GetName() == "ZoneI") {
	
	int _pdg = aTrack->GetDefinition()->GetPDGEncoding();
	
	if(_pdg != 22) _pdg -= 1000000000;
	
	BxLog(trace) << "Isotope: " 
             << aTrack->GetDefinition()->GetParticleName() << " "
             << aTrack->GetVolume()->GetName() << " "
             << _pdg<< " "
             << aTrack->GetStep()->GetPostStepPoint ()->GetProcessDefinedStep ()->GetProcessName () << " "
             << aTrack->GetStep()->GetPostStepPoint ()->GetProcessDefinedStep ()->GetProcessType () << " "
             << aTrack->GetLocalTime()/ns<< " "
             << aTrack->GetTrackID()<<  " " 
	     << aTrack->GetStep()->GetPostStepPoint ()->GetPosition()
             << endlog ;

	BxOutputVertex::Get()->SetDId(BxOutputVertex::Get()->GetDId()+10000);
	BxOutputVertex::Get()->SetDPDG(_pdg);
	BxOutputVertex::Get()->SetDEnergy(ene);
	BxOutputVertex::Get()->SetDTime(aTrack->GetLocalTime()/ns); 
	BxOutputVertex::Get()->SetDPosition(aTrack->GetVertexPosition());       
	BxOutputVertex::Get()->SetDDirection(aTrack->GetVertexMomentumDirection());	
	BxOutputVertex::Get()->SetNDaughters(BxOutputVertex::Get()->GetNDaughters() + 1);	
	BxOutputVertex::Get()->SetDaughters();

      }
  }

BxDataCollection::Get()->GetPhotonData().push_back(thePhotonData);

  if(aTrack->GetDefinition()->GetPDGEncoding() != 50) {

    if(BxLogger::GetSeverity() == BxLogger::debugging) 
      BxLog(debugging) 
       <<   aTrack->GetDefinition()->GetParticleName() << " "
       <<  " gtime: " << aTrack->GetGlobalTime()/ns << " "
       <<  " ltime: " << aTrack->GetLocalTime()/ns  << " "
       <<  " steps: " << aTrack->GetCurrentStepNumber() << " "
       <<  " E: " <<aTrack->GetKineticEnergy()/eV<< " "
       <<  " ID: " <<aTrack->GetTrackID() << " "
       <<  " PID: " <<aTrack->GetParentID() << " "
       << endlog ;
  }

  
   if(aTrack->GetTrackStatus() == fStopAndKill)	{

    if (aTrack->GetDefinition()->GetPDGEncoding() == 50) { //opticalphoton
        
      thePhotonEnergy = (aTrack->GetTotalEnergy());
      thePhysicalVolume = aTrack->GetStep()->GetPostStepPoint()->GetPhysicalVolume();

G4ThreeVector posizione=aTrack->GetPosition();

      if(thePhysicalVolume)	{
	theLogicalVolume = thePhysicalVolume->GetLogicalVolume();
	
	if(theLogicalVolume->GetName() == "PMT_logic") { 

          PMTQE = PMTQEVector->Value(thePhotonEnergy);

          BxPMTNumber = (thePhysicalVolume->GetCopyNo());

//Here we read the relative PMTQE

          RelPMTQE = BxPropertyCollection::Get()->GetChannelRelQE()[BxPMTNumber];

          PMTQE = RelPMTQE*PMTQE;

//The subtraction of neutron capture time from time of the detected photon (to make the neutron capture event visible by bx_elec)
          G4double Tneutron=0;
          if(BxDataCollection::Get()->GetNeutronCaptureFlag()==1)
		Tneutron=BxDataCollection::Get()->GetNeutronTime(); 

	    	HistoManager::Get()->FillHisto(4, aTrack->GetGlobalTime()/ns,1);

	  if (aTrack->GetGlobalTime()/ns < 35) HistoManager::Get()->FillHisto(2,aTrack->GetTotalEnergy()/eV);
	HistoManager::Get()->FillNtuple(aTrack->GetTotalEnergy()/eV,aTrack->GetGlobalTime()/ns);

	  if(G4UniformRand() < PMTQE) {
	    
	        G4double hit_time = aTrack->GetGlobalTime()-Tneutron;
		
		if (hit_time < 0) // Can occur if the neutron is captured just near the PMT face
		                  // or due to errors in Geant4 neutron capture database 
		{                 
		hit_time = 0.1 * ns;
		G4cout << "hit_time is less then ZERO -attention !" << G4endl;
		}
	    
	    
	      thePhotonID=thePhotonID+1;
	      BxDataCollection::Get()->SetNumberOfPhotons(thePhotonID);
	      BxDataCollection::Get()->GetPhotonData()[thePhotonID-1].SetID(thePhotonID-1);
	      BxDataCollection::Get()->GetPhotonData()[thePhotonID-1].SetFlightTime(hit_time);
	      
	      BxDataCollection::Get()->GetPhotonData()[thePhotonID-1].SetPMTNumber(BxPMTNumber); 
	      BxDataCollection::Get()->SetNumberOfHitPMT(BxDataCollection::Get()->GetNumberOfHitPMT()+1);
	      BxDataCollection::Get()->GetPhotonData()[thePhotonID-1].SetNbOfPhReflections(NbOfPhReflections/2);
	      
	      BxOutputVertex::Get()->SetTOF(hit_time);
	      BxOutputVertex::Get()->SetHole(BxPropertyCollection::Get()->GetChannelNumber()[BxPMTNumber]);
	      BxOutputVertex::Get()->SetPhotons();
	     

	  }
	}

	if(theLogicalVolume->GetName() == "PMTmu_logic") {
	  //theWave = (h_Planck*c_light*1E3)/thePhotonEnergy;
	  //the same (and the maximum one) QE for all the PMTs
	  PMTQE = PMTQEVector->Value(thePhotonEnergy);

	  if(G4UniformRand() < PMTQE) {

	    thePhotonID=thePhotonID+1;
	    BxDataCollection::Get()->SetNumberOfPhotons(thePhotonID);
	    BxDataCollection::Get()->GetPhotonData()[thePhotonID-1].SetID(thePhotonID-1);
	    BxDataCollection::Get()->GetPhotonData()[thePhotonID-1].SetFlightTime(aTrack->GetGlobalTime());	   			       	    
	    BxPMTNumber = (thePhysicalVolume->GetCopyNo());
	    BxDataCollection::Get()->GetPhotonData()[thePhotonID-1].SetPMT_OD_Number(BxPMTNumber); 
	    BxDataCollection::Get()->SetNumberOfHitPMTmu(BxDataCollection::Get()->GetNumberOfHitPMTmu()+1);
	    BxDataCollection::Get()->GetPhotonData()[thePhotonID-1].SetNbOfPhReflections(NbOfPhReflections/2);
	    
	    BxOutputVertex::Get()->SetMuTOF(aTrack->GetGlobalTime());
	    BxOutputVertex::Get()->SetMuHole(BxPMTNumber);
  	    BxOutputVertex::Get()->SetMuPhotons();
	  } 
	}
      }
      
      NbOfPhReflections=0;		
    }
  }   
}

void BxTrackingAction::ReadPMTQE() {

  G4int	    theNum0, theNum1;

  G4double*	pointer1;
  G4double*	pointer2;

  theNum0 = BxPropertyCollection::Get()->GetPhotonEnergyQE().size();
  theNum1 = BxPropertyCollection::Get()->GetPMTQuantumEfficiency().size();


  if(theNum0 != theNum1)	{
    BxLog(warning) << "Different numbers of entries in Wave and Efficiency data :";
    BxLog(warning) << "for Waves table -> " << theNum0 << ", for Efficiency table ";
    BxLog(warning) << "-> " << theNum1 << endlog;
    BxLog(warning) << "  minimum value of entries used!" << endlog;
    theNum0 = min(theNum0,theNum1);
  }

  pointer1 = BxPropertyCollection::Get()->GetPMTQuantumEfficiency();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyQE();

  PMTQETable = new G4MaterialPropertiesTable();
  PMTQETable->AddProperty("PMTQE",pointer2, pointer1, theNum0);
}
