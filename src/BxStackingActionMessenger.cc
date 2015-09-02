//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "BxReadParameters.hh"
#include "BxStackingActionMessenger.hh"
#include "BxStackingAction.hh"
//#include "BxStackingEB.hh"
//#include "BxStackingMuon.hh"
#include "BxStackingRDM.hh"
#include "BxStackingRDMChain.hh"
#include "BxStackingNeutron.hh"
#include "BxOutputVertex.hh"

using namespace std;

BxStackingActionMessenger::BxStackingActionMessenger(BxStackingAction* stack){
  fStacking = stack ;
  fDirectory = new G4UIdirectory("/bx/stack/");
  fDirectory->SetGuidance("Control commands for stack:");
  // /Bx/generator/select command
  fSelectCmd = new G4UIcmdWithAString("/bx/stack/select", this);
  fSelectCmd->SetGuidance("Selects stack options");
  G4String candidates = "Default EB Muon RDM RDMChain Neutron";
  fSelectCmd->SetCandidates(candidates);
 }


BxStackingActionMessenger::~BxStackingActionMessenger() {

  delete fDirectory;
  delete fSelectCmd;
//  if(fStacking)  delete fStacking ;
}


void BxStackingActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue) { 
   if(command == fSelectCmd) { /*
     if(newValue == "EB")        fStacking->SetBxStackingAction( new BxStackingEB );
     else if(newValue == "Muon") fStacking->SetBxStackingAction( new BxStackingMuon) ;
    else*/ if(newValue == "RDM")  fStacking->SetBxStackingAction( new BxStackingRDM) ;
     else if(newValue == "Neutron")  fStacking->SetBxStackingAction( new BxStackingNeutron) ;
     else if(newValue == "RDMChain") {
        fStacking->SetBxStackingAction( new BxStackingRDMChain) ;
        BxOutputVertex::Get()->SetChain(1);
     }
   }
} 
