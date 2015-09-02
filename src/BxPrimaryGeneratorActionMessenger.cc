//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "BxReadParameters.hh"
#include "BxGeneratorG4Gun.hh"
#include "BxLogger.hh"
#include "BxGeneratorMultiEvent.hh"
#include "BxGeneratorSolarNeutrino.hh"
#include "BxGeneratorSolarNeutrino2.hh"
/*#include "BxGeneratorCosmicRayMuons.hh"
#include "BxGeneratorCERNRays.hh"
#include "BxGeneratorNeutronsAtGS.hh"
#include "BxGeneratorLaser.hh"
#include "BxGeneratorTiming.hh"
*/
#include "BxGeneratorSterileAntiNu.hh"
#include "BxPrimaryGeneratorActionMessenger.hh"
#include "BxPrimaryGeneratorAction.hh"

/*#include "BxGeneratorGeneb.hh"
#include "BxGeneratorGenebDecay.hh"
#include "BxGeneratorG4Bx.hh"
*/
#include "BxGeneratorSCS.hh"
#include "BxGeneratorRDMDecayGun.hh"
#include "BxGeneratorAntiNeutrino.hh"
#include "BxGeneratorLi9He8.hh"
#include "BxGeneratorAmBeSource.hh"
#include "BxGeneratorNeutrinoC13.hh"
//#include "BxGeneratorSterileNeutrino.hh"

#include "BxGeneratorPositronium.hh"
#include "BxOutputVertex.hh"

using namespace std;

BxPrimaryGeneratorActionMessenger::BxPrimaryGeneratorActionMessenger(BxPrimaryGeneratorAction *act){
  fGeneratorPrimary = act;
  fDirectory = new G4UIdirectory("/bx/generator/");
  fDirectory->SetGuidance("Control commands for generators:");
  fDirectory->SetGuidance("/bx/generator/select: Select generator.");

  // /Bx/generator/select command
  fSelectCmd = new G4UIcmdWithAString("/bx/generator/select", this);
  fSelectCmd->SetGuidance("Selects generator for events.");
 /* fSelectCmd->SetGuidance("Options are:");
  fSelectCmd->SetGuidance("G4gun: Standard G4 gun.");
  fSelectCmd->SetGuidance("CosmicRayMuons: cosmic-ray muons generator");
  fSelectCmd->SetGuidance("CERNRays: CERN muons generator");  
  G4String candidates 
  = "G4Gun Multi SNU SNU2 CosmicRayMuons NeutronsAtGS CERNRays Laser Timing Geneb GenebDecay G4Bx RDM SCS AntiNu AmBe NeutrinoC13 Positronium Li9He8 Sterile";
  */
  G4String candidates="G4Gun RDM AntiNu AmBe SNU SNU2 Positronium NeutrinoC13 Li9He8 SCS Multi SterileAntiNu"; 
fSelectCmd->SetCandidates(candidates);
 }


BxPrimaryGeneratorActionMessenger::~BxPrimaryGeneratorActionMessenger() {
  delete fSelectCmd;
  delete fDirectory;
  //delete fGeneratorPrimary;

}


void BxPrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue) { 
   if(command == fSelectCmd) {   
   
     BxLog(routine) << "Generator: " << newValue << endlog;
    
     if(newValue == "G4Gun") {
       BxOutputVertex::Get()->SetGenerator(3);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorG4Gun);            
} else if (newValue == "Multi") {     
       BxOutputVertex::Get()->SetGenerator(11);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorMultiEvent);        
     }/* else if (newValue == "CosmicRayMuons") {     
       BxOutputVertex::Get()->SetGenerator(6);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorCosmicRayMuons);        
     } else if (newValue == "CERNRays") {     
       BxOutputVertex::Get()->SetGenerator(8);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorCERNRays);        
     } else if (newValue == "NeutronsAtGS") {     
      BxOutputVertex::Get()->SetGenerator(7);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorNeutronsAtGS);        
     } else if (newValue == "Laser") {     
       BxOutputVertex::Get()->SetGenerator(4);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorLaser);        
     } else if (newValue == "Timing") {     
       BxOutputVertex::Get()->SetGenerator(5);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorTiming);        
     } else if (newValue == "Geneb") {     
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorGeneb);     
       BxOutputVertex::Get()->SetGenebFlag(true);   
       BxOutputVertex::Get()->SetGenerator(0);   
     } else if (newValue == "GenebDecay") {     
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorGenebDecay);        
       BxOutputVertex::Get()->SetGenerator(1);   
     }*/ else if (newValue == "RDM") {     
       BxOutputVertex::Get()->SetGenerator(2);   
       BxReadParameters::Get()->SetRDMDecay(true);
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorRDMDecayGun);
     } else if (newValue == "SCS") {     
       BxOutputVertex::Get()->SetGenerator(9);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorSCS);
     } else if (newValue == "SNU") {     
       BxOutputVertex::Get()->SetGenerator(13);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorSolarNeutrino);
     } else if (newValue == "SNU2") {     
       BxOutputVertex::Get()->SetGenerator(10);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorSolarNeutrino2);
     } else if (newValue == "AntiNu") {     
       BxOutputVertex::Get()->SetGenerator(14);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorAntiNeutrino);
     } else if (newValue == "AmBe") {     
       BxOutputVertex::Get()->SetGenerator(12);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorAmBeSource);
     }/* else if (newValue == "G4Bx") {     
       BxOutputVertex::Get()->SetGenerator(15);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorG4Bx);
     }*/ else if (newValue == "NeutrinoC13") {     
       BxOutputVertex::Get()->SetGenerator(16);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorNeutrinoC13);
     } else if (newValue == "Positronium") {     
       BxOutputVertex::Get()->SetGenerator(17);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorPositronium);
     } else if (newValue == "Li9He8") {     
       BxOutputVertex::Get()->SetGenerator(18);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorLi9He8);
     } else if (newValue == "SterileAntiNu") {     
       BxOutputVertex::Get()->SetGenerator(19);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorSterileAntiNu);
     }/* else if (newValue == "Sterile") {     
       BxOutputVertex::Get()->SetGenerator(19);   
       fGeneratorPrimary->SetBxGenerator(new BxGeneratorSterileNeutrino);
     }
     */
   }

}
