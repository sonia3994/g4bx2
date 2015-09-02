//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#include "BxPrimaryGeneratorAction.hh"
#include "BxPrimaryGeneratorActionMessenger.hh"



BxPrimaryGeneratorAction::BxPrimaryGeneratorAction() {
   fMessenger = new BxPrimaryGeneratorActionMessenger(this);
}

BxPrimaryGeneratorAction::~BxPrimaryGeneratorAction(){
   delete fMessenger;
}

void BxPrimaryGeneratorAction::GeneratePrimaries(G4Event *event) {
  
  generator->BxGeneratePrimaries(event);  
  

}
