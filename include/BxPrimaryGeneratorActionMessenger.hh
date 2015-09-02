//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014

#ifndef BxPrimaryGeneratorActionMessenger_h
#define BxPrimaryGeneratorActionMessenger_h


#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

class BxPrimaryGeneratorAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;

///Messenger for BxPrimaryGeneratorAction
class BxPrimaryGeneratorActionMessenger: public G4UImessenger {

  public:
    BxPrimaryGeneratorActionMessenger(BxPrimaryGeneratorAction*);
   ~BxPrimaryGeneratorActionMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    BxPrimaryGeneratorAction*   fGeneratorPrimary;   
    G4UIdirectory*              fDirectory;
    G4UIcmdWithAString*         fSelectCmd;  
    
};

#endif
