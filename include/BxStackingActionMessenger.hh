//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef BxStackingActionMessenger_h
#define BxStackingActionMessenger_h

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

#include "BxStackingAction.hh"

class G4UserStackingAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class BxStackingAction;
class G4StackManager;

///Messenger of BxStackingAction
class BxStackingActionMessenger: public G4UImessenger {

  public:
    BxStackingActionMessenger(BxStackingAction*);
   ~BxStackingActionMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    G4UIdirectory*              fDirectory;
    G4UIcmdWithAString*         fSelectCmd;  
    BxStackingAction*           fStacking;
};

#endif
