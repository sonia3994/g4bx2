//Created by I. Machulin
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef BxSteppingAction_h
#define BxSteppingAction_h 1

#include "G4UserSteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
///It manages what happens at each step of the simulation
class BxSteppingAction : public G4UserSteppingAction
{
  public:
    BxSteppingAction();
   ~BxSteppingAction(){};

    void UserSteppingAction(const G4Step*);
  private:
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
