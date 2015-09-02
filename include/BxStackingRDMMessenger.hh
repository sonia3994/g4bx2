//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef BxStackingRDMMessenger_h
#define BxStackingRDMMessenger_h 1

#include "BxIO.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4ParticleGun.hh"
#include "G4UIcmdWithAnInteger.hh"

class G4UIcmdWithADoubleAndUnit;
class BxStackingRDM;

using namespace std;

class BxStackingRDMMessenger: public G4UImessenger {

  public:
    BxStackingRDMMessenger(BxStackingRDM*  );
   ~BxStackingRDMMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    BxStackingRDM*                       stacking;   
    G4UIdirectory*                       fDirectory;
    G4UIcmdWithADoubleAndUnit*           fRangeCmd;
    G4UIcmdWithADoubleAndUnit*           fKillAngleCmd;
    G4UIcmdWithABool*                    fBulkCmd;
    G4UIcmdWithAnInteger*                fKillParticleCmd;
    G4UIcommand*                         fKillLEParticleCmd; 
   
    
};

#endif
