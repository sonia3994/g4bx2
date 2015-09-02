//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef BxStackingEBMessenger_h
#define BxStackingEBMessenger_h 1

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
class BxStackingEB;

using namespace std;

class BxStackingEBMessenger: public G4UImessenger {

  public:
    BxStackingEBMessenger(BxStackingEB*  );
   ~BxStackingEBMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
/*    BxStackingEB*                       stacking;   
    G4UIdirectory*                       fDirectory;
    G4UIcmdWithADoubleAndUnit*           fRangeCmd;
    G4UIcmdWithADoubleAndUnit*           fKillAngleCmd;
    G4UIcmdWithABool*                    fBulkCmd;
    G4UIcmdWithAnInteger*                fKillParticleCmd;
    G4UIcommand*                         fKillLEParticleCmd; 
    G4UIcmdWithAString*                  fVesselSupDef; 
    G4UIcmdWithAString*                  fBufferDef;
    G4UIcmdWithADoubleAndUnit*           fVesselThickness;
    G4UIcmdWithADoubleAndUnit*           fRhoCut; 
    G4UIcmdWithADoubleAndUnit*           fZCut;
  */  
};

#endif
