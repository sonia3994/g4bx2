//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef BxStackingRDMChainMessenger_h
#define BxStackingRDMChainMessenger_h 1
#include "BxIO.hh"

#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4ParticleGun.hh"


class G4UIcmdWithADoubleAndUnit;
class BxStackingRDMChain;

using namespace std;

class BxStackingRDMChainMessenger: public G4UImessenger {

  public:
    BxStackingRDMChainMessenger(BxStackingRDMChain*  );
   ~BxStackingRDMChainMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    BxStackingRDMChain*                       stacking;   
    G4UIdirectory*                       fDirectory;
    G4UIcmdWithADoubleAndUnit*           fLifeTimeCmd;
    G4UIcmdWithAString*                  fVesselSupDef; 
    G4UIcmdWithAString*                  fBufferDef;
    G4UIcmdWithADoubleAndUnit*           fVesselThickness;
    G4UIcmdWithADoubleAndUnit*           fRhoCut; 
    G4UIcmdWithADoubleAndUnit*           fZCut;
};

#endif
