//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef BxGeneratorG4GunMessenger_h
#define BxGeneratorG4GunMessenger_h 1
#include "BxIO.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4ParticleGun.hh"

class BxPrimaryGeneratorAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class BxGeneratorG4Gun;

using namespace std;

///Messenger of BxGeneratorG4Gun
class BxGeneratorG4GunMessenger: public G4UImessenger {

  public:
///ctor
    BxGeneratorG4GunMessenger(BxGeneratorG4Gun*  );
///dtor
   ~BxGeneratorG4GunMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    BxGeneratorG4Gun*                    generator;   
    G4UIdirectory*                       fDirectory;
    G4UIcmdWith3VectorAndUnit*           fPositionCmd; 
    G4UIcmdWith3Vector*                  fDirectionCmd;
    G4UIcmdWithADoubleAndUnit*           fEnergyCmd;
    G4UIcmdWithAString*                  fParticleCmd;
    G4UIcmdWithADoubleAndUnit*           fSphereBulkCmd;
    G4UIcmdWithADoubleAndUnit*           fSphereSurfCmd;
    G4UIcmdWith3VectorAndUnit*           fSphereCentreCmd; 
    G4UIcmdWithADoubleAndUnit*           fSphereBulkRadMinCmd;
    G4UIcmdWithABool*                    fVesselCmd;
    G4UIcmdWithABool*                    fBulkCmd;
    G4UIcmdWithABool*                    fBufferCmd;
    G4UIcmdWithAString*                  fConfineCmd ;
    G4UIcmdWithAString*                  fEnergyDisTypeCmd;
    G4UIcmdWithADoubleAndUnit*           fSetEminCmd;
    G4UIcmdWithADoubleAndUnit*           fSetEmaxCmd;   
    G4UIcmdWithADouble*                  fSetAlphaCmd;  
    G4UIcmdWithADoubleAndUnit*           fSetTempCmd;   
    G4UIcmdWithADoubleAndUnit*           fSetEzeroCmd;  
    G4UIcmdWithADouble*                  fSetGradientCmd;
    G4UIcmdWithADouble*                  fSetInterCeptCmd;
};

#endif
