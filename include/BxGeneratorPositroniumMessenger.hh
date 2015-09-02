// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#ifndef BxGeneratorPositroniumMessenger_h
#define BxGeneratorPositroniumMessenger_h 1
#include "BxIO.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
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
class BxGeneratorPositronium;

using namespace std;
///Messenger for BxGeneratorPositronium
class BxGeneratorPositroniumMessenger: public G4UImessenger {

  public:
    BxGeneratorPositroniumMessenger(BxGeneratorPositronium*  );
   ~BxGeneratorPositroniumMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    BxGeneratorPositronium*              generator;   
    G4UIdirectory*                       fDirectory;
    G4UIcmdWith3VectorAndUnit*           fPositionCmd; 
    G4UIcmdWith3Vector*                  fDirectionCmd;
    G4UIcmdWithADoubleAndUnit*           fSphereBulkCmd;
    G4UIcmdWithADoubleAndUnit*           fSphereSurfCmd;
    G4UIcmdWith3VectorAndUnit*           fSphereCentreCmd; 
    G4UIcmdWithADoubleAndUnit*           fSphereBulkRadMinCmd;
    G4UIcmdWithABool*                    fVesselCmd;
    G4UIcmdWithABool*                    fBulkCmd;
    G4UIcmdWithABool*                    fBufferCmd;
    G4UIcmdWithAString*                  fConfineCmd ;
    G4UIcmdWithADoubleAndUnit*           fEnergyCmd;
    G4UIcmdWithADouble*                  fMeanlifeCmd;
    G4UIcmdWithADouble*                  fProbabilityCmd;
    G4UIcmdWithABool*                    fSpectrumCmd;
    
};

#endif
