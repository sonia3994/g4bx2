// --------------------------------------------------------------------------//
/** 
 * AUTHOR: I.Machulin
 * CONTACT: machulin@lngs.infn.it
* Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#ifndef BxGeneratorAntiNeutrinoMessenger_h
#define BxGeneratorAntiNeutrinoMessenger_h 1
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
class BxGeneratorAntiNeutrino;

using namespace std;
/**Geo - to simulate both positrons and neutrons of Geo origin
*
*Reactor - to simulate both positrons and neutrons of Reactor origin
*
*GeoPositron - to simulate only positrons of Geo origin
*
*ReactorPositron - to simulate only positrons of Reactor origin
*
*neutron - to simulate only neutrons
*
*Flat - to simulate flat positron distribution
*
*Sun - to simulate antineutrino from Sun (supposed to have the same spectre as neutrino)
*/

///Messenger of BxGeneratorAntiNeutrino
class BxGeneratorAntiNeutrinoMessenger: public G4UImessenger {

  public:
    BxGeneratorAntiNeutrinoMessenger(BxGeneratorAntiNeutrino*  );
   ~BxGeneratorAntiNeutrinoMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    BxGeneratorAntiNeutrino*                    generator;   
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
    G4UIcmdWithAString*                  fNeutrinoCmd ;

    G4UIcmdWithAnInteger*                fNumberOfStepCmd;
};

#endif
