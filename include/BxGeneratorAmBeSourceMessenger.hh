// --------------------------------------------------------------------------//
/** 
 * AUTHOR: I.Machulin
 * CONTACT: machulin@lngs.infn.it
* Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#ifndef BxGeneratorAmBeSourceMessenger_h
#define BxGeneratorAmBeSourceMessenger_h 1
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
class BxGeneratorAmBeSource;

using namespace std;
 /** all - to simulate the full source
 *
 *  neutron0G - to simulate neutrons witout gamma
 *
 *   neutron1G -to simulate neutrons with 1 gamma (group n1)
 *
 *   neutron2G -to simulate neutrons of group n2
*/

///Messenger of BxGeneratorAmBeSource
class BxGeneratorAmBeSourceMessenger: public G4UImessenger {

  public:
    BxGeneratorAmBeSourceMessenger(BxGeneratorAmBeSource*  );
   ~BxGeneratorAmBeSourceMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    BxGeneratorAmBeSource*                    generator;   
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
