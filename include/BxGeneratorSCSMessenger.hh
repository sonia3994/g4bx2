//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef BxGeneratorSCSMessenger_h
#define BxGeneratorSCSMessenger_h 1
#include "BxIO.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4ParticleGun.hh"

class BxPrimaryGeneratorAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class BxGeneratorSCS;

using namespace std;
///Messenger of BxGeneratorSCS
/**With the command /bx/generator/scs/c14reference xxx it is possible to define the space factor from references in literature:
 * – 1: (1 - 0.38 * E) from F. P. Calaprice and B. R. Holstein, Nucl. Phys. A 273, 301 (1976);
 *
 * – 2: (1 - 0.72 * E) from G. Alimonti, G. Angloher, C. Arpesell, et al., Phys. Lett. B 422, 349 (1998);
 *
 * – 3: (1 - 1.179 * E) from H. Genz, G. Kuhner, A. Richter, and H. Behrens, Z. Phys. A 341, 9 (1991);
 *
 * – 4: (1 - 0.37 * E) from A. Garcia and B. A. Brown, Phys. Rev. C 52, 3416 (1995);
 *
 * – 5: (1 - 0.45 * E) from F. E. Weitfeldt, E. B. Norman, Y. D. Chan, et al., Phys. Rev. C 52, 1028 (1995);
 *
 * – 6: (1 + 1.1 * (E0 - E)) from B. Sur, E. B. Norman, K. T. Lesko, et al., Phys. Rev. Lett. 66 2444 (1991);
 *
 *– 7: (1 + 1.24 * (E0 - E)) from V. V. Kuzminov and N. J. Osetrova, Exp. Res. Meth. And Fac. 63 (2000) 1365;
*
– 8: (1 - 4.67*E + 3/E +2*E*E) from Ch. Sonntag, H. Rebel, B. Ribbat, et al., Lett. Nuovo Cimento 4, 717 (1970);
 */

class BxGeneratorSCSMessenger: public G4UImessenger {

  public:
    BxGeneratorSCSMessenger(BxGeneratorSCS*  );
   ~BxGeneratorSCSMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    BxGeneratorSCS*                      generator;   
    G4UIdirectory*                       fDirectory;
    G4UIcmdWith3VectorAndUnit*           fPositionCmd; 
    G4UIcmdWith3Vector*                  fDirectionCmd;
    //G4UIcmdWithADoubleAndUnit*           fEnergyCmd;
    G4UIcmdWithAString*                  fParticleCmd;
    G4UIcmdWithADoubleAndUnit*           fSphereBulkCmd;
    G4UIcmdWithADoubleAndUnit*           fSphereSurfCmd;
    G4UIcmdWith3VectorAndUnit*           fSphereCentreCmd; 
    G4UIcmdWithADoubleAndUnit*           fSphereBulkRadMinCmd;
    G4UIcmdWithABool*                    fVesselCmd;
    G4UIcmdWithABool*                    fBulkCmd;
    G4UIcmdWithABool*                    fBufferCmd;
    G4UIcmdWithAString*                  fConfineCmd ;

    G4UIcmdWithAnInteger*                fNumberOfBinsCmd ;
    G4UIcmdWithADouble*                  fShapeFactorCmd;
    G4UIcmdWithAnInteger*                fReferenceCmd ;
    G4UIcmdWithAnInteger*                fShapeFactorTermCmd ;
};

#endif
