#ifndef BxGeneratorSupernovaAntiNuMessenger_h
#define BxGeneratorSupernovaAntiNuMessenger_h 1
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
class BxGeneratorSupernovaAntiNu;

using namespace std;

class BxGeneratorSupernovaAntiNuMessenger: public G4UImessenger {

  public:
    	BxGeneratorSupernovaAntiNuMessenger(BxGeneratorSupernovaAntiNu*  );
   	~BxGeneratorSupernovaAntiNuMessenger();

  	void SetNewValue(G4UIcommand*, G4String);

  private:
    	BxGeneratorSupernovaAntiNu*              generator;   
    	G4UIdirectory*                           fDirectory;

        G4UIcmdWith3VectorAndUnit*           fPositionCmd;
        G4UIcmdWithAString*                  fNeutrinoCmd ;
        G4UIcmdWithAnInteger*                fGenInScintCmd;

};
#endif
