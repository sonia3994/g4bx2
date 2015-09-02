//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef BxGeneratorMultiEventMessenger_h
#define BxGeneratorMultiEventMessenger_h 1
#include "BxIO.hh"
#include "BxOutputStructure.hh"

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
class BxGeneratorMultiEvent;

using namespace std;

class BxGeneratorMultiEventMessenger: public G4UImessenger {

  public:
    BxGeneratorMultiEventMessenger(BxGeneratorMultiEvent*  );
   ~BxGeneratorMultiEventMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    BxGeneratorMultiEvent*               generator;   
    G4UIdirectory*                       fDirectory;
    G4UIcmdWith3VectorAndUnit*           fPositionCmd; 
    //G4UIcmdWith3Vector*                  fDirectionCmd;
    //G4UIcmdWithADoubleAndUnit*           fEnergyCmd;
    //G4UIcmdWithAString*                  fParticleCmd;
    G4UIcmdWithADoubleAndUnit*           fSphereBulkCmd;
    G4UIcmdWithADoubleAndUnit*           fSphereSurfCmd;
    G4UIcmdWith3VectorAndUnit*           fSphereCentreCmd; 
    G4UIcmdWithADoubleAndUnit*           fSphereBulkRadMinCmd;
    G4UIcmdWithABool*                    fVesselCmd;
    G4UIcmdWithABool*                    fBulkCmd;
    G4UIcmdWithABool*                    fBufferCmd;
    G4UIcmdWithAString* 		 fConfineCmd;    
    G4UIcmdWithABool*                    fOnlyGammaCmd;
    G4UIcmdWithAnInteger*                fSkipEventCmd;
    G4UIcmdWithAString*                  fPDGCmd;

    MultiEvent                           fMultiParticle;
    MultiEvent&                          ConvertFromString(string);

};

#endif
