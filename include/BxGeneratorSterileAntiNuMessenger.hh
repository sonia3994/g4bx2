#ifndef BxGeneratorSterileAntiNuMessenger_h
#define BxGeneratorSterileAntiNuMessenger_h 1
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
#include "BxGeneratorSterileAntiNu.hh"

class BxPrimaryGeneratorAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class BxGeneratorSterileAntiNu;

using namespace std;

class BxGeneratorSterileAntiNuMessenger: public G4UImessenger {

  public:
    BxGeneratorSterileAntiNuMessenger(BxGeneratorSterileAntiNu*  );
   ~BxGeneratorSterileAntiNuMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    BxGeneratorSterileAntiNu*              generator;   
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
    /// New:
   G4UIcmdWithADouble*                  fDM2Cmd;
   G4UIcmdWithADouble*                  fSIN2T2Cmd ;
   G4UIcmdWithADouble*                  fCOEFFICIENTACmd;
   G4UIcmdWithADouble*                  fCOEFFICIENTBCmd ;
   G4UIcmdWithADouble*                  fCOEFFICIENTCCmd ;
};

#endif
/*
 * $Log: BxGeneratorSterileAntiNuMessenger.hh,v $
 * Revision 1.2  2015/02/25 18:17:27  acaminata
 * Shape factor added
 *
 * Revision 1.1  2015/02/12 14:32:07  acaminata
 * Sterile generator added
 *
 * Revision 1.7  2007-11-12 12:09:45  dfranco
 * added to g4gun, the following  energy distributions:
 * Lin (linear), Pow (power-law), Exp (exponential), Gauss (gaussian),
 * Brem (bremsstrahlung), BBody (black-body), Cdg (cosmic diffuse gamma-ray)
 *
 * Revision 1.6  2007-03-22 14:48:49  dfranco
 * Stable version
 *
 */
