//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//Created by D. Franco 
//Revised by A. Caminata and S. Marcocci, Sept. 2014

#ifndef BxPhysicsListMessenger_h
#define BxPhysicsListMessenger_h 1


#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIdirectory.hh"

using namespace std;

class BxPhysicsList;
///Messenger of BxPhysicsList; It allows to activate or deactivate certain interactions during the simulation.
class BxPhysicsListMessenger: public G4UImessenger {

  public:
    BxPhysicsListMessenger(BxPhysicsList*);
   ~BxPhysicsListMessenger();

    void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    
    BxPhysicsList*        fPhysicsList;
    G4UIcmdWithABool*     fOPCmd;
    G4UIcmdWithABool*     fOPcherCmd;    
    G4UIcmdWithABool*     fOPscintCmd;    
    G4UIcmdWithAString*   fHadronicCmd;    
    G4UIcmdWithABool*     fAltHadronicCmd;    
    G4UIdirectory*        fDirectory;
    G4UIcmdWithAnInteger* fMaxNumPhotonCmd;
    G4UIcmdWithAnInteger* fMaxParentIdForQuenching;
    G4UIcmdWithADouble*   fReduceLightYieldCmd; 
    G4UIcmdWithADouble*   fReduceCerenkovLightCmd; 
    G4UIcmdWithABool*     fARCmd;
    G4UIcmdWithABool*     fRayleighCmd;
    G4UIcmdWithABool*     fReemission;
    G4UIcmdWithABool*     fDecayCmd;
    G4UIcmdWithADoubleAndUnit*  fCutCmd;
    G4UIcmdWithADoubleAndUnit*  fElectronCutCmd;
    G4UIcmdWithADoubleAndUnit*  fGammaCutCmd;
    G4UIcmdWithAnInteger* fLightYieldCmd;
    G4UIcmdWithABool*     fNewEMCmd;    
    G4UIcmdWithABool*     fNucIntCmd;    
    G4UIcmdWithABool*     fIonPhysCmd;    
    G4UIcmdWithABool*     fTrackSecondariesFirstCmd;    
};

#endif
