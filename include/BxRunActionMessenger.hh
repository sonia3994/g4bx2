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

#ifndef BxRunActionMessenger_h
#define BxRunActionMessenger_h 1
#include "BxIO.hh"

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include <sstream>
#include "globals.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

class BxRunAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
///Messenger for BxRunAction
class BxRunActionMessenger: public G4UImessenger {

  public:
///Constructor
    BxRunActionMessenger(BxRunAction*);
///Destructor
   ~BxRunActionMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    BxRunAction*	       runAction;	    
///Autoseed on/off    
G4UIcmdWithABool*	       fSetAutoSeedCmd;     
///Sets random number generator seed
    G4UIcmdWithAnInteger*      fHEPRandomSeedCmd;   
///Name for output files 
   G4UIcmdWithAString*	       fSetFileNameCmd;     
///Set binary name
    G4UIcmdWithAString*	       fSetG4BxNameCmd;    
    G4UIcmdWithAString*	       fSetCommentCmd;   
    G4UIdirectory*	       fDirectory;	    
    G4UIcmdWithABool*	       fSetIsBinaryCmd;     
    G4UIcmdWithAnInteger*      fRunCmd; 	    
    G4UIcmdWithADoubleAndUnit* fRateCmd;
    G4UIcmdWithABool*	       fWriteDepositsCmd;
    G4UIcmdWithABool*	       fWriteIsotopesCmd;
    G4UIcmdWithABool*	       fWriteEBCmd;
};

#endif
