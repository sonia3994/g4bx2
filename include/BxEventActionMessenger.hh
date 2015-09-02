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
//
//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef BxEventActionMessenger_h
#define BxEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class BxEventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
///Messenger for BxEventAction
class BxEventActionMessenger: public G4UImessenger
{
  public:
    BxEventActionMessenger(BxEventAction*);
   ~BxEventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:

    BxEventAction*              eventAction;   
///It selects which tracks to draw: none, charged, neutral, all
    G4UIcmdWithAString*         fDrawCmd;
///Print event statistics if the event is multiple of
    G4UIcmdWithAnInteger*       fPrintCmd;   
    G4UIcmdWithAnInteger*       fCounterCmd;
///Discard events with npe < 
    G4UIcmdWithAnInteger*       fEmptyCmd;
   ///Set Verbosity
 G4UIcmdWithAnInteger*       fVerbosityCmd; 
///Discard events without any deposit with r < 
    G4UIcmdWithADoubleAndUnit*  fSetDepoRadiusCmd ;
///Discard events with visible energy < 
    G4UIcmdWithADoubleAndUnit*  fSetVisEneCmd ;
///Kill particles with global time > 
    G4UIcmdWithADoubleAndUnit*  fSetTimeCutCmd ;
///If false the output is written at 32bit, if true writes at 64bit
    G4UIcmdWithABool*           fBitCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
