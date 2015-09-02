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

#ifndef BxManagerMessenger_h
#define BxManagerMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "globals.hh"
#include "G4UIcmdWithAString.hh"

using namespace std;

class BxManager;
class G4UIcmdWithAString;

///class BxManagerMessenger inherits from G4UImessenger 
class BxManagerMessenger: public G4UImessenger {

  public:
///Constructor
    BxManagerMessenger(BxManager*);
///Destructor  
 ~BxManagerMessenger();
///It is called when the command "/bxlog" is written in batch mode or read in a .mac file.
  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    G4UIcmdWithAString     *fBxLogCmd;
};

#endif
