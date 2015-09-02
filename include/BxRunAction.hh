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

#ifndef BxRunAction_h
#define BxRunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UImanager.hh"

class G4Timer;
class G4Run;
class BxRunActionMessenger;
/**
 * It contains the actions to be done at the beginning and at the end of each run
*/

 class BxRunAction : public G4UserRunAction
{
  public:

    BxRunAction();
    virtual ~BxRunAction();

  public:

    virtual void BeginOfRunAction(const G4Run* );
    virtual void EndOfRunAction(const G4Run* aRun);
    
  public:
    inline void SetAutoSeed (const G4bool val)    {autoSeed    = val;}

  private:

    BxRunActionMessenger* runMessenger;
    G4bool                autoSeed;
    G4Timer*              timer;

};

#endif 
