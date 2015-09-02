//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014

#ifndef BxVStackingAction_h
#define BxVStackingAction_h 1
#include "globals.hh"
//#include "G4ClassificationOfNewTrack.hh"
#include "G4UserStackingAction.hh"
#include "BxStackingAction.hh"
#include "G4UImanager.hh"
#include "G4StackingMessenger.hh"

class G4StackManager ;
/**
 * Purely virtual class to manage the Stacking Actions
 */

class BxVStackingAction  {
  public:

      BxVStackingAction();
  
      virtual ~BxVStackingAction();

  public: // virtual methods which are then inherited by other stacking classes

      virtual G4ClassificationOfNewTrack  BxClassifyNewTrack(const G4Track* aTrack) = 0;

      virtual void BxNewStage() = 0;

      virtual void BxPrepareNewEvent() = 0;

 
  public: // concrete methods

      inline void BxStackAbort()                 { UImanager->ApplyCommand("/event/abort");         }
      inline void BxStackClearAll()              { UImanager->ApplyCommand("/event/stack/clear 2"); }
      inline void BxStackClearUrgentAndWaiting() { UImanager->ApplyCommand("/event/stack/clear 1"); }
      inline void BxStackClearWaiting()          { UImanager->ApplyCommand("/event/stack/clear 0"); }
      inline void BxStackClearUrgent()           { UImanager->ApplyCommand("/event/stack/clear -1");}
      inline void BxStackClearPostponed()        { UImanager->ApplyCommand("/event/stack/clear -2");}
      inline void BxStackCheckStatus()           { UImanager->ApplyCommand("/event/stack/status");}
      inline void BxStackReClassify()            { UImanager->ApplyCommand("/event/stack/clear -3");}
    

   private:
   
      //G4StackManager *fManager;  
      G4UImanager * UImanager ;
     

};

#endif

