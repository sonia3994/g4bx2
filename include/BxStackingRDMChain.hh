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

#ifndef BxStackingRDMChain_h
#define BxStackingRDMChain_h 1
#include "G4ClassificationOfNewTrack.hh"
#include "G4UserStackingAction.hh"
#include "G4ParticleChange.hh"
#include "BxVStackingAction.hh"
#include "BxStackingAction.hh"

class G4StackManager;
class G4Track;
class BxStackingRDMChainMessenger;

/**
*  This is the base class of one of the user's optional action classes.
* This class gives the hooks for G4StackManager which controls the stacks
* of G4Track objects.
*/

class BxStackingRDMChain : public BxVStackingAction
{
  public:
///costructor
      BxStackingRDMChain();
///destructor
      virtual ~BxStackingRDMChain();

  public: // with description
//---------------------------------------------------------------
// vitual methods to be implemented by user
//---------------------------------------------------------------
//
///It classifies the track: possible state fUrgent, fWaiting, fPostpone, fKill
      virtual G4ClassificationOfNewTrack 
        BxClassifyNewTrack(const G4Track* aTrack);
//
//    Reply G4ClassificationOfNewTrack determined by the
//  newly coming G4Track.
//
//    enum G4ClassificationOfNewTrack
//    {
//      fUrgent,    // put into the urgent stack
//      fWaiting,   // put into the waiting stack
//      fPostpone,  // postpone to the next event
//      fKill       // kill without stacking
//    };
//
//    The parent_ID of the track indicates the origin of it.
//                
//    G4int parent_ID = aTrack->get_parentID();
//   
//      parent_ID = 0 : primary particle
//                > 0 : secondary particle
//                < 0 : postponed from the previous event
//
//---------------------------------------------------------------


/**    This method is called by G4StackManager when the urgentStack
*  becomes empty and contents in the waitingStack are transtered
*  to the urgentStack.
*    Note that this method is not called at the begining of each
*  event, but "PrepareNewEvent" is called.
*
*    In case re-classification of the stacked tracks is needed,
*  use the following method to request to G4StackManager.
*
*    stackManager->ReClassify();
*
*  All of the stacked tracks in the waitingStack will be re-classified 
*  by "ClassifyNewTrack" method.
*    To abort current event, use the following method.
*
*    stackManager->clear();
*
*  Note that this way is valid and safe only for the case it is called
*  from this user class. The more global way of event abortion is
*
*    G4UImanager * UImanager = G4UImanager::GetUIpointer();
*    UImanager->ApplyCommand("/event/abort");
**/
      virtual void BxNewStage();
//---------------------------------------------------------------
  
/**
*   This method is called by G4StackManager at the begining of
*  each event.
*    Be careful that the urgentStack and the waitingStack of 
*  G4StackManager are empty at this moment, because this method
*  is called before accepting primary particles. Also, note that
*  the postponeStack of G4StackManager may have some postponed
*  tracks.
**/
    virtual void BxPrepareNewEvent();
//---------------------------------------------------------------
      void   SetIsReclassify(G4bool val) { IsReclassify = val;  }      
      G4bool GetIsReclassify()           { return IsReclassify; }      
      
      void   SetfSupDef(G4bool val) { fSupDef = val;  }      
      G4bool GetfSupDef()           { return fSupDef; } 
 
      ///Set the vessel thickness     
      void   SetfVesThick(G4double val) { fVthick = val;  }   

      ///FIXME used on cuts: specify
      void   SetfRhoCut(G4double val) { frhoCut = val;  }  
      void   SetffRhocut(G4bool val)  { ffRhocut = val;  }     
      G4bool GetffRhocut()            { return ffRhocut; }          

      ///FIXME used on cuts: specify
      void   SetfZCut(G4double val) { fzCut = val;  }  
      void   SetffZcut(G4bool val)  { ffZcut = val;  }  
      G4bool GetffZcut()            { return ffZcut; }          

      ///FIXME used on cuts: specify
      void   SetfBufDef(G4bool val) { fBufDef = val;  }      
      G4bool GetfBufDef()           { return fBufDef; }  

      ///Set the isotope max lifetime. if lifetime is greater than fMaxLifeTime, the decay chain is stopped     
      void   SetMaxLifeTime(G4double ltime) {  fMaxLifeTime = ltime;  }      

  private:

      BxStackingRDMChainMessenger  *fMessenger ;

      G4bool            IsReclassify ;
      G4bool            IsShort ;
      G4int             fCounter; 
      G4int             stage;
      G4bool            IsValid ;
      G4bool            isDaughter; 
      G4bool            isSecondDaughter; 
      G4bool            isFirst; 
      G4int             fAlpha;
      G4int             fBeta; 
      G4int             fGamma; 
      G4int             fDaughters;
      G4double          fRadius;
      ///Max lifetime, 10 days
      G4double          fMaxLifeTime;
      G4ParticleChange* changedParticle;
      G4int             fNSequence;
      G4bool            fSupDef;
      G4bool            fBufDef;
      G4bool            ffZcut;
      G4bool            ffRhocut;
      G4double          rr, zz;
      G4double          RhoV[200], ZV[200];
      G4double          fVthick; 
      ///rho cut value   
      G4double          frhoCut;
      ///Z cut value 
      G4double          fzCut; 
      G4bool            VSloaded;
      G4int             RZDim;
};

#endif
