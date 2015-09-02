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
// Created by D. Franco
// Revised by A. Caminata and S. Marcocci, Sept. 2014

#ifndef BxStackingEB_h
#define BxStackingEB_h 1
#include "G4ClassificationOfNewTrack.hh"
#include "G4UserStackingAction.hh"
#include "BxVStackingAction.hh"
#include "BxStackingAction.hh"

class G4StackManager;
class G4Track;
class BxStackingEBMessenger;

/**
*  This is the base class of one of the user's optional action classes.
* This class gives the hooks for G4StackManager which controls the stacks
* of G4Track objects.
*/

class BxStackingEB : public BxVStackingAction
{
  public:
      BxStackingEB();
      virtual ~BxStackingEB();

  public: // with description
//---------------------------------------------------------------
// vitual methods to be implemented by user
//---------------------------------------------------------------
//
      virtual G4ClassificationOfNewTrack 
        BxClassifyNewTrack(const G4Track* aTrack);
//---------------------------------------------------------------
//
/**
*    This method is called by G4StackManager when the urgentStack
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
*/
      virtual void BxNewStage();
//---------------------------------------------------------------
//
/**
*    This method is called by G4StackManager at the begining of
*  each event.
*    Be careful that the urgentStack and the waitingStack of 
*  G4StackManager are empty at this moment, because this method
*  is called before accepting primary particles. Also, note that
*  the postponeStack of G4StackManager may have some postponed
*  tracks.
*/
      virtual void BxPrepareNewEvent();
//---------------------------------------------------------------
    /*  void   SetIsReclassify(G4bool val) { IsReclassify = val;  }      
      G4bool GetIsReclassify()           { return IsReclassify; }      
      
      void   SetfSupDef(G4bool val) { fSupDef = val;  }      
      G4bool GetfSupDef()           { return fSupDef; } 
      
      void   SetfVesThick(G4double val) { fVthick = val;  } 
      
      void   SetfRhoCut(G4double val) { frhoCut = val;  }  
      void   SetffRhocut(G4bool val)   { ffRhocut = val;  }     
      G4bool GetffRhocut()            { return ffRhocut; }      
      
      void   SetfZCut(G4double val) { fzCut = val;  }   
      void   SetffZcut(G4bool val)  { ffZcut = val;  }       
      G4bool GetffZcut()            { return ffZcut; }      
      
      void   SetfBufDef(G4bool val) { fBufDef = val;  }      
      G4bool GetfBufDef()           { return fBufDef; }      
           
      void   SetRadius(G4double radius) { fRadius = radius;  }      
 
      void   KillParticles( G4int val ) { fPDGToBeKilled.push_back(val) ;}    
      void   KillLEParticles( G4int val, G4double ene ) { fLEPDGToBeKilled.push_back(val); fLEnergyToBeKilled.push_back(ene);}    
      void   SetKillingAngle( G4double val) {fKillOuterDirected = true ; fAngle = val ; } 
  */ 
   
    private: 
     /* BxStackingEBMessenger  *fMessenger ;
      
      G4bool IsReclassify ;
      G4int fCounter, stage;
      G4bool IsValid ;
      G4int fDaughters;  
      G4int fAlpha, fBeta, fGamma;
      G4double fRadius;
      G4bool IsTheFirstDaughter;
      vector<G4int>     fPDGToBeKilled ;
      vector<G4int>     fLEPDGToBeKilled;
      vector<G4double>  fLEnergyToBeKilled;
      G4bool            fKillOuterDirected;
      G4double          fAngle;
      G4double          rr, zz;
      G4bool            fKillVessel; 
      G4bool            fSupDef;
      G4bool            fBufDef;
      G4bool            isRead;
      G4double          RhoV[200], ZV[200];
      G4int             RZDim;
      G4double          fVthick;
      G4double          frhoCut; 
      G4double          fzCut; 
      G4bool            VSloaded;
      G4bool            ffZcut;
      G4bool            ffRhocut;
     
*/
};

#endif

