#ifndef BxStackingRDM_h
#define BxStackingRDM_h 1
#include "G4ClassificationOfNewTrack.hh"
#include "G4UserStackingAction.hh"
#include "BxVStackingAction.hh"
#include "BxStackingAction.hh"

class G4StackManager;
class G4Track;
class BxStackingRDMMessenger;
// class description:
//
//  This is the base class of one of the user's optional action classes.
// This class gives the hooks for G4StackManager which controls the stacks
// of G4Track objects.
//

class BxStackingRDM : public BxVStackingAction
{
  public:
      BxStackingRDM();
      virtual ~BxStackingRDM();

  public: // with description
//---------------------------------------------------------------
// vitual methods to be implemented by user
//---------------------------------------------------------------
//
/**
*    Reply G4ClassificationOfNewTrack determined by the
*  newly coming G4Track.
*
*    enum G4ClassificationOfNewTrack
*    {
*      fUrgent,    // put into the urgent stack
*      fWaiting,   // put into the waiting stack
*      fPostpone,  // postpone to the next event
*      fKill       // kill without stacking
*    };
*
*    The parent_ID of the track indicates the origin of it.
*                
*    G4int parent_ID = aTrack->get_parentID();
*   
*      parent_ID = 0 : primary particle
*                > 0 : secondary particle
*                < 0 : postponed from the previous event
**/
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
*
**/
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
**/
      virtual void BxPrepareNewEvent();
//---------------------------------------------------------------
///Set IsReclassify
      void   SetIsReclassify(G4bool val) { IsReclassify = val;  }      
///Get IsReclassify
      G4bool GetIsReclassify()           { return IsReclassify; }      
      
   ///Kill aprticle 
      void   KillParticles( G4int val ) { fPDGToBeKilled.push_back(val) ;}    
///Kill particle by energy 
      void   KillLEParticles( G4int val, G4double ene ) { fLEPDGToBeKilled.push_back(val); fLEnergyToBeKilled.push_back(ene);}    
///Kill particle by angle
      void   SetKillingAngle( G4double val) {fKillOuterDirected = true ; fAngle = val ; }    
      
    private:
    ///The messenger  
      BxStackingRDMMessenger  *fMessenger ;
      
      G4bool IsReclassify ;
      G4int fCounter, stage;
      G4bool IsValid ;
      G4int fDaughters;
///FIXME never used?  
      G4int fAlpha, fBeta, fGamma;
///Fiducial volume radius
      G4double fRadius;
///Boolean variable
      //G4bool IsTheFirstDaughter;
      vector<G4int>     fPDGToBeKilled ;
      vector<G4int>     fLEPDGToBeKilled;
      vector<G4double>  fLEnergyToBeKilled;
///kill if directed outside the detector;
      G4bool            fKillOuterDirected;
      G4double          fAngle;
     
     
      
};

#endif
/*
 * $Log: BxStackingRDM.hh,v $
 * Revision 1.2  2014/12/03 15:55:28  marcocci
 * Vessels' EndCaps added
 *
 * Revision 1.1  2014/10/20 14:53:08  marcocci
 * g4bx2 added to repository
 *
 * Revision 1.2  2014/06/03 13:30:23  acaminata
 * *** empty log message ***
 *
 * Revision 1.1  2014/05/03 10:33:31  marcocci
 * CMakeLists.txt
 *
 * Revision 1.14  2011/07/13 13:57:08  buizza
 * moved option for ext bkg from StackingRDM to StackingEB
 *
 * Revision 1.13  2011-01-20 17:07:35  buizza
 * added some minor modification for z and rho cut
 *
 * Revision 1.11  2010-10-13 17:39:55  buizza
 * Adding the possibility to simulate events IN deformed buffers
 *
 * Revision 1.10  2010-10-13 13:44:41  buizza
 * Adding the possibility to simulate background ON deformed vessels
 *
 * Revision 1.9  2009-10-16 15:41:41  davini
 * float, double, int -> G4float G4double G4int;
 *
 * Revision 1.8  2008-04-10 15:14:01  dfranco
 * added the possibility to kill particles when the angle between direction and position (oriented toward the center) vectors is larger than the user defined one. The command is: /bx/stack/rdm/killAngle xxx.
 *
 * Revision 1.7  2008-04-10 13:06:56  dfranco
 * Added the possibility to kill particles with energy below a certain threshold in keV. The command is /bx/stack/rdm/killLE xxx yyy. Look in the manual for a more detail explanation
 *
 * Revision 1.6  2008-03-14 11:23:18  dfranco
 * added a method to kill specific particle types in a decay.
 *
 * Revision 1.5  2007-03-30 15:58:20  dfranco
 * Added the cvs log at the end of file
 *
 */

