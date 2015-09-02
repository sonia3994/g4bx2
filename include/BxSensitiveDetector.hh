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

#ifndef BxSensitiveDetector_h
#define BxSensitiveDetector_h 1


#include "G4VHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4VReadOutGeometry.hh"
#include "G4TouchableHistory.hh"
#include "G4CollectionNameVector.hh"
#include "G4VSDFilter.hh"
#include "G4VSensitiveDetector.hh"

#include  <vector>

/**
*  This is the abstract base class of the sensitive detector. The user's
* sensitive detector which generates hits must be derived from this
* class.
*  In the derived class constructor, name(s) of hits collection(s) which
* are made by the sensitive detector must be set to "collectionName" string
* vector.
 */

class BxSensitiveDetector : public G4VSensitiveDetector
{

  public: // with description
      BxSensitiveDetector(G4String name);
      /// Constructors. The user's concrete class must use one of these constructors
      /// by the constructor initializer of the derived class. The name of
      /// the sensitive detector must be unique.

  public:
      virtual ~BxSensitiveDetector();

 
  public: // with description
      ///  These two methods are invoked at the begining and at the end of each
      /// event. The hits collection(s) created by this sensitive detector must
      /// be set to the G4HCofThisEvent object at one of these two methods.
      virtual void Initialize(G4HCofThisEvent*);
      virtual void EndOfEvent(G4HCofThisEvent*);
      ///  This method is invoked if the event abortion is occured. Hits collections
      /// created but not beibg set to G4HCofThisEvent at the event should be deleted.
      /// Collection(s) which have already set to G4HCofThisEvent will be deleted 
      /// automatically.
      virtual void clear();


  protected: // with description
      /**  The user MUST implement this method for generating hit(s) using the 
      * information of G4Step object. Note that the volume and the position
      * information is kept in PreStepPoint of G4Step.
      *  Be aware that this method is a protected method and it sill be invoked 
      * by Hit() method of Base class after Readout geometry associated to the
      * sensitive detector is handled.
      *  "ROhist" will be given only is a Readout geometry is defined to this
      * sensitive detector. The G4TouchableHistory object of the tracking geometry
      * is stored in the PreStepPoint object of G4Step.
      *virtual G4int GetCollectionID(G4int i);
      *  This is a utility method which returns the hits collection ID of the
      * "i"-th collection. "i" is the order (starting with zero) of the collection
      * whose name is stored to the collectionName protected vector.
 */
      virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist) ;
  
  protected:
      G4String SensitiveDetectorName; // detector name
      
  private:
  
      std::vector<G4int>         trackID;
  
};




#endif

