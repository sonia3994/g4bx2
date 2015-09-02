//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef BXTRACKINGACTION_H
#define BXTRACKINGACTION_H

#include "globals.hh"
#include "G4UserTrackingAction.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"

///It manages the tracking of each generated particle
class BxTrackingAction : public G4UserTrackingAction {

  public:

	BxTrackingAction();
	virtual ~BxTrackingAction(){};
	virtual void PostUserTrackingAction(const G4Track*);

  private:

	G4MaterialPropertiesTable 	*PMTQETable;
  	G4MaterialPropertyVector*	PMTQEVector;
 	void				ReadPMTQE();
};

#endif
