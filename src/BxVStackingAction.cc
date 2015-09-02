//---------------------------------------------------------------------------//
/**
 *
 * Pure virtual base class for Bx StackingActions. 
 * 
 */
// End class description
/** 
 * AUTHOR: davide.franco@mi.infn.it
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
 */
// --------------------------------------------------------------------------//

#include "BxVStackingAction.hh"
#include "G4StackManager.hh"

//---------------------------------------------------------------------------//

BxVStackingAction::BxVStackingAction() {
   UImanager = G4UImanager::GetUIpointer();
}

//---------------------------------------------------------------------------//


BxVStackingAction::~BxVStackingAction(){}
