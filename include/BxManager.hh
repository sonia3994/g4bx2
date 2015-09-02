//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef _BXMANAGER_HH
#define _BXMANAGER_HH

#include "G4RunManager.hh"
//---------------------------------------------------------------------------//

class  BxManagerMessenger;


//---------------------------------------------------------------------------//

class BxManager : public G4RunManager 
/**Class BxManager inherits from class G4RunManager. It is friend with BxManagerMessenger.
* This is the concrete class of the G4RunManager
*/

  {
private:
  ///Default Constructor
  BxManager();

public:
///Singleton
  static BxManager* Get();
     ///Destructor
  virtual ~BxManager();

  //private  members
private:
   static BxManager *Manager;
 /// Messenger to select the log severity level.
  BxManagerMessenger              *fBxMessenger;
};
#endif
