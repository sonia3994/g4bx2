//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef BxPrimaryGeneratorAction_h
#define BxPrimaryGeneratorAction_h
#include "G4VUserPrimaryGeneratorAction.hh"
#include "BxVGenerator.hh"
#include <iostream>
using namespace std;
class G4Event;
class BxPrimaryGeneratorActionMessenger;


///Base class for the generator
class BxPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {							   
  public:						   

  	  BxPrimaryGeneratorAction();			   
  	  ~BxPrimaryGeneratorAction();  		   

  public:
         virtual void GeneratePrimaries(G4Event*);
         void SetBxGenerator(BxVGenerator *gen) {generator = gen;}
	 


  private:
         BxPrimaryGeneratorActionMessenger* fMessenger;	
	 BxVGenerator*                      generator ;
	 
};
#endif 
