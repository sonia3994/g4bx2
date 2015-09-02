//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#include "BxDataCollection.hh"
#include <iostream>
#include "stdlib.h"

using namespace std;

BxDataCollection* BxDataCollection::me=0;

BxDataCollection* BxDataCollection::Get(){
	if(!me)
		me=new BxDataCollection();

return me;
}


BxDataCollection::BxDataCollection() {
//For standart Event at one point
	fEID = 0;
	fEPosition = G4ThreeVector(0,0,0);
	fEDirection = G4ThreeVector(0,0,0);
/*        fShootDist = 0.; */
	fEParticle = 0;
	fEEnergy = 0.0;
	fNumberOfPhotons = 0;
	fNumberOfHitPMT = 0;
	fNumberOfHitPMTmu = 0;
	fVerboseFlag = 0;
	fNeutronTime =0;
}

BxDataCollection::~BxDataCollection(){
	Clear();
}

void BxDataCollection::Clear() {
//For standart Event at one point
	fEID = 0;
	fEPosition = G4ThreeVector(0,0,0);
	fEDirection = G4ThreeVector(0,0,0);
//      fShootDist = 0.;
	fEParticle = 0;
	fEEnergy = 0.0;
	fNumberOfPhotons = 0;
	fNumberOfHitPMT = 0;

	fNumberOfHitPMTmu = 0;

	fVerboseFlag = 0;
	
	fNeutronTime =0;
	fNeutronCaptureFlag =0;
	fPhData.clear();
}


PhotonData::PhotonData()  {
	fID = 0;
//	fStartWavelength = 0.0;
//	fEndWavelength = 0.0;
	fTrackLength = 0.0;
	fMeanStepLength = 0.0;
  	fNumberOfSteps = 0;
	fAbsZoneStatus = 0;
	fAbsMatStatus = 0;	
  	fPMTNumber = -1;		
  	fPMT_OD_Number = -1;
	fFlightTime = 0.0;
//	fPMTSpreadTime = 0.0;
//	fScintillatorTime = 0.0;
//	fPPOReemissionTime = 0.0;
//	fPCReemissionTime = 0.0;
  //Addition to control the reemission processes
  //      fNumPPOReemissions = 0;
    //    fNumPCReemissions = 0;
      //  fNumPCRelays = 0;
  //Addition to control the reflection processes	
	fNbOfPhReflections=0;

}
