//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef BxGeneratorRDMDecayGunMessenger_h
#define BxGeneratorRDMDecayGunMessenger_h 1
#include "BxIO.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4ParticleGun.hh"
#include "BxGeneratorRDMUIcmdWithNucleusAndUnit.hh"

class BxPrimaryGeneratorAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class BxGeneratorRDMDecayGun;

using namespace std;
///Messenger of BxGeneratorRDMDecayGun
class BxGeneratorRDMDecayGunMessenger: public G4UImessenger {

	public:
		///Ctor
		BxGeneratorRDMDecayGunMessenger(BxGeneratorRDMDecayGun*  );
		///Dtor
		~BxGeneratorRDMDecayGunMessenger();
		///Sets the new value
		void SetNewValue(G4UIcommand*, G4String);


	private:
		///The related generator
		BxGeneratorRDMDecayGun*              generator;   
		/// Control of BxrdmDecayGun event generator
		G4UIdirectory*                       fDirectory;
		///Set the gun position
		G4UIcmdWith3VectorAndUnit*           fPositionCmd; 
		/// Set the gun direction  
		G4UIcmdWith3Vector*                  fDirectionCmd;
		///Set the gun energy
		G4UIcmdWithADoubleAndUnit*           fEnergyCmd;
		///Set the gun type (i.e. wich nucleus)
		G4UIcmdWithAString*                  fParticleCmd;
		/// Bulk radius   
		G4UIcmdWithADoubleAndUnit*           fSphereBulkCmd;
		/// Surface radius
		G4UIcmdWithADoubleAndUnit*           fSphereSurfCmd;
		///Set the sphere position
		G4UIcmdWith3VectorAndUnit*           fSphereCentreCmd; 
		/// Bulk minimum radius 
		G4UIcmdWithADoubleAndUnit*           fSphereBulkRadMinCmd;
		//confine the particle generation on a vessel
		G4UIcmdWithABool*                    fVesselCmd;
		//confine the particle generation in the bulk 
		G4UIcmdWithABool*                    fBulkCmd;
		//confine the particle generation on the buffer
		G4UIcmdWithABool*                    fBufferCmd;
		///confine the particle generation in a given volume
		G4UIcmdWithAString*                  fConfineCmd ;
		///define the primary ion (a,z,e)   
		BxGeneratorRDMUIcmdWithNucleusAndUnit*      fIonCmd;
};

#endif
