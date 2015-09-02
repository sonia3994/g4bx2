// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#ifndef _BXGENERATORRDMDecayGun_HH
#define _BXGENERATORRDMDecayGun_HH

//---------------------------------------------------------------------------//

#include "G4GeneralParticleSource.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "BxVGenerator.hh"
#include "BxReadParameters.hh"
#include "G4Event.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"
#include "BxGeneratorRDMNucleus.hh"
#include "G4SPSRandomGenerator.hh"
//---------------------------------------------------------------------------//

class BxGeneratorRDMDecayGunMessenger ;

///It shoots the desired nucleus
/**When you select the RDM generator, you MUST switch on radioactive decay process with: /bx/physics/decay true
*If you want to generate a decay chain segment, starting from the source isotope, activate the /bx/stack/select RDMChain.
*/
class BxGeneratorRDMDecayGun : public BxVGenerator {
	public:

		///default constructor
		BxGeneratorRDMDecayGun();

		///destructor
		virtual ~BxGeneratorRDMDecayGun();

		/// Sets the isotope.
		void  SetNucleus(BxGeneratorRDMNucleus theIon1);
		/// Returns the specified isotope.
		BxGeneratorRDMNucleus GetNucleus() {return theIon;}

		///public interface
		virtual void BxGeneratePrimaries(G4Event *event);


		///Sets the particle position
		void SetParticlePosition(G4ThreeVector pos)               { fParticleGun->SetParticlePosition(pos); }
		///Sets the particle: i.e. it is the ion 
		void SetParticleDefinition(G4ParticleDefinition* part)    { fParticleGun->SetParticleDefinition(part); }
		///Sets the particle momentum distribution
		void SetParticleMomentumDirection(G4ParticleMomentum mom) { fParticleGun->SetParticleMomentumDirection(mom); }
		///Sets the particle energy 
		void SetParticleEnergy(G4double ene)                      { fParticleGun->SetParticleEnergy(ene); }
		///Sets if the background is generated inside a volume
		void SetIfVolDist(G4bool a)         { fVolumeFlag = a; }

		///Allows user to choose Point, Plane, Surface or Volume source position distributions.
		void SetPosDisType(G4String string) { fSPSPos->SetPosDisType(string); }

		///Allows the user to choose the particular shape they wish for the position distribution. Choices are Square, Circle, Ellipse, Rectangle, Sphere, Ellipsoid, Cylinder, Parallelepiped.
		void SetPosDisShape(G4String string) { fSPSPos->SetPosDisShape(string); }

		///Sets the co-ordinates of the centre of the position distribution.
		void SetCentreCoords(G4ThreeVector pos) { fSPSPos->SetCentreCoords(pos); }

		///Used to specify the co-ordinate system for the position distribution along with SetPosRot2. SetPosRot1 sets the vector x' and need not be a unit vector.
		void SetPosRot1(G4ThreeVector rot) { fSPSPos->SetPosRot1(rot); }

		///Used in connection with SetPosRot1. This sets a vector in the plane x'y'. By a series of cross products x', y', z' are generated. Again need not be a unit vector.
		void SetPosRot2(G4ThreeVector rot){ fSPSPos->SetPosRot2(rot); }

		///Sets the radius where appropriate for source distribution shapes.
		void SetRadius(G4double radius){ fSPSPos->SetRadius(radius); }

		///Sets the inner radius where appropriate for source distribution shapes.
		void SetRadius0(G4double radius){ fSPSPos->SetRadius0(radius); }

		///Used to confine the start positions to a particular volume.
		void ConfineSourceToVolume(G4String vol){ fSPSPos->ConfineSourceToVolume(vol); }



		//private  members
	private:
		///The messenger
		BxGeneratorRDMDecayGunMessenger*   fTheMessenger;
		G4ParticleGun*               fParticleGun;
		///Position distribution
		G4SPSPosDistribution*        fSPSPos;
		///Angular distribution
		G4SPSAngDistribution*        fSPSAng;
		G4bool                       fVolumeFlag ;
		///This class contains the characteristics of the nucleus
		BxGeneratorRDMNucleus        theIon;
		///The mass number
		G4int                        fA ; 
		///The atomic nu,ber
		G4int                        fZ;
		///The nucleus energy
		G4double                     fE ;
		G4bool                       IsFirst;
		G4SPSRandomGenerator*        fRndGen;
};
#endif
