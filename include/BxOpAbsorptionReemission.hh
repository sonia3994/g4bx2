/** //////////////////////////////////////////////////////////////////////
* Bx Optical Photon Absorption Class Definition
* ////////////////////////////////////////////////////////////////////////
*
* File:        BxOpAbsorptionReemission.hh
* Description: Discrete Process -- Absorption and reemission of Optical Photons
* Version:     1.0
* Created:     1998-02-06
* Author:      Stefano Magni
* mail:        Stefano.Magni@mi.infn.it
*
* Revised by A. Caminata and S. Marcocci, Sept. 2014
* //////////////////////////////////////////////////////////////////////
*/
#ifndef BxOPAbsorptionREEMISSION_H
#define BxOPAbsorptionREEMISSION_H

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4OpticalPhoton.hh"


#include "G4PhysicsTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "BxLightSource.hh"

/////////////////////
// Class Definition
/////////////////////

///It simulates the absorption and reemission of light in the scintillator
class BxOpAbsorptionReemission : public G4VDiscreteProcess 
{

  public:

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

        BxOpAbsorptionReemission(const G4String& processName =
						"AbsorptionReemission");

        // G4OpAbsorption(const G4OpAbsorption &right);

	~BxOpAbsorptionReemission();

	////////////
	// Methods
        ////////////

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);

	G4double GetMeanFreePath(const G4Track& aTrack,
				G4double ,
				G4ForceCondition* );

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
					const G4Step&  aStep);


 
  private:

  	//const G4DynamicParticle* 	aParticle;
	G4int 				fAbsorber;
        
	G4int                           fMaterialIndexPC;
	G4int                           fMaterialIndexPPO;
	G4int                           fMaterialIndexDMP;

	size_t                          fScintillatorIndex;
	size_t                          fBufferIndex;
	
	G4MaterialPropertyVector*	AbsorptionLengthVector[10];

	G4MaterialPropertyVector* 	ReemissionProbabilityVector1;
	G4MaterialPropertyVector* 	ReemissionProbabilityVector2;

	G4bool 				G4BooleanRand(const G4double prob)
									const;
	void 				DoReemission(BxLightSource::SourceType type,
					const G4Track& aTrack);
	void 				DoAbsorption ();
	void 				DoScattering (const G4Track& aTrack);
};


////////////////////
// Inline methods
////////////////////

inline
G4bool BxOpAbsorptionReemission::IsApplicable(const G4ParticleDefinition& aParticleType)
{
   return ( &aParticleType == G4OpticalPhoton::OpticalPhoton() );
}

inline
G4bool BxOpAbsorptionReemission::G4BooleanRand(const G4double prob) const
{
  /* Returns a random boolean variable with the specified probability */

  return (G4UniformRand() < prob);
}

#endif 
