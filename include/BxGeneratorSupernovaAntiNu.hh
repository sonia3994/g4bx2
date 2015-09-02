#ifndef _BXGENERATORSUPERNOVAANTINU_HH
#define _BXGENERATORSUPERNOVAANTINU_HH

//---------------------------------------------------------------------------//


#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "BxVGenerator.hh"
#include "BxGeneratorSupernovaAntiNuMessenger.hh"
#include "BxReadParameters.hh"
#include "G4Event.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"

//---------------------------------------------------------------------------//
/**Please use the stacking action command in macro
 * file /bx/stack/select to postpone the neutron capture
 * event. This is necessary for bx-elec simultion to
 * process the neutron capture as the new event.
 */
///Antineutrino + Proton reaction simulation
class BxGeneratorSupernovaAntiNu : public BxVGenerator {

 public:
	///default constructor
  	BxGeneratorSupernovaAntiNu();
	///destructor
  	virtual ~BxGeneratorSupernovaAntiNu();

  	///public interface
  	virtual void BxGeneratePrimaries(G4Event *event);

 private:
	G4double ShootAnglePositron();
	G4double GetKinEnergyPositron(G4double angle);
	G4double GetKinEnergyNeutron(G4double angle, G4double energyPositron);
	G4double GetAngleNeutron(G4double angle);

	BxGeneratorSupernovaAntiNuMessenger*   fTheMessenger;
	
	G4ParticleTable*             fParticleTable;
	G4ParticleDefinition*        fParticle;
	G4ThreeVector                fPosition;
	G4ThreeVector                fDirection;

	G4bool                       isFirstTime;

};
#endif

	
