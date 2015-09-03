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
        G4double const pi  = 3.14159265358979323846;
        G4double const f  = 1.0;
        G4double const g  = 1.26;
        G4double const f2 = 3.706;
        G4double const delt = 939.565378 - 938.272046;
        G4double const M    = (939.565378 + 938.272046)/2.;
        G4double const me   = 0.510999;
        G4double const mn   = 939.565378;
        G4double const mp   = 938.272046;
        G4double const cosUCab   = 0.974;
        G4double const deltInRad = 0.024;
        G4double const Gfermi    = 1.166 /100000000000. ;
        G4double const y         = (delt*delt-me*me)/2.;
        G4double const sigma0    = ((Gfermi*Gfermi)*cosUCab*cosUCab)*(1+deltInRad)/pi;

        void     initFunc(G4double eNu);
	G4double ShootAnglePositron();
	G4double GetKinEnergyPositron(G4double angle);
	G4double GetKinEnergyNeutron(G4double angle, G4double energyPositron);
	G4double GetAngleNeutron(G4double angle);

        G4ThreeVector  GetVParticlePosition(); //замутить ивент рандомно в определенном объеме


	BxGeneratorSupernovaAntiNuMessenger*   fTheMessenger;
	
	G4ParticleTable*             fParticleTable;
	G4ParticleDefinition*        fParticle;
	G4ThreeVector                fPosition;
	G4ThreeVector                fDirection;

	G4bool                       isFirstTime;
        G4SPSPosDistribution*        fSPSPos;  //для GetVParticlePosition()

};
#endif

	
