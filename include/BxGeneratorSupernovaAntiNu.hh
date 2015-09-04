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

#include <math>
#include <vector>




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
        //вспомогательные константы и переменные
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

        G4double eNu;
        G4double e0;
        G4double ve0;
        //std::vector<G4double> cosTeta;
        std::vector<G4double> dS_dc;
        std::vector<G4double> pos_probability;
        std::vector<G4double> pos_energyBin; //pos_energyBin <=> cosTeta

        // вспомогательные функции
        void     initFunc(G4double);
        G4double getPosEnergy(G4double cosTeta) {
            return e0*(1-(eNu/M)*(1-ve0*cosTeta)) - y*y/M;
        }
        G4double dSigma_dcos(G, posE, cosTeta) {
            return (sigma0/2.)*((f*f+3*g*g)+(f*f-g*g)*(pow((posE*posE - me*me), 1./2.) /
                             (posE))*cosTeta)*posE*pow((posE*posE - me*me), 1./2.) -
                    (sigma0/2.)*(G/M)*e0*pow((e0*e0-me*me),1./2.);
        }

        G4double getBigGamma(G4double cosTeta) {
            return 2*(f+f2)*g*((2*e0 + delt)*(1-ve0*cosTeta) - (me**2)/e0) +
                    (f**2 + g**2)*(delt*(1+ve0*cosTeta) + (me**2)/e0) +
                    (f**2 + 3*g*g)*((e0 + delt)*(1 - (1/ve0)*cosTeta) - delt) +
                    (f**2 - g**2)*((e0 + delt)*(1 - (1/ve0)*cosTeta) - delt)*ve0*cosTeta;
        }
        G4double getKinNeutron(G4double cosTeta) {
            return (eNu*e0/M)*(1-ve0*cosTeta) + y**2/M;
        }
        G4double getAngleNeutron(G4double cosTeta, enPos, kinNeutr) {
            return (mp*mp + me*me - 2*mp*enPos - mn**2 +2*eNu*(kinNeutr+mn))/
                    (2*pow(((kinNeutr+mn)*(kinNeutr+mn)-mn*mn),(1./2.))*eNu);
        }

        G4double ShootAnglePositron();


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

	
