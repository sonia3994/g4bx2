#include "BxEventAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4Neutron.hh"
#include "BxLogger.hh"
#include "BxDataCollection.hh"
#include "BxOutputVertex.hh"
#include "BxGeneratorSupernovaAntiNu.hh"
#include "G4SystemOfUnits.hh"
#include "BxGeneratorSupernovaAntiNuMessenger.hh"
//---------------------------------------------------------------------------//
/**Please use the stacking action command in macro
 * file /bx/stack/select to postpone the neutron capture
 * event. This is necessary for bx-elec simultion to
 * process the neutron capture as the new event.
 */
///Antineutrino + Proton reaction simulation


BxGeneratorSupernovaAntiNu::BxGeneratorSupernovaAntiNu(): BxVGenerator("BxGeneratorSupernovaAntiNu") {


  isFirstTime = true ; 

  fPosition = G4ThreeVector(0.,0.,0.) ;
  fDirection = G4ThreeVector(1.,0.,0.);
  fParticleTable = G4ParticleTable::GetParticleTable();
  fSPSPos = new G4SPSPosDistribution ;
  G4SPSRandomGenerator *RndGen = new G4SPSRandomGenerator;
  fSPSPos->SetBiasRndm(RndGen);

  fTheMessenger = new BxGeneratorSupernovaAntiNuMessenger(this);
 
}
//---------------------------------------------------------------------------//


BxGeneratorSupernovaAntiNu::~BxGeneratorSupernovaAntiNu()
{
  delete fTheMessenger;
  //delete fSPSPos;
  //delete fSPSAng;
  
}

//---------------------------------------------------------------------------//

void BxGeneratorSupernovaAntiNu::BxGeneratePrimaries(G4Event *event) {


	if(isFirstTime) {
		if(fNeutrinoType < 0) {
      	    BxLog(error) << "Set the antineutrino source type"<<endlog;
      	    BxLog(fatal) << endlog;
    	}

    	  BxLog(routine) << "AntiNeutrino source " << fNeutrinoType << endlog ;

    	  isFirstTime = false ;
	}

	if ((event->GetEventID())%2 == 0) {

		fParticle = fParticleTable->FindParticle(-11); // -11 for the positron
                G4double anglePos = ShootAnglePositron();
	}

}

void BxGeneratorSupernovaAntiNu::initFunc(G4double eN) {
    //вспомогрательные переменные
    G4double summ=0.;
    G4double norma=0.;
    G4double init=-1.;
    G4double step=0.01;

    eNu = eN;
    e0  = eNu - delt;
    ve0 = pow((e0*e0-me*me),1./2.)/e0;
    for (int i=0; i < 201; i++) {
        pos_energyBin.push_back(init + step*i); //pos_energyBin <=> cosTeta
        G4double posE = getPosEnergy(pos_energyBin[i]);
        G4double gamma = getBigGamma(pos_energyBin[i]);
        dS_dc.push_back(dSigma_dcos(gamma,posE,pos_energyBin[i]));
        if (i !=0) //dS_dc[0] = 0.
            norma +=dS_dc[i];
    }
    dS_dc[0] = 0.; //для правильнной интерполяции
    for (int i=0; i < dS_dc.size(); i++) {
        summ += dS_dc[i]/norma;
        pos_probability.push_back(summ);
    }
}

G4double BxGeneratorSupernovaAntiNu::ShootAnglePositron() {
    G4double val = G4UniformRand();   //узнать, можно ли так делать
    G4double deltaX,x,y;
    for (int i=0; i < dS_dc.size(); i++) {
        if(pos_probability[i] >= val) {
            deltaX = val - pos_probability[i];
            y =pos_energyBin[i]-pos_energyBin[i-1];
            x= pos_probability[i]-pos.probability[i-1];
            return deltaX*y/x + pos_energyBin[i];
        }
    }
    return -2.; //подумай

}


