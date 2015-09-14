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

#include "G4Material.hh"
#include "G4TransportationManager.hh"
//---------------------------------------------------------------------------//
/**Please use the stacking action command in macro
 * file /bx/stack/select to postpone the neutron capture
 * event. This is necessary for bx-elec simultion to
 * process the neutron capture as the new event.
 */
///Antineutrino + Proton reaction simulation


BxGeneratorSupernovaAntiNu::BxGeneratorSupernovaAntiNu(): BxVGenerator("BxGeneratorSupernovaAntiNu") {

  initFunc(5);

  isFirstTime = true ;
  fScintFlag  = false;
  fNeutrinoType = -1;
  gNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

  fPosition = G4ThreeVector(0.,0.,0.) ;
  fDirection = G4ThreeVector(1.,0.,0.);
  fParticleTable = G4ParticleTable::GetParticleTable();
  fSPSPos = new G4SPSPosDistribution ;
  G4SPSRandomGenerator *RndGen = new G4SPSRandomGenerator;
  fSPSPos->SetBiasRndm(RndGen); //для разыгрования по объему
  fSPSPos->SetPosDisType("Volume");
  fSPSPos->SetPosDisShape("Sphere");
  fSPSPos->SetCentreCoords(G4ThreeVector(0.,0.,0.));
  fSPSPos->SetRadius(4255.*mm);



  fTheMessenger = new BxGeneratorSupernovaAntiNuMessenger(this);
 
}
//---------------------------------------------------------------------------//


BxGeneratorSupernovaAntiNu::~BxGeneratorSupernovaAntiNu()
{
  delete fTheMessenger;
  delete fSPSPos;
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
        if(fScintFlag)
            fPosition = GetVParticlePositionInScint();
        G4cout << G4endl << G4endl << "Tuta ia!!!!" << G4endl << G4endl;


	if ((event->GetEventID())%2 == 0) {
            G4cout << G4endl << G4endl << "Suda ia!!!!" << G4endl << G4endl;

            G4double cosThetaPos = ShootAnglePositron();
            G4double enPos      = getPosEnergy(cosThetaPos)*MeV;
            G4cout << G4endl << G4endl << "enPos = " << enPos << G4endl << G4endl;
            G4double phiPos     = G4UniformRand()*360*degree;
            G4double phiNeutr   = 0.0*deg;
            if (phiPos >= 180.*deg)
                phiNeutr = phiPos - 180.*deg;
            else
                phiNeutr = phiPos + 180.*deg;


            if(GetNeutrinoType() == both || GetNeutrinoType() == positron) {  //разыгрываем позитрон

                fParticle = fParticleTable->FindParticle(-11); // -11 for the positron
                fDirection.setX(cos(phiPos)*cosThetaPos);
                fDirection.setY(sin(phiPos)*cosThetaPos);
                fDirection.setZ(-1.0*sqrt(1.-cosThetaPos*cosThetaPos));

                G4double pmom = sqrt(enPos*enPos - (fParticle->GetPDGMass())*(fParticle->GetPDGMass()));    //Здесь могут быть проблемы с размерностями
                G4double px = pmom*fDirection.x();
                G4double py = pmom*fDirection.y();
                G4double pz = pmom*fDirection.z();

                G4PrimaryVertex*   vertex   = new G4PrimaryVertex(fPosition,0);
                G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle,px,py,pz);
                vertex->SetPrimary( particle );
                event->AddPrimaryVertex( vertex );

                 BxOutputVertex::Get()->SetIsotopeCoinc(0);
                BxOutputVertex::Get()->SetEnergy(enPos - (fParticle->GetPDGMass()));
                BxOutputVertex::Get()->SetPosition(fPosition);
                BxOutputVertex::Get()->SetDirection(fDirection);
                BxOutputVertex::Get()->SetTime(0);
                BxOutputVertex::Get()->SetPDG(fParticle->GetPDGEncoding());
            }

            if(GetNeutrinoType() == both || GetNeutrinoType() == neutron) {  //разыгрываем нейтрон

                G4double kinNeutr = getKinNeutron(cosThetaPos)*MeV;
                G4double cosThetaNeutr = getAngleNeutron(enPos, kinNeutr);
                 G4cout << G4endl << G4endl << "cosThetaNeutr = " << cosThetaNeutr << G4endl << G4endl;
                 //cosThetaNeutr = 0.8;


                fParticle = fParticleTable->FindParticle(2112); // -11 for the positron
                fDirection.setX(cos(phiNeutr)*cosThetaNeutr);
                fDirection.setY(sin(phiNeutr)*cosThetaNeutr);
                fDirection.setZ(-1.0*sqrt(1.-cosThetaNeutr*cosThetaNeutr));

                G4double enNeutr = kinNeutr + fParticle->GetPDGMass();

                G4double pmom = sqrt(enNeutr*enNeutr - (fParticle->GetPDGMass())*(fParticle->GetPDGMass()));
                G4double px = pmom*fDirection.x();
                G4double py = pmom*fDirection.y();
                G4double pz = pmom*fDirection.z();

                G4PrimaryVertex*   vertex   = new G4PrimaryVertex(fPosition,0);
                G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle,px,py,pz);
                vertex->SetPrimary( particle );
                event->AddPrimaryVertex( vertex );

                BxOutputVertex::Get()->SetDId(0);
                BxOutputVertex::Get()->SetDPosition(fPosition);
                BxOutputVertex::Get()->SetDPDG(2112);
                BxOutputVertex::Get()->SetDTime(0.);
                BxOutputVertex::Get()->SetDEnergy(kinNeutr);
                BxOutputVertex::Get()->SetDaughters();
            }



	}

}

G4ThreeVector BxGeneratorSupernovaAntiNu::GetVParticlePositionInScint() {
    G4ThreeVector myPos(0.,0.,0.);
    myPos = fSPSPos->GenerateOne();
    while( !CheckMaterial( myPos, "Scintillator" ) ) myPos = fSPSPos->GenerateOne();

    return myPos;

}

G4bool BxGeneratorSupernovaAntiNu::CheckMaterial(G4ThreeVector pos, G4String MatName) {
    G4bool found = false ;

    G4ThreeVector null(0.,0.,0.);
    G4ThreeVector *ptr;
    ptr = &null;

    G4VPhysicalVolume *theVolume;
    theVolume = gNavigator->LocateGlobalPointAndSetup(pos,ptr,true);
    G4Material *amat = theVolume->GetLogicalVolume()->GetMaterial();
    G4String theMatName = amat->GetName();
    if(theMatName == MatName) {

        found = true ;
        BxLog(routine) << theMatName << endlog;
    }

    return found ;
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
    G4double deltaX, x, yy;
    for (int i=0; i < dS_dc.size(); i++) {
        if(pos_probability[i] >= val) {
            deltaX = val - pos_probability[i];
            yy =pos_energyBin[i]-pos_energyBin[i-1];
            x= pos_probability[i]-pos_probability[i-1];
            return deltaX*yy/x + pos_energyBin[i];
        }
    }
    return -2.; //подумай

}


