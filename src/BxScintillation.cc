//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "BxOutputVertex.hh"
#include "G4Material.hh"
#include "BxScintillation.hh" 
#include "G4Electron.hh" 
#include "G4Alpha.hh" 
#include "G4Proton.hh"  
#include "G4SystemOfUnits.hh"
#include "BxLightSource.hh"
#include "BxIO.hh"
using namespace std ;

BxScintillation::BxScintillation(const G4String& processName, const G4int& particlePDG)
                  : G4VRestDiscreteProcess(processName) {
  
  /**IMPORTANT 
  * Quenching formalization requires that secondaries are tracked 
  * at the beginning!! 
  * Leave this flag in the false status!!!!
  * For further explanation aske to DF
  * davide.franco@mi.infn.it
  *fTrackSecondariesFirst = false;
  *
  * false (default) for standard physics with quenching
  * true  for muon simulations
  */
  fTrackSecondariesFirst = BxOutputVertex::Get()->GetTrackSecondariesFirst() ;
   
  BxLog(trace) << "Track Secondaries First: " << fTrackSecondariesFirst << endlog ;
  //-------------------------------
   
  fParticlePDG = particlePDG;
 
  fMeanExcitationEnergy  = 6.04e-5 ;

  emass = G4Electron::Definition()->GetPDGMass ()/MeV;
  amass = G4Alpha::Definition()->GetPDGMass ()/MeV;
  pmass = G4Proton::Definition()->GetPDGMass ()/MeV;
  // Electron quenching: coefficient in the Bethe equation
  // fBetheCoeff = K/2 * ro * Z/A
  //  ro  = 0.877
  //  K   = 0.307075  (PDG)
  //  Z/A = 0.55
  //  ref: PDG
  fBetheCoeff = G4double(0.877)*G4double(0.307075)*G4double(0.55/2.);

  fResolutionScale = 1.0;
  fYieldFactor     = 1.0;
  fExcitationRatio = 1.0;
  fEMinusEPlusEndPoint = 5.*MeV ;
  fAlphaEndPoint = 10.*MeV ;
  fProtonEndPoint = 20.*MeV ;  
  fMaxParentId    =  BxReadParameters::Get()->GetMaxParentIdForQuenching();
  
  fNumberOfStep = 20000 ;
  NumAlphaSteps = 0;
  
  fParticleCode=0;
  
  BxOutputVertex::Get()->SetNPhotons(0);
    
  fScintillationYieldBulk = (BxReadParameters::Get()->GetLightYield())*(BxReadParameters::Get()->GetLightYieldScale());

  // Quenching factor in buffer
  fScintillationYieldBuffer = fScintillationYield*0.10*G4double(0.40/0.82);
  fBufferQuenching          = 0.10*G4double(0.40/0.82);
  // Quenching variables
  fBirksConstantBeta  = BxReadParameters::Get()->GetBirksBeta();
  fBirksConstantBeta2  = BxReadParameters::Get()->GetBirksSecondOrderBeta();
  fBirksConstantAlpha  = BxReadParameters::Get()->GetBirksAlpha();
  fBirksConstantAlpha2 = BxReadParameters::Get()->GetBirksSecondOrderAlpha();
  fBirksConstantProton  = BxReadParameters::Get()->GetBirksProton();
  fBirksConstantProton2 = BxReadParameters::Get()->GetBirksSecondOrderProton();  
 
  fComputeBetaBirks = false ;
  fComputeAlphaBirks   = false ;
  fComputeProtonBirks = false ; 
  
  fScintillatorIndex = BxOutputVertex::Get()->GetScintillatorIndex();
  fDMPBufferIndex    = BxOutputVertex::Get()->GetDMPBufferIndex();

}

BxScintillation::~BxScintillation()  {}

G4VParticleChange* BxScintillation::AtRestDoIt(const G4Track& aTrack, const G4Step& aStep) {
  return BxScintillation::PostStepDoIt(aTrack, aStep);

}

G4VParticleChange* BxScintillation::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) {
  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle* sParticle = aTrack.GetDynamicParticle();
  const G4Material* aMaterial = aTrack.GetMaterial();

  if((aMaterial->GetName() != "Scintillator") &&(aMaterial->GetName() != "DMPbuffer"))	{
    aParticleChange.SetNumberOfSecondaries(0);
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);

  } 

  if(!fTrackSecondariesFirst) {
  
    // Compute beta quenching array for positrons and electrons
    if(!fComputeBetaBirks) BetaQuenchingFactorVector(fEMinusEPlusEndPoint/MeV, fBirksConstantBeta, fBirksConstantBeta2);

    // Compute alpha quenching array up
    if(!fComputeAlphaBirks) AlphaQuenchingFactorVector(fBirksConstantAlpha, fBirksConstantAlpha2);

    // Compute proton quenching array up
    if(!fComputeProtonBirks) ProtonQuenchingFactorVector(fProtonEndPoint/MeV, fBirksConstantProton, fBirksConstantProton2); 
  
  }
  
  G4double EnergyDeposit = aStep.GetTotalEnergyDeposit();
  if(EnergyDeposit == 0.) {
    aParticleChange.SetNumberOfSecondaries(0);
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
 

  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

  G4ThreeVector x0 = pPreStepPoint->GetPosition();
  G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
  G4double      t0 = pPreStepPoint->GetGlobalTime();

   //
  //Quenching for alpha and and proton and beta particles
  //
  if(!fTrackSecondariesFirst && aTrack.GetParentID () < fMaxParentId ) {
  
    //positron
    if(aTrack.GetDefinition()->GetPDGEncoding() == -11) 
			fQuenchingFactor = Interpolator(aTrack.GetVertexKineticEnergy(), fEMinusEPlusEnergy, fBetaPlusFactor) ;
		 
    //electron
    else if(aTrack.GetDefinition()->GetPDGEncoding() == 11) 
			fQuenchingFactor = Interpolator(aTrack.GetVertexKineticEnergy(), fEMinusEPlusEnergy, fBetaFactor) ;
    
    //Alpha
    else if(aTrack.GetDefinition()->GetPDGEncoding() == 1000020040) 
      fQuenchingFactor = Interpolator(aTrack.GetVertexKineticEnergy(), fAlphaEnergy, fAlphaFactor) ; 

    //Proton
    else if(aTrack.GetDefinition()->GetPDGEncoding() == 2212)    
    fQuenchingFactor = Interpolator(aTrack.GetVertexKineticEnergy(), fProtonEnergy, fProtonFactor) ;   
  }

 //If it is a positron, electron, alpha or proton 
  if(((!fTrackSecondariesFirst) && (abs(aTrack.GetDefinition()->GetPDGEncoding()) == 11) )|| (aTrack.GetDefinition()->GetPDGEncoding() == 1000020040 ) || (aTrack.GetDefinition()->GetPDGEncoding() == 2212))     
     EnergyDeposit *= fQuenchingFactor ;
  //  cout << fQuenchingFactor << " " <<  aTrack.GetTrackID () << " " << aTrack.GetParentID () << " " << endl ;
  //  cout << aTrack.GetStep ()->GetPostStepPoint ()->GetProcessDefinedStep ()->GetProcessName () << " " 
  //     <<  aTrack.GetStep ()->GetPostStepPoint ()->GetProcessDefinedStep ()-> GetProcessType () << endl ;

  //  Quenching for muon pdg= 13 - can be put here in future
  //  Now the muon crossing the Scintillator do not produce the scintillation
/*cout << aTrack.GetTrackID () << " " 
       << aTrack.GetParentID () << " " 
       << aTrack.GetDefinition()->GetParticleName() << " "
       <<pPreStepPoint->GetKineticEnergy()<< " " 
       << pPreStepPoint->GetLocalTime()/ns << " " 
       << fQuenchingFactor 
       << endl ;
  */

      
  fScintillationYield  =  fScintillationYieldBulk;

  //Here is the part for Light Yield reduction in DMPbuffer
  if(aMaterial->GetName() == "DMPbuffer") EnergyDeposit *= fBufferQuenching;

  G4double MeanNumPhotons = fScintillationYield * EnergyDeposit;
  
 // G4cout << "aTrack.GetDefinition()->GetPDGEncoding() =" << aTrack.GetDefinition()->GetPDGEncoding() << G4endl;
 // G4cout << MeanNumPhotons << " MeanNumPhotons" << G4endl ;
  
  BxOutputVertex::Get()->SetVisEnergy(EnergyDeposit + BxOutputVertex::Get()->GetVisEnergy());

  //bias of the simulation. SM: Fluctuations a la Davide F. are better (fluctuate after you divide by maxQE) 
  MeanNumPhotons=MeanNumPhotons*BxReadParameters::Get()->GetPMTQEMaximum();
  G4int NumPhotons = 0;

  if (MeanNumPhotons > 30.) {
    G4double sigma = fResolutionScale * sqrt(MeanNumPhotons);
    NumPhotons = G4int(G4RandGauss::shoot(MeanNumPhotons,sigma)+0.5);
  } else {
    NumPhotons = G4int(G4Poisson(MeanNumPhotons));
  }
  if (NumPhotons <= 0) {
    aParticleChange.SetNumberOfSecondaries(0);
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  
  BxOutputVertex::Get()->SetNPhotons(G4int(NumPhotons) + BxOutputVertex::Get()->GetNPhotons());
  G4double thePhotonEnergy;

  aParticleChange.SetNumberOfSecondaries(NumPhotons);
 	BxLog(development) << "Number of generated photons: " << NumPhotons << endlog;
  
  if (fTrackSecondariesFirst) {
     if (aTrack.GetTrackStatus() == fAlive )
	     aParticleChange.ProposeTrackStatus(fSuspend);
  }

  for (G4int i = 0; i < NumPhotons; i++) {

    if (aMaterial->GetName() == "Scintillator") thePhotonEnergy = BxLightSource::Get()->DefineEmissionEnergy(BxLightSource::PPOEmission);
    else thePhotonEnergy = BxLightSource::Get()->DefineEmissionEnergy(BxLightSource::PCEmission); 

    // Generate random photon direction

    G4double cost = 1. - 2.*G4UniformRand();
    G4double sint = sqrt((1.-cost)*(1.+cost));

    G4double phi = 2*M_PI*G4UniformRand();
    G4double sinp = sin(phi);
    G4double cosp = cos(phi);

    G4double px = sint*cosp;
    G4double py = sint*sinp;
    G4double pz = cost;

    // Create photon momentum direction vector 

    G4ParticleMomentum thePhotonMomentum(px, py, pz);

    // Determine polarization of new photon 

    G4double sx = cost*cosp;
    G4double sy = cost*sinp; 
    G4double sz = -sint;

    G4ThreeVector thePhotonPolarization(sx, sy, sz);

    G4ThreeVector perp = thePhotonMomentum.cross(thePhotonPolarization);

    phi = 2*M_PI*G4UniformRand();
    sinp = sin(phi);
    cosp = cos(phi);

    thePhotonPolarization = cosp * thePhotonPolarization + sinp * perp;

    thePhotonPolarization = thePhotonPolarization.unit();

    // Generate a new photon:

    G4DynamicParticle* aScintillationPhoton =
      new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), 
  					     thePhotonMomentum);
    aScintillationPhoton->SetPolarization
			 (thePhotonPolarization.x(),
			  thePhotonPolarization.y(),
			  thePhotonPolarization.z());

    aScintillationPhoton->SetKineticEnergy(thePhotonEnergy);

    // Generate new G4Track object:

    G4double rand;

    if (sParticle->GetDefinition()->GetPDGCharge() != 0)
       rand = G4UniformRand();
    else 
       rand = 1.0;


    G4double delta = aStep.GetStepLength();
    G4double deltaTime = delta /
	   ((pPreStepPoint->GetVelocity()+
             pPostStepPoint->GetVelocity())/2.);


  //HERE THE EMISSION TIME is defined

//if alpha or proton
    if (fParticlePDG == 1000020040 || fParticlePDG == 2212) fParticleCode=1;
    else fParticleCode=0;

  //Here is the part for different scintillationtime  in Scintillator or DMPbuffer

    G4double theScintillationTime =0;

    if((aMaterial->GetName() == "Scintillator"))	{
      theScintillationTime =  BxLightSource::Get()->DefineEmissionTime(BxLightSource::ScintillatorEmission,fParticleCode);
    } else {
      theScintillationTime =  BxLightSource::Get()->DefineEmissionTime(BxLightSource::PCinDMPEmission,fParticleCode);
    }

    G4double aSecondaryTime = t0 + deltaTime + theScintillationTime;

    G4ThreeVector aSecondaryPosition = x0 + rand * aStep.GetDeltaPosition();

    G4Track* aSecondaryTrack = 
	    new G4Track(aScintillationPhoton,aSecondaryTime,aSecondaryPosition);

    aSecondaryTrack->SetTouchableHandle((G4VTouchable*)0);

    aSecondaryTrack->SetParentID(aTrack.GetTrackID());

    aParticleChange.AddSecondary(aSecondaryTrack);
  }

  if (verboseLevel>0) {
    G4cout << "\n Exiting from scintillation::DoIt -- NumberOfSecondaries = " 
       << aParticleChange.GetNumberOfSecondaries() << G4endl;
  }

  return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);			
}

G4double BxScintillation::GetMeanFreePath(const G4Track& , G4double , G4ForceCondition* condition) {
  *condition = StronglyForced;
  return DBL_MAX;
}

G4double BxScintillation::GetMeanLifeTime(const G4Track& , G4ForceCondition* condition) {
   *condition = StronglyForced;
   return DBL_MAX;
}


// Analitical stopping power for the quenching model: 
// electrons and alphas: Sternheimer, Berger, Seltzer - Density effect for Ionization Loss
// positrons: F. Rohrlich and B. C. Carlson, 
//            Positron-electron differences in energy loss and multiple scattering
//            Phys. Rew. 93 (1954) 38-44
 

void BxScintillation::BetaQuenchingFactorVector(G4double E, G4double kB, G4double kB2) {
  G4double beta2, Factor, bin, DEDB,gamma,FactorPlus,DEDBPlus ;


  G4double step = E /fNumberOfStep ;
  G4double qeneMinus = 0. ;
  G4double qenePlus  = 0. ;
  for (G4int i = 0; i<fNumberOfStep; i++) {

    bin    = G4double(i)*step + step/3.;			 
    beta2  = 1. - pow(emass/(bin+emass),2.);
    gamma  = 1./sqrt(1 - beta2) ;
    
      // BetaMinus
    Factor =   log(emass*bin*beta2/(1-beta2)) 
							- log(2.)*(2*sqrt(1-beta2)-1+beta2)
							+ 1 - beta2 + 1./8.*(1-sqrt(1-beta2))*(1-sqrt(1-beta2))
							- 2.*log(fMeanExcitationEnergy) - log(2.);

    DEDB = fBetheCoeff/beta2*Factor;
    
    
    
    if(DEDB <= 0) {
	    fBetaFactor.push_back(0.01) ;

    } else {
	    qeneMinus += step / (1. + kB*DEDB+kB2*pow(DEDB,2.));
	    fBetaFactor.push_back(qeneMinus/bin) ;
	    fEMinusEPlusEnergy.push_back(bin) ;
    }


		
  
      // BetaPlus
    FactorPlus =   log(pow(bin,2.)*(gamma+1.)/2.)
									+ 2.*log(2.)  - beta2/12.*  
									(23.+14./(gamma+1)+10./pow(gamma+1,2.)+4./pow(gamma+1,3.) )
									- 2.*log(fMeanExcitationEnergy);

    DEDBPlus = fBetheCoeff/beta2*FactorPlus;
    if(DEDBPlus <= 0) {
      fBetaPlusFactor.push_back(0.01) ;
     
    } else {
      qenePlus += step / (1. + kB*DEDBPlus+kB2*pow(DEDBPlus,2.));
			fBetaPlusFactor.push_back(qenePlus/bin);
    }
  }
	
  fComputeBetaBirks = true ;
}  

//void BxScintillation::AlphaQuenchingFactorVector(G4double E, G4double kB, G4double kB2) {
void BxScintillation::AlphaQuenchingFactorVector(G4double kB, G4double kB2) {
	G4double bin = 0., DEDB = 0., step = 0., qeneAlpha = 0.;
	
	while(BxIO::Get()->GetStreamPCAlphadEdx() >> bin >> DEDB){
		if(step == 0){
			step = bin;
		}

		if(DEDB <= 0){
			fAlphaFactor.push_back(0.0001);
		} else {
			qeneAlpha += step / (1. + kB*DEDB + kB2*pow(DEDB,2));
			fAlphaFactor.push_back(qeneAlpha/bin);
			fAlphaEnergy.push_back(bin);
		}
		
	}
  fComputeAlphaBirks = true ;
}  

void BxScintillation::ProtonQuenchingFactorVector(G4double E, G4double kB, G4double kB2) {
  G4double beta2, Factor, bin, DEDB;

  G4double step = E /fNumberOfStep ;
  G4double qeneProton = 0. ;

  for (G4int i = 0; i<fNumberOfStep; i++) {

    bin    = G4double(i)*step + step/3.;			 
    beta2  = 1. - pow(pmass/(bin+pmass),2.);
 
    Factor =   2.*log(2.*emass*beta2/(1-beta2)) 
             - 2.*log(fMeanExcitationEnergy);
    
    DEDB = fBetheCoeff/beta2*Factor;
    if(DEDB <= 0) {
      fProtonFactor.push_back(0.0001) ;
    } else {
      qeneProton += step / (1. + kB*1.*DEDB+ kB2*pow(DEDB,2));
      fProtonFactor.push_back(qeneProton/bin);
      fProtonEnergy.push_back(bin) ;

    }

  }
  fComputeProtonBirks = true ;
}  


G4double BxScintillation::Interpolator(G4double Ene, vector<G4double>& Nrg, vector<G4double>& QuenchingFactor){



	G4double QuenchingFactorInterp = 0.;
	G4double dQ = 0.;
	G4double dNrg = 0.;


	if(Ene < Nrg.front()){
		dQ   = QuenchingFactor[0] - 1;
		dNrg = Nrg[0];
		QuenchingFactorInterp = 1 + Ene*dQ/dNrg;
	}

	else if(Ene >= Nrg.back())
		QuenchingFactorInterp = QuenchingFactor.back();




	else{
		G4int j = 0;
		for(G4int i = 0; i < G4int(Nrg.size()); i++){
			if(Ene < Nrg[i]){
				j = i;
				break;
			}	
		}
	
		dQ   = QuenchingFactor[j] - QuenchingFactor[j-1];

		dNrg = Nrg[j] - Nrg[j-1];

		QuenchingFactorInterp = QuenchingFactor[j-1] + (Ene - Nrg[j-1])*dQ/dNrg;

	}

	return QuenchingFactorInterp;
	
}
