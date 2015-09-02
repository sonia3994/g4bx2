// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//
#include <TF1.h>
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4Electron.hh"

#include "BxLogger.hh"
#include "BxDataCollection.hh"
#include "BxOutputVertex.hh"
#include "BxGeneratorSolarNeutrino.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>
//---------------------------------------------------------------------------//


BxGeneratorSolarNeutrino::BxGeneratorSolarNeutrino(): BxVGenerator("BxGeneratorSolarNeutrino") {
  //pi = 3.14159265;
  emass  = G4Electron::Definition()->GetPDGMass ()/MeV;
  DM2    = 7.59 ;
  TG2T   = 0.42 ;
  isFirstTime = true ;
  fNeutrinoType = -1;
  fOldSigma=false;
  fNumberOfSteps = 2000;
  
  AlphaPee[pp]  = 4.68  ;
  AlphaPee[pep] = 5.13  ;
  AlphaPee[hep] = 3.96  ;
  AlphaPee[be7] = 6.16  ;
  AlphaPee[b8]  = 6.81  ;
  AlphaPee[n13] = 6.22  ;
  AlphaPee[o15] = 6.69  ;
  AlphaPee[f17] = 6.74  ;

  BetaPee[pp]  =  0.109  ;
  BetaPee[pep] =  0.079  ;
  BetaPee[hep] =  0.165 ;
  BetaPee[be7] =  0.029 ;
  BetaPee[b8]  =  0.01 ;
  BetaPee[n13] =  0.054 ;
  BetaPee[o15] =  0.013 ;
  BetaPee[f17] =  0.012 ;

  NuEndPoint[pp]  =  0.423  ;
  NuEndPoint[pep] =  1.44  ;
  NuEndPoint[hep] =  19.0 ;
  NuEndPoint[be7] =  9.0 ;
  NuEndPoint[b8]  =  16.56;
  NuEndPoint[n13] =  1.199 ;
  NuEndPoint[o15] =  1.732 ;
  NuEndPoint[f17] =  1.74 ;


  fVolumeFlag = false ;
  G4ThreeVector zero(0., 0., 0.) ;
  fSPSAng = new G4SPSAngDistribution ;
  fSPSPos = new G4SPSPosDistribution ;
  G4SPSRandomGenerator *RndGen = new G4SPSRandomGenerator;
  fSPSPos->SetBiasRndm(RndGen);
  fSPSAng->SetBiasRndm(RndGen);
  SetRadius0(0.);
  SetCentreCoords(zero);    
  
  fParticleGun  = new G4ParticleGun ;
  fParticleGun->SetParticleDefinition(G4Electron::Definition());

  fTheMessenger = new BxGeneratorSolarNeutrinoMessenger(this);
 
  if(!fParticleGun) {
    BxLog(error) << "Could not allocate G4ParticleGun! Out of memory?"<<endlog;
    BxLog(fatal) << endlog;
  }
  BxLog(routine) << "G4ParticleGun Constructed." << endlog;
  

}

//---------------------------------------------------------------------------//

BxGeneratorSolarNeutrino::~BxGeneratorSolarNeutrino()
{
  delete fTheMessenger;
  delete fParticleGun;
  delete fSPSPos;
  delete fSPSAng;
  
}

//---------------------------------------------------------------------------//


void BxGeneratorSolarNeutrino::BxGeneratePrimaries(G4Event *event) {
  if(isFirstTime) {
    if(fNeutrinoType < 0) {
      BxLog(error) << "Set the neutrino source type"<<endlog;
      BxLog(fatal) << endlog;
    }  
    DM2   *= 1.E-5;
    BxLog(routine) << "deltaM^2 = " << DM2 << " eV" << endlog ;
    BxLog(routine) << "tan(theta)^2 = " << TG2T << endlog ;
    BxLog(routine) << "Neutrino source " << fNeutrinoType << endlog ;
    CumulativeDistribution(fNeutrinoType);
    isFirstTime = false ;
  }
  if(fVolumeFlag)  {
	  fSPSAng->SetVerbosity(0);
	  SetParticlePosition(fSPSPos->GenerateOne());  
	  fSPSAng->SetAngDistType("iso");    
	  fSPSAng->SetPosDistribution(fSPSPos);
	  SetParticleMomentumDirection(fSPSAng->GenerateOne());    
  }

  fParticleGun->SetNumberOfParticles(1);
  fParticleGun->SetParticleEnergy(ShootEnergy()/MeV);
  fParticleGun->GeneratePrimaryVertex(event);
  


  BxOutputVertex::Get()->SetPosition(fParticleGun->GetParticlePosition());  
  BxOutputVertex::Get()->SetDirection(fParticleGun->GetParticleMomentumDirection()); 
  BxOutputVertex::Get()->SetPDG(fParticleGun->GetParticleDefinition()->GetPDGEncoding());
  BxOutputVertex::Get()->SetTime(0.);
  BxOutputVertex::Get()->SetEnergy(fParticleGun->GetParticleEnergy()/MeV);		

}
//-------------------------------------------------------------------------
//     Survival probability Pee
//     from: P.C. de Holanda, Wei Liao and A. Yu. Smirnov (hep-ph/0404042)
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino::survive(G4double E, G4int k) {
  
	G4double theta  =  atan(sqrt(TG2T));
	G4double sin22  =  pow(sin(2*theta),2.);
	G4double cos2   =  sqrt(1.-sin22) ;
	G4double gamma  =  2.*AlphaPee[k]*E*1.E6/DM2*1.E-12 ;

  G4double delta  = 3./2.*pow(gamma,2.)*sin22*BetaPee[k]
	  /pow(pow(cos2 - gamma,2.)+sin22,2.);

  G4double cos2m  = (cos2 - gamma)/sqrt(pow(cos2-gamma,2.)+sin22);

  return 0.5 + 0.5*(1.-delta)*cos2m*cos2 ;  
}


G4double BxGeneratorSolarNeutrino::I_f(G4double E) {
	G4double x=sqrt(1+2*emass/E);
	return 1./6.*(1./3.+(3.-x*x)*(1./2.*x*log((x+1.)/(x-1.))-1.));
}

G4double BxGeneratorSolarNeutrino::k_e(G4double E) {
	return 0.9786+0.0097*I_f(E);
}
G4double BxGeneratorSolarNeutrino::k_mu(G4double E) {
	return 0.9965+0.00037*I_f(E);
}

G4double BxGeneratorSolarNeutrino::fMinus(G4double z, G4double q){
	G4double E=z*q+emass;
	G4double l=sqrt(E*E-emass*emass);
	TF1 fFun("fFun","log(abs(1-x))/x",0.0,E);
	G4double val= (E/l*log((E+l)/emass)-1.)*(2*log(1.-z-emass/(E+l))-log(1.-z)-1./2.*log(z)-5./12.)+1./2.*(fFun.Integral(0,z)-fFun.Integral(0,l/E))-1./2.*log(1.-z)*log(1.-z)-(11./12.+z/2.)*log(1-z)+z*(log(z)+1/2.*log(2*q/emass))-(31./18.+1./12.*log(z))*l/E-11./12.*z+z*z/24.;
if(val!=val)
return 0;
else
return val;
}


G4double BxGeneratorSolarNeutrino::fPlus(G4double z, G4double q){
	G4double E=z*q+emass;
	G4double l=sqrt(E*E-emass*emass);
	TF1 fFun("fFun","log(abs(1-x))/x",0.0,E);
	G4double val= 1./(1.-z)/(1.-z)*(E/l*log((E+l)/emass)-1.)*((1-z)*(1-z)*(2*log(1-z-emass/(E+l))-log(1-z)-log(z)/2.-2./3.)-(z*z*log(z)+1-z)/2)-(1-z)*(1-z)/2.*(log(1-z)*log(1-z)+l/E*(fFun.Integral(0,1.-z)-log(z)*log(1-z)))+log(1-z)*(z*z/2.*log(z)+(1-z)/3.*(2*z-1./2.))-z*z/2.*fFun.Integral(0,1-z)-z*(1-2*z)/3.*log(z)-z*(1-z)/6.-l/E/12.*(log(z)+(1-z)*(115-109*z)/6.);
if(val==val)
return val;
else
return 0;
}


G4double BxGeneratorSolarNeutrino::fPlusMinus(G4double z, G4double q){
  G4double E=z*q+emass;
  G4double l=sqrt(E*E-emass*emass);
  TF1 fFun("fFun","log(abs(1-x))/x",0.0,E);
  G4double val= (E/l*log((E+l)/emass)-1.)*2.*log(1.-z-emass/(E+l));
  if(val==val)
	  return val;
  else
	  return 0;
}


//-------------------------------------------------------------------------
//     Electroweak interaction nu_e + e-
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino::cross_nue(G4double E, G4double bin) {
	G4double endpoint = bin;
	G4double sigma = 0;

	if(!fOldSigma){
		//from K. Nakamura et al. (Particle Data Group), J. Phys. G, 37, 075021 (2010)
		G4double rhoNC=1.0127;
		G4double sintw=0.23116;
		//from Bahcall et al. PhysRevD. 51, 11 1995
		G4double gl=rhoNC*(1./2.-k_e(E)*sintw*sintw)-1;
		G4double gr=-rhoNC*k_e(E)*sintw*sintw;
		G4double z= E/endpoint;
		G4double alpha=7.297e-3;
		if(E < endpoint/(1.+emass/2./endpoint))
			sigma=gl*gl*(1.+alpha/M_PI*fMinus(z,endpoint))+gr*gr*(1.-z)*(1.-z)*(1.+alpha/M_PI*fPlus(z,endpoint))-gr*gl*emass*z/endpoint*(1+alpha/M_PI*fPlusMinus(z,endpoint));
		if(sigma < 0) sigma = 0 ;
	}else{
		G4double gr = 0.23 ;
		G4double gl = 0.5 + gr ;
		if(E < endpoint/(1+emass/2./endpoint))
			sigma = gl*gl+gr*gr*pow(1.-E/endpoint,2.)
				- gl*gr*E*emass/pow(endpoint,2) ;
		if(sigma < 0) sigma = 0 ;

	}
	return sigma ;
}

//-------------------------------------------------------------------------
//     Electroweak interaction nu_mu + e- or nu_tau + e-
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino::cross_nux(G4double E, G4double bin) {
	G4double endpoint = bin;
	G4double sigma = 0;
	if(!fOldSigma){
		//from K. Nakamura et al. (Particle Data Group), J. Phys. G, 37, 075021 (2010) 
		G4double rhoNC=1.0127;
		G4double sintw=0.23116;
		//from Bahcall et al. PhysRevD. 51, 11 1995
		G4double gl=rhoNC*(1./2.-k_mu(E)*sintw*sintw);
		G4double gr=-rhoNC*k_mu(E)*sintw*sintw;
		G4double z= E/endpoint;
		G4double alpha=7.297e-3;

		if(E < endpoint/(1+emass/2./endpoint))
			sigma=gl*gl*(1.+alpha/M_PI*fMinus(z,endpoint))+gr*gr*(1.-z)*(1.-z)*(1.+alpha/M_PI*fPlus(z,endpoint))-gr*gl*emass*z/endpoint*(1+alpha/M_PI*fPlusMinus(z,endpoint));
		if(sigma < 0) sigma = 0 ;
	}else{
		G4double gr = 0.23 ;
		G4double gl = - 0.5 + gr ;
		if(E < endpoint/(1+emass/2./endpoint))
			sigma = gl*gl+gr*gr*pow(1.-E/endpoint,2.)
				- gl*gr*E*emass/pow(endpoint,2) ;
		if(sigma < 0) sigma = 0 ;
	}
	return sigma ;
}

//-------------------------------------------------------------------------
//     N13
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino::N13NuSpectrum(G4double E) {
	G4double value = 0.;
	if(E < NuEndPoint[n13])  
		value =  3.1436*E*E*(1.71-E)*sqrt(pow(1.71-E,2.)-0.26112);
	return value ;
}


//-------------------------------------------------------------------------
//     O15
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino::O15NuSpectrum(G4double E) {
  G4double value = 0.;
  if(E < NuEndPoint[o15])  
   value =  0.677966*E*E*(2.243-E)*sqrt(pow(2.243-E,2)-0.26112);
  return value ;
}

//-------------------------------------------------------------------------
//     F17
//-------------------------------------------------------------------------
// Only 0.04% of the main CNO branch
G4double BxGeneratorSolarNeutrino::F17NuSpectrum(G4double E) {
  G4double value = 0.;
  if(E < NuEndPoint[f17]) 
    value =  0.6680926*E*E*(2.251-E)*sqrt(pow(2.251-E,2)-0.26112);
  return value ;
}

//-------------------------------------------------------------------------
//     B8
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino::B8NuSpectrum(G4double E) {
  G4double value = 0.;
  if(E < NuEndPoint[b8])  
    value = 8.0496e-7*E*E*pow(17.071-E,3.4899);
  return value ;
}

//-------------------------------------------------------------------------
//     PP
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino::PPNuSpectrum(G4double E) {
  G4double value = 0.;
  if(E < NuEndPoint[pp])  
    value = 188.981*E*E*(0.9339-E)*sqrt(pow(0.9339-E,2)-0.26112);
  return value ;
}
//-------------------------------------------------------------------------
//    Energy Spectrum
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino::EnergySpectrum(G4double E, G4int k) {
  G4double endpoint = NuEndPoint[k] ;
  G4double step = endpoint/double(fNumberOfSteps);
  G4double sum = 0;
  if(k == pep) {
    if(E < endpoint ) 
      sum =  survive(E,k)*cross_nue(E,endpoint)+(1.-survive(E,k))*cross_nux(E,endpoint); 
    return sum;
  }  else if(k == be7) {
  
//    G4double endpoint2 = 0.3843;
    
//    G4double sum1 =  survive(E,k)*cross_nue(E,endpoint)
//       +(1.-survive(E,k))*cross_nux(E,endpoint); 
       
      G4double sum1 =  cross_nue(E,endpoint);      
       
       
 //   G4double sum2 =  survive(E,k)*cross_nue(E,endpoint2)
 //      +(1.-survive(E,k))*cross_nux(E,endpoint2);  
    
 //   return  0.897*sum1 + 0.103*sum2 ;
 
      return  sum1 ;
  }

  G4double cross_section = 0;
  for(G4int i = 0; i< fNumberOfSteps ; i++) {
    G4double bin = G4double(i)*step + step/3. ;
    if(k == pp) {
      cross_section = PPNuSpectrum(bin)*
	  (survive(bin,k)*cross_nue(E,bin) + (1. - survive(bin,k))*cross_nux(E,bin));                   
    } else if(k == b8) {
       cross_section = B8NuSpectrum(bin)*
	  (survive(bin,k)*cross_nue(E,bin) + (1. - survive(bin,k))*cross_nux(E,bin));                   
    } else if(k == f17) {
       cross_section = F17NuSpectrum(bin)*
	  (survive(bin,k)*cross_nue(E,bin) + (1. - survive(bin,k))*cross_nux(E,bin));                   
    } else if(k == o15) {
       cross_section = O15NuSpectrum(bin)*
	  (survive(bin,k)*cross_nue(E,bin) + (1. - survive(bin,k))*cross_nux(E,bin));                   
    } else if(k == n13) {
       cross_section = N13NuSpectrum(bin)*
	  (survive(bin,k)*cross_nue(E,bin) + (1. - survive(bin,k))*cross_nux(E,bin));                   
    }  
    sum += cross_section*step ;
  }

  return sum ;
}
//-------------------------------------------------------------------------
//    Energy Spectrum Normalizer
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino::Normalizer(G4int k) {
  G4double endpoint = NuEndPoint[k]/(1+emass/2./NuEndPoint[k]) ;
  G4double step = endpoint/G4double(fNumberOfSteps);
  G4double sum = 0;
  for(G4int i = 0; i< fNumberOfSteps ; i++) {
    G4double bin = G4double(i)*step + step/3. ;
    sum += EnergySpectrum(bin,k)*step ;
  }
  return sum ;
}
//-------------------------------------------------------------------------
//    Cumulative Distribution
//-------------------------------------------------------------------------
void BxGeneratorSolarNeutrino::CumulativeDistribution(G4int k) {
  G4double sum = 0;
  G4double norma, endpoint, normaN13, normaO15, normaF17, step ;
  if(k != cno) {
    norma = Normalizer(k);
    endpoint = NuEndPoint[k]/(1+emass/2./NuEndPoint[k]) ;
    step = endpoint/G4double(fNumberOfSteps);
    for(G4int i = 0; i< fNumberOfSteps ; i++) {
      G4double bin = G4double(i)*step + step/3. ;
      sum += EnergySpectrum(bin,k)*step/norma ;
      fEnergyBin.push_back(bin);
      fProbability.push_back(sum);
    }
  } else {
    normaN13 = Normalizer(n13);
    normaO15 = Normalizer(o15);
    normaF17 = Normalizer(f17);
    endpoint = NuEndPoint[f17]/(1+emass/2./NuEndPoint[f17]) ;
    step = endpoint/G4double(fNumberOfSteps);
    for(G4int i = 0; i< fNumberOfSteps ; i++) {
      G4double bin = G4double(i)*step + step/3. ;
      sum += 0.9996*(2.00*EnergySpectrum(bin,n13)/normaN13 
             +1.44*EnergySpectrum(bin,o15)/normaO15)/3.44 
             +4.E-4*EnergySpectrum(bin,f17)/normaF17; 	     
      fEnergyBin.push_back(bin);
      fProbability.push_back(sum*step);
    }
  }
}
//-------------------------------------------------------------------------
//    Energy Shooter with linear interpolation
//-------------------------------------------------------------------------

G4double BxGeneratorSolarNeutrino::ShootEnergy() {
  G4double val = G4UniformRand()  ;
  for(G4int i=0;i<G4int(fProbability.size());i++) {
    if(fProbability[i] >= val) {
      if(i == 0) return fEnergyBin[0];
      G4double deltaX = val - fProbability[i] ;
      G4double y = fEnergyBin[i] - fEnergyBin[i-1] ;
      G4double x = fProbability[i] - fProbability[i-1] ;
      return deltaX*y/x + fEnergyBin[i];
    }
  }
  return 0;
}
