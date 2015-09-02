// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//
#include "TF1.h"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4Electron.hh"
#include <CLHEP/Random/RandFlat.h>
#include "BxLogger.hh"
#include "BxDataCollection.hh"
#include "BxOutputVertex.hh"
#include "BxGeneratorSolarNeutrino2.hh"
#include "G4SystemOfUnits.hh"
//---------------------------------------------------------------------------//


BxGeneratorSolarNeutrino2::BxGeneratorSolarNeutrino2(): BxVGenerator("BxGeneratorSolarNeutrino2") {
  //pi = 3.14159265;
  emass  = G4Electron::Definition()->GetPDGMass ()/MeV;
//according to pp-nu internal note April 16, 2014
  DM2    = 7.54 ;
  TG2T   = 0.44 ;
  isFirstTime = true ;
  is_read     = false ;
  fOldSigma=false; //neutrino cross-sections at tree-level
  fNeutrinoType = -1;

  fNumberOfSteps = 10000;

  AlphaPee[pp]  = 4.68  ;
  AlphaPee[pep] = 5.13  ;
  AlphaPee[hep] = 3.96  ;
  AlphaPee[be7] = 6.16  ;
  AlphaPee[b8]  = 6.81  ;
  AlphaPee[n13] = 6.22  ;
  AlphaPee[o15] = 6.69  ;
  AlphaPee[f17] = 6.74  ;
  AlphaPee[be7_862] = AlphaPee[be7]  ;
  AlphaPee[be7_384] = AlphaPee[be7]  ;

  BetaPee[pp]  =  0.109  ;
  BetaPee[pep] =  0.079  ;
  BetaPee[hep] =  0.165 ;
  BetaPee[be7] =  0.029 ;
  BetaPee[b8]  =  0.01 ;
  BetaPee[n13] =  0.054 ;
  BetaPee[o15] =  0.013 ;
  BetaPee[f17] =  0.012 ;
  BetaPee[be7_862] = BetaPee[be7] ;
  BetaPee[be7_384] = BetaPee[be7] ;

  NuEndPoint[pp]  =  0.423  ;
  NuEndPoint[pep] =  1.44  ;
  NuEndPoint[hep] =  19.0 ;
  NuEndPoint[be7] =  0.862 ;
  NuEndPoint[b8]  =  16.56;
  NuEndPoint[n13] =  1.199 ;
  NuEndPoint[o15] =  1.732 ;
  NuEndPoint[f17] =  1.74 ;
  NuEndPoint[be7_862] =  0.862 ;
  NuEndPoint[be7_384] =  0.3843 ;
  
  fnbin  = 4000;
  fNormaNue = 0;
  fNormaNux = 0;
  
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

  fTheMessenger = new BxGeneratorSolarNeutrino2Messenger(this);
 
  if(!fParticleGun) {
    BxLog(error) << "Could not allocate G4ParticleGun! Out of memory?"<<endlog;
    BxLog(fatal) << endlog;
  }
  BxLog(routine) << "G4ParticleGun Constructed." << endlog;
  

}

//---------------------------------------------------------------------------//

BxGeneratorSolarNeutrino2::~BxGeneratorSolarNeutrino2()
{
  delete fTheMessenger;
  delete fParticleGun;
  delete fSPSPos;
  delete fSPSAng;
  
}

//---------------------------------------------------------------------------//


void BxGeneratorSolarNeutrino2::BxGeneratePrimaries(G4Event *event) {
  if(isFirstTime) {
    if(fNeutrinoType < 0) {
      BxLog(error) << "Set the neutrino source type"<<endlog;
      BxLog(fatal) << endlog;
    }  
    DM2   *= 1.E-5;
    BxLog(routine) << "deltaM^2 = " << DM2 << " eV^2" << endlog ;
    BxLog(routine) << "tan(theta)^2 = " << TG2T << endlog ;
    BxLog(routine) << "Neutrino source " << fNeutrinoType << endlog ;
    CumulativeDistribution(fNeutrinoType);
    fstep =  (NuEndPoint[fNeutrinoType]*1.05)/G4double(fnbin);
    for(G4int k=0; k<fnbin;k++) {
G4cout<<"k = "<<k<<G4endl;
      G4double norma_nue = 0; 
      G4double norma_nux = 0; 
      for(G4int i=0;i<fnbin;i++) {
      	norma_nue += cross_nue(fstep*i+fstep/3.,fstep*k+fstep/3.)*fstep;
      	norma_nux += cross_nux(fstep*i+fstep/3.,fstep*k+fstep/3.)*fstep;
	}
      if(fNormaNue < norma_nue) fNormaNue = norma_nue;
      if(fNormaNux < norma_nux) fNormaNux = norma_nux;
    }
G4cout<<"fNormaNue = "<<fNormaNue<<G4endl;
G4cout<<"fNormaNux = "<<fNormaNux<<G4endl;

    isFirstTime = false ;
  }
  if(fVolumeFlag)  {
    fSPSAng->SetVerbosity(0);
    SetParticlePosition(fSPSPos->GenerateOne());  
    fSPSAng->SetAngDistType("iso");    
    fSPSAng->SetPosDistribution(fSPSPos);
    SetParticleMomentumDirection(fSPSAng->GenerateOne());    
  }
  G4double EEnergy = 0;
  G4double NuEnergy ; 

  G4int nutype = 0; // nutype = 0 for nu_e and 1 for nu_mu or nu_tau
  bool isGenerated = false;
  while(!isGenerated) {
    NuEnergy  =  ShootEnergy(fNeutrinoType);
    G4double probx = 1-survive(NuEnergy, fNeutrinoType) ;
    // here I choose if it is a nue or nux
    if(G4UniformRand() < probx) {  // nux
      nutype = 1; 
      EEnergy = 0;
      G4double sum = 0;
      G4double value = CLHEP::RandFlat::shoot(0.,1.);
      if(value < cross_nux(0,NuEnergy)*fstep/fNormaNux) {
	  isGenerated = false ; 
	  break ;        
      } 
      for(G4int i=0;i<fnbin;i++) {  // generation of the recoiled electron pdf
	EEnergy = fstep*i+fstep/3.;
	sum += cross_nux(EEnergy,NuEnergy)*fstep/fNormaNux;
	if(cross_nux(EEnergy,NuEnergy) == 0) {
	  isGenerated = false ; 
	  break ;
	}
	if(sum > value) {
	  isGenerated = true ;
	  if(G4UniformRand()> fNormaNux/fNormaNue) isGenerated = false;	   
	  break;	  
        }
      }
    } else {    // nue
      nutype = 0; 
      EEnergy = 0;
      G4double sum = 0;
      G4double value = CLHEP::RandFlat::shoot(0.,1.);
      if(value < cross_nue(0,NuEnergy)*fstep/fNormaNue) {
	  isGenerated = false ; 
	  break ;        
      }
      for(G4int i=0;i<fnbin;i++) { // generation of the recoiled electron pdf
	EEnergy = fstep*i+fstep/3.;
	sum += cross_nue(EEnergy,NuEnergy)*fstep/fNormaNue;
	if(cross_nue(EEnergy,NuEnergy) == 0) {
	  isGenerated = false ; 
	  break ;
	}
	if(sum > value) {
	  isGenerated = true ; 
	  break;	  
        }      
      }
    }        
  } 
  
  BxOutputVertex::Get()->SetPosition(fParticleGun->GetParticlePosition());  
  BxOutputVertex::Get()->SetDirection(fParticleGun->GetParticleMomentumDirection()); 
  BxOutputVertex::Get()->SetPDG(fParticleGun->GetParticleDefinition()->GetPDGEncoding());
  BxOutputVertex::Get()->SetTime(0.);
  BxOutputVertex::Get()->SetEnergy(NuEnergy/MeV);		
  BxOutputVertex::Get()->SetDEnergy(EEnergy/MeV);		
  if(!nutype) BxOutputVertex::Get()->SetPDG(12);
  else        BxOutputVertex::Get()->SetPDG(14);
  BxOutputVertex::Get()->SetDPDG(fParticleGun->GetParticleDefinition()->GetPDGEncoding());		
  BxOutputVertex::Get()->SetDaughters();


  fParticleGun->SetNumberOfParticles(1);
  fParticleGun->SetParticleEnergy(EEnergy/MeV);
  fParticleGun->GeneratePrimaryVertex(event);
}
//-------------------------------------------------------------------------
//     Survival probability Pee
//     from: P.C. de Holanda, Wei Liao and A. Yu. Smirnov (hep-ph/0404042)
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino2::survive(G4double E, G4int k) {
  if(TG2T == 0)  return 1 ;
  if(DM2  == 0)  return 1 ;
  G4double theta  =  atan(sqrt(TG2T));
  G4double sin22  =  pow(sin(2*theta),2.);
  G4double cos2   =  sqrt(1.-sin22) ;
  G4double gamma  =  2.*AlphaPee[k]*E*1.E6/DM2*1.E-12 ;
  
  G4double delta  = 3./2.*pow(gamma,2.)*sin22*BetaPee[k]
                  /pow(pow(cos2 - gamma,2.)+sin22,2.);
		  
  G4double cos2m  = (cos2 - gamma)/sqrt(pow(cos2-gamma,2.)+sin22);
  
  return 0.5 + 0.5*(1.-delta)*cos2m*cos2 ;  
}

G4double BxGeneratorSolarNeutrino2::I_f(G4double E) {
        G4double x=sqrt(1+2*emass/E);
        return 1./6.*(1./3.+(3.-x*x)*(1./2.*x*log((x+1.)/(x-1.))-1.));
}

G4double BxGeneratorSolarNeutrino2::k_e(G4double E) {
        return 0.9786+0.0097*I_f(E);
}
G4double BxGeneratorSolarNeutrino2::k_mu(G4double E) {
        return 0.9965+0.00037*I_f(E);
}

G4double BxGeneratorSolarNeutrino2::fMinus(G4double z, G4double q){
	G4double E=z*q+emass;
	G4double l=sqrt(E*E-emass*emass);
	TF1 fFun("fFun","log(abs(1-x))/x",0.0,E);
	G4double val= (E/l*log((E+l)/emass)-1.)*(2*log(1.-z-emass/(E+l))-log(1.-z)-1./2.*log(z)-5./12.)+1./2.*(fFun.Integral(0,z)-fFun.Integral(0,l/E))-1./2.*log(1.-z)*log(1.-z)-(11./12.+z/2.)*log(1-z)+z*(log(z)+1/2.*log(2*q/emass))-(31./18.+1./12.*log(z))*l/E-11./12.*z+z*z/24.;
	if(val!=val)
		return 0;
	else
		return val;
}


G4double BxGeneratorSolarNeutrino2::fPlus(G4double z, G4double q){
        G4double E=z*q+emass;
	G4double l=sqrt(E*E-emass*emass);
	TF1 fFun("fFun","log(abs(1-x))/x",0.0,E);
	G4double val= 1./(1.-z)/(1.-z)*(E/l*log((E+l)/emass)-1.)*((1-z)*(1-z)*(2*log(1-z-emass/(E+l))-log(1-z)-log(z)/2.-2./3.)-(z*z*log(z)+1-z)/2)-(1-z)*(1-z)/2.*(log(1-z)*log(1-z)+l/E*(fFun.Integral(0,1.-z)-log(z)*log(1-z)))+log(1-z)*(z*z/2.*log(z)+(1-z)/3.*(2*z-1./2.))-z*z/2.*fFun.Integral(0,1-z)-z*(1-2*z)/3.*log(z)-z*(1-z)/6.-l/E/12.*(log(z)+(1-z)*(115-109*z)/6.);
	if(val==val)
		return val;
	else
		return 0;
}

G4double BxGeneratorSolarNeutrino2::fPlusMinus(G4double z, G4double q){
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
G4double BxGeneratorSolarNeutrino2::cross_nue(G4double E, G4double nuE) {
	//from K. Nakamura et al. (Particle Data Group), J. Phys. G, 37, 075021 (2010)
	G4double sigma = 0;
	if(!fOldSigma){
		G4double rhoNC=1.0127;
		G4double sintw=0.23116;
		//from Bahcall et al. PhysRevD. 51, 11 1995
		G4double gl=rhoNC*(1./2.-k_e(E)*sintw*sintw)-1;
		G4double gr=-rhoNC*k_e(E)*sintw*sintw;
		G4double z= E/nuE;
		G4double alpha=7.297e-3;
		if(E < nuE/(1.+emass/2./nuE))
			sigma=gl*gl*(1.+alpha/M_PI*fMinus(z,nuE))+gr*gr*(1.-z)*(1.-z)*(1.+alpha/M_PI*fPlus(z,nuE))-gr*gl*emass*z/nuE*(1+alpha/M_PI*fPlusMinus(z,nuE));
		if(sigma < 0) sigma = 0 ;
	}else{
		G4double gr = 0.23 ;
		G4double gl = 0.5 + gr ;
		if(E > nuE/(1+emass/2./nuE)) return 0 ;
		else  sigma = gl*gl+gr*gr*pow(1.-E/nuE,2.) - gl*gr*E*emass/nuE/nuE ;

		if(sigma < 0) sigma = 0 ;
	}
	return sigma ;
}

//-------------------------------------------------------------------------
//     Electroweak interaction nu_mu + e- or nu_tau + e-
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino2::cross_nux(G4double E, G4double nuE) {
	G4double sigma = 0;
	if(!fOldSigma){
		//from K. Nakamura et al. (Particle Data Group), J. Phys. G, 37, 075021 (2010) 
		G4double rhoNC=1.0127;
		G4double sintw=0.23116;
		//                        //from Bahcall et al. PhysRevD. 51, 11 1995
		G4double gl=rhoNC*(1./2.-k_mu(E)*sintw*sintw);
		G4double gr=-rhoNC*k_mu(E)*sintw*sintw;
		G4double z= E/nuE;
		G4double alpha=7.297e-3;

	if(E < nuE/(1+emass/2./nuE))
		sigma=gl*gl*(1.+alpha/M_PI*fMinus(z,nuE))+gr*gr*(1.-z)*(1.-z)*(1.+alpha/M_PI*fPlus(z,nuE))-gr*gl*emass*z/nuE*(1+alpha/M_PI*fPlusMinus(z,nuE));
	if(sigma < 0) sigma = 0 ;
	}else{
		G4double gr = 0.23 ;
		G4double gl = - 0.5 + gr ;
		if(E > nuE/(1+emass/2./nuE)) return 0 ;
		else   sigma = gl*gl+gr*gr*pow(1.-E/nuE,2.) - gl*gr*E*emass/nuE/nuE ;

		if(sigma < 0) sigma = 0 ;
	}
	return sigma ;

}

//-------------------------------------------------------------------------
//     N13
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino2::N13NuSpectrum(G4double E) {
  if(!is_read) {
    is_read = true;
    G4double ene, prob;
//from Bahcall:http://www.sns.ias.edu/~jnb/SNdata/Export/CNOspectra/n13.dat
    ifstream fil("../g4bx2/data/dat/n13spectrum.dat");
    while(!fil.eof()) {
      fil >> ene >> prob;
      if(fil.eof()) break ;
      fNuEnergy.push_back(ene);
      fNuProbability.push_back(prob);
    }
    fil.close();
  }
  if(E>NuEndPoint[n13]) return 0.;
  G4double bin = fNuEnergy[1] - fNuEnergy[0];
  G4int    k = G4int(E/bin) ;
  return fNuProbability[k];
}


//-------------------------------------------------------------------------
//     O15
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino2::O15NuSpectrum(G4double E) {
  if(!is_read) {
    is_read = true;
    G4double ene, prob;
//from Bahcall: http://www.sns.ias.edu/~jnb/SNdata/Export/CNOspectra/o15.dat
    ifstream fil("../g4bx2/data/dat/o15spectrum.dat");
    while(!fil.eof()) {
      fil >> ene >> prob;
      if(fil.eof()) break ;
      fNuEnergy.push_back(ene);
      fNuProbability.push_back(prob);
    }
    fil.close();
  }

  if(E>NuEndPoint[o15]) return 0.;  
  G4double bin = fNuEnergy[1] - fNuEnergy[0];
  G4int    k = G4int(E/bin) ;
  return fNuProbability[k];
}

//-------------------------------------------------------------------------
//     F17
//-------------------------------------------------------------------------
// Only 0.04% of the main CNO branch
G4double BxGeneratorSolarNeutrino2::F17NuSpectrum(G4double E) {
  if(!is_read) {
    is_read = true;
    G4double ene, prob;
//from Bahcall: http://www.sns.ias.edu/~jnb/SNdata/Export/CNOspectra/f17.dat
    ifstream fil("../g4bx2/data/dat/f17spectrum.dat");
    while(!fil.eof()) {
      fil >> ene >> prob;
      if(fil.eof()) break ;
      fNuEnergy.push_back(ene);
      fNuProbability.push_back(prob);
    }
    fil.close();
  }
  
  if(E>NuEndPoint[f17]) return 0.;  
  G4double bin = fNuEnergy[1] - fNuEnergy[0];
  G4int    k = G4int(E/bin) ;
  return fNuProbability[k];
}

//-------------------------------------------------------------------------
//     B8
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino2::B8NuSpectrum(G4double E) {
  if(!is_read) {
    is_read = true;
    G4double ene, prob;
    //from:W. Winter et al., PHYSICAL REVIEW C 73, 025503 (2006)
    ifstream fil("../g4bx2/data/dat/b8spectrum.dat");
    while(!fil.eof()) {
	    fil >> ene >> prob;
	    if(fil.eof()) break ;
	    fNuEnergy.push_back(ene);
	    fNuProbability.push_back(prob);
    }
    fil.close();
  }

  if(E>NuEndPoint[b8]) return 0.;  
  G4double bin = fNuEnergy[1] - fNuEnergy[0];
  G4int    k = G4int(E/bin) ;
  return fNuProbability[k];
}

//-------------------------------------------------------------------------
//     PP
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino2::PPNuSpectrum(G4double E) {
  if(!is_read) {
    is_read = true;
    G4double ene, prob;
//from Bahcall: http://www.sns.ias.edu/~jnb/SNdata/Export/PPenergyspectrum/ppenergytab
    ifstream fil("../g4bx2/data/dat/ppspectrum.dat");
    while(!fil.eof()) {
      fil >> ene >> prob;
G4cout<<" ene= "<<ene<<" prob = "<<prob<<G4endl;
      if(fil.eof()) break ;
      fNuEnergy.push_back(ene);
      fNuProbability.push_back(prob);
    }
    fil.close();
  }
  
  if(E>NuEndPoint[pp]) return 0.;  
  G4double bin = fNuEnergy[1] - fNuEnergy[0];
  G4int    k = G4int(E/bin) ;
  return fNuProbability[k];
}
//-------------------------------------------------------------------------
//     Be7
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino2::BE7NuSpectrum(G4double E) {
  if(!is_read) {
    is_read = true;
    G4double ene, prob;
    ifstream fil("../g4bx2/data/dat/be7spectrum.dat");
    while(!fil.eof()) {
      fil >> ene >> prob;
      if(fil.eof()) break ;
      fNuEnergy.push_back(ene);
      fNuProbability.push_back(prob);
    }
    fil.close();
  }
  
  if(E>NuEndPoint[be7]) return 0.;  
  G4double bin = fNuEnergy[1] - fNuEnergy[0];
  G4int    k = G4int(E/bin) ;
  return fNuProbability[k];
}
//-------------------------------------------------------------------------
//     HEP
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino2::HEPNuSpectrum(G4double E) {
  if(!is_read) {
    is_read = true;
    G4double ene, prob;
//from Bahcall: http://www.sns.ias.edu/~jnb/SNdata/Export/Hepspectrum/hepspectrum.dat
    ifstream fil("../g4bx2/data/dat/hepspectrum.dat");
    while(!fil.eof()) {
      fil >> ene >> prob;
      if(fil.eof()) break ;
      fNuEnergy.push_back(ene);
      fNuProbability.push_back(prob);
    }
    fil.close();
  }
  
  if(E>NuEndPoint[hep]) return 0.;  
  G4double bin = fNuEnergy[1] - fNuEnergy[0];
  G4int    k = G4int(E/bin) ;
  return fNuProbability[k];
}
//-------------------------------------------------------------------------
//    Energy Spectrum
//-------------------------------------------------------------------------
G4double BxGeneratorSolarNeutrino2::EnergySpectrum(G4double E, G4int k) {
  G4double value = 0;
  if(k == pp) {
    value = PPNuSpectrum(E);                   
  } else if(k == b8) {
     value = B8NuSpectrum(E);                   
  } else if(k == f17) {
     value = F17NuSpectrum(E);                   
  } else if(k == o15) {
     value = O15NuSpectrum(E);                   
  } else if(k == n13) {
     value = N13NuSpectrum(E);                   
  }  else if(k == be7) {
     //value = BE7NuSpectrum(E);                        
     if(E == 0.3843) value =  0.103;
     else if(E == NuEndPoint[be7]) value =  0.897;
     else value = 0;
  }  else if(k == be7_862) {                      
     if(E == NuEndPoint[be7_862]) value =  1;
     else value = 0;
  }  else if(k == be7_384) {                      
     if(E == NuEndPoint[be7_384]) value =  1;
     else value = 0;
  }  else if(k == hep) {
     value = HEPNuSpectrum(E);                        
  }  else if(k == pep) {
     if(E == NuEndPoint[pep]) value =  1;  
     else value = 0;
  }

  return value ;
}
//-------------------------------------------------------------------------
//    Cumulative Distribution
//-------------------------------------------------------------------------
void BxGeneratorSolarNeutrino2::CumulativeDistribution(G4int k) {
  G4double sum = 0;
  G4double endpoint = NuEndPoint[k] ;
  G4double step = endpoint/G4double(fNumberOfSteps);
 for(G4int i = 0; i< fNumberOfSteps ; i++) {
    G4double bin = G4float(i)*step + step/3. ;
    sum += EnergySpectrum(bin,k)*step ;
    fEnergyBin.push_back(bin);
    fProbability.push_back(sum);
  }
}
//-------------------------------------------------------------------------
//    Energy Shooter with linear interpolation
//-------------------------------------------------------------------------

G4double BxGeneratorSolarNeutrino2::ShootEnergy(G4int k) {
  G4double val = G4UniformRand()  ;
  if((k != pep) && (k != be7) && (k != be7_862) && (k != be7_384)) {
    for(G4int i=0;i<G4int(fProbability.size());i++) {
      if(fProbability[i] >= val) {
	if(i == 0) return fEnergyBin[0];
	G4double deltaX = val - fProbability[i] ;
	G4double y = fEnergyBin[i] - fEnergyBin[i-1] ;
	G4double x = fProbability[i] - fProbability[i-1] ;
	return deltaX*y/x + fEnergyBin[i];
      }
    }
  } else if(k == be7) {
    if(val < 0.897) return NuEndPoint[be7];
    else return 0.3843;
  } else if(k == be7_862) return NuEndPoint[be7_862];
  else if(k == be7_384) return NuEndPoint[be7_384];
  else if(k == pep) return NuEndPoint[pep];
  
  return 0;
}
//-------------------------------------------------------------------------
//    Pure interpolation
//-------------------------------------------------------------------------

G4double BxGeneratorSolarNeutrino2::Interpolate(G4double prob1, G4double prob2, G4double eE1, G4double eE2, G4double val) {
  if(eE1 < 0) eE1 = 0;
  G4double deltaX = val -  prob1;
  G4double y = eE2 - eE1 ;
  G4double x = prob2 - prob1 ;
  return deltaX*y/x + eE1;
}

