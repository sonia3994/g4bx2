// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//
#include "G4Poisson.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4SPSAngDistribution.hh"
#include "BxLogger.hh"
#include "BxReadParameters.hh"
#include "BxDataCollection.hh"
#include "G4SPSEneDistribution.hh"
#include "BxGeneratorSCS.hh"
#include "BxDataCollection.hh"
#include "BxGeneratorSCSMessenger.hh"
#include "BxOutputVertex.hh"
#include "BxPropertyCollection.hh"
#include <math.h>
#include <fstream>
#include <math.h>
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
using namespace std;
//---------------------------------------------------------------------------//


BxGeneratorSCS::BxGeneratorSCS(): BxVGenerator("BxGeneratorSCS"){
  

  ffp[0] =  0.99999999999980993;
  ffp[1] =  676.5203681218851;
  ffp[2] =  -1259.1392167224028; 
  ffp[3] =  771.32342877765313; 
  ffp[4] =  -176.61502916214059;
  ffp[5] =  12.507343278686905; 
  ffp[6] =  -0.13857109526572012; 
  ffp[7] =  9.9843695780195716e-6;
  ffp[8] =  1.5056327351493116e-7;

  fSize = 200 ;
  fC14ShapeFactor = 0.3 ;
  fShape = 1;
  fCrossSection = "bi210";
  fFirstTime = true ;
  theParticleTable = G4ParticleTable::GetParticleTable();  
  thePosition  = G4ThreeVector(0.0,0.0,0.0);
  theDirection = G4ThreeVector(0.0,0.0,1.0);
  fVolumeFlag    = false ;
  G4ThreeVector zero(0., 0., 0.) ;
  G4SPSRandomGenerator *RndGen = new G4SPSRandomGenerator;
  fSPSPos = new G4SPSPosDistribution ;
  fSPSAng = new G4SPSAngDistribution ;
  fSPSPos->SetBiasRndm(RndGen);
  fSPSAng->SetBiasRndm(RndGen);
  SetRadius0(0.);
  SetCentreCoords(zero);  
 
  fMessenger  = new  BxGeneratorSCSMessenger(this);

  BxLog(routine) << "Special Cross Section Generator Built" << endlog;
}

//---------------------------------------------------------------------------//

BxGeneratorSCS::~BxGeneratorSCS()
{
  delete fSPSPos;
  delete fSPSAng;
  delete fMessenger;
}

//---------------------------------------------------------------------------//

void BxGeneratorSCS::BxGeneratePrimaries(G4Event *event) {
   if(fFirstTime){
     if(fCrossSection == "c14") { 
       particle_definition = theParticleTable->FindParticle("e-");
       C14CrossSection(); 
     } else if(fCrossSection == "bi210") {
       particle_definition = theParticleTable->FindParticle("e-");
       Bi210CrossSection(); 
     } else if(fCrossSection == "c14oleg") {
       particle_definition = theParticleTable->FindParticle("e-");
       C14OlegCrossSection(); 
     } else  {
       BxLog(error) << "Wrong cross section" << endl ;
       BxLog(fatal) << endl ;
     }

     fFirstTime = false ;
   }

   if(fVolumeFlag)  {
     fSPSAng->SetAngDistType("iso"); 
     thePosition  = fSPSPos->GenerateOne();
     fSPSAng->SetPosDistribution(fSPSPos);
     theDirection = fSPSAng->GenerateOne();
   }
   
   G4double ene = ShootEnergy() ;
   G4double mass = particle_definition->GetPDGMass();
   G4double energy = ene + mass; 
   G4double pmom = std::sqrt(energy*energy-mass*mass);   
   G4double px = pmom*theDirection.x();
   G4double py = pmom*theDirection.y();
   G4double pz = pmom*theDirection.z();
   G4PrimaryVertex*   vertex   = new G4PrimaryVertex(thePosition,0.0);
   G4PrimaryParticle* particle = new G4PrimaryParticle(particle_definition,px,py,pz);
   vertex->SetPrimary( particle );
//G4cout<<ene<<G4endl;
   event->AddPrimaryVertex( vertex );
   BxOutputVertex::Get()->SetPosition(thePosition);  
   BxOutputVertex::Get()->SetDirection(theDirection); 
   BxOutputVertex::Get()->SetPDG(particle_definition->GetPDGEncoding());
   BxOutputVertex::Get()->SetTime(0.);
   BxOutputVertex::Get()->SetEnergy(ene);		


} 

//---------------------------------------------------------------------------//

void BxGeneratorSCS::C14CrossSection() {
   fFirstTime = false ;
   const G4double e_mass = 0.510998903;
   const G4double qvalue = 0.156475 ;
   G4double shape =  fC14ShapeFactor;
   G4double norma = 0; 
   G4double shape_term = 0;
   
   G4int nbin = fSize;
   for(G4int i=0;i<nbin;i++) {
     G4double xt = G4double(i+0.5)/G4double(nbin)*qvalue;
     G4double factor = 6.*(e_mass+xt)/(137.*sqrt(pow(xt,2)+2*e_mass*xt));
     G4double fermi2 = 2*M_PI*factor/(1-exp(-2.*M_PI*factor));
     if(fShape == 1)      shape_term = 1-shape*(xt+e_mass);     // Shape term type 0  
     else if(fShape == 2) shape_term = 1+shape*(qvalue - xt);	// Shape term type 1		        
     else if(fShape == 3) shape_term = 1-4.67*(xt+e_mass) + 3./(xt+e_mass) +	        
                          	       2.*pow(xt+e_mass,2);	// Shape term type 2
     else if(fShape == 4){ 
	     G4double correction=BxPropertyCollection::Get()->GetC14CorrectionFactor(i);
		     shape_term = correction*(1+shape*(qvalue-(0.5+i*(qvalue-0.5)/(nbin-1))));
     }
     else  BxLog(fatal) <<"Wrong shape factor term" << endl ;

     norma +=  sqrt(2*e_mass*xt+xt*xt)*(xt+e_mass)*pow(qvalue - xt,2.)
	     *shape_term*fermi2*qvalue/nbin;
   }

   for(G4int k=0;k<nbin;k++) {
	   G4double _fShootEnergy = 0;
	   G4double _fShootRandom = 0;

	   for(G4int i=0;i<k+1;i++) {
		   G4double xt = G4double(i+0.5)/G4double(nbin)*qvalue;
		   G4double factor = 6.*(e_mass+xt)/(137.*sqrt(pow(xt,2)+2*e_mass*xt));
		   G4double fermi2 = 2*M_PI*factor/(1-exp(-2.*M_PI*factor));
		   if(fShape == 1)      shape_term = 1-shape*(xt+e_mass);     // Shape term type 0  
		   else if(fShape == 2) shape_term = 1+shape*(qvalue - xt);	// Shape term type 1		        
		   else if(fShape == 3) shape_term = 1-4.67*(xt+e_mass) + 3./(xt+e_mass) +	        
			   2.*pow(xt+e_mass,2);	// Shape term type 2		        
     else if(fShape == 4){ 
	     G4double correction=BxPropertyCollection::Get()->GetC14CorrectionFactor(i);
		     shape_term = correction*(1+shape*(qvalue-(0.5+i*(qvalue-0.5)/(nbin-1))));
     }
		   _fShootEnergy = xt*MeV ;
		   _fShootRandom +=  sqrt(2*e_mass*xt+xt*xt)*(xt+e_mass)*pow(qvalue - xt,2.)
			   *shape_term*fermi2*qvalue/nbin/norma;
	   }

	   fShootEnergy.push_back(_fShootEnergy);
	   fShootRandom.push_back(_fShootRandom);
   } 
}
//---------------------------------------------------------------------------//

void BxGeneratorSCS::C14OlegCrossSection() {
  G4double array[60] = {
    9.82, 10, 10.5, 10.8, 11, 11.4, 11.5, 11.6, 11.7, 11.7,
    11.7, 11.6, 11.5, 11.4, 11.2, 11.1, 10.9, 10.7, 10.4, 10.2,
    9.94, 9.67, 9.39, 9.1, 8.8, 8.49, 8.17, 7.85, 7.53, 7.2,
    6.87, 6.54, 6.21, 5.87, 5.54, 5.21, 4.88, 4.56, 4.24, 3.92,
    3.61, 3.31, 3.02, 2.73, 2.46, 2.19, 1.94, 1.69, 1.46, 1.25,
    1.04, 0.857, 0.686, 0.532, 0.396, 0.279, 0.181, 0.103,
    0.0463, 0.0117
  };
  
  G4double norma = 0;
  for(G4int i=0;i<60;i++) norma += array[i];
  
  G4double cum = 0;
  for(G4int i=0;i<60;i++) {
    cum += array[i];
    G4double ene = (i*2.6 + 0.5)*keV;
    fShootEnergy.push_back(ene);    
    fShootRandom.push_back(cum/norma);
  }
}

//---------------------------------------------------------------------------//
void BxGeneratorSCS::Bi210CrossSection() {
  const G4double znucl = 83. ;
  const G4double ANUCL = 210.;
  const G4double AALPHA = 0.0072973531;
  const G4double ELMAS  = 0.5110034;
  const G4double zeta   = znucl+1 ;
  const G4double RNUCL  = 0.426* AALPHA*pow(ANUCL,0.333);
  const G4double QBETA  = 1.1621 ;
  const G4double GG0 = sqrt(1.-(AALPHA*zeta)*(AALPHA*zeta) );
  G4double norm = 0. ;
  for(G4int i=0; i<fSize; i++) {
    G4double T = G4double(i+0.5)/G4double(fSize)*QBETA;
    G4double WBETA = (T+ELMAS)/ELMAS;
    G4double YNU   = AALPHA*zeta*WBETA/sqrt(WBETA*WBETA -1.);
    G4double RSPC  = 2.*RNUCL*sqrt(WBETA*WBETA-1.);
    complex<G4double> ARGCOF1 = CGAMMA(complex<G4double>(GG0,YNU));
    complex<G4double> ARGCOF2 = CGAMMA(complex<G4double>(2*GG0+1.,0.));		
    G4double COFFERMI1 = pow(real(ARGCOF1),2)+pow(imag(ARGCOF1),2);
    G4double COFFERMI2 = pow(real(ARGCOF2),2)+pow(imag(ARGCOF2),2);
    G4double FFERMI = 2*(1+GG0)*pow(RSPC,2.*(GG0-1)) * exp(3.1415926*YNU) * COFFERMI1 / COFFERMI2 ;
    G4double PBETA= sqrt(T*T+2.*T*ELMAS);
    G4double FUNCNEW = FFERMI*pow(QBETA-T,2.)*(T+ELMAS)*PBETA;
    G4double BISHAPE  =1.+9.9*(T/ELMAS+1.)+383/((T/ELMAS)+1.)-8.4*pow(T/ELMAS+1.,2.); //from Aldo's note about Bi210
    FUNCNEW *= BISHAPE;
    norm += FUNCNEW*QBETA/G4double(fSize);
  }
  for(G4int k=0;k<fSize;k++) {
    G4double _fShootEnergy = 0;
    G4double _fShootRandom = 0;
    
    for(G4int i=0; i<k+1; i++) {
      G4double T = G4double(i+0.5)/G4double(fSize)*QBETA;
      G4double WBETA = (T+ELMAS)/ELMAS;
      G4double YNU   = AALPHA*zeta*WBETA/sqrt(WBETA*WBETA -1.);
      G4double RSPC  = 2.*RNUCL*sqrt(WBETA*WBETA-1.);
      complex<G4double> ARGCOF1 = CGAMMA(complex<G4double>(GG0,YNU));
      complex<G4double> ARGCOF2 = CGAMMA(complex<G4double>(2*GG0+1.,0.));		
      G4double COFFERMI1 = pow(real(ARGCOF1),2)+pow(imag(ARGCOF1),2);
      G4double COFFERMI2 = pow(real(ARGCOF2),2)+pow(imag(ARGCOF2),2);
      G4double FFERMI = 2*(1+GG0)*pow(RSPC,2.*(GG0-1)) * exp(3.1415926*YNU) * COFFERMI1 / COFFERMI2 ;
      G4double PBETA= sqrt(T*T+2.*T*ELMAS);
      G4double FUNCNEW = FFERMI*pow(QBETA-T,2.)*(T+ELMAS)*PBETA;
      G4double BISHAPE  =1.+9.9*(T/ELMAS+1.)+383/((T/ELMAS)+1.)-8.4*pow(T/ELMAS+1.,2.);
      FUNCNEW *= BISHAPE;
      _fShootEnergy = T*MeV ;
      _fShootRandom += FUNCNEW*QBETA/G4double(fSize)/norm;
    }     
    fShootEnergy.push_back(_fShootEnergy);
    fShootRandom.push_back(_fShootRandom);
    
  }  
  
}
//---------------------------------------------------------------------------//

G4double BxGeneratorSCS::ShootEnergy() {
  G4double val = 0;
  while(val==0) {
    G4double _val = G4UniformRand()  ;
    if(_val > fShootRandom[0]) val = _val;
  }
  for(G4int i=0;i<G4int(fShootRandom.size());i++) {
    if(fShootRandom[i] >= val) {
      if(i == 0) return fShootEnergy[0];
      //G4double deltaX = val - fShootRandom[i] ;
      G4double y = fShootEnergy[i] - fShootEnergy[i-1] ;
      G4double x = fShootRandom[i] - fShootRandom[i-1] ;
      G4double q = fShootEnergy[i] - y/x*fShootRandom[i];
      
      return y/x*val+q;
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------//
complex<G4double> BxGeneratorSCS::CGAMMA( complex<G4double> z) {
  const G4int par = 7 ;
  if ( real(z)<0.5 ) return M_PI / (sin(M_PI*z)*CGAMMA(1.0-z));
  z -= 1.0;
  complex<G4double> x = ffp[0];
  for (G4int i=1; i<par+2; i++) x += ffp[i]/(z+complex<G4double>(i,0));
  complex<G4double> t = z + (par + 0.5);
  return sqrt(2.*M_PI) * pow(t,z+0.5) * exp(-t) * x;
}
