// --------------------------------------------------------------------------//
/** 
 * AUTHOR: I.Machulin
 * CONTACT: machulin@lngs.infn.it
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
 *
 * The model for Li-He beta-neutron decay is based mainly on the
 * Nuclear Data Sheets and recent article
 * "Low-lying resonance states in the Be9 continium" Ph.Let. B 618 2005 43-50
 * Parameters for 8He decay are taken from 
 * Nuclear Physics A 366 (1981) 461-468
 * Nuclear Physics A 487 (1988) 269-278
 * Parameters for the 9Li decay taken from the internal note "Detection of the cosmogenic radioactive isotopes 8He and 9Li via the (beta-n) decay channel in Borexino" by V. Kobychev and Y. Suvorov (December 2, 2009)
 *A small decay channel 8He->6He+d->5He+t is neglected. The branching ratio is small (about 0.9%).
*/
// --------------------------------------------------------------------------//
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
#include "BxGeneratorLi9He8.hh"
#include "G4SystemOfUnits.hh"
//---------------------------------------------------------------------------//


BxGeneratorLi9He8::BxGeneratorLi9He8(): BxVGenerator("BxGeneratorLi9He8") {


  isFirstTime = true ;
  fNeutrinoType = -1;
  


  fVolumeFlag = false ;
  G4ThreeVector zero(0., 0., 0.) ;
  fSPSAng = new G4SPSAngDistribution ;
  fSPSPos = new G4SPSPosDistribution ;
  G4SPSRandomGenerator *RndGen = new G4SPSRandomGenerator;
  fSPSPos->SetBiasRndm(RndGen);
  fSPSAng->SetBiasRndm(RndGen);
  SetRadius0(0.);
  SetCentreCoords(zero);    
  fPosition = G4ThreeVector(0.,0.,0.) ;
  fParticleTable = G4ParticleTable::GetParticleTable();    

  fTheMessenger = new BxGeneratorLi9He8Messenger(this);
 
}
//---------------------------------------------------------------------------//


BxGeneratorLi9He8::~BxGeneratorLi9He8()
{
  delete fTheMessenger;
  delete fSPSPos;
  delete fSPSAng;
  
}

//---------------------------------------------------------------------------//

void BxGeneratorLi9He8::BxGeneratePrimaries(G4Event *event) {

  if(isFirstTime) {
    if(fNeutrinoType < 0) {
      BxLog(error) << "Set the antineutrino source type"<<endlog;
      BxLog(fatal) << endlog;
    }  

    BxLog(routine) << "Li9He8 source " << fNeutrinoType << endlog ;

    isFirstTime = false ;
  }

   if ((event->GetEventID())%2 == 0){

  if(fVolumeFlag)          
  {
  fSPSAng->SetAngDistType("iso");   
  fSPSAng->SetPosDistribution(fSPSPos); 
  fPosition = fSPSPos->GenerateOne();
  }  

    if ((event->GetEventID())%2 == 0)
    {         
             
    fParticle = fParticleTable->FindParticle(11);     // 11 for the electron
    
    
    G4double KinE=0, E0=0;
    G4double KinEalpha1=0,  KinEalpha2=0, KinEneutron=0;        

    G4int FlagIso=0;
    

 //  Here i generate the beta-particles

   G4int ChHe = 0, ChLi = 0, ChLiBeta=0, ChHeBeta = 0;
   
   G4double   HeW =0, LiW=0;
   
   if (fNeutrinoType == 0) // He production
    
    {

  ChHeBeta = ChooseBranchHe(); 
    
   if (ChHeBeta == 1) 
   {
   ChHe = 1;
   E0 = 7.44;   

   G4int D=0;

   while (D==0)
{
   HeW = G4RandGauss::shoot(0,1.0/2.35) ;  // The width  of the level is 1.0  MeV  

   if ( (HeW > -0.7) && (HeW < + 0.7 ) ) D=1;
}

   E0= E0 + HeW; 

//  //G4cout << "  E0 7.44  = " << E0 << G4endl;
   }
   
   else 
   {   
   E0 = 5.25; 
   ChHe=2;      

   G4int D=0;

   while (D==0)
{
   HeW = G4RandGauss::shoot(0 ,0.65/2.35) ;  // The width  of the level is 0.65  MeV  

   if ( (HeW > -2.0) && (HeW <  2.0))  D=1;

}

   E0= E0 + HeW; 

}

    KinE = (ShootEnergyElectronHe(E0)) /MeV;
 
//   //G4cout << "  KinE = " << KinE << G4endl;
//    //G4cout << "  Ch = " << ChHe << G4endl;
//    //G4cout << "  E0 = " << E0 << G4endl;    
    
    FlagIso=2;
    
    }
    

    
   if (fNeutrinoType == 1)   //Li production, total delta - 13.61 MeV

// Branch-E_level         31.9% - 2.43 Mev; 11.6% - 2.78 MeV; 3.2 % - 5.0 MeV; 1.5 % - 7.94 MeV ; 2.7 % - 11.81 MeV;   Total 50.9 %
// E_beta                 	 11.18 MeV;        10.83 MeV;         8.61 MeV;        5.67 MeV;           1.80 MeV;
// Width                  	 0.77  keV;         1.08 MeV;         2.0 MeV;         1.0  MeV;           0.4  MeV
// ChLiBeta			  1		    2                 3                4                   5
    {
    
    FlagIso=3;
       
   ChLiBeta = ChooseBranchLi(); 

   if (ChLiBeta == 1) 
   {
   E0 = 11.18;
   ChLi = 1; //Low energy neutron emission
   ChLiBeta =1;
 // The width  of the level is 0.77 keV  -no need to widen
}

         else  if (ChLiBeta == 2) 
   {
   E0 = 10.83;
   ChLi = 1; 
   G4int D=0;
   ChLiBeta =2;
 // The width  of the level is 1080  keV
    while (D==0)
{
   LiW = G4RandGauss::shoot(0,1.08/2.35) ;  // The width  of the level is 1.08  MeV  
   if ( (LiW > -1.08) && (LiW <  1.08))  D=1;
}
   E0= E0 + LiW; 
   }

   else  if (ChLiBeta == 3) 
   {
   E0 = 8.6;
   ChLi = 1;
   ChLiBeta = 3;
   G4int D=0;
 // The width  of the level is 2  MeV
    while (D==0)
{
   LiW = G4RandGauss::shoot(0,2.0/2.35) ;  // The width  of the level is 2.0  MeV  
   if ( (LiW > -2.0) && (LiW <  2.0))  D=1;
}
   E0= E0 + LiW; 
   }

   else  if (ChLiBeta == 4 ) 
   {
   E0 = 5.66;
   ChLi = 2;
   ChLiBeta = 4;
   G4int D=0;
 // The width  of the level is 1.0  MeV
    while (D==0)
{
   LiW = G4RandGauss::shoot(0,1.0/2.35) ;  // The width  of the level is 1.0  MeV  
   if ( (LiW > -1.0) && (LiW <  1.0))  D=1;
}
   E0= E0 + LiW; 
   }

   
   else  
   {
   E0 = 1.8;       
   ChLi = 2;
   ChLiBeta = 5;
   G4int D=0;
 // The width  of the level is 400  keV
    while (D==0)
{
   LiW = G4RandGauss::shoot(0,0.400/2.35) ;  // The width  of the level is 0.4  MeV  
   if ( (LiW > -0.4) && (LiW <  0.4))  D=1;
}
   E0= E0 + LiW; 
   }
    
    KinE = ShootEnergyElectronLi(E0)/MeV;  //This is electron Energy
    
    }
    
         
    G4double energy = KinE + fParticle->GetPDGMass();	
    
    fDirection = fSPSAng->GenerateOne();
    
    G4double pmom = std::sqrt(pow(energy,2.) - pow(fParticle->GetPDGMass(),2.));
    G4double px = pmom*fDirection.x();
    G4double py = pmom*fDirection.y();
    G4double pz = pmom*fDirection.z();
    G4PrimaryVertex*   vertex   = new G4PrimaryVertex(fPosition,0);
    G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle,px,py,pz);
    vertex->SetPrimary( particle );
    event->AddPrimaryVertex( vertex );    

    BxOutputVertex::Get()->SetIsotopeCoinc(0);
    
    
    BxOutputVertex::Get()->SetEnergy(KinE);		
    BxOutputVertex::Get()->SetPosition(fPosition);		
    BxOutputVertex::Get()->SetDirection(fDirection);
    BxOutputVertex::Get()->SetTime(0);
//    BxOutputVertex::Get()->SetPDG(fParticle->GetPDGEncoding());
    
    
//    BxOutputVertex::Get()->SetUserInt1 (ChLiBeta);   
    
    BxOutputVertex::Get()->SetPDG(ChLiBeta);    
    
         
// Here i generate the gamma 478 keV for He-8 -> Li7* (Li7 + 478 keV) + n    
// supposed to be 33 % case for neutron decay    
    
   if (FlagIso == 2) 
    
    {
	    G4double    branch = G4UniformRand();
	    if(branch<0.33){
		    KinE =  0.478 /MeV;

		    fParticle = fParticleTable->FindParticle(22);     

		    energy = KinE + fParticle->GetPDGMass();

		    fDirection = fSPSAng->GenerateOne();

		    pmom = std::sqrt(pow(energy,2.) - pow(fParticle->GetPDGMass(),2.));
		    px = pmom*fDirection.x();
		    py = pmom*fDirection.y();
		    pz = pmom*fDirection.z();

		    vertex   = new G4PrimaryVertex(fPosition,0);
		    particle = new G4PrimaryParticle(fParticle,px,py,pz);

		    vertex->SetPrimary( particle );
		    event->AddPrimaryVertex( vertex );    

		    BxOutputVertex::Get()->SetDId(1);     
		    BxOutputVertex::Get()->SetDPosition(fPosition);  
		    //    BxOutputVertex::Get()->SetEnergy(KinE);	
		    //   BxOutputVertex::Get()->SetDDirection(pmom); 
		    BxOutputVertex::Get()->SetDPDG(22);
		    BxOutputVertex::Get()->SetDTime(0.);
		    BxOutputVertex::Get()->SetDEnergy(KinE);		
		    BxOutputVertex::Get()->SetDaughters();        
	    }
    }    

// the generation of the neutrons energy from Li decay    

   KinEneutron = ShootEnergyNeutronLi(ChLi)/MeV; 

// the generation of the alpha particles from Li decay.       
   if (FlagIso == 3) 
    
    {
    KinEalpha1 = ShootEnergyAlpha1Li(ChLiBeta)/MeV;

      //G4cout << "KinEalpha =" << KinEalpha1 << G4endl;
      //G4cout << "EnergyNeutron= " << KinEneutron << G4endl;

    KinE =  KinEalpha1;
    
     fParticle = fParticleTable->FindParticle(1000020040);     
    
     energy = KinE + fParticle->GetPDGMass();

     fDirection = fSPSAng->GenerateOne();
    
     pmom = std::sqrt(pow(energy,2.) - pow(fParticle->GetPDGMass(),2.));
     px = pmom*fDirection.x();
     py = pmom*fDirection.y();
     pz = pmom*fDirection.z();

     vertex   = new G4PrimaryVertex(fPosition,0);
     particle = new G4PrimaryParticle(fParticle,px,py,pz);

    vertex->SetPrimary( particle );
    event->AddPrimaryVertex( vertex );    
 
    BxOutputVertex::Get()->SetDId(2);     
    BxOutputVertex::Get()->SetDPosition(fPosition);  
//    BxOutputVertex::Get()->SetEnergy(KinE);	
//    BxOutputVertex::Get()->SetDDirection(pmom); 
    BxOutputVertex::Get()->SetDPDG(1000020040);
    BxOutputVertex::Get()->SetDTime(0.);
    BxOutputVertex::Get()->SetDEnergy(KinE);		
    BxOutputVertex::Get()->SetDaughters();     
    
    G4double Sum_alpha_n = KinEalpha1 + KinEneutron;

    KinEalpha2 = ShootEnergyAlpha2Li(ChLiBeta, Sum_alpha_n)/MeV;

fParticle = fParticleTable->FindParticle(1000020040);

     energy = KinEalpha2 + fParticle->GetPDGMass();

     fDirection = fSPSAng->GenerateOne();
//the momentum direction is randomly chosen. The correlation between the two gamma directions is not taken into account.
     pmom = std::sqrt(pow(energy,2.) - pow(fParticle->GetPDGMass(),2.));
     px = pmom*fDirection.x();
     py = pmom*fDirection.y();
     pz = pmom*fDirection.z();

     vertex   = new G4PrimaryVertex(fPosition,0);
     particle = new G4PrimaryParticle(fParticle,px,py,pz);

    vertex->SetPrimary( particle );
    event->AddPrimaryVertex( vertex );
    BxOutputVertex::Get()->SetDId(3);
    BxOutputVertex::Get()->SetDPosition(fPosition);
    //    BxOutputVertex::Get()->SetEnergy(KinEalpha2);   
    //    BxOutputVertex::Get()->SetDDirection(pmom); 
    BxOutputVertex::Get()->SetDPDG(1000020040);
    BxOutputVertex::Get()->SetDTime(0.);
    BxOutputVertex::Get()->SetDEnergy(KinEalpha2);
    BxOutputVertex::Get()->SetDaughters();

    }



// the generation of the neutrons from Li and He decay

    fParticle = fParticleTable->FindParticle(2112);     

     if (FlagIso == 2)
  {
     KinE = ShootEnergyNeutronHe(ChHe)/MeV;
}

     else KinE = KinEneutron;     
	  	  
     energy = KinE + fParticle->GetPDGMass();

     fDirection = fSPSAng->GenerateOne();
    
     pmom = std::sqrt(pow(energy,2.) - pow(fParticle->GetPDGMass(),2.));
     px = pmom*fDirection.x();
     py = pmom*fDirection.y();
     pz = pmom*fDirection.z();

     vertex   = new G4PrimaryVertex(fPosition,0);
     particle = new G4PrimaryParticle(fParticle,px,py,pz);

    vertex->SetPrimary( particle );
    event->AddPrimaryVertex( vertex );    
    BxOutputVertex::Get()->SetDId(0);     
    BxOutputVertex::Get()->SetDPosition(fPosition);  
    BxOutputVertex::Get()->SetEnergy(KinEalpha1+KinE+KinEalpha2);
 //   BxOutputVertex::Get()->SetDDirection(pmom); 
    BxOutputVertex::Get()->SetDPDG(2112);
    BxOutputVertex::Get()->SetDTime(0.);
    BxOutputVertex::Get()->SetDEnergy(KinE);		
    BxOutputVertex::Get()->SetDaughters(); 
     
}
           	
}

}

//-------------------------------------------------------------------------
//   Branch for beta decay for Li spectrum
//-------------------------------------------------------------------------

G4int BxGeneratorLi9He8::ChooseBranchLi() {

G4double    Branch = G4UniformRand(); 
G4int    Ch =0;

   if (Branch < (31.9/50.9)) { Ch = 1; }      
   else  if (Branch < ((31.9+11.6)/50.9)) { Ch = 2;}
   else  if (Branch < ((31.9+11.6+3.2)/50.9))  { Ch = 3;}
   else  if (Branch <((31.9+11.6+3.2+1.5)/50.9) ) { Ch = 4;} 
   else   { Ch = 5;}
   
      return Ch;
}

//-------------------------------------------------------------------------
//   Branch for beta decay for He spectrum
//-------------------------------------------------------------------------
//The branching ratio of the channel 8He->8Li*(3.21) and 8He->8Li*(5.4)
//is not well known in literature. The brookhaven web site 
//(http://www.nndc.bnl.gov/) reports that the two branching ratio are similar.
//Other authors (Nuclear Physics A 487 (1988) 269-278) report P(8Li*(3.21))/P(8Li*(3.21))~3
G4int BxGeneratorLi9He8::ChooseBranchHe() {

G4double    Branch = G4UniformRand(); 
G4int   Ch =0;
   if (Branch < 0.5) { Ch = 1; }      
   else   { Ch = 2;}   
      return Ch;
}

//-------------------------------------------------------------------------
//    Energy Shooter for first alpha for Li spectrum
//-------------------------------------------------------------------------
G4double BxGeneratorLi9He8::ShootEnergyAlpha1Li(G4int Ch) {

   G4double AlphaEnergySp1 [11] =    
 {    0.,
      0.,
      2.,
      40.,
      90.,
      60.,
      15.,
      2.,
      0.,
      0.,
      0.      }; // End point for alpha = 1 MeV (0.7 MeV in exp data) bin = 100 keV -2.43 MeV level

   G4double AlphaEnergySp2_3 [21] =    
 {    0.,
      40.,
      100.,
      130.,
      180.,
      130.,
      110.,
      90.,
      80.,
      60.,
      50.,      
      30.,
      20.,
      10.,      
      7.,
      2.,
      1.,
      0.,      
      0.,      
      0.,
      0.         }; // End point for alpha = 4 MeV (3.2 in exp data) bin = 200 keV - 2.78 MeV level (and by hand - 5.0 MeV Level)

   G4double AlphaEnergySp4 [22] =    
 {    0.,
      1.5,
      3.,
      8.,
      11.,
      15.,
      25.,
      28.,
      31.,
      32.,
      35.,      
      39.,
      40.,
      50.,      
      50.,
      48.,
      25.,
      22.,      
      10.,      
      3.,
      2.,
      0.8	}; // End point for alpha = 4.2 MeV (4.2 in exp data) bin = 200 keV -  7.94 MeV level

   G4double AlphaEnergySp5 [31] =    
 {    0.,
      8.,
      40.,
      50.,
      55.,
      60.,
      70.,
      90.,
      100.,
      100.,
      90.,     
      80.,
      80.,
      75.,      
      70.,
      70.,
      65.,
      65.,      
      60.,      
      50.,
      39.,
      40.,
      35.,
      33.,
      24.,
      20.,
      18.,
      12.,
      10.,
      0.,
      0.      }; // End point for alpha = 6.0 MeV (5.8 in exp data) bin = 200 keV - 11.81 MeV level

   G4int N = 0;
   G4double sum = 0.;
   G4double norma = 0.;
   G4double xE;
   G4double deltaX, x, y;

if (Ch == 1)
{
N=11;
   G4double Probability [11];
   G4double EnergyBin [11];
   G4double E0 = 1.0; //end Total Energy of the spectre

   G4double scale = E0/N;

     for(G4int i = 0; i< N ; i++) {
      norma += AlphaEnergySp1[i];
     }   
        
   for (G4int i = 0; i < N; i++) {
     xE = scale*G4float(i);
     EnergyBin[i] = xE;     
}

   for(G4int i = 0; i < N; i++) {
     sum += AlphaEnergySp1[i]/norma;
     Probability[i] = sum;
   }

   G4double val = G4UniformRand();

   for(G4int i = 0; i < N; i++) {
     if(Probability[i] >= val) {
       if(i == 0) return EnergyBin[0];
       deltaX = val - Probability[i] ;
       y = EnergyBin[i] - EnergyBin[i-1] ;
       x = Probability[i] - Probability[i-1] ;
       return (deltaX*y/x + EnergyBin[i]) ;
     }
   }
}
if (Ch == 2)
{
N=21;
   G4double Probability [21];
   G4double EnergyBin [21];
   G4double E0 = 4.0; //end Total Energy of the spectre

   G4double scale = E0/N;

     for(G4int i = 0; i< N ; i++) {
      norma += AlphaEnergySp2_3[i];
     }   
        
   for (G4int i = 0; i < N; i++) {
     xE = scale*G4float(i);
     EnergyBin[i] = xE;     
}

   for(G4int i = 0; i < N; i++) {
     sum += AlphaEnergySp2_3[i]/norma;
     Probability[i] = sum;
   }

   G4double val = G4UniformRand();

   for(G4int i = 0; i < N; i++) {
     if(Probability[i] >= val) {
       if(i == 0) return EnergyBin[0];
       deltaX = val - Probability[i] ;
       y = EnergyBin[i] - EnergyBin[i-1] ;
       x = Probability[i] - Probability[i-1] ;
       return (deltaX*y/x + EnergyBin[i]) ;
     }
   }
}
if (Ch == 3)
{
N=21;
   G4double Probability [21];
   G4double EnergyBin [21];
   G4double E0 = 4.0; //end Total Energy of the spectre

   G4double scale = E0/(N-1);

     for(G4int i = 0; i< N ; i++) {
      norma += AlphaEnergySp2_3[i];
     }   
        
   for (G4int i = 0; i < N; i++) {
     xE = scale*G4float(i);
     EnergyBin[i] = xE;     
}

   for(G4int i = 0; i < N; i++) {
     sum += AlphaEnergySp2_3[i]/norma;
     Probability[i] = sum;
   }

   G4double val = G4UniformRand();

   for(G4int i = 0; i < N; i++) {
     if(Probability[i] >= val) {
       if(i == 0) return EnergyBin[0];
       deltaX = val - Probability[i] ;
       y = EnergyBin[i] - EnergyBin[i-1] ;
       x = Probability[i] - Probability[i-1] ;
       return (deltaX*y/x + EnergyBin[i]) ;
     }
   }
}
if (Ch == 4)
{
N=22;
   G4double Probability [22];
   G4double EnergyBin [22];
   G4double E0 = 4.2; //end Total Energy of the spectre

   G4double scale = E0/N;

     for(G4int i = 0; i< N ; i++) {
      norma += AlphaEnergySp4[i];
     }   
        
   for (G4int i = 0; i < N; i++) {
     xE = scale*G4float(i);
     EnergyBin[i] = xE;     
}

   for(G4int i = 0; i < N; i++) {
     sum += AlphaEnergySp4[i]/norma;
     Probability[i] = sum;
   }

   G4double val = G4UniformRand();

   for(G4int i = 0; i < N; i++) {
     if(Probability[i] >= val) {
       if(i == 0) return EnergyBin[0];
       deltaX = val - Probability[i] ;
       y = EnergyBin[i] - EnergyBin[i-1] ;
       x = Probability[i] - Probability[i-1] ;
       return (deltaX*y/x + EnergyBin[i]) ;
     }
   }
}
if (Ch == 5)
{
N=31;
   G4double Probability [31];
   G4double EnergyBin [31];
   G4double E0 = 6.0; //end Total Energy of the spectre

   G4double scale = E0/N;

     for(G4int i = 0; i< N ; i++) {
      norma += AlphaEnergySp5[i];
     }   
        
   for (G4int i = 0; i < N; i++) {
     xE = scale*G4float(i);
     EnergyBin[i] = xE;     
}

   for(G4int i = 0; i < N; i++) {
     sum += AlphaEnergySp5[i]/norma;
     Probability[i] = sum;
   }

   G4double val = G4UniformRand();

   for(G4int i = 0; i < N; i++) {
     if(Probability[i] >= val) {
       if(i == 0) return EnergyBin[0];
       deltaX = val - Probability[i] ;
       y = EnergyBin[i] - EnergyBin[i-1] ;
       x = Probability[i] - Probability[i-1] ;
       return (deltaX*y/x + EnergyBin[i]) ;
     }
   }
}
return 0;

}

//-------------------------------------------------------------------------
//    Energy Shooter for second alpha for Li spectrum
//-------------------------------------------------------------------------
G4double BxGeneratorLi9He8::ShootEnergyAlpha2Li(G4int Ch, G4double E2Sum) {
G4double ESum;
G4double EAlpha2;

//G4double AlphaQ=1.57;

if (Ch == 1)
{
ESum=2.43-1.57;
}
if (Ch == 2)
{
ESum=2.78-1.57;
}
if (Ch == 3)
{
ESum=5.0-1.57;
}
if (Ch == 4)
{
ESum=7.94-1.57;
}
if (Ch == 5)
{
ESum=11.81-1.57;
}


EAlpha2= ESum-E2Sum;

if (EAlpha2 < 0) EAlpha2=0;

return EAlpha2;

}

//-------------------------------------------------------------------------
//    Energy Shooter for neutrons for He spectrum
//-------------------------------------------------------------------------

G4double BxGeneratorLi9He8::ShootEnergyNeutronHe(G4int Ch) {
       G4double EnergyBin; 
      if (Ch == 1 ) EnergyBin=0.7;
      else EnergyBin=2.9;
      return EnergyBin;
}

//-------------------------------------------------------------------------
//    Energy Shooter for neutrons for Li spectrum
//-------------------------------------------------------------------------

G4double BxGeneratorLi9He8::ShootEnergyNeutronLi(G4int Ch) {

//This is the high energy neutron spectrum - which should correspond to the beta decay to high energy levels
//from 0 up to 3 MeV

   G4double NeutronEnergySpHigh [16] =    
 { 
      0.0,
      2.24,
      5.0,
      7.94,
      11.2,
      12.6,
      11.0,
      10.0,
      7.9,
      5.0,
      4.0,
      3.2,
      2.0,
      0.5,
      0.1,
      0.0
};

//This is the low energy neutron spectrum - which should correspond to the beta decay on low energy levels
//From 0 up to 1 MeV

   G4double NeutronEnergySpLow [21] = 
{
      0.0,
      2100.0,
      5447.74,
      6500,
      7306.69,
      7000.,
      6625.06,
      5000.,
      3322.82,
      2200.,  
      1149.,
      900.,
      911.57,
      2100,
      955.64,
      100.0,
      0.36,
      0.08,
      0.08,
      0.0,
      0.0,
  };

   G4double ProbabilityHigh[16];
   G4double ProbabilityLow[21];   
   G4double EnergyBinHigh[16];  
   G4double EnergyBinLow[21];    

   G4double sum = 0;
   G4double norma =0; 
   G4double deltaX,x,y;
   G4int Num;
   G4double scale;

if (Ch==2) { Num =16; 
	     scale=200;

   norma= 0;

     for(G4int i = 0; i< Num; i++) {
      norma += NeutronEnergySpHigh[i];
     }

    for(G4int i = 0; i< Num ; i++) {
      sum += NeutronEnergySpHigh[i]/norma;
      ProbabilityHigh[i]=sum;
      EnergyBinHigh[i]=G4float(i)*scale;
    }

  G4double val = G4UniformRand()  ;

  for(G4int i=0;i<Num;i++) {
    if(ProbabilityHigh[i] >= val) {

      if(i == 0) {  
      return EnergyBinHigh[0]/1000;
}

      deltaX = val - ProbabilityHigh[i] ;
      y = EnergyBinHigh[i] - EnergyBinHigh[i-1] ;
      x = ProbabilityHigh[i] - ProbabilityHigh[i-1] ;
      
       
      return ((deltaX*y/x + EnergyBinHigh[i])/1000) ;

    }
  }


}

if (Ch==1) { Num =21; 
	     scale=50;
   norma= 0;
   
     for(G4int i = 0; i< Num; i++) {
      norma += NeutronEnergySpLow[i];
     }

    for(G4int i = 0; i< Num ; i++) {
      sum += NeutronEnergySpLow[i]/norma;
      ProbabilityLow[i]=sum;
      EnergyBinLow[i]=G4float(i)*scale;
    }

  G4double val = G4UniformRand()  ;

  for(G4int i=0;i<Num;i++) {
    if(ProbabilityLow[i] >= val) {

      if(i == 0) {  
      return EnergyBinLow[0]/1000;
}

      deltaX = val - ProbabilityLow[i] ;
      y = EnergyBinLow[i] - EnergyBinLow[i-1] ;
      x = ProbabilityLow[i] - ProbabilityLow[i-1] ;
      
      return ((deltaX*y/x + EnergyBinLow[i])/1000) ;

    }
  }


}   



  return 0;
}

//-------------------------------------------------------------------------
//    Energy Shooter for electrons from Li spectrum
//-------------------------------------------------------------------------

G4double BxGeneratorLi9He8::ShootEnergyElectronLi(G4double lll) {
   
   G4double E0 = lll; //  end-point kinetic energy


   G4double sum = 0.;
   G4double norma = 0.;
   G4double Probability [1000];
   G4double EnergyBin [1000];
              
   G4double Z0 = 4.;                  // 7. - electron (N), 5. = positron (B): doughter nucleus charge 
   G4double xE;                       // beta kinetic energy 
   G4double me = 0.511;
   G4double kk = 4.586E-2, FF = 0., xx = 0.;   // Fermi function, kk=1/137*2*pi
   G4double a1 = -1./2., a2 = -1./24.;         // "-" - electron, "+" - positron
   G4double scale = E0/1000.;
   G4double deltaX, x, y;
   G4double PositronEnergySp[1000]; 



   for (G4int i = 0; i < 1000; i++) {

     xE = scale*(G4double(i) + 0.5);
     EnergyBin[i] = xE;
     PositronEnergySp[i] = 0.;

     xx = kk*Z0*(xE+me)/(sqrt(xE*xE + 2.*me*xE));
     FF = xx/(xx + a1*xx*xx + xx*xx*xx/6. + a2*xx*xx*xx*xx);//Fermi function calculation, the exponential is expanded in Taylor's series 
     if (xE < E0) 
         PositronEnergySp[i] = FF*(E0 - xE)*(E0 - xE)*(xE + me)*sqrt(xE*xE + 2.*me*xE); 
	
     norma += PositronEnergySp[i];	
   }

   for(G4int i = 0; i < 1000; i++) {
     sum += PositronEnergySp[i]/norma;
     Probability[i] = sum;
   }

   G4double val = G4UniformRand();

   for(G4int i = 0; i < 1000; i++) {
     if(Probability[i] >= val) {
       if(i == 0) return EnergyBin[0];
       deltaX = val - Probability[i] ;
       y = EnergyBin[i] - EnergyBin[i-1] ;
       x = Probability[i] - Probability[i-1] ;
       return (deltaX*y/x + EnergyBin[i]) ;
     }
   }
   
   return 0;

}

//-------------------------------------------------------------------------
//    Energy Shooter for electrons from He spectrum
//-------------------------------------------------------------------------

G4double BxGeneratorLi9He8::ShootEnergyElectronHe(G4double Egr) {


   G4double E0 = Egr;  //  end-point kinetic energy

   G4double sum = 0.;
   G4double norma = 0.;
   G4double Probability [1000];
   G4double EnergyBin [1000];           
   G4double Z0 = 3.;                  // 7. - electron (N), 5. = positron (B): doughter nucleus charge 
   G4double xE;                       // beta kinetic energy 
   G4double me = 0.511;
   G4double kk = 4.586E-2, FF = 0., xx = 0.;   // Fermi function, kk=1/137*2*pi
   G4double a1 = -1./2., a2 = -1./24.;         // "-" - electron, "+" - positron
   G4double scale = E0/1000.;
   G4double deltaX, x, y;
   G4double PositronEnergySp[1000]; 



   for (G4int i = 0; i < 1000; i++) {

     xE = scale*(G4double(i) + 0.5);
     EnergyBin[i] = xE;
     PositronEnergySp[i] = 0.;

     xx = kk*Z0*(xE+me)/(sqrt(xE*xE + 2.*me*xE));
     FF = xx/(xx + a1*xx*xx + xx*xx*xx/6. + a2*xx*xx*xx*xx);//Fermi function calculation, the exponential is expanded in Taylor's series  
     if (xE < E0) 

     PositronEnergySp[i] = FF*(E0 - xE)*(E0 - xE)*(xE + me)*sqrt(xE*xE + 2.*me*xE); 
	
     norma += PositronEnergySp[i];	

   }

   for(G4int i = 0; i < 1000; i++) {

     sum += PositronEnergySp[i]/norma;
     Probability[i] = sum;

   }

   G4double val = G4UniformRand();

   for(G4int i = 0; i < 1000; i++) {
   
     if(Probability[i] >= val) {
       if(i == 0) return EnergyBin[0];
       deltaX = val - Probability[i] ;
       y = EnergyBin[i] - EnergyBin[i-1] ;
       x = Probability[i] - Probability[i-1] ;
       return (deltaX*y/x + EnergyBin[i]) ;
     }

}
   return 0;

}



