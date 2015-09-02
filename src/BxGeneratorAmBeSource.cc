// --------------------------------------------------------------------------//
/** 
 * AUTHOR: I.Machulin
 * CONTACT: machulin@lngs.infn.it
* Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

//Useful note: Calibration of the Borexino using the Am-Be source
//E. Litvinovich, I. Machulin, S. Silaeva, S. Sukhotin, Yu. Suvorov

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
#include "BxGeneratorAmBeSource.hh"
#include "G4SystemOfUnits.hh"
#include "BxGeneratorAmBeSourceMessenger.hh"
//---------------------------------------------------------------------------//

BxGeneratorAmBeSource::BxGeneratorAmBeSource(): BxVGenerator("BxGeneratorAmBeSource") {


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

  fTheMessenger = new BxGeneratorAmBeSourceMessenger(this);
}
//---------------------------------------------------------------------------//


BxGeneratorAmBeSource::~BxGeneratorAmBeSource()
{
  delete fTheMessenger;
  delete fSPSPos;
  delete fSPSAng;
  
}

//---------------------------------------------------------------------------//


void BxGeneratorAmBeSource::BxGeneratePrimaries(G4Event *event) {

	// We escape the second event by choosing the false e- generation to process the capture gamma from stack
	if ((event->GetEventID())%2 == 0) {

		if(isFirstTime) {
			if(fNeutrinoType < 0) {
				BxLog(error) << "Set the AmBeSource source type"<<endlog;
				BxLog(fatal) << endlog;
			}  

			BxLog(routine) << "AmBeSource  " << fNeutrinoType << endlog ;

			isFirstTime = false ;
		}

		if((event->GetEventID())%2 ==0){
			if(fVolumeFlag)          
			{
				fSPSAng->SetAngDistType("iso");   
				fSPSAng->SetPosDistribution(fSPSPos); 
				fPosition = fSPSPos->GenerateOne();
			}  

			// We choose one of the channels to simulate 0 Gamma, 1 Gamma, 2 Gammas 

			G4int    Scheme=0;
			//0 -NoGamma
			//1 -1Gamma
			//2 -2Gamma

//			G4double ScRnd = G4UniformRand();

//			if (ScRnd < 0.36)  Scheme=0;

//			else if(ScRnd < 0.97) Scheme=1;

//			else Scheme=2;


			//Here I choose different branches of the alpha + Be-9 reactions
			//full simulation 0
			//neutrons without Gamma (group n0) 
			//neutrons with 1 Gamma  (group n1)
			//neutrons with 0 Gammas (group n2)
			//see articles reported below

//			if (fNeutrinoType == 1) Scheme=0;
//			if (fNeutrinoType == 2) Scheme=1;
//			if (fNeutrinoType == 3) Scheme=2;

                        Scheme=0;

  G4ThreeVector RandPos;
  G4double z=0.;
  G4double r=0,phi=0, theta=0;
  G4double rxy = 0;
  G4double twopi = 2*3.1414;

      r     = 9200*mm;                                   //!!!
      theta= twopi/4 * G4UniformRand();
      phi   = twopi * G4UniformRand();

      rxy = r* sin (theta);    
      
      
      RandPos.setX(rxy*cos(phi));
      RandPos.setY(rxy*sin(phi));
      RandPos.setZ(1037+r*cos(theta));

      fPosition = RandPos;



//Angular distribution the same as for muon
  
      //According to Phys. Rev. D 44 (1991) 3543, the angular spectrum underground 
      //is proportional to 1/cos(theta) --> with respect to the horizontal
      G4double costhetamin = 0.00001;
      G4double costheta = pow(costhetamin,(1.0-G4UniformRand()));
      phi = G4UniformRand()*360*degree;
      G4double sintheta = sqrt(1.0-costheta*costheta);
      fDirection.setX(cos(phi)*costheta);
      fDirection.setY(sin(phi)*costheta);
      fDirection.setZ(-1.0*sintheta);



			//Here we generate the Energy of the neutron
			G4double KinE = 0.0;

			if (Scheme==0){ KinE = ShootEnergyNeutron0G()/MeV;}

			if (Scheme==1){ KinE = ShootEnergyNeutron1G()/MeV;}

			if (Scheme==2){ KinE = ShootEnergyNeutron2G()/MeV;}

			fParticle = fParticleTable->FindParticle(2112);     

			G4double mass = fParticle->GetPDGMass();
			G4double energy = KinE + mass;	
G4cout << "KinE = "<< KinE << G4endl;
//			fDirection = fSPSAng->GenerateOne();



			G4double pmom = std::sqrt(energy*energy - mass*mass);
			G4double px = pmom*fDirection.x();
			G4double py = pmom*fDirection.y();
			G4double pz = pmom*fDirection.z();

			G4PrimaryVertex*   vertex   = new G4PrimaryVertex(fPosition,0);
			G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle,px,py,pz);
			vertex->SetPrimary( particle );
			event->AddPrimaryVertex( vertex );    

			BxOutputVertex::Get()->SetEnergy(KinE);		
			BxOutputVertex::Get()->SetPosition(fPosition);		
//			BxOutputVertex::Get()->SetDirection(pmom); 
			BxOutputVertex::Get()->SetPDG(2112);
			BxOutputVertex::Get()->SetTime(0.);


			if (Scheme==1){ 

				fParticle = fParticleTable->FindParticle(22);     

				G4double pmass = fParticle->GetPDGMass();
				G4double penergy = 4.439/MeV;

				fDirection = fSPSAng->GenerateOne();

				G4double partmom = std::sqrt(penergy*penergy - pmass*pmass);
				G4double p_x = partmom*fDirection.x();
				G4double p_y = partmom*fDirection.y();
				G4double p_z = partmom*fDirection.z();

				G4PrimaryVertex*   newvertex   = new G4PrimaryVertex(fPosition,0);
				G4PrimaryParticle* newparticle = new G4PrimaryParticle(fParticle,p_x,p_y,p_z);
				newvertex->SetPrimary( newparticle );
				event->AddPrimaryVertex( newvertex );    

				BxOutputVertex::Get()->SetDId(1);     
				BxOutputVertex::Get()->SetDPosition(fPosition);  
				BxOutputVertex::Get()->SetDDirection(G4ThreeVector(p_x,p_y,p_z)); 
				BxOutputVertex::Get()->SetDPDG(22);
				BxOutputVertex::Get()->SetDTime(0.);
				BxOutputVertex::Get()->SetDEnergy(KinE);		
				BxOutputVertex::Get()->SetDaughters();  
			}
		}
	}

}

//Spectra from A. D. Vijaya and Arun Kumar, The neutron spectrum of Am-Be neutron source, nuclear instruments
//and methods 3 (1973) 435-440 and 
//T. Vilaithong et al., J. Sci. Soc. Thailand, 9 (1983) 129-142.

//-------------------------------------------------------------------------
//    Energy Shooter for Am-Be neutron without Gamma: group 0
//-------------------------------------------------------------------------

G4double BxGeneratorAmBeSource::ShootEnergyNeutron0G() {
   G4double sum = 0;
   G4double norma =0; 
   G4double Probability [351];
   G4double EnergyBin [351];
   G4double deltaX,x,y;
   G4double En=0;

   G4double NeutronEnergySp [351];
   
   NeutronEnergySp [0]=0;
   NeutronEnergySp [1]=0;
      
    


	for(G4int i = 2; i< 351 ; i++)
        
{

En=i*0.01;

NeutronEnergySp [i] = (exp(-7.828*En))/En + 0.40094* exp(-2.23*En) - 7.7 *
pow(10.0,-15.0) * pow(En,-2.831) ;

}



   norma= 0;

     for(G4int i = 0; i< 351 ; i++) {
      norma += NeutronEnergySp[i];
     }

    for(G4int i = 0; i< 351 ; i++) {
      sum += NeutronEnergySp[i]/norma;
      Probability[i]=sum;
      EnergyBin[i]=G4float(i)*10;
    }

  G4double val = G4UniformRand()  ;

  for(G4int i=0;i<351;i++) {
    if(Probability[i] >= val) {

      if(i == 0) {  
      return EnergyBin[0];
}

      deltaX = val - Probability[i] ;
      y = EnergyBin[i] - EnergyBin[i-1] ;
      x = Probability[i] - Probability[i-1] ;	  
	  
         return ((deltaX*y/x + EnergyBin[i])) ;                      //!!!
	 
  //      return (100) ;                                              //!!!
		  
	  
    }
  }

    return 0;
}

//-------------------------------------------------------------------------
//    Energy Shooter for Am-Be neutron with 1 Gamma
//-------------------------------------------------------------------------

G4double BxGeneratorAmBeSource::ShootEnergyNeutron1G() {

	G4double sum = 0;
	G4double norma =0; 
	G4double Probability [53];
	G4double EnergyBin [53];
	G4double deltaX,x,y;

   G4double neutronEnergySp [53] =    
   {
	   0.0,0.0,0.0,0.0,0.0, //End 1 MeV
	   0.00 , 
	   0.4 ,
	   1.2 ,            
	   2.4 ,
	   4.0, //End 2 Mev 
	   5.6,
	   6.6,
	   7.1,
	   7.8,
	   8.8,
	   10,
	   11.25,
	   10.75,
	   10,
	   9.0,
	   8.4,
	   8.0,
	   8.0,
	   8.7,
	   9.0,
	   9.3,
	   8.7,
	   8.1,
	   7,
	   6.2,
	   4.5,
	   2.7,
	   0.2,
	   0.0,0.0,0.0,0.0,0.0,
	   0.0,0.0,0.0,0.0,0.0,
	   0.0,0.0,0.0,0.0,0.0,
	   0.0,0.0,0.0,0.0,0.0
   };  

   norma= 0;

     for(G4int i = 0; i< 53 ; i++) {
      norma += neutronEnergySp[i];
     }

     for(G4int i = 0; i< 53 ; i++) {
	     sum += neutronEnergySp[i]/norma;
	     Probability[i]=sum;
	     EnergyBin[i]=G4double(i)*200;
     }

     G4double val = G4UniformRand()  ;

     for(G4int i=0;i<53;i++) {
	     if(Probability[i] >= val) {

		     if(i == 0) {
			     return EnergyBin[0]/1000;
		     }
		     deltaX = val - Probability[i] ;
		     y = EnergyBin[i] - EnergyBin[i-1] ;
		     x = Probability[i] - Probability[i-1] ;

		     return ((deltaX*y/x + EnergyBin[i])/1000) ;

	     }
     }
     return 0;
}

//-------------------------------------------------------------------------
//    Energy Shooter for AmBe neutron: group 2
//-------------------------------------------------------------------------

G4double BxGeneratorAmBeSource::ShootEnergyNeutron2G() {

	G4double sum = 0;
	G4double norma =0; 
	G4double Probability [58];
	G4double EnergyBin [58];
	G4double deltaX,x,y;

   G4double neutronEnergySp [58] =    
   {
	   0.0,0.0,
	   0.0,
	   0.3,
	   1.3,
	   2.2,
	   2.0,
	   2.06,
	   2.25,
	   1.9,
	   1.6,
	   1.3,
	   1.25,
	   1.0,
	   0.5,
	   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
	   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
	   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
	   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
	   0.0,0.0,0.0
   };  

   norma= 0;

     for(G4int i = 0; i< 58 ; i++) {
      norma += neutronEnergySp[i];
     }

     for(G4int i = 0; i< 58 ; i++) {
	     sum += neutronEnergySp[i]/norma;
	     Probability[i]=sum;
	     EnergyBin[i]=G4double(i)*200;
     }

     G4double val = G4UniformRand()  ;

     for(G4int i=0;i<58;i++) {
	     if(Probability[i] >= val) {

		     if(i == 0) {
			     return EnergyBin[0]/1000;
		     }

		     deltaX = val - Probability[i] ;
		     y = EnergyBin[i] - EnergyBin[i-1] ;
		     x = Probability[i] - Probability[i-1] ;

		     return ((deltaX*y/x + EnergyBin[i])/1000) ;

	     }
     }
     return 0;
}
