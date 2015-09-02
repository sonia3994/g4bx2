// --------------------------------------------------------------------------//
/** 
 * AUTHOR: M. Meyer (UHH)
 * Anti-Neutrino Generator 
 */
// --------------------------------------------------------------------------//

#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TROOT.h"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4Electron.hh"
#include "G4Neutron.hh"
#include "G4SystemOfUnits.hh"

#include "BxLogger.hh"
#include "BxDataCollection.hh"
#include "BxOutputVertex.hh"
#include "BxGeneratorSterileAntiNu.hh"
#include "BxGeneratorSterileAntiNuMessenger.hh"

#include "HistoManager.hh"
#include "TSystem.h"
#include <string>
//---------------------------------------------------------------------------//


BxGeneratorSterileAntiNu::BxGeneratorSterileAntiNu(): BxVGenerator("BxGeneratorSterileAntiNu") {
	fVolumeFlag = false ;
	fEnergyDistribution = false ;
	G4ThreeVector zero(0., 0., 0.) ;
	fSPSAng = new G4SPSAngDistribution ;
	fSPSPos = new G4SPSPosDistribution ;
	fSPSEne = new G4SPSEneDistribution ;
	G4SPSRandomGenerator *RndGen = new G4SPSRandomGenerator;
	fSPSPos->SetBiasRndm(RndGen);
	fSPSAng->SetBiasRndm(RndGen);
	fSPSEne->SetBiasRndm(RndGen);
	SetRadius0(0.);
	SetCentreCoords(zero);   
	fshapeCoefficientA=0;
	fshapeCoefficientB=0;
	fshapeCoefficientC=0;

	fParticleTable = G4ParticleTable::GetParticleTable();    


	fTheMessenger = new BxGeneratorSterileAntiNuMessenger(this);
	fParticleGun  = new G4ParticleGun ;  


	string CrossSectioSpectrumRootfile ="../data/dat/File_CrossSection_Spectrum.root"; 
	string Pr144SpectrumRootfile = "../data/dat/Pr144Spectrum_and_crossSection_Feb2015_keV.root";    

	/// PDF building of Neutrino Spectrum and Cross Section (IBD, first order)	
	fCrossSectionFile = new TFile (CrossSectioSpectrumRootfile.c_str(), "RECREATE");	// currently not active (might be changed in the future)
	fSpectrumRootFile	= new TFile(Pr144SpectrumRootfile.c_str(),"readonly");		// antineutrino spectrum x cross section


	fSpectrumRootFile->GetObject("hist_spectrum",hist_AntiNuSpec);


	if (hist_AntiNuSpec==0) {

		BxLog(fatal) << "WARNING: Neutrino Source Spectrum Histogram NOT FOUND!" << endlog;
	}

	/// Histogram including shape deformation (added February 2015):
	shapeFactorInit = false;
	// cout << "TOP: shapeFactorInit " << shapeFactorInit << endl;
	hist_AntiNuSpec_ShapeFactor = (TH1F*)hist_AntiNuSpec->Clone("hist_AntiNuSpec_ShapeFactor");
	hist_AntiNuSpec_ShapeFactor->Reset();

	
	

	fCrossSectionFile->cd();

	hist_AntiNuSpec->Write();
	//hist_AntiNuSpec_ShapeFactor->Write();

	BxLog(routine) << "PDF read-in of antineutrino spectrum and cross-section complete!" << endlog;

	/// END of PDF building

	if(!fParticleGun) {
		BxLog(error) << "Could not allocate G4ParticleGun! Out of memory?"<<endlog;
		BxLog(fatal) << endlog;
	}
	BxLog(routine) << "G4ParticleGun Constructed." << endlog;
}

//---------------------------------------------------------------------------//

BxGeneratorSterileAntiNu::~BxGeneratorSterileAntiNu()
{
  delete fTheMessenger;
  delete fParticleGun;
  delete fSPSPos;
  delete fSPSAng;
  delete fSPSEne;
  delete fCrossSectionFile;
  delete fSpectrumRootFile; 
}

//---------------------------------------------------------------------------//


void BxGeneratorSterileAntiNu::BxGeneratePrimaries(G4Event *event) {

	/** Comments: (not active, only antiNu)
	  SimulationType = 1 --> Neutrinos
	  SimulationType = 2 --> Anti-Neutrinos

	 **/
	BxLog(debugging) << "DOWN-A: shapeFactorInit " << shapeFactorInit << endlog;
	initShapeFactor(shapeFactorInit);
	BxLog(debugging) << "DOWN-B: shapeFactorInit " << shapeFactorInit << endlog;
	/*
	G4cout << endl;
	G4cout << "fshapeCoefficientA = " << fshapeCoefficientA << endl;
	G4cout << "fshapeCoefficientB = " << fshapeCoefficientB << endl;
	G4cout << "fshapeCoefficientC = " << fshapeCoefficientC << endl;
	G4cout << endl; 
	*/

	
	G4int 	SimulationType 	= 2;		
	survivalProb(3., 3.);



	if(fVolumeFlag)  {
		fSPSAng->SetVerbosity(0);
		fNuPosition = fSPSPos->GenerateOne(); 	
		SetParticlePosition(fSPSPos->GenerateOne());  
		fSPSAng->SetAngDistType("iso");    
		fSPSAng->SetPosDistribution(fSPSPos);
		SetParticleMomentumDirection(fSPSAng->GenerateOne());    
	}
	if(fEnergyDistribution) {
		fSPSEne->SetVerbosity(3);
		SetParticleEnergy(fSPSEne->GenerateOne(fParticleGun->GetParticleDefinition()));
	}



	/// **********************************************************************************************************************
	/// random generator:

	TRandom *rand = new TRandom3(0);	
	gRandom = rand;	

	/// ######################################################################################################################




	/// Generator...
	if ((event->GetEventID())%2 == 0){ 
		fParticle = fParticleTable->FindParticle(-11);	

		G4double nuEnergy 		= 0;
		//bool interaction		= true;	
		G4double oscillationBaseline 	= 0.;	

		G4bool saveEvent		= false;
		fPosition			= GetInteractionPositionIntSource(fNuPosition);
		G4ThreeVector distance		= fNuPosition - fPosition;				
		oscillationBaseline		= distance.mag()/m;					
		// G4cout << "oscillationBaseline: " << oscillationBaseline << endl;
		nuEnergy 			= shootNuEnergy(SimulationType);
		if(G4UniformRand() < survivalProb(nuEnergy/MeV, distance.mag()/m)) saveEvent=true;	

		/// for testing
		G4double	eEnergy 		= shootEnergy(nuEnergy/MeV, SimulationType);				
		G4ThreeVector 	eDirection 		= fSPSAng->GenerateOne(); 
		G4ThreeVector 	nDirection 		= fSPSAng->GenerateOne();    

		/// Events will only be stored, if oscillation probability allows it...:

		if (saveEvent)  {

			G4double kinE	=	eEnergy;	


			G4ParticleDefinition* positron;
			positron 			= fParticleTable->FindParticle(-11);	
			G4ThreeVector eMom;
			G4double px 			= eMom.x();
			G4double py 			= eMom.y();
			G4double pz 			= eMom.z();	

			G4double energy			= eEnergy + positron->GetPDGMass(); 	
			G4double pmom 			= (std::sqrt(pow(energy,2.) - pow(positron->GetPDGMass(),2.))); 
			px 				= pmom*eDirection.x();
			py 				= pmom*eDirection.y();
			pz 				= pmom*eDirection.z();

			G4PrimaryVertex* vertex 	= new G4PrimaryVertex(fPosition,0);
			G4PrimaryParticle* particle	= new G4PrimaryParticle(positron,px,py,pz);			

			vertex->SetPrimary(particle);
			event->AddPrimaryVertex(vertex);

			BxOutputVertex::Get()->SetIsotopeCoinc(0);	


			G4double  positronKineticEnergy = eEnergy;


			/// May 2014:
			/// Attention:
			/*
			   For this generator, the position of the events corresponds
			   to the position of the neutrino source. You can also gener-
			   ate events spatially distributed in a source. Neutrino event
			   characteristics are stored in the vertex. Electrons are au-
			   tomatically generated within the specified FV, and their
			   properties are stored in the daughter structure
			   */	

			fParticleGun->SetParticleEnergy(positronKineticEnergy*MeV);	
			fParticleGun->SetParticlePosition(fPosition);

			fParticleGun->SetParticleMomentumDirection(eDirection);
			fParticleGun->GeneratePrimaryVertex(event);

			// neutrino
			BxOutputVertex::Get()->SetPDG(-12);
			BxOutputVertex::Get()->SetTime(0.);
			BxOutputVertex::Get()->SetEnergy(nuEnergy/MeV);
			BxOutputVertex::Get()->SetPosition(fNuPosition);
			BxOutputVertex::Get()->SetDirection(fDirection);	
			// electron/positron
			BxOutputVertex::Get()->SetDEnergy(eEnergy/MeV);		
			BxOutputVertex::Get()->SetDPDG(fParticleGun->GetParticleDefinition()->GetPDGEncoding());		
			BxOutputVertex::Get()->SetDPosition(fParticleGun->GetParticlePosition());  
			BxOutputVertex::Get()->SetDDirection(fParticleGun->GetParticleMomentumDirection()); 
			BxOutputVertex::Get()->SetDaughters();		



			/// Neutron:
			G4ParticleDefinition* particleNeutron;
			G4ParticleDefinition* particleProton;
			G4ParticleDefinition* particlePositron;

			particleNeutron 	= fParticleTable->FindParticle(2112);
			particleProton		= fParticleTable->FindParticle(2212);
			particlePositron	= fParticleTable->FindParticle(-11);

			G4double neutronMass	= particleNeutron->GetPDGMass();
			G4double protonMass	= particleProton->GetPDGMass();

			fParticle = fParticleTable->FindParticle(2112);     


			kinE = nuEnergy + protonMass - neutronMass - energy;
			energy = kinE + fParticle->GetPDGMass();

			fDirection = fSPSAng->GenerateOne();

			pmom = std::sqrt(pow(energy,2.) - pow(fParticle->GetPDGMass(),2.));
			px = pmom*fDirection.x();
			py = pmom*fDirection.y();
			pz = pmom*fDirection.z();

			vertex   = new G4PrimaryVertex(fPosition,0);
			particle = new G4PrimaryParticle(fParticle,px,py,pz);

			vertex->SetPrimary( particle );
			event->AddPrimaryVertex( vertex );    
	//HistoManager::Get()->FillNtuple(eEnergy,fPosition.x(),fPosition.y(),fPosition.z(),fNuPosition.x(),fNuPosition.y(),fNuPosition.z());
	
	}
	
	
	
   }

}

//---------------------------------------------------------------------------//

/// Additional methods:

/*
G4double BxGeneratorSterileAntiNu::ShootEnergyNeutron() {
	// adapted from geo-neutrino generator
	G4double sum = 0;
	G4double norma =0; 
	G4double Probability [51];
	G4double EnergyBin [51];
	G4double deltaX,x,y;

	G4double NeutronEnergySp [51] =    
	{ 0., 274.4119, 342.3826, 257.9844, 187.8962, 141.1938,
		107.2338,  78.0118,  62.4729,  45.4234,  37.3674,
		28.1552,  27.4770,  17.0363,  17.4592,   9.9503,
		9.9544,   7.9960,   4.7464,   5.1454,   3.6949,
		2.1576,   1.7950,   1.1628,   1.1344,   0.9440,
		0.5982,   0.4545,   0.3037,   0.3071,   0.2996,
		0.1526,   0.1976,   0.0576,   0.1063,   0.0565,
		0.0597,   0.0334,   0.0270,   0.0326,   0.0113,
		0.0121,   0.0025,   0.0063,   0.0022,   0.0011,
		0.0000,   0.0000,   0.0000,   0.0000,   0.0000};
	//     l'energie est au milieu du bin : 2kev est au milieu du bin de
	//                             0kev a 4 kev   qui contient 274.4119  

	norma= 0;

	for(G4int i = 0; i< 51 ; i++) {
		norma += NeutronEnergySp[i];
	}

	for(G4int i = 0; i< 51 ; i++) {
		sum += NeutronEnergySp[i]/norma;
		Probability[i]=sum;
		EnergyBin[i]=G4float(i)*4;
	}

	G4double val = G4UniformRand()  ;

	for(G4int i=0;i<51;i++) {
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
*/
//---------------------------------------------------------------------------//

G4double BxGeneratorSterileAntiNu::survivalProb(G4double E, G4double L) {
	if(fSin2theta2 == 0)  return 1. ;	
	if(fDeltaM2  == 0)    return 1. ;	

	G4double survivalProbability = 1 - fSin2theta2*pow(sin(1.276 * fDeltaM2*L/E),2);


	return survivalProbability; //Test

}

//---------------------------------------------------------------------------//

/*
   G4ThreeVector BxGeneratorSterileAntiNu::shootE_Direction(G4ThreeVector nuPos, G4ThreeVector interactionPos, G4double eEnergy, G4double nuEnergy) {

   G4double emass  = G4Electron::Definition()->GetPDGMass ()/MeV;
// Shoot Electron/Positron direction:
G4ThreeVector rr(0,0,0);	// Vector (source Point and interaction Point)
rr = interactionPos - nuPos;

// "Recoil angle":
// For definition see http://www-sk.icrr.u-tokyo.ac.jp/sk/ykphd/chap3-2.html
//G4double theta = acos ( (1.+emass/nuEnergy) / (sqrt(1+2.*emass/eEnergy) ) );


}
*/

G4double BxGeneratorSterileAntiNu::shootNuEnergy(G4int k) {

	AntiNuEnergyMax	= 2.999;
	AntiNuThreshold 	= 1.806;
	if (k==1) {
		BxLog(debugging) << "Selected Generator: Cr-51" << endlog;
	}
	if (k==2) {
		BxLog(debugging) << "Selected Generator: Anti-Neutrino (Ce-144 / Pr-144)" << endlog;
	}
	/// Cr-51:
	G4double val = G4UniformRand()  ;
	// G4cout << "val = " << val << endl;
	if(k == 1) { // vorher Cr51
		if(val < 0.09) return 0.751; 
		if(val < 0.90) return 0.746; 
		if(val < 0.91) return 0.431; 
		if(val < 1.00) return 0.426;     
	} 
	/// Anti-Neutrino:
	G4double nuEnergy=100.;
	if (k==2) {


		TRandom *rand = new TRandom3(0);	
		gRandom = rand;


		nuEnergy = hist_AntiNuSpec_ShapeFactor->GetRandom() / 1000. *MeV;	// keV --> MeV	


		return nuEnergy;


	}


	return 0;

}

//---------------------------------------------------------------------------//

G4double BxGeneratorSterileAntiNu::shootEnergy(G4double NuEnergy, G4int simType) {

	if (simType==1){
		BxLog(debugging) << "Positron Energy calculation for scattering ve" << endlog;
	}
	if (simType==2){
		BxLog(debugging) << "Positron Energy calculation for IBD" << endlog;
	}

	if(simType==1){

		G4double 	NuEndPoint 	= 0.75273;
		G4int		numberOfSteps 	= 500;	

		G4double fStep = (NuEndPoint*1.01)/numberOfSteps;
		G4double fNorma = 0 ;
		for(G4int i=0;i<numberOfSteps;++i) {
			if(fStep*i+fStep/3. > NuEndPoint) break ;
			fNorma += cross_nue(fStep*i+fStep/3.,NuEnergy)*fStep;
		}

		G4double value = G4UniformRand();

		G4double sum = 0 ;
		for(G4int i=0;i<numberOfSteps;++i) {
			G4double bin = fStep*i+fStep/3. ;
			sum += cross_nue(bin,NuEnergy)*fStep/fNorma;
			if(sum>value) return bin ;
		}
	}

	if(simType==2){
		G4ParticleDefinition* particleNeutron;
		G4ParticleDefinition* particleProton;
		G4ParticleDefinition* particlePositron;

		particleNeutron 	= fParticleTable->FindParticle(2112);
		particleProton		= fParticleTable->FindParticle(2212);
		particlePositron	= fParticleTable->FindParticle(-11);

		G4double neutronMass	= particleNeutron->GetPDGMass();
		G4double protonMass	= particleProton->GetPDGMass();
		G4double positronMass	= particlePositron->GetPDGMass();

		G4double Delta		= neutronMass-protonMass;

		G4double positronEnergy	= 0.;	
		G4double positronKineticEnergy = 0.;

		/// FIXME 22.5: positronEnergy 		= NuEnergy-Delta;
		positronEnergy 		= 1./2. * (std::sqrt(neutronMass*neutronMass - 4.*protonMass*(-NuEnergy+Delta+(Delta*Delta-positronMass*positronMass)/(2.*protonMass))) - neutronMass);

		positronKineticEnergy = positronEnergy - positronMass;

		// pow(energy,2.)

		G4bool CoutON = false;

		if (CoutON) {
			// Ausgabe:
			G4cout << endl;
			G4cout << "neutronMass     = " << neutronMass<< endl;
			G4cout << "protonMass  	   = " << protonMass << endl;
			G4cout << "Delta      	   = " << Delta << endl;
			G4cout << "NuEnergy        = " << NuEnergy << endl;
			G4cout << "positronEnergy  = " << positronEnergy << endl;
			G4cout << endl;
			// Histogramm:
			hist_truePositronEnergy = new TH1F("hist_truePositronEnergy","hist_truePositronEnergy",100,0.,3.5);
			fCrossSectionFile->cd();
			hist_truePositronEnergy->Fill(positronEnergy);
			//hist_truePositronEnergy->Write(); // FIXME
		}
		return positronKineticEnergy;

	}



	return 0;

}

//---------------------------------------------------------------------------//

G4double BxGeneratorSterileAntiNu::cross_nue(G4double E, G4double nuE) {

	G4double emass  = G4Electron::Definition()->GetPDGMass ()/MeV;
	G4double gr = 0.23 ;
	G4double gl = 0.5 + gr ;
	G4double sigma = 0;
	if(E > nuE/(1+emass/2./nuE)) return 0 ;
	else  sigma = gl*gl+gr*gr*pow(1.-E/nuE,2.) - gl*gr*E*emass/nuE/nuE ;

	if(sigma < 0) sigma = 0 ;

  return sigma ;
}

//---------------------------------------------------------------------------//

G4ThreeVector BxGeneratorSterileAntiNu::GetInteractionPositionIntSource(G4ThreeVector nupos) {

	G4double NuPosx0 = nupos.x();						// Position of emmision
	G4double NuPosy0 = nupos.y(); 
	G4double NuPosz0 = nupos.z();

	G4double fFVRadius = 4.25*m;		

	G4double R = nupos.mag() + fFVRadius ;

	G4double r = nupos.mag() + fFVRadius + 1;

	G4ThreeVector X(nupos.mag() + fFVRadius + 1,0,0);

	while(X.mag() > fFVRadius) {						// dice till inside FV

		G4double phi 		= 2.*CLHEP::pi* G4UniformRand();	// between 0° und 360°
		G4double v		= 2.0*G4UniformRand()-1.0;	
		G4double theta		= acos(v);

		r  			= R*G4UniformRand() ;			// point source

		X.setX(r*sin(theta) * cos(phi) + NuPosx0);
		X.setY(r*sin(theta) * sin(phi) + NuPosy0);
		X.setZ(r*cos(theta) + NuPosz0);

	}

	return X;
}


G4bool BxGeneratorSterileAntiNu::initShapeFactor(G4bool isSet) {
	
	G4bool makeHistoForCheck  = true;
	if (isSet==true) return true;
	else { 
		G4int numberOfBinsForClone = 0;
		numberOfBinsForClone = hist_AntiNuSpec->GetSize()-2; 
	
		// Set shape coefficients: 
		G4float spectrumShapeCoefficient_A = 0.;// .01;// -0.058;
		G4float spectrumShapeCoefficient_B = 0.;//.35;//+0.389;
		G4float spectrumShapeCoefficient_C = 0.;//.0;// 0.;//0.;		

		spectrumShapeCoefficient_A = fshapeCoefficientA;
		spectrumShapeCoefficient_B = fshapeCoefficientB;
		spectrumShapeCoefficient_C = fshapeCoefficientC;	
	
		// Building deformation:
		for (G4int iter=0; iter<numberOfBinsForClone; iter++) {
	
			G4float energy_kev 			= (G4float)hist_AntiNuSpec->GetBinCenter(iter+1);
			G4float spectrumPDF			= (G4float)hist_AntiNuSpec->GetBinContent(iter+1);
			G4float shapeFactor			= (G4float)(1. + spectrumShapeCoefficient_A*energy_kev 
								+ spectrumShapeCoefficient_B/energy_kev 
								+ spectrumShapeCoefficient_C*energy_kev*energy_kev);
			G4float spectrumPDFwithShapeFactor	= spectrumPDF*shapeFactor;		
	
			hist_AntiNuSpec_ShapeFactor->SetBinContent(iter+1, spectrumPDFwithShapeFactor);	
			
		}
		shapeFactorInit=true;  
		normalizeHisto(hist_AntiNuSpec);		// should already be normalized
		normalizeHisto(hist_AntiNuSpec_ShapeFactor);	// will be normalized for optical inspection, otherwise not necessary
		if (makeHistoForCheck) {
			TFile *fout = new TFile("checkAntiNuSpectrum.root","RECREATE");		
			hist_AntiNuSpec->Write();
			hist_AntiNuSpec_ShapeFactor->Write();
			fout->Close();
		}
  		return true;
	}


  	
}

G4bool BxGeneratorSterileAntiNu::normalizeHisto(TH1F* hist) {

	G4int binnr = hist->GetSize()-2;

	BxLog(debugging)  << endl;
	BxLog(debugging)  <<"Normalize Histo..." << endl;
	BxLog(debugging)  << "Bins: " << binnr << endl;
	G4double sum = 0;
	G4double content;
	
	for(G4int i=1; i<binnr+1; i++) {
		content = hist->GetBinContent(i);
		sum += content;
	}
	BxLog(debugging)  << "Factor for Normalization: " << sum << endl;
	hist->Scale(1./sum);

	return true;
}





/*
 * $Log: BxGeneratorSterileAntiNu.cc,v $
 * Revision 1.2  2015/02/25 18:17:18  acaminata
 * Shape factor added
 *
 * Revision 1.1  2015/02/12 14:31:35  acaminata
 * Sterile generator added
 *
 * Revision 1.12  2009-10-20 12:37:37  dfranco
 * Added a new generator in order to read a g4bx output file and begin a simulation
 * from the energy deposits. Useful for external background.
 *
 * Revision 1.11  2007-11-12 12:09:40  dfranco
 * added to g4gun, the following  energy distributions:
 * Lin (linear), Pow (power-law), Exp (exponential), Gauss (gaussian),
 * Brem (bremsstrahlung), BBody (black-body), Cdg (cosmic diffuse gamma-ray)
 *
 * Revision 1.10  2007-04-26 16:58:26  dfranco
 * Development of the new output format
 *
 * Revision 1.9  2007-04-26 08:51:25  dfranco
 * Added the structure for the new ouput format. It is not yet active!!
 *
 * Revision 1.8  2007-04-12 12:08:33  dfranco
 * Quenching formalization modified: now, for each charged particles, at each step, dE/dx is evaluated. The Birks equation is applied to the light response. The remaining kinetic energy of the particle is not affected by the quenching effect.
 * New commands for modifying Birks (and second order Birks) parameter:
 * /bx/detector/birks xxx
 * /bx/detector/birks2 xxx
 * Both are already in cm/MeV. By default, birks = 0.0085 and birks2 = 0. They must be tuned.
 * Photon yield is set to 18800, based on the last Borexino results. Fine tuning of these parameters is required.
 *
 * Revision 1.7  2007-03-22 14:48:55  dfranco
 * Stable version
 *
 */
