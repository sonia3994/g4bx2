//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Created by A. Caminata and S. Marcocci, Sept. 2014
//
//
////////////////////////////////////////////////////////////////////////
// Cerenkov Radiation Class Implementation
////////////////////////////////////////////////////////////////////////
#include "HistoManager.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Poisson.hh"
#include "G4EmProcessSubType.hh"
#include "G4LossTableManager.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"
#include "BxReadParameters.hh"
#include "BxCerenkov.hh"
#include <vector>
/////////////////////////
// Class Implementation  
////////////////////////

BxCerenkov::BxCerenkov(const G4String& processName, G4ProcessType type)
           : G4VProcess(processName, type)
{
        SetProcessSubType(fCerenkov);

	fTrackSecondariesFirst = false;
        fMaxBetaChange = 0.;
	fMaxPhotons = 0;

        thePhysicsTable = NULL;

	if (verboseLevel>0) {
           G4cout << GetProcessName() << " is created " << G4endl;
	}
}

// BxCerenkov::BxCerenkov(const BxCerenkov &right)
// {
// }

        ////////////////
        // Destructors
        ////////////////

BxCerenkov::~BxCerenkov() 
{
	if (thePhysicsTable != NULL) {
	   thePhysicsTable->clearAndDestroy();
           delete thePhysicsTable;
	}
}

        ////////////
        // Methods
        ////////////

G4bool BxCerenkov::IsApplicable(const G4ParticleDefinition& aParticleType)
{
    G4bool result = false;
    if (aParticleType.GetPDGCharge() != 0.0 && 
	aParticleType.GetPDGMass() != 0.0 &&
	aParticleType.GetParticleName() != "chargedgeantino" &&
	!aParticleType.IsShortLived() ) { result = true; }

    return result;
}

void BxCerenkov::BuildPhysicsTable(const G4ParticleDefinition&)
{
    if (!thePhysicsTable) BuildThePhysicsTable();
}

// PostStepDoIt
// -------------
//
G4VParticleChange*
BxCerenkov::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)

// This routine is called for each tracking Step of a charged particle
// in a radiator. A Poisson-distributed number of photons is generated
// according to the Cerenkov formula, distributed evenly along the track
// segment and uniformly azimuth w.r.t. the particle direction. The 
// parameters are then transformed into the Master Reference System, and 
// they are added to the particle change. 

{
	//////////////////////////////////////////////////////
	// Should we ensure that the material is dispersive?
	//////////////////////////////////////////////////////

        aParticleChange.Initialize(aTrack);
        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
        const G4Material* aMaterial = aTrack.GetMaterial();

	G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
	G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

	G4ThreeVector x0 = pPreStepPoint->GetPosition();
        G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
	G4double t0 = pPreStepPoint->GetGlobalTime();

        G4MaterialPropertiesTable* aMaterialPropertiesTable =
                               aMaterial->GetMaterialPropertiesTable();
        if (!aMaterialPropertiesTable) return pParticleChange;

	G4MaterialPropertyVector* Rindex = 
                aMaterialPropertiesTable->GetProperty("RINDEX"); 
        if (!Rindex) return pParticleChange;

	G4bool RindIsGot;

	//here the maximum value of the refraction index is calculated
	G4double nMax=Rindex->GetValue(Rindex->GetLowEdgeEnergy(0),RindIsGot);
	for(size_t i=0 ; i<Rindex->GetVectorLength()-1;i++){
		if(Rindex->GetValue(Rindex->GetLowEdgeEnergy(i),RindIsGot)>nMax)
			nMax=Rindex->GetValue(Rindex->GetLowEdgeEnergy(i),RindIsGot);
	}


	// particle charge
	const G4double charge = aParticle->GetDefinition()->GetPDGCharge();

	// particle beta
	const G4double beta = (pPreStepPoint ->GetBeta() +
			pPostStepPoint->GetBeta())/2.;

	vector <G4double> MeanNumberOfPh;
	vector <G4double> MomentumRange;
	vector <G4int> NumPhotons;
	GetAverageNumberOfPhotons(charge,beta,aMaterial,Rindex, MeanNumberOfPh, MomentumRange);

	G4double TotalMeanNumberOfPhotons=0;
	for(size_t i=0;i<MeanNumberOfPh.size();i++)
		TotalMeanNumberOfPhotons+=MeanNumberOfPh.at(i);

	if (TotalMeanNumberOfPhotons<= 0.0) {
		// return unchanged particle and no secondaries
		aParticleChange.SetNumberOfSecondaries(0);
		return pParticleChange;
	}

	G4double step_length;
	step_length = aStep.GetStepLength();

	G4double MeanNumberOfPhotons=0; 
	for(size_t i=0;i<MeanNumberOfPh.size();i++)
		NumPhotons.push_back((G4int) G4Poisson(MeanNumberOfPh.at(i) * step_length*BxReadParameters::Get()->GetPMTQEMaximum()));

	for(size_t i=0;i<MeanNumberOfPh.size();i++)
		MeanNumberOfPhotons+=NumPhotons.at(i);

	if (MeanNumberOfPhotons <= 0) {

		// return unchanged particle and no secondaries  

		aParticleChange.SetNumberOfSecondaries(0);

		return pParticleChange;
	}

		////////////////////////////////////////////////////////////////

	aParticleChange.SetNumberOfSecondaries(MeanNumberOfPhotons);

	if (fTrackSecondariesFirst) {
		if (aTrack.GetTrackStatus() == fAlive )
			aParticleChange.ProposeTrackStatus(fSuspend);
	}




	for(size_t IntegratingRegionIterator=0;IntegratingRegionIterator<MeanNumberOfPh.size();IntegratingRegionIterator++){
//G4cout<<G4endl;
//G4cout<<"------------------------------------------------------"<<G4endl;
//G4cout<<" CICLO DI ITERAZIONE DEL CHERENKOV NUMERO "<<IntegratingRegionIterator<<G4endl;
//G4cout<<"------------------------------------------------------"<<G4endl;

//		G4cout<<" Beta = "<<beta<<G4endl;
//		G4cout<<"STEP LENGTH = "<<aStep.GetStepLength()<<G4endl;

//		G4cout<<" NUMBER OF CHERENKOV PHOTONS GENERATED "<< NumPhotons <<G4endl;
//		G4cout<<" Pmin "<< MomentumRange.at(2*IntegratingRegionIterator)<<" Pmax "<<MomentumRange.at(2*IntegratingRegionIterator+1)<<" dp "<<MomentumRange.at(2*IntegratingRegionIterator+1)-MomentumRange.at(2*IntegratingRegionIterator) <<G4endl;

		////////////////////////////////////////////////////////////////
		G4double Pmin = MomentumRange.at(2*IntegratingRegionIterator);
		G4double Pmax = MomentumRange.at(2*IntegratingRegionIterator+1);
		G4double dp = Pmax - Pmin;

		G4double BetaInverse = 1./beta;

		G4double maxCos = BetaInverse / nMax; 
		G4double maxSin2 = (1.0 - maxCos) * (1.0 + maxCos);

		const G4double beta1 = pPreStepPoint ->GetBeta();
		const G4double beta2 = pPostStepPoint->GetBeta();

		vector<G4double> MeanNumberOfPh1;
		vector<G4double> dummy;
		G4double MeanNumberOfPhotons1=0;
		GetAverageNumberOfPhotons(charge,beta1,aMaterial,Rindex, MeanNumberOfPh1, dummy);
		for(size_t i=0;i<MeanNumberOfPh1.size();i++){
			MeanNumberOfPhotons1+=MeanNumberOfPh1.at(i);
		}
		vector<G4double> MeanNumberOfPh2;
		G4double MeanNumberOfPhotons2=0; 
		GetAverageNumberOfPhotons(charge,beta2,aMaterial,Rindex,MeanNumberOfPh2, dummy);
		for(size_t i=0;i<MeanNumberOfPh2.size();i++){
			MeanNumberOfPhotons2+=MeanNumberOfPh2.at(i);
		}
		for (G4int i = 0; i < NumPhotons.at(IntegratingRegionIterator); i++) {

			// Determine photon energy

			G4double rand;
			G4double sampledEnergy, sampledRI; 
			G4double cosTheta, sin2Theta;

			// sample an energy
			do {
				rand = G4UniformRand();	
				sampledEnergy = Pmin + rand * dp; 
				sampledRI = Rindex->Value(sampledEnergy);
				cosTheta = BetaInverse / sampledRI;  

				sin2Theta = (1.0 - cosTheta)*(1.0 + cosTheta);
				rand = G4UniformRand();	

			} while (rand*maxSin2 > sin2Theta);

			// Generate random position of photon on cone surface 
			// defined by Theta 

			rand = G4UniformRand();

			G4double phi = twopi*rand;
			G4double sinPhi = std::sin(phi);
			G4double cosPhi = std::cos(phi);

			// calculate x,y, and z components of photon energy
			// (in coord system with primary particle direction 
			//  aligned with the z axis)

			G4double sinTheta = std::sqrt(sin2Theta); 
			G4double px = sinTheta*cosPhi;
			G4double py = sinTheta*sinPhi;
			G4double pz = cosTheta;

			// Create photon momentum direction vector 
			// The momentum direction is still with respect
			// to the coordinate system where the primary
			// particle direction is aligned with the z axis  

			G4ParticleMomentum photonMomentum(px, py, pz);

			// Rotate momentum direction back to global reference
			// system 

			photonMomentum.rotateUz(p0);

			// Determine polarization of new photon 

			G4double sx = cosTheta*cosPhi;
			G4double sy = cosTheta*sinPhi; 
			G4double sz = -sinTheta;

			G4ThreeVector photonPolarization(sx, sy, sz);

			// Rotate back to original coord system 

			photonPolarization.rotateUz(p0);

			// Generate a new photon:

			G4DynamicParticle* aCerenkovPhoton =
				new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), 
						photonMomentum);
			aCerenkovPhoton->SetPolarization
				(photonPolarization.x(),
				 photonPolarization.y(),
				 photonPolarization.z());
			HistoManager::Get()->FillHisto(5, sampledEnergy/eV);
			aCerenkovPhoton->SetKineticEnergy(sampledEnergy);

			// Generate new G4Track object:

			G4double delta, NumberOfPhotons, N;

			do {
				rand = G4UniformRand();
				delta = rand * aStep.GetStepLength();
				NumberOfPhotons = MeanNumberOfPhotons1 - delta *
					(MeanNumberOfPhotons1-MeanNumberOfPhotons2)/
					aStep.GetStepLength();
				N = G4UniformRand() *
					std::max(MeanNumberOfPhotons1,MeanNumberOfPhotons2);
			} while (N > NumberOfPhotons);

			G4double deltaTime = delta /
				((pPreStepPoint->GetVelocity()+
				  pPostStepPoint->GetVelocity())/2.);

			G4double aSecondaryTime =t0 + deltaTime;

			G4ThreeVector aSecondaryPosition =x0 + rand * aStep.GetDeltaPosition();

			G4Track* aSecondaryTrack = 
				new G4Track(aCerenkovPhoton,aSecondaryTime,aSecondaryPosition);

			aSecondaryTrack->SetTouchableHandle(
					aStep.GetPreStepPoint()->GetTouchableHandle());

			aSecondaryTrack->SetParentID(aTrack.GetTrackID());

			aParticleChange.AddSecondary(aSecondaryTrack);
		}

		if (verboseLevel>0) {
			G4cout <<"\n Exiting from BxCerenkov::DoIt -- NumberOfSecondaries = "
				<< aParticleChange.GetNumberOfSecondaries() << G4endl;
		}
	}
	return pParticleChange;

}
// BuildThePhysicsTable for the Cerenkov process
// ---------------------------------------------
//

void BxCerenkov::BuildThePhysicsTable()
{
	if (thePhysicsTable) return;

	const G4MaterialTable* theMaterialTable=
	 		       G4Material::GetMaterialTable();
	G4int numOfMaterials = G4Material::GetNumberOfMaterials();

	// create new physics table
	
	thePhysicsTable = new G4PhysicsTable(numOfMaterials);

	// loop for materials

	for (G4int i=0 ; i < numOfMaterials; i++)
	{
	        G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector = 0;

		// Retrieve vector of refraction indices for the material
		// from the material's optical properties table 

		G4Material* aMaterial = (*theMaterialTable)[i];

		G4MaterialPropertiesTable* aMaterialPropertiesTable =
				aMaterial->GetMaterialPropertiesTable();

		if (aMaterialPropertiesTable) {

		   aPhysicsOrderedFreeVector = new G4PhysicsOrderedFreeVector();
		   G4MaterialPropertyVector* theRefractionIndexVector = 
		    	   aMaterialPropertiesTable->GetProperty("RINDEX");

		   if (theRefractionIndexVector) {
		
		      // Retrieve the first refraction index in vector
		      // of (photon energy, refraction index) pairs 

                      G4double currentRI = (*theRefractionIndexVector)[0];

		      if (currentRI > 1.0) {

			 // Create first (photon energy, Cerenkov Integral)
			 // pair  

                         G4double currentPM = theRefractionIndexVector->
                                                 Energy(0);
			 G4double currentCAI = 0.0;

			 aPhysicsOrderedFreeVector->
			 	 InsertValues(currentPM , currentCAI);

			 // Set previous values to current ones prior to loop

			 G4double prevPM  = currentPM;
			 G4double prevCAI = currentCAI;
                	 G4double prevRI  = currentRI;

			 // loop over all (photon energy, refraction index)
			 // pairs stored for this material  

                         for (size_t ii = 1;
                              ii < theRefractionIndexVector->GetVectorLength();
                              ++ii)
			 {
                                currentRI = (*theRefractionIndexVector)[ii];
                                currentPM = theRefractionIndexVector->Energy(ii);

				currentCAI = 0.5*(1.0/(prevRI*prevRI) +
					          1.0/(currentRI*currentRI));

				currentCAI = prevCAI + 
					     (currentPM - prevPM) * currentCAI;

				aPhysicsOrderedFreeVector->
				    InsertValues(currentPM, currentCAI);

				prevPM  = currentPM;
				prevCAI = currentCAI;
				prevRI  = currentRI;
			 }

		      }
		   }
		}

	// The Cerenkov integral for a given material
	// will be inserted in thePhysicsTable
	// according to the position of the material in
	// the material table. 

	thePhysicsTable->insertAt(i,aPhysicsOrderedFreeVector); 

	}
}

// GetMeanFreePath
// ---------------
//

G4double BxCerenkov::GetMeanFreePath(const G4Track&,
                                           G4double,
                                           G4ForceCondition*)
{
        return 1.;
}

G4double BxCerenkov::PostStepGetPhysicalInteractionLength(
                                           const G4Track& aTrack,
                                           G4double,
                                           G4ForceCondition* condition)
{
        *condition = NotForced;
        G4double StepLimit = DBL_MAX;

        const G4Material* aMaterial = aTrack.GetMaterial();
	G4int materialIndex = aMaterial->GetIndex();

	// If Physics Vector is not defined no Cerenkov photons
	//    this check avoid string comparison below
	if(!(*thePhysicsTable)[materialIndex]) { return StepLimit; }

        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
        const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();

        G4double kineticEnergy = aParticle->GetKineticEnergy();
        const G4ParticleDefinition* particleType = aParticle->GetDefinition();
        G4double mass = particleType->GetPDGMass();

        // particle beta
        G4double beta = aParticle->GetTotalMomentum() /
	                aParticle->GetTotalEnergy();
        // particle gamma
        G4double gamma = aParticle->GetTotalEnergy()/mass;

        G4MaterialPropertiesTable* aMaterialPropertiesTable =
                            aMaterial->GetMaterialPropertiesTable();

        G4MaterialPropertyVector* Rindex = NULL;

        if (aMaterialPropertiesTable)
		Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");

	G4double nMax;
	G4bool RindIsGot;
	if (Rindex) {
		nMax=Rindex->GetValue(Rindex->GetLowEdgeEnergy(0),RindIsGot);
		for(size_t i=0 ; i<Rindex->GetVectorLength()-1;i++){
			if(Rindex->GetValue(Rindex->GetLowEdgeEnergy(i),RindIsGot)>nMax)
				nMax=Rindex->GetValue(Rindex->GetLowEdgeEnergy(i),RindIsGot);
		}
	} else {
		return StepLimit;
	}

        G4double BetaMin = 1./nMax;
        if ( BetaMin >= 1. ) return StepLimit;

        G4double GammaMin = 1./std::sqrt(1.-BetaMin*BetaMin);

        if (gamma < GammaMin ) return StepLimit;

        G4double kinEmin = mass*(GammaMin-1.);

        G4double RangeMin = G4LossTableManager::Instance()->
                                                   GetRange(particleType,
                                                            kinEmin,
                                                            couple);
        G4double Range    = G4LossTableManager::Instance()->
                                                   GetRange(particleType,
                                                            kineticEnergy,
                                                            couple);

        G4double Step = Range - RangeMin;
        if (Step < 1.*um ) return StepLimit;

        if (Step > 0. && Step < StepLimit) StepLimit = Step; 

        // If user has defined an average maximum number of photons to
        // be generated in a Step, then calculate the Step length for
        // that number of photons. 
 
        if (fMaxPhotons > 0) {

           // particle charge
           const G4double charge = aParticle->
                                   GetDefinition()->GetPDGCharge();

	   G4double MeanNumberOfPhotons =0;
	   vector<G4double> MeanNumberOfPh;
	   vector<G4double> dummy; 
	   GetAverageNumberOfPhotons(charge,beta,aMaterial,Rindex,MeanNumberOfPh, dummy);
	   for(size_t i=0;i<MeanNumberOfPh.size();i++)
		   MeanNumberOfPhotons+=MeanNumberOfPh.at(i);

           Step = 0.;
           if (MeanNumberOfPhotons > 0.0) Step = fMaxPhotons /
                                                 MeanNumberOfPhotons;

           if (Step > 0. && Step < StepLimit) StepLimit = Step;
        }

        // If user has defined an maximum allowed change in beta per step
        if (fMaxBetaChange > 0.) {

           G4double dedx = G4LossTableManager::Instance()->
                                                   GetDEDX(particleType,
                                                           kineticEnergy,
                                                           couple);

           G4double deltaGamma = gamma - 
                                 1./std::sqrt(1.-beta*beta*
                                                 (1.-fMaxBetaChange)*
                                                 (1.-fMaxBetaChange));

           Step = mass * deltaGamma / dedx;

           if (Step > 0. && Step < StepLimit) StepLimit = Step;

        }

        *condition = StronglyForced;
        return StepLimit;
}

// GetAverageNumberOfPhotons
// -------------------------
// This routine computes the number of Cerenkov photons produced per
// GEANT-unit (millimeter) in the current medium. 
//             ^^^^^^^^^^

void 
BxCerenkov::GetAverageNumberOfPhotons(const G4double charge,
		const G4double beta, 
		const G4Material* aMaterial,
		G4MaterialPropertyVector* Rindex, vector <G4double> &NumPh, vector <G4double> &range) const
{
	const G4double Rfact = 369.81/(eV * cm);

	if(beta <= 0.0)
		NumPh.push_back(0.0);

	G4double BetaInverse = 1./beta;

	// Vectors used in computation of Cerenkov Angle Integral:
	// 	- Refraction Indices for the current material
	//	- new G4PhysicsOrderedFreeVector allocated to hold CAI's

	G4int materialIndex = aMaterial->GetIndex();

	// Retrieve the Cerenkov Angle Integrals for this material  

	G4PhysicsOrderedFreeVector* CerenkovAngleIntegrals =
		(G4PhysicsOrderedFreeVector*)((*thePhysicsTable)(materialIndex));

	if(!(CerenkovAngleIntegrals->IsFilledVectorExist()))
		NumPh.push_back(0.0);


	G4double dp, ge;
	G4double beta1=beta;
	FindZeros(beta1, Rindex, range);

	if((range.size()%2)!=0){
		BxLog(warning)<<"ERROR IN GENERATING CHERENKOV PHOTONS "<<endlog;
	}else{
		for(size_t i=0;i<range.size();i++){
			dp = range.at(i+1) - range.at(i);

			// Max Cerenkov Angle Integral
			G4double CAImin = CerenkovAngleIntegrals->Value(range.at(i));
			G4double CAImax=CerenkovAngleIntegrals->Value(range.at(i+1));
			ge = CAImax - CAImin;

			if (verboseLevel>0) {
				G4cout << "CAImin = " << CAImin << G4endl;
				G4cout << "ge = " << ge << G4endl;
			}
			// Calculate number of photons 

			G4double NumPhotons = Rfact * charge/eplus * charge/eplus *
				(dp - ge * BetaInverse*BetaInverse);
			NumPh.push_back(NumPhotons);
			i++;
		}
	}

}

void BxCerenkov::FindZeros( const G4double beta,
		G4MaterialPropertyVector* Rindex, vector <G4double> &range)const{
	range.clear();
	G4double BetaInverse = 1./beta;

	G4bool RindIsGot;
	if(Rindex->GetValue(Rindex->GetLowEdgeEnergy(0),RindIsGot)>BetaInverse){
		range.push_back(Rindex->GetLowEdgeEnergy(0));
	}

	for(size_t i=0 ; i<Rindex->GetVectorLength()-1;i++){
		if(((Rindex->GetValue(Rindex->GetLowEdgeEnergy(i),RindIsGot)>=BetaInverse)&&(Rindex->GetValue(Rindex->GetLowEdgeEnergy(i+1),RindIsGot)<=BetaInverse))||((Rindex->GetValue(Rindex->GetLowEdgeEnergy(i),RindIsGot)<=BetaInverse)&&(Rindex->GetValue(Rindex->GetLowEdgeEnergy(i+1),RindIsGot)>=BetaInverse))){
			range.push_back(Rindex->GetLowEdgeEnergy(i));
		}
	}
	if(range.size()%2!=0)
		range.push_back(Rindex->GetLowEdgeEnergy(Rindex->GetVectorLength()));

}

