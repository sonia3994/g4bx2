//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014

#include "G4ios.hh"

#include "globals.hh"
#include "BxPhysicsList.hh"
#include "BxPhysicsListMessenger.hh"
#include "BxReadParameters.hh"
#include "BxLogger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTable.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"

#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

BxPhysicsList::BxPhysicsList() :  G4VUserPhysicsList()
{ 
 electronCutValue = 5.0*mm;
//  defaultCutValue=0.01*mm;
 gammaCutValue = 1*mm;

 OPIsActivate       = true;
 
 OPcherIsActivate   = true;
 OPscintIsActivate  = true;
 AbsorptionIsActive = true;
 RayleighIsActive   = true;
 fHadronic          = -1 ;
 IsNucInt           = false ;
 IsNewEM            = false ;
 IsIonPhys          = false ;
 fHadronicPhysicsListFlag = true ;

  fMaxNumberOfPhoton = 300;

 BxReadParameters::Get()->SetScattering(true);
 BxReadParameters::Get()->SetReemission(true);
 BxPhys = new BxPhysicsListMessenger(this);

 SetVerboseLevel(1);
}

BxPhysicsList::~BxPhysicsList() {
  delete BxPhys;
}


void BxPhysicsList::ConstructParticle()
{

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  ConstructIons();
  ConstructShortLived();
}

///////////////////////////////////////////////////////////////////////////////
// Construct particles ////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void BxPhysicsList::ConstructBosons() {
  G4BosonConstructor consBos; 
  consBos.ConstructParticle(); 
}

void BxPhysicsList::ConstructLeptons() {
  G4LeptonConstructor consLep; 
  consLep.ConstructParticle();
}

void BxPhysicsList::ConstructMesons() {
  G4MesonConstructor consMes; 
  consMes.ConstructParticle();
}

void BxPhysicsList::ConstructBaryons(){
  G4BaryonConstructor theBaryonConstructor;
  theBaryonConstructor.ConstructParticle();
}

void BxPhysicsList::ConstructIons() {
  G4IonConstructor consIon; 
  consIon.ConstructParticle();  
}

void BxPhysicsList::ConstructShortLived(){
  G4ShortLivedConstructor consShL; 
  consShL.ConstructParticle();
}

//#include "G4DecayPhysics.hh"
//////////////////////////////////////////////////////////////////////////////
// Processes          ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void BxPhysicsList::ConstructProcess() {
  BxLog(routine) << "Default physics cut value = "<< defaultCutValue/mm << " mm"  << endlog;
  AddTransportation();


  if(IsNewEM) ConstructEM2();
  else        ConstructEM();
 
 
  if(GetOPIsActivate())                    ConstructOp();
 
  if(fHadronic == 1)                       ConstructHad_QGSP_BERT_HP();
  //else if(fHadronic == 2)                  ConstructHad_QSGP_BIC_HP();
  //else if(fHadronic == 3)                  ConstructHad_FTF_BIC_HP();
  //else if(fHadronic == 4)                  ConstructHad_FTFP_BERT_HP();
  else if(fHadronic == 0)                  ConstructHad();
  
  if(IsNucInt)                             ConstructNuclearReactions();
  
  if(BxReadParameters::Get()->IsRDMDecay())ConstructGeneral();
  
  if(IsIonPhys)                            ConstructIonPhysics();

}

//////////////////////////////////////////////////////////////////////////////
// EM1                ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCapture.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
// The results can be different for G4hIonisation.hh and
// G4hLowEnergyIonisation.hh
// alpha and GenericIon and deuterons, triton, He3:
//#include "G4hLowEnergyIonisation.hh"
#include "G4EnergyLossTables.hh"
// hLowEnergyIonisation uses Ziegler 1988 as the default
#include "G4hBremsstrahlung.hh"
#include "G4IonParametrisedLossModel.hh"

void BxPhysicsList::ConstructEM()
{
  BxLog(routine) << "EM Physics Active"      << endlog;

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    G4double charge = particle->GetPDGCharge();

    if (particleName == "gamma") {
    // gamma
      // Construct processes for gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
      

    } else if (particleName == "e-") {
    //electron
      G4VProcess* theeminusMultipleScattering = new G4eMultipleScattering();
      G4VProcess* theeminusIonisation         = new G4eIonisation();
      G4VProcess* theeminusBremsstrahlung     = new G4eBremsstrahlung();
      //
      // add processes
      pmanager->AddProcess(theeminusMultipleScattering);
      pmanager->AddProcess(theeminusIonisation);
      pmanager->AddProcess(theeminusBremsstrahlung);
      //      
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxAlongStep,1);
      pmanager->SetProcessOrdering(theeminusIonisation,         idxAlongStep,2);
      pmanager->SetProcessOrdering(theeminusBremsstrahlung,     idxAlongStep,3);
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(theeminusIonisation,         idxPostStep,2);
      pmanager->SetProcessOrdering(theeminusBremsstrahlung,     idxPostStep,3);
     
    } else if (particleName == "e+") {
    //positron
      G4VProcess* theeplusMultipleScattering = new G4eMultipleScattering();
      G4VProcess* theeplusIonisation         = new G4eIonisation();
      G4VProcess* theeplusBremsstrahlung     = new G4eBremsstrahlung();
      G4VProcess* theeplusAnnihilation       = new G4eplusAnnihilation();
      //
      // add processes
      pmanager->AddProcess(theeplusMultipleScattering);
      pmanager->AddProcess(theeplusIonisation);
      pmanager->AddProcess(theeplusBremsstrahlung);
      pmanager->AddProcess(theeplusAnnihilation);
      //
      // set ordering for AtRestDoIt
      pmanager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);
      //
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxAlongStep,1);
      pmanager->SetProcessOrdering(theeplusIonisation,         idxAlongStep,2);
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(theeplusIonisation,         idxPostStep,2);
      pmanager->SetProcessOrdering(theeplusBremsstrahlung,     idxPostStep,3);
      pmanager->SetProcessOrdering(theeplusAnnihilation,       idxPostStep,4);
  

    } else if( particleName == "mu+" ||
               particleName == "mu-"    ) {
    //muon
      G4VProcess* aMultipleScattering = new G4MuMultipleScattering();
      G4VProcess* aBremsstrahlung     = new G4MuBremsstrahlung();
      G4VProcess* aPairProduction     = new G4MuPairProduction();
      G4VProcess* anIonisation        = new G4MuIonisation();
 if( particleName == "mu-" )
          pmanager->AddProcess(new G4MuonMinusCapture(), 0,-1,-1);

      //
      // add processes
      pmanager->AddProcess(anIonisation);
      pmanager->AddProcess(aMultipleScattering);
      pmanager->AddProcess(aBremsstrahlung);
      pmanager->AddProcess(aPairProduction);
      //
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,1);
      pmanager->SetProcessOrdering(anIonisation,        idxAlongStep,2);
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(anIonisation,        idxPostStep,2);
      pmanager->SetProcessOrdering(aBremsstrahlung,     idxPostStep,3);
      pmanager->SetProcessOrdering(aPairProduction,     idxPostStep,4);

 }else if (particleName == "proton" ||
             particleName == "pi+" ||
             particleName == "pi-")
      {
	      //multiple scattering
	      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);

	      //ionisation
	      G4hIonisation* hIonisation = new G4hIonisation();
	      hIonisation->SetStepFunction(0.2, 50*um);
	      pmanager->AddProcess(hIonisation, -1, 2, 2);
	                                                      //bremmstrahlung
	      pmanager->AddProcess(new G4hBremsstrahlung,     -1,-3, 3);
      }
else if(particleName == "alpha"      ||
             particleName == "deuteron"   ||
             particleName == "triton"     ||
             particleName == "He3")
      {
	      //multiple scattering
	      pmanager->AddProcess(new G4hMultipleScattering,-1,1,1);

	      //ionisation
	      G4ionIonisation* ionIoni = new G4ionIonisation();
	      ionIoni->SetStepFunction(0.1, 20*um);
	      pmanager->AddProcess(ionIoni,                   -1, 2, 2);
      }
 else if (particleName == "GenericIon")
 {
	 // OBJECT may be dynamically created as either a GenericIon or nucleus
	 // G4Nucleus exists and therefore has particle type nucleus
	 // genericIon:
	 //multiple scattering
	 pmanager->AddProcess(new G4hMultipleScattering,-1,1,1);

	 //ionisation
	 G4ionIonisation* ionIoni = new G4ionIonisation();
	 ionIoni->SetEmModel(new G4IonParametrisedLossModel());
	 ionIoni->SetStepFunction(0.1, 20*um);
	 pmanager->AddProcess(ionIoni,                   -1, 2, 2);
 }
else if ((!particle->IsShortLived()) &&
             (charge != 0.0) &&
             (particle->GetParticleName() != "chargedgeantino"))
{
	//all others charged particles except geantino
	G4hMultipleScattering* aMultipleScattering = new G4hMultipleScattering();
	G4hIonisation* ahadronIon = new G4hIonisation();

	//multiple scattering
	pmanager->AddProcess(aMultipleScattering,-1,1,1);

	//ionisation
	pmanager->AddProcess(ahadronIon,       -1,2,2);
}









/*
// HERE FINISHES THE PART MODIFIED BY ME    
    } else if (particleName == "proton"     ||
	       particleName == "alpha"      ||
	       particleName == "deuteron"   ||
	       particleName == "triton"     ||
	       particleName == "He3"        ||
	       particleName == "GenericIon" || 
	      (particleName == "nucleus" && charge != 0)) {

      G4hMultipleScattering*	theIonMultipleScattering = new G4hMultipleScattering();
      //G4hLowEnergyIonisation*	theIonIonisation 	 = new G4hLowEnergyIonisation();
      G4ionIonisation* theIonIonisation = new G4ionIonisation();    
      pmanager->AddProcess(theIonMultipleScattering,-1,1,1);
      pmanager->AddProcess(theIonIonisation,-1,2,2);                      
    } else if  ((!particle->IsShortLived()) &&
		    (particle->GetPDGCharge() != 0.0) && 
		    (particle->GetParticleName() != "chargedgeantino")) {
	    // all others charged particles except geantino     
	    G4VProcess* aMultipleScattering = new G4hMultipleScattering();
	    G4VProcess* anIonisation        = new G4hIonisation();
	    //
	    // add processes
	    pmanager->AddProcess(anIonisation);
	    pmanager->AddProcess(aMultipleScattering);
	    //
	    // set ordering for AlongStepDoIt
	    pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,1);
	    pmanager->SetProcessOrdering(anIonisation,        idxAlongStep,2);
	    //
	    // set ordering for PostStepDoIt
	    pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep,1);
	    pmanager->SetProcessOrdering(anIonisation,        idxPostStep,2);
    }*/
  }
}

//////////////////////////////////////////////////////////////////////////////
// EM2                ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include "G4EmStandardPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LossTableManager.hh"
#include "G4EmProcessOptions.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4eMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"


void BxPhysicsList::ConstructEM2()  {
  BxLog(routine) << "EM2 Physics Active"      << endlog;
  
//  G4int verbose = 0;
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    G4double charge = particle->GetPDGCharge();
 

    if (particleName == "gamma") {

      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);

    } else if (particleName == "e-") {

      pmanager->AddProcess(new G4eMultipleScattering,-1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung(),  -1,-3, 3);

    } else if (particleName == "e+") {

      pmanager->AddProcess(new G4eMultipleScattering,-1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1,-3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1, 4);

    } else if (particleName == "mu+" ||
               particleName == "mu-"    ) {

      pmanager->AddProcess(new G4MuMultipleScattering,-1, 1, 1);
      pmanager->AddProcess(new G4MuIonisation,       -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung,   -1,-3, 3);
      pmanager->AddProcess(new G4MuPairProduction,   -1,-4, 4);
           
 
    } else if (particleName == "alpha" ||
               particleName == "He3" ||
               particleName == "GenericIon") {

      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4ionIonisation,       -1, 2, 2);
 //     pmanager->AddProcess(new G4hLowEnergyIonisation  -1, 2, 2);
  
    } else if (particleName == "pi+" ||
               particleName == "pi-" ||
               particleName == "proton" ) {

      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
      pmanager->AddProcess(new G4hBremsstrahlung,     -1,-3, 3);
      pmanager->AddProcess(new G4hPairProduction,     -1,-4, 4);

    } else if (particleName == "B+" ||
	       particleName == "B-" ||
	       particleName == "D+" ||
	       particleName == "D-" ||
	       particleName == "Ds+" ||
	       particleName == "Ds-" ||
               particleName == "anti_lambda_c+" ||
               particleName == "anti_omega-" ||
               particleName == "anti_proton" ||
               particleName == "anti_sigma_c+" ||
               particleName == "anti_sigma_c++" ||
               particleName == "anti_sigma+" ||
               particleName == "anti_sigma-" ||
               particleName == "anti_xi_c+" ||
               particleName == "anti_xi-" ||
               particleName == "deuteron" ||
	       particleName == "kaon+" ||
               particleName == "kaon-" ||
	       particleName == "lambda_c+" ||
               particleName == "omega-" ||
               particleName == "sigma_c+" ||
               particleName == "sigma_c++" ||
               particleName == "sigma+" ||
               particleName == "sigma-" ||
               particleName == "tau+" ||
               particleName == "tau-" ||
               particleName == "triton" ||
               particleName == "xi_c+" ||
               particleName == "xi-" ||
	      (particleName == "nucleus" && charge != 0)) {

      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
    
    } else if  ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {

 	       pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
	       pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
    
    }
  }

 // G4EmProcessOptions opt;
 // opt.SetVerbose(verbose);

 

}

//////////////////////////////////////////////////////////////////////////////
// Optics             ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include "BxCerenkov.hh"
#include "G4OpRayleigh.hh"

#include "G4OpBoundaryProcess.hh"

#include "BxOpAbsorptionReemission.hh"
#include "BxScintillation.hh"


void BxPhysicsList::ConstructOp(){

  BxLog(routine) << "Cherenkov                                     = " << GetOPcherIsActivate()      << endlog;
  BxLog(routine) << "Scintillation                                 = " << GetOPscintIsActivate()    << endlog;  
  BxLog(routine) << "Rayleigh Scattering in Water                  = " <<  RayleighIsActive << endlog;  
  BxLog(routine) << "Complete Absorption and Re-emission Processes = " << AbsorptionIsActive << endlog;  
  BxLog(routine) << "  Sub-process: PC Scattering                  = " <<  BxReadParameters::Get()->GetScattering() << endlog;  
  BxLog(routine) << "  Sub-process: PPO Absorption and Re-emission = " << BxReadParameters::Get()->GetReemission() << endlog;  

  BxCerenkov*     theCerenkovProcess = new BxCerenkov("Cerenkov");
  G4OpRayleigh*   theRayleighScatteringProcess = new G4OpRayleigh();

  BxOpAbsorptionReemission* theAbsorptionProcess = new BxOpAbsorptionReemission();

  G4OpBoundaryProcess* theBoundaryProcess = new G4OpBoundaryProcess();

  BxScintillation* theBetaScintillationProcess = new BxScintillation("Scintillation",11);
  BxScintillation* theMuScintillationProcess = new BxScintillation("Scintillation",13);
  BxScintillation* theProtonScintillationProcess = new BxScintillation("Scintillation",2212);
  BxScintillation* theAlphaScintillationProcess = new BxScintillation("Scintillation",1000020040);

/*
  BxLog(debugging) << "Dumping of the Cherenkov physics table: " << endlog;
  if(BxLogger::GetSeverity() == BxLogger::debugging) theCerenkovProcess->DumpPhysicsTable();
  BxLog(debugging) << "Dumping of the Ryleigh physics table: " << endlog;
  if(BxLogger::GetSeverity() <= BxLogger::debugging) theRayleighScatteringProcess->DumpPhysicsTable();
*/
    
  theAbsorptionProcess->SetVerboseLevel(0);
  theRayleighScatteringProcess->SetVerboseLevel(0);
 
 theBoundaryProcess->SetVerboseLevel(0);
  theCerenkovProcess->SetVerboseLevel(0);
  theCerenkovProcess->SetTrackSecondariesFirst(true);
  theCerenkovProcess->SetMaxNumPhotonsPerStep(fMaxNumberOfPhoton);

  theParticleIterator->reset();
  while((*theParticleIterator)()){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    
    if (GetOPscintIsActivate()){
	    if (particleName == "e+" || particleName == "e-"){
		    pmanager->AddProcess(theBetaScintillationProcess);
		    pmanager->SetProcessOrdering(theBetaScintillationProcess, idxAtRest);
		    pmanager->SetProcessOrdering(theBetaScintillationProcess, idxPostStep);
	    }else if(particleName == "mu+" || particleName == "mu-"){
		    pmanager->AddProcess(theMuScintillationProcess);
		    pmanager->SetProcessOrdering(theMuScintillationProcess, idxAtRest);
		    pmanager->SetProcessOrdering(theMuScintillationProcess, idxPostStep);
	    }else if(particleName == "proton"){
		    pmanager->AddProcess(theProtonScintillationProcess);
		    pmanager->SetProcessOrdering(theProtonScintillationProcess, idxAtRest);
		    pmanager->SetProcessOrdering(theProtonScintillationProcess, idxPostStep);
	    }else if(particleName == "alpha"){
		    pmanager->AddProcess(theAlphaScintillationProcess);
		    pmanager->SetProcessOrdering(theAlphaScintillationProcess, idxAtRest);
		    pmanager->SetProcessOrdering(theAlphaScintillationProcess, idxPostStep);
	    }
    } 

    if (GetOPcherIsActivate()){
	    if (theCerenkovProcess->IsApplicable(*particle)) {
		    pmanager->AddProcess(theCerenkovProcess);
		    pmanager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
	    }
    }

    if (particleName == "opticalphoton"){
	    BxLog(trace) << " AddDiscreteProcess to OpticalPhoton " << endlog;
	    if(AbsorptionIsActive) pmanager->AddDiscreteProcess(theAbsorptionProcess);
	    if(RayleighIsActive)   pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
	    pmanager->AddDiscreteProcess(theBoundaryProcess);
    }
  }
}


//////////////////////////////////////////////////////////////////////////////
// Nuclear Reactions  ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


#include "G4EmExtraPhysics.hh"
void BxPhysicsList::ConstructNuclearReactions() {   
   BxLog(routine) << "Nuclear EM Physics Active" << endlog;
   G4EmExtraPhysics *p = new G4EmExtraPhysics;
   G4String on = "on";
   p->GammaNuclear (on);
   p->MuonNuclear (on);
   p->Synch(on);
   p->ConstructProcess();
  
}
//////////////////////////////////////////////////////////////////////////////
// Ion processes ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include "G4IonPhysics.hh"
void BxPhysicsList::ConstructIonPhysics() {   
   BxLog(routine) << "Ion Physics Active" << endlog;
   G4IonPhysics *p = new G4IonPhysics("ion");
   p->ConstructProcess();
  
}

//////////////////////////////////////////////////////////////////////////////
// Hadronic processes ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//QGSP BERTINI cascade NEUTRONHP model
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
void BxPhysicsList::ConstructHad_QGSP_BERT_HP() {
  BxLog(routine) << "Hadronic Physics Active: QSGP_BERT_HP model" << endlog;
  G4HadronPhysicsQGSP_BERT_HP *hadPhysicsList = new G4HadronPhysicsQGSP_BERT_HP;
  hadPhysicsList->ConstructProcess();
}

/*
//QGSP BINARY cascade NEUTRONHP model
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
void BxPhysicsList::ConstructHad_QSGP_BIC_HP() {
  BxLog(routine) << "Hadronic Physics Active: QSGP_BIC_HP model" << endlog;
  G4HadronPhysicsQGSP_BIC_HP *hadPhysicsList = new G4HadronPhysicsQGSP_BIC_HP;
  hadPhysicsList->ConstructProcess();
}

//FTF BINARY cascade NEUTRONHP model
#include "G4HadronPhysicsFTF_BIC.hh"
void BxPhysicsList::ConstructHad_FTF_BIC_HP() {
  BxLog(routine) << "Hadronic Physics Active: FTF_BIC_HP model" << endlog;
  G4HadronPhysicsFTF_BIC *hadPhysicsList = new G4HadronPhysicsFTF_BIC;
  hadPhysicsList->ConstructProcess();
}

//FTFP BERTINI cascade NEUTRONHP model
#include "G4HadronPhysicsFTFP_BERT.hh"
void BxPhysicsList::ConstructHad_FTFP_BERT_HP() {
  BxLog(routine) << "Hadronic Physics Active: FTFP_BERT_HP model" << endlog;
  G4HadronPhysicsFTFP_BERT *hadPhysicsList = new G4HadronPhysicsFTFP_BERT;
  hadPhysicsList->ConstructProcess();
}
*/
// Hadronic processes ////////////////////////////////////////////////////////
// Taken from Geant4 10.0 examples: Advanced Underground physics
// Elastic processes:
#include "G4HadronElasticProcess.hh"
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"

// Inelastic processes:
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

// added for revision (it is needed in the alternative proton)
#include "G4ProtonInelasticCrossSection.hh"

// Low-energy Models: < 20GeV
 /*#include "G4LElastic.hh"
 #include "G4LEPionPlusInelastic.hh"
 #include "G4LEPionMinusInelastic.hh"
 #include "G4LEKaonPlusInelastic.hh"
 #include "G4LEKaonZeroSInelastic.hh"
 #include "G4LEKaonZeroLInelastic.hh"
 #include "G4LEKaonMinusInelastic.hh"
 #include "G4LEProtonInelastic.hh"
 #include "G4LEAntiProtonInelastic.hh"
 #include "G4LENeutronInelastic.hh"
 #include "G4LEAntiNeutronInelastic.hh"
 #include "G4LEDeuteronInelastic.hh"
 #include "G4LETritonInelastic.hh"
 #include "G4LEAlphaInelastic.hh"
*/

// High energy FTFP model and Bertini cascade
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4CascadeInterface.hh"

// Cross sections
#include "G4VCrossSectionDataSet.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4CrossSectionElastic.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4AntiNuclElastic.hh"

#include "G4CrossSectionInelastic.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4GGNuclNuclCrossSection.hh"

#include "G4HadronElastic.hh"
#include "G4HadronCaptureProcess.hh"

// Neutron high-precision models: <20 MeV
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4LFission.hh"

// Stopping processes
#include "G4PiMinusAbsorptionBertini.hh"
#include "G4KaonMinusAbsorptionBertini.hh"
#include "G4AntiProtonAbsorptionFritiof.hh"

void BxPhysicsList::ConstructHad() {
  BxLog(routine) << "Hadronic Physics Active" << endlog;

//Elastic models
  const G4double elastic_elimitPi = 1.0*GeV;

  G4HadronElastic* elastic_lhep0 = new G4HadronElastic();
  G4HadronElastic* elastic_lhep1 = new G4HadronElastic();
  elastic_lhep1->SetMaxEnergy( elastic_elimitPi );
  G4ChipsElasticModel* elastic_chip = new G4ChipsElasticModel();
  G4ElasticHadrNucleusHE* elastic_he = new G4ElasticHadrNucleusHE();
  elastic_he->SetMinEnergy( elastic_elimitPi );


  // Inelastic scattering
  const G4double theFTFMin0 =    0.0*GeV;
  const G4double theFTFMin1 =    4.0*GeV;
  const G4double theFTFMax =   100.0*TeV;
  const G4double theBERTMin0 =   0.0*GeV;
  const G4double theBERTMin1 =  19.0*MeV;
  const G4double theBERTMax =    5.0*GeV;
  const G4double theHPMin =      0.0*GeV;
  const G4double theHPMax =     20.0*MeV;

  G4FTFModel * theStringModel = new G4FTFModel;
  G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay( new G4LundStringFragmentation );
  theStringModel->SetFragmentationModel( theStringDecay );
  G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel( new G4ExcitationHandler );
  G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface( thePreEquilib );

  G4TheoFSGenerator * theFTFModel0 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel0->SetHighEnergyGenerator( theStringModel );
  theFTFModel0->SetTransport( theCascade );
  theFTFModel0->SetMinEnergy( theFTFMin0 );
  theFTFModel0->SetMaxEnergy( theFTFMax );

  G4TheoFSGenerator * theFTFModel1 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel1->SetHighEnergyGenerator( theStringModel );
  theFTFModel1->SetTransport( theCascade );
  theFTFModel1->SetMinEnergy( theFTFMin1 );
  theFTFModel1->SetMaxEnergy( theFTFMax );

  G4CascadeInterface * theBERTModel0 = new G4CascadeInterface;
  theBERTModel0->SetMinEnergy( theBERTMin0 );
  theBERTModel0->SetMaxEnergy( theBERTMax );

  G4CascadeInterface * theBERTModel1 = new G4CascadeInterface;
  theBERTModel1->SetMinEnergy( theBERTMin1 );
  theBERTModel1->SetMaxEnergy( theBERTMax );

  G4VCrossSectionDataSet * thePiData = new G4CrossSectionPairGG( new G4PiNuclearCrossSection, 91*GeV );
  G4VCrossSectionDataSet * theAntiNucleonData = new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );
  G4VCrossSectionDataSet * theGGNuclNuclData = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4GGNuclNuclCrossSection::Default_Name());

  theParticleIterator->reset();
  while ((*theParticleIterator)())
  { 
	  G4ParticleDefinition* particle = theParticleIterator->value();
	  G4ProcessManager* pmanager = particle->GetProcessManager();
	  G4String particleName = particle->GetParticleName();

	  if (particleName == "pi+")
	  { 
		  //Elastic scattering
			  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
		  theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
		  theElasticProcess->RegisterMe( elastic_lhep1 );
		  theElasticProcess->RegisterMe( elastic_he );
		  pmanager->AddDiscreteProcess( theElasticProcess );
		  //Inelastic scattering
		  G4PionPlusInelasticProcess* theInelasticProcess =
			  new G4PionPlusInelasticProcess("inelastic");
		  theInelasticProcess->AddDataSet( thePiData );
		  theInelasticProcess->RegisterMe( theFTFModel1 );
		  theInelasticProcess->RegisterMe( theBERTModel0 );
		  pmanager->AddDiscreteProcess( theInelasticProcess );
	  }

	  else if (particleName == "pi-")
	  {
		  //Elastic scattering
			  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
		  theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
		  theElasticProcess->RegisterMe( elastic_lhep1 );
		  theElasticProcess->RegisterMe( elastic_he );
		  pmanager->AddDiscreteProcess( theElasticProcess );
		  //Inelastic scattering
		  G4PionMinusInelasticProcess* theInelasticProcess =
			  new G4PionMinusInelasticProcess("inelastic");
		  theInelasticProcess->AddDataSet( thePiData );
		  theInelasticProcess->RegisterMe( theFTFModel1 );
		  theInelasticProcess->RegisterMe( theBERTModel0 );
		  pmanager->AddDiscreteProcess( theInelasticProcess );
		  //Absorption
		  pmanager->AddRestProcess(new G4PiMinusAbsorptionBertini, ordDefault);
	  }

	  else if (particleName == "kaon+")
	  {
		  // Elastic scattering
		  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
		  theElasticProcess->RegisterMe( elastic_lhep0 );
		  pmanager->AddDiscreteProcess( theElasticProcess );
		  // Inelastic scattering       
		  G4KaonPlusInelasticProcess* theInelasticProcess =
			  new G4KaonPlusInelasticProcess("inelastic");
		  theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
				  GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name()));
		  theInelasticProcess->RegisterMe( theFTFModel1 );
		  theInelasticProcess->RegisterMe( theBERTModel0 );
		  pmanager->AddDiscreteProcess( theInelasticProcess );
	  }

	  else if (particleName == "kaon0S")
	  {
		  // Elastic scattering
		  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
		  theElasticProcess->RegisterMe( elastic_lhep0 );
		  pmanager->AddDiscreteProcess( theElasticProcess );
		  // Inelastic scattering   
		  G4KaonZeroSInelasticProcess* theInelasticProcess =
			  new G4KaonZeroSInelasticProcess("inelastic");
		  theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
				  GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
		  theInelasticProcess->RegisterMe( theFTFModel1 );
		  theInelasticProcess->RegisterMe( theBERTModel0 );
		  pmanager->AddDiscreteProcess( theInelasticProcess );
	  }

	  else if (particleName == "kaon0L")
	  { 
		  // Elastic scattering
		  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
		  theElasticProcess->RegisterMe( elastic_lhep0 );
		  pmanager->AddDiscreteProcess( theElasticProcess );
		  // Inelastic scattering
		  G4KaonZeroLInelasticProcess* theInelasticProcess =
			  new G4KaonZeroLInelasticProcess("inelastic");
		  theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
				  GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
		  theInelasticProcess->RegisterMe( theFTFModel1 );
		  theInelasticProcess->RegisterMe( theBERTModel0 );
		  pmanager->AddDiscreteProcess( theInelasticProcess );
	  }

	  else if (particleName == "kaon-")
	  {
		  // Elastic scattering
		  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
		  theElasticProcess->RegisterMe( elastic_lhep0 );
		  pmanager->AddDiscreteProcess( theElasticProcess );
		  // Inelastic scattering
		  G4KaonMinusInelasticProcess* theInelasticProcess =
			  new G4KaonMinusInelasticProcess("inelastic");
		  theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
				  GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name()));
		  theInelasticProcess->RegisterMe( theFTFModel1 );
		  theInelasticProcess->RegisterMe( theBERTModel0 );
		  pmanager->AddDiscreteProcess( theInelasticProcess );
		  pmanager->AddRestProcess(new G4KaonMinusAbsorptionBertini, ordDefault);
	  }

	  else if (particleName == "proton")
	  {

		  //Elastic scattering
			  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
		  theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->
				  GetCrossSectionDataSet(G4ChipsProtonElasticXS::Default_Name()));
		  theElasticProcess->RegisterMe( elastic_chip );
		  pmanager->AddDiscreteProcess( theElasticProcess );
		  // Inelastic scattering
		  G4ProtonInelasticProcess* theInelasticProcess =
			  new G4ProtonInelasticProcess("inelastic");
		  theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Proton::Proton() ) );
		  theInelasticProcess->RegisterMe( theFTFModel1 );
		  theInelasticProcess->RegisterMe( theBERTModel0 );
		  pmanager->AddDiscreteProcess( theInelasticProcess );
/*
		  G4ProtonInelasticProcess* theProtonInelasticProcess = new G4ProtonInelasticProcess("pInelastic");
		  G4CascadeInterface* theBertiniModel = new G4CascadeInterface();
		  theBertiniModel->SetMinEnergy(0.0*MeV);
		  theBertiniModel->SetMaxEnergy(9.9*GeV);
		  theProtonInelasticProcess->RegisterMe(theBertiniModel);
*/
	  }
	  else if (particleName == "anti_proton")
	  {
		  // Elastic scattering
		  const G4double elastic_elimitAntiNuc = 100.0*CLHEP::MeV;
		  G4AntiNuclElastic* elastic_anuc = new G4AntiNuclElastic();
		  elastic_anuc->SetMinEnergy( elastic_elimitAntiNuc );
		  G4CrossSectionElastic* elastic_anucxs = new G4CrossSectionElastic( elastic_anuc->GetComponentCrossSection() );
		  G4HadronElastic* elastic_lhep2 = new G4HadronElastic();
		  elastic_lhep2->SetMaxEnergy( elastic_elimitAntiNuc );
		  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
		  theElasticProcess->AddDataSet( elastic_anucxs );
		  theElasticProcess->RegisterMe( elastic_lhep2 );
		  theElasticProcess->RegisterMe( elastic_anuc );
		  pmanager->AddDiscreteProcess( theElasticProcess );
		  // Inelastic scattering
		  G4AntiProtonInelasticProcess* theInelasticProcess =
			  new G4AntiProtonInelasticProcess("inelastic");
		  theInelasticProcess->AddDataSet( theAntiNucleonData );
		  theInelasticProcess->RegisterMe( theFTFModel0 );
		  pmanager->AddDiscreteProcess( theInelasticProcess );
		  // Absorption
		  pmanager->AddRestProcess(new G4AntiProtonAbsorptionFritiof, ordDefault);
	  }

	  else if (particleName == "neutron") {
//G4HadronElasticProcess* theNeutronElasticProcess = new G4HadronElasticProcess("nElastic");
//G4NeutronInelasticProcess* theNeutronInelasticProcess = new G4NeutronInelasticProcess("nInelastic");
G4HadronFissionProcess* theNeutronFissionProcess = new G4HadronFissionProcess();
//G4HadronCaptureProcess* theNeutronCaptureProcess = new G4HadronCaptureProcess();

		  // elastic scattering
		  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
		  theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsNeutronElasticXS::Default_Name()));
		  G4HadronElastic* elastic_neutronChipsModel = new G4ChipsElasticModel();
		  elastic_neutronChipsModel->SetMinEnergy( 19.0*CLHEP::MeV );
		  theElasticProcess->RegisterMe( elastic_neutronChipsModel );
		  G4NeutronHPElastic * theElasticNeutronHP = new G4NeutronHPElastic;
		  theElasticNeutronHP->SetMinEnergy( theHPMin );
		  theElasticNeutronHP->SetMaxEnergy( theHPMax );
		  theElasticProcess->RegisterMe( theElasticNeutronHP );
		  theElasticProcess->AddDataSet( new G4NeutronHPElasticData );
		  pmanager->AddDiscreteProcess( theElasticProcess );
		  // inelastic scattering         
		  G4NeutronInelasticProcess* theInelasticProcess =
			  new G4NeutronInelasticProcess("inelastic");
		  theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Neutron::Neutron() ) );
		  theInelasticProcess->RegisterMe( theFTFModel1 );
		  theInelasticProcess->RegisterMe( theBERTModel1 );
		  G4NeutronHPInelastic * theNeutronInelasticHPModel = new G4NeutronHPInelastic;
		  theNeutronInelasticHPModel->SetMinEnergy( theHPMin );
		  theNeutronInelasticHPModel->SetMaxEnergy( theHPMax );
		  theInelasticProcess->RegisterMe( theNeutronInelasticHPModel );
		  theInelasticProcess->AddDataSet( new G4NeutronHPInelasticData );
		  pmanager->AddDiscreteProcess(theInelasticProcess);
		  // capture
		  G4HadronCaptureProcess* theCaptureProcess =
			  new G4HadronCaptureProcess;
		  G4NeutronHPCapture * theLENeutronCaptureModel = new G4NeutronHPCapture;
		  theLENeutronCaptureModel->SetMinEnergy(theHPMin);
		  theLENeutronCaptureModel->SetMaxEnergy(theHPMax);
		  theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
		  theCaptureProcess->AddDataSet( new G4NeutronHPCaptureData);
		  pmanager->AddDiscreteProcess(theCaptureProcess);

		  //fission
		  G4LFission* theNeutronFissionModel = new G4LFission();
		  theNeutronFissionModel->SetMinEnergy(20.0*MeV);
		  theNeutronFissionModel->SetMaxEnergy(20.0*TeV);
		  theNeutronFissionProcess->RegisterMe(theNeutronFissionModel);
/*
		  G4NeutronHPElastic* theHPElastic = new G4NeutronHPElastic();
		  G4NeutronHPElasticData* theHPElasticData = new G4NeutronHPElasticData();
		  theNeutronElasticProcess->AddDataSet(theHPElasticData);
		  theNeutronElasticProcess->RegisterMe(theHPElastic);

		  G4NeutronHPInelastic* theHPInelastic = new G4NeutronHPInelastic();
		  G4NeutronHPInelasticData* theHPInelasticData = new G4NeutronHPInelasticData();
		  theNeutronInelasticProcess->AddDataSet(theHPInelasticData);
		  theNeutronInelasticProcess->RegisterMe(theHPInelastic);

		  G4NeutronHPFission* theHPFission = new G4NeutronHPFission();
		  G4NeutronHPFissionData* theHPFissionData = new G4NeutronHPFissionData();
		  theNeutronFissionProcess->AddDataSet(theHPFissionData);
		  theNeutronFissionProcess->RegisterMe(theHPFission);

		  G4NeutronHPCapture* theHPCapture = new G4NeutronHPCapture();
		  G4NeutronHPCaptureData* theHPCaptureData = new G4NeutronHPCaptureData();
		  theNeutronCaptureProcess->AddDataSet(theHPCaptureData);
		  theNeutronCaptureProcess->RegisterMe(theHPCapture);
		  pmanager->AddDiscreteProcess(theNeutronElasticProcess);
		  pmanager->AddDiscreteProcess(theNeutronInelasticProcess);
		  pmanager->AddDiscreteProcess(theNeutronFissionProcess);
		  pmanager->AddDiscreteProcess(theNeutronCaptureProcess);

*/
	  }
	  else if (particleName == "anti_neutron")
	  {

		  // Elastic scattering
		  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
		  theElasticProcess->RegisterMe( elastic_lhep0 );
		  pmanager->AddDiscreteProcess( theElasticProcess );
		  // Inelastic scattering (include annihilation on-fly)
		  G4AntiNeutronInelasticProcess* theInelasticProcess =
			  new G4AntiNeutronInelasticProcess("inelastic");
		  theInelasticProcess->AddDataSet( theAntiNucleonData );
		  theInelasticProcess->RegisterMe( theFTFModel0 );
		  pmanager->AddDiscreteProcess( theInelasticProcess );
	  }

	  else if (particleName == "deuteron")
	  { 
		  // Elastic scattering
		  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
		  theElasticProcess->RegisterMe( elastic_lhep0 );
		  pmanager->AddDiscreteProcess( theElasticProcess );
		  // Inelastic scattering
		  G4DeuteronInelasticProcess* theInelasticProcess =
			  new G4DeuteronInelasticProcess("inelastic");
		  theInelasticProcess->AddDataSet( theGGNuclNuclData );
		  theInelasticProcess->RegisterMe( theFTFModel1 );
		  theInelasticProcess->RegisterMe( theBERTModel0 );
		  pmanager->AddDiscreteProcess( theInelasticProcess );
	  }

	  else if (particleName == "triton")
	  { 
		  // Elastic scattering
		  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
		  theElasticProcess->RegisterMe( elastic_lhep0 );
		  pmanager->AddDiscreteProcess( theElasticProcess );
		  // Inelastic scattering
		  G4TritonInelasticProcess* theInelasticProcess =
			  new G4TritonInelasticProcess("inelastic");
		  theInelasticProcess->AddDataSet( theGGNuclNuclData );
		  theInelasticProcess->RegisterMe( theFTFModel1 );
		  theInelasticProcess->RegisterMe( theBERTModel0 );
		  pmanager->AddDiscreteProcess( theInelasticProcess );
	  }
	  else if (particleName == "alpha")
	  {
		  // Elastic scattering
		  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
		  theElasticProcess->RegisterMe( elastic_lhep0 );
		  pmanager->AddDiscreteProcess( theElasticProcess );
		  // Inelastic scattering
		  G4AlphaInelasticProcess* theInelasticProcess =
			  new G4AlphaInelasticProcess("inelastic");
		  theInelasticProcess->AddDataSet( theGGNuclNuclData );
		  theInelasticProcess->RegisterMe( theFTFModel1 );
		  theInelasticProcess->RegisterMe( theBERTModel0 );
		  pmanager->AddDiscreteProcess( theInelasticProcess );
	  }

  }
}

// Decays ///////////////////////////////////////////////////////////////////
#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
void BxPhysicsList::ConstructGeneral()  {
  BxLog(routine) << "Radioactive Decay Physics Active" << endlog;

  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager->AddProcess(theDecayProcess);
      pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
  
  const G4IonTable *theIonTable = G4ParticleTable::GetParticleTable()->GetIonTable();
  G4RadioactiveDecay *theRadioactiveDecay = new G4RadioactiveDecay();
  theRadioactiveDecay->SetICM(true);                //Internal Conversion
  theRadioactiveDecay->SetARM(true);               //Atomic Rearangement
  for (G4int i=0; i<theIonTable->Entries(); i++) {
	  G4String particleName = theIonTable->GetParticle(i)->GetParticleName();
	  if (particleName == "GenericIon") {
		  G4ProcessManager* pmanager = theIonTable->GetParticle(i)->GetProcessManager();
		  pmanager->SetVerboseLevel(0);
		  pmanager->AddProcess(theRadioactiveDecay);
		  pmanager->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
		  pmanager->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
	  } 
  }
}

void BxPhysicsList::SetCuts()
{

  //special for low energy physics
  //G4double lowlimit=250*eV;  
  //G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit,100.*GeV);

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  // 
  //double cut_em = defaultCutValue;
  G4double cut_had= defaultCutValue;
  G4double cut_ion= 1000*mm;

  SetCutValue(gammaCutValue, "gamma");
  SetCutValue(electronCutValue, "e-");
  SetCutValue(electronCutValue, "e+");

  //SetCutValue(10*keV, "proton");
  SetCutValue(0.6*mm, "proton");
  SetCutValue(cut_had, "anti_proton");
  SetCutValue(0.6*mm, "neutron");
  
  SetCutValue(cut_ion, "alpha");
  SetCutValue(cut_ion, "GenericIon");

   
}
