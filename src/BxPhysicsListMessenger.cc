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
#include "BxPhysicsListMessenger.hh"
#include "BxReadParameters.hh"
#include "BxPhysicsList.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"
#include "BxLogger.hh"
#include "BxOutputVertex.hh"
#include "BxIO.hh"

class BxPhysicsListMessenger;

using namespace std;


BxPhysicsListMessenger::BxPhysicsListMessenger(BxPhysicsList *phys){
  fPhysicsList = phys;
  fDirectory = new G4UIdirectory("/bx/physics/");
  fDirectory->SetGuidance("Control of physical processes");

  fOPCmd = new G4UIcmdWithABool("/bx/physics/optical", this);
  fOPCmd->SetGuidance("Activate or disactivate the optical physics processes");
  fOPCmd->SetGuidance("Default: true");

  fOPcherCmd = new G4UIcmdWithABool("/bx/physics/cherenkov", this);
  fOPcherCmd->SetGuidance("Activate or deactivate the cherenkov physics processes");
  fOPcherCmd->SetGuidance("Default: true");

  fOPscintCmd = new G4UIcmdWithABool("/bx/physics/scintillation", this);
  fOPscintCmd->SetGuidance("Activate or deactivate the scintillation physics processes");
  fOPscintCmd->SetGuidance("Default: false");

  fARCmd = new G4UIcmdWithABool("/bx/physics/scintillator_optical_processes", this);
  fARCmd->SetGuidance("Activate or deactivate absorption and re-emission");
  fARCmd->SetGuidance("Default: true");

  fRayleighCmd = new G4UIcmdWithABool("/bx/physics/rayleigh", this);
  fRayleighCmd->SetGuidance("Activate or deactivate the rayleigh scattering");
  fRayleighCmd->SetGuidance("Default: true");

  fLightYieldCmd = new G4UIcmdWithAnInteger("/bx/physics/lightyield", this);
  fLightYieldCmd->SetGuidance("Light Yield");
  fLightYieldCmd->SetGuidance("Default: 14000");

  fReduceLightYieldCmd = new G4UIcmdWithADouble("/bx/physics/lightyieldscale", this);
  fReduceLightYieldCmd->SetGuidance("Light Yield Reduction Factor (> 1)");
  fReduceLightYieldCmd->SetGuidance("Default: 1.");
 
  fReduceCerenkovLightCmd = new G4UIcmdWithADouble("/bx/physics/cerenkovscale", this);
  fReduceCerenkovLightCmd->SetGuidance("Cerenkov Light  Reduction Factor (> 1)");
  fReduceCerenkovLightCmd->SetGuidance("Default: 1.");

  fMaxNumPhotonCmd = new G4UIcmdWithAnInteger("/bx/physics/max_photons_per_step", this);
  fMaxNumPhotonCmd->SetGuidance("Set the max number of photons per step (it isi an integer!)");
  fMaxNumPhotonCmd->SetGuidance("Default: 300");


  fReemission = new G4UIcmdWithABool("/bx/physics/reemission", this);
  fReemission->SetGuidance("Able or disable reemission process");
  fReemission->SetGuidance("Default: true");
  
  
  fMaxParentIdForQuenching = new G4UIcmdWithAnInteger("/bx/physics/maxparentid", this);

  fHadronicCmd = new G4UIcmdWithAString("/bx/physics/hadronic", this);
  G4String candidates = "standard QGSP_BERT_HP QGSP_BIC_HP FTF_BIC_HP FTF_BERT_HP";
  fHadronicCmd->SetCandidates(candidates);

  fAltHadronicCmd = new G4UIcmdWithABool("/bx/physics/alternative_hadronic_physics", this);
  fAltHadronicCmd->SetGuidance("Turn on/off alternative hadronic physics");
  fAltHadronicCmd->SetGuidance("Default: false");

  fDecayCmd = new G4UIcmdWithABool("/bx/physics/decay", this);
  fDecayCmd->SetGuidance("Turn on/off radioactive decays");
  fDecayCmd->SetGuidance("Default: false");
  
  fCutCmd = new G4UIcmdWithADoubleAndUnit("/bx/physics/cut", this);
  fCutCmd->SetDefaultUnit("mm");
  fCutCmd->SetUnitCandidates("m cm mm");

  fElectronCutCmd = new G4UIcmdWithADoubleAndUnit("/bx/physics/ecut", this);
  fElectronCutCmd->SetDefaultUnit("mm");
  fElectronCutCmd->SetUnitCandidates("m cm mm");
  
  fGammaCutCmd = new G4UIcmdWithADoubleAndUnit("/bx/physics/gammacut", this);
  fGammaCutCmd->SetDefaultUnit("mm");
  fGammaCutCmd->SetUnitCandidates("m cm mm");

  fNewEMCmd = new G4UIcmdWithABool("/bx/physics/new_em", this);

  fNucIntCmd = new G4UIcmdWithABool("/bx/physics/nuclear", this);

  fIonPhysCmd = new G4UIcmdWithABool("/bx/physics/ionphysics", this);

  fTrackSecondariesFirstCmd = new G4UIcmdWithABool("/bx/physics/tracksecondariesfirst", this);
}


BxPhysicsListMessenger::~BxPhysicsListMessenger() {
  delete fOPCmd;
  delete fOPcherCmd;
  delete fOPscintCmd;
  delete fDirectory;
  delete fMaxNumPhotonCmd;
  delete fReduceLightYieldCmd;
  delete fLightYieldCmd;
  delete fReduceCerenkovLightCmd;
  delete fARCmd;  
  delete fRayleighCmd;
  delete fNewEMCmd;
  delete fNucIntCmd;  
  delete fReemission;
  delete fHadronicCmd;
  delete fAltHadronicCmd;
  delete fDecayCmd;
  delete fCutCmd;
  delete fElectronCutCmd;
  delete fGammaCutCmd;
  delete fMaxParentIdForQuenching;
  delete fIonPhysCmd;
  delete fTrackSecondariesFirstCmd;
}


void BxPhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue) { 
  if(command == fOPCmd) {
     fPhysicsList->SetOPIsActivate(fOPCmd->ConvertToBool(newValue));
  } else if(command == fOPcherCmd) {
     fPhysicsList->SetOPcherIsActivate(fOPcherCmd->ConvertToBool(newValue));
  } else if(command == fIonPhysCmd) {
     fPhysicsList->SetIonPhysics(fIonPhysCmd->ConvertToBool(newValue));
     BxLog(routine) << "Ion Physics: " << newValue << endlog;
  } else if(command == fReemission) {
     BxReadParameters::Get()->SetReemission(fReemission->ConvertToBool(newValue));
     BxLog(routine) << "Reemission: " << newValue << endlog;
  } else if(command == fHadronicCmd) {
     if(newValue == "standard") {
       fPhysicsList->SetHadronic(0);
       BxLog(routine) << "Hadronic Physics: " << newValue << endlog;
     } else if(newValue == "QGSP_BERT_HP") {
       fPhysicsList->SetHadronic(1);     
       BxLog(routine) << "Hadronic Physics: " << newValue << endlog;
     } else if(newValue == "QGSP_BIC_HP") {
       fPhysicsList->SetHadronic(2);     
       BxLog(routine) << "Hadronic Physics: " << newValue << endlog;
     } else if(newValue == "FTF_BIC_HP") {
       fPhysicsList->SetHadronic(3);     
       BxLog(routine) << "Hadronic Physics: " << newValue << endlog;
     } else if(newValue == "FTF_BERT_HP") {
       fPhysicsList->SetHadronic(4);     
       BxLog(routine) << "Hadronic Physics: " << newValue << endlog;
     }
  } else if(command == fDecayCmd) {
     BxReadParameters::Get()->SetRDMDecay(fDecayCmd->ConvertToDouble(newValue));
     BxLog(routine) << "Radioactive Decays:  " << newValue << endlog;
  } else if(command == fMaxParentIdForQuenching) {
     BxReadParameters::Get()->SetMaxParentIdForQuenching(fMaxParentIdForQuenching->ConvertToInt(newValue));
     BxLog(routine) << "Quenching: max parent id " << newValue << endlog;
  } else if(command == fAltHadronicCmd) {
    fPhysicsList->SetAltHadIsActivate(fAltHadronicCmd->ConvertToBool(newValue));
  } else if(command == fOPscintCmd) {
    BxLog(routine) << "Scintillation:  " << newValue << endlog;
    fPhysicsList->SetOPscintIsActivate(fOPscintCmd->ConvertToBool(newValue));
  } else if(command == fARCmd) {
    fPhysicsList->SetAbsorptionIsActive(fOPscintCmd->ConvertToBool(newValue));
  } else if(command == fRayleighCmd) {
    fPhysicsList->SetRayleighIsActive(fOPscintCmd->ConvertToBool(newValue));
    BxLog(routine) << "Rayleigh sacttering:  " << newValue << endlog;
  }  else if(command == fMaxNumPhotonCmd) {
    fPhysicsList->SetMaxNumberOfPhotons(fMaxNumPhotonCmd->ConvertToInt(newValue));
  }  else if(command == fCutCmd) {
    fPhysicsList->SetDefCutValue(fCutCmd->ConvertToDimensionedDouble(newValue));
  }  else if(command == fElectronCutCmd) {
    fPhysicsList->SetElectronCutValue(fElectronCutCmd->ConvertToDimensionedDouble(newValue));
  }  else if(command == fGammaCutCmd) {
    fPhysicsList->SetGammaCutValue(fGammaCutCmd->ConvertToDimensionedDouble(newValue));
  } else if(command == fLightYieldCmd) {
    BxReadParameters::Get()->SetLightYield(fLightYieldCmd->ConvertToInt(newValue));
    BxLog(routine) << "Light Yield Re-defined to " << newValue << endlog;
  } else if(command == fReduceLightYieldCmd) {
    BxReadParameters::Get()->SetLightYieldScale(fReduceLightYieldCmd->ConvertToDouble(newValue));
    if(fReduceLightYieldCmd->ConvertToDouble(newValue) < 1.) { 
      BxLog(error) << "Light Yield reduction factor < 1: you are increasing the light yield!!" << endlog;
      BxLog(fatal)  << endlog;    
    }
    BxLog(routine) << "Light Yield prescaling activated: reduction factor " << newValue << endlog;
    BxIO::Get()->GetStreamLogFile() <<  "Light Yield prescaling activated: reduction factor " << newValue << endl; 
  } else if(command == fReduceCerenkovLightCmd) {
    BxReadParameters::Get()->SetLightCerenkovScale(fReduceCerenkovLightCmd->ConvertToDouble(newValue));    
    if(fReduceCerenkovLightCmd->ConvertToDouble(newValue) < 1.)  {
      BxLog(error) << "Cerenkov light reduction factor < 1: you are increasing the Cerenkov light yield!!" << endlog;
      BxLog(fatal) << endl ;
    }
    BxLog(routine) << "Cerenkov light prescaling activated: reduction factor " << newValue << endlog;
    BxIO::Get()->GetStreamLogFile() << "Cerenkov light prescaling activated: reduction factor " << newValue << endl;
  } else if(command == fNewEMCmd) {
    fPhysicsList->SetNewEM(fNewEMCmd->ConvertToBool(newValue));
    if(fNewEMCmd->ConvertToBool(newValue) == false)
      BxLog(routine) << "Old EM Model"  << endlog;
    else 
      BxLog(routine) << "New EM Model"  << endlog;
  } else if(command == fNucIntCmd) {
    fPhysicsList->SetNucInt(fNucIntCmd->ConvertToBool(newValue));
    BxLog(routine) << "Nuclear Production Activation: " << newValue  << endlog;
  } else if(command == fTrackSecondariesFirstCmd) {
    BxOutputVertex::Get()->SetTrackSecondariesFirst(fTrackSecondariesFirstCmd->ConvertToBool(newValue));
    BxLog(routine) << "Track Secondaries First (And Quenching On): " << newValue  << endlog;
  }
}
