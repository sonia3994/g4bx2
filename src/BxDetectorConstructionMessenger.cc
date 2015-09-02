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
#include "BxDetectorConstructionMessenger.hh"
#include "BxDetectorConstruction.hh"
#include "BxPropertyCollection.hh"
#include "BxMaterial.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIdirectory.hh"
#include "BxLogger.hh"
#include "BxIO.hh"
#include "BxReadParameters.hh"
#include "G4SystemOfUnits.hh"
class BxDetectorConstructionMessenger;

using namespace std;


BxDetectorConstructionMessenger::BxDetectorConstructionMessenger(BxDetectorConstruction *det){
  fDetectorConstruction = det;
  fDirectory = new G4UIdirectory("/bx/detector/");
  fDirectory->SetGuidance("Controls detector gemetry and materials");

  fTankCmd = new G4UIcmdWithABool("/bx/detector/tank", this);
  fTankCmd->SetGuidance("Able or disable the tank");
  fTankCmd->SetGuidance("Default: true");

/*  fRockCmd = new G4UIcmdWithABool("/bx/detector/rock", this);
  fRockCmd->SetGuidance("Build rock and concrete");
  fRockCmd->SetGuidance("Default: true");
*/  
  fConfCmd = new G4UIcmdWithAString("/bx/detector/configuration", this);
  fConfCmd->SetGuidance("Choose the detector configuration");
  fConfCmd->SetGuidance("Default: inScint");
  fConfCmd->SetCandidates("full fullWater inWater inScint outer ctf");
/*
  fOPERACmd = new G4UIcmdWithABool("/bx/detector/opera", this);
  fOPERACmd->SetGuidance("Activate or disactivate OPERA detector skeleton");
  fOPERACmd->SetGuidance("Default: false");
*/
  fNoPMTCmd = new G4UIcmdWithABool("/bx/detector/nopmts", this);
  fNoPMTCmd->SetGuidance("Default: false");

  fPMTCmd = new G4UIcmdWithAString("/bx/detector/PMTconfiguration", this);
  fPMTCmd->SetGuidance("Choose the PMT configuration");
  fPMTCmd->SetGuidance("Default: simplePMT");
  fPMTCmd->SetCandidates("fullPMT diskPMT diskPMTwithConcentrator simplePMT");
  
  fNoConcCmd = new G4UIcmdWithABool("/bx/detector/noConcentrators", this);
  fNoConcCmd->SetGuidance("True= no Concentrators at all in any PMT configuration");
  fNoConcCmd->SetGuidance("Default: false");

  fPMTdistribCmd = new G4UIcmdWithAString("/bx/detector/PMTdistribution", this);
  fPMTdistribCmd->SetGuidance("Choose the PMT distribution on the sphere");
  fPMTdistribCmd->SetGuidance("Default: real");
  fPMTdistribCmd->SetCandidates("real uniform");
  
  fPMTnumCmd = new G4UIcmdWithAnInteger("/bx/detector/PMTnumber", this);
  fPMTnumCmd->SetGuidance("Choose the number of PMTs uniformely distributed on the sphere");
  fPMTnumCmd->SetGuidance("Default: 1000");
  
  fSSSReflectivityCmd = new G4UIcmdWithADouble("/bx/detector/SSSReflectivity", this);
  fSSSReflectivityCmd->SetGuidance("Default: 0.68");
  
  fSSSspecLOBECmd = new G4UIcmdWithADouble("/bx/detector/SSSspecLOBE", this);
  fSSSspecLOBECmd->SetGuidance("Default: 0.0");
  
  fSSSspecSPIKECmd = new G4UIcmdWithADouble("/bx/detector/SSSspecSPIKE", this);
  fSSSspecSPIKECmd->SetGuidance("Default: 0.02");
  
  fSSSbackSCATTCmd = new G4UIcmdWithADouble("/bx/detector/SSSbackSCATT", this);
  fSSSbackSCATTCmd->SetGuidance("Default: 0.0");
  
  fConcentratorExternalReflectivityCmd = new G4UIcmdWithADouble("/bx/detector/externalConcentratorReflectivity", this);
  fConcentratorExternalReflectivityCmd->SetGuidance("Default: 0.94");
  
  fConcentratorInternalReflectivityCmd = new G4UIcmdWithADouble("/bx/detector/internalConcentratorReflectivity", this);
  fConcentratorInternalReflectivityCmd->SetGuidance("Default: 0.95");
  
  fConcentratorExternalspecSPIKECmd = new G4UIcmdWithADouble("/bx/detector/externalConcentratorSpike", this);
  fConcentratorExternalspecSPIKECmd->SetGuidance("Default: 0.9");
  
  fConcentratorInternalspecSPIKECmd = new G4UIcmdWithADouble("/bx/detector/internalConcentratorSpike", this);
  fConcentratorInternalspecSPIKECmd->SetGuidance("Default: 0.98");
  
  fPMTReflectivityCmd = new G4UIcmdWithADouble("/bx/detector/cathodeReflectivity", this);
  fPMTReflectivityCmd->SetGuidance("Default: 0.0");
  
  fPMTringReflectivityCmd = new G4UIcmdWithADouble("/bx/detector/PMTringReflectivity", this);
  fPMTringReflectivityCmd->SetGuidance("Default: 0.76");
  
  fPMTringspecSPIKECmd = new G4UIcmdWithADouble("/bx/detector/PMTringspecSPIKE", this);
  fPMTringspecSPIKECmd->SetGuidance("Default: 0.9");
  
  fShieldReflectivityCmd = new G4UIcmdWithADouble("/bx/detector/shieldReflectivity", this);
  fShieldReflectivityCmd->SetGuidance("Default: 0.4");
  
  fBirksACmd = new G4UIcmdWithADouble("/bx/detector/birksAlpha", this);
  fBirksACmd->SetGuidance("Default: 0.0075 cm/MeV");

  fBirksA2Cmd = new G4UIcmdWithADouble("/bx/detector/birksAlpha2", this);
  fBirksA2Cmd->SetGuidance("Default: 0.0 (cm/MeV)**2");
 
  fBirksBCmd = new G4UIcmdWithADouble("/bx/detector/birksBeta", this);
  fBirksBCmd->SetGuidance("Default: 0.0144 cm/MeV");

  fBirksB2Cmd = new G4UIcmdWithADouble("/bx/detector/birksBeta2", this);
  fBirksB2Cmd->SetGuidance("Default: 0.00 (cm/MeV)**2");

  fVesselOriginCmd = new G4UIcmdWith3VectorAndUnit("/bx/detector/vesselorigin", this);
  fVesselOriginCmd->SetGuidance("Default: 0. 0. 0. cm");
  fVesselOriginCmd->SetDefaultUnit("cm");
  fVesselOriginCmd->SetUnitCandidates("m cm mm");

  fRZDimCmd= new G4UIcmdWithADoubleAndUnit("/bx/detector/vessel_step",this);
  fRZDimCmd->SetGuidance("It defines the length of the step in the vessel's polycone. Default: ~1cm.");
  fRZDimCmd->SetDefaultUnit("cm");
  fRZDimCmd->SetUnitCandidates("m cm mm");

  fSourceDirectory = new G4UIdirectory("/bx/detector/source/");
  fSourceDirectory->SetGuidance("Control source gemetry and materials");
  
  fSourceTypeCmd   = new G4UIcmdWithAString("/bx/detector/source/type", this);	    
  fSourceTypeCmd->SetCandidates("sphere cylinder box Am-Be");
  
  fSourceOriginCmd = new G4UIcmdWith3VectorAndUnit("/bx/detector/source/origin", this);     
  fSourceOriginCmd->SetUnitCandidates("m cm mm");
  
  fSourceRadiusCmd = new G4UIcmdWithADoubleAndUnit("/bx/detector/source/radius", this);     
  fSourceRadiusCmd->SetUnitCandidates("m cm mm");
  
  fSourceLongCmd = new G4UIcmdWithADoubleAndUnit("/bx/detector/source/length", this);     
  fSourceLongCmd->SetUnitCandidates("m cm mm");
  
  fSourceThickCmd = new G4UIcmdWithADoubleAndUnit("/bx/detector/source/thickness", this);     
  fSourceThickCmd->SetUnitCandidates("m cm mm");
  
  fSourceXCmd = new  G4UIcmdWithADoubleAndUnit("/bx/detector/source/x", this);     
  fSourceXCmd->SetUnitCandidates("m cm mm");
  
  fSourceYCmd = new  G4UIcmdWithADoubleAndUnit("/bx/detector/source/y", this);     
  fSourceYCmd->SetUnitCandidates("m cm mm");
  
  fSourceZCmd = new  G4UIcmdWithADoubleAndUnit("/bx/detector/source/z", this);     
  fSourceZCmd->SetUnitCandidates("m cm mm");
  
  fSourceMaterialCmd  = new G4UIcmdWithAString("/bx/detector/source/material", this);	    
  fSourceMaterialCmd->SetCandidates("pc scintillator water nylon quartz buffer air (steel vacuum)");

  fSourceVialMaterialCmd  = new G4UIcmdWithAString("/bx/detector/source/vialmaterial", this);	    
  fSourceVialMaterialCmd->SetCandidates("pc scintillator water nylon quartz buffer air (steel vacuum) ");

  fCheckOverlapSourceCmd = new G4UIcmdWithABool("/bx/detector/check_source_overlaps", this);
  fCheckOverlapSourceCmd->SetGuidance("Default: false");
  
  fCheckOverlapCmd = new G4UIcmdWithABool("/bx/detector/check_overlaps", this);
  fCheckOverlapCmd->SetGuidance("Default: false");

  fUpdateGeometryCmd = new G4UIcmdWithoutParameter("/bx/detector/update",this); 

  //fIsSensitiveCmd  = new G4UIcmdWithABool("/bx/detector/sensitive", this);

  fDefVesselShapeCmd = new G4UIcmdWithABool("/bx/detector/deformed_vessel", this);
  fDefVesselShapeCmd->SetGuidance("Default: false");

  fRunNumberCmd= new G4UIcmdWithAnInteger("/bx/detector/run_number",this);
  fRunNumberCmd->SetGuidance("Default: 0, i.e. no run number specified.");

  fIsRingCmd    = new G4UIcmdWithABool("/bx/detector/setEndCaps", this);
  
  fScintillatorDirectory = new G4UIdirectory("/bx/scintillator/");
  fScintillatorDirectory->SetGuidance("Controls scintillator's properties (time constants and exponential weights)");

  fAlphaDecayCmd  = new G4UIcmdWithAString("/bx/scintillator/alphadecay", this);
  fAlphaDecayCmd->SetGuidance("Put in sequence the taus of the exponentials (in nanoseconds)");

  fBetaDecayCmd   = new G4UIcmdWithAString("/bx/scintillator/betadecay", this);
  fBetaDecayCmd->SetGuidance("Put in sequence the taus of the exponentials (in nanoseconds)");
  
  fAlphaWeightCmd = new G4UIcmdWithAString("/bx/scintillator/alphaweight", this);
  fAlphaWeightCmd->SetGuidance("Put in sequence the weigths of the exponentials (be sure their sum is 1!)");
  
  fBetaWeightCmd  = new G4UIcmdWithAString("/bx/scintillator/betaweight", this);
  fBetaWeightCmd->SetGuidance("Put in sequence the weigths of the exponentials (be sure their sum is 1!)");

  fTimePPOEmissionCmd = new G4UIcmdWithADouble ("/bx/scintillator/PPOEmissionTime",this);
  fTimePPOEmissionCmd->SetGuidance("PPO emission time in absorption and reemission processes (ns)");
  
  fTimePCEmissionCmd = new G4UIcmdWithADouble ("/bx/scintillator/PCEmissionTime",this);
  fTimePCEmissionCmd->SetGuidance("PC emission time in absorption and reemission processes (ns)");
  
  fTimePCtoPPOTransferCmd = new G4UIcmdWithADouble ("/bx/scintillator/PCtoPPOTransferTime",this);
  fTimePCtoPPOTransferCmd->SetGuidance("PC to PPO energy transfer time in absorption and reemission processes (ns)");

  fDMPAttenuationLengthFactorCmd = new G4UIcmdWithADouble ("/bx/scintillator/DMPattenuationLengthFactor",this);
  fDMPAttenuationLengthFactorCmd->SetGuidance("Scaling factor for the DMP attenuation length");

  fNylonAttenuationLengthFactorCmd = new G4UIcmdWithADouble ("/bx/scintillator/NylonattenuationLengthFactor",this);
  fNylonAttenuationLengthFactorCmd->SetGuidance("Scaling factor for the nylon attenuation length");
  
  fNylonWlShiftCmd = new G4UIcmdWithADouble ("/bx/scintillator/NylonWlshift",this);
  fNylonWlShiftCmd->SetGuidance("Additive constant for the position of the step function describing nylon transmittance (nm)(");
  
  fPPOAbsReemProbXRaysCmd = new G4UIcmdWithADouble ("/bx/scintillator/PPOAbsReemProbXRays",this);
  fPPOAbsReemProbXRaysCmd->SetGuidance("Absorption reemission probability at high optical-photon energies");
  
  fPPOAbsReemProbXRaysThCmd = new G4UIcmdWithADouble ("/bx/scintillator/PPOAbsReemProbThXRays",this);
  fPPOAbsReemProbXRaysThCmd->SetGuidance("Absorption reemission probability Threshold at high optical-photon energies");

  fPPOAttenuationLengthFactorCmd = new G4UIcmdWithADouble ("/bx/scintillator/PPOattenuationLengthFactor",this);
  fPPOAttenuationLengthFactorCmd->SetGuidance("Scaling factor for the PPO attenuation length");

  fPCAttenuationLengthFactorCmd = new G4UIcmdWithADouble ("/bx/scintillator/PCattenuationLengthFactor",this);
  fPCAttenuationLengthFactorCmd->SetGuidance("Scaling factor for the PC attenuation length");

  fPMTShieldShiftCmd = new G4UIcmdWithADouble ("/bx/detector/shieldShift",this);
  fPMTShieldShiftCmd->SetGuidance("Shift along the PMT axis of the shift (towards the SSS centre). Only mm are allowed.");

  fNylonReflectivityCmd = new G4UIcmdWithADouble ("/bx/detector/nylonReflectivity",this);
  fNylonReflectivityCmd->SetGuidance("Reflectivity of Nylon in Vessels' endcaps");
  
}


BxDetectorConstructionMessenger::~BxDetectorConstructionMessenger() {
 // delete fOPERACmd;
  delete fNoPMTCmd;
  delete fDirectory;
  delete fConfCmd;
 /* delete fRockCmd;
  delete fTankCmd;*/
  delete fPMTCmd;
  delete fPMTdistribCmd;
  delete fPMTnumCmd;
  delete fSSSReflectivityCmd;
  delete fSSSspecLOBECmd;
  delete fSSSspecSPIKECmd;
  delete fSSSbackSCATTCmd;
  delete fConcentratorExternalReflectivityCmd;
  delete fConcentratorInternalReflectivityCmd;
  delete fConcentratorExternalspecSPIKECmd;
  delete fConcentratorInternalspecSPIKECmd;
  delete fPMTReflectivityCmd;
  delete fPMTringReflectivityCmd;
  delete fPMTringspecSPIKECmd;
  delete fShieldReflectivityCmd;
  delete fNylonReflectivityCmd;
  delete fBirksACmd;
  delete fBirksA2Cmd;
  delete fBirksBCmd;
  delete fBirksB2Cmd;
 /* delete fVesselOriginCmd;*/
  delete fCheckOverlapCmd;
  delete fCheckOverlapSourceCmd;
  delete fUpdateGeometryCmd;
  delete fSourceDirectory;
  delete fSourceTypeCmd;
  delete fSourceOriginCmd;
  delete fSourceRadiusCmd;
  delete fSourceLongCmd;
  delete fSourceThickCmd;
  delete fSourceXCmd;
  delete fSourceYCmd;
  delete fSourceZCmd;
  delete fSourceMaterialCmd;
  delete fSourceVialMaterialCmd;
 // delete fIsSensitiveCmd;
  delete fDefVesselShapeCmd;
  delete fRunNumberCmd;
  delete fIsRingCmd;
  delete fAlphaDecayCmd;
  delete fBetaDecayCmd;
  delete fAlphaWeightCmd;
  delete fBetaWeightCmd;
  delete fTimePPOEmissionCmd;
  delete fTimePCEmissionCmd;
  delete fTimePCtoPPOTransferCmd;
  delete fPCAttenuationLengthFactorCmd;
  delete fPPOAttenuationLengthFactorCmd;
  delete fDMPAttenuationLengthFactorCmd;
  delete fNylonAttenuationLengthFactorCmd;
  delete fNylonWlShiftCmd;
  delete fPMTShieldShiftCmd;
  delete fPPOAbsReemProbXRaysCmd;
  delete fPPOAbsReemProbXRaysThCmd;
}


void BxDetectorConstructionMessenger::SetNewValue(G4UIcommand* command, G4String newValue) { 
	/*if(command == fOPERACmd) {
    fDetectorConstruction->SetOperaDetector(fOPERACmd->ConvertToBool(newValue));
  } else if(command == fIsSensitiveCmd) {
    fDetectorConstruction->SetSensitive(fIsSensitiveCmd->ConvertToBool(newValue));
    BxLog(routine) << "Zone I (bulk) is the sensitive detector" <<  endlog;   
  } else*/ if(command == fNoPMTCmd) {
    fDetectorConstruction->SetNoPMTs(fNoPMTCmd->ConvertToBool(newValue));
    if(fNoPMTCmd->ConvertToBool(newValue))
    BxLog(routine) << "All PMTs removed " <<  endlog;   
  } else if(command == fCheckOverlapSourceCmd) {
    fDetectorConstruction->SetOverlapSource(fCheckOverlapCmd->ConvertToBool(newValue));
    BxLog(routine) << "Check Source volume overlapping " << newValue << endlog;   
  } else if(command == fCheckOverlapCmd) {
    fDetectorConstruction->SetOverlap(fCheckOverlapCmd->ConvertToBool(newValue));
    BxLog(routine) << "Check volume overlapping " << newValue << endlog;   
  } else if(command == fUpdateGeometryCmd) {
    fDetectorConstruction->UpdateGeometry();
       BxLog(routine) << "Geometry parameters updated" << endlog;    
  }  else if(command == fTankCmd) {
    fDetectorConstruction->SetTank(fTankCmd->ConvertToBool(newValue));
    if(!fTankCmd->ConvertToBool(newValue)) 
       BxLog(routine) << "Outer tank disabled" << endlog;    
  } /*else if(command == fRockCmd) {
    fDetectorConstruction->SetRock(fRockCmd->ConvertToBool(newValue)); 
    if(fRockCmd->ConvertToBool(newValue))   
       BxLog(routine) << "Build rock and concrete" << endlog;    
  } */else if(command == fConfCmd) {
    if(newValue == "full") {
      fDetectorConstruction->SetDetectorFlag(0);
      BxLog(routine) << "Detector Configuration 0: " << newValue << endlog;
    } else if(newValue == "fullWater") {
      fDetectorConstruction->SetDetectorFlag(1);
      BxLog(routine) << "Detector Configuration 1: " << newValue << endlog;
    } else if(newValue == "inWater") {
      fDetectorConstruction->SetDetectorFlag(2);
      fDetectorConstruction->SetTank(0);
      BxLog(routine) << "Detector Configuration 2: " << newValue << endlog;
    } else if(newValue == "inScint") {
      fDetectorConstruction->SetDetectorFlag(4);
      fDetectorConstruction->SetTank(0);
      BxLog(routine) << "Detector Configuration 4: " << newValue << endlog;
    } else if(newValue == "outer") {
      fDetectorConstruction->SetDetectorFlag(3);
      BxLog(routine) << "Detector Configuration 3: " << newValue << endlog;
    } else if(newValue == "ctf") {
      fDetectorConstruction->SetDetectorFlag(5);
      BxLog(routine) << "Detector Configuration 5: " << newValue << endlog;
      fDetectorConstruction->SetTank(0);
      BxReadParameters::Get()->SetLightYieldScale(21./30.);
    }
  } else if(command == fPMTCmd) {
    if(newValue == "fullPMT") {
      fDetectorConstruction->SetPMTFlag(0);
      BxLog(routine) << "Full PMT Configuration 0: " << newValue << endlog;
    } else if(newValue == "diskPMT") {
      fDetectorConstruction->SetPMTFlag(1);
      BxLog(routine) << "Disk PMT Configuration 1: " << newValue << endlog;
    } else if(newValue == "diskPMTwithConcentrator") {
      fDetectorConstruction->SetPMTFlag(2);
      BxLog(routine) << "Disk PMT Configuration 2: " << newValue << endlog;
    } else if(newValue == "simplePMT") {
      fDetectorConstruction->SetPMTFlag(3);
      BxLog(routine) << "Spherical Cathode with Concentrator PMT Configuration 3: " << newValue << endlog;
    } 
  } else if(command == fPMTdistribCmd) {
    if(newValue == "uniform") {
      fDetectorConstruction->SetPMTDistributionFlag(0);
      BxLog(routine) << "Distribution of PMTs on the sphere is uniform and not the real one" << endlog;
    } else if(newValue == "real") {
      fDetectorConstruction->SetPMTDistributionFlag(1);
      BxLog(routine) << "Distribution of PMTs on the sphere is the real one" << endlog;}
  } else if(command == fPMTnumCmd) {
      fDetectorConstruction->SetNumberOfPMT(atoi(newValue.c_str()));
      BxLog(routine) << "Number of uniformly distributed PMTs is " << newValue <<  endlog;
  }  else if(command == fNoConcCmd) {
    fDetectorConstruction->SetNoConcentrators(fNoConcCmd->ConvertToBool(newValue));
    if(fTankCmd->ConvertToBool(newValue)) 
  	BxLog(routine) << "NO concentrators!" << endlog;
	} else if(command == fSSSReflectivityCmd) {
      BxPropertyCollection::Get()->SetSSSReflectivity(fSSSReflectivityCmd->ConvertToDouble(newValue));
      BxLog(routine) << "SSS Reflectivity = " << newValue << endlog;
  } else if(command == fSSSspecLOBECmd) {
      BxPropertyCollection::Get()->SetSSSspecLOBE(fSSSspecLOBECmd->ConvertToDouble(newValue));
      BxLog(routine) << "SSSspecLOBE = " << newValue << endlog;
    } else if(command == fSSSspecSPIKECmd) {
      BxPropertyCollection::Get()->SetSSSspecSPIKE(fSSSspecSPIKECmd->ConvertToDouble(newValue));
      BxLog(routine) << "SSSspecSPIKE = " << newValue << endlog;
    	} else if(command == fSSSbackSCATTCmd) {
      BxPropertyCollection::Get()->SetSSSbackSCATT(fSSSbackSCATTCmd->ConvertToDouble(newValue));
      BxLog(routine) << "SSSbackSCATT = " << newValue << endlog;
    	} else if(command == fConcentratorExternalReflectivityCmd) {
      BxPropertyCollection::Get()->SetConcentratorExternalReflectivity(fConcentratorExternalReflectivityCmd->ConvertToDouble(newValue));
      BxLog(routine) << "Concentrator External Reflectivity = " << newValue << endlog;
    	} else if(command == fConcentratorInternalReflectivityCmd) {
      BxPropertyCollection::Get()->SetConcentratorInternalReflectivity(fConcentratorInternalReflectivityCmd->ConvertToDouble(newValue));
      BxLog(routine) << "Concentrator Internal Reflectivity = " << newValue << endlog;
    	} else if(command == fConcentratorExternalspecSPIKECmd) {
      BxPropertyCollection::Get()->SetConcentratorExternalspecSPIKE(fConcentratorExternalspecSPIKECmd->ConvertToDouble(newValue));
      BxLog(routine) << "Concentrator External specular reflectivity = " << newValue << " (diffusive=" <<1-fConcentratorExternalspecSPIKECmd->ConvertToDouble(newValue) << ")" << endlog;
    	} else if(command == fConcentratorInternalspecSPIKECmd) {
      BxPropertyCollection::Get()->SetConcentratorInternalspecSPIKE(fConcentratorInternalspecSPIKECmd->ConvertToDouble(newValue));
      BxLog(routine) << "Concentrator Internal specular reflectivity = " << newValue << " (diffusive=" <<1-fConcentratorInternalspecSPIKECmd->ConvertToDouble(newValue) << ")" << endlog;
	} else if(command == fPMTReflectivityCmd) {
      BxPropertyCollection::Get()->SetCathodeReflectivity(fPMTReflectivityCmd->ConvertToDouble(newValue));
      BxLog(routine) << "Cathode Reflectivity = " << newValue << endlog;
	} else if(command == fPMTringReflectivityCmd) {
      BxPropertyCollection::Get()->SetPMTringReflectivity(fPMTringReflectivityCmd->ConvertToDouble(newValue));
      BxLog(routine) << "PMT ring Reflectivity = " << newValue << endlog;
	} else if(command == fPMTringspecSPIKECmd) {
      BxPropertyCollection::Get()->SetPMTringspecSPIKE(fPMTringspecSPIKECmd->ConvertToDouble(newValue));
      BxLog(routine) << "PMT ring specular reflectivity = " << newValue << " (diffusive=" <<1-fPMTringspecSPIKECmd->ConvertToDouble(newValue) << ")" << endlog;
	} else if(command == fShieldReflectivityCmd) {
      BxPropertyCollection::Get()->SetShieldReflectivity(fShieldReflectivityCmd->ConvertToDouble(newValue));
      BxLog(routine) << "PMT's Shield Reflectivity = " << newValue << endlog;
	} else if(command == fBirksACmd) {
      BxReadParameters::Get()->SetBirksAlpha(fBirksACmd->ConvertToDouble(newValue));
      BxLog(routine) << "Alpha Birks parameter = " << newValue << " cm/MeV" << endlog;
  } else if(command == fBirksA2Cmd) {
      BxReadParameters::Get()->SetBirksSecondOrderAlpha(fBirksA2Cmd->ConvertToDouble(newValue));
      BxLog(routine) << "Alpha Second Order Birks parameter = " << newValue << " (cm/MeV)^2" << endlog;
  } else if(command == fBirksBCmd) {
      BxReadParameters::Get()->SetBirksBeta(fBirksBCmd->ConvertToDouble(newValue));
      BxLog(routine) << "Beta Birks parameter = " << newValue << " cm/MeV" << endlog;
  } else if(command == fBirksB2Cmd) {
      BxReadParameters::Get()->SetBirksSecondOrderBeta(fBirksB2Cmd->ConvertToDouble(newValue));
      BxLog(routine) << "Beta Second Order Birks parameter = " << newValue << " (cm/MeV)^2" << endlog;
  } else if(command == fVesselOriginCmd) {
      BxReadParameters::Get()->SetVesselOrigin(fVesselOriginCmd->ConvertToDimensioned3Vector(newValue));
      BxLog(routine) << "Vessel Origin = " << fVesselOriginCmd->ConvertToDimensioned3Vector(newValue)/cm << " cm" << endlog;
  } else if (command == fRZDimCmd){
      BxReadParameters::Get()->SetRZDim(fRZDimCmd->ConvertToDimensionedDouble(newValue));
      BxLog(routine) << "Step in the vessel's polycone = " << fRZDimCmd->ConvertToDimensionedDouble(newValue)/cm << " cm" << endlog;
  }else if(command == fSourceTypeCmd) {
           if(newValue == "sphere")   fDetectorConstruction->SetSourceType(1);
      else if(newValue == "cylinder") fDetectorConstruction->SetSourceType(2);
      else if(newValue == "box")      fDetectorConstruction->SetSourceType(3);
      else if(newValue == "Am-Be")    fDetectorConstruction->SetSourceType(4);
      else {
        BxLog(error) << newValue << " is not an accepted source geometry " << endlog; 
	BxLog(fatal) << endlog; 
      }
      BxLog(routine) << "Source Type: " << newValue << endlog;
  } else if(command == fSourceOriginCmd) {
      fDetectorConstruction->SetSourceOrigin(fVesselOriginCmd->ConvertToDimensioned3Vector(newValue));
      BxLog(routine) << "Source Origin = " << fVesselOriginCmd->ConvertToDimensioned3Vector(newValue)/cm << " cm" << endlog;
  } else if(command == fSourceRadiusCmd) {
      fDetectorConstruction->SetSourceRadius(fSourceRadiusCmd->ConvertToDimensionedDouble(newValue));
  } else if(command == fSourceLongCmd) {
      fDetectorConstruction->SetSourceLong(fSourceLongCmd->ConvertToDimensionedDouble(newValue));
  } else if(command == fSourceThickCmd) {
      fDetectorConstruction->SetSourceVialThick(fSourceThickCmd->ConvertToDimensionedDouble(newValue));
  } else if(command == fSourceXCmd) {
      fDetectorConstruction->SetSourceX(fSourceXCmd->ConvertToDimensionedDouble(newValue));
  } else if(command == fSourceYCmd) {
      fDetectorConstruction->SetSourceY(fSourceYCmd->ConvertToDimensionedDouble(newValue));
  } else if(command == fSourceZCmd) {
      fDetectorConstruction->SetSourceZ(fSourceZCmd->ConvertToDimensionedDouble(newValue));
  } else if(command == fSourceMaterialCmd) {
      if(newValue == "pc")            fDetectorConstruction->SetSourceMaterial(BxMaterial::Get()->GetPC());
      if(newValue == "scintillator")  fDetectorConstruction->SetSourceMaterial(BxMaterial::Get()->GetScintillator());
      if(newValue == "water")         fDetectorConstruction->SetSourceMaterial(BxMaterial::Get()->GetWater());
      if(newValue == "nylon")         fDetectorConstruction->SetSourceMaterial(BxMaterial::Get()->GetNylon());
      if(newValue == "quartz")        fDetectorConstruction->SetSourceMaterial(BxMaterial::Get()->GetQuartz());
      if(newValue == "buffer")        fDetectorConstruction->SetSourceMaterial(BxMaterial::Get()->GetDMPbuffer());
      if(newValue == "steel")  BxLog(fatal) << "source made of steel not yet implemented" << endlog;
	      //fDetectorConstruction->SetSourceMaterial(BxMaterial::Get()->GetSteel());
      if(newValue == "air")           fDetectorConstruction->SetSourceMaterial(BxMaterial::Get()->GetAir());
      if(newValue == "vacuum")    BxLog(fatal) << "source made of vacuum not yet implemented" << endlog;    
	      //fDetectorConstruction->SetSourceMaterial(BxMaterial::Get()->GetVacuum());
      BxLog(routine) << "Source Material: " << newValue << endlog;
  } else if(command == fSourceVialMaterialCmd) {
      if(newValue == "pc")            fDetectorConstruction->SetSourceVialMaterial(BxMaterial::Get()->GetPC());
      if(newValue == "scintillator")  fDetectorConstruction->SetSourceVialMaterial(BxMaterial::Get()->GetScintillator());
      if(newValue == "water")         fDetectorConstruction->SetSourceVialMaterial(BxMaterial::Get()->GetWater());
      if(newValue == "nylon")         fDetectorConstruction->SetSourceVialMaterial(BxMaterial::Get()->GetNylon());
      if(newValue == "quartz")        fDetectorConstruction->SetSourceVialMaterial(BxMaterial::Get()->GetQuartz());
      if(newValue == "buffer")        fDetectorConstruction->SetSourceVialMaterial(BxMaterial::Get()->GetDMPbuffer());
      if(newValue == "steel")   BxLog(fatal) << "source vial made of steel not yet implemented" << endlog;
	      //fDetectorConstruction->SetSourceVialMaterial(BxMaterial::Get()->GetSteel());
      if(newValue == "air")           fDetectorConstruction->SetSourceVialMaterial(BxMaterial::Get()->GetAir());
      if(newValue == "vacuum")      BxLog(fatal) << "source vial made of vacuum not yet implemented" << endlog;
	      //fDetectorConstruction->SetSourceVialMaterial(BxMaterial::Get()->GetVacuum());
      BxLog(routine) << "Source Vial Material: " << newValue << endlog;
  } else if(command == fDefVesselShapeCmd) {
    BxReadParameters::Get()->SetIsDVessel(fDefVesselShapeCmd->ConvertToBool(newValue));
    BxLog(routine) << "Deformed vessel activation: " << newValue  <<  endlog;   
  } else if(command == fRunNumberCmd) {
    BxReadParameters::Get()->SetRunNumber(atoi(newValue.c_str()));
    BxLog(routine) << "Run Number being simulated " << newValue  <<  endlog;   
  } else if(command == fIsRingCmd) {
    fDetectorConstruction->SetEndCaps(fIsRingCmd->ConvertToBool(newValue));
    BxLog(routine) << "EndCaps activation: " << newValue  <<  endlog;   
  } else if(command == fAlphaDecayCmd) {
    BxReadParameters::Get()->SetAlphaDecayTimeConstant(BxReadParameters::Get()->ConvertTo4DoubleVector(newValue));
    BxLog(routine) << "Set Alpha Decay Constants to: " << newValue  <<  endlog;   
  } else if(command == fBetaDecayCmd) {
    BxReadParameters::Get()->SetBetaDecayTimeConstant(BxReadParameters::Get()->ConvertTo4DoubleVector(newValue));
    BxLog(routine) << "Set Beta Decay Constants to: " << newValue  <<  endlog;   
  } else if(command == fAlphaWeightCmd) {
    BxReadParameters::Get()->SetAlphaDecayWeight(BxReadParameters::Get()->ConvertTo4DoubleVector(newValue));
    BxLog(routine) << "Set Alpha Weight Constants to: " << newValue  <<  endlog;   
  } else if(command == fBetaWeightCmd) {
    BxReadParameters::Get()->SetBetaDecayWeight(BxReadParameters::Get()->ConvertTo4DoubleVector(newValue));
    BxLog(routine) << "Set Beta Weight Constants to: " << newValue  <<  endlog;   
  } else if(command == fTimePPOEmissionCmd) {
    BxReadParameters::Get()->SetTimePPOEmission(fTimePPOEmissionCmd->ConvertToDouble(newValue));
    BxLog(routine) << "Set PPO Emission Time to: " << newValue  <<  endlog;   
  } else if(command == fTimePCEmissionCmd) {
    BxReadParameters::Get()->SetTimePCEmission(fTimePCEmissionCmd->ConvertToDouble(newValue));
    BxLog(routine) << "Set PC Emission Time to: " << newValue  <<  endlog;   
  } else if(command == fTimePCtoPPOTransferCmd) {
    BxReadParameters::Get()->SetTimePCtoPPOTransfer(fTimePCtoPPOTransferCmd->ConvertToDouble(newValue));
    BxLog(routine) << "Set PC to PPO energy transfer time to: " << newValue  <<  endlog;   
  } else if(command == fPPOAttenuationLengthFactorCmd) {
    //BxReadParameters::Get()->SetPPOAttenuationLengthFactor(fPPOAttenuationLengthFactorCmd->ConvertToDouble(newValue));
    //No need to call this method now since the parameter has already been updated in g4bx2.cc
    BxLog(routine) << "Set PPO attenuation length scaling factor to: " << newValue  <<  endlog;
  } else if(command == fPCAttenuationLengthFactorCmd) {
    //BxReadParameters::Get()->SetPCAttenuationLengthFactor(fPCAttenuationLengthFactorCmd->ConvertToDouble(newValue));
    //No need to call this method now since the parameter has already been updated in g4bx2.cc
    BxLog(routine) << "Set PC attenuation length scaling factor to: " << newValue  <<  endlog;
  } else if(command == fDMPAttenuationLengthFactorCmd) {
    //BxReadParameters::Get()->SetDMPAttenuationLengthFactor(fDMPAttenuationLengthFactorCmd->ConvertToDouble(newValue));
    //No need to call this method now since the parameter has already been updated in g4bx2.cc
    BxLog(routine) << "Set DMP attenuation length scaling factor to: " << newValue  <<  endlog;
  } else if(command == fNylonAttenuationLengthFactorCmd) {
    //BxReadParameters::Get()->SetNylonAttenuationLengthFactor(fNylonAttenuationLengthFactorCmd->ConvertToDouble(newValue));
    //No need to call this method now since the parameter has already been updated in g4bx2.cc
    BxLog(routine) << "Set Nylon attenuation length scaling factor to: " << newValue  <<  endlog;
  } else if(command == fNylonWlShiftCmd) {
    //BxReadParameters::Get()->SetNylonWlShift(fNylonWlShiftCmd->ConvertToDouble(newValue));
    //No need to call this method now since the parameter has already been updated in g4bx2.cc
    BxLog(routine) << "Set Nylon wavelength shift to (nm): " << newValue  <<  endlog;
  } else if(command == fPPOAbsReemProbXRaysCmd) {
    //This parameter has already been updated in g4bx2.cc
    BxLog(routine) << "Set PPO reemission probability for UV photons to: " << newValue  <<  endlog;
  } else if(command == fPPOAbsReemProbXRaysThCmd) {
    //This parameter has already been updated in g4bx2.cc
    BxLog(routine) << "Set wavelength threshold above which the reemission probability decreases to: " << newValue  <<  endlog;
  } else if(command == fPMTShieldShiftCmd) {
    BxReadParameters::Get()->SetPMTShieldShift(fPMTShieldShiftCmd->ConvertToDouble(newValue));
    BxLog(routine) << "Set PMTs' shield distance from SSS to: " << BxReadParameters::Get()->GetPMTShieldZPosition()/cm  << " cm" <<  endlog;   
  } else if(command == fNylonReflectivityCmd) {
    BxPropertyCollection::Get()->SetNylonReflectivity(fNylonReflectivityCmd->ConvertToDouble(newValue));
    BxLog(routine) << "Endcaps' Nylon Reflectivity set to " << fNylonReflectivityCmd->ConvertToDouble(newValue) <<  endlog;   
}

}

