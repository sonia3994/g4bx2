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
// *********--***********************************************************
//
// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "BxGeneratorSolarNeutrino2Messenger.hh"
#include "G4UIcmdWithADouble.hh"
#include "BxGeneratorSolarNeutrino2.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "BxReadParameters.hh"
#include "BxOutputVertex.hh"
#include "BxLogger.hh"
#include "G4SystemOfUnits.hh"
using namespace std;

BxGeneratorSolarNeutrino2Messenger::BxGeneratorSolarNeutrino2Messenger(BxGeneratorSolarNeutrino2* gen){
  generator = gen;
  fDirectory = new G4UIdirectory("/bx/generator/snu2/");
  fDirectory->SetGuidance("Control of BxSolarNeutrino2 event generator");


  fPositionCmd = new G4UIcmdWith3VectorAndUnit("/bx/generator/snu2/position",this);
  fPositionCmd->SetGuidance("Set the gun position");
  fPositionCmd->SetUnitCategory("Length");
  fPositionCmd->SetDefaultUnit("cm");
  fPositionCmd->SetUnitCandidates("mm cm m");
  
  fSphereBulkCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/snu2/sphere_radius",this);
  fSphereBulkCmd->SetGuidance("Bulk radius");
  fSphereBulkCmd->SetUnitCategory("Length");
  fSphereBulkCmd->SetDefaultUnit("cm");
  fSphereBulkCmd->SetUnitCandidates("mm cm m");
  
  fSphereBulkRadMinCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/snu2/sphere_radius_min",this);
  fSphereBulkRadMinCmd->SetGuidance("Bulk minimum radius");
  fSphereBulkRadMinCmd->SetGuidance("Default: 0");
  fSphereBulkRadMinCmd->SetUnitCategory("Length");
  fSphereBulkRadMinCmd->SetDefaultUnit("cm");
  fSphereBulkRadMinCmd->SetUnitCandidates("mm cm m");

  fSphereSurfCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/snu2/surface_radius",this);
  fSphereSurfCmd->SetGuidance("Surface radius");
  fSphereSurfCmd->SetUnitCategory("Length");
  fSphereSurfCmd->SetDefaultUnit("cm");
  fSphereSurfCmd->SetUnitCandidates("mm cm m");

  fSphereCentreCmd = new G4UIcmdWith3VectorAndUnit("/bx/generator/snu2/sphere_origin",this);
  fSphereCentreCmd->SetGuidance("Set the sphere position");
  fSphereCentreCmd->SetGuidance("Default: 0. 0. 0. cm");
  fSphereCentreCmd->SetUnitCategory("Length");
  fSphereCentreCmd->SetDefaultUnit("cm");
  fSphereCentreCmd->SetUnitCandidates("mm cm m");

  fVesselCmd = new G4UIcmdWithABool("/bx/generator/snu2/vessel",this);
  fVesselCmd->SetGuidance("Set events on the vessel inner surface");

  fOldSigmaCmd = new G4UIcmdWithABool("/bx/generator/snu2/oldsigma",this);
  fOldSigmaCmd->SetGuidance("Set tree level sigma in snu");

  fBulkCmd = new G4UIcmdWithABool("/bx/generator/snu2/bulk",this);
  fBulkCmd->SetGuidance("Set events on the vessel inner surface");

  fBufferCmd = new G4UIcmdWithABool("/bx/generator/snu2/buffer",this);
  fBufferCmd->SetGuidance("Set events on the vessel inner surface");
 
  fDirectionCmd = new G4UIcmdWith3Vector("/bx/generator/snu2/direction",this);
  fDirectionCmd->SetGuidance("Set the gun direction");
  

  // confine to volume
  fConfineCmd = new G4UIcmdWithAString("/bx/generator/snu2/confine",this);
  fConfineCmd->SetGuidance("Confine source to volume (NULL to unset).");
  fConfineCmd->SetGuidance("usage: confine VolName");
  fConfineCmd->SetParameterName("VolName",true,true);
  fConfineCmd->SetDefaultValue("NULL");
  
  fNeutrinoCmd = new G4UIcmdWithAString("/bx/generator/snu2/source",this);
  G4String candidates = "Be7 pep pp CNO B8 N13 O15 F17 Be7_862 Be7_384 hep";  
  fNeutrinoCmd->SetCandidates(candidates);
 
  fDM2Cmd  =  new G4UIcmdWithADouble("/bx/generator/snu2/deltam2",this);
  
  fTG2TCmd =  new G4UIcmdWithADouble("/bx/generator/snu2/tan2t",this);
  
  fNumberOfStepCmd =  new G4UIcmdWithAnInteger("/bx/generator/snu2/stepnumber",this);
}


BxGeneratorSolarNeutrino2Messenger::~BxGeneratorSolarNeutrino2Messenger() {

  delete fDirectory;
  delete fPositionCmd; 
  delete fDirectionCmd;
  delete fSphereBulkCmd; 
  delete fSphereSurfCmd;
  delete fSphereCentreCmd;
  delete fVesselCmd;
  delete fOldSigmaCmd;
  delete fBulkCmd;
  delete fBufferCmd;
  delete fConfineCmd;
  delete fNeutrinoCmd;
  delete fDM2Cmd;
  delete fTG2TCmd;
  delete fNumberOfStepCmd;
}


void BxGeneratorSolarNeutrino2Messenger::SetNewValue(G4UIcommand* cmd, G4String newValue) { 
   if (cmd == fPositionCmd) {         
	 generator->SetParticlePosition(fPositionCmd->ConvertToDimensioned3Vector(newValue));
         BxOutputVertex::Get()->SetSpatialDist(0);	 
   } else if (cmd == fDirectionCmd){
   	 generator->SetParticleMomentumDirection(fDirectionCmd->ConvertTo3Vector(newValue));
   } else if (cmd == fSphereBulkCmd){
         generator->SetIfVolDist(true);
	 generator->SetPosDisType("Volume");
         generator->SetPosDisShape("Sphere");
         generator->SetRadius(fSphereBulkCmd->ConvertToDimensionedDouble(newValue));	 
 	 BxOutputVertex::Get()->SetRagMin(0);
	 BxOutputVertex::Get()->SetRagMax(fSphereBulkCmd->ConvertToDimensionedDouble(newValue));
         BxOutputVertex::Get()->SetSpatialDist(1);	 
   } else if (cmd == fSphereSurfCmd){
         generator->SetIfVolDist(true);
	 generator->SetPosDisType("Surface");
         generator->SetPosDisShape("Sphere");
         generator->SetRadius(fSphereSurfCmd->ConvertToDimensionedDouble(newValue));	 
         BxOutputVertex::Get()->SetSpatialDist(1);	 
   } else if (cmd == fSphereCentreCmd){
   	 generator->SetCentreCoords(fSphereCentreCmd->ConvertToDimensioned3Vector(newValue));
   } else if (cmd == fSphereBulkRadMinCmd){
   	 generator->SetRadius0(fSphereBulkRadMinCmd->ConvertToDimensionedDouble(newValue));
 	 BxOutputVertex::Get()->SetRagMin(fSphereBulkRadMinCmd->ConvertToDimensionedDouble(newValue));
   } else if (cmd == fVesselCmd){
	 generator->SetIfVolDist(true);
	 generator->SetPosDisType("Surface");
         generator->SetPosDisShape("Sphere");
         generator->SetRadius(BxReadParameters::Get()->GetZoneIExternalRadius() - 0.001*mm); 
         BxOutputVertex::Get()->SetSpatialDist(1);	 
	 BxOutputVertex::Get()->SetRagMin(BxReadParameters::Get()->GetZoneIExternalRadius() - 0.001*mm);
	 BxOutputVertex::Get()->SetRagMax(BxReadParameters::Get()->GetZoneIExternalRadius() - 0.001*mm);
   }else if(cmd == fOldSigmaCmd){
           generator->SetOldSigma(newValue);
} else if (cmd == fBulkCmd){
	 generator->SetIfVolDist(true);
	 generator->SetPosDisType("Volume");
         generator->SetPosDisShape("Sphere");
         generator->SetRadius(BxReadParameters::Get()->GetZoneIExternalRadius() - 0.001*mm); 
         BxOutputVertex::Get()->SetSpatialDist(1);	 
	 BxOutputVertex::Get()->SetRagMin(0);
	 BxOutputVertex::Get()->SetRagMax(BxReadParameters::Get()->GetZoneIExternalRadius() - 0.001*mm);
   } else if (cmd == fBufferCmd){
	 generator->SetIfVolDist(true);
	 generator->SetPosDisType("Volume");
         generator->SetPosDisShape("Sphere");
         generator->SetRadius0(BxReadParameters::Get()->GetZoneIExternalRadius()); 
         generator->SetRadius (BxReadParameters::Get()->GetZoneIIIExternalRadius()); 
         BxOutputVertex::Get()->SetSpatialDist(1);	 
 	 BxOutputVertex::Get()->SetRagMin(BxReadParameters::Get()->GetZoneIExternalRadius());
	 BxOutputVertex::Get()->SetRagMax(BxReadParameters::Get()->GetZoneIIIExternalRadius());
   } else if (cmd == fConfineCmd){
	   generator->SetIfVolDist(true);
	   generator->ConfineSourceToVolume(newValue);
   } else if (cmd == fNeutrinoCmd){
	   BxLog(routine) << "Neutrino source: " << newValue << endlog ;
	   if(newValue == "pp") {
		   generator->SetNeutrinoType(BxGeneratorSolarNeutrino2::pp);
	   } else if(newValue == "Be7") {
		   generator->SetNeutrinoType(BxGeneratorSolarNeutrino2::be7);
	   } else if(newValue == "Be7_862") {
		   generator->SetNeutrinoType(BxGeneratorSolarNeutrino2::be7_862);
	   } else if(newValue == "Be7_384") {
		   generator->SetNeutrinoType(BxGeneratorSolarNeutrino2::be7_384);
	   } else if(newValue == "B8") {
		   generator->SetNeutrinoType(BxGeneratorSolarNeutrino2::b8);
	   } else if(newValue == "N13") {
		   generator->SetNeutrinoType(BxGeneratorSolarNeutrino2::n13);
	   } else if(newValue == "O15") {
		   generator->SetNeutrinoType(BxGeneratorSolarNeutrino2::o15);
	   } else if(newValue == "F17") {
		   generator->SetNeutrinoType(BxGeneratorSolarNeutrino2::f17);
	   } else if(newValue == "pep") {
		   generator->SetNeutrinoType(BxGeneratorSolarNeutrino2::pep);
	   } else if(newValue == "hep") {
		   generator->SetNeutrinoType(BxGeneratorSolarNeutrino2::hep);
	   } else if(newValue == "CNO") {
		   generator->SetNeutrinoType(BxGeneratorSolarNeutrino2::cno);
	   }
   } else if (cmd == fDM2Cmd){
	   generator->SetDM2(fDM2Cmd->ConvertToDouble(newValue));
   } else if (cmd == fTG2TCmd){
	   generator->SetTG2T(fTG2TCmd->ConvertToDouble(newValue));
   } else if (cmd == fNumberOfStepCmd){
	   generator->SetBinning(fNumberOfStepCmd->ConvertToInt(newValue));
   }       
}
