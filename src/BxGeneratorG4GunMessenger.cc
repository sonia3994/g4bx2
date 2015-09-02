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
//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "BxGeneratorG4GunMessenger.hh"
#include "G4UIcmdWithADouble.hh"
#include "BxGeneratorG4Gun.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "BxReadParameters.hh"
#include "BxOutputVertex.hh"
#include "G4SystemOfUnits.hh"
using namespace std;

BxGeneratorG4GunMessenger::BxGeneratorG4GunMessenger(BxGeneratorG4Gun* gen){
	generator = gen;
	fDirectory = new G4UIdirectory("/bx/generator/g4gun/");
	fDirectory->SetGuidance("Control of BxG4Gun event generator");


	fPositionCmd = new G4UIcmdWith3VectorAndUnit("/bx/generator/g4gun/position",this);
	fPositionCmd->SetGuidance("Set the gun position");
	fPositionCmd->SetUnitCategory("Length");
	fPositionCmd->SetDefaultUnit("cm");
	fPositionCmd->SetUnitCandidates("mm cm m");

	fSphereBulkCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/g4gun/sphere_radius",this);
	fSphereBulkCmd->SetGuidance("Bulk radius");
	fSphereBulkCmd->SetUnitCategory("Length");
	fSphereBulkCmd->SetDefaultUnit("cm");
	fSphereBulkCmd->SetUnitCandidates("mm cm m");

	fSphereBulkRadMinCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/g4gun/sphere_radius_min",this);
	fSphereBulkRadMinCmd->SetGuidance("Bulk minimum radius");
	fSphereBulkRadMinCmd->SetGuidance("Default: 0");
	fSphereBulkRadMinCmd->SetUnitCategory("Length");
	fSphereBulkRadMinCmd->SetDefaultUnit("cm");
	fSphereBulkRadMinCmd->SetUnitCandidates("mm cm m");

	fSphereSurfCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/g4gun/surface_radius",this);
	fSphereSurfCmd->SetGuidance("Surface radius");
	fSphereSurfCmd->SetUnitCategory("Length");
	fSphereSurfCmd->SetDefaultUnit("cm");
	fSphereSurfCmd->SetUnitCandidates("mm cm m");

	fSphereCentreCmd = new G4UIcmdWith3VectorAndUnit("/bx/generator/g4gun/sphere_origin",this);
	fSphereCentreCmd->SetGuidance("Set the sphere position");
	fSphereCentreCmd->SetGuidance("Default: 0. 0. 0. cm");
	fSphereCentreCmd->SetUnitCategory("Length");
	fSphereCentreCmd->SetDefaultUnit("cm");
	fSphereCentreCmd->SetUnitCandidates("mm cm m");
	fVesselCmd = new G4UIcmdWithABool("/bx/generator/g4gun/vessel",this);
	fVesselCmd->SetGuidance("Set events on the vessel inner surface");

	fBulkCmd = new G4UIcmdWithABool("/bx/generator/g4gun/bulk",this);
	fBulkCmd->SetGuidance("Set events on the vessel inner surface");

	fBufferCmd = new G4UIcmdWithABool("/bx/generator/g4gun/buffer",this);
	fBufferCmd->SetGuidance("Set events on the vessel inner surface");
	fDirectionCmd = new G4UIcmdWith3Vector("/bx/generator/g4gun/direction",this);
	fDirectionCmd->SetGuidance("Set the gun direction");

	fEnergyCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/g4gun/energy",this);
	fEnergyCmd->SetGuidance("Set the gun energy");
	fEnergyCmd->SetUnitCategory("Energy");
	fEnergyCmd->SetDefaultUnit("MeV");
	fEnergyCmd->SetUnitCandidates("eV keV MeV GeV");


	
	fEnergyDisTypeCmd = new G4UIcmdWithAString("/bx/generator/g4gun/dist_energy",this);
	G4String candidates = "Lin Pow Exp Gauss Brem BBody Cdg";  
	fEnergyDisTypeCmd->SetCandidates(candidates);

	fSetEminCmd      = new G4UIcmdWithADoubleAndUnit("/bx/generator/g4gun/emin",this);
	fSetEminCmd->SetDefaultUnit("MeV");
	fSetEminCmd->SetUnitCandidates("eV keV MeV GeV");

	fSetEmaxCmd      = new G4UIcmdWithADoubleAndUnit("/bx/generator/g4gun/emax",this);
	fSetEmaxCmd->SetDefaultUnit("MeV");
	fSetEmaxCmd->SetUnitCandidates("eV keV MeV GeV");

	fSetAlphaCmd     = new G4UIcmdWithADouble("/bx/generator/g4gun/alpha",this);

	fSetTempCmd      = new G4UIcmdWithADoubleAndUnit("/bx/generator/g4gun/temp",this);
	fSetTempCmd->SetDefaultUnit("kelvin");
	fSetTempCmd->SetUnitCandidates("kelvin");

	fSetEzeroCmd     = new G4UIcmdWithADoubleAndUnit("/bx/generator/g4gun/ezero",this);
	fSetEzeroCmd->SetDefaultUnit("MeV");
	fSetEzeroCmd->SetUnitCandidates("eV keV MeV GeV");

	fSetGradientCmd  = new G4UIcmdWithADouble("/bx/generator/g4gun/gradient",this);

	fSetInterCeptCmd = new G4UIcmdWithADouble("/bx/generator/g4gun/intercept",this);


	fParticleCmd = new G4UIcmdWithAString("/bx/generator/g4gun/particle",this);
	fParticleCmd->SetGuidance("Set the gun type");

	// confine to volume
	fConfineCmd = new G4UIcmdWithAString("/bx/generator/g4gun/confine",this);
	fConfineCmd->SetGuidance("Confine source to volume (NULL to unset).");
	fConfineCmd->SetGuidance("usage: confine VolName");
	fConfineCmd->SetParameterName("VolName",true,true);
	fConfineCmd->SetDefaultValue("NULL");

}


BxGeneratorG4GunMessenger::~BxGeneratorG4GunMessenger() {
	delete fDirectory;
	delete fPositionCmd; 
	delete fDirectionCmd;
	delete fEnergyCmd;
	delete fParticleCmd; 
	delete fSphereBulkCmd; 
	delete fSphereSurfCmd;
	delete fSphereCentreCmd;
	delete fVesselCmd;
	delete fBulkCmd;
	delete fBufferCmd;
	delete fConfineCmd;
	delete fEnergyDisTypeCmd;
	delete fSetEminCmd;
	delete fSetEmaxCmd;	
	delete fSetAlphaCmd;	
	delete fSetTempCmd;	
	delete fSetEzeroCmd;	
	delete fSetGradientCmd;
	delete fSetInterCeptCmd;

}

void BxGeneratorG4GunMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) { 
	if (cmd == fPositionCmd) {         
		generator->SetParticlePosition(fPositionCmd->ConvertToDimensioned3Vector(newValue));
		BxOutputVertex::Get()->SetSpatialDist(0);	 
	} else if (cmd == fDirectionCmd){
		generator->SetParticleMomentumDirection(fDirectionCmd->ConvertTo3Vector(newValue));
	} else if (cmd == fEnergyCmd){
		generator->SetParticleEnergy(fEnergyCmd->ConvertToDimensionedDouble(newValue));
	} else if (cmd == fParticleCmd){
		generator->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(newValue));
	} else if (cmd == fSphereBulkCmd){
		generator->SetIfVolDist(true);
		generator->SetPosDisType("Volume");
		generator->SetPosDisShape("Sphere");
		generator->SetRadius(fSphereBulkCmd->ConvertToDimensionedDouble(newValue));	 
		BxOutputVertex::Get()->SetRagMin(0);
		BxOutputVertex::Get()->SetRagMax(fSphereBulkCmd->ConvertToDimensionedDouble(newValue));
	} else if (cmd == fSphereSurfCmd){
		generator->SetIfVolDist(true);
		generator->SetPosDisType("Surface");
		generator->SetPosDisShape("Sphere");
		generator->SetRadius(fSphereSurfCmd->ConvertToDimensionedDouble(newValue));	 
		BxOutputVertex::Get()->SetRagMin(fSphereBulkCmd->ConvertToDimensionedDouble(newValue));
		BxOutputVertex::Get()->SetRagMax(fSphereBulkCmd->ConvertToDimensionedDouble(newValue));
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
	} else if (cmd == fEnergyDisTypeCmd){ 
		BxLog(routine) << "G4Gun Energy distribution " <<newValue  << endlog ;
		generator->SetEnergyDistribution(true);
		generator->SetEnergyDisType(newValue);
	} else if (cmd == fSetEminCmd){ 
		BxLog(routine) << "G4Gun Emin = " <<newValue  << endlog ;
		generator->SetEmin(fSetEminCmd->ConvertToDimensionedDouble(newValue));
	} else if (cmd == fSetEmaxCmd){ 
		BxLog(routine) << "G4Gun Emax = " <<newValue  << endlog ;
		generator->SetEmax(fSetEmaxCmd->ConvertToDimensionedDouble(newValue));
	} else if (cmd == fSetAlphaCmd){ 
		BxLog(routine) << "G4Gun Alpha = " <<newValue  << endlog ;
		generator->SetAlpha(fSetAlphaCmd->ConvertToDouble(newValue));
	} else if (cmd == fSetTempCmd){ 
		BxLog(routine) << "G4Gun Temperature = " <<newValue  << endlog ;
		generator->SetTemp(fSetTempCmd->ConvertToDimensionedDouble(newValue));
	} else if (cmd == fSetEzeroCmd){ 
		BxLog(routine) << "G4Gun Ezero = " <<newValue  << endlog ;
		generator->SetEzero(fSetEzeroCmd->ConvertToDimensionedDouble(newValue));
	} else if (cmd == fSetGradientCmd){ 
		BxLog(routine) << "G4Gun Gradient = " <<newValue  << endlog ;
		generator->SetGradient(fSetGradientCmd->ConvertToDouble(newValue));
	} else if (cmd == fSetInterCeptCmd){ 
		BxLog(routine) << "G4Gun Intercept = " <<newValue  << endlog ;
		generator->SetInterCept(fSetInterCeptCmd->ConvertToDouble(newValue));
	}

}
