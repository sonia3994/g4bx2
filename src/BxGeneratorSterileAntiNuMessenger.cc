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

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "BxGeneratorSterileAntiNu.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "BxReadParameters.hh"
#include "BxOutputVertex.hh"
#include "BxLogger.hh"
#include "G4SystemOfUnits.hh"
using namespace std;

BxGeneratorSterileAntiNuMessenger::BxGeneratorSterileAntiNuMessenger(BxGeneratorSterileAntiNu* gen){
  generator = gen;
  fDirectory = new G4UIdirectory("/bx/generator/SterileAntiNu/");
  fDirectory->SetGuidance("Control of Bx SterileAntiNu event generator");


  fPositionCmd = new G4UIcmdWith3VectorAndUnit("/bx/generator/SterileAntiNu/position",this);
  fPositionCmd->SetGuidance("Set the gun position");
  fPositionCmd->SetUnitCategory("Length");
  fPositionCmd->SetDefaultUnit("cm");
  fPositionCmd->SetUnitCandidates("mm cm m");
  
  fSphereBulkCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/SterileAntiNu/sphere_radius",this);
  fSphereBulkCmd->SetGuidance("Bulk radius");
  fSphereBulkCmd->SetUnitCategory("Length");
  fSphereBulkCmd->SetDefaultUnit("cm");
  fSphereBulkCmd->SetUnitCandidates("mm cm m");
  
  fSphereBulkRadMinCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/SterileAntiNu/sphere_radius_min",this);
  fSphereBulkRadMinCmd->SetGuidance("Bulk minimum radius");
  fSphereBulkRadMinCmd->SetGuidance("Default: 0");
  fSphereBulkRadMinCmd->SetUnitCategory("Length");
  fSphereBulkRadMinCmd->SetDefaultUnit("cm");
  fSphereBulkRadMinCmd->SetUnitCandidates("mm cm m");

  fSphereSurfCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/SterileAntiNu/surface_radius",this);
  fSphereSurfCmd->SetGuidance("Surface radius");
  fSphereSurfCmd->SetUnitCategory("Length");
  fSphereSurfCmd->SetDefaultUnit("cm");
  fSphereSurfCmd->SetUnitCandidates("mm cm m");

  fSphereCentreCmd = new G4UIcmdWith3VectorAndUnit("/bx/generator/SterileAntiNu/sphere_origin",this);
  fSphereCentreCmd->SetGuidance("Set the sphere position");
  fSphereCentreCmd->SetGuidance("Default: 0. 0. 0. cm");
  fSphereCentreCmd->SetUnitCategory("Length");
  fSphereCentreCmd->SetDefaultUnit("cm");
  fSphereCentreCmd->SetUnitCandidates("mm cm m");

  fVesselCmd = new G4UIcmdWithABool("/bx/generator/SterileAntiNu/vessel",this);
  fVesselCmd->SetGuidance("Set events on the vessel inner surface");

  fBulkCmd = new G4UIcmdWithABool("/bx/generator/SterileAntiNu/bulk",this);
  fBulkCmd->SetGuidance("Set events on the vessel inner surface");

  fBufferCmd = new G4UIcmdWithABool("/bx/generator/SterileAntiNu/buffer",this);
  fBufferCmd->SetGuidance("Set events on the vessel inner surface");
 
  fDirectionCmd = new G4UIcmdWith3Vector("/bx/generator/SterileAntiNu/direction",this);
  fDirectionCmd->SetGuidance("Set the gun direction");
  
  fEnergyCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/SterileAntiNu/energy",this);
  fEnergyCmd->SetGuidance("Set the gun energy");
  fEnergyCmd->SetUnitCategory("Energy");
  fEnergyCmd->SetDefaultUnit("MeV");
  fEnergyCmd->SetUnitCandidates("eV keV MeV GeV");

  fEnergyDisTypeCmd = new G4UIcmdWithAString("/bx/generator/SterileAntiNu/dist_energy",this);
  G4String candidates = "Lin Pow Exp Gauss Brem BBody Cdg";  
  fEnergyDisTypeCmd->SetCandidates(candidates);
 
  fSetEminCmd      = new G4UIcmdWithADoubleAndUnit("/bx/generator/SterileAntiNu/emin",this);
  fSetEminCmd->SetDefaultUnit("MeV");
  fSetEminCmd->SetUnitCandidates("eV keV MeV GeV");
  
  fSetEmaxCmd      = new G4UIcmdWithADoubleAndUnit("/bx/generator/SterileAntiNu/emax",this);
  fSetEmaxCmd->SetDefaultUnit("MeV");
  fSetEmaxCmd->SetUnitCandidates("eV keV MeV GeV");
  
  fSetAlphaCmd     = new G4UIcmdWithADouble("/bx/generator/SterileAntiNu/alpha",this);
  
  fSetTempCmd      = new G4UIcmdWithADoubleAndUnit("/bx/generator/SterileAntiNu/temp",this);
  fSetTempCmd->SetDefaultUnit("kelvin");
  fSetTempCmd->SetUnitCandidates("kelvin");
  
  fSetEzeroCmd     = new G4UIcmdWithADoubleAndUnit("/bx/generator/SterileAntiNu/ezero",this);
  fSetEzeroCmd->SetDefaultUnit("MeV");
  fSetEzeroCmd->SetUnitCandidates("eV keV MeV GeV");
  
  fSetGradientCmd  = new G4UIcmdWithADouble("/bx/generator/SterileAntiNu/gradient",this);
  
  fSetInterCeptCmd = new G4UIcmdWithADouble("/bx/generator/SterileAntiNu/intercept",this);
  
 
  fParticleCmd = new G4UIcmdWithAString("/bx/generator/SterileAntiNu/particle",this);
  fParticleCmd->SetGuidance("Set the gun type");

  // confine to volume
  fConfineCmd = new G4UIcmdWithAString("/bx/generator/SterileAntiNu/confine",this);
  fConfineCmd->SetGuidance("Confine source to volume (NULL to unset).");
  fConfineCmd->SetGuidance("usage: confine VolName");
  fConfineCmd->SetParameterName("VolName",true,true);
  fConfineCmd->SetDefaultValue("NULL");

  /// for oscillation and spectral deformation (shape factor):
  fDM2Cmd  =  new G4UIcmdWithADouble("/bx/generator/SterileAntiNu/deltam2",this);

  fSIN2T2Cmd =  new G4UIcmdWithADouble("/bx/generator/SterileAntiNu/sin2t2",this);

  fCOEFFICIENTACmd =  new G4UIcmdWithADouble("/bx/generator/SterileAntiNu/shapeFactorCoefficientA",this);

  fCOEFFICIENTBCmd =  new G4UIcmdWithADouble("/bx/generator/SterileAntiNu/shapeFactorCoefficientB",this);

  fCOEFFICIENTCCmd =  new G4UIcmdWithADouble("/bx/generator/SterileAntiNu/shapeFactorCoefficientC",this);



}


BxGeneratorSterileAntiNuMessenger::~BxGeneratorSterileAntiNuMessenger() {

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
  /// new:
  delete fDM2Cmd;
  delete fSIN2T2Cmd;
  delete fCOEFFICIENTACmd;
  delete fCOEFFICIENTBCmd;
  delete fCOEFFICIENTCCmd;
}


void BxGeneratorSterileAntiNuMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) { 
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
	   //}
	   /// Neu:
} else if (cmd == fDM2Cmd){
	generator->SetDM2(fDM2Cmd->ConvertToDouble(newValue));   
} else if (cmd == fSIN2T2Cmd){
	generator->SetSIN2T2(fSIN2T2Cmd->ConvertToDouble(newValue));
} else if (cmd == fCOEFFICIENTACmd){ // FIXME   
	generator->SetCOEFFICIENTA(fCOEFFICIENTACmd->ConvertToDouble(newValue));
} else if (cmd == fCOEFFICIENTBCmd){ 
	generator->SetCOEFFICIENTB(fCOEFFICIENTBCmd->ConvertToDouble(newValue));
} else if (cmd == fCOEFFICIENTCCmd){ 
	generator->SetCOEFFICIENTC(fCOEFFICIENTCCmd->ConvertToDouble(newValue));
}



       
}
/*
 * $Log: BxGeneratorSterileAntiNuMessenger.cc,v $
 * Revision 1.3  2015/04/22 13:05:48  acaminata
 * bug-fixing in generator confine command
 *
 * Revision 1.2  2015/02/25 18:17:18  acaminata
 * Shape factor added
 *
 * Revision 1.1  2015/02/12 14:31:35  acaminata
 * Sterile generator added
 *
 * Revision 1.12  2009-01-13 14:06:54  dfranco
 * Fixed a bug in the positioning of the centre of the generation spherical distribution
 *
 * Revision 1.11  2007-11-12 12:09:41  dfranco
 * added to g4gun, the following  energy distributions:
 * Lin (linear), Pow (power-law), Exp (exponential), Gauss (gaussian),
 * Brem (bremsstrahlung), BBody (black-body), Cdg (cosmic diffuse gamma-ray)
 *
 * Revision 1.10  2007-05-07 12:57:12  dfranco
 * Parameter tuning:
 * photonyield = 11600
 * kB  = 0.0082
 * kB2 = 0
 *
 * Revision 1.9  2007-04-26 08:51:25  dfranco
 * Added the structure for the new ouput format. It is not yet active!!
 *
 * Revision 1.8  2007-03-22 14:48:55  dfranco
 * Stable version
 *
 */
