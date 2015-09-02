// --------------------------------------------------------------------------//
/** 
 * AUTHOR: Igor Machulin
 * CONTACT: machulin@lngs.infn.it
* Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "BxGeneratorLi9He8Messenger.hh"
#include "G4UIcmdWithADouble.hh"
#include "BxGeneratorLi9He8.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "BxReadParameters.hh"
#include "BxOutputVertex.hh"
#include "BxLogger.hh"
#include "G4SystemOfUnits.hh"
using namespace std;

BxGeneratorLi9He8Messenger::BxGeneratorLi9He8Messenger(BxGeneratorLi9He8* gen){
  generator = gen;
  fDirectory = new G4UIdirectory("/bx/generator/Li9He8/");
  fDirectory->SetGuidance("Control of BxLi9He8 event generator");

  fPositionCmd = new G4UIcmdWith3VectorAndUnit("/bx/generator/Li9He8/position",this);
  fPositionCmd->SetGuidance("Set the gun position");
  fPositionCmd->SetUnitCategory("Length");
  fPositionCmd->SetDefaultUnit("cm");
  fPositionCmd->SetUnitCandidates("mm cm m");
  
  fSphereBulkCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/Li9He8/sphere_radius",this);
  fSphereBulkCmd->SetGuidance("Bulk radius");
  fSphereBulkCmd->SetUnitCategory("Length");
  fSphereBulkCmd->SetDefaultUnit("cm");
  fSphereBulkCmd->SetUnitCandidates("mm cm m");
  
  fSphereBulkRadMinCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/Li9He8/sphere_radius_min",this);
  fSphereBulkRadMinCmd->SetGuidance("Bulk minimum radius");
  fSphereBulkRadMinCmd->SetGuidance("Default: 0");
  fSphereBulkRadMinCmd->SetUnitCategory("Length");
  fSphereBulkRadMinCmd->SetDefaultUnit("cm");
  fSphereBulkRadMinCmd->SetUnitCandidates("mm cm m");

  fSphereSurfCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/Li9He8/surface_radius",this);
  fSphereSurfCmd->SetGuidance("Surface radius");
  fSphereSurfCmd->SetUnitCategory("Length");
  fSphereSurfCmd->SetDefaultUnit("cm");
  fSphereSurfCmd->SetUnitCandidates("mm cm m");

  fSphereCentreCmd = new G4UIcmdWith3VectorAndUnit("/bx/generator/Li9He8/sphere_origin",this);
  fSphereCentreCmd->SetGuidance("Set the sphere position");
  fSphereCentreCmd->SetGuidance("Default: 0. 0. 0. cm");
  fSphereCentreCmd->SetUnitCategory("Length");
  fSphereCentreCmd->SetDefaultUnit("cm");
  fSphereCentreCmd->SetUnitCandidates("mm cm m");

  fVesselCmd = new G4UIcmdWithABool("/bx/generator/Li9He8/vessel",this);
  fVesselCmd->SetGuidance("Set events on the vessel inner surface");

  fBulkCmd = new G4UIcmdWithABool("/bx/generator/Li9He8/bulk",this);
  fBulkCmd->SetGuidance("Set events on the vessel inner surface");

  fBufferCmd = new G4UIcmdWithABool("/bx/generator/Li9He8/buffer",this);
  fBufferCmd->SetGuidance("Set events on the vessel inner surface");
 

  // confine to volume
  fConfineCmd = new G4UIcmdWithAString("/bx/generator/Li9He8/confine",this);
  fConfineCmd->SetGuidance("Confine source to volume (NULL to unset).");
  fConfineCmd->SetGuidance("usage: confine VolName");
  fConfineCmd->SetParameterName("VolName",true,true);
  fConfineCmd->SetDefaultValue("NULL");
  
  fNeutrinoCmd = new G4UIcmdWithAString("/bx/generator/Li9He8/source",this);
  G4String candidates = "He8 Li9";  
  fNeutrinoCmd->SetCandidates(candidates);
 
}


BxGeneratorLi9He8Messenger::~BxGeneratorLi9He8Messenger() {

  delete fDirectory;
  delete fPositionCmd; 
  delete fDirectionCmd;
  delete fSphereBulkCmd; 
  delete fSphereSurfCmd;
  delete fSphereCentreCmd;
  delete fVesselCmd;
  delete fBulkCmd;
  delete fBufferCmd;
  delete fConfineCmd;
  delete fNeutrinoCmd;
  delete fNumberOfStepCmd;
}


void BxGeneratorLi9He8Messenger::SetNewValue(G4UIcommand* cmd, G4String newValue) { 
   if (cmd == fPositionCmd) {         
	 generator->SetParticlePosition(fPositionCmd->ConvertToDimensioned3Vector(newValue));
         BxOutputVertex::Get()->SetSpatialDist(0);	 
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
   } else if (cmd == fNeutrinoCmd){
	   BxLog(routine) << "Li9He8 source: " << newValue << endlog ;
	   if(newValue == "He8") {
		   generator->SetNeutrinoType(BxGeneratorLi9He8::He8N);
	   } else if(newValue == "Li9") {
		   generator->SetNeutrinoType(BxGeneratorLi9He8::Li9N);
	   } 
   }   
}
