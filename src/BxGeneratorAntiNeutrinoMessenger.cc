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
#include "BxGeneratorAntiNeutrinoMessenger.hh"
#include "G4UIcmdWithADouble.hh"
#include "BxGeneratorAntiNeutrino.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "BxReadParameters.hh"
#include "BxOutputVertex.hh"
#include "BxLogger.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
using namespace std;

BxGeneratorAntiNeutrinoMessenger::BxGeneratorAntiNeutrinoMessenger(BxGeneratorAntiNeutrino* gen){
  generator = gen;
  fDirectory = new G4UIdirectory("/bx/generator/AntiNu/");
  fDirectory->SetGuidance("Control of BxAntiNeutrino event generator");


  fPositionCmd = new G4UIcmdWith3VectorAndUnit("/bx/generator/AntiNu/position",this);
  fPositionCmd->SetGuidance("Set the gun position");
  fPositionCmd->SetUnitCategory("Length");
  fPositionCmd->SetDefaultUnit("cm");
  fPositionCmd->SetUnitCandidates("mm cm m");
  
  fSphereBulkCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/AntiNu/sphere_radius",this);
  fSphereBulkCmd->SetGuidance("Bulk radius");
  fSphereBulkCmd->SetUnitCategory("Length");
  fSphereBulkCmd->SetDefaultUnit("cm");
  fSphereBulkCmd->SetUnitCandidates("mm cm m");
  
  fSphereBulkRadMinCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/AntiNu/sphere_radius_min",this);
  fSphereBulkRadMinCmd->SetGuidance("Bulk minimum radius");
  fSphereBulkRadMinCmd->SetGuidance("Default: 0");
  fSphereBulkRadMinCmd->SetUnitCategory("Length");
  fSphereBulkRadMinCmd->SetDefaultUnit("cm");
  fSphereBulkRadMinCmd->SetUnitCandidates("mm cm m");

  fSphereSurfCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/AntiNu/surface_radius",this);
  fSphereSurfCmd->SetGuidance("Surface radius");
  fSphereSurfCmd->SetUnitCategory("Length");
  fSphereSurfCmd->SetDefaultUnit("cm");
  fSphereSurfCmd->SetUnitCandidates("mm cm m");

  fSphereCentreCmd = new G4UIcmdWith3VectorAndUnit("/bx/generator/AntiNu/sphere_origin",this);
  fSphereCentreCmd->SetGuidance("Set the sphere position");
  fSphereCentreCmd->SetGuidance("Default: 0. 0. 0. cm");
  fSphereCentreCmd->SetUnitCategory("Length");
  fSphereCentreCmd->SetDefaultUnit("cm");
  fSphereCentreCmd->SetUnitCandidates("mm cm m");

  fVesselCmd = new G4UIcmdWithABool("/bx/generator/AntiNu/vessel",this);
  fVesselCmd->SetGuidance("Set events on the vessel inner surface");

  fBulkCmd = new G4UIcmdWithABool("/bx/generator/AntiNu/bulk",this);
  fBulkCmd->SetGuidance("Set events on the vessel inner surface");

  fBufferCmd = new G4UIcmdWithABool("/bx/generator/AntiNu/buffer",this);
  fBufferCmd->SetGuidance("Set events on the vessel inner surface");
 

  // confine to volume
  fConfineCmd = new G4UIcmdWithAString("/bx/generator/AntiNu/confine",this);
  fConfineCmd->SetGuidance("Confine source to volume (NULL to unset).");
  fConfineCmd->SetGuidance("usage: confine VolName");
  fConfineCmd->SetParameterName("VolName",true,true);
  fConfineCmd->SetDefaultValue("NULL");
  
  fNeutrinoCmd = new G4UIcmdWithAString("/bx/generator/AntiNu/source",this);
  G4String candidates = "Geo Reactor GeoPositron ReactorPositron neutron U Th Sun Flat";  
  fNeutrinoCmd->SetCandidates(candidates);
 
}


BxGeneratorAntiNeutrinoMessenger::~BxGeneratorAntiNeutrinoMessenger() {

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


void BxGeneratorAntiNeutrinoMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) { 
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
	   BxLog(routine) << "AntiNeutrino source: " << newValue << endlog ;
	   if(newValue == "Geo") {
		   generator->SetNeutrinoType(BxGeneratorAntiNeutrino::GeoN);
	   } else if(newValue == "Reactor") {
		   generator->SetNeutrinoType(BxGeneratorAntiNeutrino::ReactorN);
	   } else if(newValue == "GeoPositron") {
		   generator->SetNeutrinoType(BxGeneratorAntiNeutrino::GeoPositroN);
	   } else if(newValue == "ReactorPositron") {
		   generator->SetNeutrinoType(BxGeneratorAntiNeutrino::ReactorPositroN);
	   } else if(newValue == "neutron") {
		   generator->SetNeutrinoType(BxGeneratorAntiNeutrino::NeutroN);
	   } else if(newValue == "U") {
		   generator->SetNeutrinoType(BxGeneratorAntiNeutrino::UrN);
	   } else if(newValue == "Th") {
		   generator->SetNeutrinoType(BxGeneratorAntiNeutrino::ThN);
	   } else if(newValue == "Sun") {
		   generator->SetNeutrinoType(BxGeneratorAntiNeutrino::SunN);
	   } else if(newValue == "Flat") {
		   generator->SetNeutrinoType(BxGeneratorAntiNeutrino::FlatN);
	   } 
   }   
}
