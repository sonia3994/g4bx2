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
#include "BxGeneratorSCSMessenger.hh"
#include "BxGeneratorSCS.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "BxReadParameters.hh"
#include "BxPropertyCollection.hh"
#include "BxOutputVertex.hh"
#include "G4SystemOfUnits.hh"
using namespace std;

BxGeneratorSCSMessenger::BxGeneratorSCSMessenger(BxGeneratorSCS* gen){
  generator = gen;

  fDirectory = new G4UIdirectory("/bx/generator/scs/");
  fDirectory->SetGuidance("Control of Bxscs event generator");

  fParticleCmd = new G4UIcmdWithAString("/bx/generator/scs/isotope",this);
  fParticleCmd->SetCandidates("c14 bi210 c14oleg");

  fPositionCmd = new G4UIcmdWith3VectorAndUnit("/bx/generator/scs/position",this);
  fPositionCmd->SetGuidance("Set the gun position");
  fPositionCmd->SetUnitCategory("Length");
  fPositionCmd->SetDefaultUnit("cm");
  fPositionCmd->SetUnitCandidates("mm cm m");
  
  fSphereBulkCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/scs/sphere_radius",this);
  fSphereBulkCmd->SetGuidance("Bulk radius");
  fSphereBulkCmd->SetUnitCategory("Length");
  fSphereBulkCmd->SetDefaultUnit("cm");
  fSphereBulkCmd->SetUnitCandidates("mm cm m");
  
  fSphereBulkRadMinCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/scs/sphere_radius_min",this);
  fSphereBulkRadMinCmd->SetGuidance("Bulk minimum radius");
  fSphereBulkRadMinCmd->SetGuidance("Default: 0");
  fSphereBulkRadMinCmd->SetUnitCategory("Length");
  fSphereBulkRadMinCmd->SetDefaultUnit("cm");
  fSphereBulkRadMinCmd->SetUnitCandidates("mm cm m");

  fSphereSurfCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/scs/surface_radius",this);
  fSphereSurfCmd->SetGuidance("Surface radius");
  fSphereSurfCmd->SetUnitCategory("Length");
  fSphereSurfCmd->SetDefaultUnit("cm");
  fSphereSurfCmd->SetUnitCandidates("mm cm m");

  fSphereCentreCmd = new G4UIcmdWith3VectorAndUnit("/bx/generator/scs/sphere_origin",this);
  fSphereCentreCmd->SetGuidance("Set the sphere position");
  fSphereCentreCmd->SetGuidance("Default: 0. 0. 0. cm");
  fSphereCentreCmd->SetUnitCategory("Length");
  fSphereCentreCmd->SetDefaultUnit("cm");
  fSphereCentreCmd->SetUnitCandidates("mm cm m");

  fVesselCmd = new G4UIcmdWithABool("/bx/generator/scs/vessel",this);
  fVesselCmd->SetGuidance("Set events on the vessel inner surface");

  fBulkCmd = new G4UIcmdWithABool("/bx/generator/scs/bulk",this);
  fBulkCmd->SetGuidance("Set events on the vessel inner surface");

  fBufferCmd = new G4UIcmdWithABool("/bx/generator/scs/buffer",this);
  fBufferCmd->SetGuidance("Set events on the vessel inner surface");
 
  fDirectionCmd = new G4UIcmdWith3Vector("/bx/generator/scs/direction",this);
  fDirectionCmd->SetGuidance("Set the gun direction");
  

  // confine to volume
  fConfineCmd = new G4UIcmdWithAString("/bx/generator/scs/confine",this);
  fConfineCmd->SetGuidance("Confine source to volume (NULL to unset).");
  fConfineCmd->SetGuidance("usage: confine VolName");
  fConfineCmd->SetParameterName("VolName",true,true);
  fConfineCmd->SetDefaultValue("NULL");


  fNumberOfBinsCmd = new G4UIcmdWithAnInteger("/bx/generator/scs/number_of_bins",this);
  fNumberOfBinsCmd->SetGuidance("default: 200");

  fReferenceCmd = new G4UIcmdWithAnInteger("/bx/generator/scs/c14reference",this) ;
  fReferenceCmd->SetGuidance("Set shape factor reference");
  
  fShapeFactorCmd = new G4UIcmdWithADouble("/bx/generator/scs/c14shapefactor",this);
 
  fShapeFactorTermCmd = new G4UIcmdWithAnInteger("/bx/generator/scs/c14shapefactorterm",this) ;

}


BxGeneratorSCSMessenger::~BxGeneratorSCSMessenger() {

  delete fDirectory;
  delete fPositionCmd; 
  delete fDirectionCmd;
  delete fParticleCmd; 
  delete fSphereBulkCmd; 
  delete fSphereSurfCmd;
  delete fSphereCentreCmd;
  delete fVesselCmd;
  delete fBulkCmd;
  delete fBufferCmd;
  delete fConfineCmd;
  
  delete fNumberOfBinsCmd ;
  delete fShapeFactorCmd;
  delete fReferenceCmd ;
  delete fShapeFactorTermCmd ;
}


void BxGeneratorSCSMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) { 
   if (cmd == fPositionCmd) {         
	 generator->SetParticlePosition(fPositionCmd->ConvertToDimensioned3Vector(newValue));
         BxOutputVertex::Get()->SetSpatialDist(0);	 
   } else if (cmd == fDirectionCmd){
   	 generator->SetParticleMomentumDirection(fDirectionCmd->ConvertTo3Vector(newValue));
   } else if (cmd == fParticleCmd){
	 BxLog(routine) << "Isotope: " << newValue << endlog ;
         generator->SetCrossSection(newValue);
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
   } else if (cmd == fReferenceCmd){ 
	   if(fReferenceCmd->ConvertToInt(newValue) == 1) {
		   generator->SetShapeFactorTerm(1);
		   generator->SetShapeFactor(0.38);
		   BxLog(routine) << "Shape Factor: 1 - 0.38 * E"  << endlog ;
		   BxLog(routine) << "From reference: F. P. Calaprice and B. R. Holstein, Nucl. Phys. A 273, 301 (1976)"  << endlog ;	 
	   } else if(fReferenceCmd->ConvertToInt(newValue) == 2) {
		   generator->SetShapeFactorTerm(1);
		   generator->SetShapeFactor(0.72);
		   BxLog(routine) << "Shape Factor: 1 - 0.72 * E"  << endlog ;
		   BxLog(routine) << "G. Alimonti, G. Angloher, C. Arpesell, et al., Phys. Lett. B 422, 349 (1998)."  << endlog ;	 
	   } else if(fReferenceCmd->ConvertToInt(newValue) == 3) {
		   generator->SetShapeFactorTerm(1);
		   generator->SetShapeFactor(1.179);
		   BxLog(routine) << "Shape Factor: 1 - 1.179 * E"  << endlog ;
		   BxLog(routine) << "H. Genz, G. Kuhner, A. Richter, and H. Behrens, Z. Phys. A 341, 9 (1991)."  << endlog ;	 
	   } else if(fReferenceCmd->ConvertToInt(newValue) == 4) {
		   generator->SetShapeFactorTerm(1);
		   generator->SetShapeFactor(0.37);
		   BxLog(routine) << "Shape Factor: 1 - 0.37 * E"  << endlog ;
		   BxLog(routine) << "A. Garcia and B. A. Brown, Phys. Rev. C 52, 3416 (1995)."  << endlog ;	 
	   } else if(fReferenceCmd->ConvertToInt(newValue) == 5) {
		   generator->SetShapeFactorTerm(1);
		   generator->SetShapeFactor(0.45);
		   BxLog(routine) << "Shape Factor: 1 - 0.45 * E"  << endlog ;
		   BxLog(routine) << "F. E. Weitfeldt, E. B. Norman, Y. D. Chan, et al., Phys. Rev. C 52, 1028 (1995)."  << endlog ;	 
	   } else if(fReferenceCmd->ConvertToInt(newValue) == 6) {
		   generator->SetShapeFactorTerm(2);
		   generator->SetShapeFactor(1.1);
		   BxLog(routine) << "Shape Factor: 1 + 1.1 * (E0 - E)"  << endlog ;
		   BxLog(routine) << "B. Sur, E. B. Norman, K. T. Lesko, et al., Phys. Rev. Lett. 66, 2444 (1991)."  << endlog ;	 
	   } else if(fReferenceCmd->ConvertToInt(newValue) == 7) {
		   generator->SetShapeFactorTerm(2);
		   generator->SetShapeFactor(1.24);
		   BxLog(routine) << "Shape Factor: 1 + 1.24 * (E0 - E)"  << endlog ;
		   BxLog(routine) << "V. V. Kuzminov and N. J. Osetrova, Exp. Res. Meth. And Fac. 63 (2000)  1365"  << endlog ;	 
	   } else if(fReferenceCmd->ConvertToInt(newValue) == 8) {
		   generator->SetShapeFactorTerm(3);
		   BxLog(routine) << "Shape Factor: 1 - 4.67*E + 3/E +2*E*E"  << endlog ;
		   BxLog(routine) << "Ch. Sonntag, H. Rebel, B. Ribbat, et al., Lett. Nuovo Cimento 4, 717 (1970)."  << endlog ;	 
	   } else if(fReferenceCmd->ConvertToInt(newValue) == 9) {
		   generator->SetShapeFactorTerm(4);
		   generator->SetShapeFactor(1.24);
		   generator->SetNumberOfBins(401);
		   BxPropertyCollection::Get()->ReadC14CorrectionFactor();
		   BxLog(routine) << "Shape Factor: alpha(E) (1 + 1.24 * (E0 - E_corrected))"  << endlog ;
		   BxLog(routine) << "From V. V. Kuzminov and N. J. Osetrova, Exp. Res. Meth. And Fac. 63 (2000)  1365 with corrections provided by Oleg"  << endlog ;	 
	   } 
   } else if (cmd == fShapeFactorCmd){
	   generator->SetShapeFactor(fShapeFactorCmd->ConvertToDouble(newValue));
	   BxLog(routine) << "Shape Factor: " << newValue  << endlog ;
   } else if (cmd == fShapeFactorTermCmd){
	   generator->SetShapeFactorTerm(fShapeFactorTermCmd->ConvertToInt(newValue));
	   if(fShapeFactorTermCmd->ConvertToInt(newValue) == 1) {
		   BxLog(routine) << "Shape Factor Term: 1 - SF * E"  << endlog ;
	   }  else if(fShapeFactorTermCmd->ConvertToInt(newValue) == 2) {
		   BxLog(routine) << "Shape Factor Term: 1 + SF * ( E0 - E )"  << endlog ;	 
	   }  else {
		   BxLog(error) << "No Shape Factor Term corresponding to " << newValue   << endlog ;	 
		   BxLog(fatal) << endl ;
	   } 
   }  else if (cmd == fNumberOfBinsCmd){
	   generator->SetNumberOfBins(fNumberOfBinsCmd->ConvertToInt(newValue));
	   BxLog(routine) << "Cross Section has been binned with " << newValue << " steps"  << endlog ;
   }
}
