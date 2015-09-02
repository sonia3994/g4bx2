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
// *********--***********************************************************
//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "BxGeneratorMultiEventMessenger.hh"
#include "BxGeneratorMultiEvent.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "BxReadParameters.hh"
#include "BxLogger.hh"
#include "BxOutputVertex.hh"
#include "G4SystemOfUnits.hh"
#include <stdio.h>

using namespace std;

BxGeneratorMultiEventMessenger::BxGeneratorMultiEventMessenger(BxGeneratorMultiEvent* gen){
  generator = gen;
  fDirectory = new G4UIdirectory("/bx/generator/multi/");
  fDirectory->SetGuidance("Control of BxMultiEvent generator");
  
  fPDGCmd  = new G4UIcmdWithAString("/bx/generator/multi/event",this);
 

  fPositionCmd = new G4UIcmdWith3VectorAndUnit("/bx/generator/multi/position",this);
  fPositionCmd->SetGuidance("Set the gun position");
  fPositionCmd->SetUnitCategory("Length");
  fPositionCmd->SetDefaultUnit("cm");
  fPositionCmd->SetUnitCandidates("mm cm m");
  
  fSphereBulkCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/multi/sphere_radius",this);
  fSphereBulkCmd->SetGuidance("Bulk radius");
  fSphereBulkCmd->SetUnitCategory("Length");
  fSphereBulkCmd->SetDefaultUnit("cm");
  fSphereBulkCmd->SetUnitCandidates("mm cm m");
  
  fSphereBulkRadMinCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/multi/sphere_radius_min",this);
  fSphereBulkRadMinCmd->SetGuidance("Bulk minimum radius");
  fSphereBulkRadMinCmd->SetGuidance("Default: 0");
  fSphereBulkRadMinCmd->SetUnitCategory("Length");
  fSphereBulkRadMinCmd->SetDefaultUnit("cm");
  fSphereBulkRadMinCmd->SetUnitCandidates("mm cm m");

  fSphereSurfCmd = new G4UIcmdWithADoubleAndUnit("/bx/generator/multi/surface_radius",this);
  fSphereSurfCmd->SetGuidance("Surface radius");
  fSphereSurfCmd->SetUnitCategory("Length");
  fSphereSurfCmd->SetDefaultUnit("cm");
  fSphereSurfCmd->SetUnitCandidates("mm cm m");

  fSphereCentreCmd = new G4UIcmdWith3VectorAndUnit("/bx/generator/multi/sphere_origin",this);
  fSphereCentreCmd->SetGuidance("Set the sphere position");
  fSphereCentreCmd->SetGuidance("Default: 0. 0. 0. cm");
  fSphereCentreCmd->SetUnitCategory("Length");
  fSphereCentreCmd->SetDefaultUnit("cm");
  fSphereCentreCmd->SetUnitCandidates("mm cm m");

  fVesselCmd = new G4UIcmdWithABool("/bx/generator/multi/vessel",this);
  fVesselCmd->SetGuidance("Set events on the vessel inner surface");

  fBulkCmd = new G4UIcmdWithABool("/bx/generator/multi/bulk",this);
  fBulkCmd->SetGuidance("Set events on the vessel inner surface");

  fBufferCmd = new G4UIcmdWithABool("/bx/generator/multi/buffer",this);
  fBufferCmd->SetGuidance("Set events on the vessel inner surface");
 
  // confine to volume
  fConfineCmd = new G4UIcmdWithAString("/bx/generator/multi/confine",this);
  fConfineCmd->SetGuidance("Confine source to volume (NULL to unset).");
  fConfineCmd->SetGuidance("usage: confine VolName");
  fConfineCmd->SetParameterName("VolName",true,true);
  fConfineCmd->SetDefaultValue("NULL");
}


BxGeneratorMultiEventMessenger::~BxGeneratorMultiEventMessenger() {

  delete fDirectory;
  delete fPositionCmd; 
  delete fSphereBulkCmd; 
  delete fSphereSurfCmd;
  delete fSphereCentreCmd;
  delete fVesselCmd;
  delete fBulkCmd;
  delete fBufferCmd;
  delete fOnlyGammaCmd;
  delete fSkipEventCmd;
  delete fPDGCmd;
  delete fVesselCmd;
  delete fBulkCmd;
  delete fBufferCmd;
  delete fConfineCmd;

}


void BxGeneratorMultiEventMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) { 
   if (cmd == fPositionCmd) { 
	 generator->SetParticlePosition(fPositionCmd->ConvertToDimensioned3Vector(newValue));
         BxOutputVertex::Get()->SetSpatialDist(0);	 
   }  else if (cmd == fPDGCmd){
         ConvertFromString(newValue);
	 generator->SetListOfParticles(fMultiParticle);
   }  else if (cmd == fSphereBulkCmd){
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
   } else if (cmd == fConfineCmd){
	   generator->SetIfVolDist(true);
	   generator->ConfineSourceToVolume(newValue);
   } else if (cmd == fVesselCmd){
	   generator->SetIfVolDist(true);
	   generator->SetPosDisType("Surface");
	   generator->SetPosDisShape("Sphere");
	   generator->SetRadius(BxReadParameters::Get()->GetZoneIExternalRadius() - 0.001*mm); 
	   BxOutputVertex::Get()->SetRagMin(BxReadParameters::Get()->GetZoneIExternalRadius() - 0.001*mm);
	   BxOutputVertex::Get()->SetRagMax(BxReadParameters::Get()->GetZoneIExternalRadius() - 0.001*mm);
	   BxOutputVertex::Get()->SetSpatialDist(1);	 
   } else if (cmd == fBulkCmd){
	   generator->SetIfVolDist(true);
	   generator->SetPosDisType("Volume");
	   generator->SetPosDisShape("Sphere");
	   generator->SetRadius(BxReadParameters::Get()->GetZoneIExternalRadius() - 0.001*mm); 
	   BxOutputVertex::Get()->SetRagMin(0);
	   BxOutputVertex::Get()->SetRagMax(BxReadParameters::Get()->GetZoneIExternalRadius() - 0.001*mm);
	   BxOutputVertex::Get()->SetSpatialDist(1);	 
   } else if (cmd == fBufferCmd){
	   generator->SetIfVolDist(true);
	   generator->SetPosDisType("Volume");
	   generator->SetPosDisShape("Sphere");	 
	   generator->SetRadius0(BxReadParameters::Get()->GetZoneIExternalRadius()); 
	   generator->SetRadius (BxReadParameters::Get()->GetZoneIIIExternalRadius()); 
	   BxOutputVertex::Get()->SetRagMin(BxReadParameters::Get()->GetZoneIExternalRadius());
	   BxOutputVertex::Get()->SetRagMax(BxReadParameters::Get()->GetZoneIIIExternalRadius());
	   BxOutputVertex::Get()->SetSpatialDist(1);	 
	   BxLog(routine) << BxReadParameters::Get()->GetZoneIExternalRadius() << " " 
		   << BxReadParameters::Get()->GetZoneIIIExternalRadius() << endlog ;
   }
}


MultiEvent& BxGeneratorMultiEventMessenger::ConvertFromString(string text) {  
  sscanf (text.c_str(),"%d %f %d %f %f",
          &fMultiParticle.Id,
          &fMultiParticle.BRTOT,
          &fMultiParticle.PDG,
          &fMultiParticle.Energy,	  
          &fMultiParticle.BR
	 );
  return fMultiParticle;

}
