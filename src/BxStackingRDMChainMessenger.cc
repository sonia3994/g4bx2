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
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "BxStackingRDMChainMessenger.hh"
#include "BxStackingRDMChain.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "BxReadParameters.hh"
#include "BxLogger.hh"
#include "G4SystemOfUnits.hh"
using namespace std;

BxStackingRDMChainMessenger::BxStackingRDMChainMessenger(BxStackingRDMChain* stack){
  stacking = stack;
  fDirectory = new G4UIdirectory("/bx/stack/rdmchain/");
  fDirectory->SetGuidance("Control of BxGeneb event generator");
  
  fVesselSupDef = new G4UIcmdWithAString("/bx/stack/rdmchain/VesselSupDef", this);
 
  fVesselThickness = new G4UIcmdWithADoubleAndUnit("/bx/stack/rdmchain/VesselThickness", this);
  fVesselThickness->SetDefaultUnit("mm");
  fVesselThickness->SetUnitCandidates("mm cm m");
  
  fRhoCut = new G4UIcmdWithADoubleAndUnit("/bx/stack/rdmchain/RhoCut", this);
  fRhoCut->SetDefaultUnit("mm");
  fRhoCut->SetUnitCandidates("mm cm m");
  
  fZCut = new G4UIcmdWithADoubleAndUnit("/bx/stack/rdmchain/ZCut", this);
  fZCut->SetDefaultUnit("mm");
  fZCut->SetUnitCandidates("mm cm m");
  
  fBufferDef = new G4UIcmdWithAString("/bx/stack/rdmchain/BufferDef", this);
  
  fLifeTimeCmd = new G4UIcmdWithADoubleAndUnit("/bx/stack/rdmchain/maxlifetime",this);
  fLifeTimeCmd->SetGuidance("Set the max life time acceptable for a decay");
  fLifeTimeCmd->SetDefaultUnit("s");
  fLifeTimeCmd->SetUnitCandidates("ps ns mus ms s");
  
}


BxStackingRDMChainMessenger::~BxStackingRDMChainMessenger() {
  delete fVesselSupDef;
  delete fBufferDef;
  delete fDirectory;
  delete fLifeTimeCmd;
  delete fVesselThickness;
  delete fRhoCut;
  delete fZCut;
}


void BxStackingRDMChainMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) { 
   if (cmd == fLifeTimeCmd){
         stacking->SetMaxLifeTime(fLifeTimeCmd->GetNewDoubleValue(newValue));
	 BxLog(routine) << "Max Decay Mean Life < " << newValue <<  endlog ;
   //} else if(cmd == fVesselSupDef) {
	    //BxReadParameters::Get()->SetVesselFileName(newValue);
	    //BxReadParameters::Get()->ReadVesselGeometry();
	    //stacking->SetfSupDef(true);
	    //BxLog(routine) << "Deformation surface active?  " << stacking->GetfSupDef() <<  endlog; 
	    //BxLog(routine) << "Deformed vessel filename: " << newValue  <<  endlog; 
    } else if(cmd == fVesselThickness) {
	    stacking->SetfVesThick(fVesselThickness->GetNewDoubleValue(newValue)); 
	    BxLog(routine) << "Deformed vessel thickness: " << fVesselThickness->GetNewDoubleValue(newValue)/mm  << " mm"   <<  endlog;
    } else if(cmd == fRhoCut) {
	    stacking->SetfRhoCut(fRhoCut->GetNewDoubleValue(newValue)); 
	    stacking->SetffRhocut(true);
	    BxLog(routine) << "Setting rho < " << fRhoCut->GetNewDoubleValue(newValue)/mm  << " mm"   <<  endlog;
    } else if(cmd == fZCut) {
	    stacking->SetfZCut(fZCut->GetNewDoubleValue(newValue)); 
	    stacking->SetffZcut(true);
	    BxLog(routine) << "Setting z < " << fZCut->GetNewDoubleValue(newValue)/mm  << " mm"   <<  endlog;
    //} else if(cmd == fBufferDef) {
//	    BxReadParameters::Get()->SetVesselFileName(newValue);
//	    BxReadParameters::Get()->ReadVesselGeometry();
//	    stacking->SetfBufDef(true);
//	    BxLog(routine) << "Deformation buffer active?  " << stacking->GetfBufDef() <<  endlog; 
//	    BxLog(routine) << "Deformed vessel filename: " << newValue  <<  endlog; 
    }
}
