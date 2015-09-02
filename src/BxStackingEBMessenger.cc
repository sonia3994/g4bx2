//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "BxStackingEBMessenger.hh"
#include "BxStackingEB.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "BxReadParameters.hh"
#include "BxOutputVertex.hh"
#include "BxLogger.hh"
#include "G4Tokenizer.hh"

BxStackingEBMessenger::BxStackingEBMessenger(BxStackingEB* stack){
  BxStackingEB* stacking = stack;
  stack=stacking;
  /*  
  fDirectory = new G4UIdirectory("/bx/stack/eb/");
  fDirectory->SetGuidance("Control of BxGeneb event generator");
  
  fVesselSupDef = new G4UIcmdWithAString("/bx/stack/eb/VesselSupDef", this);
 
  fVesselThickness = new G4UIcmdWithADoubleAndUnit("/bx/stack/eb/VesselThickness", this);
  fVesselThickness->SetDefaultUnit("mm");
  fVesselThickness->SetUnitCandidates("mm cm m");
 
  fRhoCut = new G4UIcmdWithADoubleAndUnit("/bx/stack/eb/RhoCut", this);
  fRhoCut->SetDefaultUnit("mm");
  fRhoCut->SetUnitCandidates("mm cm m");
  
  fZCut = new G4UIcmdWithADoubleAndUnit("/bx/stack/eb/ZCut", this);
  fZCut->SetDefaultUnit("mm");
  fZCut->SetUnitCandidates("mm cm m");

  fBufferDef = new G4UIcmdWithAString("/bx/stack/eb/BufferDef", this);
  
  fRangeCmd = new G4UIcmdWithADoubleAndUnit("/bx/stack/eb/killRadius",this);
  fRangeCmd->SetDefaultUnit("cm");
  fRangeCmd->SetUnitCandidates("mm cm m");

  fKillAngleCmd = new G4UIcmdWithADoubleAndUnit("/bx/stack/eb/killAngle",this);
  fKillAngleCmd->SetGuidance("kill particles not directed toward the center with a given angle");
  fKillAngleCmd->SetDefaultUnit("deg");
  fKillAngleCmd->SetUnitCandidates("rad deg");
  
  fBulkCmd = new G4UIcmdWithABool("/bx/stack/eb/bulk",this);

  fKillParticleCmd = new G4UIcmdWithAnInteger("/bx/stack/eb/kill",this);

  fKillLEParticleCmd = new G4UIcommand("/bx/stack/eb/killLE",this);
  fKillLEParticleCmd->SetGuidance("Set properties of ion to be generated.");
  fKillLEParticleCmd->SetGuidance("[usage]/bx/stack/eb/killLE  pdg E");
  fKillLEParticleCmd->SetGuidance("	 pdg:(int) particle code");
  fKillLEParticleCmd->SetGuidance("	 E:(double) kinetic energy (in keV)");

  G4UIparameter* param;
  param = new G4UIparameter("pdg",'i',false);
  param->SetDefaultValue("-1000");
  fKillLEParticleCmd->SetParameter(param);
  param = new G4UIparameter("E",'d',true);
  param->SetDefaultValue("0.0");
  fKillLEParticleCmd->SetParameter(param);
 */
}


BxStackingEBMessenger::~BxStackingEBMessenger() {
/*
  delete fDirectory;
  delete fRangeCmd;
  delete fBulkCmd;
  delete fKillParticleCmd;
  delete fKillLEParticleCmd;
  delete fKillAngleCmd;
  delete fVesselSupDef;
  delete fBufferDef;
  delete fVesselThickness;
  delete fRhoCut;
  delete fZCut;
*/}


void BxStackingEBMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) { 
//dummy operation to avoid warning
	if (cmd) G4cout << "EB command " << newValue << endl;
	/* if (cmd == fRangeCmd){
         BxOutputVertex::Get()->SetPostponeFlag(true);
         stacking->SetIsReclassify(true);
         //stacking->SetRadius(fRangeCmd->GetNewDoubleValue(newValue));
	 BxOutputVertex::Get()->SetKillRadius(fRangeCmd->GetNewDoubleValue(newValue));
	 BxLog(routine) << "Only events with at least a deposit with r < " << newValue << " are selected"  << endlog ;
   }  else if (cmd == fBulkCmd){
         stacking->SetIsReclassify(fBulkCmd->ConvertToBool(newValue));
	 BxLog(routine) << "Only events with at least a deposit in the bulk are selected"  << endlog ;
   }  else if (cmd == fKillParticleCmd){
         stacking->KillParticles(fKillParticleCmd->ConvertToInt(newValue));
	 BxLog(routine) << "Kill particle with PDG code: " << newValue   << endlog ;
   }  else if( cmd == fKillLEParticleCmd) {
      G4Tokenizer next( newValue );
      G4int PDG   = StoI(next());
      G4double ene = StoD(next())*keV ;
      stacking->KillLEParticles(PDG,ene); 
      BxLog(routine) << "Kill particle with PDG code: " <<  PDG 
       << " and kinetic energy below " << ene/keV << " keV"   << endlog ;
    } else if(cmd == fVesselSupDef) {
         
	    BxReadParameters::Get()->SetVesselFileName(newValue);
	    BxReadParameters::Get()->ReadVesselGeometry();
	    stacking->SetfSupDef(true);
	    BxLog(routine) << "Deformation surface active?  " << stacking->GetfSupDef() <<  endlog; 
	    BxLog(routine) << "Deformed vessel filename: " << newValue  <<  endlog; 
    } else if(cmd == fVesselThickness) {
	    stacking->SetfVesThick(fVesselThickness->GetNewDoubleValue(newValue)); 
	    BxLog(routine) << "Deformed vessel thickness: " << fVesselThickness->GetNewDoubleValue(newValue)/mm  << " mm"   <<  endlog;
    } else if(cmd == fRhoCut) {
            stacking->SetffRhocut(true);
	    stacking->SetfRhoCut(fRhoCut->GetNewDoubleValue(newValue)); 
	    BxLog(routine) << "Setting rho < " << fRhoCut->GetNewDoubleValue(newValue)/mm  << " mm"   <<  endlog;
    } else if(cmd == fZCut) {
            stacking->SetffZcut(true);
	    stacking->SetfZCut(fZCut->GetNewDoubleValue(newValue)); 
	    BxLog(routine) << "Setting z < " << fZCut->GetNewDoubleValue(newValue)/mm  << " mm"   <<  endlog;
    } else if(cmd == fBufferDef) {
         
	    BxReadParameters::Get()->SetVesselFileName(newValue);
	    BxReadParameters::Get()->ReadVesselGeometry();
	    stacking->SetfBufDef(true);
	    BxLog(routine) << "Deformation buffer active?  " << stacking->GetfBufDef() <<  endlog; 
	    BxLog(routine) << "Deformed vessel filename: " << newValue  <<  endlog; 
    } if (cmd == fKillAngleCmd){
         stacking->SetKillingAngle(fKillAngleCmd->GetNewDoubleValue(newValue));
	 BxLog(routine) << "Kill particles not directed toward the center within  < " << fKillAngleCmd->GetNewDoubleValue(newValue)/deg
	 << " degree"  << endlog ;
    }
*/}
