//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "BxStackingRDMMessenger.hh"
#include "BxStackingRDM.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "BxReadParameters.hh"
#include "BxOutputVertex.hh"
#include "BxLogger.hh"
#include "G4Tokenizer.hh"

BxStackingRDMMessenger::BxStackingRDMMessenger(BxStackingRDM* stack){
  stacking = stack;
  
  fDirectory = new G4UIdirectory("/bx/stack/rdm/");
  fDirectory->SetGuidance("Control of BxGeneb event generator");
  

  fKillAngleCmd = new G4UIcmdWithADoubleAndUnit("/bx/stack/rdm/killAngle",this);
  fKillAngleCmd->SetGuidance("kill particles not directed toward the center with a given angle");
  fKillAngleCmd->SetDefaultUnit("deg");
  fKillAngleCmd->SetUnitCandidates("rad deg");
  
  fBulkCmd = new G4UIcmdWithABool("/bx/stack/rdm/bulk",this);

  fKillParticleCmd = new G4UIcmdWithAnInteger("/bx/stack/rdm/kill",this);

  fKillLEParticleCmd = new G4UIcommand("/bx/stack/rdm/killLE",this);
  fKillLEParticleCmd->SetGuidance("Set properties of ion to be generated.");
  fKillLEParticleCmd->SetGuidance("[usage]/bx/stack/rdm/killLE  pdg E");
  fKillLEParticleCmd->SetGuidance("	 pdg:(int) particle code");
  fKillLEParticleCmd->SetGuidance("	 E:(double) kinetic energy (in keV)");

  G4UIparameter* param;
  param = new G4UIparameter("pdg",'i',false);
  param->SetDefaultValue("-1000");
  fKillLEParticleCmd->SetParameter(param);
delete param; 
  param = new G4UIparameter("E",'d',true);
  param->SetDefaultValue("0.0");
  fKillLEParticleCmd->SetParameter(param);
delete param; 
}


BxStackingRDMMessenger::~BxStackingRDMMessenger() {

  delete fDirectory;
  delete fRangeCmd;
  delete fBulkCmd;
  delete fKillParticleCmd;
  delete fKillLEParticleCmd;
  delete fKillAngleCmd;
 
}


void BxStackingRDMMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) { 
   if (cmd == fRangeCmd){
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
    } if (cmd == fKillAngleCmd){
         stacking->SetKillingAngle(fKillAngleCmd->GetNewDoubleValue(newValue));
	 BxLog(routine) << "Kill particles not directed toward the center within  < " << fKillAngleCmd->GetNewDoubleValue(newValue)/deg
	 << " degree"  << endlog ;
    }
}
