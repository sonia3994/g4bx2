#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "BxGeneratorSupernovaAntiNuMessenger.hh"
#include "G4UIcmdWithADouble.hh"
#include "BxGeneratorSupernovaAntiNu.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "BxReadParameters.hh"
#include "BxOutputVertex.hh"
#include "BxLogger.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
using namespace std;

BxGeneratorSupernovaAntiNuMessenger::BxGeneratorSupernovaAntiNuMessenger (BxGeneratorSupernovaAntiNu* gen){
	generator = gen;
	fDirectory = new G4UIdirectory("/bx/generator/SupernovaAntiNu/");
        fDirectory->SetGuidance("Control of BxSupernovaAntiNu event generator");

        fPositionCmd = new G4UIcmdWith3VectorAndUnit("/bx/generator/SupernovaAntiNu/position",this);
        fPositionCmd->SetGuidance("Set the gun position");
        fPositionCmd->SetUnitCategory("Length");
        fPositionCmd->SetDefaultUnit("cm");
        fPositionCmd->SetUnitCandidates("mm cm m");

        fGenInScintCmd = new G4UIcmdWithAnInteger("/bx/generator/SupernovaAntiNu/inScintillator",this);

        fNeutrinoCmd = new G4UIcmdWithAString("/bx/generator/SupernovaAntiNu/source",this);
        G4String candidates = "both neutron positron";
        fNeutrinoCmd->SetCandidates(candidates);
}

BxGeneratorSupernovaAntiNuMessenger::~BxGeneratorSupernovaAntiNuMessenger() {
    delete fDirectory;
    delete fPositionCmd;
    delete fNeutrinoCmd;
    delete fGenInScintCmd;
}

void BxGeneratorSupernovaAntiNuMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) {

    if (cmd == fPositionCmd) {
          generator->SetParticlePosition(fPositionCmd->ConvertToDimensioned3Vector(newValue));
          BxOutputVertex::Get()->SetSpatialDist(0);  // chto eto?????
    } else if (cmd == fNeutrinoCmd) {
        BxLog(routine) << "SupernovaAntiNU source: " << newValue << endlog ;
        if(newValue == "both") {
                generator->SetNeutrinoType(BxGeneratorSupernovaAntiNu::both);
        } else if(newValue == "neutron") {
                generator->SetNeutrinoType(BxGeneratorSupernovaAntiNu::neutron);
        } else if(newValue == "positron") {
                generator->SetNeutrinoType(BxGeneratorSupernovaAntiNu::positron);
        }

    } else if (cmd == fGenInScintCmd){
        generator->SetScintFlag(true);
    }
}
