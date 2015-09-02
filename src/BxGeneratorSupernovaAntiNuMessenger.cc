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

BxGeneratorSupernovaAntiNuMessenger::BxGeneratorSupernovaAntiNuMessenge (BxGeneratorSupernovaAntiNu* gen){
	generator = gen;
	fDirectory = new G4UIdirectory("/bx/generator/SupernovaAntiNu/");
	fDirectory->SetGuidance("Control of BxAntiNeutrino event generator");
