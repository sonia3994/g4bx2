//
// --------------------------------------------------------------------------//
/**
* AUTHOR: 
* CONTACT: ydc@lbl.gov
* FIRST SUBMISSION: Wed Mar 10 14:42:07 PST 2004
*
* REVISION:Revised by A. Caminata and S. Marcocci, Sept. 2014
*
* 06-02-2004, Initial port to MG
* mm-dd-yyyy, What is changed, Whoami
*/
// --------------------------------------------------------------------------//
#include "G4UIcommand.hh"

#include "BxGeneratorRDMUIcmdWithNucleusAndUnit.hh"

#include <sstream>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//
BxGeneratorRDMUIcmdWithNucleusAndUnit::BxGeneratorRDMUIcmdWithNucleusAndUnit
(const char * theCommandPath,G4UImessenger * theMessenger)
:G4UIcommand(theCommandPath,theMessenger)
{
  G4UIparameter * intParamA = new G4UIparameter('i');
  SetParameter(intParamA);
  G4UIparameter * intParamZ = new G4UIparameter('i');
  SetParameter(intParamZ);
  G4UIparameter * dblParamE = new G4UIparameter('d');
  SetParameter(dblParamE);
  G4UIparameter * untParam = new G4UIparameter('s');
  SetParameter(untParam);
  untParam->SetParameterName("Unit");
}
////////////////////////////////////////////////////////////////////////////////
//
BxGeneratorRDMUIcmdWithNucleusAndUnit::~BxGeneratorRDMUIcmdWithNucleusAndUnit()
{
  ;
}
////////////////////////////////////////////////////////////////////////////////
//
BxGeneratorRDMNucleus BxGeneratorRDMUIcmdWithNucleusAndUnit::GetNewNucleusValue(G4String paramString)
{
  G4int a;
  G4int z;
  G4double e;
  char unts[30];

  istringstream is(paramString);
  is >> a >> z >> e >>unts;
  G4String unt = unts;

  return BxGeneratorRDMNucleus(a,z,e*ValueOf(unt));
}

G4double BxGeneratorRDMUIcmdWithNucleusAndUnit::GetNewUnitValue(G4String paramString)
{
  G4int a;
  G4int z;
  G4double e;  

  char unts[30];
  
  istringstream is(paramString);
  is >> a >> z >> e  >> unts;

  G4String unt = unts;
  
  return ValueOf(unt);
}

////////////////////////////////////////////////////////////////////////////////
//
G4String BxGeneratorRDMUIcmdWithNucleusAndUnit::ConvertToString(BxGeneratorRDMNucleus def, 
						    const char *unitName)
{
  G4double uv = ValueOf(unitName);

  ostringstream os;
  os << def.GetA() << " " << def.GetZ()
     << " "<< def.GetE()/uv<<" "<< unitName <<  '\0';
  return G4String(os.str());
}                         
////////////////////////////////////////////////////////////////////////////////
//
void BxGeneratorRDMUIcmdWithNucleusAndUnit::SetParameterName
(const char * theNameA, const char * theNameZ,
const char * theNameE,G4bool omittable,G4bool currentAsDefault)
{
  G4UIparameter * theParamA = GetParameter(0);
  theParamA->SetParameterName(theNameA);
  theParamA->SetOmittable(omittable);
  theParamA->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter * theParamZ = GetParameter(1);
  theParamZ->SetParameterName(theNameZ);
  theParamZ->SetOmittable(omittable);
  theParamZ->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter * theParamE = GetParameter(2);
  theParamE->SetParameterName(theNameE);
  theParamE->SetOmittable(omittable);
  theParamE->SetCurrentAsDefault(currentAsDefault);
}
////////////////////////////////////////////////////////////////////////////////
//
void BxGeneratorRDMUIcmdWithNucleusAndUnit::SetDefaultValue(BxGeneratorRDMNucleus def)
{
  G4UIparameter * theParamA = GetParameter(0);
  theParamA->SetDefaultValue(def.GetA());
  G4UIparameter * theParamZ = GetParameter(1);
  theParamZ->SetDefaultValue(def.GetZ());
  G4UIparameter * theParamE = GetParameter(2);
  theParamE->SetDefaultValue(def.GetE());
}


void BxGeneratorRDMUIcmdWithNucleusAndUnit::SetUnitCandidates(const char * candidateList)
{
  G4UIparameter * untParam = GetParameter(3);
  G4String canList = candidateList;
  untParam->SetParameterCandidates(canList);
}

void BxGeneratorRDMUIcmdWithNucleusAndUnit::SetDefaultUnit(const char * defUnit)
{
  G4UIparameter * untParam = GetParameter(3);
  untParam->SetOmittable(true);
  untParam->SetDefaultValue(defUnit);
}












