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
#ifndef BxGeneratorRDMUIcmdWithNucleusAndUnit_h
#define BXGeneratorRDMUIcmdWithNucleusAndUnit_h 1
/** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
* MODULE:              MGGeneratorRDMUIcmdWithNucleusAndUnit.hh
*
* Version:             0.b.3
* Date:                29/02/00
* Author:              F Lei & P R Truscott
* Organisation:        DERA UK
* Customer:            ESA/ESTEC, NOORDWIJK
* Contract:            12115/96/JG/NL Work Order No. 3
*
* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
* Class Description
*
* The G4MGGeneratorRDMUIcmdWithNucleusAndUnit permits input of the nucleus definition
* in terms of its (atomic weight, atomic number, excitation energy).
* Input is expected in the form:
*
*			A, Z, E (energy unit)
*
* where A, Z, E are respectively the atomic weight, atomic number and
* excitation energy.  The energy unit can be any of the geant4 defined energy
* units.  The default is "keV"
*
* class description - end
*/
#include "G4UIcommand.hh"
#include "globals.hh"

#include "BxGeneratorRDMNucleus.hh"
////////////////////////////////////////////////////////////////////////////////
//
///It permits input of the nucleus definitionin terms of its atomic weight, atomic number and excitation energy.
class BxGeneratorRDMUIcmdWithNucleusAndUnit : public G4UIcommand
{
public:  // with description
  ///    Constructor identifying the command path in the User Interface and the associated G4UImessenger which will use this G4UIcommand object.
  BxGeneratorRDMUIcmdWithNucleusAndUnit
  (const char * theCommandPath,G4UImessenger * theMessenger);
  //
  ///    Desctructor
  ~BxGeneratorRDMUIcmdWithNucleusAndUnit();
  //
  ///    Extracts the values A, Z, E and unit from paramString.
  BxGeneratorRDMNucleus GetNewNucleusValue(G4String paramString);
  //
  ///    Returns the value of the unit (paramString) as defined in geant4
  G4double GetNewUnitValue(G4String paramString);
  //
  ///    Converts the Nucleus defined by nuc and the associated units of energy *unit into a G4String.
  G4String ConvertToString(BxGeneratorRDMNucleus nuc, const char * unit);
  

  ///    Identifies the parameter names associated with A, Z, E
void SetParameterName(const char * theNameA,const char * theNameZ,
                        const char * theNameE,
                        G4bool omittable, G4bool currentAsDefault=true);
  //
  ///    Sets the default Nucleus if the command is invoked without any parameters.
  void SetDefaultValue(BxGeneratorRDMNucleus defVal);
  ///    Sets the list of unit candidates
  void SetUnitCandidates(const char * candidateList);
  //
  //    Sets the default unit if the command is invoked without any parameters.
  void SetDefaultUnit(const char * defUnit);
};
////////////////////////////////////////////////////////////////////////////////
#endif

