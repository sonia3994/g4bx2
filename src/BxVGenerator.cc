//---------------------------------------------------------------------------//
/**                                                            
 *      
 * CLASS DECLARATION:  BxVGenerator.cc
 *
 *---------------------------------------------------------------------------//
 *
 * DESCRIPTION: 
 *
 */ 
// Begin description of class here
/**
 *
 * Pure virtual base class for Bx generators. 
 * 
 */
// End class description
//
/**  
 * SPECIAL NOTES:
 *
 */
// 
// --------------------------------------------------------------------------//
/** 
 * AUTHOR: davide.franco@mi.infn.it
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
 */
// --------------------------------------------------------------------------//

#include "BxVGenerator.hh"

//---------------------------------------------------------------------------//

BxVGenerator::BxVGenerator(const G4String &myname):
  fGeneratorName(myname), fG4Messenger(0), fReportingFrequency(1000)
{;}

//---------------------------------------------------------------------------//

BxVGenerator::~BxVGenerator()
{;}

//---------------------------------------------------------------------------//
