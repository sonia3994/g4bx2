//---------------------------------------------------------------------------//
/**                                                            
 *      
 * CLASS DECLARATION:  BxVGenerator.hh
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

#ifndef _BXVGENERATOR_HH
#define _BXVGENERATOR_HH

//---------------------------------------------------------------------------//

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleGun.hh"
#include "G4RadioactiveDecay.hh"


//---------------------------------------------------------------------------//

class G4Event;
class G4UImessenger;
class G4Run;

//---------------------------------------------------------------------------//

class BxVGenerator {
  public:

    ///default constructor
     BxVGenerator(const G4String &myname);


    ///destructor
    virtual ~BxVGenerator();

    //public interface
    /// Called at beginning of Run

    G4String GetGeneratorName() { return fGeneratorName; }
    void SetReportingFrequency(G4int freq) { fReportingFrequency = freq; }

    /// Public interface to BxGeneratePrimaries
    virtual void BxGeneratePrimaries(G4Event *event) = 0;

  private:
    G4String        fGeneratorName;  // Name of Generator.

    //protected members
  protected:
    G4UImessenger*   fG4Messenger;   // G4Messenger for setting up generator.
    G4int            fReportingFrequency; // Generate report every fReportingFrequency events;  
};
#endif
