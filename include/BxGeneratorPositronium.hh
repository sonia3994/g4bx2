// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#ifndef _BXGENERATORPositronium_HH
#define _BXGENERATORPositronium_HH

//---------------------------------------------------------------------------//


#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "BxVGenerator.hh"
#include "BxGeneratorPositroniumMessenger.hh"
#include "BxReadParameters.hh"
#include "G4Event.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"
#include "G4SystemOfUnits.hh"
//---------------------------------------------------------------------------//

///C11 decay simulation
/**This generator allows to simulate 11 C events taking into
 * account the possible formation of orthoPositronium and
 * its decay dalayed in time. Probability of oPs formation
 * and mean live are settable. It is possible to generate oPs
 * events for monocromatic positrons from 11 C or following
 * their β + energy spectrum.
 */
class BxGeneratorPositronium : public BxVGenerator {
public:

  ///default constructor
  BxGeneratorPositronium();

  ///destructor
  virtual ~BxGeneratorPositronium();

  ///public interface
  virtual void BxGeneratePrimaries(G4Event *event);

///Set the particle energy
  void SetParticleEnergy(G4double a){ fEnergy = a; }
///Get the particle energy
  G4double GetParticleEnergy()      { return fEnergy/MeV; }
///Set the orthopositronium probability 
  void SetPsProbability(G4double a) { fOrthoProb = a; }
///Set the orthopositronium mean life 
  void SetPsMeanLife(G4double a)    { fTau = a; }
///Set the particle position
  void SetParticlePosition(G4ThreeVector pos){ fPosition = pos; }
///Set the particle momentum 
  void SetParticleMomentumDirection(G4ParticleMomentum mom){ fDirection = mom; }
  
  void SetIfVolDist(G4bool a) { fVolumeFlag = a; }  
 
    ///Allows user to choose Point, Plane, Surface or Volume source position distributions.
  void SetPosDisType(const G4String& fstring) { fSPSPos->SetPosDisType(fstring); }
 
   //Allows the user to choose the particular shape they wish for the position distribution. Choices are Square, Circle, Ellipse, Rectangle, Sphere, Ellipsoid, Cylinder, Parallelepiped.
  void SetPosDisShape(const G4String& fstring ) { fSPSPos->SetPosDisShape(fstring); }

    ///Sets the co-ordinates of the centre of the position distribution.
  void SetCentreCoords(G4ThreeVector pos) { fSPSPos->SetCentreCoords(pos); }

    ///Used to specify the co-ordinate system for the position distribution along with SetPosRot2. SetPosRot1 sets the vector x' and need not be a unit vector.
  void SetPosRot1(G4ThreeVector rot) { fSPSPos->SetPosRot1(rot); }
  
     ///Used in connection with SetPosRot1. This sets a vector in the plane x'y'. By a series of cross products x', y', z' are generated. Again need not be a unit vector.
   void SetPosRot2(G4ThreeVector rot){ fSPSPos->SetPosRot2(rot); }
   
      ///Sets the radius where appropriate for source distribution shapes.
   void SetRadius(G4double radius){ fSPSPos->SetRadius(radius); }

     ///Sets the inner radius where appropriate for source distribution shapes.
   void SetRadius0(G4double radius){ fSPSPos->SetRadius0(radius); }
  
     ///Used to confine the start positions to a particular volume.
   void ConfineSourceToVolume(const G4String& vol){ fSPSPos->ConfineSourceToVolume(vol); }
  
///Set the generation of the energy following the spectrum of 11 C
   void SetSpectrum(G4bool val) { fSpectrum = val; }  
///Get the C11 spectrum
   G4bool GetSpectrum()         { return fSpectrum; }      

  //private  members
private:
    
  BxGeneratorPositroniumMessenger*   fTheMessenger;
  G4double                     ShootEnergy();
  G4SPSPosDistribution*        fSPSPos;
  G4SPSAngDistribution*        fSPSAng;
  //G4SPSEneDistribution*        fSPSEne;
  G4bool                       fVolumeFlag ;
  //G4bool                       fEnergyDistribution;
  
  G4ParticleTable*             fParticleTable;
  G4double                     fTau;
  G4double                     fOrthoProb;
  G4double                     fEnergy;
  G4ThreeVector                fPosition ;
  G4ThreeVector                fDirection ;
  G4bool                       isFirstTime;
  G4ParticleDefinition         *fElectron ;
  G4ParticleDefinition         *fGamma ;
  vector<G4double>             fShootRandom;
  vector<G4double>             fShootEnergy;
  G4bool                       isRead;
  G4bool                       fSpectrum;
  
};
#endif
