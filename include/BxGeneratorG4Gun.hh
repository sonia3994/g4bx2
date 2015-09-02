// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
* Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#ifndef _BXGENERATORG4GUN_HH
#define _BXGENERATORG4GUN_HH

//---------------------------------------------------------------------------//
/**the most general and flexible generator;
 * it allows to shoot single particles (e-, gamma, alpha,
 * etc.) in pointlike positions or with spatial distribu-
 * tions; moreover, you can generate energy distribu-
 * tions;*/
///the most general and flexible generator. It allows to shoot single particles
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "BxVGenerator.hh"
#include "BxGeneratorG4GunMessenger.hh"
#include "BxReadParameters.hh"
#include "G4Event.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"

//---------------------------------------------------------------------------//

class BxGeneratorG4Gun : public BxVGenerator {
public:

  ///default constructor
  BxGeneratorG4Gun();

  ///destructor
  virtual ~BxGeneratorG4Gun();

  ///public interface
  virtual void BxGeneratePrimaries(G4Event *event);

  void SetParticlePosition(G4ThreeVector pos){fParticleGun->SetParticlePosition(pos); }
  void SetParticleDefinition(G4ParticleDefinition* part) { fParticleGun->SetParticleDefinition(part); }
  void SetParticleMomentumDirection(G4ParticleMomentum mom){ fParticleGun->SetParticleMomentumDirection(mom); }
  void SetParticleEnergy(G4double ene){ fParticleGun->SetParticleEnergy(ene); }
  
  void SetIfVolDist(G4bool a) { fVolumeFlag = a; }
 
///Sets the energy distribution 
  /** Energy Distribution
  *    Allows the user to choose the energy distribution type. The arguments
  *    are Mono (mono-energetic), Lin (linear), Pow (power-law), Exp 
  *    (exponential), Gauss (gaussian), Brem (bremsstrahlung), BBody (black-body),
  *    Cdg (cosmic diffuse gamma-ray), User (user-defined), Arb (arbitrary point-wise), Epn (energy per nucleon).
  */
  void SetEnergyDisType(G4String fstring) { fSPSEne->SetEnergyDisType(fstring); }  
  void SetEmin(G4double val)  { fSPSEne->SetEmin(val); }
  void SetEmax(G4double val)  { fSPSEne->SetEmax(val); }
   ///Sets alpha for a power-law distribution.
  void SetAlpha(G4double val) { fSPSEne->SetAlpha(val); }
   ///Sets Temperature for a Brem or BBody distributions.
  void SetTemp(G4double val)  { fSPSEne->SetTemp(val); }
   ///Sets Ezero for an exponential distribution.
  void SetEzero(G4double val) { fSPSEne->SetEzero(val); }
  /// Sets gradient for a linear distribution.
  void SetGradient(G4double val){ fSPSEne->SetGradient(val); }
  /// Sets intercept for a linear distribution.	
  void SetInterCept(G4double val){ fSPSEne->SetInterCept(val); }
  
 
    ///Allows user to choose Point, Plane, Surface or Volume source position distributions.
  void SetPosDisType(G4String fstring) { fSPSPos->SetPosDisType(fstring); }
 
   //Allows the user to choose the particular shape they wish for the position distribution. Choices are Square, Circle, Ellipse, Rectangle, Sphere, Ellipsoid, Cylinder, Parallelepiped.
  void SetPosDisShape(G4String fstring ) { fSPSPos->SetPosDisShape(fstring); }

    ///Sets the co-ordinates of the centre of the position distribution.
  void SetCentreCoords(G4ThreeVector pos) { fSPSPos->SetCentreCoords(pos); }


    //Used to specify the co-ordinate system for the position distribution along with SetPosRot2. SetPosRot1 sets the vector x' and need not be a unit vector.
  void SetPosRot1(G4ThreeVector rot) { fSPSPos->SetPosRot1(rot); }
  
     //Used in connection with SetPosRot1. This sets a vector in the plane
     //x'y'. By a series of cross products x', y', z' are generated. Again
     //need not be a unit vector.
   void SetPosRot2(G4ThreeVector rot){ fSPSPos->SetPosRot2(rot); }
  
      ///Sets the radius where appropriate for source distribution shapes.
   void SetRadius(G4double r){ fSPSPos->SetRadius(r); }

     ///Sets the inner radius where appropriate for source distribution shapes.
   void SetRadius0(G4double r){ fSPSPos->SetRadius0(r);}
  
     ///Used to confine the start positions to a particular volume.
   void ConfineSourceToVolume(G4String vol){ fSPSPos->ConfineSourceToVolume(vol); }
  
///To set the energy distribution of the generated particle using the mac file.
    void SetEnergyDistribution(G4bool val){ fEnergyDistribution = val ; }
  
  //private  members
private:
  BxGeneratorG4GunMessenger*   fTheMessenger;
  G4ParticleGun*               fParticleGun;
  G4SPSPosDistribution*        fSPSPos;
  G4SPSAngDistribution*        fSPSAng;
  G4SPSEneDistribution*        fSPSEne;
  G4bool                       fVolumeFlag ;
  G4bool                       fEnergyDistribution;
  G4bool                       fDebugGenerator;
};
#endif
