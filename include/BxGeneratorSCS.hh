// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#ifndef _BxGeneratorSCS_HH
#define _BxGeneratorSCS_HH

//---------------------------------------------------------------------------//


#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "BxVGenerator.hh"
#include "G4Event.hh"
#include  <vector>
#include "G4SPSRandomGenerator.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"
#include <complex>

/**Special cross section generator (C14, Bi210)
*It can generate the shape according to the literature
*/

class BxGeneratorSCSMessenger;

//---------------------------------------------------------------------------//

class BxGeneratorSCS : public BxVGenerator {
public:

  ///default constructor
  BxGeneratorSCS();

  ///destructor
  virtual ~BxGeneratorSCS();

  ///public interface
  virtual void BxGeneratePrimaries(G4Event *event);
  
  void SetCrossSection(const G4String& a) {fCrossSection = a;}
  
  void SetIsNewPosition(G4bool a) { IsNewPosition = a; }

  void SetParticlePosition(G4ThreeVector pos)         {  thePosition = pos; }

  void SetParticleMomentumDirection(G4ThreeVector dir){  theDirection = dir; }
  
  void SetIfVolDist(G4bool a) { fVolumeFlag = a; }
 
    ///Allows user to choose Point, Plane, Surface or Volume source position distributions.
  void SetPosDisType(const G4String& string) { fSPSPos->SetPosDisType(string); }

    ///Allows the user to choose the particular shape they wish for the position distribution. Choices are Square, Circle, Ellipse, Rectangle, Sphere, Ellipsoid, Cylinder, Parallelepiped.
  void SetPosDisShape(const G4String& string) { fSPSPos->SetPosDisShape(string); }

    ///Sets the co-ordinates of the centre of the position distribution.
  void SetCentreCoords(const G4ThreeVector& pos) { fSPSPos->SetCentreCoords(pos); }

    ///Used to specify the co-ordinate system for the position distribution along with SetPosRot2. SetPosRot1 sets the vector x' and need not be a unit vector.
  void SetPosRot1(G4ThreeVector rot) { fSPSPos->SetPosRot1(rot); }
  
    ///Used in connection with SetPosRot1. This sets a vector in the plane x'y'. By a series of cross products x', y', z' are generated. Again need not be a unit vector.
  void SetPosRot2(G4ThreeVector rot){ fSPSPos->SetPosRot2(rot); }

     ///Sets the radius where appropriate for source distribution shapes.
  void SetRadius(G4double rad){ fSPSPos->SetRadius(rad); }

    ///Sets the inner radius where appropriate for source distribution shapes.
  void SetRadius0(G4double rad){ fSPSPos->SetRadius0(rad); }

    ///Used to confine the start positions to a particular volume.
  void ConfineSourceToVolume(const G4String& vol){ fSPSPos->ConfineSourceToVolume(vol); }
  
  void SetShapeFactorTerm(G4int val){ fShape = val; }

  void SetShapeFactor(G4double  val){ fC14ShapeFactor = val; }
 
  void SetNumberOfBins(G4int  val){ fSize = val; }
 ///Calculates the gamma function
  complex<G4double>  CGAMMA (complex<G4double>) ;

  
  //private  members
private:

    void C14CrossSection();
    void C14OlegCrossSection();
    void Bi210CrossSection();
    
    G4double ShootEnergy();
    
    
    G4String               fCrossSection;
    //G4ParticleGun*         fParticleGun  ;
    G4int                  fSize ;
    G4ThreeVector          thePosition  ;
    G4ThreeVector          theDirection ;
    //G4double*              theEnergy;
    G4double               fC14ShapeFactor ;
    G4ParticleDefinition*  particle_definition;
    G4ParticleTable*       theParticleTable;
    G4ThreeVector	   theNewPosition;
    //G4int		   fTipo;
    G4SPSPosDistribution*  fSPSPos;
    G4SPSAngDistribution*  fSPSAng;
    G4bool		   fVolumeFlag ;
    G4bool		   IsNewPosition;
    G4int                  fShape;
    G4bool                 fFirstTime;
    BxGeneratorSCSMessenger  *fMessenger;
    vector<G4double>         fShootRandom;
    vector<G4double>         fShootEnergy;
    
    G4double                 ffp[9];
    
  
};
#endif
