// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#ifndef _BXGENERATORMultiEvent_HH
#define _BXGENERATORMultiEvent_HH

//---------------------------------------------------------------------------//

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "BxVGenerator.hh"
#include "BxGeneratorMultiEventMessenger.hh"
#include "BxReadParameters.hh"
#include "G4Event.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"
#include "BxOutputStructure.hh"

//---------------------------------------------------------------------------//


class BxGeneratorMultiEvent : public BxVGenerator {
public:

  //default constructor
  BxGeneratorMultiEvent();

  //destructor
  virtual ~BxGeneratorMultiEvent();

  //public interface
  virtual void BxGeneratePrimaries(G4Event *event);

  void SetParticlePosition(G4ThreeVector pos) { fPosition  = pos; }

  
  void SetIfVolDist(G4bool a) { fVolumeFlag = a; }
 


  G4int GetNumberOfParticles() const { return fPDG.size(); }
  void SetfParticles( G4ThreeVector fPart);   


    ///Allows user to choose Point, Plane, Surface or Volume source position distributions.
  void SetPosDisType(const G4String& fstring) { fSPSPos->SetPosDisType(fstring); }
 
   ///Allows the user to choose the particular shape they wish for the position distribution. Choices are Square, Circle, Ellipse, Rectangle,
   ///Sphere, Ellipsoid, Cylinder, Parallelepiped.
  void SetPosDisShape(const G4String& fstring ) { fSPSPos->SetPosDisShape(fstring); }

    ///Sets the co-ordinates of the centre of the position distribution.
  void SetCentreCoords(G4ThreeVector pos) { fSPSPos->SetCentreCoords(pos); }

    ///Used to specify the co-ordinate system for the position distribution
    ///along with SetPosRot2. SetPosRot1 sets the vector x' and need not be
    ///a unit vector.
  void SetPosRot1(G4ThreeVector rot) { fSPSPos->SetPosRot1(rot); }
  
     ///Used in connection with SetPosRot1. This sets a vector in the plane
     ///x'y'. By a series of cross products x', y', z' are generated. Again
     ///need not be a unit vector.
   void SetPosRot2(G4ThreeVector rot){ fSPSPos->SetPosRot2(rot); }
   
      ///Sets the radius where appropriate for source distribution shapes.
   void SetRadius(G4double radius){ fSPSPos->SetRadius(radius); }

     ///Sets the inner radius where appropriate for source distribution shapes.
   void SetRadius0(G4double radius){ fSPSPos->SetRadius0(radius); }
  
     ///Used to confine the start positions to a particular volume.
   void ConfineSourceToVolume(const G4String& vol){ fSPSPos->ConfineSourceToVolume(vol); }
  
   void SetListOfParticles (MultiEvent val)   { fListOfParticles.push_back(val) ; }
  
  //private  members
private:

  BxGeneratorMultiEventMessenger*   fMessenger;

  G4ParticleTable*             fParticleTable;
  G4ParticleDefinition*        fParticle;
  G4ThreeVector                fPosition;
  G4ThreeVector                fDirection;
  

  vector<G4int>                fPDG;
  vector<G4double>             fPDGEnergy;
  vector<G4double>             fPDGBR;
  
  G4SPSPosDistribution*        fSPSPos;
  G4SPSAngDistribution*        fSPSAng;

  G4bool                       fVolumeFlag ;
  G4bool                       fRead ;
  
  vector<MultiEvent>           fListOfParticles;
  G4int                        fNVertex;
  vector<G4double>             fBR;

};
#endif
