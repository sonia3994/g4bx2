// --------------------------------------------------------------------------//
/** 
 * AUTHOR: I.Machulin
 * CONTACT: machulin@lngs.infn.it
* Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#ifndef _BXGENERATORLI9HE8_HH
#define _BXGENERATORLI9HE8_HH

//---------------------------------------------------------------------------//

#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "BxVGenerator.hh"
#include "BxGeneratorLi9He8Messenger.hh"
#include "BxReadParameters.hh"
#include "G4Event.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"

//---------------------------------------------------------------------------//
///Simulate the beta-neutron branches for Li9 and He8

/**Please use the stacking action command in macro
 * file /bx/stack/select to postpone the neutron capture
 * event. This is necessary for bx-elec simultion to pro-
 * cess the neutron capture as the new event.
*/
class BxGeneratorLi9He8 : public BxVGenerator {
public:

  ///default constructor
  BxGeneratorLi9He8();

  //copy constructor
  //BxGeneratorLi9He8(const BxGeneratorLi9He8 &);

  ///destructor
  virtual ~BxGeneratorLi9He8();

  ///public interface
  virtual void BxGeneratePrimaries(G4Event *event);



  void SetParticlePosition(G4ThreeVector pos){ fPosition  = pos; }
 // void SetParticleMomentumDirection(G4ParticleMomentum mom){ fParticleGun->SetParticleMomentumDirection(mom); }
  
  void SetIfVolDist(G4bool a) { fVolumeFlag = a; }  
 
    ///Allows user to choose Point, Plane, Surface or Volume source position distributions.
  void SetPosDisType(G4String fstring) { fSPSPos->SetPosDisType(fstring); }
 
   ///Allows the user to choose the particular shape they wish for the position distribution. Choices are Square, Circle, Ellipse, Rectangle, Sphere, Ellipsoid, Cylinder, Parallelepiped.
  void SetPosDisShape(G4String fstring ) { fSPSPos->SetPosDisShape(fstring); }

    ///Sets the co-ordinates of the centre of the position distribution.
  void SetCentreCoords(G4ThreeVector pos) { fSPSPos->SetCentreCoords(pos); }

    ///Used to specify the co-ordinate system for the position distribution along with SetPosRot2. SetPosRot1 sets the vector x' and need not be a unit vector.
  void SetPosRot1(G4ThreeVector rot) { fSPSPos->SetPosRot1(rot); }
  
     ///Used in connection with SetPosRot1. This sets a vector in the plane x'y'. By a series of cross products x', y', z' are generated. Again need not be a unit vector.
   void SetPosRot2(G4ThreeVector rot){ fSPSPos->SetPosRot2(rot); }
   
      ///Sets the radius where appropriate for source distribution shapes.
   void SetRadius(G4double rad){ fSPSPos->SetRadius(rad); }

     ///Sets the inner radius where appropriate for source distribution shapes.
   void SetRadius0(G4double rad){ fSPSPos->SetRadius0(rad); }
  
     ///Used to confine the start positions to a particular volume.
   void ConfineSourceToVolume(G4String vol){ fSPSPos->ConfineSourceToVolume(vol); }
  
///Used to choose the generator type (9Li=1 or 8He=0)
   void     SetNeutrinoType(G4int k)  { fNeutrinoType = k ;}
///Used to get the generator type (9Li=1 or 8He=0)
   G4int    GetNeutrinoType()         { return fNeutrinoType ;}

  //protected members

  enum Neutrinos
  {He8N=0,Li9N=1};

  //private  members
private:
   
  G4int ChooseBranchLi();
  G4int ChooseBranchHe();
  
  G4double ShootEnergyNeutronLi(G4int);
  G4double ShootEnergyNeutronHe(G4int);

  G4double ShootEnergyAlpha1Li(G4int);
  G4double ShootEnergyAlpha2Li(G4int,G4double);
  
  G4double ShootEnergyElectronHe(G4double);
  G4double ShootEnergyElectronLi(G4double);

  BxGeneratorLi9He8Messenger*   fTheMessenger;

  G4ParticleTable*             fParticleTable;
  G4ParticleDefinition*        fParticle;
  G4ThreeVector                fPosition;
  G4ThreeVector                fDirection;

  G4SPSPosDistribution*        fSPSPos;
  G4SPSAngDistribution*        fSPSAng;
  //G4SPSEneDistribution*        fSPSEne;
  G4bool                       fVolumeFlag ;
  //G4bool                       fEnergyDistribution;

  vector<G4double>             fEnergyBin;
  vector<G4double>             fProbability ;

  G4int                        fNeutrinoType;
  G4bool                       isFirstTime;
};
#endif
