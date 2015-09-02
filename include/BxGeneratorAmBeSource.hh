// --------------------------------------------------------------------------//
/** 
 * AUTHOR: I.Machulin
 * CONTACT: machulin@lngs.infn.it
* Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#ifndef _BXGENERATORAMBESOURCE_HH
#define _BXGENERATORAMBESOURCE_HH

//---------------------------------------------------------------------------//


#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "BxVGenerator.hh"
#include "BxGeneratorAmBeSourceMessenger.hh"
#include "BxReadParameters.hh"
#include "G4Event.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"

//---------------------------------------------------------------------------//
/**Please use the stacking action command in macro
file /bx/stack/select to postpone the neutron capture
event. 
*This is necessary for bx-elec simultion to
process the neutron capture as the new event.
*/
///define the branch in AmBe calibrationi source products simulation
class BxGeneratorAmBeSource : public BxVGenerator {
public:

  ///default constructor
  BxGeneratorAmBeSource();

  ///destructor
  virtual ~BxGeneratorAmBeSource();

  ///public interface
  virtual void BxGeneratePrimaries(G4Event *event);



  void SetParticlePosition(G4ThreeVector pos){ fPosition  = pos; }
  
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
  

   void     SetNeutrinoType(G4int k)  { fNeutrinoType = k ;}
   G4int    GetNeutrinoType()         { return fNeutrinoType ;}

  
  ///protected members

  enum Neutrinos
  {AllPN=0,NeutronNoG=1,Neutron1G=2,Neutron2G=3};

  ///private  members
private:
   
  G4double ShootEnergyNeutron0G()  ;
  G4double ShootEnergyNeutron1G()  ;
  G4double ShootEnergyNeutron2G()  ;
  
  G4double NeutronEnergySp (G4int );

 
  BxGeneratorAmBeSourceMessenger*   fTheMessenger;

  G4ParticleTable*             fParticleTable;
  G4ParticleDefinition*        fParticle;
  G4ThreeVector                fPosition;
  G4ThreeVector                fDirection;

  G4SPSPosDistribution*        fSPSPos;
  G4SPSAngDistribution*        fSPSAng;
  G4SPSEneDistribution*        fSPSEne;
  G4bool                       fVolumeFlag ;
  G4bool                       fEnergyDistribution;

  vector<G4double>             fEnergyBin;
  vector<G4double>             fProbability ;

  G4int                        fNeutrinoType;
  G4bool                       isFirstTime;
};
#endif
