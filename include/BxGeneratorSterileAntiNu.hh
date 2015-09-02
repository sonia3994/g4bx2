// --------------------------------------------------------------------------//
/** 
 * AUTHOR: M. Meyer
 * CONTACT: mikko.meyer@desy.de
*/
// --------------------------------------------------------------------------//

#ifndef _BXGENERATORSterileAntiNu_HH
#define _BXGENERATORSterileAntiNu_HH 1

//---------------------------------------------------------------------------//


#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "BxVGenerator.hh"
#include "BxGeneratorSterileAntiNuMessenger.hh"
#include "BxReadParameters.hh"
#include "G4Event.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"


#include "TFile.h"
#include "TH1F.h"

//---------------------------------------------------------------------------//

class BxGeneratorSterileAntiNu : public BxVGenerator {
public:

  //default constructor
  BxGeneratorSterileAntiNu();

  //copy constructor
  //BxGeneratorSterileAntiNu(const BxGeneratorSterileAntiNu &);

  //destructor
  virtual ~BxGeneratorSterileAntiNu();

  //public interface
  virtual void BxGeneratePrimaries(G4Event *event);



  void SetParticlePosition(G4ThreeVector pos){ fParticleGun->SetParticlePosition(pos); }
  void SetParticleDefinition(G4ParticleDefinition* part) { fParticleGun->SetParticleDefinition(part); }
  void SetParticleMomentumDirection(G4ParticleMomentum mom){ fParticleGun->SetParticleMomentumDirection(mom); }
  void SetParticleEnergy(G4double ene){ fParticleGun->SetParticleEnergy(ene); }
  
  void SetIfVolDist(G4bool a) { fVolumeFlag = a; }
 

  
  // Energy Distribution
  //    Allows the user to choose the energy distribution type. The arguments
  //    are Mono (mono-energetic), Lin (linear), Pow (power-law), Exp 
  //    (exponential), Gauss (gaussian), Brem (bremsstrahlung), BBody (black-body), Cdg
  //    (cosmic diffuse gamma-ray), User (user-defined), Arb (arbitrary
  //    point-wise), Epn (energy per nucleon).
  void SetEnergyDisType(G4String fstring) { fSPSEne->SetEnergyDisType(fstring); }  
  void SetEmin(G4double val)  { fSPSEne->SetEmin(val); }
  void SetEmax(G4double val)  { fSPSEne->SetEmax(val); }
   //Sets alpha for a power-law distribution.
  void SetAlpha(G4double val) { fSPSEne->SetAlpha(val); }
   //Sets Temperature for a Brem or BBody distributions.
  void SetTemp(G4double val)  { fSPSEne->SetTemp(val); }
   //Sets Ezero for an exponential distribution.
  void SetEzero(G4double val) { fSPSEne->SetEzero(val); }
  // Sets gradient for a linear distribution.
  void SetGradient(G4double val){ fSPSEne->SetGradient(val); }
  // Sets intercept for a linear distribution.	
  void SetInterCept(G4double val){ fSPSEne->SetInterCept(val); }
  
 
    //Allows user to choose Point, Plane, Surface or Volume source
    //position distributions.
  void SetPosDisType(G4String fstring) { fSPSPos->SetPosDisType(fstring); }
 
   //Allows the user to choose the particular shape they wish for the
    //osition distribution. Choices are Square, Circle, Ellipse, Rectangle,
    //Sphere, Ellipsoid, Cylinder, Parallelepiped.
  void SetPosDisShape(G4String fstring ) { fSPSPos->SetPosDisShape(fstring); }

    //Sets the co-ordinates of the centre of the position distribution.
  void SetCentreCoords(G4ThreeVector pos) { fSPSPos->SetCentreCoords(pos); }

    //Used to specify the co-ordinate system for the position distribution
    //along with SetPosRot2. SetPosRot1 sets the vector x' and need not be
    //a unit vector.
  void SetPosRot1(G4ThreeVector rot) { fSPSPos->SetPosRot1(rot); }
  
     //Used in connection with SetPosRot1. This sets a vector in the plane
     //x'y'. By a series of cross products x', y', z' are generated. Again
     //need not be a unit vector.
   void SetPosRot2(G4ThreeVector rot){ fSPSPos->SetPosRot2(rot); }
   
      //Sets the radius where appropriate for source distribution shapes.
   void SetRadius(G4double radius){ fSPSPos->SetRadius(radius); }

     //Sets the inner radius where appropriate for source distribution shapes.
   void SetRadius0(G4double radius){ fSPSPos->SetRadius0(radius); }
  
     //Used to confine the start positions to a particular volume.
   void ConfineSourceToVolume(G4String vol){ fSPSPos->ConfineSourceToVolume(vol); }
  
   void SetEnergyDistribution(G4bool val){ fEnergyDistribution = val ; }

  /// Neu:
  void     SetDM2(G4double val)      { fDeltaM2 =  val ;}
  G4double GetDM2()                  { return fDeltaM2 ;}

  void     SetSIN2T2(G4double val)   { fSin2theta2 = val;}
  G4double GetSIN2T2()               { return  fSin2theta2;}

  /// spectrum deformation:
  void     SetCOEFFICIENTA(G4double val)      { fshapeCoefficientA =  val ;}
  G4double GetCOEFFICIENTA()                  { return fshapeCoefficientA ;}

  void     SetCOEFFICIENTB(G4double val)      { fshapeCoefficientB =  val ;}
  G4double GetCOEFFICIENTB()                  { return fshapeCoefficientB ;}

  void     SetCOEFFICIENTC(G4double val)      { fshapeCoefficientC =  val ;}
  G4double GetCOEFFICIENTC()                  { return fshapeCoefficientC ;}
  // ENDE-Neu

  
  //protected members


  //private  members
private:

  G4double survivalProb(G4double ,  G4double) ;
  G4double shootNuEnergy(G4int);					// creates Neutrino Energy
  G4double shootEnergy(G4double, G4int);				// creates Positron/Electron Energy  
  G4double cross_nue(G4double, G4double) ;				// cross section for Neutrino Interaction
  G4ThreeVector GetInteractionPositionIntSource(G4ThreeVector) ;
  G4ThreeVector shootE_Direction(G4ThreeVector, G4ThreeVector, G4double, G4double);				// shoots Direction of electron (including kinematics)
  G4bool initShapeFactor(G4bool);
  G4bool normalizeHisto(TH1F*);


  BxGeneratorSterileAntiNuMessenger* fTheMessenger;
  G4ParticleGun*               fParticleGun;
  G4SPSPosDistribution*        fSPSPos;
  G4SPSAngDistribution*        fSPSAng;
  G4SPSEneDistribution*        fSPSEne;
  G4bool                       fVolumeFlag ;
  G4bool                       fEnergyDistribution;

  /// new:
  G4double                     fDeltaM2;
  G4double                     fSin2theta2;
  G4double                     fshapeCoefficientA;
  G4double                     fshapeCoefficientB;
  G4double                     fshapeCoefficientC;

  G4bool                       shapeFactorInit;

  G4ThreeVector                fNuPosition;
  G4ThreeVector                fPosition;

  // Globale Variablen:
/*
  G4int 			SimulationType;	
  G4double 			NuEndPoint;	
  G4int				numberOfSteps;
*/
 G4ParticleTable*            	fParticleTable;
 G4DynamicParticle* 		dynamicParticleLastEvent;
 G4ParticleDefinition*        	fParticle;

 // Neutron:
 G4double ShootEnergyNeutron();
 // G4ThreeVector                fPosition;
 G4ThreeVector                fDirection;


  /// For Rootifile IO and loading histograms:
  TFile *fSpectrumRootFile;
  TFile *fCrossSectionFile;
  TH1F  *hist_AntiNuSpec;
  TH1F  *hist_CrossSection_firstOrder;
  TH1F  *hist_Convolution_CrossSection_AntiNuSpec;
  TH1F  *hist_truePositronEnergy;
  G4double AntiNuEnergyMax;
  G4double AntiNuThreshold;
  // for antiNu spectral deformation (added Feb. 2015):
  TH1F  *hist_AntiNuSpec_ShapeFactor;


};
#endif
/*
 * $Log: BxGeneratorSterileAntiNu.hh,v $
 * Revision 1.2  2015/02/25 18:17:26  acaminata
 * Shape factor added
 *
 * Revision 1.1  2015/02/12 14:32:07  acaminata
 * Sterile generator added
 *
 * Revision 1.8  2007-11-12 12:09:44  dfranco
 * added to g4gun, the following  energy distributions:
 * Lin (linear), Pow (power-law), Exp (exponential), Gauss (gaussian),
 * Brem (bremsstrahlung), BBody (black-body), Cdg (cosmic diffuse gamma-ray)
 *
 * Revision 1.7  2007-03-22 14:48:49  dfranco
 * Stable version
 *
 */
