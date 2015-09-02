// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#ifndef _BXGENERATORSOLARNEUTRINO_HH
#define _BXGENERATORSOLARNEUTRINO_HH

//---------------------------------------------------------------------------//


#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "BxVGenerator.hh"
#include "BxGeneratorSolarNeutrinoMessenger.hh"
#include "BxReadParameters.hh"
#include "G4Event.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"

//---------------------------------------------------------------------------//
///Solar Neutrino generator: generates solar neutrino spectra, including distortions due to oscillations; stores information only about recoiled electron
class BxGeneratorSolarNeutrino : public BxVGenerator {
public:

  ///default constructor
  BxGeneratorSolarNeutrino();

  ///destructor
  virtual ~BxGeneratorSolarNeutrino();

  ///public interface
  virtual void BxGeneratePrimaries(G4Event *event);
///Set the particle (e-) position
  void SetParticlePosition(G4ThreeVector pos){ fParticleGun->SetParticlePosition(pos); }
///Set the particle momentum
  void SetParticleMomentumDirection(G4ParticleMomentum mom){ fParticleGun->SetParticleMomentumDirection(mom); }
  
  void SetIfVolDist(G4bool a) { fVolumeFlag = a; }  
 
    ///Allows user to choose Point, Plane, Surface or Volume source position distributions.
  void SetPosDisType(const G4String& fstring) { fSPSPos->SetPosDisType(fstring); }
 
   ///Allows the user to choose the particular shape they wish for the position distribution. Choices are Square, Circle, Ellipse, Rectangle, Sphere, Ellipsoid, Cylinder, Parallelepiped.
  void SetPosDisShape(const G4String& fstring ) { fSPSPos->SetPosDisShape(fstring); }

    ///Sets the co-ordinates of the centre of the position distribution.
  void SetCentreCoords(G4ThreeVector pos) { fSPSPos->SetCentreCoords(pos); }

    ///Used to specify the co-ordinate system for the position distribution along with SetPosRot2. SetPosRot1 sets the vector x' and need not be a unit vector.
  void SetPosRot1(G4ThreeVector rot) { fSPSPos->SetPosRot1(rot); }
  
     //Used in connection with SetPosRot1. This sets a vector in the plane x'y'. By a series of cross products x', y', z' are generated. Again need not be a unit vector.
   void SetPosRot2(G4ThreeVector rot){ fSPSPos->SetPosRot2(rot); }
   
      ///Sets the radius where appropriate for source distribution shapes.
   void SetRadius(G4double rad){ fSPSPos->SetRadius(rad); }

     ///Sets the inner radius where appropriate for source distribution shapes.
   void SetRadius0(G4double rad){ fSPSPos->SetRadius0(rad); }
  
     ///Used to confine the start positions to a particular volume.
   void ConfineSourceToVolume(const G4String& vol){ fSPSPos->ConfineSourceToVolume(vol); }
  
///Set the neutrino type
   void     SetNeutrinoType(G4int k)  { fNeutrinoType = k ;}
///Get the neutrino type
   G4int    GetNeutrinoType() const   { return fNeutrinoType ;}

///Set Delta m^2
   void     SetDM2(G4double val)      { DM2 =  val ;}
///Get Delta m^2
   G4double GetDM2()          const   { return DM2 ;}
   
///Set TG2T
   void     SetTG2T(G4double val)     { TG2T = val;}
///Get TG2T
   G4double GetTG2T()         const   { return TG2T ;}
  
   void     SetBinning(G4int k)       { fNumberOfSteps = k;}
   G4int    GetBinning()      const   { return fNumberOfSteps ;}
  //protected members

///Neutrino types
  enum Neutrinos {pp=0,pep=1,hep=2,be7=3,b8=4,n13=5,o15=6,f17=7,cno=8};

///Set Tree level cross-sections
void SetOldSigma(G4bool val) {fOldSigma=val;}

  //private  members
private:
   
  ///Calculates the survive probability
  G4double survive(G4double ,  G4int ) ;
  G4double I_f (G4double ) ;
  G4double k_e (G4double ) ;
  G4double k_mu (G4double ) ;
  ///Calculates the nu_e - e  cross section 
  G4double cross_nue(G4double, G4double ) ;
  ///Calculates the nu_mu - e  cross section 
  G4double cross_nux(G4double, G4double ) ;

  G4double EnergySpectrum  (G4double, G4int ) ;
  G4double Normalizer  (G4int ) ;
  void CumulativeDistribution  (G4int ) ;
  G4double ShootEnergy()  ;
  
  G4double N13NuSpectrum (G4double ) ;
  G4double O15NuSpectrum (G4double ) ;
  G4double F17NuSpectrum (G4double ) ;
  G4double B8NuSpectrum  (G4double ) ;
  G4double PPNuSpectrum  (G4double ) ;
  G4double fMinus(G4double , G4double );
  G4double fPlus(G4double , G4double );
  G4double fPlusMinus(G4double , G4double );
//Double_t fFunct(Double_t* x, Double_t* par);
  BxGeneratorSolarNeutrinoMessenger*   fTheMessenger;
  G4ParticleGun*               fParticleGun;
  G4SPSPosDistribution*        fSPSPos;
  G4SPSAngDistribution*        fSPSAng;
  //G4SPSEneDistribution*        fSPSEne;
  G4bool                       fVolumeFlag ;
  //G4bool                       fEnergyDistribution;
  G4bool 		       fOldSigma;
  G4int                        fNumberOfSteps;
  G4double                     emass;

  G4double                     DM2;
  G4double                     TG2T;

  G4double                     AlphaPee[8];
  G4double                     BetaPee[8];
  G4double                     NuEndPoint[8];
  
  vector<G4double>             fEnergyBin;
  vector<G4double>             fProbability ;

  G4int                        fNeutrinoType;
  G4bool                       isFirstTime;
};
#endif
