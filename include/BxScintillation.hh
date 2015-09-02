//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef BXSCINTILLATION_HH
#define BXSCINTILLATION_HH

#include <vector>
#include "globals.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4OpticalPhoton.hh"

#include "BxReadParameters.hh"

class BxLightSource;
///It manages the scintillation process
class BxScintillation : public G4VRestDiscreteProcess
{
public:
///Constructor
	BxScintillation(const G4String& processName = "Scintillation", 
				const G4int& particlePDG = 11);
	
///Destructor
	~BxScintillation();
	
///Scintillation can be applied? Always true unless PDGEncoding==50
        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);	
	///FIXME: unused?	
	G4double GetMeanFreePath(const G4Track& aTrack, G4double , G4ForceCondition* );
	///FIXME: unused?	
	G4double GetMeanLifeTime(const G4Track& aTrack, G4ForceCondition* );
	///It performs the scintillation process
	G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step&  aStep);
	///It returns PostStepDoIt	
	G4VParticleChange* AtRestDoIt (const G4Track& aTrack,  const G4Step& aStep);

	void SetTrackSecondariesFirst(const G4bool state);
	G4bool GetTrackSecondariesFirst() const;				       					

	///FIXME: unused?	
	void			SetScintillationYieldFactor(const G4double YieldFactor);
	///FIXME: unused?	
	G4double		GetScintillationYieldFactor() const;

	///FIXME: unused?	
	void			SetScintillationExcitationRatio(const G4double ExcitationRatio);
	///FIXME: unused?	
	G4double		GetScintillationExcitationRatio() const;	

	///FIXME: unused?	
	void  SetMaxParentID(G4int val) { fMaxParentId = val;}
	///FIXME: unused?	
	G4int GetMaxParentId() const    {return fMaxParentId ;}			       					
	
private:
	/// It computes the quenching array for positrons and electrons
	void            BetaQuenchingFactorVector  (G4double, G4double, G4double) ;
	/// It computes the proton quenching array 
	void            ProtonQuenchingFactorVector (G4double , G4double, G4double ) ;
	/// It computes the proton quenching array 
	void            AlphaQuenchingFactorVector (G4double, G4double) ;
	/// It interpolates the values of the previous arrays 
	G4double	Interpolator (G4double, vector<G4double>&, vector<G4double>&) ;
	/// Max parent Id for quenching
	G4int           fMaxParentId;
	G4bool 		fTrackSecondariesFirst;

        G4double 	fScintillationYield;
        G4double 	fScintillationYieldBulk;
        G4double 	fScintillationYieldBuffer;
	/// Quenching in the buffer liquid
	G4double 	fBufferQuenching;
	///FIXME: unused?	
	G4double 	fResolutionScale;
	///FIXME: unused?	
	G4double	fYieldFactor;
	///FIXME: unused?	
	G4double	fExcitationRatio;
	///Birks constant for electron and positrons	
	G4double	fBirksConstantBeta;
	///It is set to zero in BorexProperty file
	G4double	fBirksConstantBeta2;
	///Birks constant for alphas	
	G4double	fBirksConstantAlpha;
	///It is set to zero in BorexProperty file
	G4double	fBirksConstantAlpha2;
	///Birks constant for protons	
	G4double	fBirksConstantProton;
	///It is set to zero in BorexProperty file
	G4double	fBirksConstantProton2;	
	///Endpoint for quenching calculation	
	G4double        fEMinusEPlusEndPoint;
	///Endpoint for quenching calculation	
	G4double        fAlphaEndPoint;
	///Endpoint for quenching calculation	
	G4double        fProtonEndPoint;

	///Vector of energies for quenching calculation
	vector<G4double>	fEMinusEPlusEnergy;
	///Vector of energies for quenching calculation
	vector<G4double>	fAlphaEnergy;
	///Vector of energies for quenching calculation
	vector<G4double>	fProtonEnergy;
	///The quenching for electrons
	vector<G4double>	fBetaFactor;
	///The quenching for positrons 
	vector<G4double>	fBetaPlusFactor;
	///The quenching for alphas 
	vector<G4double>	fAlphaFactor;
	///The quenching for protons 
	vector<G4double>	fProtonFactor;
	///The quenching calculated at the particle energy
	G4double	fQuenchingFactor;

	///FIXME: unused?	
	G4int		NumAlphaSteps;
///Zero unless alpha or proton
	G4int		fParticleCode;

	G4bool          fComputeBetaBirks;
	G4bool          fComputeAlphaBirks;
	G4bool          fComputeProtonBirks;

	//BxLightSource*	fLightSource;

	G4int   	fParticlePDG;
	G4int           fNumberOfStep ;
	///Mean energy necessary to excitate the scintillator: FIXME: references?
	G4double        fMeanExcitationEnergy;
	///Electron mass
	G4double        emass  ;
	///Alpha mass
	G4double        amass  ;
	///Proton mass
	G4double        pmass  ;	
	/// coefficient in the Bethe equation: fBetheCoeff = K/2 * ro * Z/A
	G4double        fBetheCoeff ;

	///FIXME: unused?	
	G4int           fScintillatorIndex ;
	///FIXME: unused?	
	G4int           fDMPBufferIndex ;

};

inline
G4bool BxScintillation::IsApplicable(const G4ParticleDefinition& aParticleType)
{
        if (aParticleType.GetPDGEncoding() == 50)
           return false;
        else
           return true;
           
}

inline 
void BxScintillation::SetTrackSecondariesFirst(const G4bool state) 
{ 
	fTrackSecondariesFirst = state;
}

inline
G4bool BxScintillation::GetTrackSecondariesFirst() const
{
        return fTrackSecondariesFirst;
}

inline 	void BxScintillation::SetScintillationYieldFactor(const G4double YieldFactor)
{
	fYieldFactor = YieldFactor;
}

inline	G4double BxScintillation::GetScintillationYieldFactor() const
{
	return fYieldFactor;
}

inline 	void BxScintillation::SetScintillationExcitationRatio(const G4double ExcitationRatio)
{
	fExcitationRatio = ExcitationRatio;
}

inline	G4double BxScintillation::GetScintillationExcitationRatio() const
{
	return fExcitationRatio;
}					

#endif
