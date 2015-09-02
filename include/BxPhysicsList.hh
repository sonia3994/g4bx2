//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014

#ifndef BxPhysicsList_h
#define BxPhysicsList_h 1

#include "globals.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4SystemOfUnits.hh"
class BxPhysicsListMessenger;

class BxOpAbsorptionReemission;

class Scintillation;
/**
 * This class is responsible to activate all the physics taken into account during the simulation
 */

class BxPhysicsList : public G4VUserPhysicsList
{
  public:
    BxPhysicsList();
    virtual ~BxPhysicsList();

    G4bool GetOPIsActivate()      const  { return OPIsActivate ;}
    G4bool GetOPcherIsActivate()  const  { return OPcherIsActivate ;}    
    G4bool GetOPscintIsActivate() const  { return OPscintIsActivate ;}
    G4bool SetAbsorptionIsActive()const  { return AbsorptionIsActive;}
    G4bool SetRayleighIsActive()  const  { return RayleighIsActive;}


    G4int  GetHadronicModel()     const  { return fHadronic ;}
    G4bool GetAltHadIsActivate()  const  { return fHadronicPhysicsListFlag; }
    G4bool GetDecayIsActivate()   const  { return DecayIsActivate; }
    G4bool GetNewEM()             const  { return IsNewEM; }
    G4bool GetNucInt()            const  { return IsNucInt; }
    

    void SetOPIsActivate(G4bool a)       { OPIsActivate = a ;}
    void SetOPcherIsActivate(G4bool a)   { OPcherIsActivate = a ;}
    void SetOPscintIsActivate(G4bool a)  { OPscintIsActivate = a ;}
    void SetAbsorptionIsActive(G4bool a) { AbsorptionIsActive = a ;}
    void SetRayleighIsActive(G4bool a)   { RayleighIsActive= a ;}
  
    void SetAltHadIsActivate(G4bool a)   { fHadronicPhysicsListFlag = a ;}
    void SetDecayIsActivate(G4bool a)    { DecayIsActivate = a ;}
    void SetDefCutValue(G4double a)      { defaultCutValue = a ;}
    void SetElectronCutValue(G4double a) { electronCutValue = a ;}
    void SetGammaCutValue(G4double a)    { gammaCutValue = a ;}
    void SetNewEM(G4bool a)              { IsNewEM = a; }
    void SetNucInt(G4bool a)             { IsNucInt = a ; }

    void SetMaxNumberOfPhotons(G4int a) { fMaxNumberOfPhoton = a ;}
    void SetNoReemission()              { IsReemission = false; }
    void SetNoScattering()              { IsScattering = false; }

    void SetHadronic(G4int a)           { fHadronic = a ;}
    void SetIonPhysics(G4bool a)         { IsIonPhys = a ;}

  protected:
    // Construct particles and processes

///It calls methods to construct particles 
virtual void ConstructParticle();

///It constructs the physical processes
    virtual void ConstructProcess();
    virtual void SetCuts();

  protected:
    /// this method construct bosons 
    void ConstructBosons();
    /// this method construct leptons 
    void ConstructLeptons();
    /// this method construct mesons 
    void ConstructMesons();
    /// this method construct barions 
    void ConstructBaryons();
    /// this method construct ions 
    void ConstructIons();
    /// this method construct short lived particles (particles with very short life time that immediately and are never tracked) 
    void ConstructShortLived();

  protected:
    // these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();
    void ConstructEM2();
///It constructs the opticals processes    
    void ConstructOp();
    void ConstructHad();
    void ConstructHad_QGSP_BERT_HP();
    //void ConstructHad_QSGP_BIC_HP();
    //void ConstructHad_FTF_BIC_HP();
    //void ConstructHad_FTFP_BERT_HP();
    void ConstructNuclearReactions();
    void ConstructIonPhysics();

    
  private:  
    
    G4bool IsScattering;
    G4bool IsReemission;
    G4bool OPIsActivate;
    G4bool OPcherIsActivate;
    G4bool OPscintIsActivate;
    G4bool AbsorptionIsActive;
    G4bool RayleighIsActive;
    G4bool DecayIsActivate;
    G4int  fHadronic;    
    
    G4int  fMaxNumberOfPhoton;
    BxPhysicsListMessenger* BxPhys;
    
    //BxOpAbsorptionReemission* theAbsorptionProcess;
   // BxOpBoundaryProcess* 	 theBoundaryProcess;

    //Scintillation*	theBetaScintillationProcess;
    //Scintillation*	theProtonScintillationProcess;    
    //Scintillation*	theAlphaScintillationProcess;    
    
    G4bool fHadronicPhysicsListFlag;
    G4double       electronCutValue ;
    G4double       gammaCutValue ; 
    
    G4bool         IsNucInt;
    G4bool         IsNewEM;
    G4bool         IsIonPhys;
};

#endif 
