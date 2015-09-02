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
// Created by D. Franco
// Revised by A. Caminata and S. Marcocci, Sept. 2014

#ifndef BxDetectorConstructionMessenger_h
#define BxDetectorConstructionMessenger_h 1


#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"			
#include "globals.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

using namespace std;

class G4UIcmdWith3VectorAndUnit;
class BxDetectorConstruction;
class G4UIcmdWithADoubleAndUnit;
///Messenger of BxDetectorConstruction
class BxDetectorConstructionMessenger: public G4UImessenger {

  public:
    BxDetectorConstructionMessenger(BxDetectorConstruction*);
   ~BxDetectorConstructionMessenger();

    void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    
    BxDetectorConstruction*        fDetectorConstruction;
    G4UIcmdWithABool*              fTankCmd;
  //  G4UIcmdWithABool*              fOPERACmd;
  //  G4UIcmdWithABool*              fRockCmd;
    G4UIcmdWithABool*              fNoPMTCmd;
    G4UIdirectory*                 fDirectory;
    G4UIcmdWithAString*            fConfCmd;
    G4UIcmdWithAString*            fPMTCmd;
    G4UIcmdWithAString*            fPMTdistribCmd;
    G4UIcmdWithAnInteger*            fPMTnumCmd;
    G4UIcmdWithADouble*            fBirksACmd;
    G4UIcmdWithADouble*            fBirksA2Cmd;
    G4UIcmdWithADouble*            fBirksBCmd;
    G4UIcmdWithADouble*            fBirksB2Cmd;
    G4UIcmdWithABool*              fNoConcCmd;
    G4UIcmdWith3VectorAndUnit*     fVesselOriginCmd;
    G4UIcmdWithABool*              fCheckOverlapCmd;
    G4UIcmdWithABool*              fCheckOverlapSourceCmd;
    //G4UIcmdWithABool*              fIsSensitiveCmd;
    G4UIcmdWithADouble*            fSSSReflectivityCmd;
    G4UIcmdWithADouble*            fSSSspecLOBECmd;
    G4UIcmdWithADouble*            fSSSspecSPIKECmd;
    G4UIcmdWithADouble*            fSSSbackSCATTCmd;
    
    G4UIcmdWithADouble*            fPMTringReflectivityCmd;
    G4UIcmdWithADouble*            fPMTringspecSPIKECmd;
    G4UIcmdWithADouble*            fConcentratorExternalReflectivityCmd;
    G4UIcmdWithADouble*            fConcentratorInternalReflectivityCmd;
    G4UIcmdWithADouble*            fConcentratorExternalspecSPIKECmd;
    G4UIcmdWithADouble*            fConcentratorInternalspecSPIKECmd;
    G4UIcmdWithADouble*            fPMTReflectivityCmd;
    G4UIcmdWithADouble*            fShieldReflectivityCmd;
    G4UIcmdWithADouble*            fNylonReflectivityCmd;

    G4UIdirectory*                 fScintillatorDirectory;
    G4UIcmdWithAString*            fAlphaDecayCmd;
    G4UIcmdWithAString*            fBetaDecayCmd;
    G4UIcmdWithAString*            fAlphaWeightCmd;
    G4UIcmdWithAString*            fBetaWeightCmd;
    G4UIcmdWithADouble*            fTimePCEmissionCmd;
    G4UIcmdWithADouble*            fTimePPOEmissionCmd;
    G4UIcmdWithADouble*            fTimePCtoPPOTransferCmd;
    
    G4UIcmdWithADouble*            fPPOAttenuationLengthFactorCmd;
    G4UIcmdWithADouble*            fPCAttenuationLengthFactorCmd;
    G4UIcmdWithADouble*            fDMPAttenuationLengthFactorCmd;
    G4UIcmdWithADouble*            fNylonAttenuationLengthFactorCmd;
    G4UIcmdWithADouble*            fNylonWlShiftCmd;
    G4UIcmdWithADouble*            fPMTShieldShiftCmd;
    G4UIcmdWithADouble*		   fPPOAbsReemProbXRaysCmd;    
    G4UIcmdWithADouble*		   fPPOAbsReemProbXRaysThCmd;    
    G4UIdirectory*                 fSourceDirectory;
    G4UIcmdWithAString*            fSourceTypeCmd;
    G4UIcmdWith3VectorAndUnit*     fSourceOriginCmd;
    G4UIcmdWithADoubleAndUnit*     fSourceRadiusCmd;
    G4UIcmdWithADoubleAndUnit*     fSourceLongCmd;
    G4UIcmdWithADoubleAndUnit*     fSourceThickCmd;
    G4UIcmdWithADoubleAndUnit*     fSourceXCmd;
    G4UIcmdWithADoubleAndUnit*     fSourceYCmd;
    G4UIcmdWithADoubleAndUnit*     fSourceZCmd;
    G4UIcmdWithAString*            fSourceMaterialCmd;
    G4UIcmdWithAString*            fSourceVialMaterialCmd;
    G4UIcmdWithABool*              fDefVesselShapeCmd;
    G4UIcmdWithAnInteger*          fRunNumberCmd;
    
    G4UIcmdWithADoubleAndUnit*     fRZDimCmd;
    G4UIcmdWithABool*              fIsRingCmd;
    G4UIcmdWithoutParameter*       fUpdateGeometryCmd;   
    
};

#endif
