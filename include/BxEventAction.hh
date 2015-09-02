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
//
// Created by D. Franco
// Revised by A. Caminata and S. Marcocci, Sept. 2014
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#ifndef BxEventAction_h
#define BxEventAction_h 1
#include "BxOutputStructure.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;
class G4Timer;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class BxEventActionMessenger;

/**This class manages the actions performed at the beginning and at the end of each event: e.g. writing on the binary file.
*/
class BxEventAction : public G4UserEventAction
{
  public:
    BxEventAction();
    virtual ~BxEventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    void SetDrawFlag   (G4String val)  {drawFlag = val;};
    void SetPrintModulo(G4int    val)  {printModulo = val;};

    void SetLogFileName(G4String  val)     {fLogFileName = val;};
    void SetBigLogFileName(G4String  val)  {fBigLogFileName = val;};
    void SetBinFileName(G4String  val)     {fBinFileName = val;};

    void SetCounter(G4int  val)        {Counter = val;};

    void SetWriteLimit(G4int  val)     { fWriteLimit = val;};

    void SetDepoRadius(G4double  val)     { fCheckDepositPosition = true ; fMaxDepositRadius = val;};
    
    void SetVisEnergyCut(G4double  val)     {  fVisibleEnergy = val;};
    
    void Set64Bit(G4bool val)  {f64bit = val;};
    

  private:
    G4String                    drawFlag;
    G4int                       printModulo;                         
    BxEventActionMessenger*     eventMessenger;
    G4String                    fBigLogFileName;
    G4String                    fLogFileName;   
    G4String                    fBinFileName ;  
  
    G4int                       Counter;                         
    G4Timer*                    timer;
    G4int                       totnpe;
    G4int                       fWriteLimit;
    G4bool                      fCheckDepositPosition ;
    G4double                    fMaxDepositRadius ;
    G4double                    fVisibleEnergy ;
    G4bool                      f64bit;
                 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
