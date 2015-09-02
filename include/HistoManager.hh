//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//Created by A. Caminata and S. Marcocci, Sept. 2014
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
//#include "TAxis.h"
//#include "TGraph.h"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 class TFile;
 class TTree;
 class TH1D;
 class TGraph;
 class TAxis;
 const G4int MaxHisto = 6;
 const G4int MaxGraph = 20;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/** This class provides an easy way to produce output rootfile to perform quick cross checks on the code
 */

class HistoManager
{
  public:
  static HistoManager* Get();
   virtual ~HistoManager();
   
    void book();
    void save();

    void FillHisto(G4int id, G4double bin, G4double weight = 1.0);
    void CreateGraph(G4int id, G4double *x, G4double *y,G4int Nentries);
    void Normalize(G4int id, G4double fac);    

    void FillNtuple(G4double fE, G4double ftime);
    
    void PrintStatistic();
    void SetFileName(G4String astr){fileName=astr;}
    G4String GetFileName(){return fileName;}

  private:
    HistoManager();
    static HistoManager *me;
    TFile*   fRootFile;
    TH1D*    fHisto[MaxHisto];            
    TTree*   fNtuple1;    
    //TTree*   fNtuple2;    
    TGraph*  graph[MaxGraph];

    TAxis *xaxis;
    TAxis *yaxis; 
    G4double fE;
    G4double ftime;

G4String fileName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

