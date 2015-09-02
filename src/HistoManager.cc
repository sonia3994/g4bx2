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

#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include <TGraph.h>
#include <TAxis.h>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager* HistoManager::me=0;


HistoManager* HistoManager::Get(){
	if(!me)
		me= new HistoManager();
	return me;
}

HistoManager::HistoManager()
:fRootFile(0), 
 fNtuple1(0), 
 fE(0), ftime(0)
{
      
  // histograms
  for (G4int k=0; k<MaxHisto; k++) fHisto[k] = 0;
    
  // ntuple
  fNtuple1 = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  if ( fRootFile ) delete fRootFile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{ 
 // Creating a tree container to handle histograms and ntuples.
 // This tree is associated to an output file.
 //
 //G4String fileName = "AnaEx02.root";
// fRootFile = new TFile(fileName,"RECREATE");
// if(!fRootFile) {
//   G4cout << " HistoManager::book :" 
  //        << " problem creating the ROOT TFile "
   //       << G4endl;
 //  return;
// }
   
 fHisto[1] = new TH1D("1", "# photons in PMTS", 1000, 0., 2000.);
 if (!fHisto[1]) G4cout << "\n can't create histo 1" << G4endl;
 fHisto[2] = new TH1D("2", "photon dies at r=...", 1000, 0, 800);
 if (!fHisto[2]) G4cout << "\n can't create histo 2" << G4endl;
 fHisto[3] = new TH1D("3", "Number of generated photons", 20000, 0, 20000);
 if (!fHisto[3]) G4cout << "\n can't create histo 3" << G4endl;
 fHisto[4] = new TH1D("time", "Time of flight (ns)", 1000, 0, 100);
 if (!fHisto[4]) G4cout << "\n can't create histo 4" << G4endl;  
 fHisto[5] = new TH1D("5", "Cerenkov energy", 100, 0, 30);
 if (!fHisto[5]) G4cout << "\n can't create histo 4" << G4endl;  

 // create 1st ntuple in subdirectory "tuples"
 //
 //fNtuple1 = new TTree("T","Time and Energy tree");
// fNtuple1->Branch("E", &fE, "E/D");
// fNtuple1->Branch("t", &ftime, "t/D");

 // create 2nd ntuple in subdirectory "tuples"
 //
 //fNtuple2 = new TTree("102", "TrackL");
 //fNtuple2->Branch("Labs", &fLabs, "Labs/D");
 //fNtuple2->Branch("Lgap", &fLgap, "Lgap/D");

 
 G4cout << "\n----> Histogram file is opened in " << fileName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{ 
  if (fRootFile) {
    fRootFile->Write();       // Writing the histograms to the file
    fRootFile->Close();        // and closing the tree (and the file)
    G4cout << "\n----> Histogram Tree is saved \n" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::CreateGraph(G4int id, G4double *xx, G4double *yy, G4int Nentries){
	
	G4double *x=new G4double [ Nentries];
	G4double *y=new G4double [ Nentries];
	for(int i=0;i<Nentries;i++){
		x[i]=xx[i];//4.13e-15*3e8/xx[i]/1e6;
		y[i]=yy[i];
	}
	if (id==0){
		graph[id]=new TGraph(Nentries,x,y);
		graph[id]->Draw("AP");
		graph[id]->SetTitle("PC Abs Length");
		xaxis=graph[id]->GetXaxis();
		xaxis->SetTitle("Lambda (m)");
		yaxis=graph[id]->GetYaxis();
		yaxis->SetTitle("Abs Length (mm)");
		graph[id]->Write("ABS_LENGTH_PC");
	}
	if (id==1){
		graph[id]=new TGraph(Nentries,x,y);
		graph[id]->Draw("AP");
		graph[id]->SetTitle("Refraction index PC");
		xaxis=graph[id]->GetXaxis();
		xaxis->SetTitle("lambda (m)");
		yaxis=graph[id]->GetYaxis();
		yaxis->SetTitle("Refraction index PC");
		graph[id]->Write("RINDEX_PC");
	}
	if (id==2){
		graph[id]=new TGraph(Nentries,x,y);
		graph[id]->Draw("AP");
		graph[id]->SetTitle("PC emission");
		xaxis=graph[id]->GetXaxis();
		xaxis->SetTitle("lambda (m)");
		yaxis=graph[id]->GetYaxis();
		yaxis->SetTitle("PC emission spectrum");
		graph[id]->Write("EMISSION_PC");
	}
	if (id==3){
		graph[id]=new TGraph(Nentries,x,y);
		graph[id]->Draw("AP");
		graph[id]->SetTitle("Re-emission probability PC");
		xaxis=graph[id]->GetXaxis();
		xaxis->SetTitle("lambda (m)");
		yaxis=graph[id]->GetYaxis();
		yaxis->SetTitle("Re-emission probability PC");
		graph[id]->Write("REEMISSION_PC");
	}
	if (id==4){
		graph[id]=new TGraph(Nentries,x,y);
		graph[id]->Draw("AP");
		graph[id]->SetTitle("Rayleigh PC");
		xaxis=graph[id]->GetXaxis();
		xaxis->SetTitle("lambda (m)");
		yaxis=graph[id]->GetYaxis();
		yaxis->SetTitle("Rayleigh");
		graph[id]->Write("RAYLEIGH_PC");
	}
	if (id==5){
		graph[id]=new TGraph(Nentries,x,y);
		graph[id]->Draw("AP");
		graph[id]->SetTitle("PPO Abs Length");
		xaxis=graph[id]->GetXaxis();
		xaxis->SetTitle("lambda (m)");
		yaxis=graph[id]->GetYaxis();
		yaxis->SetTitle("PPO Abs Length (mm)");
		graph[id]->Write("ABS_LENGTH_PPO");
	}
	if (id==6){
		graph[id]=new TGraph(Nentries,x,y);
		graph[id]->Draw("AP");
		graph[id]->SetTitle("PPO emission");
		xaxis=graph[id]->GetXaxis();
		xaxis->SetTitle("lambda (m)");
		yaxis=graph[id]->GetYaxis();
		yaxis->SetTitle("PPO emission spectrum");
		graph[id]->Write("EMISSION_PPO");
	}
	if (id==7){
		graph[id]=new TGraph(Nentries,x,y);
		graph[id]->Draw("AP");
		graph[id]->SetTitle("Re-emission probability PPO");
		xaxis=graph[id]->GetXaxis();
		xaxis->SetTitle("lambda (m)");
		yaxis=graph[id]->GetYaxis();
		yaxis->SetTitle("Re-emission probability PPO");
		graph[id]->Write("REEMISSION_PPO");
	}
	if (id==8){
		graph[id]=new TGraph(Nentries,x,y);
		graph[id]->Draw("AP");
		graph[id]->SetTitle("Scintillator RINDEX");
		xaxis=graph[id]->GetXaxis();
		xaxis->SetTitle("lambda (m)");
		yaxis=graph[id]->GetYaxis();
		yaxis->SetTitle("RINDEX");
		graph[id]->Write("RINDEX_SCINTILLATOR");
	}
	if (id==9){
		graph[id]=new TGraph(Nentries,x,y);
		graph[id]->Draw("AP");
		graph[id]->SetTitle("Scintillator Rayleigh");
		xaxis=graph[id]->GetXaxis();
		xaxis->SetTitle("lambda (m)");
		yaxis=graph[id]->GetYaxis();
		yaxis->SetTitle("Rayleigh");
		graph[id]->Write("SCINTILLATOR_RAYLEIGH");
	}
	if (id==10){
		graph[id]=new TGraph(Nentries,x,y);
		graph[id]->Draw("AP");
		graph[id]->SetTitle("Glass abslength");
		xaxis=graph[id]->GetXaxis();
		xaxis->SetTitle("lambda (m)");
		yaxis=graph[id]->GetYaxis();
		yaxis->SetTitle("abslength");
		graph[id]->Write("GLASS_ABSLENGTH");
	}
	if (id==11){
		graph[id]=new TGraph(Nentries,x,y);
		graph[id]->Draw("AP");
		graph[id]->SetTitle("Glass refraction index");
		xaxis=graph[id]->GetXaxis();
		xaxis->SetTitle("lambda (m)");
		yaxis=graph[id]->GetYaxis();
		yaxis->SetTitle("Glass refraction index");
		graph[id]->Write("GLASS_REFRACTION_INDEX");
	}
if (id==12){
                graph[id]=new TGraph(Nentries,x,y);
                graph[id]->Draw("AP");
                graph[id]->SetTitle("Nylon abslength");
                xaxis=graph[id]->GetXaxis();
                xaxis->SetTitle("lambda (m)");
                yaxis=graph[id]->GetYaxis();
                yaxis->SetTitle("abslength");
                graph[id]->Write("NYLON_ABSLENGTH");
        }
if (id==13){
                graph[id]=new TGraph(Nentries,x,y);
                graph[id]->Draw("AP");
                graph[id]->SetTitle("Nylon refraction index");
                xaxis=graph[id]->GetXaxis();
                xaxis->SetTitle("lambda (m)");
                yaxis=graph[id]->GetYaxis();
                yaxis->SetTitle("Nylon refraction index");
                graph[id]->Write("NYLON_REFRACTION_INDEX");
        }
if (id==14){
                graph[id]=new TGraph(Nentries,x,y);
                graph[id]->Draw("AP");
                graph[id]->SetTitle("Tyvek refraction index VS photon energy");
                xaxis=graph[id]->GetXaxis();
                xaxis->SetTitle("Energy (eV)");
                yaxis=graph[id]->GetYaxis();
                yaxis->SetTitle("Tyvek refraction index");
                graph[id]->Write("TYVEK_REFRACTION_INDEX");
        }


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
           << " does not exist. (xbin=" << xbin << " weight=" << weight << ")"
           << G4endl;
    return;
  }
 if  (fHisto[ih]) { fHisto[ih]->Fill(xbin, weight); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double fac)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::Normalize() : histo " << ih
           << " does not exist. (fac=" << fac << ")" << G4endl;
    return;
  }
  if (fHisto[ih]) fHisto[ih]->Scale(fac);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillNtuple(G4double energy,G4double time)
{
 fE = energy;
 ftime = time;

  if (fNtuple1) fNtuple1->Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::PrintStatistic()
{
  if(fHisto[1]) {
     G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " EAbs : mean = " << G4BestUnit(fHisto[1]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(fHisto[1]->GetRMS(),  "Energy") << G4endl;
    G4cout                
    << " EGap : mean = " << G4BestUnit(fHisto[2]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(fHisto[2]->GetRMS(),  "Energy") << G4endl;
    G4cout 
    << " LAbs : mean = " << G4BestUnit(fHisto[3]->GetMean(), "Length") 
            << " rms = " << G4BestUnit(fHisto[3]->GetRMS(),  "Length") << G4endl;
    G4cout 
    << " LGap : mean = " << G4BestUnit(fHisto[4]->GetMean(), "Length") 
            << " rms = " << G4BestUnit(fHisto[4]->GetRMS(),  "Length") << G4endl;

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


