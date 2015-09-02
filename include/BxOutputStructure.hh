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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#ifndef BxOutputStructure_h
#define BxOutputStructure_h 1

#include <vector>
#include "G4ThreeVector.hh"
using namespace std;

struct BinaryStructure {
  int NRUN;
  int NEVENTS;
  float TEMPOSTAT;
  int NEVENTO;
  int NUMDO1;
  int JDECAD1;
  int JDECAD2;
  int ISOTOPO1;
  int ISOTOPO2;
  int ICOINC1;
  int ICOINC2;
  double ABSTIME1;
  double ABSTIME2;
  float EIN1;
  float EIN2;
  float PSIN1[3];
  float PSIN2[3];
  int ITIPO1;
  int ITIPO2;
  int NMAX1;
  int NMAX2;
  int NUMCER1;
  int NUMCER2;
  float PLIGHT[3];
  int NPH;
  vector<int>   IPHS;
  vector<float> TPHS; 
};  

struct BinaryStructureFake {
  int NRUN;
  int NEVENTS;
  float TEMPOSTAT;
  int NEVENTO;
  int NUMDO1;
  int JDECAD1;
  int JDECAD2;
  int ISOTOPO1;
  int ISOTOPO2;
  int ICOINC1;
  int ICOINC2;
  double ABSTIME1;
  double ABSTIME2;
  float EIN1;
  float EIN2;
  float PSIN1[3];
  float PSIN2[3];
  int ITIPO1;
  int ITIPO2;
  int NMAX1;
  int NMAX2;
  int NUMCER1;
  int NUMCER2;
  float PLIGHT[3];
  int NPH;
}; 


struct HeaderStructure {
  int     Events;
  int     Run;
  int     PDG; 
  int     Chain;  
  int     Generator;   
  float   Rate;
  int     DetectorFlag;
  int     PMTFlag;
  int     PhotonYield;
  float   KB;
  float   KB2;
  int     SpatialDist; 
  float   RagMin;
  float   RagMax;
  float   NuSin2;
  float   NuDeltaMass;
  int     CommentLength;
  char    Comment[10000];
}; 

struct UserStructure {
  int     UserInt1;
  int     UserInt2;
  float   UserFloat1;
  float   UserFloat2;
  double  UserDouble;
}; 

struct PhotonStructure {
  int     Hole;
  float   TOF;
}; 

struct MuPhotonStructure {
  int     MuHole;
  float   MuTOF;
}; 


struct DepositStructure {
  int    PDGParent;
  float  Energy;
  float  Position[3];
};  

struct DaughterStructure {
  int    Id;
  int    PDG;
  double Time;
  float  Energy;
  float  Position[3];
  float  Direction[3];
}; 

struct VertexStructureDiskFormat {
  int    EventID;                   // Event ID 
  int    NSequence;                 // Sequence number of isotope in the chain
  int    IsotopeCoinc;              // 1 if there is a coincidence
  int    PDG;                       // Particle Data Group Code 
  double Time;                      // Time
  float  Energy;                    // Kinetic energy 
  float  VisEnergy;                 // Visible Kinetic energy 
  float  Position[3];               // Position
  float  Barycenter[3];             // Barycenter
  float  Direction[3];              // Direction
  int    NPE;                       // Number of photoelectrons in the ID (size of the vector thePhotons)
  int    MuNPE;                     // Number of photoelectrons in the OD (size of the vector theMuPhotons)
  int    NDaughters;                // Number of daughters (size of the vector theDaughters)
  int    NDeposits;                 // Number of deposits  (size of the vector theDeposits)
  int    NUsers;                    // Size of the vector theUsers
  int    NPhotons;                  // Number of generated photons in the ID
  int    NMuPhotons;                // Number of generated photons in the OD
}; 

struct VertexStructure: public VertexStructureDiskFormat {
  vector<DaughterStructure>   theDaughters ;
  vector<DepositStructure>    theDeposits  ;
  vector<UserStructure>       theUsers     ;
  vector<PhotonStructure>     thePhotons   ;
  vector<MuPhotonStructure>   theMuPhotons ;
};

/// Structure for the MultiEvent generator
struct MultiEvent {
  int    Id;
  float  BRTOT;
  int    PDG;
  float  Energy;
  float  BR;
}; 

#endif
