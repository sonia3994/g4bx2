//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef BxG4BxReader_h
#define BxG4BxReader_h 1
#include "BxInputVertex.hh"
#include "BxOutputStructure.hh"

#include "globals.hh"
#include <stdio.h>
#include <iostream>

using namespace std;

class BxG4BxReader {
 private:
   BxG4BxReader();
 public:
  ~BxG4BxReader();
  static BxG4BxReader* Get();
  void    ReadHeader();
  G4bool  ReadEvent();
  void DumpHeader();
  void DumpEvent();
  void SetG4BxFile();

private:

  static BxG4BxReader *me;

  void                 SkipEvents();
  
  VertexStructureDiskFormat fVertex;
  DaughterStructure fDaughter;
  DepositStructure fDeposit;
  UserStructure fUser;
  PhotonStructure fPhoton;
  MuPhotonStructure fMuPhoton;
  
};

#endif
