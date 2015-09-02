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
#include "BxG4BxReader.hh"
#include "BxIO.hh"
#include "BxLogger.hh"
#include "BxOutputVertex.hh"
using namespace std;


BxG4BxReader* BxG4BxReader::me = 0;

// singleton
BxG4BxReader::BxG4BxReader() {
}

BxG4BxReader* BxG4BxReader::Get() {
  if (!me)
    me = new BxG4BxReader();
  return me;
}


G4bool BxG4BxReader::ReadEvent() { 
  int BuffDimension;
  BxIO::Get()->GetG4BxFile().read(reinterpret_cast<char*>(&BuffDimension),  sizeof ( int ));

  
  BxIO::Get()->GetG4BxFile().read(reinterpret_cast<char*>(&fVertex),  sizeof ( VertexStructureDiskFormat ));  
  BxOutputVertex::Get()->SetVertexStructureDiskFormat(fVertex);

  BxOutputVertex::Get()->SetEventID(fVertex.EventID) ;		 
  BxOutputVertex::Get()->SetPDG(fVertex.PDG)	 ; 
  BxOutputVertex::Get()->SetNSequence(fVertex.NSequence) ;  	 
  BxOutputVertex::Get()->SetIsotopeCoinc(fVertex.IsotopeCoinc) ;	
  BxOutputVertex::Get()->SetTime(fVertex.Time)	 ;	 
  BxOutputVertex::Get()->SetEnergy(fVertex.Energy)	 ;	 
  BxOutputVertex::Get()->SetVisEnergy(0) ;	
  BxOutputVertex::Get()->SetPosition(BxOutputVertex::Get()->CopyArrayToVector(fVertex.Position)); 
  BxOutputVertex::Get()->SetBarycenter(BxOutputVertex::Get()->CopyArrayToVector(fVertex.Barycenter) ); 
  BxOutputVertex::Get()->SetDirection(BxOutputVertex::Get()->CopyArrayToVector(fVertex.Direction ) ) ;
  BxOutputVertex::Get()->SetNDaughters(fVertex.NDaughters) ; 	
  BxOutputVertex::Get()->SetNDeposits(0)   ;	 
  BxOutputVertex::Get()->SetNPE(0) ;	 
  BxOutputVertex::Get()->SetMuNPE(0) ;		 
  BxOutputVertex::Get()->SetNPhotons(0) ;		 
  BxOutputVertex::Get()->SetNMuPhotons(0) ; 	
  BxOutputVertex::Get()->SetNUsers(fVertex.NUsers) ;		 
		       
  for(int i=0; i<fVertex.NDaughters;i++) {
    BxIO::Get()->GetG4BxFile().read(reinterpret_cast<char*>(&fDaughter),  sizeof ( DaughterStructure ));
    BxOutputVertex::Get()->SetDaughter(fDaughter);  
    BxOutputVertex::Get()->SetDaughters();  
  }

  for(int i = 0; i <fVertex.NDeposits; i++ ) { 
    BxIO::Get()->GetG4BxFile().read(reinterpret_cast<char*>(&fDeposit),  sizeof ( DepositStructure ));
    BxInputVertex::Get()->SetDeposit(fDeposit);   
    BxInputVertex::Get()->SetDeposits();   
  }
  
  for(int i = 0; i <fVertex.NUsers; i++ ) { 
    BxIO::Get()->GetG4BxFile().read(reinterpret_cast<char*>(&fUser),  sizeof ( UserStructure ));
    BxOutputVertex::Get()->SetUser(fUser);   
    BxOutputVertex::Get()->SetUsers();   
  }

  for(int i = 0; i <fVertex.NPE; i++ ) { 
    BxIO::Get()->GetG4BxFile().read(reinterpret_cast<char*>(&fPhoton),  sizeof ( PhotonStructure ));
    BxOutputVertex::Get()->SetPhoton(fPhoton);   
    BxOutputVertex::Get()->SetPhotons();   
  }

  for(int i = 0; i <fVertex.MuNPE; i++ ) { 
    BxIO::Get()->GetG4BxFile().read(reinterpret_cast<char*>(&fMuPhoton),  sizeof ( MuPhotonStructure ));
    BxOutputVertex::Get()->SetMuPhoton(fMuPhoton);   
    BxOutputVertex::Get()->SetMuPhotons();   
  }
  
  int BuffDimension2;
  BxIO::Get()->GetG4BxFile().read(reinterpret_cast<char*>(&BuffDimension2),  sizeof ( int ));

  if(BuffDimension == BuffDimension2) return true;
  return false;
}
void BxG4BxReader::ReadHeader() {
  int BuffDimension;
  BxIO::Get()->GetG4BxFile().read(reinterpret_cast<char*>(&BuffDimension),  sizeof ( int ));
  BxIO::Get()->GetG4BxFile().ignore(BuffDimension);
  int BuffDimension2;
  BxIO::Get()->GetG4BxFile().read(reinterpret_cast<char*>(&BuffDimension2),  sizeof ( int ));
  if(BuffDimension != BuffDimension2)
    BxLog(error) << "Error in reading the G4Bx file header!" << endl;
}


void BxG4BxReader::DumpHeader() {


}


void BxG4BxReader::DumpEvent() {


}
