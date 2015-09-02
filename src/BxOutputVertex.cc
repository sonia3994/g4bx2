//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#include "BxOutputVertex.hh"
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

BxOutputVertex* BxOutputVertex::me = 0;

// singleton
BxOutputVertex::BxOutputVertex(){
  fWriteDeposits         = false ;
  fWriteEB               = false ;
  fWriteEBevent          = false ;
  fGeneb                 = false ;
  fPostpone              = false ;
  fIsPostponed           = false ;
  fKillOutRadius         = 0 ;
  fVerbosity             = 0;
  fSaveIsotopes          = false;
  fTrackSecondariesFirst = false ;
  ClearAll();
}

BxOutputVertex* BxOutputVertex::Get() {
  if (!me) 
    me = new BxOutputVertex();
  return me;
}

G4ThreeVector BxOutputVertex::CopyArrayToVector( float *val) {
   G4ThreeVector vector;
   vector[0] = val[0] ;
   vector[1] = val[1] ;
   vector[2] = val[2] ;
  return vector;
}

void   BxOutputVertex::SetComment(string val) {
   fComment = val;
   theHeaderStructure.CommentLength = strlen (fComment); 
   strcpy(theHeaderStructure.Comment, fComment.c_str());
}



void  BxOutputVertex::SetPosition(G4ThreeVector  val)	   { 
  for(int i = 0; i < 3; i++) theVertexStructure.Position[i]    = val[i]; 
}

void  BxOutputVertex::SetBarycenter(G4ThreeVector val)     { 
  for(int i = 0; i < 3; i++) theVertexStructure.Barycenter[i]    = val[i]; 
}

void  BxOutputVertex::SetDirection(G4ThreeVector val)	   { 
  for(int i = 0; i < 3; i++) theVertexStructure.Direction[i]    = val[i]; 
}

void  BxOutputVertex::SetDPosition(G4ThreeVector val)       { 
  for(int i = 0; i < 3; i++) theDaughterStructure.Position[i] = val[i]  ; 
}

void  BxOutputVertex::SetDDirection(G4ThreeVector val)      { 
  for(int i = 0; i < 3; i++) theDaughterStructure.Direction[i] = val[i]  ; 
}

void  BxOutputVertex::SetDepPosition(G4ThreeVector val) { 
  for(int i = 0; i < 3; i++) theDepositStructure.Position[i]     = val[i]; 
}
void BxOutputVertex::DumpHeader() {			    
  BxLog(trace) << "Header: N    " << theHeaderStructure.Events	        << endlog;
  BxLog(trace) << "Header: Run  " << theHeaderStructure.Run 	        << endlog;
  BxLog(trace) << "Header: Rate " << theHeaderStructure.Rate	        << endlog; 
  BxLog(trace) << "Header: Det  " << theHeaderStructure.DetectorFlag    << endlog;   
  BxLog(trace) << "Header: Gen  " << theHeaderStructure.Generator	<< endlog;   
  BxLog(trace) << "Header: Ch   " << theHeaderStructure.Chain	        << endlog;   
  BxLog(trace) << "Header: PMT  " << theHeaderStructure.PMTFlag	        << endlog;   
  BxLog(trace) << "Header: PY   " << theHeaderStructure.PhotonYield     << endlog;   
  BxLog(trace) << "Header: KB   " << theHeaderStructure.KB  	        << endlog;   
  BxLog(trace) << "Header: KB2  " << theHeaderStructure.KB2 	        << endlog;   
  BxLog(trace) << "Header: SP   " << theHeaderStructure.SpatialDist     << endlog;   
  BxLog(trace) << "Header: Rmin " << theHeaderStructure.RagMin	        << endlog;  
  BxLog(trace) << "Header: Rmax " << theHeaderStructure.RagMax	        << endlog;
  BxLog(trace) << "Header: NuS  " << theHeaderStructure.NuSin2	        << endlog;
  BxLog(trace) << "Header: NuD  " << theHeaderStructure.NuDeltaMass     << endlog; 
  BxLog(trace) << "Header: CL   " << theHeaderStructure.CommentLength   << endlog; 
  BxLog(trace) << "Header: Com  " << theHeaderStructure.Comment         << endlog; 
  BxLog(trace) << endlog;      
}

void BxOutputVertex::DumpVertex(){
  BxLog(trace) << "Vertex: ID  " << theVertexStructure.EventID	        << endlog;
  BxLog(trace) << "Vertex: PDG " << theVertexStructure.PDG  	        << endlog;     
  BxLog(trace) << "Vertex: N   " << theVertexStructure.NSequence        << endlog;     
  BxLog(trace) << "Vertex: IC  " << theVertexStructure.IsotopeCoinc    << endlog;     
  BxLog(trace) << "Vertex: T   " << theVertexStructure.Time 	        << endlog;     
  BxLog(trace) << "Vertex: E   " << theVertexStructure.Energy	        << endlog; 
  BxLog(trace) << "Vertex: QE  " << theVertexStructure.VisEnergy	<< endlog; 
  BxLog(trace) << "Vertex: X   " << theVertexStructure.Position[0]      << endlog;	 
  BxLog(trace) << "Vertex: Y   " << theVertexStructure.Position[1]      << endlog;	 
  BxLog(trace) << "Vertex: Z   " << theVertexStructure.Position[2]      << endlog;	 
  BxLog(trace) << "Vertex: BX  " << theVertexStructure.Barycenter[0]    << endlog;	 
  BxLog(trace) << "Vertex: BY  " << theVertexStructure.Barycenter[1]    << endlog;	 
  BxLog(trace) << "Vertex: BZ  " << theVertexStructure.Barycenter[2]    << endlog;	 
  BxLog(trace) << "Vertex: pX  " << theVertexStructure.Direction[0]     << endlog;	 
  BxLog(trace) << "Vertex: pY  " << theVertexStructure.Direction[1]     << endlog;	 
  BxLog(trace) << "Vertex: pZ  " << theVertexStructure.Direction[2]     << endlog;	 
  BxLog(trace) << "Vertex: ND  " << theVertexStructure.NDaughters	<< endlog;		   
  BxLog(trace) << "Vertex: NDep" << theVertexStructure.NDeposits	<< endlog;		   
  BxLog(trace) << "Vertex: NUs " << theVertexStructure.NUsers   	<< endlog;		   
  BxLog(trace) << "Vertex: NPE " << theVertexStructure.NPE  	        << endlog;		   
  BxLog(trace) << "Vertex: NPH " << theVertexStructure.NPhotons	        << endlog;		   
  BxLog(trace) << endlog;      
}

void BxOutputVertex::DumpPhoton(){
  for(int i = 0; i < int(theVertexStructure.thePhotons.size()); i++) {
    BxLog(trace) << "Photon: Hole " << theVertexStructure.thePhotons[i].Hole << endlog;
    BxLog(trace) << "Photon: TOF  " << theVertexStructure.thePhotons[i].TOF  << endlog;
    BxLog(trace) << endlog;     
  }
}

void BxOutputVertex::DumpMuPhoton(){
  for(int i = 0; i < int(theVertexStructure.theMuPhotons.size()); i++) {
    BxLog(trace) << "MuPhoton: Hole " << theVertexStructure.theMuPhotons[i].MuHole << endlog;
    BxLog(trace) << "MuPhoton: TOF  " << theVertexStructure.theMuPhotons[i].MuTOF  << endlog;
    BxLog(trace) << endlog;     
  }
}

void BxOutputVertex::DumpDaughter(){
  for(int i = 0; i < int(theVertexStructure.theDaughters.size()); i++) {
    BxLog(trace) << "Daughter: Id  " << theVertexStructure.theDaughters[i].Id  << endlog;
    BxLog(trace) << "Daughter: PDG " << theVertexStructure.theDaughters[i].PDG << endlog;  
    BxLog(trace) << "Daughter: T   " << theVertexStructure.theDaughters[i].Time	  << endlog;  
    BxLog(trace) << "Daughter: E   " << theVertexStructure.theDaughters[i].Energy      << endlog; 
    BxLog(trace) << "Daughter: X   " << theVertexStructure.theDaughters[i].Position[0]  << endlog;     
    BxLog(trace) << "Daughter: Y   " << theVertexStructure.theDaughters[i].Position[1]  << endlog;     
    BxLog(trace) << "Daughter: Z   " << theVertexStructure.theDaughters[i].Position[2]  << endlog;     
    BxLog(trace) << "Daughter: pX  " << theVertexStructure.theDaughters[i].Direction[0] << endlog;     
    BxLog(trace) << "Daughter: pY  " << theVertexStructure.theDaughters[i].Direction[1] << endlog;     
    BxLog(trace) << "Daughter: pZ  " << theVertexStructure.theDaughters[i].Direction[2] << endlog;     
    BxLog(trace) << endlog;     
  }
}

void BxOutputVertex::DumpDeposit(){
  for(int i = 0; i < int(theVertexStructure.theDeposits.size()); i++) {
    BxLog(trace) << "Deposit: PDG " << theVertexStructure.theDeposits[i].PDGParent	   << endlog;
    BxLog(trace) << "Deposit: E   " << theVertexStructure.theDeposits[i].Energy	   << endlog;    
    BxLog(trace) << "Deposit: X   " << theVertexStructure.theDeposits[i].Position[0]   << endlog;	
    BxLog(trace) << "Deposit: Y   " << theVertexStructure.theDeposits[i].Position[1]   << endlog; 
    BxLog(trace) << "Deposit: Z   " << theVertexStructure.theDeposits[i].Position[2]   << endlog;	
    BxLog(trace) << endlog;     
  }
}

void BxOutputVertex::DumpUser(){
  for(int i = 0; i < int(theVertexStructure.theUsers.size()); i++) {
    BxLog(trace) << "User: I1  " << theVertexStructure.theUsers[i].UserInt1	   << endlog;
    BxLog(trace) << "User: I2  " << theVertexStructure.theUsers[i].UserInt2	   << endlog;  
    BxLog(trace) << "User: F1  " << theVertexStructure.theUsers[i].UserFloat1  << endlog;	 
    BxLog(trace) << "User: F2  " << theVertexStructure.theUsers[i].UserFloat2  << endlog; 
    BxLog(trace) << "User: D   " << theVertexStructure.theUsers[i].UserDouble  << endlog;	 
    BxLog(trace) << endlog;     
  }
}

void BxOutputVertex::DumpAll(){
  DumpHeader();
  DumpVertex();
  DumpPhoton();
  DumpMuPhoton();
  DumpDaughter();
  DumpDeposit();
  DumpUser();
}

void BxOutputVertex::ClearHeader() {
  theHeaderStructure.Events	     = 0;
  theHeaderStructure.Run	     = 0;
  theHeaderStructure.Rate	     = 0;
  theHeaderStructure.DetectorFlag    = 0;
  theHeaderStructure.Generator	     = 0;
  theHeaderStructure.Chain           = 0;
  theHeaderStructure.PMTFlag	     = 0;
  theHeaderStructure.PhotonYield     = 0;
  theHeaderStructure.KB 	     = 0;
  theHeaderStructure.KB2	     = 0;
  theHeaderStructure.SpatialDist     = 0;
  theHeaderStructure.RagMin          = 0;
  theHeaderStructure.RagMax          = 0;
  theHeaderStructure.NuSin2          = 0;
  theHeaderStructure.NuDeltaMass     = 0;
  theHeaderStructure.CommentLength   = 0;
  sprintf(theHeaderStructure.Comment, "KEIN KOMMENTAR");
}

void BxOutputVertex::ClearVertex()  {
  theVertexStructure.EventID	     = 0; 
  theVertexStructure.NSequence       = 0;
  theVertexStructure.IsotopeCoinc    = 0;
  theVertexStructure.PDG	     = 0;
  theVertexStructure.Time	     = 0;
  theVertexStructure.Energy	     = 0;
  theVertexStructure.VisEnergy	     = 0;
  theVertexStructure.Position[0]     = 0;
  theVertexStructure.Position[1]     = 0;
  theVertexStructure.Position[2]     = 0;
  theVertexStructure.Barycenter[0]   = 0;
  theVertexStructure.Barycenter[1]   = 0;
  theVertexStructure.Barycenter[2]   = 0;
  theVertexStructure.Direction[0]    = 0;
  theVertexStructure.Direction[1]    = 0;
  theVertexStructure.Direction[2]    = 0;
  theVertexStructure.NDaughters      = 0;
  theVertexStructure.NPE	     = 0;
  theVertexStructure.NPhotons	     = 0;
  (theVertexStructure.theDaughters).clear();
  (theVertexStructure.thePhotons).clear();    
  (theVertexStructure.theMuPhotons).clear();    
  (theVertexStructure.theDeposits).clear();  
}   

void BxOutputVertex::ClearPhoton()    {
  thePhotonStructure.Hole 	 = 0;
  thePhotonStructure.TOF  	 = 0;
}   

void BxOutputVertex::ClearMuPhoton()    {
  theMuPhotonStructure.MuHole 	 = 0;
  theMuPhotonStructure.MuTOF  	 = 0;
}   

void BxOutputVertex::ClearDaughter() {
  theDaughterStructure.Id	       = 0;
  theDaughterStructure.PDG	       = 0;
  theDaughterStructure.Time	       = 0;
  theDaughterStructure.Energy	       = 0;
  theDaughterStructure.Position[0]     = 0; 
  theDaughterStructure.Position[1]     = 0; 
  theDaughterStructure.Position[2]     = 0; 
  theDaughterStructure.Direction[0]    = 0;
  theDaughterStructure.Direction[1]    = 0;
  theDaughterStructure.Direction[2]    = 0;
}   

void BxOutputVertex::ClearDeposit()  {
  theDepositStructure.PDGParent	   = 0; 
  theDepositStructure.Energy	   = 0; 
  theDepositStructure.Position[0]  = 0;   
  theDepositStructure.Position[1]  = 0;   
  theDepositStructure.Position[2]  = 0;   
}    
void BxOutputVertex::ClearUser()  {
  theUserStructure.UserInt1        = 0;
  theUserStructure.UserInt2        = 0;
  theUserStructure.UserFloat1      = 0;
  theUserStructure.UserFloat2      = 0;
  theUserStructure.UserDouble      = 0;
}    

void BxOutputVertex::ClearAll()  {
  ClearVertex();
  ClearPhoton();
  ClearMuPhoton();
  ClearDaughter();
  ClearUser();
  ClearDeposit();
}
