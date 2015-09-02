
#ifndef _BXInputVertex_HH
#define _BXInputVertex_HH 1

#include "BxOutputStructure.hh"
#include "G4ThreeVector.hh"

#include "BxLogger.hh"
#include "G4String.hh"
#include <iostream>
#include <fstream>

using namespace std;

class BxInputVertex  {
  private:
    BxInputVertex();
  public:
    static BxInputVertex* Get();
    
    virtual ~BxInputVertex() {}
    
     
  private:
  
    static BxInputVertex *me;
    
  
  public:
  
    HeaderStructure& GetHeader()        { return theHeaderStructure               ; }
    
	   

    int&     GetEvents()	        { return theHeaderStructure.Events	  ; }	  
    int&     GetRun() 	                { return theHeaderStructure.Run	          ; }
    int&     GetParticle()	        { return theHeaderStructure.PDG	          ; }
    int&     GetGenerator()	        { return theHeaderStructure.Generator     ; }
    int&     GetChain()		        { return theHeaderStructure.Chain	  ; }
    float&   GetRate()  	        { return theHeaderStructure.Rate	  ; }	  
    int&     GetDetectorFlag()	        { return theHeaderStructure.DetectorFlag  ; }	  
    int&     GetPMTFlag() 	        { return theHeaderStructure.PMTFlag	  ; }	  
    int&     GetPhotonYield() 	        { return theHeaderStructure.PhotonYield   ; }	  
    float&   GetKB() 		        { return theHeaderStructure.KB  	  ; }	  
    float&   GetKB2() 	                { return theHeaderStructure.KB2 	  ; }	  
    int&     GetSpatialDist()           { return theHeaderStructure.SpatialDist   ; }	  
    float&   GetRagMin()  		{ return theHeaderStructure.RagMin  	  ; }	  
    float&   GetRagMax()  		{ return theHeaderStructure.RagMax  	  ; }	  
    float&   GetNuSin2()  		{ return theHeaderStructure.NuSin2  	  ; }	  
    float&   GetNuDeltaMass() 		{ return theHeaderStructure.NuDeltaMass   ; }	  
    int&     GetCommentLength()         { return theHeaderStructure.CommentLength ; }	  
    string   GetComment()               { return fComment                         ; }	
    
    
    
    VertexStructure&             GetVertex()               { return theVertexStructure             ; }   
    VertexStructureDiskFormat&   GetVertexDiskFormat()     { return theVertexStructureDiskFormat   ; }   
    
    
    int&            GetEventID()	     	 { return theVertexStructure.EventID	  ; }		
    int&            GetPDG()		     	 { return theVertexStructure.PDG	  ; }		
    int&            GetNSequence()		 { return theVertexStructure.NSequence	  ; }		
    int&            GetIsotopeCoinc()	         { return theVertexStructure.IsotopeCoinc ; }
    double&         GetTime()   		 { return theVertexStructure.Time	  ; }		
    float&          GetEnergy() 	     	 { return theVertexStructure.Energy	  ; }		
    float&          GetVisEnergy() 	     	 { return theVertexStructure.VisEnergy	  ; }		
    G4ThreeVector   GetPosition()		 { return CopyArrayToVector(theVertexStructure.Position)     ; }		
    G4ThreeVector   GetBarycenter()		 { return CopyArrayToVector(theVertexStructure.Barycenter)   ; }		   
    G4ThreeVector   GetDirection()		 { return CopyArrayToVector(theVertexStructure.Direction)    ; }		   
    int&            GetNDaughters()	         { return theVertexStructure.NDaughters   ; }		
    int&            GetNDeposits()	         { return theVertexStructure.NDeposits    ; }		
    int&            GetNPE()		         { return theVertexStructure.NPE	  ; }		
    int&            GetMuNPE()		         { return theVertexStructure.MuNPE	  ; }		
    int&            GetNPhotons()	         { return theVertexStructure.NPhotons	  ; }		    
    int&            GetNMuPhotons()	         { return theVertexStructure.NMuPhotons   ; }		    
    int&            GetNUsers()	        	 { return theVertexStructure.NUsers	  ; }		    
    
    vector<DaughterStructure>&   GetVDaughters() { return theVertexStructure.theDaughters ; }		
    vector<PhotonStructure>&     GetVPhotons()	 { return theVertexStructure.thePhotons   ; }		
    vector<MuPhotonStructure>&   GetVMuPhotons() { return theVertexStructure.theMuPhotons ; }		
    vector<DepositStructure>&    GetVDeposits()	 { return theVertexStructure.theDeposits  ; }		
    vector<UserStructure>&       GetVUsers()	 { return theVertexStructure.theUsers	  ; }		
     
    
    
    
    DaughterStructure& GetDaughters()            { return theDaughterStructure           ; }
    

    int&    GetDId()	                         { return theDaughterStructure.Id	 ; } 
    int&    GetDPDG()                            { return theDaughterStructure.PDG	 ; } 
    double& GetDTime()                           { return theDaughterStructure.Time	 ; } 
    float&  GetDEnergy()                         { return theDaughterStructure.Energy	 ; } 
    G4ThreeVector   GetDDirection()              { return CopyArrayToVector(theDaughterStructure.Position)  ; } 
    G4ThreeVector   GetDPosition()               { return CopyArrayToVector(theDaughterStructure.Direction) ; } 
   
    

    PhotonStructure&   GetPhotons()   { return thePhotonStructure       ;  }
    

    int&    GetHole()                    { return theMuPhotonStructure.MuHole    ; }
    float& GetTOF()                      { return theMuPhotonStructure.MuTOF	  ; }
   
    
    
    MuPhotonStructure&   GetMuPhotons()   { return theMuPhotonStructure           ; }
    

    int&    GetMuHole()                    { return theMuPhotonStructure.MuHole    ; }
    float& GetMuTOF()                      { return theMuPhotonStructure.MuTOF	  ; }

    
    
    DepositStructure&  GetDeposits()         { return theDepositStructure            ; }

    
    int&   GetDepPDG()                       { return theDepositStructure.PDGParent ; }    
    float& GetDepEnergy()                    { return theDepositStructure.Energy    ; }
    G4ThreeVector  GetDepPosition()          { return CopyArrayToVector(theDepositStructure.Position)  ; }

    void SetVertexStructureDiskFormat(VertexStructureDiskFormat val) {theVertexStructureDiskFormat = val;}
    void SetDaughter(DaughterStructure val)                          { theDaughterStructure = val;}   
    void SetPhoton(PhotonStructure val)                              { thePhotonStructure = val;} 
    void SetMuPhoton(MuPhotonStructure val)                          { theMuPhotonStructure = val;} 
    void SetDeposit(DepositStructure val)                            { theDepositStructure  = val;}    
    void SetUser(UserStructure val)                                  { theUserStructure = val;} 

    void    SetDeposits()	                 { theVertexStructure.theDeposits.push_back(theDepositStructure)  ; }

    UserStructure&  GetUsers()             { return theUserStructure           ; }

 
    int&     GetUserInt1()                 { return  theUserStructure.UserInt1   ; }	 
    int&     GetUserInt2()                 { return  theUserStructure.UserInt2   ; }	 
    float&   GetUserFloat1()               { return  theUserStructure.UserFloat1 ; }	 
    float&   GetUserFloat2()               { return  theUserStructure.UserFloat2 ; }	 
    double&  GetUserDouble()               { return  theUserStructure.UserDouble ; }	 

    
    
    G4bool GetWriteDeposits()            { return fWriteDeposits; }

    G4bool GetGenebFlag()                { return fGeneb        ; }

    void ClearHeader();
    void ClearVertex();
    void ClearPhoton();
    void ClearMuPhoton();
    void ClearDaughter();
    void ClearDeposit();
    void ClearUser();
    void ClearAll();

    void DumpHeader();
    void DumpVertex();
    void DumpPhoton();
    void DumpMuPhoton();
    void DumpDaughter();
    void DumpDeposit();
    void DumpUser();
    void DumpAll();

    void SetVerbosity(G4int val)  { fVerbosity = val;  }
    G4int GetVerbosity()          { return fVerbosity; }           
    
  private:
    VertexStructureDiskFormat     theVertexStructureDiskFormat;
    VertexStructure               theVertexStructure;
    PhotonStructure               thePhotonStructure;
    MuPhotonStructure             theMuPhotonStructure;
    DaughterStructure             theDaughterStructure;
    DepositStructure              theDepositStructure;
    HeaderStructure               theHeaderStructure;
    UserStructure                 theUserStructure;
    G4String                      fComment;
    G4bool                        fWriteDeposits;
    G4bool                        fGeneb;
    G4ThreeVector                 CopyArrayToVector(G4float*);
    G4int                         fVerbosity;
     
};

#endif
