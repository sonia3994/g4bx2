//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef _BXOutputVertex_HH
#define _BXOutputVertex_HH 1

#include "BxOutputStructure.hh"
#include "G4ThreeVector.hh"

#include "BxLogger.hh"
#include "G4String.hh"
#include <iostream>
#include <fstream>

using namespace std;
/** This class manages the storage of the event infos before they are stored in the binary file at
 * the end of the event.
 */
 
class BxOutputVertex  {
  private:
///Constructor
    BxOutputVertex();
  public:
///It gets the singleton
    static BxOutputVertex* Get();
    ///Destructor
    virtual ~BxOutputVertex() {}
    
     
  private:
  ///the singleton
    static BxOutputVertex *me;
    
  
  public:
  ///It returns the pointer to the header structure
    HeaderStructure& GetHeader()        { return theHeaderStructure               ; }
    
    void   SetEvents(int val)   	{ theHeaderStructure.Events          = val; }
    void   SetRun(int val)	        { theHeaderStructure.Run	     = val; }
    void   SetParticle(int val)	        { theHeaderStructure.PDG	     = val; }
    void   SetRate(float val)		{ theHeaderStructure.Rate            = val; }  
    void   SetGenerator(int val)   	{ theHeaderStructure.Generator       = val; }
    void   SetChain(int val)	        { theHeaderStructure.Chain	     = val; }
    void   SetDetectorFlag(int val)  	{ theHeaderStructure.DetectorFlag    = val; }
    void   SetPMTFlag(int val)  	{ theHeaderStructure.PMTFlag         = val; }
    void   SetPhotonYield(int val)  	{ theHeaderStructure.PhotonYield     = val; }
    void   SetKB(float val)		{ theHeaderStructure.KB              = val; }
    void   SetKB2(float val)		{ theHeaderStructure.KB2             = val; }
    void   SetSpatialDist(int val) 	{ theHeaderStructure.SpatialDist     = val; }
    void   SetRagMin(float val)		{ theHeaderStructure.RagMin	     = val; } 
    void   SetRagMax(float val)		{ theHeaderStructure.RagMax	     = val; } 
    void   SetNuSin2(float val)		{ theHeaderStructure.NuSin2	     = val; } 
    void   SetNuDeltaMass(float val)	{ theHeaderStructure.NuDeltaMass     = val; } 
    void   SetCommentLength(int val)   	{ theHeaderStructure.CommentLength   = val; }
    void   SetComment(string );
	   

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
    
    
  ///It returns the pointer to the vertex structure
    VertexStructure&             GetVertex()               { return theVertexStructure             ; }   
    VertexStructureDiskFormat&   GetVertexDiskFormat()     { return theVertexStructureDiskFormat   ; }   
    
    void    SetEventID(int val)                  { theVertexStructure.EventID	    = val; }
    void    SetPDG(int val)	                 { theVertexStructure.PDG	    = val; }
    void    SetNSequence(int val)	         { theVertexStructure.NSequence	    = val; }
    void    SetIsotopeCoinc(int val)	         { theVertexStructure.IsotopeCoinc  = val; }
    void    SetTime(double val)                  { theVertexStructure.Time	    = val; }
    void    SetEnergy(float val)                 { theVertexStructure.Energy	    = val; }
    void    SetVisEnergy(float val)              { theVertexStructure.VisEnergy	    = val; }
    void    SetPosition(G4ThreeVector   );       
    void    SetBarycenter(G4ThreeVector );     
    void    SetDirection(G4ThreeVector  ) ;     
    void    SetNDaughters(int val)               { theVertexStructure.NDaughters    = val; }
    void    SetNDeposits(int val)                { theVertexStructure.NDeposits     = val; }
    void    SetNPE(int val)	                 { theVertexStructure.NPE	    = val; }
    void    SetMuNPE(int val)	                 { theVertexStructure.MuNPE	    = val; }
    void    SetNPhotons(int val)                 { theVertexStructure.NPhotons      = val; }
    void    SetNMuPhotons(int val)               { theVertexStructure.NMuPhotons    = val; }
    void    SetNUsers(int val)                   { theVertexStructure.NUsers	    = val; }
    void    SetDaughters()	                 { theVertexStructure.theDaughters.push_back(theDaughterStructure); }
    void    SetPhotons()	                 { theVertexStructure.thePhotons.push_back(thePhotonStructure)    ; }
    void    SetMuPhotons()	                 { theVertexStructure.theMuPhotons.push_back(theMuPhotonStructure); }
    void    SetDeposits()	                 { theVertexStructure.theDeposits.push_back(theDepositStructure)  ; }
    void    SetUsers()	                         { theVertexStructure.theUsers.push_back(theUserStructure)  ; }
    
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
    
    void   SetDId(int val)	                 { theDaughterStructure.Id	  = val  ; }
    void   SetDPDG(int val)                      { theDaughterStructure.PDG	  = val  ; }
    void   SetDTime(double val)                  { theDaughterStructure.Time	  = val  ; }
    void   SetDEnergy(float val)                 { theDaughterStructure.Energy    = val  ; }
    void   SetDPosition(G4ThreeVector );       
    void   SetDDirection(G4ThreeVector );      

    int&    GetDId()	                         { return theDaughterStructure.Id	 ; } 
    int&    GetDPDG()                            { return theDaughterStructure.PDG	 ; } 
    double& GetDTime()                           { return theDaughterStructure.Time	 ; } 
    float&  GetDEnergy()                         { return theDaughterStructure.Energy	 ; } 
    G4ThreeVector   GetDDirection()              { return CopyArrayToVector(theDaughterStructure.Position)  ; } 
    G4ThreeVector   GetDPosition()               { return CopyArrayToVector(theDaughterStructure.Direction) ; } 
   
    

    PhotonStructure&   GetPhotons()   { return thePhotonStructure       ;  }
    
    void   SetHole(int val)           { thePhotonStructure.Hole     = val; }
    void   SetTOF(float val)          { thePhotonStructure.TOF      = val; }

    int&    GetHole()                    { return thePhotonStructure.Hole    ; }
    float& GetTOF()                      { return thePhotonStructure.TOF	  ; }
   
    
    
    MuPhotonStructure&   GetMuPhotons()   { return theMuPhotonStructure           ; }
    
    void   SetMuHole(int val)             { theMuPhotonStructure.MuHole      = val; }
    void   SetMuTOF(float val)            { theMuPhotonStructure.MuTOF       = val; }

    int&    GetMuHole()                    { return theMuPhotonStructure.MuHole    ; }
    float& GetMuTOF()                      { return theMuPhotonStructure.MuTOF	  ; }

    
    
    DepositStructure&  GetDeposits()         { return theDepositStructure            ; }

    void   SetDepPDG(int val)                { theDepositStructure.PDGParent    = val; }    
    void   SetDepEnergy(float val)           { theDepositStructure.Energy       = val; }
    void   SetDepPosition(G4ThreeVector) ;
    
    int&   GetDepPDG()                       { return theDepositStructure.PDGParent ; }    
    float& GetDepEnergy()                    { return theDepositStructure.Energy    ; }
    G4ThreeVector  GetDepPosition()          { return CopyArrayToVector(theDepositStructure.Position)  ; }
    


    UserStructure&  GetUsers()             { return theUserStructure           ; }

    void  SetUserInt1(int val) 	   	   { theUserStructure.UserInt1 = val   ; }   
    void  SetUserInt2(int val) 	   	   { theUserStructure.UserInt2 = val   ; }   
    void  SetUserFloat1(float val) 	   { theUserStructure.UserFloat1 = val ; }   
    void  SetUserFloat2(float val) 	   { theUserStructure.UserFloat2 = val ; }   
    void  SetUserDouble(double val)	   { theUserStructure.UserDouble = val ; }   
 
    int&     GetUserInt1()                 { return  theUserStructure.UserInt1   ; }	 
    int&     GetUserInt2()                 { return  theUserStructure.UserInt2   ; }	 
    float&   GetUserFloat1()               { return  theUserStructure.UserFloat1 ; }	 
    float&   GetUserFloat2()               { return  theUserStructure.UserFloat2 ; }	 
    double&  GetUserDouble()               { return  theUserStructure.UserDouble ; }	 

    void SetVertexStructureDiskFormat(VertexStructureDiskFormat val) {theVertexStructureDiskFormat = val;}
    void SetDaughter(DaughterStructure val)                          { theDaughterStructure = val;}   
    void SetPhoton(PhotonStructure val)                              { thePhotonStructure = val;} 
    void SetMuPhoton(MuPhotonStructure val)                          { theMuPhotonStructure = val;} 
    void SetDeposit(DepositStructure val)                            { theDepositStructure  = val;}    
    void SetUser(UserStructure val)                                  { theUserStructure = val;} 
   
    
    G4bool GetWriteDeposits()            { return fWriteDeposits; }
    void   SetWriteDeposits(G4bool val)  { fWriteDeposits = val ; }

    G4bool GetGenebFlag()                { return fGeneb        ; }
    void   SetGenebFlag(G4bool val)      { fGeneb = val         ; }

    G4bool PostponeFlag()               { return fPostpone; }
    void   SetPostponeFlag(G4bool val)  { fPostpone = val ; }

    G4bool IsPostponed()                 { return fIsPostponed; }
    void   SetIsPostponed(G4bool val)    { fIsPostponed = val ; }

    G4double GetKillRadius()                { return fKillOutRadius; }
    void     SetKillRadius(G4double val)    { fKillOutRadius = val ; }

    G4bool GetWriteIsotope()                 { return fSaveIsotopes; }
    void   SetWriteIsotope(G4bool val)       { fSaveIsotopes = val ; }

    G4bool GetWriteEB()            { return fWriteEB; }
    void   SetWriteEB(G4bool val)  { fWriteEB = val ; }

    G4bool GetWriteEBevent()            { return fWriteEBevent; }
    void   SetWriteEBevent(G4bool val)  { fWriteEBevent = val ; }
 
    G4bool GetTrackSecondariesFirst()            { return fTrackSecondariesFirst; }
    void   SetTrackSecondariesFirst(G4bool val)  { fTrackSecondariesFirst = val ; }

    G4int GetScintillatorIndex()            { return fScintillatorIndex ; }
    void  SetScintillatorIndex(G4int val)   { fScintillatorIndex = val ; }

    G4int GetDMPBufferIndex()            { return fDMPBufferIndex ; }
    void  SetDMPBufferIndex(G4int val)   { fDMPBufferIndex = val ; }

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
    G4ThreeVector                 CopyArrayToVector(G4float*);
   
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
    G4bool                        fWriteEB;
    G4bool                        fWriteEBevent;
    G4bool                        fGeneb;
    G4int                         fVerbosity;
    G4bool                        fPostpone;
    G4bool                        fIsPostponed;
    G4double                      fKillOutRadius;
    G4bool                        fSaveIsotopes;
    G4bool                        fTrackSecondariesFirst;
    G4int                         fScintillatorIndex ;
    G4int                         fDMPBufferIndex ;
};

#endif
