//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef _BXIO_HH
#define _BXIO_HH 1

#include "G4String.hh"
#include <iostream>
#include <fstream>

using namespace std;
/**
 * This class manages the I/O operations with data files
 * */
class TFile;
///Input-output from and to file
class BxIO  {
  private:
    BxIO();
  public:
///It gets the singleton
    static BxIO* Get();
    virtual ~BxIO() {}
///It checks the filename (G4String) and if it is already used adds _v1, _v2 and so on
    G4String CheckFileName(G4String); 
    G4String GetFileName()                    { return fFileName ;} 
    void SetFileName(G4String a )             { fFileName = a ;}   
 
      /// Output Binary File getter   
    ofstream&  GetBinaryFile()                { return fBinaryFile; }
    G4String GetBinaryFileName()              { return fBinaryFileName;}
    void OpenBinaryFile(); 
    void CloseBinaryFile(); 
    G4bool GetIsBinary()                      { return fIsBinary ;} 
    void SetIsBinary(G4bool a)                { fIsBinary = a ;} 
   
      /// Log Files    
    ofstream&  GetStreamLogFile()        { return fStreamLogFile;}
    void OpenLogFile(); 
    void CloseLogFile(); 

    /// G4Bx File  (the one used to start the simulation from energy deposits for external background)  
    ifstream&  GetG4BxFile()                 { return fG4BxFile; }
    G4String GetG4BxFileName()               { return fG4BxFileName;}
    void OpenG4BxFile(); 
    void CloseG4BxFile(); 
    void   SetIsG4Bx(G4bool a)               { fIsG4Bx = a;}    
    G4bool IsG4Bx()                          { return fIsG4Bx ;} 
    void SetG4BxFileName(G4String a )        { fG4BxFileName = a ;}   

    ///Get and Set of BxProperty's file name
    void    SetBorexProperty(G4String val)      { fStreamBxPropertyFileName = val ;}
    string  GetBorexProperty()                { return fStreamBxPropertyFileName ;}

    ///Get and Set of dst's file name
    void  SetDstFileName(G4String name)         {fDstFileName=name;}       
    string GetDstFileName() 		     {return fDstFileName;}
 
    ///Get and Set of the Concentrators Shape's file name
    void SetLightGuideShapeFileName(G4String val) {fLightGuideShapeFileName=val;}
    string GetLightGuideShapeFileName()  {return fLightGuideShapeFileName;}


  private:
///the singleton
    static BxIO *me;

    G4bool fIsBinary ;
    G4bool fIsG4Bx;

    ofstream  fBinaryFile;    
    ofstream  fStreamLogFile;
    ifstream  fG4BxFile;
    
    G4String  fFileName;
    G4String  fLogFileName;
    G4String  fBinaryFileName;
    G4String  fG4BxFileName;
    G4String  fBorexPropertyFileName;
    G4String fBorexGeometryFileName;
    G4String fDstFileName;
    G4String fLightGuideShapeFileName;

    void SetLogFileName(G4String name)   { fLogFileName = name;}
    void SetBinaryFileName(G4String name)     { fBinaryFileName = name;}
    
    void ChangeName(G4String); 
    
    // Input files
  
  public:
  

    void  SetStreamPCFileName(G4String name)	         {  fStreamPCFileName = name; }
    void  SetStreamPCspectrFileName(G4String name)	 {  fStreamPCspectrFileName = name; }

    void  SetStreamPPOFileName(G4String name)	         {  fStreamPPOFileName = name; }
    void  SetStreamPPOspectrFileName(G4String name)	 {  fStreamPPOspectrFileName = name; }

    void  SetStreamDMPFileName(G4String name)	         {  fStreamDMPFileName = name; }

    void  SetStreamNylonFileName(G4String name)	         {  fStreamNylonFileName = name; }
    void  SetStreamBxPropertyFileName(G4String name)	 {  fStreamBxPropertyFileName = name; }
    void SetStreamBxGeometryFileName(G4String name)      {  fStreamBxGeometryFileName=name;   }
    void  SetStreamBxPMTFileName(G4String name)	         {  fStreamBxPMTFileName = name; }
    void  SetStreamBxPMTVetoFileName(G4String name)	 {  fStreamBxPMTVetoFileName = name; }

    void  SetStreamBxPMTRelQEFileName(G4String name)	 {  fStreamBxPMTRelQEFileName = name; }

    void  SetStreamBxQEFileName(G4String name) 	         {  fStreamBxQEFileName = name; }
    void  SetStreamBxQEWaveFileName(G4String name)	 {  fStreamBxQEWaveFileName = name; }

    void  SetStreamLGreflFileName(G4String name)	 {  fStreamLGreflFileName = name; }

    void  SetStreamBxSRIFileName(G4String name)          {  fStreamBxSRIFileName = name; }
    void  SetStreamTyvekFileName(G4String name)		{fStreamTyvekFileName=name;}
    
    ifstream&  GetStreamPC();
    ifstream&  GetStreamPCspectr();
    ifstream&  GetStreamPCAlphadEdx();
    ifstream&  GetStreamPPO();
    ifstream&  GetStreamPPOspectr();
    ifstream&  GetStreamDMP();
    ifstream&  GetStreamTyvek();
    ifstream&  GetStreamNylon();  
    ifstream&  GetStreamBxProperty();  
    ifstream& GetStreamBxGeometry();
    
    ifstream&  GetStreamBxPMT();	  
    ifstream&  GetStreamBxPMTVeto();

    ifstream&  GetStreamBxPMTRelQE();

    ifstream&  GetStreamBxQE();            
    ifstream&  GetStreamBxQEWave();	  
    ifstream&  GetStreamLGrefl(); 
 
    ifstream&  GetStreamBxSRI();            
    ifstream&  GetStreamQuartzRI();            
    
    TFile* GetStreamDst();
    void QueryDB_Dst_RelQE();

    ifstream& GetStreamLightGuideShape();

    ifstream& GetStreamC14CorrectionFactor(); 
    private:

    //    ifstream&  fCheckStreamer(ifstream*, G4String);
    ifstream  fStreamPC;    
    ifstream  fStreamPCspectr;
    ifstream  fStreamPCAlphadEdx;
    ifstream  fStreamPPO;
    ifstream  fStreamPPOspectr;    
    ifstream  fStreamDMP;    
    TFile*    fStreamDst;
    ifstream  fStreamNylon;
    ifstream  fStreamBxProperty;
    ifstream  fStreamBxGeometry;
    ifstream  fStreamBxPMT;
    ifstream  fStreamBxPMTVeto;
    ifstream  fStreamLightGuideShape;
    ifstream  fStreamBxPMTRelQE;
    ifstream  fStreamTyvek;
    ifstream  fStreamBxQE;
    ifstream  fStreamBxQEWave;

    ifstream  fStreamLGrefl;

    ifstream  fStreamBxSRI;
    ifstream  fStreamQuartzRI;
    ifstream  fStreamC14CorrectionFactor;


    G4String  fStreamPCFileName;    
    G4String  fStreamPCAlphadEdxFileName;    
    G4String  fStreamPCspectrFileName;
    G4String  fStreamPPOFileName;
    G4String  fStreamPPOspectrFileName;
    G4String  fStreamDMPFileName;

    G4String  fStreamNylonFileName;
    G4String  fStreamBxPropertyFileName;
    G4String fStreamBxGeometryFileName;
    G4String  fStreamBxPMTFileName; 
    G4String  fStreamBxPMTVetoFileName; 
    G4String  fStreamTyvekFileName;
    G4String  fStreamBxPMTRelQEFileName;

    G4String  fStreamBxQEFileName;
    G4String  fStreamBxQEWaveFileName;
    
    G4String  fStreamLGreflFileName;
    
    G4String  fStreamBxSRIFileName;
    G4String  fStreamQuartzRIFileName;
    G4String  fStreamRIFileName;
    G4String  fStreamC14CorrectionFactorFileName;
	

};
#endif
