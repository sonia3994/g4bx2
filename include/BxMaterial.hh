//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef BxMaterial_h
#define BxMaterial_h

#include "G4Material.hh"

/**
 * This class stores the properties of the simulated materials
*/
 class BxMaterial  {

  public:

   ~BxMaterial();

///It defines the basic atomic elements and their compounds
    void DefineMaterials();

///It defines the absorption lengths, the refraction indexes and the scintillation properties in function of the wavelength
    void DefineProperties();

    void SetMaterialIndexes(); 

    G4Material*		GetPPO()          { return fPPO ;}
    G4Material*		GetDMP()          { return fDMP ;}
    G4Material*		GetPC()           { return fPC ;}
    G4Material*		GetWater()        { return fWater;}
    G4Material*		GetGlass()        { return fGlass ;}
    
    G4Material*		GetSteel()        { return fSteel ;}
    G4Material*		GetCopper()        { return fCopper ;}
    G4Material*		GetMuMetal()        { return fMuMetal ;}

    G4Material* 	GetTyvek()	{return fTyvek;}
    G4Material*		GetNylon()        { return fNylon ;}
  //G4Material*		GetNylonL()       { return fNylonL ;}

    G4Material*		GetAir()          { return fAir ;}
    G4Material*		GetVacuum()       { return fVacuum ;}  
    G4Material*		GetBialkali()     { return fBialkali ;}
   // G4Material*		GetPaint()        { return fPaint ;}
    G4Material*		GetScintillator() { return fScintillator ;}
    G4Material*		GetDMPbuffer()    { return fDMPbuffer ;}
    G4Material*		GetQuartz()       { return fQuartz ;}
    G4Material*		GetAmBe()           { return fAmBe ;}
    G4Material*		GetPb()           { return fPb ;}
    G4Material*		GetDerlin()       { return fDerlin ;}
  /*  G4Material*		GetRock()         { return fRock ;}
    G4Material*		GetConcrete()     { return fConcrete ;}*/
    G4Material*		GetAluminum()     { return fAluminum ;}
  
///Getter of the singleton      
   static BxMaterial* Get();

 private:
///singleton
    static BxMaterial * me;
    BxMaterial();
    
    G4Material*		fPPO;
    G4Material*		fCopper;
    G4Material*		fDMP;
    G4Material*		fPC;
    G4Material*		fWater;
    G4Material*		fNylon;
//    G4Material*		fNylonL;

   G4Material*         fSteel;
   G4Material*         fMuMetal;
   G4Material*		fAir;
   G4Material*		fVacuum;  
    G4Material*		fGlass;
 G4Material* 	fTyvek; 
   // G4Material*		fPaint;
  G4Material*		fBialkali;
    G4Material*		fScintillator;
    G4Material*		fDMPbuffer;
    G4Material*         fQuartz;
    G4Material*         fAmBe;
    G4Material*         fPb;    
    G4Material*         fDerlin;    
  /*  G4Material*         fRock;    
    G4Material*         fConcrete; */   
    G4Material*         fAluminum;    
};

#endif
