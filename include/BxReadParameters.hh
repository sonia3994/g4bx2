// --------------------------------------------------------------------------//
/** 
 * AUTHOR: Davide Franco
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
 */
// --------------------------------------------------------------------------//

#ifndef BXREADPARAMETERSH_H
#define BXREADPARAMETERSH_H

#include "G4ThreeVector.hh"
#include <vector>
#include "math.h"
#include <iostream>
#include "BxLogger.hh"
class TF1;
class TDirectory;
class TFile;
using namespace std;

/** This class reads parameters from files and stores them
 */

class BxReadParameters  {
  private:
    BxReadParameters();
  public:
    static BxReadParameters* Get();
    
    virtual ~BxReadParameters() {}
///It reads the detector geometry sizes
    void	ReadDetectorGeometry();
////It reads the optical properties
    void        ReadBorexProperty();
    void	DefineDetectorVariables();
    ///It reads the ID PMT positions
    void	ReadPMTPositions();
    ///It reads the OD PMT positions
    void	ReadPMTVetoPositions();
    ///It reads the shape of light concentrators
    void	ReadLGParam();
    ///It reads the vessel shape
    G4bool	ReadVesselGeometry();
  ///It gets the flag describing the detector configurations being simulated
    G4int      GetDetectorFlag()                 {return fDetectorFlag;}  
    ///It gets the number of PMTs in OD 
    G4int      GetNumberOfVetoPMT()              {return fNumberOfVetoPMT;}	  

//World
#ifdef BX_USE_OPERA
    G4double   GetWorldSizeX()                   {return fWorldSizeXOpera;}        
#else
    G4double   GetWorldSizeX()                   {return fWorldSizeX;}        
#endif 
    G4double   GetWorldSizeY()                   {return fWorldSizeY;}        
    G4double   GetWorldSizeZ()                   {return fWorldSizeZ;}        

//TankUp
    G4double   GetTankUpSizeRmax()               {return fTankUpSizeRmax;}        

//TankDawn
    G4double   GetTankDawnSizeZ()                {return fTankDawnSizeZ;}	   
    G4double   GetTankSteelThickness()		{return fTankSteelThickness;}
    G4double   GetDisalignment()                 {return fDisalignment;}

//Steel plate on the ground
   G4double   GetUpPlateSteelR()			{return fUpPlateSteelR;}
   G4double   GetUpPlateSteelZ()			{return fUpPlateSteelZ;}
   G4double   GetDownPlateSteelR()			{return fDownPlateSteelR;}
   G4double   GetDownPlateSteelZ()			{return fDownPlateSteelZ;}

   //Tyvek Thickness
    G4double  GetTyvekThickness()		{return fTyvekThickness;}

//Tyvek height w.r.t. the SSS
   G4double   GetTyvekSphereThickness()		{return fTyvekSphereThickness;}

//OD platform and SSS legs
    G4double   GetPlatformZ()                    {return fPlatformZ;}	      
    G4double   GetPlatformOuterRadius()		{return fPlatformOuterRadius;}
    G4double   GetLegCentre()			{return fLegCentre;}
    G4double   GetLegExternalRadius()		{return fLegExternalRadius;}
    G4double   GetLegThickness()		{return fLegThickness;}
    G4int      GetNlegs()			{return fNlegs;}

 //parameters for tank polycone (water and tyvek)   
    G4int GetNumOfZplanesTankTyvek()	{return fNumOfZplanesTankTyvek;}
    G4double* GetZplanesTankTyvek()		{return fZplanesTankTyvek;}
    G4double* GetRinTankTyvek()		{return fRinTankTyvek;}
    G4double* GetRoutTankTyvek()		{return fRoutTankTyvek;}

    G4int GetNumOfZplanesTank()	{return fNumOfZplanesTank;}
    G4double* GetZplanesTank()		{return fZplanesTank;}
    G4double* GetRinTank()		{return fRinTank;}
    G4double* GetRoutTank()	{return fRoutTank;}
    
 //SSS
    G4double   GetSSSExternalRadius()            {return fSSSExternalRadius;}	
    G4double   GetSSSThickness()                 {return fSSSThickness;}	

//ZoneIII
    G4double   GetZoneIIIExternalRadius()        {return fZoneIIIExternalRadius;}


//Shroud
    G4double   GetShroudExternalRadius()         {return fShroudExternalRadius;}
    G4double   GetShroudThickness()              {return fShroudThickness;}    

//ZoneII
    G4double   GetZoneIIExternalRadius()         {return fZoneIIExternalRadius;}
    G4double   GetZoneIIThickness()              {return fZoneIIThickness;}    

//Vessel
    G4double         GetVesselExternalRadius()   {return fVesselExternalRadius;}
    G4double         GetVesselThickness()        {return fVesselThickness;}    
    G4ThreeVector&   GetVesselOrigin()           {return fVesselCenter;}    
    void             SetVesselOrigin(G4ThreeVector val) {fVesselCenter = val;}    

//ZoneI
    G4double   GetZoneIExternalRadius()          {return fZoneIExternalRadius;}


///Vessels' Endcaps
   
   G4int 	GetNumOfZplanesNylonTube()	{return fNumOfZplanesNylonTube;}
   const G4double*    GetZplanesNylonTube()	{return fZplanesNylonTube;}
   const G4double*    GetRinNylonTube()		{return fRinNylonTube;}
   const G4double*    GetRoutNylonTube()	{return fRoutNylonTube;}

   G4int 	GetNumOfZplanesIVTube()		{return fNumOfZplanesIVTube;}
   const G4double*    GetZplanesIVTube()	{return fZplanesIVTube;}
   const G4double*    GetRinIVTube()		{return fRinIVTube;}
   const G4double*    GetRoutIVTube()		{return fRoutIVTube;}
 
   G4int 	GetNumOfZplanesOVTube()		{return fNumOfZplanesOVTube;}
   const G4double*    GetZplanesOVTube()	{return fZplanesOVTube;}
   const G4double*    GetRinOVTube()		{return fRinOVTube;}
   const G4double*    GetRoutOVTube()		{return fRoutOVTube;}
   
   G4double GetNylonRingRadiusOut()		{return fNylonRingRadiusOut;}
   G4double GetNylonRingRadiusIn()		{return fNylonRingRadiusIn;}
   G4double GetNylonRingHeight()		{return fNylonRingHeight;}
   
   G4double GetNylonBridgeRadius()		{return fNylonBridgeRadius;}
   G4double GetNylonBridgeHeight()		{return fNylonBridgeHeight;}

   G4double GetCopperStrutLength()		{return fCopperStrutLength;}
   G4double GetCopperStrutRadius()		{return	fCopperStrutRadius;}
   
   G4double GetIVTubeZposition()		{return fIVTubeZposition;}

//PMT (Cathode)    
    G4double   GetPMTQEMaximum()                 {return fPMTQEMaximum;}	
    
    G4double   GetPMTSizeRmin()                  {return  fPMTSizeRmin;} 
    G4double   GetPMTSizeRmax()                  {return  fPMTSizeRmax;} 
    G4double   GetPMTSizeThetaDelta()            {return fPMTSizeThetaDelta;}  
    G4double   GetPMTCathodeProjection()         {return fPMTCathodeProjection;}
    G4double   GetPMTThickness()                 {return fPMTThickness;}       
    
    //PMT
    
    G4double   GetPMTLength()          		 {return fPMTLength;}    
    G4double   GetPMTZShift()          		 {return fPMTZShift;}    

    //PMT Shield
    G4double   GetPMTShieldSizeR1min()           {return fPMTShieldSizeR1min;}    
    G4double   GetPMTShieldSizeR1max()           {return fPMTShieldSizeR1max;}    
    G4double   GetPMTShieldSizeR2min()           {return fPMTShieldSizeR2min;}    
    G4double   GetPMTShieldSizeR2max()           {return fPMTShieldSizeR2max;}    
    G4double   GetPMTShieldSizeZ()               {return fPMTShieldSizeZ;}	  
    G4double   GetPMTShieldZPosition()               {return fPMTShieldZPosition;}
    
    void   SetPMTShieldShift(G4double shift)      {
	BxLog(warning) << "You are setting an additional " << 
		shift << "mm shift of the PMT shield shift: be aware that a default 20mm shift is already implemented!" << endlog; 
	    fPMTShieldZPosition+=shift;
    }		 
    G4double   GetPMTShieldThickness()           {return fPMTShieldThickness;}    

   G4int 	GetNumOfZplanesShield()		{return fNumOfZplanesShield;}
   const G4double*    GetZplanesShield()		{return fZplanesShield;}
   const G4double*    GetRinShield()			{return fRinShield;}
   const G4double*    GetRoutShield()			{return fRoutShield;}
   
//PMT Housing
    G4double   GetPMTHousingSizeRmin()          { return fPMTHousingSizeRmin; }	     
    G4double   GetPMTHousingSizeRmax()          { return fPMTHousingSizeRmax; }	     
    G4double   GetPMTHousingShift() 	         { return fPMTHousingShift; }	     
    G4double   GetPMTHousingSizeZ()             { return fPMTHousingSizeZ; }
    G4double   GetPMTHousingThickness()            { return fPMTHousingThickness;}			    

   G4int 	GetNumOfZplanesHousing()	{return fNumOfZplanesHousing;}
   const G4double*    GetZplanesHousing()		{return fZplanesHousing;}
   const G4double*    GetRinHousing()		  	{return fRinHousing;}
   const G4double*    GetRoutHousing()		{return fRoutHousing;}
 
//PMT inox steel ring for PMTs without concentrator
   
   G4double    GetPMTringZShift()			{return fPMTringSizeShift;}
   G4int 	GetNumOfZplanesPMTring()	        {return fNumOfZplanesPMTring;}
   const G4double*    GetZplanesPMTring()		{return fZplanesPMTring;}
   const G4double*    GetRinPMTring()		  	{return fRinPMTring;}
   const G4double*    GetRoutPMTring()	        	{return fRoutPMTring;}
   

//PMT Light Concentrators
    G4double   GetPMTLightGuideZShift()                  {return  fPMTLightGuideZShift;}

//PMT simple disk radius 

    G4double   GetPMTSizeDiskR()                  {return  fPMTSizeDiskR;}
    G4double   GetPMTSizeDiskH()                  {return  fPMTSizeDiskH;}

// OD PMTs Flanges
    G4double	GetPMTmuFlangeSizeRmin()	{return fPMTmuFlangeSizeRmin;}
    G4double	GetPMTmuFlangeSizeRmax()	{return fPMTmuFlangeSizeRmax;}
    G4double	GetPMTmuFlangeSizeZ()		{return fPMTmuFlangeSizeZ;}
    
//OD PMTs Shieldings for the ones on the ground
   G4int 	GetNumOfZplanesShieldmu()	{return fNumOfZplanesShieldmu;}
 const  G4double*    GetZplanesShieldmu()		{return fZplanesShieldmu;}
 const  G4double*    GetRinShieldmu()		{return fRinShieldmu;}
 const  G4double*    GetRoutShieldmu()		{return fRoutShieldmu;}
   
// AmBe source holder construction

    G4int GetNumOfZplanesDerlin()	{return fNumOfZplanesDerlin;}
  const  G4double* GetZplanesDerlin()		{return fZplanesDerlin;}
   const  G4double* GetRinDerlin()		{return fRinDerlin;}
    const G4double* GetRoutDerlin()		{return fRoutDerlin;}
  
    //Screws for AmBe source holder
    
    G4double GetScrewHeadHeight()		{return fScrewHeadHeight;}
    G4double GetScrewHeadRadius()		{return fScrewHeadRadius;}
    G4double GetScrewBaseHeight()		{return fScrewBaseHeight;}
    G4double GetScrewBaseRadius()		{return fScrewBaseRadius;}
   
    //AmBe Pb container

    G4double	GetAmBePbRadius()		{return fAmBePbRadius;}
    G4double	GetAmBePbHeight()		{return fAmBePbHeight;}
    
   //AmBe source

    G4double	GetAmBeRadius()		{return fAmBeRadius;}
    G4double	GetAmBeHeight()		{return fAmBeHeight;}
    
    
    //Opera
/*
    G4double   GetOperaZoneSizeX()               {return fOperaZoneSizeX;}        
    G4double   GetOperaZoneSizeY()               {return fOperaZoneSizeY;}        
    G4double   GetOperaZoneSizeZ()               {return fOperaZoneSizeZ;}        

    G4double   GetOperaSMSizeX()                 {return fOperaSMSizeX;}        
    G4double   GetOperaSMSizeY()                 {return fOperaSMSizeY;}        
    G4double   GetOperaSMSizeZ()                 {return fOperaSMSizeZ;}   

    G4double   GetOperaTTSizeX()                 {return fOperaTTSizeX;}        
    G4double   GetOperaTTSizeY()                 {return fOperaTTSizeY;}        
    G4double   GetOperaTTSizeZ()                 {return fOperaTTSizeZ;}   

    G4double   GetOperaMSSizeX()                 {return fOperaMSSizeX;}        
    G4double   GetOperaMSSizeY()                 {return fOperaMSSizeY;}        
    G4double   GetOperaMSSizeZ()                 {return fOperaMSSizeZ;}   

    G4double   GetOperaTTWallsSpacing()          {return fOperaTTWallsSpacing;}				     
    G4double   GetOperaMSWallsSpacing()          {return fOperaMSWallsSpacing;}
    G4double   GetOperaTTWallsFirstPos()         {return fOperaTTWallsFirstPos;}				     
    G4double   GetOperaMSWallsFirstPos()         {return fOperaMSWallsFirstPos;}    

    G4ThreeVector GetOperaZonePos()              {return fOperaZonePos;}
    G4ThreeVector GetOperaSM1Pos()               {return fOperaSM1Pos;}
    G4ThreeVector GetOperaSM2Pos()               {return fOperaSM2Pos;}
*/
// Scintillator light yield and Birks constant
    void     SetLightYield(G4int  Yield)            {  fLightYield = Yield ; }
    G4int    GetLightYield()                        {  return fLightYield  ; }
   
    void     SetBirksAlpha(G4double BirksAlpha)          {  fBirksAlpha = BirksAlpha ; }
    G4double GetBirksAlpha()                             {  return fBirksAlpha  ; }

    void     SetBirksSecondOrderAlpha(G4double BirksAlpha)  {  fBirksSecondOrderAlpha = BirksAlpha ; }
    G4double GetBirksSecondOrderAlpha()                     {  return fBirksSecondOrderAlpha  ; }

    void     SetBirksBeta(G4double BirksBeta)           {  fBirksBeta = BirksBeta ; }
    G4double GetBirksBeta()                             {  return fBirksBeta  ; }

    void     SetBirksSecondOrderBeta(G4double BirksBeta)  {  fBirksSecondOrderBeta = BirksBeta ; }
    G4double GetBirksSecondOrderBeta()                     {  return fBirksSecondOrderBeta  ; }

    void     SetBirksProton(G4double BirksProton)          {  fBirksProton = BirksProton ; }
    G4double GetBirksProton()                             {  return fBirksProton  ; }

    void     SetBirksSecondOrderProton(G4double BirksProton)  {  fBirksSecondOrderProton = BirksProton ; }
    G4double GetBirksSecondOrderProton()                     {  return fBirksSecondOrderProton  ; }



/// Light Yield and cerenkov mean free number of photons per step scale reduction    
    void     SetLightYieldScale(G4double scale)     {  fLightYieldScale = scale ; }
    G4double GetLightYieldScale()                   {  return fLightYieldScale  ; }
   
    void     SetLightCerenkovScale(G4double scale)  {  fLightCerenkovScale = scale ; }
    G4double GetLightCerenkovScale()                {  return fLightCerenkovScale  ; }

    std::vector<G4ThreeVector>&  GetPMTPosition()         { return fPMTPosition;}
    std::vector<G4int>&          GetVetoPmt()             { return fVetoPmt;}
    std::vector<G4ThreeVector>&  GetPMT_VetoOD_Position() { return fPMT_VetoOD_Position;}
    std::vector<G4int>&		 GetVetoOD_Pmt()          { return fVetoOD_Pmt;}
    G4double* 	 GetLGZValue()              { return fLGZValue;}
    G4double* 	 GetLGInnerRadius()            { return fLGInnerRadius;}
    G4double* 	 GetLGMidRadius()            { return fLGMidRadius;}
    G4double* 	 GetLGOuterRadius()            { return fLGOuterRadius;}
    
    void SetReemission(G4bool val) { IsReemission = val; }
    void SetScattering(G4bool val) { IsScattering = val; }

    G4bool GetReemission() { return IsReemission; }
    G4bool GetScattering() { return IsScattering; }

    void SetPCAttenuationLengthFactor(G4double val) { fPCAttenuationLengthFactor = val; }
    G4double GetPCAttenuationLengthFactor()        { return fPCAttenuationLengthFactor; }

    void SetPPOAttenuationLengthFactor(G4double val) { fPPOAttenuationLengthFactor = val; }
    G4double GetPPOAttenuationLengthFactor()        { return fPPOAttenuationLengthFactor; }

    void SetDMPAttenuationLengthFactor(G4double val) { fDMPAttenuationLengthFactor = val; }
    G4double GetDMPAttenuationLengthFactor()        { return fDMPAttenuationLengthFactor; }
   
    void SetNylonAttenuationLengthFactor(G4double val) { fNylonAttenuationLengthFactor = val; }
    G4double GetNylonAttenuationLengthFactor()        { return fNylonAttenuationLengthFactor; }

    void SetPPOAbsReemProbXRays(G4double val) { fPPOAbsReemProbXRays = val;}
    G4double GetPPOAbsReemProbXRays()        { return fPPOAbsReemProbXRays; }

    void SetPPOAbsReemProbXRaysTh(G4double val) {fPPOAbsReemProbXRaysTh = val;} 
    G4double GetPPOAbsReemProbXRaysTh()        {return fPPOAbsReemProbXRaysTh; }

    void SetRelPMTQEflag(G4int val) { fRelPMTQE = val; }
    G4int GetRelPMTQEflag()         { return fRelPMTQE; }
    
    void SetRDMDecay(G4bool val) { fRDMDecay = val; }
    G4bool IsRDMDecay()         { return fRDMDecay; }
    
    void SetRDMChain(G4bool val) { fRDMChain = val; }
    G4bool IsRDMChain()         { return fRDMChain; }

    void SetMaxParentIdForQuenching(G4int val) { fMaxParentIdForQuenching = val; }
    G4int  GetMaxParentIdForQuenching()        { return fMaxParentIdForQuenching; }

    void SetRealPDGMeanLife(G4double val) { fRealPDGMeanLife = val; }
    G4double GetRealPDGMeanLife()        { return fRealPDGMeanLife ; }

    void SetPreAbsTime(G4double val) {  fPreAbsTime = val; }
    G4double GetPreAbsTime()        { return fPreAbsTime ; }
    
    void SetIsTimeCut(G4bool val) { fIsTimeCut = val ;}
    G4bool IsTimeCut()            { return fIsTimeCut ; }

    void SetTimeCut(G4double val)      { fTimeCut = val ;}
    G4double GetTimeCut()            { return fTimeCut ; }
    
    G4bool GetIsDVessel()                  { return fDVessel; }
    void   SetIsDVessel(G4bool val)        { fDVessel = val; }


    void SetRZDim(int num)		{fRZDim=num;}
    G4int GetRZDim()			{return fRZDim;}

    G4int  GetRunNumber()		{return fRunNumber;}
    void SetRunNumber(int num)          {fRunNumber=num;}

    std::vector<G4double>& 	GetAlphaDecayTimeConstant()  { return fAlphaDecayTimeConstant; }
    std::vector<G4double>& 	GetBetaDecayTimeConstant()   { return fBetaDecayTimeConstant; }
    std::vector<G4double>& 	GetAlphaDecayWeight()        { return fAlphaDecayWeight; }
    std::vector<G4double>& 	GetBetaDecayWeight()         { return fBetaDecayWeight; }

    void   SetAlphaDecayTimeConstant(std::vector<G4double> val)  { fAlphaDecayTimeConstant = val; }	      
    void   SetBetaDecayTimeConstant(std::vector<G4double> val)   { fBetaDecayTimeConstant = val; }	      
    void   SetAlphaDecayWeight(std::vector<G4double> val) 	 { fAlphaDecayWeight = val; } 	      
    void   SetBetaDecayWeight(std::vector<G4double> val) 	 { fBetaDecayWeight = val; }    	      

    void SetTimePPOEmission(G4double val)      { fTimePPOEmission = val ;}
    G4double GetTimePPOEmission()            { return fTimePPOEmission ; }

    void SetTimePCEmission(G4double val)      { fTimePCEmission = val ;}
    G4double GetTimePCEmission()            { return fTimePCEmission ; }

    void SetTimePCtoPPOTransfer(G4double val)  	{ fTimePCtoPPOTransfer = val ;}
    G4double GetTimePCtoPPOTransfer()            { return fTimePCtoPPOTransfer ; }

    const G4double* GetVessel_Rin() {return fVessel_Rin;}
    const G4double* GetVessel_Rout() {return fVessel_Rout;}
    const G4double* GetZoneI_Z() {return fZoneI_Z;}
    const G4double* GetVessel_Z()  {return fVessel_Z;}

    std::vector<G4double> ConvertTo4DoubleVector(const char* st);

    void SetZCutsNumber(int num)      {fZCutsNumber=num;}
    G4int GetZCutsNumber()            {return fZCutsNumber;}
  
    G4int GetLGZCutsNumber()		{return fLGZCutsNumber;}

    private:
  ///The singleton
    static BxReadParameters *me;       
    
    G4int 	fLGZCutsNumber;    
    G4double*   fVessel_Rin;
    G4double*   fVessel_Rout;
    G4double*   fVessel_Z;
    G4double*   fZoneI_Z;
    TF1* VesselShape;
    TDirectory* VesselDir;
    G4bool      fRDMDecay;
    G4bool      fRDMChain;
    G4bool      IsReemission;	    	       
    G4bool      IsScattering;     	       
    G4double    fRealPDGMeanLife;
    G4double    fPreAbsTime;
    G4double    fPCAttenuationLengthFactor;
    G4double    fPPOAttenuationLengthFactor;
    G4double    fDMPAttenuationLengthFactor;
    G4double    fNylonAttenuationLengthFactor;
    G4double    fPPOAbsReemProbXRays;
    G4double    fPPOAbsReemProbXRaysTh;
    G4int       fMaxParentIdForQuenching;
    G4double	fRZDim;
    G4int      fRelPMTQE;
    G4int       fZCutsNumber;
    G4bool      fIsTimeCut;     	       
    G4double    fTimeCut;
    G4int      fDetectorFlag;	    	       
    G4int      fNumberOfVetoPMT;     	       

    //World
    G4double   fWorldSizeX;	
    //G4double   fWorldSizeXOpera;		     
    G4double   fWorldSizeY;				     
    G4double   fWorldSizeZ;				    


//TankUp
    G4double   fTankUpSizeRmax; 			     
    

//TankDawn
    G4double   fTankDawnSizeZ;  			     
    G4double   fDisalignment;
    G4double   fTankSteelThickness;
  
//parameters for the tank down (water and tyvek) polycone    
    G4int      fNumOfZplanesTankTyvek;
    G4double*  fZplanesTankTyvek;
    G4double*  fRinTankTyvek;
    G4double*  fRoutTankTyvek;
    G4int      fNumOfZplanesTank;
    G4double*  fZplanesTank;
    G4double*  fRinTank;
    G4double*  fRoutTank;

//Plate Steel on the ground

   G4double    fUpPlateSteelR;
   G4double    fUpPlateSteelZ;
   G4double    fDownPlateSteelR;
   G4double    fDownPlateSteelZ;

//Tyvek thickness
    G4double  fTyvekThickness;

//Tyvek height w.r.t. SSS
   G4double   fTyvekSphereThickness;

//OD Platform and SSS Legs
    G4double   fPlatformZ;
    G4double   fPlatformOuterRadius;
    G4double   fLegCentre;
    G4double   fLegThickness;
    G4double   fLegExternalRadius;
    G4int fNlegs;
    
    //SSS

    G4double   fSSSExternalRadius;			     
    G4double   fSSSThickness;				    

//ZoneIII

    G4double   fZoneIIIExternalRadius;  		     


//Shroud

    G4double   fShroudExternalRadius;			     
    G4double   fShroudThickness;			    

//ZoneII

    G4double   fZoneIIExternalRadius;			     
    G4double   fZoneIIThickness;			    

//Vessel

    G4double   fVesselExternalRadius;			     
    G4double   fVesselThickness;			    
    G4ThreeVector   fVesselCenter;			    

//ZoneI

    G4double   fZoneIExternalRadius;			     

//Vessels' endcaps
   G4double	fNylonRingRadiusOut;
   G4double	fNylonRingRadiusIn;
   G4double 	fNylonRingHeight;
   G4double	fIVTubeZposition;
   G4double	fCopperStrutLength;
   G4double	fCopperStrutRadius;
   
   G4double	fNylonBridgeRadius;
   G4double	fNylonBridgeHeight;
   
   G4int	fNumOfZplanesNylonTube;
   G4double*	fZplanesNylonTube;
   G4double*	fRinNylonTube;
   G4double*	fRoutNylonTube;
  
   G4int	fNumOfZplanesIVTube;
   G4double*	fZplanesIVTube;
   G4double*	fRinIVTube;
   G4double*	fRoutIVTube;

   G4int	fNumOfZplanesOVTube;
   G4double*	fZplanesOVTube;
   G4double*	fRinOVTube;
   G4double*	fRoutOVTube;

   //PMT (Cathode)
    G4double   fPMTQEMaximum;	
    G4double   fPMTSizeRmin;				     
    G4double   fPMTSizeRmax;				     
    G4double   fPMTSizeThetaDelta;			    
    G4double   fPMTCathodeProjection;			    
    G4double   fPMTThickness;				    

    G4double   fPMTFaceSizeRmin;			     
    G4double   fPMTFaceSizeRmax;			     
    G4double   fPMTFaceSizeThetaDelta;  		    
    G4double   fPMTFaceProjection;			    

    G4double   fPMTMidSizeRmin; 			     
    G4double   fPMTMidSizeRmax; 			     
    G4double   fPMTFrameThickness;			    
    G4double   fPMTLength;
    G4double   fPMTZShift;

    //PMT Shield

    G4double   fPMTShieldSizeR1min;			    
    G4double   fPMTShieldSizeR1max;			    
    G4double   fPMTShieldSizeR2min;			    
    G4double   fPMTShieldSizeR2max;			    
    G4double   fPMTShieldSizeZ; 			    
    G4double   fPMTShieldThickness;			    
    G4double   fPMTShieldZPosition;

   G4int	fNumOfZplanesShield;
   G4double*	fZplanesShield;
   G4double*	fRinShield;
   G4double*	fRoutShield;

//PMT Housing

    G4double   fPMTHousingSizeRmin;			     
    G4double   fPMTHousingSizeRmax;			     
    G4double   fPMTHousingSizeZ;			     
    G4double   fPMTHousingThickness;			    
    G4double   fPMTHousingShift;			    
   
    G4int	fNumOfZplanesHousing;
   G4double*	fZplanesHousing;
   G4double*	fRinHousing;
   G4double*	fRoutHousing;

//PMT ring for PMTs without concentrator
    
    G4double   fPMTringSizeR1;			     
    G4double   fPMTringSizeR2;			     
    G4double   fPMTringSizeR3;			     
    G4double   fPMTringSizeZ1;			     
    G4double   fPMTringSizeZ2;			     
    G4double   fPMTringSizeThickness;			    
    G4double   fPMTringSizeShift;			    
   
   G4int	fNumOfZplanesPMTring;
   G4double*	fZplanesPMTring;
   G4double*	fRinPMTring;
   G4double*	fRoutPMTring;


   //Light Concentrators
   
    G4double   fPMTLightGuideZShift;			    

    //PMT simple disk radius 

    G4double   fPMTSizeDiskR;
    G4double   fPMTSizeDiskH;

//OD PMTs Flanges
    G4double	fPMTmuFlangeSizeRmin;
    G4double	fPMTmuFlangeSizeRmax;
    G4double	fPMTmuFlangeSizeZ;

//OD PMTs Shields for the ones on the ground
   G4int	fNumOfZplanesShieldmu;
   G4double*	fZplanesShieldmu;
   G4double*	fRinShieldmu;
   G4double*	fRoutShieldmu;

  //AmBe holder
    G4int      fNumOfZplanesDerlin;
    G4double*  fZplanesDerlin;
    G4double*  fRinDerlin;
    G4double*  fRoutDerlin;
  
  //Screws for Derlin holder
    G4double   fScrewHeadHeight;
    G4double   fScrewHeadRadius;
    G4double   fScrewBaseHeight;
    G4double   fScrewBaseRadius;
  
    //AmBe Pb container
    G4double fAmBePbRadius;
    G4double fAmBePbHeight;
   
    //AmBe source
    G4double fAmBeRadius;
    G4double fAmBeHeight;

   /*
//Opera
    G4double   fOperaZoneSizeX;				     
    G4double   fOperaZoneSizeY;				     
    G4double   fOperaZoneSizeZ;	

    G4double   fOperaSMSizeX;				     
    G4double   fOperaSMSizeY;				     
    G4double   fOperaSMSizeZ;	

    G4double   fOperaTTSizeX;				     
    G4double   fOperaTTSizeY;				     
    G4double   fOperaTTSizeZ;	

    G4double   fOperaMSSizeX;				     
    G4double   fOperaMSSizeY;				     
    G4double   fOperaMSSizeZ;

    G4double   fOperaTTWallsSpacing;				     
    G4double   fOperaMSWallsSpacing;
    G4double   fOperaTTWallsFirstPos;				     
    G4double   fOperaMSWallsFirstPos;    

    G4ThreeVector fOperaZonePos;
    G4ThreeVector fOperaSM1Pos;
    G4ThreeVector fOperaSM2Pos;
*/
    std::vector<G4ThreeVector>  fPMTPosition;
    std::vector<G4int>          fVetoPmt;
    std::vector<G4ThreeVector>  fPMT_VetoOD_Position;
    std::vector<G4int>		fVetoOD_Pmt;
    G4double* 	fLGZValue;
    G4double* 	fLGInnerRadius;
    G4double* 	fLGMidRadius;
    G4double* 	fLGOuterRadius;
    
    G4int      fLightYield;

    G4double   fBirksBeta;
    G4double   fBirksAlpha;
    G4double   fBirksProton;
    G4double   fBirksSecondOrderBeta;
    G4double   fBirksSecondOrderAlpha;
    G4double   fBirksSecondOrderProton;    
    //G4bool     fElectronQuenching;

    G4double   fLightYieldScale ;
    G4double   fLightCerenkovScale;
    G4bool                      fDVessel;
    G4String                    fVesselFileName;
    G4int			fRunNumber;
    std::vector<G4double> 	fAlphaDecayTimeConstant;
    std::vector<G4double> 	fBetaDecayTimeConstant;
    std::vector<G4double> 	fAlphaDecayWeight;
    std::vector<G4double> 	fBetaDecayWeight;

    G4double	fTimePPOEmission;
    G4double	fTimePCEmission;
    G4double 	fTimePCtoPPOTransfer;

};



#endif
