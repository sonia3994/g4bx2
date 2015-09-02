//Created by D. Franco
//Major Revision by A. Caminata and S. Marcocci, Sept. 2014
#ifndef BxDETECTORCONSTRUCTION_H
#define BxDETECTORCONSTRUCTION_H
#include <vector>

#include "G4VUserDetectorConstruction.hh" 
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

class   G4VSolid;
class	G4VPhysicalVolume;
class	G4PVPlacement;

class	G4LogicalVolume;
class	G4LogicalBorderSurface;
class	G4OpticalSurface;
class   G4LogicalSkinSurface;
class   G4GenericPolycone;
class	G4Box;
class	G4Tubs;
class	G4Sphere;
class	G4Cons;
class	G4Polycone;
class	G4BREPSolidPCone;
class	G4CSGSolid;
class	G4Material;
class   BxDetectorConstructionMessenger;
class	G4SubtractionSolid;
class	G4UnionSolid;
class   G4IntersectionSolid;
class	BxMaterial;
///It defines the geometry and the physical properties of the detector
class   BxDetectorConstruction : public G4VUserDetectorConstruction {

  public:
///Constructor. It defines some defaults variables for the geometry construction
    BxDetectorConstruction();
    virtual ~BxDetectorConstruction();

///It calls DefineSizes() and DefineMaterials()
    G4VPhysicalVolume* Construct();
///If in graphic mode you change something in the geometry, you need to invoke this method in order for the Run Manager to rebuild the detector
    void UpdateGeometry();
///Output of the sizes of the detector
    void PrintDetectorParameters();
   
///It gets the dimensions of the components from BxReadParameters
    void DefineSizes();

///It gets the material properties of the components from the singleton BxMaterial
    void DefineMaterials();

///It sets the visualization attributes of the constructed volumes
    void SetVisAttributes();
    /* 
    void DefineOpera();
    
    void SetOperaDetector(G4bool opera)        { IsOpera = opera; } 
*/
    ///This function creates the calibration source
    void AddSource();
    
///Boolean Flag for the construction of the Water Tank. true=tank is built
    void SetTank(G4bool tank)                  { IsTank  = tank; } 

//    void SetRock(G4bool rock)                  { IsRock  = rock; } 
   
///It defines the flag for the detector construction.
///0=full, 1=fullWater (outer+inner), 2=only inner filled with water, 3=only outer, 4=only inner filled with scintillator, 5=ctf like
    void SetDetectorFlag(G4int flag)           { fDetectorFlag = flag; }
    
    
///It defines the flag for the PMT construction. 1=disk PMT without concentrator, 2=disk PMT with concentrator, 3=simple PMT    
   void SetPMTFlag(G4int flag)                { fPMTFlag = flag; }

///If you set "no concentrators=true", no concentrators will be built, no matter of the PMT configuration you chose
   void SetNoConcentrators(G4bool val)        { fNoConcentrators = val; }

///It defines the solid, logical and physical volumes
    G4VPhysicalVolume* ConstructDetector();

///It sets the type of PMT distribution in the inner detector. 0 means uniform distribution, while 1 means real distribution according to data file.   
    void SetPMTDistributionFlag(G4int flag)   {fPMTDistributionFlag=flag;}
   
///Number of constructed PMTs in the case of uniform distribution inside the sphere.    
    void SetNumberOfPMT(G4int num)            {fNumberOfPMT=num;}
   
///If "noPMTs" is set to true, no PMTs will be built
    void SetNoPMTs(G4bool val)                 { noPMTs = val; } 

    ///Source Type. 1=sphere, 2=cylinder, 3=Box, 4 not yet implemented
    void SetSourceType(G4int val)                  { fSource = val; }
    ///Positioning of the centre of the source
    void SetSourceOrigin(G4ThreeVector val)        { fSourceOrigin = val; } 
    ///Source radius for sphere and cylinder cases
    void SetSourceRadius(G4double val)             { fSourceRadius = val; } 
    ///Cylindrical source height
    void SetSourceLong(G4double val)               { fSourceLong = val; } ;
    ///Thickness of the vial
    void SetSourceVialThick(G4double val)          { fSourceVialThick= val; } 
    ///X-length for the box-like source
    void SetSourceX(G4double val)                  { fSourceX = val; } 
    ///Y-length for the box-like source
    void SetSourceY(G4double val)                  { fSourceY = val; } 
    ///Z-length for the box-like source
    void SetSourceZ(G4double val)                  { fSourceZ = val; } 
    ///Vial's material
    void SetSourceMaterial(G4Material* val)        { fSourceMaterial = val; }
    ///Fluid in which sources are dissolved
    void SetSourceVialMaterial(G4Material* val)    { fSourceVialMaterial = val; }
   
    ///If IsOverlap=true, the check for the volume overlaps will be performed at the beginning of the run 
    void SetOverlap (G4bool val)                   { IsOverlap = val; }
    ///This OverlapCheck is perfomed only for the volumes of the source, to speed up the calculation
    void SetOverlapSource (G4bool val)                   { IsOverlapSource = val; }
    
    //void SetSensitive (G4bool val)                 { IsSensitive = val; }
    ///This method allows you to create (true) or not to create (false) the nylon rings supporting the IV at the poles
    void SetEndCaps (G4bool val)                   { IsEndCap = val; }
  
  private:
  BxDetectorConstructionMessenger *fMessenger;

//_____World_____________________________________________

    G4double           		fWorldSizeX;
    G4double           		fWorldSizeY;
    G4double			fWorldSizeZ;
    G4Material*			fWorldMaterial;
    G4Box*			fsolidWorld;
    G4LogicalVolume*		flogicWorld;

///It contains the whole experiment
    G4VPhysicalVolume*		fphysicWorld;
/*
//_____Rock______________________________________________
    G4Box*			fsolidRock;
    G4LogicalVolume*		flogicRock;
    G4VPhysicalVolume*		fphysicRock;
    G4bool                      IsRock;
//_____Concrete___________________________________________
    G4Box*			fsolidConcrete;
    G4LogicalVolume*		flogicConcrete;
    G4VPhysicalVolume*		fphysicConcrete;
*/
//_____TankSteel________________________________________________
    G4bool                      IsTank;  
    G4LogicalVolume*		flogicTankSteel;
    G4VPhysicalVolume*		fphysicTankSteel;
    G4Material*			fTankSteelMaterial;

    G4OpticalSurface*		fOpTankSteelSurface;
    G4LogicalBorderSurface*	fTankSteelSurface;

//_____TankUpSteel________________________________________________
 
    G4Sphere*			fsolidTankUpSphereSteel;
    G4SubtractionSolid*		fsolidTankUpSteel;


//_____TankTubeSteel________________________________________________
    
    G4Tubs*			fsolidTankTubeSteel;

//_____TankTyvek________________________________________________
    G4LogicalVolume*		flogicTankTyvek;
    G4VPhysicalVolume*		fphysicTankTyvek;
    G4Material*			fTyvekMaterial;

    G4OpticalSurface*		fOpTankTyvekSurface;
    G4LogicalBorderSurface*	fTankTyvekSurface;

//_____TankUpTyvek________________________________________________
 
    G4Sphere*			fsolidTankUpSphereTyvek;
    G4SubtractionSolid*		fsolidTankUpTyvek;


//_____TankTubeTyvek________________________________________________
    
    G4Polycone*			fsolidTankTubeTyvek;
//_____Tank________________________________________________

    G4LogicalVolume*		flogicTank;
    G4VPhysicalVolume*		fphysicTank;
    G4Material*			fTankMaterial;
 
    G4OpticalSurface*		fOpTankOutSurface;
    G4LogicalBorderSurface*	fTankOutSurface;

    G4OpticalSurface*		fOpTankInSurface;
    G4LogicalBorderSurface*	fTankInSurface;

    //_____TankUp________________________________________________
    G4double           		fTankUpSizeRmax;
    G4Sphere*			fsolidTankUpSphere;    
    G4SubtractionSolid*		fsolidTankUp;    

    //_____TankDawn______________________________________________
    G4double           		fTankDawnSizeZ;
    G4Polycone*			fsolidTankDawnTube;
    G4double			fTankSteelThickness;
    G4double			fDisalignment;
//_____Tyvek_______________________________________________
	G4double 	fTyvekThickness;
	G4double	fTyvekSphereThickness;

//_____Tank Platform____________________________________________
    G4double           		fPlatformZ;
    G4double 			fPlatformOuterRadius;
    
    G4Tubs*			fsolidExternalPlatformTank;
    G4Tubs*			fsolidInternalPlatformTank;
    G4LogicalVolume*		flogicInternalPlatformTank;
    G4VPhysicalVolume*		fphysicInternalPlatformTank;
    
    G4OpticalSurface*		fOpInternalPlatformSurface;
    G4LogicalBorderSurface*	fInternalPlatformSurface;
    G4Material*			fPlatformInternalMaterial;

//_______SteelPlates on the ground_______________________________

    G4double			fUpPlateSteelR;
    G4double			fUpPlateSteelZ;
    G4Material*			fPlateSteelMaterial;

    G4Tubs*			fsolidUpPlateSteel;
    G4LogicalVolume*		flogicUpPlateSteel;
    G4VPhysicalVolume*		fphysicUpPlateSteel;

//No surfaces with optical properties needed for now

    G4double			fDownPlateSteelR;
    G4double			fDownPlateSteelZ;

    G4Tubs*			fsolidDownPlateSteel;
    G4LogicalVolume*		flogicDownPlateSteel;
    G4VPhysicalVolume*		fphysicDownPlateSteel;

//No surfaces with optical properties needed for now

//_______SSS Legs in OD ___________________________________
   
    G4double fLegCentre;
    G4double fLegExternalRadius;
    G4double fLegThickness;
    G4double fLegInternalRadius;
    G4int fNlegs;

    G4Material*		fLegMaterial;
    G4Tubs*		fsolidLegTubs;
    G4IntersectionSolid*	fsolidLeg;
    G4LogicalVolume*		flogicLeg;
    G4LogicalVolume*		flogicLeg8;
    G4VPhysicalVolume**         fphysicLeg;	
    
    G4OpticalSurface*		fOpLegSurface;
    G4LogicalBorderSurface**	fLegSurface;

//_____AddTyvekSphere (between SSS and water)___________

    G4Sphere*			fsolidAddTyvekSphere;
    G4LogicalVolume*		flogicAddTyvek;
    G4VPhysicalVolume*		fphysicAddTyvek;

    G4OpticalSurface*		fOpAddTyvekOutSurface;
    G4LogicalBorderSurface*   	fAddTyvekOutSurface;

    G4OpticalSurface*		fOpAddTyvekInSurface;
    G4LogicalBorderSurface*   	fAddTyvekInSurface;
  //  G4Material*			fFalseWaterMaterial;


//_______INNER Detector___________________________________
//
//_____SSS________________________________________________

    G4double           		fSSSExternalRadius;
    G4double			fSSSThickness;
    G4Material*			fSSSMaterial;
    G4Sphere*			fsolidSSS;
    G4LogicalVolume*		flogicSSS;
    G4VPhysicalVolume*		fphysicSSS;
    G4OpticalSurface*		fOpSSSSurface;
    G4LogicalBorderSurface*   	fSSSSurface;

    
//_____ZoneIII____________________________________________
    G4double           		fZoneIIIExternalRadius;
    G4Material*			fZoneIIIMaterial;
    G4Sphere*			fsolidZoneIII;
    G4LogicalVolume*		flogicZoneIII;
    G4VPhysicalVolume*		fphysicZoneIII;
    G4OpticalSurface*		fOpZoneIIIInSurface;
    G4LogicalBorderSurface*   	fZoneIIIInSurface;

//_____Shroud___________________________________________

    G4double           		fShroudExternalRadius;
    G4double			fShroudThickness;
    G4Material*			fShroudMaterial;
    G4Sphere*			fsolidShroud;
    G4LogicalVolume*		flogicShroud;
    G4VPhysicalVolume*		fphysicShroud;
    G4OpticalSurface*		fOpShroudInSurface;
    G4LogicalBorderSurface*   	fShroudInSurface;
    G4OpticalSurface*		fOpShroudOutSurface;
    G4LogicalBorderSurface*   	fShroudOutSurface;


//_____ZoneII___________________________________________

    G4double           		fZoneIIExternalRadius;
    G4double			fZoneIIThickness;
    G4Material*			fZoneIIMaterial;
    G4Sphere*			fsolidZoneII;
    G4LogicalVolume*		flogicZoneII;
    G4VPhysicalVolume*		fphysicZoneII;
    G4OpticalSurface*		fOpZoneIIInSurface;
    G4LogicalBorderSurface*   	fZoneIIInSurface;
    G4OpticalSurface*		fOpZoneIIOutSurface;
    G4LogicalBorderSurface*   	fZoneIIOutSurface;

 //_____Vessel___________________________________________

    G4double           		fVesselExternalRadius;
    G4double			fVesselThickness;
    G4Material*			fVesselMaterial;
    G4Sphere*			fsolidVessel;
    G4LogicalVolume*		flogicVessel;
    G4VPhysicalVolume*		fphysicVessel;
    G4OpticalSurface*		fOpVesselInSurface;
    G4LogicalBorderSurface*   	fVesselInSurface;
    G4OpticalSurface*		fOpVesselOutSurface;
    G4LogicalBorderSurface*   	fVesselOutSurface;
    G4ThreeVector               fVesselOrigin ; 
    
//_____Vessel__from__DST_____________________

   G4GenericPolycone*   	        fsolidDVessel; 
   G4GenericPolycone*   	        fsolidDZoneI; 
   
//_____ZoneI_____________________________________________

    G4double           		fZoneIExternalRadius;
    G4Material*			fZoneIMaterial;
    G4Sphere*			fsolidZoneI;
    G4LogicalVolume*		flogicZoneI;
    G4VPhysicalVolume*		fphysicZoneI;
    G4OpticalSurface*		fOpZoneISurface;
    G4LogicalBorderSurface*   	fZoneISurface;


//_____Inner rings to support the ZoneI Vessel construction_____________
    G4bool     IsEndCap;
    G4double   fNylonRingRadiusIn;	
    G4double   fNylonRingRadiusOut;
    G4double   fNylonRingHeight;    
    G4double   fIVNylonRingNorthZ;
    G4double   fIVNylonRingSouthZ;
    G4double   fNylonTubeNorthZ;
    G4double   fNylonTubeSouthZ;
    G4double   fIVTubeNorthZ;
    G4double   fIVTubeSouthZ;
    G4double   fOVNylonRingNorthZ;
    G4double   fOVNylonRingSouthZ;
    G4double   fOVTubeNorthZ;
    G4double   fOVTubeSouthZ;
    G4double   fCopperStrutNorthZ;
    G4double   fCopperStrutSouthZ;
    G4double   fCopperStrutRadialDisplacement;
    G4int      fNumOfStruts;

    G4Material*			fNylonRingMaterial;
    G4Material*			fNylonTubeMaterial;
    G4Material*			fIVTubeMaterial;
    G4Material*			fOVTubeMaterial;
    G4Material*			fCopperStrutMaterial;
   
    //Nylon Bridge
    G4double	fNylonBridgeHeight;
    G4double 	fNylonBridgeRadius;
    G4double 	fNylonBridgeNorthZ;
    G4double	fNylonBridgeSouthZ;
    G4Material*			fNylonBridgeMaterial;

    G4Tubs*			fsolidNylonBridge;
    G4LogicalVolume*	        flogicNylonBridge;
    G4VPhysicalVolume*	        fphysicNylonBridgeNorth;
    G4VPhysicalVolume*	        fphysicNylonBridgeSouth;
    G4OpticalSurface*	   	fOpNylonBridgeNorthInSurface;
    G4LogicalBorderSurface*	fNylonBridgeNorthInSurface;
    G4OpticalSurface*	   	fOpNylonBridgeNorthOutSurface;
    G4LogicalBorderSurface*	fNylonBridgeNorthOutSurface;
    G4OpticalSurface*	   	fOpNylonBridgeSouthInSurface;
    G4LogicalBorderSurface*	fNylonBridgeSouthInSurface;
    G4OpticalSurface*	   	fOpNylonBridgeSouthOutSurface;
    G4LogicalBorderSurface*	fNylonBridgeSouthOutSurface;
    
    G4Tubs*			fsolidIVNylonRing;
    G4Tubs*			fsolidOVNylonRing;
    G4Tubs*			fsolidCopperStrut;
    G4Polycone* 		fsolidNylonTube; 
    G4Polycone* 		fsolidIVTube_raw; 
    G4Polycone* 		fsolidOVTube_raw; 
    G4IntersectionSolid*	fsolidIVTube_North; 
    G4IntersectionSolid*	fsolidOVTube; 
    G4IntersectionSolid*	fsolidIVTube_South; 
    
    G4LogicalVolume*	        flogicCopperStrut;
    G4VPhysicalVolume**	        fphysicCopperStrut;
    
    G4LogicalVolume*	        flogicIVNylonRing;
    G4VPhysicalVolume*	        fphysicIVNylonRingNorth;
    G4VPhysicalVolume*	        fphysicIVNylonRingSouth;
    G4LogicalVolume*	        flogicOVNylonRing;
    G4VPhysicalVolume*	        fphysicOVNylonRingNorth;
    G4VPhysicalVolume*	        fphysicOVNylonRingSouth;
    G4LogicalVolume*	        flogicNylonTube;
    G4VPhysicalVolume*	        fphysicNylonTubeNorth;
    G4VPhysicalVolume*	        fphysicNylonTubeSouth;
    G4LogicalVolume*	        flogicIVTube_North;
    G4LogicalVolume*	        flogicIVTube_South;
    G4VPhysicalVolume*	        fphysicIVTubeNorth;
    G4VPhysicalVolume*	        fphysicIVTubeSouth;
    G4LogicalVolume*	        flogicOVTube;
    G4VPhysicalVolume*	        fphysicOVTubeNorth;
    G4VPhysicalVolume*	        fphysicOVTubeSouth;
    
    G4OpticalSurface*	   	fOpNylonInSurface;
    G4OpticalSurface*	   	fOpNylonOutSurface;
    G4LogicalBorderSurface*	fIVNylonRingNorthInSurface;
    G4LogicalBorderSurface*	fIVNylonRingSouthInSurface;
    G4LogicalBorderSurface*	fIVNylonRingNorthOutSurface;
    G4LogicalBorderSurface*	fIVNylonRingSouthOutSurface;
    G4LogicalBorderSurface*	fOVNylonRingNorthInSurface;
    G4LogicalBorderSurface*	fOVNylonRingSouthInSurface;
    G4LogicalBorderSurface*	fOVNylonRingNorthOutSurface;
    G4LogicalBorderSurface*	fOVNylonRingSouthOutSurface;
    G4LogicalBorderSurface*	fNylonTubeNorthInSurface;
    G4LogicalBorderSurface*	fNylonTubeSouthInSurface;
    G4LogicalBorderSurface*	fNylonTubeNorthOutSurface;
    G4LogicalBorderSurface*	fNylonTubeSouthOutSurface;

    G4OpticalSurface*	   	fOpSteelTubeSurface;
    G4LogicalBorderSurface*	fIVTubeNorthSurface;
    G4LogicalBorderSurface*	fIVTubeSouthSurface;
    G4LogicalBorderSurface*	fOVTubeNorthSurface;
    G4LogicalBorderSurface*	fOVTubeSouthSurface;
    
    G4OpticalSurface*	   	fOpCopperStrutSurface;
    G4LogicalBorderSurface**	fCopperStrutSurface;

    
    //____PMT (Cathode)________________________________________

    G4double           		fPMTSizeDiskR;
    G4double 			fPMTSizeDiskH;
    G4double           		fPMTSizeRmin;
    G4double           		fPMTSizeRmax;
    G4double			fPMTSizeThetaDelta;
    G4Material*			fPMTMaterial;

    G4SubtractionSolid*		fsolidPMT;
    G4Sphere*			fsolidPMTsphere;
    
    G4IntersectionSolid*	fsolidPMTdisk;
    G4LogicalVolume*		flogicPMT;
    G4VPhysicalVolume**		fphysicPMT;
    G4OpticalSurface*		fOpPMTSurface;
    G4LogicalBorderSurface**	fPMTSurface;

//____PMT Shield_________________________________________

    G4Material*			fPMTShieldMaterial;
    G4IntersectionSolid*	fsolidPMTShielddisk;
    G4Polycone*			fsolidPMTShield;

    G4LogicalVolume*		flogicPMTShield;
    G4VPhysicalVolume**		fphysicPMTShield;
    G4OpticalSurface*		fOpPMTShieldSurface;
    G4LogicalBorderSurface**	fPMTShieldSurface;



//____PMT Housing_________________________________________

    G4Polycone*			fsolidPMTbase_cylinder;
    G4IntersectionSolid*	fsolidPMTHousing;
    
    G4LogicalVolume*		flogicPMTHousing;
    G4VPhysicalVolume**		fphysicPMTHousing;
    G4OpticalSurface*		fOpPMTHousingSurface;
    G4LogicalBorderSurface**	fPMTHousingSurface;
    G4Material*			fPMTHousingMaterial;

//____PMT ring for PMTs without concentrators_______________

    G4Polycone*			fsolidPMTring;
    
    G4LogicalVolume*		flogicPMTring;
    G4VPhysicalVolume**		fphysicPMTring;
    G4OpticalSurface*		fOpPMTringSurface;
    G4LogicalBorderSurface**	fPMTringBorderSurface;
    G4Material*			fPMTringMaterial;

//____PMT ring for PMTs without concentrators_______________


//____Light Guides________________________________________

    G4Material*			fExternalLightGuideMaterial;
    G4SubtractionSolid*        	fsolidExternalLightGuide;
    G4LogicalVolume* 		flogicExternalLightGuide;
    G4VPhysicalVolume** 	fphysicExternalLightGuide;
    G4OpticalSurface* 		fOpExternalLGSurface;
    G4LogicalBorderSurface**  	fExternalLGBorderSurface;

    G4Material*			fInternalLightGuideMaterial;
    G4SubtractionSolid*        	fsolidInternalLightGuide;
    G4LogicalVolume* 		flogicInternalLightGuide;
    G4VPhysicalVolume** 	fphysicInternalLightGuide;
    G4OpticalSurface* 		fOpInternalLGSurface;
    G4LogicalBorderSurface**  	fInternalLGBorderSurface;

//____PMT of OUTER DETECTOR _________________________________
//
    G4int     			fNumberOfPMTmu;
    G4int			fNumberOfPMTmuSSS;
    G4int			fNumberOfPMTmuGround;
    G4double			fPMTmuFlangeSizeRmin;
    G4double			fPMTmuFlangeSizeRmax;
    G4double			fPMTmuFlangeSizeZ;
//____PMTmu (Cathode)________________________________________

    G4Material*			fPMTmuMaterial;
    G4Sphere*			fsolidPMTmusphere;
    G4SubtractionSolid*		fsolidPMTmu;
    G4LogicalVolume*		flogicPMTmu;
    G4VPhysicalVolume**		fphysicPMTmu;
    G4OpticalSurface*		fOpPMTmuSurface;
    G4LogicalBorderSurface**	fPMTmuSurface;

//_____PMTmu Shield (only for PMTs on the ground)______________________________________________

    G4Material*			fPMTmuShieldMaterial;
 
    G4Polycone*			fsolidPMTmuShield;
    G4LogicalVolume*		flogicPMTmuShield;
    G4VPhysicalVolume**		fphysicPMTmuShield;
    
    G4OpticalSurface*		fOpPMTmuShieldInSurface;
    //G4OpticalSurface*		fOpPMTmuFaceOutSurface;
    G4LogicalBorderSurface**	fPMTmuShieldInSurface;
    //G4LogicalBorderSurface**	fPMTmuFaceOutSurface;

//____PMTmu Flange_________________________________________
    G4Material*			fPMTmuFlangeMaterial;

    G4Tubs*			fsolidPMTmuFlange;
    G4LogicalVolume*		flogicPMTmuFlange;
    G4VPhysicalVolume**		fphysicPMTmuFlange;

    G4OpticalSurface*		fOpPMTmuFlangeInSurface;
    G4LogicalBorderSurface**	fPMTmuFlangeInSurface;


// BX_USE_OPERA
 /*
    G4bool                      IsOpera;
 
    G4double           		fOperaZoneSizeX;
    G4double           		fOperaZoneSizeY;
    G4double			fOperaZoneSizeZ;
    G4Material*			fOperaZoneMaterial;
    G4Box*			fsolidOperaZone;
    G4LogicalVolume*		flogicOperaZone;
    G4VPhysicalVolume*		fphysicOperaZone;

    G4double           		fOperaSMSizeX;
    G4double           		fOperaSMSizeY;
    G4double			fOperaSMSizeZ;
    G4Material*			fOperaSMMaterial;
    G4Box*			fsolidOperaSM;
    G4LogicalVolume*		flogicOperaSM1;
    G4LogicalVolume*		flogicOperaSM2;
    G4VPhysicalVolume*		fphysicOperaSM1;
    G4VPhysicalVolume*		fphysicOperaSM2;
 
    G4double           		fOperaTTSizeX;
    G4double           		fOperaTTSizeY;
    G4double			fOperaTTSizeZ;
    G4Material*			fOperaTTMaterial;
    G4Box*			fsolidOperaTT;
    G4LogicalVolume*		flogicOperaTT;
    G4VPhysicalVolume*		fphysicOperaTT1;
    G4VPhysicalVolume*		fphysicOperaTT2;

    G4double           		fOperaMSSizeX;
    G4double           		fOperaMSSizeY;
    G4double			fOperaMSSizeZ;
    G4Material*			fOperaMSMaterial;
    G4Box*			fsolidOperaMS;
    G4LogicalVolume*		flogicOperaMS;
    G4VPhysicalVolume*		fphysicOperaMS11;
    G4VPhysicalVolume*		fphysicOperaMS12;
    G4VPhysicalVolume*		fphysicOperaMS21;
    G4VPhysicalVolume*		fphysicOperaMS22;

    G4double                    fOperaTTWallsSpacing;				     
    G4double                    fOperaMSWallsSpacing;
    G4double                    fOperaTTWallsFirstPos;				     
    G4double                    fOperaMSWallsFirstPos;    

    G4ThreeVector               fOperaZonePos;
    G4ThreeVector               fOperaSM1Pos;
    G4ThreeVector               fOperaSM2Pos;
// Opera ends here
*/
    G4int			fDetectorFlag;
    G4int			fPMTFlag;
    G4bool                      fNoConcentrators;

    G4int			fPMTDistributionFlag;
    G4int 			fNumberOfPMT;
    G4int 			fNumberOfVetoPMT;
    G4bool                      IsOverlap;
    G4bool                      noPMTs;
    //G4bool                      IsSensitive;
 
    //Source geometry
    G4int                       fSource ;
    G4double                    fSourceRadius ;
    G4ThreeVector               fSourceOrigin ;
    G4double                    fSourceLong;
    G4double                    fSourceVialThick;
    G4double                    fSourceX;
    G4double                    fSourceY;
    G4double                    fSourceZ;
    G4Material*                 fSourceMaterial;
    G4Material*                 fSourceVialMaterial;
    G4CSGSolid*                 fsolidVialSource; 
    G4CSGSolid*                 fsolidSource; 
    G4LogicalVolume*		flogicSource;
    G4LogicalVolume*		flogicSourceVial;
    G4VPhysicalVolume*		fphysicSource;
    G4VPhysicalVolume*		fphysicSourceVial;
    
    G4OpticalSurface*		fOpSourceVialInSurface;
    G4LogicalBorderSurface*   	fSourceVialInSurface;
    G4OpticalSurface*		fOpSourceVialOutSurface;
    G4LogicalBorderSurface*   	fSourceVialOutSurface;
    G4OpticalSurface*		fOpSourceOutSurface;
    G4LogicalBorderSurface*   	fSourceOutSurface;
    G4bool			IsOverlapSource;
    //surface between ZoneI and source Vial
    G4OpticalSurface*		fOpZoneISourceSurface;
    G4LogicalBorderSurface*   	fZoneISourceSurface;
    
    
//Am-Be Source  Derlin container
    
    G4Material*			fDerlinMaterial;
    G4Polycone*	        	fsolidDerlin;
    G4LogicalVolume*		flogicDerlin;
    G4VPhysicalVolume*		fphysicDerlin;
    
    G4OpticalSurface*		fOpDerlinSurfaceIn;
    G4LogicalBorderSurface*   	fDerlinSurfaceIn;
    G4OpticalSurface*		fOpDerlinSurfaceOut;
    G4LogicalBorderSurface*   	fDerlinSurfaceOut;

//Am-Be Source  screws for derlin container
    
    G4Material*			fScrewMaterial;
    G4UnionSolid* 		fsolidScrew;
    G4LogicalVolume*		flogicScrew;
    G4VPhysicalVolume**		fphysicScrew;
    
    G4OpticalSurface*		fOpScrewSurface;
    G4LogicalBorderSurface**   	fScrewSurface;

//Am-Be Source  Pb container  

    G4Material*			fPbMaterial;
    G4Tubs*			fsolidPb;
    G4LogicalVolume*		flogicPb;
    G4VPhysicalVolume*		fphysicPb;
    
    G4OpticalSurface*		fOpPbSurface;
    G4LogicalBorderSurface*   	fPbSurface;

//Am-Be Source   

    G4Material*			fAmBeMaterial;
    G4Tubs*			fsolidAmBe;
    G4LogicalVolume*		flogicAmBe;
    G4VPhysicalVolume*		fphysicAmBe;
    

};

#endif
