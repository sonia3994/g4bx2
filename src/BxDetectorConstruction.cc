//Created by D. Franco
//Major Revision by A. Caminata and S. Marcocci, Sept. 2014

#include "G4RunManager.hh" 
#include "G4SDManager.hh"
#include "G4String.hh"
#include <fstream>
#include <iostream>
#include "G4SDManager.hh"
//#include "BxSensitiveDetector.hh"
#include "G4Paraboloid.hh"
#include "G4GenericPolycone.hh"
#include "G4CSGSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"
#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include "BxIO.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4RotationMatrix.hh"
#include "G4Point3D.hh"
#include "G4Transform3D.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"

#ifdef BX_USE_OPERA
#include "BxChamberParameterization.hh"
#include "G4PVParameterised.hh"
#endif
#include "BxDetectorConstruction.hh"
#include "BxPropertyCollection.hh"
#include "BxLogger.hh"
#include "BxReadParameters.hh"
#include "BxMaterial.hh"
#include "BxDetectorConstructionMessenger.hh"
//#include "BxLightSource.hh" 
#include "BxOutputVertex.hh"

#include <cmath>
#include <math.h>
#include <time.h>
#include <sstream>

using namespace std;


BxDetectorConstruction::BxDetectorConstruction() {
  IsOverlap        = false;
  if(BxLogger::GetSeverity() == BxLogger::debugging) IsOverlap = true;
  fPMTDistributionFlag=1; //Default is real PMT distribution (0=uniform)
  fNumberOfPMT=1000; //Default only for the uniform PMT distribution case
  fVesselOrigin    = G4ThreeVector(0.,0.,0.);
  fMessenger = new BxDetectorConstructionMessenger(this);
  fPMTFlag=3;
  time_t rawtime;
  time ( &rawtime );
  BxLog(routine) << ctime(&rawtime) << endlog ;
  IsTank=true;
  noPMTs=false;
  fDetectorFlag=0;
  fNoConcentrators = false ;	
  IsEndCap           = true ;
  /*
  IsOpera          = false ;
  IsSensitive      = false ;
  IsRock           = false ;
  */
  fSource          = 0;
  fSourceRadius    = 3.*cm;
  fSourceOrigin    = G4ThreeVector(0.,0.,0.);
  fSourceVialThick = 1*mm ;
  fSourceLong      = 3.*cm ;
  fSourceX         = 3.*cm ;
  fSourceY         = 3.*cm ;
  fSourceZ         = 3.*cm ;
  IsOverlapSource  = false;
  fSourceMaterial      = BxMaterial::Get()->GetWater();
  fSourceVialMaterial  = BxMaterial::Get()->GetQuartz();
}

BxDetectorConstruction::~BxDetectorConstruction(){
  delete fMessenger;
}

G4VPhysicalVolume* BxDetectorConstruction::Construct() {
  DefineSizes();
  DefineMaterials(); 
 // new BxLightSource;
  return ConstructDetector();
}

void BxDetectorConstruction::DefineSizes() {

  fNumberOfVetoPMT =	 	    BxReadParameters::Get()->GetNumberOfVetoPMT();	      

//World
  fWorldSizeX =	     	            BxReadParameters::Get()->GetWorldSizeX();		
  fWorldSizeY =	     	            BxReadParameters::Get()->GetWorldSizeY();				      
  fWorldSizeZ =	     	            BxReadParameters::Get()->GetWorldSizeZ();				      

//TankUp							      
  fTankUpSizeRmax =	            BxReadParameters::Get()->GetTankUpSizeRmax();			      

//TankDawn							      
  fTankDawnSizeZ =	            BxReadParameters::Get()->GetTankDawnSizeZ();			      
  fTankSteelThickness=			BxReadParameters::Get()->GetTankSteelThickness();
//TankUpSphere and SSS centres disalignment
  fDisalignment=		BxReadParameters::Get()->GetDisalignment();

//Plate Steel on the ground
  fUpPlateSteelR=				BxReadParameters::Get()->GetUpPlateSteelR();
  fUpPlateSteelZ=				BxReadParameters::Get()->GetUpPlateSteelZ();
  fDownPlateSteelR=				BxReadParameters::Get()->GetDownPlateSteelR();
  fDownPlateSteelZ=				BxReadParameters::Get()->GetDownPlateSteelZ();

//Outer Detector Platform and SSS Legs 							      
  fPlatformZ =	     	            BxReadParameters::Get()->GetPlatformZ();				      
  fPlatformOuterRadius= 	    BxReadParameters::Get()->GetPlatformOuterRadius();
  fLegCentre=			BxReadParameters::Get()->GetLegCentre();
  fLegExternalRadius=		BxReadParameters::Get()->GetLegExternalRadius();
  fLegThickness=		BxReadParameters::Get()->GetLegThickness();
  fLegInternalRadius= 		fLegExternalRadius-fLegThickness;
  fNlegs=			BxReadParameters::Get()->GetNlegs();
//Tyvek
  fTyvekThickness=		BxReadParameters::Get()->GetTyvekThickness();
  fTyvekSphereThickness=	BxReadParameters::Get()->GetTyvekSphereThickness();

//SSS								      
  fSSSExternalRadius =  BxReadParameters::Get()->GetSSSExternalRadius();			      
  fSSSThickness =  BxReadParameters::Get()->GetSSSThickness();  			      

  //ZoneIII							      
  fZoneIIIExternalRadius =          BxReadParameters::Get()->GetZoneIIIExternalRadius();		      

  
//Shroud							      
  fShroudExternalRadius =           BxReadParameters::Get()->GetShroudExternalRadius(); 		      
  fShroudThickness =	            BxReadParameters::Get()->GetShroudThickness();			      

//ZoneII							      
  fZoneIIExternalRadius =           BxReadParameters::Get()->GetZoneIIExternalRadius(); 		      
  fZoneIIThickness =	            BxReadParameters::Get()->GetZoneIIThickness();			      

//Vessel						      
  fVesselExternalRadius =           BxReadParameters::Get()->GetVesselExternalRadius(); 		      
  if(fDetectorFlag == 5)      fVesselExternalRadius = 1.*m ;
  fVesselThickness =	            BxReadParameters::Get()->GetVesselThickness();			      
  fVesselOrigin = 	            BxReadParameters::Get()->GetVesselOrigin();	
//ZoneI 							      
  fZoneIExternalRadius =            BxReadParameters::Get()->GetZoneIExternalRadius();  	      

//Inner ring to support the ZoneI Vessel construction

  fNylonRingRadiusOut =           BxReadParameters::Get()->GetNylonRingRadiusOut();
  fNylonRingRadiusIn =            BxReadParameters::Get()->GetNylonRingRadiusIn();	
  fNylonRingHeight =  		BxReadParameters::Get()->GetNylonRingHeight();    
  
//Nylon bridge

  fNylonBridgeHeight =		BxReadParameters::Get()->GetNylonBridgeHeight();
  fNylonBridgeRadius =		BxReadParameters::Get()->GetNylonBridgeRadius();

//PMT (Cathode) 						      
  fPMTSizeRmin =	            BxReadParameters::Get()->GetPMTSizeRmin();  			      
  fPMTSizeRmax =	            BxReadParameters::Get()->GetPMTSizeRmax();  			      
  fPMTSizeThetaDelta =	            BxReadParameters::Get()->GetPMTSizeThetaDelta();			      
 
//PMT simple disk

  fPMTSizeDiskR =                  BxReadParameters::Get()->GetPMTSizeDiskR();             
  fPMTSizeDiskH=		   BxReadParameters::Get()->GetPMTSizeDiskH();  

//OD PMTs
  fPMTmuFlangeSizeRmin=			BxReadParameters::Get()->GetPMTmuFlangeSizeRmin();
  fPMTmuFlangeSizeRmax=			BxReadParameters::Get()->GetPMTmuFlangeSizeRmax();
  fPMTmuFlangeSizeZ=			BxReadParameters::Get()->GetPMTmuFlangeSizeZ();
}

void BxDetectorConstruction::DefineMaterials() {
  
  BxLog(trace) << "Detector Flag = " << fDetectorFlag << endlog;
 
  BxLog(trace) << "PMT Flag = " << fPMTFlag << endlog;
  
  fWorldMaterial 		= BxMaterial::Get()->GetAir();
 
 fTankSteelMaterial             = BxMaterial::Get()->GetSteel(); 
 fTyvekMaterial			= BxMaterial::Get()->GetTyvek(); 
 fPlatformInternalMaterial      = BxMaterial::Get()->GetSteel();     
 fLegMaterial			= BxMaterial::Get()->GetSteel();
 fPlateSteelMaterial		= BxMaterial::Get()->GetSteel();
 fPbMaterial                    = BxMaterial::Get()->GetPb();    
 fDerlinMaterial                = BxMaterial::Get()->GetDerlin(); 
 fScrewMaterial                 = BxMaterial::Get()->GetSteel(); 
 fAmBeMaterial                  = BxMaterial::Get()->GetAmBe();   

  if (fDetectorFlag == 0) {     //Full Model Scintillator (OUTER + Inner Detector)
    BxLog(routine) << "Full detector with scintillator" << endlog;
    fTankMaterial               = BxMaterial::Get()->GetWater();     
    fTankSteelMaterial          = BxMaterial::Get()->GetSteel();       
    fTyvekMaterial		= BxMaterial::Get()->GetTyvek(); 
    fPlatformInternalMaterial   = BxMaterial::Get()->GetSteel();      
    fLegMaterial		= BxMaterial::Get()->GetSteel();
    fPlateSteelMaterial		= BxMaterial::Get()->GetSteel();
    fPMTmuMaterial 		= BxMaterial::Get()->GetBialkali();
    fPMTmuFlangeMaterial   	= BxMaterial::Get()->GetSteel();  
    fPMTmuShieldMaterial   	= BxMaterial::Get()->GetSteel();  
    fSSSMaterial  		= BxMaterial::Get()->GetSteel();
    fZoneIIIMaterial 		= BxMaterial::Get()->GetDMPbuffer();
    fShroudMaterial 	      	= BxMaterial::Get()->GetNylon();
    fZoneIIMaterial 		= BxMaterial::Get()->GetDMPbuffer();
    fVesselMaterial 	      	= BxMaterial::Get()->GetNylon();
    fZoneIMaterial 	      	= BxMaterial::Get()->GetScintillator();
    fExternalLightGuideMaterial = BxMaterial::Get()->GetAluminum();    	   
    fInternalLightGuideMaterial = BxMaterial::Get()->GetAluminum();    	   
    fPMTMaterial 	        = BxMaterial::Get()->GetBialkali();        
    fPMTShieldMaterial		= BxMaterial::Get()->GetMuMetal();
    fPMTHousingMaterial		= BxMaterial::Get()->GetSteel();
    fNylonTubeMaterial          = BxMaterial::Get()->GetNylon();
    fNylonRingMaterial          = BxMaterial::Get()->GetNylon();
    fNylonBridgeMaterial        = BxMaterial::Get()->GetNylon();
    fCopperStrutMaterial        = BxMaterial::Get()->GetCopper();
    fIVTubeMaterial             = BxMaterial::Get()->GetSteel();
    fOVTubeMaterial             = BxMaterial::Get()->GetSteel();
    fPMTringMaterial		= BxMaterial::Get()->GetSteel();

  } else if (fDetectorFlag == 1) {     //Full Model Water (OUTER + Inner Detector)
	BxLog(routine) << "Full detector with water" << endlog;
    fTankMaterial               = BxMaterial::Get()->GetWater();    
    fTankSteelMaterial          = BxMaterial::Get()->GetSteel();       
    fTyvekMaterial		= BxMaterial::Get()->GetTyvek(); 
    fPlatformInternalMaterial   = BxMaterial::Get()->GetSteel();      
    fLegMaterial		= BxMaterial::Get()->GetSteel();
    fPlateSteelMaterial		= BxMaterial::Get()->GetSteel();
    fPMTmuMaterial 		= BxMaterial::Get()->GetBialkali();
    fPMTmuFlangeMaterial   	= BxMaterial::Get()->GetSteel();  
    fPMTmuShieldMaterial   	= BxMaterial::Get()->GetSteel();  
    fSSSMaterial  		= BxMaterial::Get()->GetSteel();
    fZoneIIIMaterial 		= BxMaterial::Get()->GetWater();
    fShroudMaterial 	      	= BxMaterial::Get()->GetNylon();
    fZoneIIMaterial 		= BxMaterial::Get()->GetWater();
    fVesselMaterial 	      	= BxMaterial::Get()->GetNylon();
    fZoneIMaterial 	      	= BxMaterial::Get()->GetWater();
    fExternalLightGuideMaterial = BxMaterial::Get()->GetAluminum();  	  
    fInternalLightGuideMaterial = BxMaterial::Get()->GetAluminum();  	  
    fPMTMaterial 	        = BxMaterial::Get()->GetBialkali();	  
    fPMTShieldMaterial		= BxMaterial::Get()->GetMuMetal();
    fPMTHousingMaterial		= BxMaterial::Get()->GetSteel();
    fNylonTubeMaterial          = BxMaterial::Get()->GetNylon();
    fNylonRingMaterial          = BxMaterial::Get()->GetNylon();
    fNylonBridgeMaterial        = BxMaterial::Get()->GetNylon();
    fCopperStrutMaterial        = BxMaterial::Get()->GetCopper();
    fIVTubeMaterial             = BxMaterial::Get()->GetSteel();
    fOVTubeMaterial             = BxMaterial::Get()->GetSteel();
    fPMTringMaterial		= BxMaterial::Get()->GetSteel();

  } else if (fDetectorFlag == 2) { //Only  Inner Detector
    BxLog(routine) << "Inner detector with water" << endlog;    
    fTankMaterial               = BxMaterial::Get()->GetWater();	  
    fTankSteelMaterial          = BxMaterial::Get()->GetSteel();	  
    fTyvekMaterial		= BxMaterial::Get()->GetTyvek(); 
    fPlatformInternalMaterial   = BxMaterial::Get()->GetSteel();      
    fLegMaterial		= BxMaterial::Get()->GetSteel();
    fPlateSteelMaterial		= BxMaterial::Get()->GetSteel();
    fPMTmuMaterial 		= BxMaterial::Get()->GetWater();
    fPMTmuShieldMaterial 	= BxMaterial::Get()->GetWater();
    fPMTmuFlangeMaterial   	= BxMaterial::Get()->GetWater();  
    fSSSMaterial  		= BxMaterial::Get()->GetSteel();
    fZoneIIIMaterial 		= BxMaterial::Get()->GetWater();
    fShroudMaterial 	      	= BxMaterial::Get()->GetNylon();
    fZoneIIMaterial 		= BxMaterial::Get()->GetWater();
    fVesselMaterial 	      	= BxMaterial::Get()->GetNylon();
    fZoneIMaterial 	      	= BxMaterial::Get()->GetWater();
    fExternalLightGuideMaterial = BxMaterial::Get()->GetAluminum();    	   
    fInternalLightGuideMaterial = BxMaterial::Get()->GetAluminum();    	   
    fPMTMaterial 	        = BxMaterial::Get()->GetBialkali();        
    fPMTShieldMaterial		= BxMaterial::Get()->GetMuMetal();
    fPMTHousingMaterial		= BxMaterial::Get()->GetSteel();
    fNylonTubeMaterial          = BxMaterial::Get()->GetNylon();
    fNylonRingMaterial          = BxMaterial::Get()->GetNylon();
    fNylonBridgeMaterial        = BxMaterial::Get()->GetNylon();
    fCopperStrutMaterial        = BxMaterial::Get()->GetCopper();
    fIVTubeMaterial             = BxMaterial::Get()->GetSteel();
    fOVTubeMaterial             = BxMaterial::Get()->GetSteel();
    fPMTringMaterial		= BxMaterial::Get()->GetSteel();

  } else if (fDetectorFlag == 3){ //Only Outer Detector
    BxLog(routine) << "Outer detector" << endlog;    
    fTankMaterial               = BxMaterial::Get()->GetWater();    
    fTankSteelMaterial          = BxMaterial::Get()->GetSteel();       
    fTyvekMaterial		= BxMaterial::Get()->GetTyvek(); 
    fPlatformInternalMaterial   = BxMaterial::Get()->GetSteel();      
    fLegMaterial		= BxMaterial::Get()->GetSteel();
    fPlateSteelMaterial		= BxMaterial::Get()->GetSteel();
    fPMTmuMaterial 		= BxMaterial::Get()->GetBialkali();
    fPMTmuFlangeMaterial   	= BxMaterial::Get()->GetSteel();
    fPMTmuShieldMaterial   	= BxMaterial::Get()->GetSteel();
    fSSSMaterial  		= BxMaterial::Get()->GetSteel();
    fZoneIIIMaterial 		= BxMaterial::Get()->GetWater();
    fShroudMaterial 	      	= BxMaterial::Get()->GetWater();
    fZoneIIMaterial 		= BxMaterial::Get()->GetWater();
    fVesselMaterial 	      	= BxMaterial::Get()->GetWater();
    fZoneIMaterial 	      	= BxMaterial::Get()->GetWater();
    fExternalLightGuideMaterial = BxMaterial::Get()->GetWater();	
    fInternalLightGuideMaterial = BxMaterial::Get()->GetWater();	
    fPMTMaterial 	        = BxMaterial::Get()->GetWater();	
    fPMTShieldMaterial		= BxMaterial::Get()->GetWater();
    fPMTHousingMaterial		= BxMaterial::Get()->GetWater();
    fNylonTubeMaterial          = BxMaterial::Get()->GetWater();
    fNylonRingMaterial          = BxMaterial::Get()->GetWater();
    fNylonBridgeMaterial        = BxMaterial::Get()->GetWater();
    fCopperStrutMaterial        = BxMaterial::Get()->GetWater();
    fIVTubeMaterial             = BxMaterial::Get()->GetWater();
    fOVTubeMaterial             = BxMaterial::Get()->GetWater();
    fPMTringMaterial		= BxMaterial::Get()->GetWater();
    
  } else if (fDetectorFlag == 4) { //Only  Inner Detector with Scintillator
    BxLog(routine) << "Inner detector with scintillator" << endlog;    
    fTankMaterial               = BxMaterial::Get()->GetWater();	  
    fTankSteelMaterial          = BxMaterial::Get()->GetSteel();	  
    fTyvekMaterial		= BxMaterial::Get()->GetTyvek(); 
    fPlatformInternalMaterial   = BxMaterial::Get()->GetSteel();      
    fLegMaterial		= BxMaterial::Get()->GetSteel();
    fPlateSteelMaterial		= BxMaterial::Get()->GetSteel();
    fPMTmuMaterial 		= BxMaterial::Get()->GetWater();
    fPMTmuFlangeMaterial   	= BxMaterial::Get()->GetWater();  
    fPMTmuShieldMaterial   	= BxMaterial::Get()->GetWater();  
    fSSSMaterial  		= BxMaterial::Get()->GetSteel();
    fZoneIIIMaterial 		= BxMaterial::Get()->GetDMPbuffer();
    fShroudMaterial 	      	= BxMaterial::Get()->GetNylon();
    fZoneIIMaterial 		= BxMaterial::Get()->GetDMPbuffer();
    fVesselMaterial 	      	= BxMaterial::Get()->GetNylon();
    fZoneIMaterial 	      	= BxMaterial::Get()->GetScintillator();
    fExternalLightGuideMaterial = BxMaterial::Get()->GetAluminum();    	   
    fInternalLightGuideMaterial = BxMaterial::Get()->GetAluminum();    	   
    fPMTMaterial 	        = BxMaterial::Get()->GetBialkali();        
    fPMTShieldMaterial		= BxMaterial::Get()->GetMuMetal();
    fPMTHousingMaterial		= BxMaterial::Get()->GetSteel();
    fNylonTubeMaterial          = BxMaterial::Get()->GetNylon();
    fNylonRingMaterial          = BxMaterial::Get()->GetNylon();
    fNylonBridgeMaterial        = BxMaterial::Get()->GetNylon();
    fCopperStrutMaterial        = BxMaterial::Get()->GetCopper();
    fIVTubeMaterial             = BxMaterial::Get()->GetSteel();
    fOVTubeMaterial             = BxMaterial::Get()->GetSteel();
    fPMTringMaterial		= BxMaterial::Get()->GetSteel();
    
  } else if (fDetectorFlag == 5) { //CTF like configuration. 
    BxLog(routine) << "CTF like configuration" << endlog;    
    fTankMaterial               = BxMaterial::Get()->GetWater();	  
    fTankSteelMaterial          = BxMaterial::Get()->GetSteel();	  
    fTyvekMaterial		= BxMaterial::Get()->GetTyvek(); 
    fPlatformInternalMaterial   = BxMaterial::Get()->GetSteel();      
    fLegMaterial		= BxMaterial::Get()->GetSteel();
    fPlateSteelMaterial		= BxMaterial::Get()->GetSteel();
    fPMTmuMaterial 		= BxMaterial::Get()->GetWater();
    fPMTmuFlangeMaterial   	= BxMaterial::Get()->GetWater();  
    fPMTmuShieldMaterial   	= BxMaterial::Get()->GetWater();  
    fSSSMaterial  		= BxMaterial::Get()->GetSteel();
    fZoneIIIMaterial 		= BxMaterial::Get()->GetWater();
    fShroudMaterial 	      	= BxMaterial::Get()->GetNylon();
    fZoneIIMaterial 		= BxMaterial::Get()->GetWater();
    fVesselMaterial 	      	= BxMaterial::Get()->GetNylon();
    fZoneIMaterial 	      	= BxMaterial::Get()->GetScintillator();
    fExternalLightGuideMaterial = BxMaterial::Get()->GetAluminum();    	   
    fInternalLightGuideMaterial = BxMaterial::Get()->GetAluminum();    	   
    fPMTMaterial 	        = BxMaterial::Get()->GetBialkali();        
    fPMTShieldMaterial		= BxMaterial::Get()->GetMuMetal();
    fPMTHousingMaterial		= BxMaterial::Get()->GetSteel();
    fNylonTubeMaterial          = BxMaterial::Get()->GetNylon();
    fNylonRingMaterial          = BxMaterial::Get()->GetNylon();
    fNylonBridgeMaterial          = BxMaterial::Get()->GetNylon();
    fCopperStrutMaterial        = BxMaterial::Get()->GetCopper();
    fIVTubeMaterial             = BxMaterial::Get()->GetSteel();
    fOVTubeMaterial             = BxMaterial::Get()->GetSteel();
    fPMTringMaterial		= BxMaterial::Get()->GetSteel();
    
	}
}

 
G4VPhysicalVolume* BxDetectorConstruction::ConstructDetector() {

  BxOutputVertex::Get()->SetDetectorFlag(fDetectorFlag);
  BxOutputVertex::Get()->SetPMTFlag(fPMTFlag);    
  
  std::vector<G4ThreeVector>   fPMTPosition;
if (fPMTDistributionFlag==1)
fPMTPosition=BxReadParameters::Get()->GetPMTPosition();

  std::vector<G4int>&		fVetoPmt             =	   BxReadParameters::Get()->GetVetoPmt();	  
  std::vector<G4ThreeVector>&   fPMT_VetoOD_Position =	   BxReadParameters::Get()->GetPMT_VetoOD_Position();  
  std::vector<G4int>&		fVetoOD_Pmt	     =     BxReadParameters::Get()->GetVetoOD_Pmt();  

    // World 
  fsolidWorld  = new G4Box("World_solid",fWorldSizeX/2.0,fWorldSizeY/2.0,fWorldSizeZ/2.0);
 
  flogicWorld  = new G4LogicalVolume(fsolidWorld, fWorldMaterial, "World_logic");
  fphysicWorld = new G4PVPlacement(0,
  				 G4ThreeVector(0,0,0),
                                 "World",
                                 flogicWorld,
                                 NULL,
                                 false,
                                 0);
   

/*	
if(IsRock) {
    BxLog(routine) << "Rock and concrete builded " << endlog ;
    G4double rockX = 8.*m;
    G4double rockY = 8.*m;
    G4double rockZ = 1.5*m;
    G4double concX = 8.*m;
    G4double concY = 8.*m;
    G4double concZ = 0.15*m;

    BxLog(routine) << "Attention: put the cosmic ray generator above " 
                   << (10.*m+2*rockZ+2*concZ)/m << " m to cross rock and concrete!!!"  << endlog ;

    fsolidRock  = new G4Box("Rock_solid",rockX,rockY,rockZ);

    flogicRock  = new G4LogicalVolume(fsolidRock, BxMaterial::Get()->GetRock(), "Rock_logic");

    fphysicRock = new G4PVPlacement(0,
  				   G4ThreeVector(0,0,10.*m+rockZ+2*concZ),
                                   "Rock",
                                   flogicRock,
                                   fphysicWorld,
                                   false,
                                   0,
			           IsOverlap);

    fsolidConcrete  = new G4Box("Concrete_solid",concX,concY,concZ);

    flogicConcrete  = new G4LogicalVolume(fsolidConcrete, BxMaterial::Get()->GetConcrete(), "Concrete_logic");

    fphysicConcrete = new G4PVPlacement(0,
  				   G4ThreeVector(0,0,10.*m+concZ),
                                   "Concrete",
                                   flogicConcrete,
                                   fphysicWorld,
                                   false,
                                   0,
			           IsOverlap);
  }*/
  
  if(IsTank) {

     //2 Steel plates on the ground at the bottom of the detector 
    
     fsolidUpPlateSteel = new G4Tubs("UpPlateSteel_solid",
				0.0, fUpPlateSteelR,
				fUpPlateSteelZ/2.0,
				0.0, 2*pi);


    flogicUpPlateSteel = new G4LogicalVolume(fsolidUpPlateSteel,
                         fPlateSteelMaterial,
                         "UpPlateSteel_logic");

    fphysicUpPlateSteel = new G4PVPlacement(0,
  			 G4ThreeVector(0,0,-fUpPlateSteelZ/2.0-fDisalignment-fTankDawnSizeZ-fTankSteelThickness-1*mm), 
			 //correction for the non alignment of the SSS and TankUpSteel centres + Tank height + 1mm tolerance
                	 "UpPlateSteel",
                	 flogicUpPlateSteel,
                	 fphysicWorld,
                	 false,
                	 0,
			 IsOverlap);


     fsolidDownPlateSteel = new G4Tubs("DownPlateSteel_solid",
				0.0, fDownPlateSteelR,
				fDownPlateSteelZ/2.0,
				0.0, 2*pi);


    flogicDownPlateSteel = new G4LogicalVolume(fsolidDownPlateSteel,
                         fPlateSteelMaterial,
                         "DownPlateSteel_logic");

    fphysicDownPlateSteel = new G4PVPlacement(0,
  			 G4ThreeVector(0,0,-fDownPlateSteelZ/2.0-fUpPlateSteelZ-fDisalignment-fTankDawnSizeZ-fTankSteelThickness-2*mm), 
			 //correction for the non alignment of the SSS and TankUpSteel centres + Tank height +PlateUp height+ 1mm tolerance
                	 "DownPlateSteel",
                	 flogicDownPlateSteel,
                	 fphysicWorld,
                	 false,
                	 0,
			 IsOverlap);


  //Box to create a portion of the upper sphere on the tank 
  G4Box* SubtractionBox=new G4Box("SubtractionBox", 5*m, 5*m, 5*m);
 
  //G4ThreeVector for the translation of the SubtractionBox
  G4ThreeVector zTransBox(0.,0.,-6.*m);

  //Null rotation for boolean operations
  G4RotationMatrix *yRot = new G4RotationMatrix;
  yRot->rotateY(0.);
  
  
  //Bxon Steel SphereTank ( Upper Part) 
    fsolidTankUpSphereSteel = new G4Sphere("TankUpSphereSteel_solid",
			  0.0, fTankUpSizeRmax+fTankSteelThickness,
			  0.0, 2*pi,
			  0.0, pi);


    fsolidTankUpSteel=new G4SubtractionSolid("TankUpSteel_solid",fsolidTankUpSphereSteel,SubtractionBox,yRot,zTransBox);
    
    //Bxon Steel  TubeTank ( Dawn Part of OUTER Detector) 
    fsolidTankTubeSteel = new G4Tubs("TankTubeSteel_solid",
				0.0, fTankUpSizeRmax+fTankSteelThickness,
				(fTankDawnSizeZ+fTankSteelThickness)/2.0,
				0.0, 2*pi);


    G4ThreeVector zTrans(0. ,0. , -(fTankDawnSizeZ+fTankSteelThickness)/2.0);//6mm of steel on the floor


    G4UnionSolid *fsolidTankSteel = new G4UnionSolid("TankSteel_solid", fsolidTankUpSteel, fsolidTankTubeSteel, yRot, zTrans); 


    flogicTankSteel = new G4LogicalVolume(fsolidTankSteel,
                         fTankSteelMaterial,
                         "TankSteel_logic");

    fphysicTankSteel = new G4PVPlacement(0,
  			 G4ThreeVector(0,0,-fDisalignment), //correction for the non alignment of the SSS and TankUpSteel centres
                	 "TankSteel",
                	 flogicTankSteel,
                	 fphysicWorld,
                	 false,
                	 0,
			 IsOverlap);



     //Bxon Tyvek  SphereTank ( Upper Part) 
    fsolidTankUpSphereTyvek = new G4Sphere("TankUpSphereTyvek_solid",
			  0.0, fTankUpSizeRmax, 
			  0.0, 2*pi,
			  0.0, pi);
	

    fsolidTankUpTyvek=new G4SubtractionSolid("TankUpTyvek_solid",fsolidTankUpSphereTyvek,SubtractionBox,yRot,zTransBox);

    fsolidTankTubeTyvek=new G4Polycone("TankTubeTyvek_solid",
			 			0.0, twopi,
						BxReadParameters::Get()->GetNumOfZplanesTankTyvek(),
						BxReadParameters::Get()->GetZplanesTankTyvek(),
						BxReadParameters::Get()->GetRinTankTyvek(),
						BxReadParameters::Get()->GetRoutTankTyvek());

    G4ThreeVector zTrans1(0. ,0. , -fTankDawnSizeZ);


    G4UnionSolid *fsolidTankTyvek = new G4UnionSolid("TankTyvek_solid", fsolidTankUpTyvek, fsolidTankTubeTyvek, yRot, zTrans1); 


    flogicTankTyvek = new G4LogicalVolume(fsolidTankTyvek,
                         fTyvekMaterial,
                         "TankTyvek_logic");

    fphysicTankTyvek = new G4PVPlacement(0,
  			 G4ThreeVector(0,0,0),
                	 "TankTyvek",
                	 flogicTankTyvek,
			 fphysicTankSteel,
			 false,
			 0,
			 IsOverlap);

    
    //Bxon Veto Water Tank (Upper Part)
    fsolidTankUpSphere  = new G4Sphere("TankUpSphere_solid",
			   0, fTankUpSizeRmax-fTyvekThickness,
			   0.0, 2*pi,
			   0.0, pi);

    fsolidTankUp=new G4SubtractionSolid("TankUp_solid",fsolidTankUpSphere,SubtractionBox,yRot,zTransBox);
    
    fsolidTankDawnTube=new G4Polycone("TankDawnTube_solid",
			 			0.0, twopi,
						BxReadParameters::Get()->GetNumOfZplanesTank(),
						BxReadParameters::Get()->GetZplanesTank(),
						BxReadParameters::Get()->GetRinTank(),
						BxReadParameters::Get()->GetRoutTank());
  
	G4ThreeVector zTrans2(0. ,0. , -fTankDawnSizeZ+fTyvekThickness);

    G4UnionSolid *fsolidTank = new G4UnionSolid("Tank_solid", fsolidTankUp, fsolidTankDawnTube, yRot, zTrans2); 


    flogicTank  = new G4LogicalVolume(fsolidTank,
                          fTankMaterial,
                          "Tank_logic");

    fphysicTank = new G4PVPlacement(0,
                          G4ThreeVector(0,0,0),
                          "Tank",
                          flogicTank,
                          fphysicTankTyvek,
                          false,
                          0,
			  IsOverlap);
      
      //AddTyvekSphere (For the Tyvek on the SSS )+ the External Platform to divide upper and lower volumes of the tank (made of Tyvek)
    fsolidAddTyvekSphere = new G4Sphere("AddTyvekSphere_solid",
			   0.0, fSSSExternalRadius+fTyvekSphereThickness,
			   0.0, twopi,
			   0.0, pi);

    fsolidExternalPlatformTank = new G4Tubs("ExternalPlatformTank_solid",
                                 fSSSExternalRadius+fTyvekSphereThickness/2., fPlatformOuterRadius+fTyvekThickness,
 	                         fPlatformZ,
	                         0.0, 2*pi);

    
    G4ThreeVector zTrans3(0. ,0. , -fPlatformZ/2.0+300*mm);

    G4UnionSolid *fsolidAddTyvek = new G4UnionSolid("AddTyvek_solid", fsolidAddTyvekSphere, fsolidExternalPlatformTank, yRot, zTrans3); 
      
    flogicAddTyvek = new G4LogicalVolume(fsolidAddTyvek,
                           fTyvekMaterial,
                           "AddTyvek_logic");

    fphysicAddTyvek = new G4PVPlacement(0,	
  			   G4ThreeVector(0,0,fDisalignment),
                           "AddTyvek",
                           flogicAddTyvek,
                           fphysicTank,
                           false,
                           0, IsOverlap);

    //Inside the External Platform, I insert a steel Internal Platform to take into account the muons energy loss properly in the detector

      fsolidInternalPlatformTank = new G4Tubs("InternalPlatformTank_solid",
                                 fSSSExternalRadius+fTyvekSphereThickness+fTyvekThickness, fPlatformOuterRadius,
 	                         fPlatformZ-fTyvekThickness,
	                         0.0, 2*pi);
      
      flogicInternalPlatformTank = new G4LogicalVolume(fsolidInternalPlatformTank,
	                         fPlatformInternalMaterial,
 	                         "PlatformInternal_logic");
  
      fphysicInternalPlatformTank = new G4PVPlacement(0,
 	                         G4ThreeVector(0,0,-(fPlatformZ-fTyvekThickness)/2+300*mm),
 	                         "PlatformInternal",
 	                         flogicInternalPlatformTank,
 	                         fphysicAddTyvek,
 	                         false,
 	                         0, IsOverlap);
      
   //  SSS Legs 
   
      fsolidLegTubs=new G4Tubs("LegTubs_solid",
		               fLegInternalRadius, fLegExternalRadius,
			       fTankDawnSizeZ/2.,
			       0.0, twopi);
      G4ThreeVector zTrans4(0.0, fLegCentre , -fTankDawnSizeZ/2.);

      G4SubtractionSolid* sub=new G4SubtractionSolid("Subtraction", fsolidLegTubs,fsolidAddTyvek,yRot,-zTrans4);

      fsolidLeg=new G4IntersectionSolid("LegSolid",fsolidTank,sub,yRot,zTrans4);
      

      flogicLeg = new G4LogicalVolume(fsolidLeg,
	                         fLegMaterial,
 	                         "Leg_logic");

      fphysicLeg=new G4VPhysicalVolume*[fNlegs];
     
      //Volumes for Leg #8 which has an overlap with an OD PMT on the sphere
      G4double Leg8ExternalRadius=5*cm;
      G4Tubs* fsolidLegTub8=new G4Tubs("LegTub8_solid",
		               Leg8ExternalRadius-fLegThickness, Leg8ExternalRadius,
			       fTankDawnSizeZ/2.,
			       0.0, twopi);

      G4SubtractionSolid* sub2=new G4SubtractionSolid("Subtraction2", fsolidLegTub8,fsolidAddTyvek,yRot,-zTrans4);

      G4IntersectionSolid* fsolidLeg8=new G4IntersectionSolid("Leg8Solid",fsolidTank,sub2,yRot,zTrans4);
      
      flogicLeg8 = new G4LogicalVolume(fsolidLeg8,
	                         fLegMaterial,
 	                         "Leg8_logic");
      
      for (G4int i=0; i<fNlegs; i++){
      G4String LegName="Leg ";
      stringstream  str;
      str << i;
      LegName+=str.str();
      
      G4Transform3D tr;
      tr=G4RotateZ3D(twopi*i/fNlegs);
      
      if (i!=8){
      fphysicLeg[i] = new G4PVPlacement(tr,
				    flogicLeg,
				    LegName,
                                   flogicTank,
	                            false,
                                    i, IsOverlap);
      } else {
      fphysicLeg[i] = new G4PVPlacement(tr,
				    flogicLeg8,
				    LegName,
                                   flogicTank,
	                            false,
                                    i, IsOverlap);
      }
      }
      
      
      
      BxLog(routine) << "Tank built" << endlog ;
  }



  


//SSS
   
    fsolidSSS = new G4Sphere("SSS_solid",0, fSSSExternalRadius,
		    		0.0,twopi,
                                0.0, pi);

   
    flogicSSS=new G4LogicalVolume(fsolidSSS,fSSSMaterial,   "SSS_logic");
   

  if(IsTank) {
    fphysicSSS = new G4PVPlacement(0,
  				   G4ThreeVector(0,0,0),
                                   "SSS",
                                   flogicSSS,
                                   fphysicAddTyvek,
				   false,
                                   0, IsOverlap);
  } else {
   fphysicSSS = new G4PVPlacement(0, G4ThreeVector(0,0,0),flogicSSS, "SSS",flogicWorld,false,0,IsOverlap);
   
  }
 
   //ZoneIII
      fsolidZoneIII= new G4Sphere("ZoneIII_solid", 0,fZoneIIIExternalRadius,0,twopi,0,pi);

      flogicZoneIII=new G4LogicalVolume(fsolidZoneIII,    fZoneIIIMaterial,   "ZoneIII_logic");
     
      fphysicZoneIII = new G4PVPlacement(0, G4ThreeVector(0,0,0), flogicZoneIII,"ZoneIII",flogicSSS,false,0,IsOverlap);

  
    // Shroud 
  fsolidShroud = new G4Sphere("Shroud_solid",
			0.0, fShroudExternalRadius,
			0.0, twopi,
			0.0, pi);

  flogicShroud = new G4LogicalVolume(fsolidShroud,
                                     fShroudMaterial,
                                     "Shroud_logic");

  fphysicShroud = new G4PVPlacement(0, G4ThreeVector(0,0,0), flogicShroud, "Shroud", flogicZoneIII, false, 0, IsOverlap);

    // ZoneII 
  fsolidZoneII = new G4Sphere("ZoneII_solid",
			0.0, fZoneIIExternalRadius,
			0.0, twopi,
			0.0, pi);

  flogicZoneII = new G4LogicalVolume(fsolidZoneII,
                                     fZoneIIMaterial,
                                     "ZoneII_logic");


  fphysicZoneII = new G4PVPlacement(0,G4ThreeVector(0,0,0.0),flogicZoneII, "ZoneII", flogicShroud, false, 0, IsOverlap);

    // Vessel
    // to do: check if the vessel is centered
  BxLog(routine) << "Vessel Origin Coordinates: " << fVesselOrigin/cm << " cm" <<  endlog ;
  
  if(!BxReadParameters::Get()->GetIsDVessel() || BxReadParameters::Get()->GetRunNumber()==0) {

    //Spherical vessel
    fsolidVessel = new G4Sphere("Vessel_solid",
			  0.0, fVesselExternalRadius,
			  0.0, twopi,
			  0.0, pi);

// ZoneI 
    fsolidZoneI = new G4Sphere("ZoneI_solid",
			  0.0, fZoneIExternalRadius,
			  0.0, twopi,
			  0.0, pi);

    flogicVessel = new G4LogicalVolume(fsolidVessel,
                                       fVesselMaterial,
                                       "Vessel_logic");

    flogicZoneI = new G4LogicalVolume(fsolidZoneI,
                                       fZoneIMaterial,
                                       "ZoneI_logic");

BxLog(routine) << "Spherical Vessel constructed. " << endlog;
BxReadParameters::Get()->SetIsDVessel(false);

  } else { // Vessel Shape from dst
 bool deformed=true;
  while (!BxReadParameters::Get()->ReadVesselGeometry()){
	  if (BxReadParameters::Get()->GetRunNumber()<=0){
		  BxLog(warning) <<"No available previous runs were found!! I construct a spherical vessel." << endlog;
	          deformed=false;
	  }
    }		
if (deformed){
fsolidDVessel = new G4GenericPolycone ("Vessel_solid",
                                       0.0, twopi,
                                    BxReadParameters::Get()->GetZCutsNumber(),
                                    BxReadParameters::Get()->GetVessel_Rout(),
				    BxReadParameters::Get()->GetVessel_Z()); 

fsolidDZoneI = new G4GenericPolycone ("ZoneI_solid",
		                        0.0, twopi,
		                       BxReadParameters::Get()->GetZCutsNumber(),
		                       BxReadParameters::Get()->GetVessel_Rin(),
		                       BxReadParameters::Get()->GetZoneI_Z());

    flogicVessel = new G4LogicalVolume(fsolidDVessel,
                                       fVesselMaterial,
                                       "Vessel_logic");

    flogicZoneI = new G4LogicalVolume(fsolidDZoneI,
                                       fZoneIMaterial,
                                       "ZoneI_logic");   
   
    BxLog(routine) << "Vessel constructed according to run number " << BxReadParameters::Get()->GetRunNumber() << " ." << endlog;

} else {

//spherical vessel
    fsolidVessel = new G4Sphere("Vessel_solid",
			  0.0, fVesselExternalRadius,
			  0.0, twopi,
			  0.0, pi);

	// ZoneI 
    fsolidZoneI = new G4Sphere("ZoneI_solid",
			  0.0, fZoneIExternalRadius,
			  0.0, twopi,
			  0.0, pi);

    flogicVessel = new G4LogicalVolume(fsolidVessel,
                                       fVesselMaterial,
                                       "Vessel_logic");

    flogicZoneI = new G4LogicalVolume(fsolidZoneI,
                                       fZoneIMaterial,
                                       "ZoneI_logic");

    BxLog(routine) << "Spherical Vessel constructed. " << endlog;

BxReadParameters::Get()->SetIsDVessel(false);

  }

  }
  
  fphysicVessel = new G4PVPlacement(0,
  				   fVesselOrigin,
                                   "Vessel",
                                   flogicVessel,
                                   fphysicZoneII,
                                   false,
                                   0, IsOverlap);


  fphysicZoneI = new G4PVPlacement(0,
  				  fVesselOrigin,
				  "ZoneI",
				  flogicZoneI,
				  fphysicVessel,
				  false,
                          	  0, IsOverlap);

   //ENDCAPS holding the Vessels
 
if(IsEndCap) {

    fsolidIVNylonRing = new G4Tubs("IVNylonRing_solid",
			        	 fNylonRingRadiusIn, 
					 fNylonRingRadiusOut,
			        	 fNylonRingHeight/2,
			        	 0, 
					 twopi
					 );
   
    //1mm thick nylon shadow on the endcap
    fsolidNylonBridge = new G4Tubs("NylonBridge_solid",
			        	 0, 
					 fNylonBridgeRadius,
			        	 fNylonBridgeHeight/2.,
			        	 0, 
					 twopi
					 );

    fsolidNylonTube=new G4Polycone("NylonTube_solid",
			 			0.0, twopi,
						BxReadParameters::Get()->GetNumOfZplanesNylonTube(),
						BxReadParameters::Get()->GetZplanesNylonTube(),
						BxReadParameters::Get()->GetRinNylonTube(),
						BxReadParameters::Get()->GetRoutNylonTube());
   
    fsolidIVTube_raw=new G4Polycone("IVTube_solid",
			 			0.0, twopi,
						BxReadParameters::Get()->GetNumOfZplanesIVTube(),
						BxReadParameters::Get()->GetZplanesIVTube(),
						BxReadParameters::Get()->GetRinIVTube(),
						BxReadParameters::Get()->GetRoutIVTube());

//Checks if Vessel is deformed and takes care of it 

if (BxReadParameters::Get()->GetIsDVessel()){
G4double hmin=0, hmax=0; //minimum and maximum heights of the IV
hmin=BxReadParameters::Get()->GetVessel_Z()[0];
hmax=BxReadParameters::Get()->GetVessel_Z()[0];

for (G4int i=0; i<BxReadParameters::Get()->GetZCutsNumber(); i++){
if (BxReadParameters::Get()->GetVessel_Z()[i]>hmax) hmax=BxReadParameters::Get()->GetVessel_Z()[i];
if (BxReadParameters::Get()->GetVessel_Z()[i]<hmin) hmin=BxReadParameters::Get()->GetVessel_Z()[i];
}
	fNylonBridgeNorthZ=hmax+fNylonBridgeHeight/2.+0.1*mm;
	fNylonBridgeSouthZ=-hmin+fNylonBridgeHeight/2.+0.1*mm;
	fIVNylonRingNorthZ=fNylonBridgeNorthZ+fNylonBridgeHeight/2.+fNylonRingHeight/2.+0.1*mm;
	fIVNylonRingSouthZ=fNylonBridgeSouthZ+fNylonBridgeHeight/2.+fNylonRingHeight/2.+0.1*mm;
	fNylonTubeNorthZ=fNylonBridgeNorthZ+fNylonBridgeHeight/2.+0.1*mm;
	fNylonTubeSouthZ=fNylonBridgeSouthZ+fNylonBridgeHeight/2.+0.1*mm;

}else {
	fNylonBridgeNorthZ=fZoneIExternalRadius+fNylonBridgeHeight/2.+0.1*mm;
	fNylonBridgeSouthZ=fZoneIExternalRadius+fNylonBridgeHeight/2.+0.1*mm;
	fIVNylonRingNorthZ=fNylonBridgeNorthZ+fNylonBridgeHeight/2.+fNylonRingHeight/2.+0.1*mm;
	fIVNylonRingSouthZ=fNylonBridgeSouthZ+fNylonBridgeHeight/2.+fNylonRingHeight/2.+0.1*mm;
	fNylonTubeNorthZ=fNylonBridgeNorthZ+fNylonBridgeHeight/2.+0.1*mm;
	fNylonTubeSouthZ=fNylonBridgeSouthZ+fNylonBridgeHeight/2.+0.1*mm;

}

	fIVTubeNorthZ=fNylonTubeNorthZ+BxReadParameters::Get()->GetIVTubeZposition()+0.1*mm;
	fIVTubeSouthZ=fNylonTubeSouthZ+BxReadParameters::Get()->GetIVTubeZposition()+0.1*mm;
	fOVNylonRingNorthZ=fZoneIIExternalRadius+fNylonRingHeight/2+0.5*mm;
	fOVNylonRingSouthZ=fZoneIIExternalRadius+fNylonRingHeight/2+0.5*mm;
        fOVTubeNorthZ=fZoneIIExternalRadius+fNylonRingHeight+0.6*mm;
        fOVTubeSouthZ=fZoneIIExternalRadius+fNylonRingHeight+0.6*mm;
	fCopperStrutNorthZ=fIVNylonRingNorthZ+fNylonRingHeight+BxReadParameters::Get()->GetCopperStrutLength()/2.+1*mm;
	fCopperStrutSouthZ=fIVNylonRingSouthZ+fNylonRingHeight+BxReadParameters::Get()->GetCopperStrutLength()/2.+1*mm;
	fCopperStrutRadialDisplacement=275*mm;

	G4RotationMatrix *rot = new G4RotationMatrix;
	rot->rotateY(180*deg);
  
	
  G4ThreeVector zTrans_IVtube_North(0,0,fIVTubeNorthZ);

  G4Transform3D transform_IVtube_North(G4RotationMatrix(),zTrans_IVtube_North);
    
  fsolidIVTube_North			    = new G4IntersectionSolid("IVTube_solid",
		  						fsolidZoneII,
								fsolidIVTube_raw,
								transform_IVtube_North);
   
  G4ThreeVector zTrans_IVtube_South(0,0,-fIVTubeSouthZ);

  G4Transform3D transform_IVtube_South(*(rot),zTrans_IVtube_South);
  
  fsolidIVTube_South			    = new G4IntersectionSolid("IVTube_solid",
		  						fsolidZoneII,
								fsolidIVTube_raw,
								transform_IVtube_South);
    fsolidOVTube_raw=new G4Polycone("OVTube_solid",
			 			0.0, twopi,
						BxReadParameters::Get()->GetNumOfZplanesOVTube(),
						BxReadParameters::Get()->GetZplanesOVTube(),
						BxReadParameters::Get()->GetRinOVTube(),
						BxReadParameters::Get()->GetRoutOVTube());


  G4ThreeVector zTrans_OVtube(0, 0, fOVTubeNorthZ);

  G4Transform3D transform_OVtube(G4RotationMatrix(),zTrans_OVtube);
    
  fsolidOVTube			    = new G4IntersectionSolid("OVTube_solid",
		  						fsolidZoneIII,
								fsolidOVTube_raw,
								transform_OVtube);

    fsolidOVNylonRing = new G4Tubs("OVNylonRing_solid",
			        	 0, 
					 fNylonRingRadiusOut,
			        	 fNylonRingHeight/2,
			        	 0, 
					 twopi
					 );

    fsolidCopperStrut = new G4Tubs("CopperStrut_solid",
			        	 0, 
					 BxReadParameters::Get()->GetCopperStrutRadius(),
			        	 BxReadParameters::Get()->GetCopperStrutLength()/2.,
			        	 0, 
					 twopi
					 );

    flogicIVNylonRing= new G4LogicalVolume(fsolidIVNylonRing,
                                       fNylonRingMaterial,
                                       "IVNylonRing_logic");
    
    flogicNylonBridge= new G4LogicalVolume(fsolidNylonBridge,
                                       fNylonBridgeMaterial,
				       "NylonBridge_logic");

    flogicOVNylonRing= new G4LogicalVolume(fsolidOVNylonRing,
                                       fNylonRingMaterial,
                                       "OVNylonRing_logic");
    
    flogicNylonTube= new G4LogicalVolume(fsolidNylonTube,
                                       fNylonTubeMaterial,
                                       "NylonTube_logic");
    
    flogicIVTube_North= new G4LogicalVolume(fsolidIVTube_North,
                                       fIVTubeMaterial,
                                       "IVTube_North_logic");

    flogicIVTube_South= new G4LogicalVolume(fsolidIVTube_South,
                                       fIVTubeMaterial,
                                       "IVTube_South_logic");
    
    flogicOVTube= new G4LogicalVolume(fsolidOVTube,
                                       fOVTubeMaterial,
                                       "OVTube_logic");

    flogicCopperStrut= new G4LogicalVolume(fsolidCopperStrut,
                                       fCopperStrutMaterial,
                                       "CopperStrut_logic");
//ENDCAPS placement


      fphysicIVNylonRingNorth = new G4PVPlacement(0,
  				    G4ThreeVector(0,0,fIVNylonRingNorthZ),
				    "IVNylonRingNorth",
				    flogicIVNylonRing,
				    fphysicZoneII,
				    false,
                          	    0, IsOverlap);
      
       fphysicNylonBridgeNorth = new G4PVPlacement(0,
  				    G4ThreeVector(0,0,fNylonBridgeNorthZ),
				    "NylonBridgeNorth",
				    flogicNylonBridge,
				    fphysicZoneII,
				    false,
                          	    0, IsOverlap);
       
      fphysicNylonBridgeSouth = new G4PVPlacement(0,
  				    G4ThreeVector(0,0,-fNylonBridgeSouthZ),
				    "NylonBridgeSouth",
				    flogicNylonBridge,
				    fphysicZoneII,
				    false,
                          	    0, IsOverlap);
      
      fphysicIVNylonRingSouth = new G4PVPlacement(0,
  				    G4ThreeVector(0,0,-fIVNylonRingSouthZ),
				    "IVNylonRingSouth",
				    flogicIVNylonRing,
				    fphysicZoneII,
				    false,
                          	    0, IsOverlap);
     

      fphysicNylonTubeNorth = new G4PVPlacement(0,
  				    G4ThreeVector(0,0,fNylonTubeNorthZ),
				    "NylonTubeNorth",
				    flogicNylonTube,
				    fphysicZoneII,
				    false,
                          	    0, IsOverlap);
      
      fphysicNylonTubeSouth = new G4PVPlacement(rot,
  				    G4ThreeVector(0,0,-fNylonTubeSouthZ),
				    "NylonTubeSouth",
				    flogicNylonTube,
				    fphysicZoneII,
				    false,
                          	    0, IsOverlap);
      
      fphysicIVTubeNorth = new G4PVPlacement(0,
  				    G4ThreeVector(0,0,0),
				    "IVTubeNorth",
				    flogicIVTube_North,
				    fphysicZoneII,
				    false,
                          	    0, IsOverlap);
      
     fphysicIVTubeSouth = new G4PVPlacement(0,
  				    G4ThreeVector(0,0,0),
				    "IVTubeSouth",
				    flogicIVTube_South,
				    fphysicZoneII,
				    false,
                          	    0, IsOverlap);

      fphysicOVNylonRingNorth = new G4PVPlacement(0,
  				    G4ThreeVector(0,0,fOVNylonRingNorthZ),
				    "OVNylonRingNorth",
				    flogicOVNylonRing,
				    fphysicZoneIII,
				    false,
                          	    0, IsOverlap);
      
      fphysicOVNylonRingSouth = new G4PVPlacement(rot,
  				    G4ThreeVector(0,0,-fOVNylonRingSouthZ),
				    "OVNylonRingSouth",
				    flogicOVNylonRing,
				    fphysicZoneIII,
				    false,
                          	    0, IsOverlap);
     
      fphysicOVTubeNorth = new G4PVPlacement(0,
  				    G4ThreeVector(0,0,0),
				    "OVTubeNorth",
				    flogicOVTube,
				    fphysicZoneIII,
				    false,
                          	    0, IsOverlap);
      
      fphysicOVTubeSouth = new G4PVPlacement(rot,
  				    G4ThreeVector(0,0,0),
				    "OVTubeSouth",
				    flogicOVTube,
				    fphysicZoneIII,
				    false,
                          	    0, IsOverlap);

//Copper Struts placement
	fNumOfStruts=16; //DO NOT change this number!
	fphysicCopperStrut = new G4VPhysicalVolume*[fNumOfStruts];
  
	G4bool first=1;
	G4ThreeVector startVersor(1,0,0), versor, perp, radialdispl;
	G4Transform3D trRotateStrut, trRTraslation, trZTraslation, trStrut;
	G4String StrutName;
	G4double theta=35;
	
	for (G4int i=0; i<fNumOfStruts/2.; i++){
    	stringstream  smu;
    	smu << i;

    	StrutName   = "Copper Strut North ";
    	StrutName  += smu.str();
	
	startVersor=G4ThreeVector(1,0,0);
	versor=startVersor.rotate(theta*deg,G4ThreeVector(0,0,1));
	radialdispl=(fCopperStrutRadialDisplacement-BxReadParameters::Get()->GetCopperStrutLength()/2.*sin(15*deg))*(versor.unit());
	perp=versor.orthogonal();

	//North Pole
	trRotateStrut= G4Rotate3D(15.0*deg,perp);
	trRTraslation=G4Translate3D(radialdispl);
	trZTraslation=G4Translate3D(G4ThreeVector(0,0,fCopperStrutNorthZ-BxReadParameters::Get()->GetCopperStrutLength()/2.*(1-cos(15*deg))));

	trStrut=trZTraslation*trRTraslation*trRotateStrut;
	fphysicCopperStrut[i] = new G4PVPlacement(trStrut,
					    StrutName,
					    flogicCopperStrut,
					    fphysicZoneII,
					    false,
					    i, IsOverlap);

	//South Pole
	trRotateStrut= G4Rotate3D(-15.0*deg,perp);
	trRTraslation=G4Translate3D(radialdispl);
	trZTraslation=G4Translate3D(G4ThreeVector(0,0,-fCopperStrutSouthZ+BxReadParameters::Get()->GetCopperStrutLength()/2.*(1-cos(15*deg))));
    	
	StrutName   = "Copper Strut South ";
    	StrutName  += smu.str();
	
	trStrut=trZTraslation*trRTraslation*trRotateStrut;
	
	fphysicCopperStrut[int(fNumOfStruts/2.+i)] = new G4PVPlacement(trStrut,
					    StrutName,
					    flogicCopperStrut,
					    fphysicZoneII,
					    false,
					    fNumOfStruts+i, IsOverlap);
	if (first) theta=theta+20;
	else theta=theta+70;

	first=!first;
	}
   BxLog(routine) << "EndCaps holding the inner vessel constructed." << endlog; 
}				    
 
 /* if(IsSensitive) {
    G4SDManager* SDman = G4SDManager::GetSDMpointer();				
    BxSensitiveDetector* SensistiveDetector = new  BxSensitiveDetector( "Bulk" );
    SDman->AddNewDetector( SensistiveDetector );
    flogicZoneI->SetSensitiveDetector(SensistiveDetector); 
  }
  */ 
//------------------------------PMTs-------------------------------------------------
//Concentrators
    
if (fPMTFlag > 1 && !fNoConcentrators){
	    BxReadParameters::Get()->ReadLGParam(); 
	    
	    G4GenericPolycone* extLG=new G4GenericPolycone("extLG",0.0,twopi,
			    BxReadParameters::Get()->GetLGZCutsNumber(),
			    BxReadParameters::Get()->GetLGOuterRadius(),
			    BxReadParameters::Get()->GetLGZValue());

	    G4GenericPolycone* midLG=new G4GenericPolycone("midLG",0.0,twopi,
			    BxReadParameters::Get()->GetLGZCutsNumber(),
			    BxReadParameters::Get()->GetLGMidRadius(),
			    BxReadParameters::Get()->GetLGZValue());

	    G4GenericPolycone* innLG=new G4GenericPolycone("innLG",0.0,twopi,
			    BxReadParameters::Get()->GetLGZCutsNumber(),
			    BxReadParameters::Get()->GetLGInnerRadius(),
			    BxReadParameters::Get()->GetLGZValue());

	fsolidExternalLightGuide=new G4SubtractionSolid("PMTExternalLightGuide_solid",extLG,midLG);
	fsolidInternalLightGuide=new G4SubtractionSolid("PMTInternalLightGuide_solid",midLG,innLG);

}
  
if (fPMTFlag == 1 || fPMTFlag == 2){
  //simple disk PMT construction with and without concentrator

    //basic cylinder for simple PMT
    G4Tubs* fsolidPMTcylinder          = new G4Tubs("PMT_simple_cylinder",
			        	 0.0, 
					 fPMTSizeDiskR,
			        	 (fPMTSizeDiskH+2.*mm)/2,
			        	 0.0, 
					 twopi
					 );
    
   G4ThreeVector zTrans(0, 0, -(fZoneIIIExternalRadius-(fPMTSizeDiskH-2*mm)/2));

    G4Transform3D transform(G4RotationMatrix(),zTrans);
    
    fsolidPMTdisk=new G4IntersectionSolid("PMT_solid",
		    			  fsolidZoneIII,
					  fsolidPMTcylinder,
					  transform);

    //basic cylinder for simple disk PMT shielding
    G4Tubs* fsolidPMTcylinderShield = new G4Tubs("PMT_symple_cylinder_Shield", 
	                                           fPMTSizeDiskR+0.1*mm,
	                                           fPMTSizeDiskR+1.1*mm,
	                                           (fPMTSizeDiskH+2.*mm)/2,
	                                           0.,
	                                           twopi);

    fsolidPMTShielddisk= new G4IntersectionSolid("PMTShield_solid",
		    			  fsolidZoneIII,
					  fsolidPMTcylinderShield,
					  transform);

    
}  else if (fPMTFlag == 3){
//simple spherical PMT construction + concentrator

  fsolidPMTbase_cylinder		    = new G4Polycone("PMT_solid_base_cylinder",
		  					0.0,
							twopi,
							BxReadParameters::Get()->GetNumOfZplanesHousing(),
							BxReadParameters::Get()->GetZplanesHousing(),
							BxReadParameters::Get()->GetRinHousing(),
							BxReadParameters::Get()->GetRoutHousing()
							);
   
  G4ThreeVector zTrans(0, 0, -fZoneIIIExternalRadius-2.*mm);

  G4Transform3D transform(G4RotationMatrix(),zTrans);
    
  fsolidPMTHousing			    = new G4IntersectionSolid("PMT_base_solid",
		  						fsolidZoneIII,
								fsolidPMTbase_cylinder,
								transform);

  fsolidPMTShield		            = new G4Polycone("PMT_solid_shield",
		  					0.0,
							twopi,
							BxReadParameters::Get()->GetNumOfZplanesShield(),
							BxReadParameters::Get()->GetZplanesShield(),
							BxReadParameters::Get()->GetRinShield(),
							BxReadParameters::Get()->GetRoutShield()
							);

  fsolidPMTsphere	           = new G4Sphere("PMT_solid_sphere",
			               0.0, 
				       fPMTSizeRmax,
			               0.0,   // phi min
			               twopi, // phi max
				       0.0, 
				       twopi);
  
  //Box to create a portion of the photocatode sphere  
  G4Box* SubtractionBox=new G4Box("SubtractionBox", 1*m, 1*m, 1*m);
 
  //G4ThreeVector for the translation of the SubtractionBox
  G4ThreeVector zTransBox(0.,0.,-1.*m+BxReadParameters::Get()->GetPMTZShift());

  //Null rotation for boolean operations
  G4RotationMatrix *yRot = new G4RotationMatrix;
  yRot->rotateY(0.);
  
  fsolidPMT=new G4SubtractionSolid("PMT_solid",fsolidPMTsphere,SubtractionBox,yRot,zTransBox);


//Steel Ring for PMTs without concentrator
  fsolidPMTring= new G4Polycone("PMT_ring_solid",
		  		0.0,
				twopi,
				BxReadParameters::Get()->GetNumOfZplanesPMTring(),
				BxReadParameters::Get()->GetZplanesPMTring(),
				BxReadParameters::Get()->GetRinPMTring(),
				BxReadParameters::Get()->GetRoutPMTring()
		  );

} else {
	BxLog(fatal) << "fPMTFlag = " << fPMTFlag << " not implemented." << endlog;
	}

//Light Concentrators' logical volumes
  if (fPMTFlag >1 && !fNoConcentrators){

  flogicExternalLightGuide  = new G4LogicalVolume(fsolidExternalLightGuide,
			  fExternalLightGuideMaterial,
			  "ExternalLightGuide_logic");
  
  flogicInternalLightGuide  = new G4LogicalVolume(fsolidInternalLightGuide,
			  fInternalLightGuideMaterial,
			  "InternaLightGuide_logic");
  }

if (fPMTFlag == 1 || fPMTFlag == 2){

	 //simple disk PMT construction made of bialkali

  flogicPMT	        = new G4LogicalVolume(fsolidPMTdisk,
                           fPMTMaterial,
                           "PMT_logic");

  flogicPMTShield       = new G4LogicalVolume(fsolidPMTShielddisk,
			  fPMTShieldMaterial,
			  "PMTShield_logic");

} else if (fPMTFlag == 3 ){

	//simple PMT construction

  flogicPMT 	        = new G4LogicalVolume(fsolidPMT,
                           fPMTMaterial,
                           "PMT_logic");


  flogicPMTShield       = new G4LogicalVolume(fsolidPMTShield,
			  fPMTShieldMaterial,
			  "PMTShield_logic");

  flogicPMTHousing       = new G4LogicalVolume(fsolidPMTHousing,
			fPMTHousingMaterial,  
		  	"PMTHousing_logic");

  flogicPMTring        = new G4LogicalVolume(fsolidPMTring,
			  fPMTringMaterial,
			  "PMTring_logic");

}

if (fPMTDistributionFlag==1)
  fNumberOfPMT = fPMTPosition.size();
  //Memory allocation for concentrators

 if (fPMTFlag > 1 && !fNoConcentrators){
  fphysicExternalLightGuide = new G4VPhysicalVolume*[fNumberOfPMT-fNumberOfVetoPMT];
  fphysicInternalLightGuide = new G4VPhysicalVolume*[fNumberOfPMT-fNumberOfVetoPMT];
  
 }
 if (fPMTFlag == 1 || fPMTFlag==2) { 
//disk PMT

 fphysicPMT = new G4VPhysicalVolume*[fNumberOfPMT];
 fphysicPMTShield = new G4VPhysicalVolume*[fNumberOfPMT];

 } else if (fPMTFlag==3){

 fphysicPMT = new G4VPhysicalVolume*[fNumberOfPMT];
 fphysicPMTShield = new G4VPhysicalVolume*[fNumberOfPMT];
 fphysicPMTHousing = new G4VPhysicalVolume*[fNumberOfPMT];
//Rings for PMTs without concentrators
 if (!fNoConcentrators)
	 fphysicPMTring=new G4VPhysicalVolume*[fNumberOfVetoPMT];
 else
	 fphysicPMTring=new G4VPhysicalVolume*[fNumberOfPMT];

} 
  
  G4ThreeVector ps0,
  		ax;

  G4Transform3D tr2,
  		tr3,
		tr4,
		tr5,
		tr6,
		tr7,
		tr8,
		tr9,
		tr10;
  
  G4String	PmtN,
                PmtSN,
		PmtHN,
		PmtELGN,
		PmtILGN,
		PmtRN;
   
G4int k_conc=0,k_noconc=0;
if(!noPMTs) {
	if (fPMTFlag == 1 || fPMTFlag == 2){ 
    //simple disk PMT construction

	    if (fPMTDistributionFlag==0){
		    G4int cycles;
		    if (fNumberOfPMT+1<pow((int)sqrt(fNumberOfPMT-1),2))
			    cycles=(int)sqrt(fNumberOfPMT-2);
		    else
			    cycles=(int)sqrt(fNumberOfPMT-1);
		    G4int cont=0;
		    for (G4int i=0; i<cycles; i++){

			    for (G4int j=0; j<cycles; j++){
				    G4ThreeVector vec(1,0,0);
				    vec.setPhi(j*2.*pi*rad/cycles);
				    vec.setTheta(abs(acos(-1+2.*i/cycles)*rad));
				    vec.unit();
				    ax=vec.orthogonal();
				    tr2 = G4Rotate3D( G4Point3D ( 0.0*cm,    0.0*cm,   -1.0*cm ),
						    G4Point3D ( 0.0*cm,    1.0*cm,    0.0*cm ),
						    G4Point3D ( vec.x()*m, vec.y()*m, vec.z()*m ),
						    G4Point3D ( ax.x()*m,  ax.y()*m,  ax.z()*m) 
						    );

				    stringstream  str;
				    str << cont;
				    PmtN = "PMT Cathode ";       
				    PmtN += str.str();


 fphysicPMT[cont] = new G4PVPlacement(tr2,
				    flogicPMT,
				    PmtN,
                                   flogicZoneIII,
	                            false,
 				cont, IsOverlap);

fphysicPMTShield[cont] = new G4PVPlacement(tr2,
				    flogicPMTShield,
				    PmtSN,
                                   flogicZoneIII,
	                            false,
 				cont, IsOverlap);

if (fPMTFlag == 2){
//all PMTs have a concentrator
vec.setMag(fZoneIIIExternalRadius-fPMTSizeDiskH);//-230./2*mm);
tr3= G4Translate3D(vec);
tr3=tr3*tr2;

	       PmtELGN = "PMT External Concentrator ";
	       PmtELGN += str.str();
	 
fphysicExternalLightGuide[cont]=new G4PVPlacement(tr3,
					flogicExternalLightGuide,
					PmtELGN,
					flogicZoneIII,
					false,
					cont,
					IsOverlap);

	       PmtILGN = "PMT Internal Concentrator ";
	       PmtILGN += str.str();
	 
fphysicInternalLightGuide[cont]=new G4PVPlacement(tr3,
					flogicInternalLightGuide,
					PmtILGN,
					flogicZoneIII,
					false,
					cont,
					IsOverlap);

}
cont++;	
	}
}
} else if(fPMTDistributionFlag==1){

	for (G4int i=0; i<fNumberOfPMT; i++){
//costruzione PMT da file

          stringstream  str;
	  str << i;

	  PmtN = "PMT Cathode ";       
	  PmtN += str.str();

	  PmtSN = "PMT Shield ";
	  PmtSN += str.str();
	  	

	ps0 = fPMTPosition[i];
	ps0.unit();
	ax = ps0.orthogonal();
	tr2 = G4Rotate3D( G4Point3D ( 0.0*cm,    0.0*cm,   -1.0*cm ),
			  G4Point3D ( 0.0*cm,    1.0*cm,    0.0*cm ),
			  G4Point3D ( ps0.x()*m, ps0.y()*m, ps0.z()*m ),
			  G4Point3D ( ax.x()*m,  ax.y()*m,  ax.z()*m) 
			);
	
	    fphysicPMT[i] = new G4PVPlacement(tr2,
				    flogicPMT,
				    PmtN,
                                   flogicZoneIII,
	                            false,
                                    i, IsOverlap);

	    fphysicPMTShield[i] = new G4PVPlacement(tr2,
				    flogicPMTShield,
				    PmtSN,
                                   flogicZoneIII,
	                            false,
                                    i, IsOverlap);
	    
	    if (fPMTFlag == 2) {
ps0.setMag(fZoneIIIExternalRadius-fPMTSizeDiskH);
tr3=G4Translate3D(ps0);
tr3=tr3*tr2;

	if(!fVetoPmt.at(i)) {
               stringstream  sst;
  	       sst << k_conc;
	       PmtELGN = "PMT External Concentrator ";
	       PmtELGN += str.str();
	 
	 fphysicExternalLightGuide[k_conc]=new G4PVPlacement(tr3,
					flogicExternalLightGuide,
					PmtELGN,
					flogicZoneIII,
					false,
					k_conc,
					IsOverlap);
 
	       PmtILGN = "PMT Internal Concentrator ";
	       PmtILGN += str.str();
	 
	 fphysicInternalLightGuide[k_conc]=new G4PVPlacement(tr3,
					flogicInternalLightGuide,
					PmtILGN,
					flogicZoneIII,
					false,
					k_conc,
					IsOverlap);
	 k_conc++; 
	}	
 }
}

}


    }  else if (fPMTFlag == 3){
//Simple PMT	    
	for (G4int i=0; i<fNumberOfPMT; i++){

          stringstream  str;
	  str << i;

	  PmtN = "PMT Cathode ";       
	  PmtN += str.str();

	  PmtSN = "PMT Shield ";
	  PmtSN += str.str();
	  
	  PmtHN = "PMT Housing ";
	  PmtHN += str.str();
	  	

	ps0 = fPMTPosition[i];
	ps0.unit();
	ax = ps0.orthogonal();
	tr2 = G4Rotate3D( G4Point3D ( 0.0*cm,    0.0*cm,   -1.0*cm ),
			  G4Point3D ( 0.0*cm,    1.0*cm,    0.0*cm ),
			  G4Point3D ( ps0.x()*m, ps0.y()*m, ps0.z()*m ),
			  G4Point3D ( ax.x()*m,  ax.y()*m,  ax.z()*m) 
			);
	
	//Arrangement of PMTs Housing
	fphysicPMTHousing[i] = new G4PVPlacement(tr2,
				    flogicPMTHousing,
				    PmtHN,
                                    flogicZoneIII,
	                            false,
                                    i, IsOverlap);
      
     //Arrangement of PMTs shielding
        G4double theShieldShift=fZoneIIIExternalRadius
				-BxReadParameters::Get()->GetPMTHousingShift()
				-BxReadParameters::Get()->GetPMTShieldZPosition();
	G4ThreeVector ps1(0,0,-theShieldShift);
	
	tr3=G4Translate3D(ps1);
	tr4=tr2*tr3;

	fphysicPMTShield[i] = new G4PVPlacement(tr4,
				    flogicPMTShield,
				    PmtSN,
                                   flogicZoneIII,
	                            false,
                                    i, IsOverlap);
    
	//Arrangement of PMT photocatode	    
        G4double thePMTShift=fZoneIIIExternalRadius+BxReadParameters::Get()->GetPMTZShift()
			 	-BxReadParameters::Get()->GetPMTHousingShift()
				-BxReadParameters::Get()->GetPMTLength()
				-BxReadParameters::Get()->GetPMTHousingThickness();

	G4ThreeVector ps2(0,0,-thePMTShift);
	
	tr5=G4Translate3D(ps2);
	tr6=tr2*tr5;

	fphysicPMT[i] = new G4PVPlacement(tr6,
				    flogicPMT,
				    PmtN,
                                   flogicZoneIII,
	                            false,
                                    i, IsOverlap);
	//Arrangement of Concentrators and PMT rings for PMTs without the concentrator.
	//If concentrators are not built, then all PMTs have rings.
	if (!fNoConcentrators && !fVetoPmt.at(i)) {

               stringstream  sst;
  	       sst << k_conc;
        
	G4double theLightGuideShift=fZoneIIIExternalRadius
			 	-BxReadParameters::Get()->GetPMTHousingShift()
				-BxReadParameters::Get()->GetPMTLength()
				-BxReadParameters::Get()->GetPMTHousingThickness()
				-BxReadParameters::Get()->GetPMTLightGuideZShift();

	G4ThreeVector ps3(0,0,-theLightGuideShift);
	
	tr7=G4Translate3D(ps3);
	tr8=tr2*tr7;
	 
	       PmtELGN = "PMT External Concentrator ";
	       PmtELGN += str.str();
	 fphysicExternalLightGuide[k_conc]=new G4PVPlacement(tr8,
					flogicExternalLightGuide,
					PmtELGN,
					flogicZoneIII,
					false,
					k_conc,
					IsOverlap);
 	
	       PmtILGN = "PMT Internal Concentrator ";
	       PmtILGN += str.str();
	 fphysicInternalLightGuide[k_conc]=new G4PVPlacement(tr8,
					flogicInternalLightGuide,
					PmtILGN,
					flogicZoneIII,
					false,
					k_conc,
					IsOverlap);
	 k_conc++; 
	}else if (!fNoConcentrators && fVetoPmt.at(i)){

               stringstream  sst;
  	       sst << k_noconc;
        
	G4double thePMTringShift=fZoneIIIExternalRadius
			 	-BxReadParameters::Get()->GetPMTHousingShift()
				-BxReadParameters::Get()->GetPMTLength()
				-BxReadParameters::Get()->GetPMTHousingThickness()
				-BxReadParameters::Get()->GetPMTringZShift();

	G4ThreeVector ps4(0,0,-thePMTringShift);
	
	tr9=G4Translate3D(ps4);
	tr10=tr2*tr9;
	 
	       PmtRN = "PMT ring ";
	       PmtRN += str.str();
	 
	       fphysicPMTring[k_noconc]=new G4PVPlacement(tr10,
					flogicPMTring,
					PmtRN,
					flogicZoneIII,
					false,
					k_noconc,
					IsOverlap);
	k_noconc++;
	} else if (fNoConcentrators){
               stringstream  sst;
  	       sst << k_noconc;
        
	G4double thePMTringShift=fZoneIIIExternalRadius
			 	-BxReadParameters::Get()->GetPMTHousingShift()
				-BxReadParameters::Get()->GetPMTLength()
				-BxReadParameters::Get()->GetPMTHousingThickness()
				-BxReadParameters::Get()->GetPMTringZShift();

	G4ThreeVector ps4(0,0,-thePMTringShift);
	
	tr9=G4Translate3D(ps4);
	tr10=tr2*tr9;
	 
	       PmtRN = "PMT ring ";
	       PmtRN += str.str();
	 
	       fphysicPMTring[k_noconc]=new G4PVPlacement(tr10,
					flogicPMTring,
					PmtRN,
					flogicZoneIII,
					false,
					k_noconc,
					IsOverlap);
	k_noconc++;

		}
    	}
  }
}

//
//_______________PMT of Bx Veto Outer Detector_____________________________
//

  fsolidPMTmusphere	           = new G4Sphere("PMTmu_solid_sphere",
			               0.0, 
				       fPMTSizeRmax,
			               0.0,   // phi min
			               twopi, // phi max
				       0.0, 
				       twopi);
  
  //Box to create a portion of the photocatode sphere  
  G4Box* SubtractionBox=new G4Box("SubtractionBox", 1*m, 1*m, 1*m);
 
  //G4ThreeVector for the translation of the SubtractionBox
  G4ThreeVector zTransBox(0.,0.,-1.*m+BxReadParameters::Get()->GetPMTZShift());

  //Null rotation for boolean operations
  G4RotationMatrix *yRot = new G4RotationMatrix;
  yRot->rotateY(0.);
  
  fsolidPMTmu=new G4SubtractionSolid("PMTmu_solid",fsolidPMTmusphere,SubtractionBox,yRot,zTransBox);

  	flogicPMTmu 	  = new G4LogicalVolume(fsolidPMTmu,
                                fPMTmuMaterial,
                                "PMTmu_logic");


  	fsolidPMTmuFlange  = new G4Tubs("PMTmuFlange_solid",
				fPMTmuFlangeSizeRmin, fPMTmuFlangeSizeRmax,
				fPMTmuFlangeSizeZ/2.,
				0.0, twopi);

  	flogicPMTmuFlange = new G4LogicalVolume(fsolidPMTmuFlange,
                                   fPMTmuFlangeMaterial,
                                   "PMTFlange_logic");


  	fsolidPMTmuShield  = new G4Polycone("PMTmuShield_solid",
				0.0, twopi,
				BxReadParameters::Get()->GetNumOfZplanesShieldmu(),
				BxReadParameters::Get()->GetZplanesShieldmu(),
				BxReadParameters::Get()->GetRinShieldmu(),
				BxReadParameters::Get()->GetRoutShieldmu());

  	flogicPMTmuShield = new G4LogicalVolume(fsolidPMTmuShield,
                                   fPMTmuShieldMaterial,
                                   "PMTShield_logic");


//Loop to count number of veto PMT on the sphere and number of PMTs on the ground
	fNumberOfPMTmu = fPMT_VetoOD_Position.size();	
        fNumberOfPMTmuSSS=0;
	fNumberOfPMTmuGround=0;

	for (G4int it=0; it<fNumberOfPMTmu; it++){
         if (fVetoOD_Pmt[it] == 0) fNumberOfPMTmuSSS++;
	 else			  fNumberOfPMTmuGround++;
	}

       if (fNumberOfPMTmuSSS+fNumberOfPMTmuGround!=fNumberOfPMTmu)
	BxLog(fatal) << "Error. Mismatch in the number of OD PMTs." << endlog;

if(IsTank && !noPMTs) {
	BxLog(routine) <<"NumberOfPMTmu in Outer Detector  =" << fNumberOfPMTmu << endlog;
	BxLog(trace) <<"NumberOfPMTmu placed on the SSS  =" <<   fNumberOfPMTmuSSS << endlog;
	BxLog(trace) <<"NumberOfPMTmu placed on the ground  =" << fNumberOfPMTmuGround << endlog;
}



  //_________________This is Volumes arrangement___________________________

  	fphysicPMTmu = new G4VPhysicalVolume*[fNumberOfPMTmu];

  	fphysicPMTmuFlange = new G4VPhysicalVolume*[fNumberOfPMTmu];

	fphysicPMTmuShield = new G4VPhysicalVolume*[fNumberOfPMTmuGround];

 //_________________End of Volumes arrangement___________________________


  G4ThreeVector  ps0mu,ps0Dmu, ps0D1mu,ps1mu,ps2mu,ps3mu,psc2mu,psf2mu,pss2mu, psc3mu, pss3mu,psf3mu;
  G4ThreeVector axmu,ax2mu,theZtranslation(0.0,0.0,fDisalignment); //this correction is due to the non alignment between SSS and Tank centres

  G4Transform3D trSSSmu,trSSSFlangemu, trSSSCathodemu,tr45Cathode, tr45Flange, tr45Shield,tr45degrees, trGroundCathode, trGroundShield, trGroundFlange ;
  G4Transform3D trSSSc1mu, trSSSf1mu, tr45c1, tr45f1, tr45s1,tr45c2,tr45f2,trgc1,trgf1,trgs1, trDisalignment;

  trDisalignment=G4Translate3D(theZtranslation);

  G4String	PmtmuName,
		PmtmuFlangeName,
		PmtmuShieldName;

  for(G4int i=0;i<fNumberOfPMTmu;i++)	{

    //	PMT Hole position (from reading data file part)
    ps0mu = fPMT_VetoOD_Position[i];
    ps1mu = ps0mu;
    ps2mu = ps0mu;
    ps3mu = ps0mu;

    ps0mu=ps0mu.unit();
    axmu = ps0mu.orthogonal();

    trSSSmu = G4Rotate3D( G4Point3D ( 0*cm, 0*cm,1.0*cm ),
		    G4Point3D ( 0*cm,-1.0*cm, 0*cm ),
		    G4Point3D ( ps0mu.x()*m, ps0mu.y()*m, ps0mu.z()*m ),
		    G4Point3D (axmu.x()*m,axmu.y()*m,axmu.z()*m) );

    //	PMT Cathode and Flange position for PMTs on the sphere
    
    ps0Dmu = (-BxReadParameters::Get()->GetPMTZShift()+fSSSExternalRadius+fTyvekSphereThickness)*ps0mu.unit();
    trSSSc1mu= G4Translate3D(ps0Dmu);

    ps0D1mu=(fSSSExternalRadius+fTyvekSphereThickness+fPMTmuFlangeSizeZ/2.)*ps0mu.unit();
    trSSSf1mu= G4Translate3D(ps0D1mu);

    trSSSCathodemu=trDisalignment*trSSSc1mu*trSSSmu;
    trSSSFlangemu=trDisalignment*trSSSf1mu*trSSSmu;

    // PMT Cathode, Flange and Shield for PMTs on the 45degrees slope
    
    //shield	    
    pss2mu=ps1mu;
    pss2mu.setZ(-fTankDawnSizeZ+1*m+fTyvekThickness);
    tr45s1= G4Translate3D(pss2mu);
    
    //rotation for the 45degrees slope
    ps2mu.setZ(0);
    ax2mu=ps2mu.orthogonal();
    tr45degrees= G4Rotate3D(45.0*deg,ax2mu);

    //transformation to put the flange and the cathode on the top of the shield
    G4ThreeVector ax3mu(0,0,1);
    ps3mu=ax3mu.rotate(45.0*deg,ax2mu);
    
    //flange
    ps3mu=(1*mm+fPMTmuFlangeSizeZ/2.
		    +BxReadParameters::Get()->GetZplanesShieldmu()[BxReadParameters::Get()->GetNumOfZplanesShieldmu()-1])*ps3mu.unit();
    tr45f1=G4Translate3D(ps3mu);
    //cathode
    ps3mu=(1*mm-BxReadParameters::Get()->GetPMTZShift()+
		    BxReadParameters::Get()->GetZplanesShieldmu()[BxReadParameters::Get()->GetNumOfZplanesShieldmu()-1])*ps3mu.unit();
    tr45c1=G4Translate3D(ps3mu);
    
    //final composition
    tr45Cathode=tr45c1*tr45s1*tr45degrees;
    tr45Flange=tr45f1*tr45s1*tr45degrees;
    tr45Shield=tr45s1*tr45degrees;
   
   
    // PMT Cathode, Flange and Shield for PMTs on the Floor
    
    //shield	    
    pss3mu=ps1mu;
    pss3mu.setZ(-fTankDawnSizeZ+fTyvekThickness);
    trgs1= G4Translate3D(pss3mu);
    //flange
    psf3mu=ps1mu;
    psf3mu.setZ(1*mm+pss3mu.z()+fPMTmuFlangeSizeZ/2.+BxReadParameters::Get()->GetZplanesShieldmu()[BxReadParameters::Get()->GetNumOfZplanesShieldmu()-1]);
    trgf1= G4Translate3D(psf3mu);
    //cathode
    psc3mu=ps1mu;
    psc3mu.setZ(psf3mu.z()-fPMTmuFlangeSizeZ/2.-BxReadParameters::Get()->GetPMTZShift());
    trgc1= G4Translate3D(psc3mu);
    
    //final composition
    trGroundCathode=trgc1;
    trGroundFlange=trgf1;
    trGroundShield=trgs1;
   
    stringstream  smu;
    smu << i;

    PmtmuName   = "PMTmu-Cathode ";
    PmtmuName  += smu.str();

    PmtmuFlangeName  = "PMTmu-Flange ";
    PmtmuFlangeName += smu.str();

    PmtmuShieldName  = "PMTmu-Shield ";
    PmtmuShieldName += smu.str();

    G4int counter=0;

    if(IsTank && !noPMTs) {
      
	if(fVetoOD_Pmt[i] == 0)	{   //PMTs on the SSS	

	  fphysicPMTmu[i] = new G4PVPlacement(trSSSCathodemu,
					    PmtmuName,
					    flogicPMTmu,
					    fphysicTank,
					    false,
					    i, IsOverlap);

	  fphysicPMTmuFlange[i] = new G4PVPlacement(trSSSFlangemu,
					       PmtmuFlangeName,
					       flogicPMTmuFlange,
					       fphysicTank,
					       false,
					       i, IsOverlap);	

	} else	if (fVetoOD_Pmt[i] == 1){	//PMTs on the 45degrees slope


      fphysicPMTmu[i] = new G4PVPlacement(tr45Cathode,
			    PmtmuName,
			    flogicPMTmu,
                	    fphysicTank,
	        	    false,
                	    i, IsOverlap);
      
      fphysicPMTmuFlange[i] = new G4PVPlacement(tr45Flange,
			    PmtmuFlangeName,
			    flogicPMTmuFlange,
			    fphysicTank,
			    false,
			    i, IsOverlap);	
      
	
      fphysicPMTmuShield[counter] = new G4PVPlacement(tr45Shield,
			    PmtmuShieldName,
			    flogicPMTmuShield,
			    fphysicTank,
			    false,
			    i, IsOverlap);	
        counter++;	
 	} else	if (fVetoOD_Pmt[i] == 2){	//PMTs on the ground


      fphysicPMTmu[i] = new G4PVPlacement(trGroundCathode,
			    PmtmuName,
			    flogicPMTmu,
                	    fphysicTank,
	        	    false,
                	    i, IsOverlap);
      
      fphysicPMTmuFlange[i] = new G4PVPlacement(trGroundFlange,
			    PmtmuFlangeName,
			    flogicPMTmuFlange,
			    fphysicTank,
			    false,
			    i, IsOverlap);	
      
	
      fphysicPMTmuShield[counter] = new G4PVPlacement(trGroundShield,
			    PmtmuShieldName,
			    flogicPMTmuShield,
			    fphysicTank,
			    false,
			    i, IsOverlap);	
	
       counter++;
      } 

    }
  } 

  //--------------------------- 
  //	Surfaces
  //--------------------------- 

  G4double	*pointer1;
  G4double	*pointer2;
  G4int		theNum;
  //--------------------SSS surface--------------------------
   fOpSSSSurface = new G4OpticalSurface("OpSSSSurface");
   fOpSSSSurface->SetType(dielectric_metal);
   fOpSSSSurface->SetModel(unified);
   fOpSSSSurface->SetFinish(ground);
   
  G4MaterialPropertiesTable *OpSSSSurfaceProperty = new G4MaterialPropertiesTable();
  
  G4double vector_size=BxPropertyCollection::Get()->GetSSSReflectivity().size();
  pointer1 =BxPropertyCollection::Get()->GetPhotonEnergySSSRef();
  pointer2 =BxPropertyCollection::Get()->GetSSSReflectivity();
  OpSSSSurfaceProperty -> AddProperty("REFLECTIVITY",pointer1,pointer2,vector_size);

  pointer1 = BxPropertyCollection::Get()->GetPhotonEnergySSSEff();
  pointer2 = BxPropertyCollection::Get()->GetSSSEfficiency();
  OpSSSSurfaceProperty -> AddProperty("EFFICIENCY",pointer1,pointer2,vector_size);
 
   vector_size  = BxPropertyCollection::Get()->GetSSSspecLOBE().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergySSSspecLOBE();
   pointer2 = BxPropertyCollection::Get()->GetSSSspecLOBE();
   OpSSSSurfaceProperty->AddProperty("SPECULARLOBECONSTANT", pointer1, pointer2, vector_size);
  
  vector_size  = BxPropertyCollection::Get()->GetSSSspecSPIKE().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergySSSspecSPIKE();
   pointer2 = BxPropertyCollection::Get()->GetSSSspecSPIKE();
   OpSSSSurfaceProperty->AddProperty("SPECULARSPIKECONSTANT", pointer1, pointer2, vector_size);

  vector_size  = BxPropertyCollection::Get()->GetSSSbackSCATT().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergySSSbackSCATT();
   pointer2 = BxPropertyCollection::Get()->GetSSSbackSCATT();
   OpSSSSurfaceProperty->AddProperty("BACKSCATTERCONSTANT", pointer1, pointer2, vector_size);

  fOpSSSSurface -> SetMaterialPropertiesTable(OpSSSSurfaceProperty);
   
  fSSSSurface = new G4LogicalBorderSurface("SSSSurface",
					fphysicZoneIII	,
						fphysicSSS,
						fOpSSSSurface);


  //--------------------ZoneIII surface--------------------------
   fOpZoneIIIInSurface = new G4OpticalSurface("OpZoneIIIInSurface");
   fOpZoneIIIInSurface->SetType(dielectric_dielectric);
   fOpZoneIIIInSurface->SetModel(unified);
   fOpZoneIIIInSurface->SetFinish(polished);
   
   fZoneIIIInSurface = new G4LogicalBorderSurface("ZoneIIIInSurface",
				  fphysicZoneIII,
				  fphysicShroud,
				  fOpZoneIIIInSurface);
   

  //
  //--------------------Shroud surfaces--------------------------
   fOpShroudOutSurface = new G4OpticalSurface("OpShroudOutSurface");
   fOpShroudOutSurface->SetType(dielectric_dielectric);
   fOpShroudOutSurface->SetModel(unified);
   fOpShroudOutSurface->SetFinish(polished);
   
   fShroudOutSurface = new G4LogicalBorderSurface("ShroudOutSurface",
				  fphysicShroud,
				  fphysicZoneIII,
				  fOpShroudOutSurface);


   fOpShroudInSurface = new G4OpticalSurface("OpShroudInSurface");
   fOpShroudInSurface->SetType(dielectric_dielectric);
   fOpShroudInSurface->SetModel(unified);
   fOpShroudInSurface->SetFinish(polished);

   fShroudInSurface = new G4LogicalBorderSurface("ShroudInSurface",
						  fphysicShroud,
						  fphysicZoneII,
						  fOpShroudInSurface);

    //
    //--------------------ZoneII surfaces--------------------------

   fOpZoneIIOutSurface = new G4OpticalSurface("OpZoneIIOutSurface");
   fOpZoneIIOutSurface->SetType(dielectric_dielectric);
   fOpZoneIIOutSurface->SetModel(unified);
   fOpZoneIIOutSurface->SetFinish(polished);

   fZoneIIOutSurface = new G4LogicalBorderSurface("ZoneIIOutSurface",
						  fphysicZoneII,
						  fphysicShroud,
						  fOpZoneIIOutSurface);
   
   fOpZoneIIInSurface = new G4OpticalSurface("OpZoneIIInSurface");
   fOpZoneIIInSurface->SetType(dielectric_dielectric);
   fOpZoneIIInSurface->SetModel(unified);
   fOpZoneIIInSurface->SetFinish(polished);

   fZoneIIInSurface = new G4LogicalBorderSurface("ZoneIIInSurface",
						  fphysicZoneII,
						  fphysicVessel,
						  fOpZoneIIInSurface);

    //
    //--------------------Vessel surfaces--------------------------

   fOpVesselOutSurface = new G4OpticalSurface("OpVesselOutSurface");
   fOpVesselOutSurface->SetType(dielectric_dielectric);
   fOpVesselOutSurface->SetModel(unified);
   fOpVesselOutSurface->SetFinish(polished);
   
   fVesselOutSurface = new G4LogicalBorderSurface("VesselOutSurface",
				  fphysicVessel,
				  fphysicZoneII,
				  fOpVesselOutSurface);
  
   fOpVesselInSurface = new G4OpticalSurface("OpVesselInSurface");
   fOpVesselInSurface->SetType(dielectric_dielectric);
   fOpVesselInSurface->SetModel(unified);
   fOpVesselInSurface->SetFinish(polished);

   fVesselInSurface = new G4LogicalBorderSurface("VesselInSurface",
						  fphysicVessel,
						  fphysicZoneI,
						  fOpVesselInSurface);

    //--------------------ZoneI surface--------------------------

   fOpZoneISurface = new G4OpticalSurface("OpZoneISurface");
   fOpZoneISurface->SetType(dielectric_dielectric);
   fOpZoneISurface->SetModel(unified);
   fOpZoneISurface->SetFinish(polished);

   fZoneISurface = new G4LogicalBorderSurface("ZoneISurface",
						  fphysicZoneI,
						  fphysicVessel,
						  fOpZoneISurface);

   if(IsTank) {

   //__________________Surfaces of Outer Detector___________

//Tyvek optical properties
  G4MaterialPropertiesTable *myST10 = new G4MaterialPropertiesTable();
  
  vector_size=BxPropertyCollection::Get()->GetTyvekReflectivity().size();
  pointer1 =BxPropertyCollection::Get()->GetPhotonEnergyTyvek();
  pointer2 =BxPropertyCollection::Get()->GetTyvekReflectivity();
  myST10 -> AddProperty("REFLECTIVITY",pointer1,pointer2,vector_size);
  
  vector_size=BxPropertyCollection::Get()->GetTyvekEfficiency().size();
  pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyTyvek();
  pointer2 = BxPropertyCollection::Get()->GetTyvekEfficiency();
  myST10 -> AddProperty("EFFICIENCY",pointer1,pointer2,vector_size);
 
   vector_size  = BxPropertyCollection::Get()->GetTyvekspecLOBE().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyTyvek();
   pointer2 = BxPropertyCollection::Get()->GetTyvekspecLOBE();
   myST10->AddProperty("SPECULARLOBECONSTANT", pointer1, pointer2, vector_size);
  
  vector_size  = BxPropertyCollection::Get()->GetTyvekspecSPIKE().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyTyvek();
   pointer2 = BxPropertyCollection::Get()->GetTyvekspecSPIKE();
   myST10->AddProperty("SPECULARSPIKECONSTANT", pointer1, pointer2, vector_size);

  vector_size  = BxPropertyCollection::Get()->GetTyvekbackSCATT().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyTyvek();
   pointer2 = BxPropertyCollection::Get()->GetTyvekbackSCATT();
   myST10->AddProperty("BACKSCATTERCONSTANT", pointer1, pointer2, vector_size);
   

   //
   //---------------------TankSteel surface attached to Tyvek surfaces against TankSteel and water--
     
   fOpTankSteelSurface = new G4OpticalSurface("OpTankSteelSurface");
     fOpTankSteelSurface->SetType(dielectric_metal);
     fOpTankSteelSurface->SetModel(unified);
     fOpTankSteelSurface->SetFinish(polished); 
  
   fTankSteelSurface =   new G4LogicalBorderSurface("TankSteelSurface",
				       			fphysicTankTyvek,
							fphysicTankSteel,
				       			fOpTankSteelSurface);
    
    
   fOpTankTyvekSurface = new G4OpticalSurface("OpTankTyvekSurface");
     fOpTankTyvekSurface->SetType(dielectric_dielectric);
     fOpTankTyvekSurface->SetModel(unified);
     fOpTankTyvekSurface->SetFinish(polished); 
  
   fTankTyvekSurface =   new G4LogicalBorderSurface("TankTyvekSurface",
				       			fphysicTankTyvek,
							fphysicTank,
							fOpTankTyvekSurface);
   
   
   //--------------------Water Tank surfaces------------------------------------------------
  

   fOpTankOutSurface = new G4OpticalSurface("OpTankOutSurface");
   fOpTankOutSurface->SetType(dielectric_dielectric);
   fOpTankOutSurface->SetModel(glisur);
   fOpTankOutSurface->SetFinish(groundbackpainted);

   fOpTankOutSurface->SetMaterialPropertiesTable(myST10);
   
   fTankOutSurface = new G4LogicalBorderSurface("TankOutSurface",
				  fphysicTank,
				  fphysicTankTyvek,
				  fOpTankOutSurface);
   
   fOpTankInSurface = new G4OpticalSurface("OpTankInSurface");

   fOpTankInSurface->SetType(dielectric_dielectric);
   fOpTankInSurface->SetModel(glisur);
   fOpTankInSurface->SetFinish(groundbackpainted);
   
   fOpTankInSurface->SetMaterialPropertiesTable(myST10);
   
   fTankInSurface = new G4LogicalBorderSurface("TankInSurface",
						  fphysicTank,
						  fphysicAddTyvek,
						  fOpTankInSurface);


   //----Inner Sphere Surface of OUTER Detector "AddTyvek"----------------------------------
   //
  
   fOpAddTyvekOutSurface = new G4OpticalSurface("OpAddTyvekOutSurface");
   fOpAddTyvekOutSurface->SetType(dielectric_dielectric);
   fOpAddTyvekOutSurface->SetModel(unified);
   fOpAddTyvekOutSurface->SetFinish(polished);

   fAddTyvekOutSurface = new G4LogicalBorderSurface("AddTyvekOutSurface",
				  fphysicAddTyvek,
				  fphysicTank,
				  fOpAddTyvekOutSurface);
   

   fOpAddTyvekInSurface = new G4OpticalSurface("OpAddTyvekInSurface");
     fOpAddTyvekInSurface->SetType(dielectric_metal);
     fOpAddTyvekInSurface->SetModel(unified);
     fOpAddTyvekInSurface->SetFinish(polished); 
  
   fAddTyvekInSurface =   new G4LogicalBorderSurface("AddTyvekInSurface",
				       			fphysicAddTyvek,
							fphysicSSS,
				       			fOpAddTyvekInSurface);


//------Internal ring surface against Tyvek of AddTyvek sphere+external ring------

    fOpInternalPlatformSurface = new G4OpticalSurface("OpInternalPlatformSurface");
     fOpTankSteelSurface->SetType(dielectric_metal);
     fOpTankSteelSurface->SetModel(unified);
     fOpTankSteelSurface->SetFinish(polished); 
  
    fInternalPlatformSurface =   new G4LogicalBorderSurface("InternalPlatformSurface",
				       			fphysicAddTyvek,
							fphysicInternalPlatformTank,
				       			fOpInternalPlatformSurface);

   //-----------SSS Legs' surfaces-------------------------------------------------- 

    fOpLegSurface = new G4OpticalSurface("OpLegSurface");
     fOpTankSteelSurface->SetType(dielectric_metal);
     fOpTankSteelSurface->SetModel(unified);
     fOpTankSteelSurface->SetFinish(polished); 
  
    fLegSurface=new G4LogicalBorderSurface*[fNlegs];

    fOpLegSurface -> SetMaterialPropertiesTable(OpSSSSurfaceProperty);
    
    for (G4int i=0; i<fNlegs; i++){

    fLegSurface[i] =   new G4LogicalBorderSurface("LegSurface",
				       			fphysicTank,
							fphysicLeg[i],
				       			fOpLegSurface);
  
    }
    
    }//end IsTank

//--------------Vessels' EndCaps surfaces------------  
  if (IsEndCap){
  
   G4MaterialPropertiesTable* OpNylonProperty = new G4MaterialPropertiesTable();
   theNum  = BxPropertyCollection::Get()->GetNylonReflectivity().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyNylonReflectivity();
   pointer2 = BxPropertyCollection::Get()->GetNylonReflectivity();
   OpNylonProperty->AddProperty("REFLECTIVITY", pointer1,  pointer2, theNum);
	
   theNum  = BxPropertyCollection::Get()->GetNylonEfficiency().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyNylonEff();
   pointer2 = BxPropertyCollection::Get()->GetNylonEfficiency();
   OpNylonProperty->AddProperty("EFFICIENCY", pointer1,  pointer2, theNum);
 

//Nylon Bridge surface
   fOpNylonBridgeNorthInSurface = new G4OpticalSurface("OpNylonBridgeNorthInSurface");
   fOpNylonBridgeNorthInSurface->SetType(dielectric_dielectric);
   fOpNylonBridgeNorthInSurface->SetModel(unified);
   fOpNylonBridgeNorthInSurface->SetFinish(groundfrontpainted);

  fOpNylonBridgeNorthInSurface->SetMaterialPropertiesTable(OpNylonProperty);
  
  fNylonBridgeNorthInSurface = new G4LogicalBorderSurface("NylonBridgeNorthInSurface",
						  fphysicZoneII,
						  fphysicNylonBridgeNorth,
						  fOpNylonBridgeNorthInSurface);


   fOpNylonBridgeNorthOutSurface = new G4OpticalSurface("OpNylonBridgeNorthOutSurface");
   fOpNylonBridgeNorthOutSurface->SetType(dielectric_dielectric);
   fOpNylonBridgeNorthOutSurface->SetModel(unified);
   fOpNylonBridgeNorthOutSurface->SetFinish(groundfrontpainted);
  
   fOpNylonBridgeNorthOutSurface->SetMaterialPropertiesTable(OpNylonProperty);

   fNylonBridgeNorthOutSurface = new G4LogicalBorderSurface("NylonBridgeNorthOutSurface",
						  fphysicNylonBridgeNorth,
						  fphysicZoneII,
						  fOpNylonBridgeNorthOutSurface);

   fOpNylonBridgeSouthInSurface = new G4OpticalSurface("OpNylonBridgeSouthInSurface");
   fOpNylonBridgeSouthInSurface->SetType(dielectric_dielectric);
   fOpNylonBridgeSouthInSurface->SetModel(unified);
   fOpNylonBridgeSouthInSurface->SetFinish(groundfrontpainted);
  
   fOpNylonBridgeSouthInSurface->SetMaterialPropertiesTable(OpNylonProperty);

   fNylonBridgeSouthInSurface = new G4LogicalBorderSurface("NylonBridgeSouthInSurface",
						  fphysicZoneII,
						  fphysicNylonBridgeSouth,
						  fOpNylonBridgeSouthInSurface);


   fOpNylonBridgeSouthOutSurface = new G4OpticalSurface("OpNylonBridgeSouthOutSurface");
   fOpNylonBridgeSouthOutSurface->SetType(dielectric_dielectric);
   fOpNylonBridgeSouthOutSurface->SetModel(unified);
   fOpNylonBridgeSouthOutSurface->SetFinish(groundfrontpainted);
  
   fOpNylonBridgeSouthOutSurface->SetMaterialPropertiesTable(OpNylonProperty);

   fNylonBridgeSouthOutSurface = new G4LogicalBorderSurface("NylonBridgeSouthOutSurface",
						  fphysicNylonBridgeSouth,
						  fphysicZoneII,
						  fOpNylonBridgeSouthOutSurface);

  
	  //IV and OV Nylon Rings + Nylon Tubes
  fOpNylonInSurface = new G4OpticalSurface("OpNylonInSurface");
  fOpNylonInSurface->SetType(dielectric_dielectric);
  fOpNylonInSurface->SetModel(unified);
  fOpNylonInSurface->SetFinish(groundfrontpainted);
  //fOpNylonInSurface->SetFinish(polished);
  
  fOpNylonInSurface->SetMaterialPropertiesTable(OpNylonProperty);
  
  fIVNylonRingNorthInSurface = new G4LogicalBorderSurface("IVNylonRingNorthInSurface",
						fphysicZoneII,
						fphysicIVNylonRingNorth,
						fOpNylonInSurface);
   
  fIVNylonRingSouthInSurface = new G4LogicalBorderSurface("IVNylonRingSouthInSurface",
						fphysicZoneII,
						fphysicIVNylonRingSouth,
						fOpNylonInSurface);
  
  fOVNylonRingNorthInSurface = new G4LogicalBorderSurface("OVNylonRingNorthInSurface",
						fphysicZoneIII,
						fphysicOVNylonRingNorth,
						fOpNylonInSurface);
   
  fOVNylonRingSouthInSurface = new G4LogicalBorderSurface("OVNylonRingSouthInSurface",
						fphysicZoneIII,
						fphysicOVNylonRingSouth,
						fOpNylonInSurface);

  fNylonTubeNorthInSurface = new G4LogicalBorderSurface("NylonTubeNorthInSurface",
						fphysicZoneII,
						fphysicNylonTubeNorth,
						fOpNylonInSurface);
  
  fNylonTubeSouthInSurface = new G4LogicalBorderSurface("NylonTubeSouthInSurface",
						fphysicZoneII,
						fphysicNylonTubeSouth,
						fOpNylonInSurface);
 


  fOpNylonOutSurface = new G4OpticalSurface("OpNylonOutSurface");
  fOpNylonOutSurface->SetType(dielectric_dielectric);
  fOpNylonOutSurface->SetModel(unified);
  fOpNylonOutSurface->SetFinish(groundfrontpainted);
  //fOpNylonOutSurface->SetFinish(polished);
  
  fOpNylonOutSurface->SetMaterialPropertiesTable(OpNylonProperty);
 
  fIVNylonRingNorthOutSurface = new G4LogicalBorderSurface("IVNylonRingNorthOutSurface",
						fphysicIVNylonRingNorth,
						fphysicZoneII,
						fOpNylonOutSurface);
   
  fIVNylonRingSouthOutSurface = new G4LogicalBorderSurface("IVNylonRingSouthOutSurface",
						fphysicIVNylonRingSouth,
						fphysicZoneII,
						fOpNylonOutSurface);
 
  fOVNylonRingNorthOutSurface = new G4LogicalBorderSurface("OVNylonRingNorthOutSurface",
						fphysicOVNylonRingNorth,
						fphysicZoneIII,
						fOpNylonOutSurface);
   
  fOVNylonRingSouthOutSurface = new G4LogicalBorderSurface("OVNylonRingSouthOutSurface",
						fphysicOVNylonRingSouth,
						fphysicZoneIII,
						fOpNylonOutSurface);

  fNylonTubeNorthOutSurface = new G4LogicalBorderSurface("NylonTubeNorthOutSurface",
						fphysicNylonTubeNorth,
						fphysicZoneII,
						fOpNylonOutSurface);
  
  fNylonTubeSouthOutSurface = new G4LogicalBorderSurface("NylonTubeSouthOutSurface",
						fphysicNylonTubeSouth,
						fphysicZoneII,
						fOpNylonOutSurface);
  
  //IV and OV Steel Tubes
  fOpSteelTubeSurface = new G4OpticalSurface("SteelTubeSurface");
  fOpSteelTubeSurface->SetType(dielectric_metal);
  fOpSteelTubeSurface->SetModel(unified);
  fOpSteelTubeSurface->SetFinish(polished);
 
   G4MaterialPropertiesTable* OpTubeSteelProperty = new G4MaterialPropertiesTable();
   theNum  = BxPropertyCollection::Get()->GetTubeSteelReflectivity().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyTubeSteelRef();
   pointer2 = BxPropertyCollection::Get()->GetTubeSteelReflectivity();
   OpTubeSteelProperty->AddProperty("REFLECTIVITY", pointer1,  pointer2, theNum);

   theNum  = BxPropertyCollection::Get()->GetTubeSteelEfficiency().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyTubeSteelEff();
   pointer2 = BxPropertyCollection::Get()->GetTubeSteelEfficiency();
   OpTubeSteelProperty->AddProperty("EFFICIENCY", pointer1,  pointer2, theNum);
 
  fOpSteelTubeSurface->SetMaterialPropertiesTable(OpTubeSteelProperty);
  
  fIVTubeNorthSurface = new G4LogicalBorderSurface("IVTubeNorthSurface",
						fphysicZoneII,
						fphysicIVTubeNorth,
						fOpSteelTubeSurface);
  
  fIVTubeSouthSurface = new G4LogicalBorderSurface("IVTubeSouthSurface",
						fphysicZoneII,
						fphysicIVTubeSouth,
						fOpSteelTubeSurface);
  
  fOVTubeNorthSurface = new G4LogicalBorderSurface("OVTubeNorthSurface",
						fphysicZoneIII,
						fphysicOVTubeNorth,
						fOpSteelTubeSurface);
  
  fOVTubeSouthSurface = new G4LogicalBorderSurface("OVTubeSouthSurface",
						fphysicZoneIII,
						fphysicOVTubeSouth,
						fOpSteelTubeSurface);
 

  //Copper Struts
  fOpCopperStrutSurface = new G4OpticalSurface("CopperStrutSurface");
  fOpCopperStrutSurface->SetType(dielectric_metal);
  fOpCopperStrutSurface->SetModel(unified);
  fOpCopperStrutSurface->SetFinish(polished);
 
   G4MaterialPropertiesTable* OpCopperStrutProperty = new G4MaterialPropertiesTable();
   theNum  = BxPropertyCollection::Get()->GetCopperStrutReflectivity().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyCopperStrutRef();
   pointer2 = BxPropertyCollection::Get()->GetCopperStrutReflectivity();
   OpCopperStrutProperty->AddProperty("REFLECTIVITY", pointer1,  pointer2, theNum);

   theNum  = BxPropertyCollection::Get()->GetCopperStrutEfficiency().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyCopperStrutEff();
   pointer2 = BxPropertyCollection::Get()->GetCopperStrutEfficiency();
   OpCopperStrutProperty->AddProperty("EFFICIENCY", pointer1,  pointer2, theNum);
 
  fOpCopperStrutSurface->SetMaterialPropertiesTable(OpCopperStrutProperty);

  fCopperStrutSurface=new G4LogicalBorderSurface*[fNumOfStruts];
  
  for (G4int i=0; i<fNumOfStruts; i++){

  fCopperStrutSurface[i] = new G4LogicalBorderSurface("CopperStrutSurface",
						fphysicZoneII,
						fphysicCopperStrut[i],
						fOpCopperStrutSurface);
  
  }

  }
   
   //--------------------PMT Concentrator optical surfaces-----------------------
if (fPMTFlag != 1 && !fNoConcentrators){

   fOpExternalLGSurface = new G4OpticalSurface("OpExternalLGSurface");
   fOpExternalLGSurface->SetType(dielectric_metal);
   fOpExternalLGSurface->SetModel(unified);
   fOpExternalLGSurface->SetFinish(ground);
  
   G4MaterialPropertiesTable* OpExternalLGProperty = new G4MaterialPropertiesTable();
   theNum  = BxPropertyCollection::Get()->GetExternalLGReflectivity().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyExternalLGRef();
   pointer2 = BxPropertyCollection::Get()->GetExternalLGReflectivity();
   OpExternalLGProperty->AddProperty("REFLECTIVITY", pointer1,  pointer2, theNum);

   theNum  = BxPropertyCollection::Get()->GetExternalLGEfficiency().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyExternalLGEff();
   pointer2 = BxPropertyCollection::Get()->GetExternalLGEfficiency();
   OpExternalLGProperty->AddProperty("EFFICIENCY", pointer1,  pointer2, theNum);
 
   theNum  = BxPropertyCollection::Get()->GetExternalLGspecLOBE().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyExternalLGRef();
   pointer2 = BxPropertyCollection::Get()->GetExternalLGspecLOBE();
   OpExternalLGProperty->AddProperty("SPECULARLOBECONSTANT", pointer1, pointer2, theNum);
  
   theNum  = BxPropertyCollection::Get()->GetExternalLGspecSPIKE().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyExternalLGRef();
   pointer2 = BxPropertyCollection::Get()->GetExternalLGspecSPIKE();
   OpExternalLGProperty->AddProperty("SPECULARSPIKECONSTANT", pointer1, pointer2, theNum);

   theNum  = BxPropertyCollection::Get()->GetExternalLGbackSCATT().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyExternalLGRef();
   pointer2 = BxPropertyCollection::Get()->GetExternalLGbackSCATT();
   OpExternalLGProperty->AddProperty("BACKSCATTERCONSTANT", pointer1, pointer2, theNum);

   fOpExternalLGSurface->SetMaterialPropertiesTable(OpExternalLGProperty);


   fOpInternalLGSurface = new G4OpticalSurface("OpInternalLGSurface");
   fOpInternalLGSurface->SetType(dielectric_metal);
   fOpInternalLGSurface->SetModel(unified);
   fOpInternalLGSurface->SetFinish(ground);

   G4MaterialPropertiesTable* OpInternalLGProperty = new G4MaterialPropertiesTable();
   theNum  = BxPropertyCollection::Get()->GetInternalLGReflectivity().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyInternalLGRef();
   pointer2 = BxPropertyCollection::Get()->GetInternalLGReflectivity();
   OpInternalLGProperty->AddProperty("REFLECTIVITY", pointer1,  pointer2, theNum);

   theNum  = BxPropertyCollection::Get()->GetInternalLGEfficiency().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyInternalLGEff();
   pointer2 = BxPropertyCollection::Get()->GetInternalLGEfficiency();
   OpInternalLGProperty->AddProperty("EFFICIENCY", pointer1,  pointer2, theNum);
 
   theNum  = BxPropertyCollection::Get()->GetInternalLGspecLOBE().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyInternalLGRef();
   pointer2 = BxPropertyCollection::Get()->GetInternalLGspecLOBE();
   OpInternalLGProperty->AddProperty("SPECULARLOBECONSTANT", pointer1, pointer2, theNum);
  
   theNum  = BxPropertyCollection::Get()->GetInternalLGspecSPIKE().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyInternalLGRef();
   pointer2 = BxPropertyCollection::Get()->GetInternalLGspecSPIKE();
   OpInternalLGProperty->AddProperty("SPECULARSPIKECONSTANT", pointer1, pointer2, theNum);

   theNum  = BxPropertyCollection::Get()->GetInternalLGbackSCATT().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyInternalLGRef();
   pointer2 = BxPropertyCollection::Get()->GetInternalLGbackSCATT();
   OpInternalLGProperty->AddProperty("BACKSCATTERCONSTANT", pointer1, pointer2, theNum);

   fOpInternalLGSurface->SetMaterialPropertiesTable(OpInternalLGProperty);


}
    
if (fPMTFlag == 1 || fPMTFlag == 2){ 
// disk PMT construction

   //--------------------PMT Cathode optical surface---------------------

   fOpPMTSurface = new G4OpticalSurface("OpPMTSurface");
   fOpPMTSurface->SetType(dielectric_metal);
   fOpPMTSurface->SetModel(glisur);
   fOpPMTSurface->SetFinish(polished);

   G4MaterialPropertiesTable* MST4 = new G4MaterialPropertiesTable();
   theNum  = BxPropertyCollection::Get()->GetPMTReflectivity().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTRef();
   pointer2 = BxPropertyCollection::Get()->GetPMTReflectivity();
   MST4->AddProperty("REFLECTIVITY", pointer1, pointer2, theNum);
   theNum  = BxPropertyCollection::Get()->GetPMTEfficiency().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTEff();
   pointer2 = BxPropertyCollection::Get()->GetPMTEfficiency();
   MST4->AddProperty("EFFICIENCY", pointer1, pointer2, theNum);
   fOpPMTSurface->SetMaterialPropertiesTable (MST4);

  
   //--------------------PMT Shield optical surface-----------------------
   fOpPMTShieldSurface = new G4OpticalSurface("OpPMTShieldSurface");
   fOpPMTShieldSurface->SetType(dielectric_metal);
   fOpPMTShieldSurface->SetModel(glisur);
   fOpPMTShieldSurface->SetFinish(polished);

   G4MaterialPropertiesTable* MST6 = new G4MaterialPropertiesTable();
   theNum  = BxPropertyCollection::Get()->GetPMTShieldReflectivity().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTShieldRef();
   pointer2 = BxPropertyCollection::Get()->GetPMTShieldReflectivity();
   MST6->AddProperty("REFLECTIVITY", pointer1,  pointer2, theNum);
   theNum  = BxPropertyCollection::Get()->GetPMTShieldEfficiency().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTShieldEff();
   pointer2 = BxPropertyCollection::Get()->GetPMTShieldEfficiency();
   MST6->AddProperty("EFFICIENCY", pointer1,  pointer2, theNum);

   fOpPMTShieldSurface->SetMaterialPropertiesTable (MST6);

} else if (fPMTFlag == 3){
//simple PMT	
   
//--------------------PMT Cathode optical surface---------------------

   fOpPMTSurface = new G4OpticalSurface("OpPMTSurface");
   fOpPMTSurface->SetType(dielectric_metal);
   fOpPMTSurface->SetModel(glisur);
   fOpPMTSurface->SetFinish(polished);

   G4MaterialPropertiesTable* MST4 = new G4MaterialPropertiesTable();
   theNum  = BxPropertyCollection::Get()->GetPMTReflectivity().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTRef();
   pointer2 = BxPropertyCollection::Get()->GetPMTReflectivity();
   MST4->AddProperty("REFLECTIVITY", pointer1, pointer2, theNum);
   theNum  = BxPropertyCollection::Get()->GetPMTEfficiency().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTEff();
   pointer2 = BxPropertyCollection::Get()->GetPMTEfficiency();
   MST4->AddProperty("EFFICIENCY", pointer1, pointer2, theNum);
   fOpPMTSurface->SetMaterialPropertiesTable (MST4);

  
   //--------------------PMT Shield optical surface-----------------------
   fOpPMTShieldSurface = new G4OpticalSurface("OpPMTShieldSurface");
   fOpPMTShieldSurface->SetType(dielectric_metal);
   fOpPMTShieldSurface->SetModel(glisur);
   fOpPMTShieldSurface->SetFinish(polished);

   G4MaterialPropertiesTable* MST6 = new G4MaterialPropertiesTable();
   theNum  = BxPropertyCollection::Get()->GetPMTShieldReflectivity().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTShieldRef();
   pointer2 = BxPropertyCollection::Get()->GetPMTShieldReflectivity();
   MST6->AddProperty("REFLECTIVITY", pointer1,  pointer2, theNum);
   theNum  = BxPropertyCollection::Get()->GetPMTShieldEfficiency().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTShieldEff();
   pointer2 = BxPropertyCollection::Get()->GetPMTShieldEfficiency();
   MST6->AddProperty("EFFICIENCY", pointer1,  pointer2, theNum);

   fOpPMTShieldSurface->SetMaterialPropertiesTable (MST6);

  
   //--------------------PMT Housing optical surface-----------------------

   fOpPMTHousingSurface = new G4OpticalSurface("OpPMTHousingSurface");
   fOpPMTHousingSurface->SetType(dielectric_metal);
   fOpPMTHousingSurface->SetModel(glisur);
   fOpPMTHousingSurface->SetFinish(polished);

   G4MaterialPropertiesTable* MST7 = new G4MaterialPropertiesTable();
   theNum  = BxPropertyCollection::Get()->GetPMTHousingReflectivity().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTHousingRef();
   pointer2 = BxPropertyCollection::Get()->GetPMTHousingReflectivity();
   MST7->AddProperty("REFLECTIVITY", pointer1,  pointer2, theNum);
   theNum  = BxPropertyCollection::Get()->GetPMTHousingEfficiency().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTHousingEff();
   pointer2 = BxPropertyCollection::Get()->GetPMTHousingEfficiency();
   MST7->AddProperty("EFFICIENCY", pointer1,  pointer2, theNum);

   fOpPMTHousingSurface->SetMaterialPropertiesTable (MST7);

   //-----------------PMT ring optical surface-----------------------------
  
   fOpPMTringSurface = new G4OpticalSurface("OpPMTringSurface");
   fOpPMTringSurface->SetType(dielectric_metal);
   fOpPMTringSurface->SetModel(unified);
   fOpPMTringSurface->SetFinish(ground);

   G4MaterialPropertiesTable* MSTPMTring = new G4MaterialPropertiesTable();
   theNum  = BxPropertyCollection::Get()->GetPMTringReflectivity().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTringRef();
   pointer2 = BxPropertyCollection::Get()->GetPMTringReflectivity();
   MSTPMTring->AddProperty("REFLECTIVITY", pointer1,  pointer2, theNum);
   
   theNum  = BxPropertyCollection::Get()->GetPMTringEfficiency().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTringEff();
   pointer2 = BxPropertyCollection::Get()->GetPMTringEfficiency();
   MSTPMTring->AddProperty("EFFICIENCY", pointer1,  pointer2, theNum);
   
   theNum  = BxPropertyCollection::Get()->GetPMTringspecLOBE().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTringspecLOBE();
   pointer2 = BxPropertyCollection::Get()->GetPMTringspecLOBE();
   MSTPMTring->AddProperty("SPECULARLOBECONSTANT", pointer1, pointer2, theNum);
  
   theNum  = BxPropertyCollection::Get()->GetPMTringspecSPIKE().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTringspecSPIKE();
   pointer2 = BxPropertyCollection::Get()->GetPMTringspecSPIKE();
   MSTPMTring->AddProperty("SPECULARSPIKECONSTANT", pointer1, pointer2, theNum);

   theNum  = BxPropertyCollection::Get()->GetPMTringbackSCATT().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTringbackSCATT();
   pointer2 = BxPropertyCollection::Get()->GetPMTringbackSCATT();
   MSTPMTring->AddProperty("BACKSCATTERCONSTANT", pointer1, pointer2, theNum);

   fOpPMTringSurface->SetMaterialPropertiesTable (MSTPMTring);


}

if (!fNoConcentrators){
	if (fPMTFlag == 2){
		if (fPMTDistributionFlag==1){
			fExternalLGBorderSurface      = new G4LogicalBorderSurface*[fNumberOfPMT-fNumberOfVetoPMT];
			fInternalLGBorderSurface      = new G4LogicalBorderSurface*[fNumberOfPMT-fNumberOfVetoPMT];
		}else if (fPMTDistributionFlag ==0){
			fExternalLGBorderSurface      = new G4LogicalBorderSurface*[fNumberOfPMT];
			fInternalLGBorderSurface      = new G4LogicalBorderSurface*[fNumberOfPMT];
		}
	} else if(fPMTFlag == 3){
		fExternalLGBorderSurface      = new G4LogicalBorderSurface*[fNumberOfPMT-fNumberOfVetoPMT];
		fInternalLGBorderSurface      = new G4LogicalBorderSurface*[fNumberOfPMT-fNumberOfVetoPMT];
	}
}

if (fPMTFlag == 1 || fPMTFlag == 2){ 
//simple disk PMT construction
   fPMTSurface            = new G4LogicalBorderSurface*[fNumberOfPMT];
   fPMTShieldSurface      = new G4LogicalBorderSurface*[fNumberOfPMT];
} else if (fPMTFlag == 3){

   fPMTSurface            = new G4LogicalBorderSurface*[fNumberOfPMT];
   fPMTShieldSurface      = new G4LogicalBorderSurface*[fNumberOfPMT];
   fPMTHousingSurface      = new G4LogicalBorderSurface*[fNumberOfPMT];
if (!fNoConcentrators)
		fPMTringBorderSurface	      = new G4LogicalBorderSurface*[fNumberOfVetoPMT];
else		
		fPMTringBorderSurface	      = new G4LogicalBorderSurface*[fNumberOfPMT];

}
   //
   //--------------------PMTmu Cathode optical surface---------------------

   fOpPMTmuSurface = new G4OpticalSurface("OpPMTmuSurface");
   fOpPMTmuSurface->SetType(dielectric_metal);
   fOpPMTmuSurface->SetModel(glisur);
   fOpPMTmuSurface->SetFinish(polished);

   G4MaterialPropertiesTable* MST15 = new G4MaterialPropertiesTable();
   theNum  = BxPropertyCollection::Get()->GetPMTReflectivity().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTRef();
   pointer2 = BxPropertyCollection::Get()->GetPMTReflectivity();
   MST15->AddProperty("REFLECTIVITY", pointer1, pointer2, theNum);
   theNum  = BxPropertyCollection::Get()->GetPMTEfficiency().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTEff();
   pointer2 = BxPropertyCollection::Get()->GetPMTEfficiency();
   MST15->AddProperty("EFFICIENCY", pointer1, pointer2, theNum);
   fOpPMTmuSurface->SetMaterialPropertiesTable (MST15);


   //
   //--------------------PMTmu Flange optical surface-----------------------

   fOpPMTmuFlangeInSurface = new G4OpticalSurface("OpPMTmuFlangeInSurface");
   fOpPMTmuFlangeInSurface->SetType(dielectric_metal);
   fOpPMTmuFlangeInSurface->SetModel(glisur);
   fOpPMTmuFlangeInSurface->SetFinish(polished);

   G4MaterialPropertiesTable* MST16 = new G4MaterialPropertiesTable();
   theNum  = BxPropertyCollection::Get()->GetPMTShieldReflectivity().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTShieldRef();
   pointer2 = BxPropertyCollection::Get()->GetPMTShieldReflectivity();
   MST16->AddProperty("REFLECTIVITY", pointer1,  pointer2, theNum);
   theNum  = BxPropertyCollection::Get()->GetPMTShieldEfficiency().size();
   pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyPMTShieldEff();
   pointer2 = BxPropertyCollection::Get()->GetPMTShieldEfficiency();
   MST16->AddProperty("EFFICIENCY", pointer1,  pointer2, theNum);

   fOpPMTmuFlangeInSurface->SetMaterialPropertiesTable (MST16);

   //--------------------PMTmu Shield optical surface-----------------------

   fOpPMTmuShieldInSurface = new G4OpticalSurface("OpPMTmuShieldInSurface");
   fOpPMTmuShieldInSurface->SetType(dielectric_metal);
   fOpPMTmuShieldInSurface->SetModel(glisur);
   fOpPMTmuShieldInSurface->SetFinish(polished);

   fOpPMTmuShieldInSurface->SetMaterialPropertiesTable (MST16);

   
 // allocation of the logical border surfaces
   fPMTmuSurface = new G4LogicalBorderSurface*[fNumberOfPMTmu];
   fPMTmuFlangeInSurface = new G4LogicalBorderSurface*[fNumberOfPMTmu];
   fPMTmuShieldInSurface = new G4LogicalBorderSurface*[fNumberOfPMTmuGround];

   
   //________Placement of Inner Detector Optical Surfaces


   for(G4int i=0;i<fNumberOfPMT;i++)	{


if (fPMTFlag == 1 || fPMTFlag == 2){ 
//simple disk PMT construction

   //--------------------PMT Cathode surface-----------------------
  //Attention - the sequence is crucial for the reflection
  //fphysicZoneIII,
  //fphysicPMT[i]
  // and not vice versa !!!

   	fPMTSurface[i] = new G4LogicalBorderSurface("PMTSurface",
						fphysicZoneIII,
						fphysicPMT[i],
						fOpPMTSurface);

   //
   //--------------------PMT Shield surface-----------------------

        fPMTShieldSurface[i] = new G4LogicalBorderSurface("PMTShieldSurface",
						fphysicZoneIII,
						fphysicPMTShield[i],
						fOpPMTShieldSurface);


} else if (fPMTFlag == 3){
   //--------------------PMT Cathode surface-----------------------

   	fPMTSurface[i] = new G4LogicalBorderSurface("PMTSurface",
						fphysicZoneIII,
						fphysicPMT[i],
						fOpPMTSurface);

   //--------------------PMT Shield surface-----------------------

        fPMTShieldSurface[i] = new G4LogicalBorderSurface("PMTShieldSurface",
						fphysicZoneIII,
						fphysicPMTShield[i],
						fOpPMTShieldSurface);

	
   //--------------------PMT Housing surface-----------------------


   	fPMTHousingSurface[i] = new G4LogicalBorderSurface("PMTHousingSurface",
						fphysicZoneIII,
						fphysicPMTHousing[i],
						fOpPMTHousingSurface);

	
	}

}

//
//-------------------Light Guide surfaces--------------------------

if (fPMTFlag > 1 && !fNoConcentrators){ 
  k_conc=0;
  k_noconc=0;
  for(G4int i=0;i<fNumberOfPMT;i++)	{
    if(!fVetoPmt.at(i))	{
        fExternalLGBorderSurface[k_conc] = new G4LogicalBorderSurface("ExternalLGBorderSurface",
					      fphysicZoneIII,
					      fphysicExternalLightGuide[k_conc],
					      fOpExternalLGSurface);
        fInternalLGBorderSurface[k_conc] = new G4LogicalBorderSurface("InternalLGBorderSurface",
					      fphysicZoneIII,
					      fphysicInternalLightGuide[k_conc],
					      fOpInternalLGSurface);
      k_conc++;      
    } else if (fPMTFlag==3 && fVetoPmt.at(i)){
        fPMTringBorderSurface[k_noconc] = new G4LogicalBorderSurface("PMTringBorderSurface",
					      fphysicZoneIII,
					      fphysicPMTring[k_noconc],
					      fOpPMTringSurface);
    k_noconc++;
    
  } else if (fPMTFlag==3 && fNoConcentrators) {

        fPMTringBorderSurface[i] = new G4LogicalBorderSurface("PMTringBorderSurface",
					      fphysicZoneIII,
					      fphysicPMTring[i],
					      fOpPMTringSurface);
		}
	}
}

  
   //
   //________Placement of OUTER Detector Optical Surfaces
   //
if (IsTank && !noPMTs){
G4int counter=0;
for(G4int i=0;i<fNumberOfPMTmu;i++)	{


	if(fVetoOD_Pmt[i] == 0)	{

          fPMTmuSurface[i] = new G4LogicalBorderSurface("PMTmuSurface",
						fphysicTank,
						fphysicPMTmu[i],
						fOpPMTmuSurface);

   	fPMTmuFlangeInSurface[i] = new
			G4LogicalBorderSurface("PMTFlangeInSurface",
						fphysicTank,
						fphysicPMTmuFlange[i],
						fOpPMTmuFlangeInSurface);

	  }	else	{


                      	fPMTmuSurface[i] = new G4LogicalBorderSurface("PMTmuSurface",
						fphysicTank,
						fphysicPMTmu[i],
						fOpPMTmuSurface);


   	fPMTmuFlangeInSurface[i] = new
			G4LogicalBorderSurface("PMTmuFlangeInSurface",
						fphysicTank,
						fphysicPMTmuFlange[i],
						fOpPMTmuFlangeInSurface);

   	fPMTmuShieldInSurface[i] = new
			G4LogicalBorderSurface("PMTmuShieldInSurface",
						fphysicTank,
						fphysicPMTmuShield[counter],
						fOpPMTmuShieldInSurface);

	counter++;
	}
}
   }
                                         
//   if(IsOpera) DefineOpera();

   if(fSource) AddSource();

   PrintDetectorParameters();
   SetVisAttributes();
   return fphysicWorld;
}

void BxDetectorConstruction::SetVisAttributes(){
   // Visualization attributes

	//World
   //flogicWorld->SetVisAttributes (G4VisAttributes::Invisible);

	//Water Tank
if (IsTank){
   
   G4VisAttributes* UpPlateSteelVisAtt= new G4VisAttributes(G4Colour(192/255.,192/255.,192/255.));
   UpPlateSteelVisAtt->SetVisibility(true);
   UpPlateSteelVisAtt->SetForceWireframe(false);
   UpPlateSteelVisAtt->SetForceSolid(true);
   flogicUpPlateSteel->SetVisAttributes(UpPlateSteelVisAtt);
   
   G4VisAttributes* DownPlateSteelVisAtt= new G4VisAttributes(G4Colour(192/255.,192/255.,192/255.));
   DownPlateSteelVisAtt->SetVisibility(true);
   DownPlateSteelVisAtt->SetForceWireframe(false);
   DownPlateSteelVisAtt->SetForceSolid(true);
   flogicDownPlateSteel->SetVisAttributes(DownPlateSteelVisAtt);
   
   G4VisAttributes* TankSteelVisAtt= new G4VisAttributes(G4Colour(192/255.,192/255.,192/255.));
   TankSteelVisAtt->SetVisibility(true);
   TankSteelVisAtt->SetForceWireframe(false);
   TankSteelVisAtt->SetForceSolid(true);
   flogicTankSteel->SetVisAttributes(TankSteelVisAtt);
   
   G4VisAttributes* TankTyvekVisAtt= new G4VisAttributes(G4Colour(1.,1.,1.));
   TankTyvekVisAtt->SetVisibility(true);
   TankTyvekVisAtt->SetForceWireframe(false);
   TankTyvekVisAtt->SetForceSolid(true);
   flogicTankTyvek->SetVisAttributes(TankTyvekVisAtt);
   
   G4VisAttributes* TankVisAtt= new G4VisAttributes(G4Colour(99./255.,51./255.,208./255.));
   TankVisAtt->SetVisibility(false);
   TankVisAtt->SetForceWireframe(false);
   TankVisAtt->SetForceSolid(true);
   flogicTank->SetVisAttributes(TankVisAtt);
   
   G4VisAttributes* AddTyvekVisAtt= new G4VisAttributes(G4Colour(1.,1.,1.));
   AddTyvekVisAtt->SetVisibility(true);
   AddTyvekVisAtt->SetForceWireframe(false);
   AddTyvekVisAtt->SetForceSolid(true);
   flogicAddTyvek->SetVisAttributes(AddTyvekVisAtt);
   
   G4VisAttributes* InternalPlatformVisAtt= new G4VisAttributes(G4Colour(199./255.,190./255.,190/255.));
   InternalPlatformVisAtt->SetVisibility(false);
   InternalPlatformVisAtt->SetForceWireframe(false);
   InternalPlatformVisAtt->SetForceSolid(true);
   flogicInternalPlatformTank->SetVisAttributes(InternalPlatformVisAtt);

   G4VisAttributes* LegsVisAtt= new G4VisAttributes(G4Colour(192./255.,192./255.,192./255.));
   LegsVisAtt->SetVisibility(true);
   LegsVisAtt->SetForceWireframe(false);
   LegsVisAtt->SetForceSolid(true);
   flogicLeg->SetVisAttributes(LegsVisAtt);
   flogicLeg8->SetVisAttributes(LegsVisAtt);



}
 	//SSS
   G4VisAttributes* SSSVisAtt= new G4VisAttributes(G4Colour(192./255,192./255.,192./255.));
   SSSVisAtt->SetVisibility(true);
   SSSVisAtt->SetForceWireframe(false);
   SSSVisAtt->SetForceSolid(true);
   flogicSSS->SetVisAttributes(SSSVisAtt);

   	//ZoneIII
   G4VisAttributes* ZoneIIIVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
   ZoneIIIVisAtt->SetVisibility(false);
   ZoneIIIVisAtt->SetForceWireframe(true);
   ZoneIIIVisAtt->SetForceSolid(false);
   flogicZoneIII->SetVisAttributes(ZoneIIIVisAtt);

   //Shroud 
   G4VisAttributes* ShroudVisAtt= new G4VisAttributes(G4Colour(171./255.,205./255.,239./255.));
   ShroudVisAtt->SetVisibility(true);
   ShroudVisAtt->SetForceWireframe(false);
   ShroudVisAtt->SetForceSolid(true);
   flogicShroud->SetVisAttributes(ShroudVisAtt);

   //ZoneII
   G4VisAttributes* ZoneIIVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
   ZoneIIVisAtt->SetVisibility(false);
   ZoneIIVisAtt->SetForceWireframe(true);
   ZoneIIVisAtt->SetForceSolid(false);
   flogicZoneII->SetVisAttributes(ZoneIIVisAtt);
	
   //Vessel   
   G4VisAttributes* VesselVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0));
   VesselVisAtt->SetVisibility(true);
   VesselVisAtt->SetForceWireframe(false);
   VesselVisAtt->SetForceSolid(true);
   flogicVessel->SetVisAttributes(VesselVisAtt);
  
   //ZoneI
   G4VisAttributes* ZoneIVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
   ZoneIVisAtt->SetVisibility(false);
   ZoneIVisAtt->SetForceWireframe(true);
   ZoneIVisAtt->SetForceSolid(false);
   flogicZoneI->SetVisAttributes(ZoneIVisAtt);
   
   //PMT

   G4VisAttributes* VisPMT= new G4VisAttributes(G4Colour(1.,0.,0.));
   VisPMT->SetVisibility(true);
   VisPMT->SetForceSolid(true);
   //VisPMT->SetForceWireframe(false);
   flogicPMT->SetVisAttributes(VisPMT);

   //PMT Shield

   G4VisAttributes* PMTShieldVisAtt = new G4VisAttributes(G4Colour(238./255.,201./255.,0.));
   PMTShieldVisAtt->SetVisibility(true);
   PMTShieldVisAtt->SetForceSolid(true);
   //PMTShieldVisAtt->SetForceWireframe(false);
   flogicPMTShield->SetVisAttributes(PMTShieldVisAtt);
 
   //PMT Housing
   if (fPMTFlag>2){
   G4VisAttributes* PMTHousingVisAtt = new G4VisAttributes(G4Colour(3./255.,192./255.,60./255.));
   PMTHousingVisAtt->SetVisibility(true);
   PMTHousingVisAtt->SetForceSolid(true);
   //PMTHousingVisAtt->SetForceWireframe(false);
   flogicPMTHousing->SetVisAttributes(PMTHousingVisAtt);
   }

   if (fPMTFlag>1 && !fNoConcentrators){
   //PMT Concentrator
   G4VisAttributes* ExternalLightGuideVisAtt = new G4VisAttributes(G4Colour(205./255.,201./255.,201./255.));
   ExternalLightGuideVisAtt->SetVisibility(true);
   ExternalLightGuideVisAtt->SetForceSolid(true);
   flogicExternalLightGuide->SetVisAttributes(ExternalLightGuideVisAtt);
   
   G4VisAttributes* InternalLightGuideVisAtt = new G4VisAttributes(G4Colour(238./255.,233./255.,233./255.));
   InternalLightGuideVisAtt->SetVisibility(true);
   InternalLightGuideVisAtt->SetForceSolid(true);
   flogicInternalLightGuide->SetVisAttributes(InternalLightGuideVisAtt);
   }
   
   if (fPMTFlag==3){
	//PMT rings
   G4VisAttributes* PMTringVisAtt = new G4VisAttributes(G4Colour(238./255.,233./255.,233./255.));
   PMTringVisAtt->SetVisibility(true);
   PMTringVisAtt->SetForceSolid(true);
   flogicPMTring->SetVisAttributes(PMTringVisAtt);

   }

   //OD PMTs
   if (IsTank){

   G4VisAttributes* PmtmuVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
   PmtmuVisAtt->SetVisibility(true);
   PmtmuVisAtt->SetForceWireframe(false);
   PmtmuVisAtt->SetForceSolid(true);
   flogicPMTmu->SetVisAttributes(PmtmuVisAtt);

   G4VisAttributes* PMTmuFlangeVisAtt= new G4VisAttributes(G4Colour(0.,1.,1.));
   PMTmuFlangeVisAtt->SetVisibility(true);
   PMTmuFlangeVisAtt->SetForceWireframe(false);
   PMTmuFlangeVisAtt->SetForceSolid(true);
   flogicPMTmuFlange->SetVisAttributes(PMTmuFlangeVisAtt);
   
   G4VisAttributes* PMTmuShieldVisAtt= new G4VisAttributes(G4Colour(1.,1.,0.));
   PMTmuShieldVisAtt->SetVisibility(true);
   PMTmuShieldVisAtt->SetForceWireframe(false);
   PMTmuShieldVisAtt->SetForceSolid(true);
   flogicPMTmuShield->SetVisAttributes(PMTmuShieldVisAtt);
   }

if (IsEndCap){
   G4VisAttributes* NylonBridgeVisAtt= new G4VisAttributes(G4Colour(240./255.,248./255.,1.));
   NylonBridgeVisAtt->SetVisibility(true);
   NylonBridgeVisAtt->SetForceWireframe(false);
   NylonBridgeVisAtt->SetForceSolid(true);
   flogicNylonBridge->SetVisAttributes(NylonBridgeVisAtt);

   G4VisAttributes* NylonTubeVisAtt= new G4VisAttributes(G4Colour(240./255.,248./255.,1.));
   NylonTubeVisAtt->SetVisibility(true);
   NylonTubeVisAtt->SetForceWireframe(false);
   NylonTubeVisAtt->SetForceSolid(true);
   flogicNylonTube->SetVisAttributes(NylonTubeVisAtt);

   G4VisAttributes* IVNylonRingVisAtt= new G4VisAttributes(G4Colour(240./255.,248./255.,1.));
   IVNylonRingVisAtt->SetVisibility(true);
   IVNylonRingVisAtt->SetForceWireframe(false);
   IVNylonRingVisAtt->SetForceSolid(true);
   flogicIVNylonRing->SetVisAttributes(IVNylonRingVisAtt);
   
   G4VisAttributes* OVNylonRingVisAtt= new G4VisAttributes(G4Colour(240./255.,248./255.,1.));
   OVNylonRingVisAtt->SetVisibility(true);
   OVNylonRingVisAtt->SetForceWireframe(false);
   OVNylonRingVisAtt->SetForceSolid(true);
   flogicOVNylonRing->SetVisAttributes(OVNylonRingVisAtt);
   
   G4VisAttributes* IVTubeVisAtt= new G4VisAttributes(G4Colour(158./255,162./255.,166./255.));
   IVTubeVisAtt->SetVisibility(true);
   IVTubeVisAtt->SetForceWireframe(false);
   IVTubeVisAtt->SetForceSolid(true);
   flogicIVTube_North->SetVisAttributes(IVTubeVisAtt);
   flogicIVTube_South->SetVisAttributes(IVTubeVisAtt);

   G4VisAttributes* OVTubeVisAtt= new G4VisAttributes(G4Colour(192./255,192./255.,192./255.));
   OVTubeVisAtt->SetVisibility(true);
   OVTubeVisAtt->SetForceWireframe(false);
   OVTubeVisAtt->SetForceSolid(true);
   flogicOVTube->SetVisAttributes(OVTubeVisAtt);
   
   G4VisAttributes* CopperStrutVisAtt= new G4VisAttributes(G4Colour(217./255,144./255.,88./255.));
   CopperStrutVisAtt->SetVisibility(true);
   CopperStrutVisAtt->SetForceWireframe(false);
   CopperStrutVisAtt->SetForceSolid(true);
   flogicCopperStrut->SetVisAttributes(CopperStrutVisAtt);
}

if (fSource) {

if (fSource==1 || fSource==2 || fSource==3){
   
   G4VisAttributes* SourceVialVisAtt= new G4VisAttributes(G4Colour(1.,0.,0.));
   SourceVialVisAtt->SetVisibility(true);
   SourceVialVisAtt->SetForceWireframe(false);
   SourceVialVisAtt->SetForceSolid(true);
   flogicSourceVial->SetVisAttributes(SourceVialVisAtt);


   G4VisAttributes* SourceVisAtt= new G4VisAttributes(G4Colour(1.,1.,0.));
   SourceVisAtt->SetVisibility(true);
   SourceVisAtt->SetForceWireframe(false);
   SourceVisAtt->SetForceSolid(true);
   flogicSource->SetVisAttributes(SourceVisAtt);

} else {

   G4VisAttributes* DerlinVisAtt= new G4VisAttributes(G4Colour(1.,1.,1.));
   DerlinVisAtt->SetVisibility(true);
   DerlinVisAtt->SetForceWireframe(false);
   DerlinVisAtt->SetForceSolid(true);
   flogicDerlin->SetVisAttributes(DerlinVisAtt);

   G4VisAttributes* ScrewVisAtt= new G4VisAttributes(G4Colour(190./255.,190./255.,190./255.));
   ScrewVisAtt->SetVisibility(true);
   ScrewVisAtt->SetForceWireframe(false);
   ScrewVisAtt->SetForceSolid(true);
   flogicScrew->SetVisAttributes(ScrewVisAtt);

   G4VisAttributes* PbContainerVisAtt= new G4VisAttributes(G4Colour(190./255.,190./255.,190./255.));
   PbContainerVisAtt->SetVisibility(true);
   PbContainerVisAtt->SetForceWireframe(false);
   PbContainerVisAtt->SetForceSolid(true);
   flogicPb->SetVisAttributes(PbContainerVisAtt);

   G4VisAttributes* AmBeVisAtt= new G4VisAttributes(G4Colour(0.,1.,0.));
   AmBeVisAtt->SetVisibility(true);
   AmBeVisAtt->SetForceWireframe(false);
   AmBeVisAtt->SetForceSolid(true);
   flogicAmBe->SetVisAttributes(AmBeVisAtt);
}

	}

}


    void BxDetectorConstruction::PrintDetectorParameters() {
    BxLog(trace) << "SSS reflectivity       = " <<
		BxPropertyCollection::Get()->GetSSSReflectivity()[0] << endlog;
    BxLog(trace) << "Extern. LG reflectivity = " <<
		BxPropertyCollection::Get()->GetExternalLGReflectivity()[0] << endlog;
    BxLog(trace) << "Intern. LG reflectivity = " <<
		BxPropertyCollection::Get()->GetInternalLGReflectivity()[0] << endlog;
   BxLog(trace) << "PMT Cathode reflectivity  = " <<
		BxPropertyCollection::Get()->GetPMTReflectivity()[0] << endlog;
   BxLog(trace) << "PMT Housing reflectivity  = " <<
		BxPropertyCollection::Get()->GetPMTHousingReflectivity()[0] << endlog;
   BxLog(trace) << "PMT Shield reflectivity  = " <<
		BxPropertyCollection::Get()->GetPMTShieldReflectivity()[0] << endlog;
}

void BxDetectorConstruction::UpdateGeometry() {
   G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());
}
/*
void BxDetectorConstruction::DefineOpera() {

  fOperaZoneSizeX =                 BxReadParameters::Get()->GetOperaZoneSizeX();
  fOperaZoneSizeY =                 BxReadParameters::Get()->GetOperaZoneSizeY();
  fOperaZoneSizeZ =                 BxReadParameters::Get()->GetOperaZoneSizeZ();

  fOperaSMSizeX =                   BxReadParameters::Get()->GetOperaSMSizeX();
  fOperaSMSizeY =                   BxReadParameters::Get()->GetOperaSMSizeY();
  fOperaSMSizeZ =                   BxReadParameters::Get()->GetOperaSMSizeZ();
 
  fOperaTTSizeX =                   BxReadParameters::Get()->GetOperaTTSizeX();
  fOperaTTSizeY =                   BxReadParameters::Get()->GetOperaTTSizeY();
  fOperaTTSizeZ =                   BxReadParameters::Get()->GetOperaTTSizeZ();

  fOperaMSSizeX =                   BxReadParameters::Get()->GetOperaMSSizeX();
  fOperaMSSizeY =                   BxReadParameters::Get()->GetOperaMSSizeY();
  fOperaMSSizeZ =                   BxReadParameters::Get()->GetOperaMSSizeZ();

  fOperaTTWallsSpacing =            BxReadParameters::Get()->GetOperaTTWallsSpacing();				     
  fOperaMSWallsSpacing =            BxReadParameters::Get()->GetOperaMSWallsSpacing();
  fOperaTTWallsFirstPos =           BxReadParameters::Get()->GetOperaTTWallsFirstPos();				     
  fOperaMSWallsFirstPos =           BxReadParameters::Get()->GetOperaMSWallsFirstPos();    

  fOperaZonePos =                   BxReadParameters::Get()->GetOperaZonePos();
  fOperaSM1Pos =                    BxReadParameters::Get()->GetOperaSM1Pos();
  fOperaSM2Pos =                    BxReadParameters::Get()->GetOperaSM2Pos();

  fOperaMSMaterial 		= BxMaterial::Get()->GetWater();

    // Opera zone 
  BxLog(trace) << "WLD: " << fWorldSizeX/m << "m x " 
	       << fWorldSizeY/m << "m x " << fWorldSizeZ/m << "m" << endlog;
  BxLog(trace) << "OPR: " << fOperaZoneSizeX/m << "m x " 
	       << fOperaZoneSizeY/m << "m x " << fOperaZoneSizeZ/m << "m" << endlog;

  fsolidOperaZone  = new G4Box("OperaZone_solid", fOperaZoneSizeX/2.0, fOperaZoneSizeY/2.0, fOperaZoneSizeZ/2.0);

  flogicOperaZone  = new G4LogicalVolume(fsolidOperaZone, fOperaZoneMaterial, "OperaZone_logic");

  fphysicOperaZone = new G4PVPlacement(0,
				       fOperaZonePos,
				       "OperaZone",
				       flogicOperaZone,
				       fphysicWorld,
				       false,
				       0, IsOverlap);

  // Super-Modules' definition
  fsolidOperaSM    = new G4Box("OperaSM_solid", fOperaSMSizeX/2.0, fOperaSMSizeY/2.0, fOperaSMSizeZ/2.0);

  flogicOperaSM1   = new G4LogicalVolume(fsolidOperaSM, fOperaSMMaterial, "OperaSM1_logic");

  fphysicOperaSM1  = new G4PVPlacement(0,
				       fOperaSM1Pos,
				       "OperaSM1",
				       flogicOperaSM1,
				       fphysicOperaZone,
				       false,
				       0, IsOverlap);

  flogicOperaSM2   = new G4LogicalVolume(fsolidOperaSM, fOperaSMMaterial, "OperaSM2_logic");

  fphysicOperaSM2  = new G4PVPlacement(0,
				       fOperaSM2Pos,
				       "OperaSM2",
				       flogicOperaSM2,
				       fphysicOperaZone,
				       false,
				       0, IsOverlap);

  // Target Tracker volumes
  fsolidOperaTT    = new G4Box("OperaTT_solid", fOperaTTSizeX/2.0, fOperaTTSizeY/2.0, fOperaTTSizeZ/2.0);

  flogicOperaTT    = new G4LogicalVolume(fsolidOperaTT, fOperaTTMaterial, "OperaTT_logic");

  G4VPVParameterisation* fOperaTTParam = new BxChamberParameterization(31,
								       fOperaTTWallsFirstPos,
								       fOperaTTWallsSpacing,
								       fOperaTTSizeX,
								       fOperaTTSizeY,
								       fOperaTTSizeY);
  fphysicOperaTT1  = new G4PVParameterised("OperaTT1",
					   flogicOperaTT,
					   flogicOperaSM1,
					   kXAxis,
					   31,
					   fOperaTTParam);

  fphysicOperaTT2  = new G4PVParameterised("OperaTT2",
					   flogicOperaTT,
					   flogicOperaSM2,
					   kXAxis,
					   31,
					   fOperaTTParam);

  fsolidOperaMS    = new G4Box("OperaMS_solid", fOperaMSSizeX/2.0, fOperaMSSizeY/2.0, fOperaMSSizeZ/2.0);

  flogicOperaMS    = new G4LogicalVolume(fsolidOperaMS, fOperaMSMaterial, "OperaMS_logic");

  fphysicOperaMS11 = new G4PVPlacement(0,
				       G4ThreeVector(fOperaMSWallsFirstPos, 0, 0),
				       "OperaMS11",
				       flogicOperaMS,
				       fphysicOperaZone,
				       false,
				       0, IsOverlap);  

  fphysicOperaMS12 = new G4PVPlacement(0,
				       G4ThreeVector(fOperaMSWallsFirstPos + 
						     fOperaMSWallsSpacing, 0, 0),
				       "OperaMS12",
				       flogicOperaMS,
				       fphysicOperaZone,
				       false,
				       0, IsOverlap);  

  fphysicOperaMS21 = new G4PVPlacement(0,
				       G4ThreeVector(fOperaMSWallsFirstPos + 
						     fOperaZoneSizeX/2.0, 0, 0),
				       "OperaMS21",
				       flogicOperaMS,
				       fphysicOperaZone,
				       false,
				       0, IsOverlap);  

  fphysicOperaMS22 = new G4PVPlacement(0,
				       G4ThreeVector(fOperaMSWallsFirstPos + 
						     fOperaMSWallsSpacing + 
						     fOperaZoneSizeX/2.0, 0, 0),
				       "OperaMS22",
				       flogicOperaMS,
				       fphysicOperaZone,
				       false,
				       0, IsOverlap);  

  G4VisAttributes* OperaZoneVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));   
  OperaZoneVisAtt->SetVisibility(false); 					  
  OperaZoneVisAtt->SetForceWireframe(false);					  
  //   OperaZoneVisAtt->SetForceSolid(true);					  
  flogicOperaZone->SetVisAttributes(OperaZoneVisAtt);				  

}
*/
void BxDetectorConstruction::AddSource() {
   
  
//AmBe - source case

  if(fSource == 4) {

  //AmBeDerlin - this is the white plastic holder of the source
  
    fsolidDerlin=new G4Polycone("Derlin_solid",
			 	0.0, twopi,
				BxReadParameters::Get()->GetNumOfZplanesDerlin(),
				BxReadParameters::Get()->GetZplanesDerlin(),
				BxReadParameters::Get()->GetRinDerlin(),
				BxReadParameters::Get()->GetRoutDerlin());

  //Stainless steel M3 screws of the Derlin holder


    G4Tubs* solidScrew_head=new G4Tubs("Screwhead_solid",
					  0.0*mm, BxReadParameters::Get()->GetScrewHeadRadius(),
					  BxReadParameters::Get()->GetScrewHeadHeight()/2.+0.5*mm, //1mm more for the overlap with the base 
			 		 0.0, twopi
			 		 );
    
    G4Tubs* solidScrew_base=new G4Tubs("Screwbase_solid",
					  0.0*mm, BxReadParameters::Get()->GetScrewBaseRadius(),
					  BxReadParameters::Get()->GetScrewBaseHeight()/2.+0.5*mm, //1mm more for the overlap with the head 
			 		 0.0, twopi
			 		 );
 
    G4RotationMatrix* rot=new G4RotationMatrix;
    fsolidScrew=new G4UnionSolid("Screw_solid", solidScrew_base, solidScrew_head, 
		    		rot, G4ThreeVector(0.,0.,BxReadParameters::Get()->GetScrewBaseHeight()/2.-1*mm));
    //1mm overlap in order not to leave contact surfaces between the two solids
   
    // Pb container of the source, nested in the derlin holder

    fsolidPb = new G4Tubs("Pb_solid",
			  0.0*mm, BxReadParameters::Get()->GetAmBePbRadius(),
			  BxReadParameters::Get()->GetAmBePbHeight()/2.,
			  0.0, twopi
			  );

   // AmBe source 
    
    fsolidAmBe = new G4Tubs("AmBe_solid",
			  0.0*cm, BxReadParameters::Get()->GetAmBeRadius(),
			  BxReadParameters::Get()->GetAmBeRadius()/2.,
			  0.0, twopi
			  );

    //Logical Volumes

    flogicDerlin = new G4LogicalVolume(fsolidDerlin,
                                       fDerlinMaterial,
                                       "Derlin_logic");

    flogicPb = new G4LogicalVolume(fsolidPb,
                                       fPbMaterial,
                                       "Pb_logic");
    
    flogicScrew = new G4LogicalVolume(fsolidScrew,
                                       fScrewMaterial,
                                       "Screw_logic");

    flogicAmBe = new G4LogicalVolume(fsolidAmBe,
                                       fAmBeMaterial,
                                       "Be_logic");


    // Physical volumes/placement

    G4ThreeVector fDerlinOrigin;
    
    fDerlinOrigin = fSourceOrigin + G4ThreeVector(0.0,0.0,-12.7*mm);

   
    fphysicDerlin = new G4PVPlacement(0,
  				    G4ThreeVector(fDerlinOrigin),
				    "Derlin",
				    flogicDerlin,
				    fphysicZoneI,
				    false,
                          	    0, IsOverlapSource);

    const G4int NumberOfScrews=6; //do not change this number unless geometry of the source changed!
    fphysicScrew=new G4VPhysicalVolume*[NumberOfScrews];

   for (G4int i=0; i<NumberOfScrews; i++){
  		ostringstream stream;
		stream << "Screw_" << i;
		G4ThreeVector traslation;
		G4double radius=18*mm, z_traslation=(6.4+0.5)*mm; //0.5mm takes into account the 1mm overlap between base and head	

		z_traslation+=BxReadParameters::Get()->GetScrewBaseHeight()/2.;

		traslation.setZ(z_traslation);
		traslation.setX(radius*cos(i*60*deg));
		traslation.setY(radius*sin(i*60*deg));
	
	   fphysicScrew[i] = new G4PVPlacement(0,
  				    traslation,
				    stream.str(),
				    flogicScrew,
				    fphysicDerlin,
				    false,
                          	    i, IsOverlapSource);
	   }

    fphysicPb = new G4PVPlacement(0,
  				    G4ThreeVector(0,0,BxReadParameters::Get()->GetAmBePbHeight()/2.+6.4*mm),
				    "Pb",
				    flogicPb,
				    fphysicDerlin,
				    false,
                          	    0, IsOverlapSource);		

    fphysicAmBe = new G4PVPlacement(0,
  				    G4ThreeVector(),
				    "Be",
				    flogicAmBe,
				    fphysicPb,
				    false,
                          	    0, IsOverlapSource);

	//optical properties of the Derlin
    
    G4double	*pointer1;
    G4double	*pointer2;
    G4int 	theNum;

    fOpDerlinSurfaceIn = new G4OpticalSurface("OpDerlinSurfaceIn");
    fOpDerlinSurfaceIn->SetType(dielectric_dielectric);
    fOpDerlinSurfaceIn->SetModel(unified);
    fOpDerlinSurfaceIn->SetFinish(groundfrontpainted);

    G4MaterialPropertiesTable* DerlinST = new G4MaterialPropertiesTable();
    theNum  = BxPropertyCollection::Get()->GetDerlinReflectivity().size();
    pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyDerlinRef();
    pointer2 = BxPropertyCollection::Get()->GetDerlinReflectivity();
    DerlinST->AddProperty("REFLECTIVITY", pointer1,  pointer2, theNum);
    theNum  = BxPropertyCollection::Get()->GetDerlinEfficiency().size();
    pointer1 = BxPropertyCollection::Get()->GetPhotonEnergyDerlinEff();
    pointer2 = BxPropertyCollection::Get()->GetDerlinEfficiency();
    DerlinST->AddProperty("EFFICIENCY", pointer1,  pointer2, theNum);

    fOpDerlinSurfaceIn->SetMaterialPropertiesTable (DerlinST);
   
    fDerlinSurfaceIn = new G4LogicalBorderSurface("DerlinSurfaceIn",
				  fphysicZoneI,
				  fphysicDerlin,
				  fOpDerlinSurfaceIn);

    fOpDerlinSurfaceOut = new G4OpticalSurface("OpDerlinSurfaceOut");
    fOpDerlinSurfaceOut->SetType(dielectric_dielectric);
    fOpDerlinSurfaceOut->SetModel(unified);
    fOpDerlinSurfaceOut->SetFinish(polished);


    fDerlinSurfaceOut = new G4LogicalBorderSurface("DerlinSurfaceOut",
				  fphysicDerlin,
				  fphysicZoneI,
				  fOpDerlinSurfaceOut);


	//optical properties for the Screws

    fOpScrewSurface = new G4OpticalSurface("OpScrewSurface");
    fOpScrewSurface->SetType(dielectric_metal);
    fOpScrewSurface->SetModel(unified);
    fOpScrewSurface->SetFinish(polished);

    //same properties as SSS (this should be not really important)
    G4MaterialPropertiesTable* ScrewST = new G4MaterialPropertiesTable();
    theNum  = BxPropertyCollection::Get()->GetSSSReflectivity().size();
    pointer1 = BxPropertyCollection::Get()->GetPhotonEnergySSSRef();
    pointer2 = BxPropertyCollection::Get()->GetSSSReflectivity();
    ScrewST->AddProperty("REFLECTIVITY", pointer1,  pointer2, theNum);
    theNum  = BxPropertyCollection::Get()->GetSSSEfficiency().size();
    pointer1 = BxPropertyCollection::Get()->GetPhotonEnergySSSEff();
    pointer2 = BxPropertyCollection::Get()->GetSSSEfficiency();
    ScrewST->AddProperty("EFFICIENCY", pointer1,  pointer2, theNum);

    fOpScrewSurface->SetMaterialPropertiesTable (ScrewST);
   
    fScrewSurface=new G4LogicalBorderSurface*[NumberOfScrews];
	    
    for (G4int i=0; i<NumberOfScrews; i++){
	ostringstream stream;
	stream << "ScrewSurface_"<< i;

    fScrewSurface[i] = new G4LogicalBorderSurface(stream.str(),
				  fphysicDerlin,
				  fphysicScrew[i],
				  fOpScrewSurface);

    }



	//optical properties for the Pb container

    fOpPbSurface = new G4OpticalSurface("OpPbSurface");
    fOpPbSurface->SetType(dielectric_metal);
    fOpPbSurface->SetModel(unified);
    fOpPbSurface->SetFinish(polished);

    //same properties as SSS (this should be not really important)
    G4MaterialPropertiesTable* PbST = new G4MaterialPropertiesTable();
    theNum  = BxPropertyCollection::Get()->GetSSSReflectivity().size();
    pointer1 = BxPropertyCollection::Get()->GetPhotonEnergySSSRef();
    pointer2 = BxPropertyCollection::Get()->GetSSSReflectivity();
    PbST->AddProperty("REFLECTIVITY", pointer1,  pointer2, theNum);
    theNum  = BxPropertyCollection::Get()->GetSSSEfficiency().size();
    pointer1 = BxPropertyCollection::Get()->GetPhotonEnergySSSEff();
    pointer2 = BxPropertyCollection::Get()->GetSSSEfficiency();
    PbST->AddProperty("EFFICIENCY", pointer1,  pointer2, theNum);

    fOpPbSurface->SetMaterialPropertiesTable (PbST);
   
    fPbSurface = new G4LogicalBorderSurface("PbSurface",
				  fphysicDerlin,
				  fphysicPb,
				  fOpPbSurface);



  }  else  {

     if(fSource == 1) {
     fsolidVialSource  = new G4Sphere("VialSource_solid",
			  fSourceRadius, fSourceRadius+fSourceVialThick,
			  0.0, twopi,
			  0.0, pi);
     BxLog(routine) << "Spherical Source "  << endlog ;
     BxLog(routine) << "Source Radius:  "  << fSourceRadius/cm << " cm"<< endlog ;
      
  } else if(fSource == 2) {
     fsolidVialSource = new G4Tubs("VialSource_solid",
			        fSourceRadius, fSourceRadius+fSourceVialThick,
				fSourceLong/2.+fSourceVialThick,
				0.0, 2*pi);
     BxLog(routine) << "Cylindrical source "  << endlog ;
     BxLog(routine) << "Source Radius:  "  << fSourceRadius/cm << " cm"<< endlog ;
     BxLog(routine) << "Source Lenght:  "  << fSourceLong/cm << " cm"<< endlog ;
  } else if(fSource == 3) {
     fsolidVialSource = new G4Box("VialSource_solid",
			        fSourceX/2.+fSourceVialThick,
			        fSourceY/2.+fSourceVialThick,
			        fSourceZ/2.+fSourceVialThick
				);
     BxLog(routine) << "Box Source "  << endlog ;
     BxLog(routine) << "Source X: "  << fSourceX/cm << " cm"<< endlog ;
     BxLog(routine) << "Source Y: "  << fSourceY/cm << " cm"<< endlog ;
     BxLog(routine) << "Source Z: "  << fSourceZ/cm << " cm"<< endlog ;
   }else {
	   BxLog(fatal) << "This kind of source has not been implemented yet. "<< endlog;
   }
   
   BxLog(routine) << "Source Thickness  =  "  << fSourceVialThick/mm << " mm"<< endlog ;
   BxLog(routine) << "Source Centered in "  << fSourceOrigin/cm    << " cm"<< endlog ;
   BxLog(routine) << "Source Material = "  <<  fSourceMaterial->GetName()  << endlog ;
   BxLog(routine) << "Source Vial Material = "  <<  fSourceVialMaterial->GetName()  << endlog ;
   
   flogicSourceVial = new G4LogicalVolume(fsolidVialSource,
                                   fSourceVialMaterial,
                                   "VialSource_logic");
				   

   fphysicSourceVial= new G4PVPlacement(0,
   		    fSourceOrigin,
   		    "VialSource",
   		    flogicSourceVial,
   		    fphysicZoneI,
   		    false,
   		    0, IsOverlapSource);

   if(fSource == 1) {
     fsolidSource  = new G4Sphere("Source_solid",
			0.0, fSourceRadius,
			0.0, twopi,
			0.0, pi);
   } else if(fSource == 2) {
     fsolidSource = new G4Tubs("Source_solid",
			        0, fSourceRadius,
				fSourceLong/2.,
				0.0, 2*pi);
   } else if(fSource == 3) {
     fsolidSource = new G4Box("Source_solid",
			        fSourceX/2.,
			        fSourceY/2.,
			        fSourceZ/2.
				);
   }
   
   flogicSource = new G4LogicalVolume(fsolidSource,
                                   fSourceMaterial,
                                   "Source_logic");
				   

   fphysicSource=	new G4PVPlacement(0, G4ThreeVector(),
  		     "Source",
                     flogicSource,
                     fphysicSourceVial,
                     false,
                     0, IsOverlapSource);

//Source optical surfaces
  
//--------------------Vial surfaces--------------------------
   fOpSourceVialOutSurface = new G4OpticalSurface("OpSourceVialOutSurface");
   fOpSourceVialOutSurface->SetType(dielectric_dielectric);
   fOpSourceVialOutSurface->SetModel(glisur);
   fOpSourceVialOutSurface->SetFinish(polished);

   fSourceVialOutSurface = new G4LogicalBorderSurface("SourceVialOutSurface",
				  fphysicSourceVial,
				  fphysicZoneI,
				  fOpSourceVialOutSurface);

   fOpSourceVialInSurface = new G4OpticalSurface("OpSourceVialInSurface");
   fOpSourceVialInSurface->SetType(dielectric_dielectric);
   fOpSourceVialInSurface->SetModel(glisur);
   fOpSourceVialInSurface->SetFinish(polished);
   
   fSourceVialInSurface = new G4LogicalBorderSurface("SourceVialInSurface",
						  fphysicSourceVial,
						  fphysicSource,
						  fOpSourceVialInSurface);


//--------------------Source surfaces--------------------------
   fOpSourceOutSurface = new G4OpticalSurface("OpSourceOutSurface");
   fOpSourceOutSurface->SetType(dielectric_dielectric);
   fOpSourceOutSurface->SetModel(glisur);
   fOpSourceOutSurface->SetFinish(polished);

   fSourceOutSurface = new G4LogicalBorderSurface("SourceOutSurface",
				  fphysicSource,
				  fphysicSourceVial,
				  fOpSourceOutSurface);


//--------------------Surface between ZoneI and the Vial--------------------------
   fOpZoneISourceSurface = new G4OpticalSurface("OpZoneISourceSurface");
   fOpZoneISourceSurface->SetType(dielectric_dielectric);
   fOpZoneISourceSurface->SetModel(glisur);
   fOpZoneISourceSurface->SetFinish(polished);

   fZoneISourceSurface = new G4LogicalBorderSurface("ZoneISourceSurface",
				  fphysicZoneI,
				  fphysicSourceVial,
				  fOpZoneISourceSurface);

  	}
  }
