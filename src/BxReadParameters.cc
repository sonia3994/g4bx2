
/** 
 * AUTHOR: Davide Franco
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
 */
// --------------------------------------------------------------------------//


//---------------------------------------------------------------------------//
#include <TF1.h>
#include <TFile.h>
#include <TDirectory.h>
#include "BxReadParameters.hh"      //Present Bx Class Headers 
//---------------------------------------------------------------------------//

#include <fstream>
#include <sstream>
#include "BxLogger.hh"
#include "BxIO.hh"
#include "G4SystemOfUnits.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include <cmath>
#include <math.h>
#include "G4PhysicalConstants.hh"

BxReadParameters* BxReadParameters::me = 0;

// singleton
BxReadParameters::BxReadParameters(){
    ReadDetectorGeometry();
    DefineDetectorVariables();
    ReadPMTPositions();
    ReadBorexProperty();
    fLightYieldScale    = 1.;
    fRZDim=5*cm;
    fRunNumber=12106;
    fDVessel=true;
    fVesselCenter= G4ThreeVector(0.,0.,0.);
    ReadPMTVetoPositions();
    /*fLightCerenkovScale = 1.;
    fRDMDecay = false ;
    fIsTimeCut = false ;
    fTimeCut = 10*ns ;
  */  fMaxParentIdForQuenching = 10000000;
}

BxReadParameters* BxReadParameters::Get() {
  if (!me)
    me = new  BxReadParameters();
  return me;
}

void BxReadParameters::DefineDetectorVariables() {
  
 BxLog(routine) << "Detector Variable Defined"   <<  endlog;   
     // default parameter values of the detector
   fWorldSizeX                  = 22000.0*mm;
   fWorldSizeY                  = 22000.0*mm;
   fWorldSizeZ                  = 30000.0*mm;

   fTankUpSizeRmax              = 9000.0*mm;
   fTankDawnSizeZ		= 7700*mm; //corrected, P.Lombardi's dwg
   fPlatformZ                   = 5*cm; //replaced with 5cm from P.Lombardi's dwg, before it was 1mm
   fTankSteelThickness		= 6*mm;
   fDisalignment		= 158*mm; //P.Lombardi's dwg

   fUpPlateSteelR			= 4*m;//P.Lombardi's dwg
   fUpPlateSteelZ			= 10*cm;
   fDownPlateSteelR			= 2*m;
   fDownPlateSteelZ			= 4*cm;
   
   fPlatformOuterRadius         =8978*mm; //P.Lombardi's dwg
   fLegCentre			=7080*mm; //P.Lombardi's dwg
   fLegExternalRadius		=162*mm; //P.Lombardi's dwg
   fLegThickness		=7*mm; //P.Lombardi's dwg
   fNlegs			=20;

   //Parameters for tank-dawn polycone (both water and tyvek)
   //Tyvek part
   fNumOfZplanesTankTyvek=3;
   fZplanesTankTyvek=new G4double[fNumOfZplanesTankTyvek];
   fRinTankTyvek=new G4double[fNumOfZplanesTankTyvek];
   fRoutTankTyvek=new G4double[fNumOfZplanesTankTyvek];

   fZplanesTankTyvek[0]=0*m;
   fZplanesTankTyvek[1]=1*m;
   fZplanesTankTyvek[2]=fTankDawnSizeZ;
  
   for (G4int i=0; i<fNumOfZplanesTankTyvek; i++) fRinTankTyvek[i]=0;
   
   fRoutTankTyvek[0]=fTankUpSizeRmax-1*m;
   fRoutTankTyvek[1]=fTankUpSizeRmax;
   fRoutTankTyvek[2]=fTankUpSizeRmax;
   
   fTyvekThickness=0.2*mm; //the Borexino Detector at the LNGS
   fTyvekSphereThickness=350*mm; //equal to OLD MC and estimated to be right from OD PMTs size

   //Water part
   fNumOfZplanesTank=3;
   fZplanesTank=new G4double[fNumOfZplanesTank];
   fRinTank=new G4double[fNumOfZplanesTank];
   fRoutTank=new G4double[fNumOfZplanesTank];

   fZplanesTank[0]=0*m;
   fZplanesTank[1]=1*m;
   fZplanesTank[2]=fTankDawnSizeZ-fTyvekThickness;
  
   for (G4int i=0; i<fNumOfZplanesTank; i++) fRinTank[i]=0;
   
   fRoutTank[0]=fTankUpSizeRmax-1*m;
   fRoutTank[1]=fTankUpSizeRmax-fTyvekThickness;
   fRoutTank[2]=fTankUpSizeRmax-fTyvekThickness;

  ReadDetectorGeometry();
     //radii
   fSSSExternalRadius          += fSSSThickness;
   fZoneIIIExternalRadius       = fSSSExternalRadius - fSSSThickness;
   fShroudExternalRadius       += fShroudThickness;
   fZoneIIExternalRadius        = fShroudExternalRadius - fShroudThickness;
   fVesselExternalRadius       += fVesselThickness;
   fZoneIExternalRadius         = fVesselExternalRadius - fVesselThickness;


    //Endcaps to support the ZoneI Vessel construction: check on An. Ianni drawings
   fNylonRingRadiusOut=300.0*mm;
   fNylonRingRadiusIn  =250.0*mm;	
   fNylonRingHeight=12.7*mm;    
  
   //Nylon bridge
   fNylonBridgeRadius=300.0*mm;
   fNylonBridgeHeight=1.0*mm;

   //IV copper struts from IV steel tube to IV big Nylon Ring
   fCopperStrutLength=690*mm;
   fCopperStrutRadius=10*mm;

   //vectors for Nylon Tube's Polycone
   fNumOfZplanesNylonTube=6;
   fZplanesNylonTube=new G4double[fNumOfZplanesNylonTube];
   fRinNylonTube=new G4double[fNumOfZplanesNylonTube];
   fRoutNylonTube=new G4double[fNumOfZplanesNylonTube];

   fZplanesNylonTube[0]=0*mm;
   fZplanesNylonTube[1]=17*mm;
   fZplanesNylonTube[2]=fZplanesNylonTube[1]+0.1*mm;
   fZplanesNylonTube[3]=fZplanesNylonTube[2]+655.5*mm;
   fZplanesNylonTube[4]=fZplanesNylonTube[3]+0.1*mm;
   fZplanesNylonTube[5]=fZplanesNylonTube[4]+25*mm;
  
   fIVTubeZposition=fZplanesNylonTube[5]+0.1*mm; //Z position of IVTube w.r.t. Nylon Tube + 0.1mm tolerance

   //fRoutNylonTube[0]=80*mm;
   //fRoutNylonTube[1]=80*mm;
   fRoutNylonTube[0]=120*mm;
   fRoutNylonTube[1]=120*mm;
   fRoutNylonTube[2]=52.5*mm;
   fRoutNylonTube[3]=52.5*mm;
   fRoutNylonTube[4]=80*mm;
   fRoutNylonTube[5]=80*mm;

   for (G4int i=0; i<fNumOfZplanesNylonTube; i++) fRinNylonTube[i]=48.5*mm;
    
   //vectors for IV Steel Tube's Polycone
   fNumOfZplanesIVTube=6;
   fZplanesIVTube=new G4double[fNumOfZplanesIVTube];
   fRinIVTube=new G4double[fNumOfZplanesIVTube];
   fRoutIVTube=new G4double[fNumOfZplanesIVTube];

   fZplanesIVTube[0]=0*mm;
   fZplanesIVTube[1]=12*mm;
   fZplanesIVTube[2]=fZplanesIVTube[1]+0.1*mm;
   fZplanesIVTube[3]=fZplanesIVTube[2]+538*mm;
   fZplanesIVTube[4]=fZplanesIVTube[3]+0.1*mm;
   fZplanesIVTube[5]=fZplanesIVTube[4]+17*mm;
   
   fRoutIVTube[0]=105*mm;
   fRoutIVTube[1]=105*mm;
   fRoutIVTube[2]=50*mm;
   fRoutIVTube[3]=50*mm;
   fRoutIVTube[4]=162.5*mm;
   fRoutIVTube[5]=162.5*mm;
   
   for (G4int i=0; i<fNumOfZplanesIVTube; i++) fRinIVTube[i]=0*mm;
   
  //vectors for OV Steel Tube's Polycone
   fNumOfZplanesOVTube=6;
   fZplanesOVTube=new G4double[fNumOfZplanesOVTube];
   fRinOVTube=new G4double[fNumOfZplanesOVTube];
   fRoutOVTube=new G4double[fNumOfZplanesOVTube];

   fZplanesOVTube[0]=0*mm;
   fZplanesOVTube[1]=85*mm;
   fZplanesOVTube[2]=fZplanesOVTube[1]+0.1*mm;
   fZplanesOVTube[3]=fZplanesOVTube[2]+48*mm;
   fZplanesOVTube[4]=fZplanesOVTube[3]+0.1*mm;
   fZplanesOVTube[5]=fZplanesOVTube[4]+1350*mm;
   
   fRoutOVTube[0]=101.6*mm;
   fRoutOVTube[1]=101.6*mm;
   fRoutOVTube[2]=127.5*mm;
   fRoutOVTube[3]=127.5*mm;
   fRoutOVTube[4]=50*mm;
   fRoutOVTube[5]=50*mm;

   for (G4int i=0; i<fNumOfZplanesOVTube; i++) fRinOVTube[i]=0*mm;

     //PMT's
   fPMTFaceSizeRmin             = fPMTFaceSizeRmax - fPMTFrameThickness;
   fPMTFaceSizeThetaDelta       = asin(fPMTFaceProjection/(2.0*fPMTFaceSizeRmax));
   fPMTZShift			= sqrt(fPMTFaceSizeRmin*fPMTFaceSizeRmin-fPMTFaceProjection*fPMTFaceProjection/4.);

   fPMTMidSizeRmax              = fPMTFaceProjection/2.0;
   fPMTMidSizeRmin              = fPMTMidSizeRmax - fPMTFrameThickness;
   fPMTSizeRmax                 = fPMTFaceSizeRmin;
   fPMTSizeRmin                 = fPMTSizeRmax - fPMTThickness;
   fPMTSizeThetaDelta           = fPMTFaceSizeThetaDelta;

   fPMTShieldSizeR1max         /= 2.0;
   fPMTShieldSizeR2max         /= 2.0;
   fPMTShieldSizeR1min          = fPMTShieldSizeR1max - fPMTShieldThickness;
   fPMTShieldSizeR2min          = fPMTShieldSizeR2max - fPMTShieldThickness;

   fPMTHousingSizeRmax        /= 2.0;
   fPMTHousingSizeRmin         = fPMTHousingSizeRmax - fPMTHousingThickness;

   //Vectors for polycones

   //Shield
   fNumOfZplanesShield=2;
   fZplanesShield=new G4double[fNumOfZplanesShield];
   fRinShield=new G4double[fNumOfZplanesShield];
   fRoutShield=new G4double[fNumOfZplanesShield];

   fZplanesShield[0]=0;
   fZplanesShield[1]=fPMTShieldSizeZ;

   fRinShield[0]=fPMTShieldSizeR1min;
   fRinShield[1]=fPMTShieldSizeR2min;

   fRoutShield[0]=fPMTShieldSizeR1max;
   fRoutShield[1]=fPMTShieldSizeR2max;

   //Housing (base)
   fNumOfZplanesHousing=4;
   fZplanesHousing=new G4double[fNumOfZplanesHousing];
   fRinHousing=new G4double[fNumOfZplanesHousing];
   fRoutHousing=new G4double[fNumOfZplanesHousing];

   fZplanesHousing[0]=0;
   fZplanesHousing[1]=fPMTHousingSizeZ+fPMTHousingShift+2*mm; //2mm added to insert the PMT in the SSS volume
   fZplanesHousing[2]=fPMTHousingShift+fPMTLength+2*mm; 
   fZplanesHousing[3]=fZplanesHousing[2]+fPMTHousingThickness; 

   fRinHousing[0]=fPMTHousingSizeRmin;
   fRinHousing[1]=fPMTHousingSizeRmin;
   fRinHousing[2]=fPMTFaceProjection/2.0-fPMTHousingThickness;
   fRinHousing[3]=0;

   fRoutHousing[0]=fPMTHousingSizeRmax;
   fRoutHousing[1]=fPMTHousingSizeRmax;
   fRoutHousing[2]=fPMTFaceProjection/2.0;
   fRoutHousing[3]=fRoutHousing[2];


   //PMT rings for PMTs without concentrator
   fNumOfZplanesPMTring=3;
   fZplanesPMTring=new G4double[fNumOfZplanesPMTring];
   fRinPMTring=new G4double[fNumOfZplanesPMTring];
   fRoutPMTring=new G4double[fNumOfZplanesPMTring];

   fZplanesPMTring[0]=1*mm; //1mm tolerance to avoid overlaps
   fZplanesPMTring[1]=fZplanesPMTring[0]+fPMTringSizeZ1;
   fZplanesPMTring[2]=fZplanesPMTring[0]+fPMTringSizeZ2; 

   fRinPMTring[0]=fPMTringSizeR1+0.5*mm; //0.5mm tolerance to avoid overlaps
   fRinPMTring[1]=fPMTringSizeR2;
   fRinPMTring[2]=fPMTringSizeR3;

   fRoutPMTring[0]=fRinPMTring[0]+fPMTringSizeThickness;
   fRoutPMTring[1]=fRinPMTring[1]+fPMTringSizeThickness;
   fRoutPMTring[2]=fRinPMTring[2]+fPMTringSizeThickness;

   //Disk PMT
    
   fPMTSizeDiskR                = 93.9*mm;//to match the concentrator's diameter+shielding
   fPMTSizeDiskH		= 251*mm; //arbitrary value
  
   //OD PMTs flanges
   fPMTmuFlangeSizeRmin=        fPMTFaceProjection/2.0+1*mm;       	
   fPMTmuFlangeSizeRmax=	fPMTFaceProjection/2.0+28*mm; 
   fPMTmuFlangeSizeZ=		15*mm;                       
  
   //measured on a real PMT by S.M. and Pagani L.
   //OD PMTs Shield (only for those on the ground)
   fNumOfZplanesShieldmu=3;
   fZplanesShieldmu=new G4double[fNumOfZplanesShieldmu];
   fRinShieldmu=new G4double[fNumOfZplanesShieldmu];
   fRoutShieldmu=new G4double[fNumOfZplanesShieldmu];

   for (G4int i=0; i<fNumOfZplanesShieldmu; i++) fRinShieldmu[i]=0;

   fZplanesShieldmu[0]=0;
   fZplanesShieldmu[1]=9.5*cm;  
   fZplanesShieldmu[2]=29.5*cm; 

   fRoutShieldmu[0]=8.5/2.*cm;
   fRoutShieldmu[1]=8.5/2.*cm;
   fRoutShieldmu[2]=fPMTFaceProjection/2.0;


	//Vectors for the polycone of the AmBe holder, references for the sizes Steven Hardy's PhD Thesis page 240

   fNumOfZplanesDerlin=5;
   fZplanesDerlin=new G4double[fNumOfZplanesDerlin];
   fRinDerlin=new G4double[fNumOfZplanesDerlin];
   fRoutDerlin=new G4double[fNumOfZplanesDerlin];

   fZplanesDerlin[0]=0*mm;
   fZplanesDerlin[1]=6.4*mm;
   fZplanesDerlin[2]=fZplanesDerlin[1]+17.6*mm;
   fZplanesDerlin[3]=fZplanesDerlin[2]+1*mm;
   fZplanesDerlin[4]=fZplanesDerlin[3]+18*mm;
  
   for (G4int i=0; i<fNumOfZplanesDerlin; i++) fRinDerlin[i]=0;
   
   fRoutDerlin[0]=13.5*mm;
   fRoutDerlin[1]=22*mm;
   fRoutDerlin[2]=22*mm;
   fRoutDerlin[3]=7.8*mm;
   fRoutDerlin[4]=6*mm;
  
   
   //Stainless steel M3 screws, ref. for sizes: http://www.stepgear.com/viti_tcei_-_passo_grosso.html
   
   //head
	
   fScrewHeadRadius=2.75*mm;
   fScrewHeadHeight=3*mm;


   //base
   
   fScrewBaseRadius=1.5*mm;
   fScrewBaseHeight=1*cm;
   
	//AmBe Pb container

   fAmBePbRadius=12.3*mm;
   fAmBePbHeight=9.8*mm;

   	//AmBe source sizes, pay attention to the height: it is chosen to reproduce the little excess of Pb on the 
	//cap of the Pb holder
	
   fAmBeRadius=9.3*mm;
   fAmBeHeight=5.6*mm;

   /*
   // Opera
   fWorldSizeXOpera             = 100*m;

   fOperaZoneSizeX              = 26.0*m;
   fOperaZoneSizeY              = 10.0*m;
   fOperaZoneSizeZ              = 10.0*m;

   fOperaSMSizeX                = 9.5*m;
   fOperaSMSizeY                = 9.5*m;
   fOperaSMSizeZ                = 9.5*m;

   fOperaTTSizeX                = 5.0*cm;
   fOperaTTSizeY                = 9.0*m;
   fOperaTTSizeZ                = 9.0*m;

   fOperaMSSizeX                = 1.0*m;
   fOperaMSSizeY                = 9.0*m;
   fOperaMSSizeZ                = 9.0*m;

   fOperaTTWallsSpacing         = 0.2*m;				     
   fOperaMSWallsSpacing         = 1.5*m;
   fOperaTTWallsFirstPos        =-4.00*m;				     
   fOperaMSWallsFirstPos        =-12.0*m;    

   fOperaZonePos                = G4ThreeVector(-33.00*m, 0.0, 0.0);
   fOperaSM1Pos                 = G4ThreeVector(  8.25*m, 0.0, 0.0);
   fOperaSM2Pos                 = G4ThreeVector( -4.75*m, 0.0, 0.0);
   
   BxLog(debugging) << "---------->" << fPMTHousing1SizeRmax << " " <<fPMTHousing1SizeRmin << " " <<  fPMTHousing2SizeRmax << " " <<
   fPMTHousing2SizeZ << endlog ;
*/
} 

void BxReadParameters::ReadDetectorGeometry() {

  G4String str;
  std::ifstream&  inDetectorGeometryFile=BxIO::Get()->GetStreamBxGeometry();

  BxIO::Get()->GetStreamLogFile() << endl ;
  BxIO::Get()->GetStreamLogFile() << "########## Detector Geometry ###########" << endl ;

  while (getline (inDetectorGeometryFile, str ) ) {
  
    BxIO::Get()->GetStreamLogFile() << str << endl ;

    if (!str.find("SSSInnerRadius")) {
      fSSSExternalRadius = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("SSSThickness")) {
      fSSSThickness = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("ShroudInnerRadius")) {
      fShroudExternalRadius = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("ShroudThickness")) {
      fShroudThickness = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("IVInnerRadius")) {
      fVesselExternalRadius = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("IVThickness")) {
      fVesselThickness = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("DetectorFlag")) {
      fDetectorFlag = (int)atof(str.substr(str.find("=")+1).c_str());
    }
    if (!str.find("PMTExternalRadius")) {
      fPMTFaceSizeRmax = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTCathodeThickness")) {
      fPMTThickness = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTFrameThickness")) {
      fPMTFrameThickness = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTExternalProjection")) {
      fPMTFaceProjection = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTShieldExternalDiameter1")) {
      fPMTShieldSizeR1max = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTShieldExternalDiameter2")) {
      fPMTShieldSizeR2max = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTShieldLength")) {
      fPMTShieldSizeZ = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTShieldZPosition")) {
      fPMTShieldZPosition = atof(str.substr(str.find("=")+1).c_str()) * mm;
      //Default shift of 20mm with respect to "official" drawings (after talking to P. Lombardi) 
      fPMTShieldZPosition += 20*mm;
    }
    if (!str.find("PMTShieldThickness")) {
      fPMTShieldThickness = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTCANShift")) {
      fPMTHousingShift = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTCANExternalDiameter")) {
      fPMTHousingSizeRmax = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTCANThickness")) {
      fPMTHousingThickness = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTCANLength")) {
      fPMTHousingSizeZ = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTLength")) {
      fPMTLength = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTLightGuideZShift")) {
      fPMTLightGuideZShift = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTringRinner1")) {
      fPMTringSizeR1 = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTringRinner2")) {
      fPMTringSizeR2 = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTringRinner3")) {
      fPMTringSizeR3 = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTringZ1")) {
      fPMTringSizeZ1 = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTringZ2")) {
      fPMTringSizeZ2 = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTringShift")) {
      fPMTringSizeShift = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }
    if (!str.find("PMTringThickness")) {
      fPMTringSizeThickness = atof(str.substr(str.find("=")+1).c_str()) * mm;
    }

  }
  BxLog(trace) << "Detector Geometry loaded" << endlog;

}

void BxReadParameters::ReadBorexProperty() {
 G4String st;
 std::ifstream&  inDetectorPropertyFile=BxIO::Get()->GetStreamBxProperty();

  while (getline (inDetectorPropertyFile, st ) ) {

    if (!st.find("PPOAttenuationLengthFactor")) 
      fPPOAttenuationLengthFactor = atof(st.substr(st.find("=")+1).c_str());
    
    if (!st.find("PCAttenuationLengthFactor")) 
      fPCAttenuationLengthFactor = atof(st.substr(st.find("=")+1).c_str());
    
    if (!st.find("DMPAttenuationLengthFactor")) 
      fDMPAttenuationLengthFactor = atof(st.substr(st.find("=")+1).c_str());
    
    if (!st.find("NylonAttenuationLengthFactor")) 
      fNylonAttenuationLengthFactor = atof(st.substr(st.find("=")+1).c_str());
    
    if (!st.find("PPOAbsReemProbXRays")){ 
      fPPOAbsReemProbXRays = atof(st.substr(st.find("=")+1).c_str());
    //G4cout<<"PPOAbsReemProbXRays in ReadParam = "<<fPPOAbsReemProbXRays<<G4endl;
   } if (!st.find("PPOAbsReemProbThXRays")){
      fPPOAbsReemProbXRaysTh = atof(st.substr(st.find("=")+1).c_str());
   // G4cout<<"PPOAbsReemProbThXRays in ReadParam = "<<fPPOAbsReemProbXRaysTh<<G4endl;
    }
    if (!st.find("PMTQEMaximum")) 
	    fPMTQEMaximum = atof(st.substr(st.find("=")+1).c_str());

    if (!st.find("PMTRelativeQEflag"))
      fRelPMTQE = atoi(st.substr(st.find("=")+1).c_str());
    
     if(!st.find("BetaPhotonYield")) 
      fLightYield = (G4int)atof(st.substr(st.find("=")+1).c_str());
    
    if (!st.find("LightYieldScale")) 
      fLightYieldScale = atof(st.substr(st.find("=")+1).c_str());

    if (!st.find("BirksBeta")) 
      fBirksBeta = atof(st.substr(st.find("=")+1).c_str());
  
    if (!st.find("SecondOrderBirksBeta")) 
      fBirksSecondOrderBeta = atof(st.substr(st.find("=")+1).c_str());

    if (!st.find("BirksAlpha")) 
      fBirksAlpha = atof(st.substr(st.find("=")+1).c_str());

    if (!st.find("SecondOrderBirksAlpha")) 
      fBirksSecondOrderAlpha = atof(st.substr(st.find("=")+1).c_str());
      
    if (!st.find("BirksProton")) 
      fBirksProton = atof(st.substr(st.find("=")+1).c_str());

    if (!st.find("SecondOrderBirksProton")) 
      fBirksSecondOrderProton = atof(st.substr(st.find("=")+1).c_str());      
      
    if (!st.find("AlphaDecayTimeConstant1")) 
      fAlphaDecayTimeConstant.push_back(atof(st.substr(st.find("=")+1).c_str()));      
      
    if (!st.find("AlphaDecayTimeConstant2")) 
      fAlphaDecayTimeConstant.push_back(atof(st.substr(st.find("=")+1).c_str()));      
      
    if (!st.find("AlphaDecayTimeConstant3")) 
      fAlphaDecayTimeConstant.push_back(atof(st.substr(st.find("=")+1).c_str()));      
      
    if (!st.find("AlphaDecayTimeConstant4")) 
      fAlphaDecayTimeConstant.push_back(atof(st.substr(st.find("=")+1).c_str()));      
      
    if (!st.find("BetaDecayTimeConstant1")) 
      fBetaDecayTimeConstant.push_back(atof(st.substr(st.find("=")+1).c_str()));      
      
    if (!st.find("BetaDecayTimeConstant2")) 
      fBetaDecayTimeConstant.push_back(atof(st.substr(st.find("=")+1).c_str()));      
      
    if (!st.find("BetaDecayTimeConstant3")) 
      fBetaDecayTimeConstant.push_back(atof(st.substr(st.find("=")+1).c_str()));      
      
    if (!st.find("BetaDecayTimeConstant4")) 
      fBetaDecayTimeConstant.push_back(atof(st.substr(st.find("=")+1).c_str()));      
      
    if (!st.find("AlphaDecayWeight1")) 
      fAlphaDecayWeight.push_back(atof(st.substr(st.find("=")+1).c_str()));      
      
    if (!st.find("AlphaDecayWeight2")) 
      fAlphaDecayWeight.push_back(atof(st.substr(st.find("=")+1).c_str()));      
      
    if (!st.find("AlphaDecayWeight3")) 
      fAlphaDecayWeight.push_back(atof(st.substr(st.find("=")+1).c_str()));      
      
    if (!st.find("AlphaDecayWeight4")) 
      fAlphaDecayWeight.push_back(atof(st.substr(st.find("=")+1).c_str()));      
      
    if (!st.find("BetaDecayWeight1")) 
      fBetaDecayWeight.push_back(atof(st.substr(st.find("=")+1).c_str()));      
      
    if (!st.find("BetaDecayWeight2")) 
      fBetaDecayWeight.push_back(atof(st.substr(st.find("=")+1).c_str()));      
      
    if (!st.find("BetaDecayWeight3")) 
      fBetaDecayWeight.push_back(atof(st.substr(st.find("=")+1).c_str()));      
     
    if (!st.find("BetaDecayWeight4")) 
      fBetaDecayWeight.push_back(atof(st.substr(st.find("=")+1).c_str()));      
    
    if (!st.find("TimePPOEmission")) 
      fTimePPOEmission=atof(st.substr(st.find("=")+1).c_str());
    
    if (!st.find("TimePCEmission")) 
      fTimePCEmission=atof(st.substr(st.find("=")+1).c_str());      
     
    if (!st.find("TimePCtoPPOTransfer")) 
      fTimePCtoPPOTransfer=atof(st.substr(st.find("=")+1).c_str());      
  }
}


void BxReadParameters::ReadPMTVetoPositions() {
 
G4double     PMTXPos,
	     PMTYPos,
	     PMTZPos;

  G4int      NumbOfElChanel, 
	     VetoFlag,
             ProfileDB;  

  std::ifstream 
  inPMTVetoFile("../data/borex/bx_pmtVeto.inp",std::ios::in);

  if (!inPMTVetoFile)
	  BxLog(fatal) << "Error while reading file ../data/borex/bx_pmtVeto.inp" << endlog;
  
  int NumberOfVetoOD_Pmt = 0;
  int NumberOfAllVetoOD_Pmt = 0;

  if (fDetectorFlag != 2) {
    while(inPMTVetoFile >> ProfileDB >> NumbOfElChanel  
                        >> PMTXPos >> PMTYPos >> PMTZPos >> VetoFlag )
    
    {

        fVetoOD_Pmt.push_back(VetoFlag);

	fPMT_VetoOD_Position.push_back(G4ThreeVector(PMTXPos*m, PMTYPos*m, PMTZPos*m));
      
    if (VetoFlag==0) NumberOfVetoOD_Pmt++;

    NumberOfAllVetoOD_Pmt++;
    
    }

  }
  BxLog(trace) << "Number of All Outer Detector Veto PMT loaded = " << NumberOfAllVetoOD_Pmt << endlog;
  BxLog(trace) << "Number of Sphere Outer Detector Veto PMT = " << NumberOfVetoOD_Pmt << endlog;

}




void BxReadParameters::ReadPMTPositions()
{
G4double	PMTXPos,
		PMTYPos,
		PMTZPos,
		dummy;

G4int           HoleNumber, 
		VetoFlag;
 
std::ifstream inLGPosDataFile("../data/borex/bx_pmtfile.inp",std::ios::in);

  if (!inLGPosDataFile)
	  BxLog(fatal) << "File for ID PMTs positioning not found" << endlog;

  fNumberOfVetoPMT = 0;
  G4int NumberOfAllPMT = 0;

 if (fDetectorFlag != 3){
    while(inLGPosDataFile >> dummy >> dummy >> HoleNumber >> dummy >> PMTXPos >> PMTYPos >> PMTZPos >>   
    dummy >> dummy >> VetoFlag >> dummy)	 {
      
      G4ThreeVector tmp(PMTXPos, PMTYPos, PMTZPos);
      tmp.setMag(fZoneIIIExternalRadius-(fPMTSizeDiskH-fSSSThickness/2)/2);
      
      
      fPMTPosition.push_back(tmp);
      NumberOfAllPMT++;
      if(VetoFlag != 1) {
	fNumberOfVetoPMT  += 1;					
	fVetoPmt.push_back(1);
      } else {
	fVetoPmt.push_back(0);   	
      }

    }
    BxLog(routine) << "Number of All PMT loaded = " << NumberOfAllPMT << endlog;
    BxLog(routine) << "Veto PMT = " << fNumberOfVetoPMT << endlog;
 }
}

void BxReadParameters::ReadLGParam() {

  G4double zvalue, inradius;
  std::vector<G4double> zValue,InValue,OutValue,MidValue;
  std::ifstream inLGParametersFile(BxIO::Get()->GetLightGuideShapeFileName().c_str(),std::ios::in);

  if (!inLGParametersFile)
	BxLog(fatal) << "Error while reading file " << BxIO::Get()->GetLightGuideShapeFileName() << endlog;

  while(inLGParametersFile >> inradius >> zvalue) {
    zValue.push_back(zvalue*mm);
    InValue.push_back(inradius*mm+0.5*mm);
    MidValue.push_back(inradius*mm+1*mm);
    OutValue.push_back(inradius*mm+1.5*mm);
  }
  
  fLGZCutsNumber=zValue.size();
  fLGZValue=new G4double[zValue.size()];
  fLGInnerRadius=new G4double[zValue.size()];
  fLGMidRadius=new G4double[zValue.size()];
  fLGOuterRadius=new G4double[zValue.size()];

  for (unsigned int i=0; i<zValue.size(); i++){
  fLGZValue[i]=zValue.at(i);
  fLGInnerRadius[i]=InValue.at(i);
  fLGMidRadius[i]=MidValue.at(i);
  fLGOuterRadius[i]=OutValue.at(i);
  }
  BxLog(trace) << "LG parameters loaded" << endlog;
}

G4bool BxReadParameters::ReadVesselGeometry() {

//	if (BxIO::Get()->GetStreamDst()->IsZombie()) 
//		return false;

//	VesselDir=(TDirectory*)BxIO::Get()->GetStreamDst()->Get("vessel_shape");
//	VesselShape=(TF1*)VesselDir->Get("Vessel_2D");

	//I assume that the radius given by the file is the external one. S.M.
	//what is the precision of this vessel shape measurement?

	//I compute the number of Z steps for the polycone
/*
        fZCutsNumber=(int)(pi*4250/fRZDim*mm);
	fVessel_Rout=new G4double[fZCutsNumber+1];
	fVessel_Rin=new G4double[fZCutsNumber+1];
	fVessel_Z=new G4double[fZCutsNumber+1];
	fZoneI_Z=new G4double[fZCutsNumber+1];

	fVesselCenter.setX(0.);
	fVesselCenter.setY(0.);
	fVesselCenter.setZ(0.);

	//Loop for Vessel's shape
	for (int i=1; i<fZCutsNumber; i++){
		fVessel_Z[i]=fVesselCenter.z()*cm+VesselShape->Eval(i*pi/fZCutsNumber)*cos(i*pi/fZCutsNumber)*m;
		fVessel_Rout[i]=VesselShape->Eval(i*pi/fZCutsNumber)*sin(i*pi/fZCutsNumber)*m;
	}

	//Loop for ZoneI's shape
	for (int i=1; i<fZCutsNumber; i++){
		fZoneI_Z[i]=fVesselCenter.z()*cm+(VesselShape->Eval(i*pi/fZCutsNumber)*m-fVesselThickness)*cos(i*pi/fZCutsNumber);
		fVessel_Rin[i]=(VesselShape->Eval(i*pi/fZCutsNumber)*m-fVesselThickness)*sin(i*pi/fZCutsNumber);

		if (fVessel_Rin[i]<0)
			fVessel_Rin[i]=0;
	}

	fVessel_Z[0]=fVesselCenter.z()*cm+VesselShape->Eval(0)*m;
	fVessel_Z[fZCutsNumber]=fVesselCenter.z()*cm-VesselShape->Eval(pi-0.000001)*m;

	fZoneI_Z[0]=fVesselCenter.z()*cm+VesselShape->Eval(0)*m-fVesselThickness;
	fZoneI_Z[fZCutsNumber]=fVesselCenter.z()*cm-VesselShape->Eval(pi-0.000001)*m+fVesselThickness;

	if (fVessel_Z[fZCutsNumber]==0 && fVessel_Z[0]==0)
		BxLog(fatal) << "Boundaries of the Vessel Shape function are zero. Impossible to create vessel's shape. " << endlog;

	if (fVessel_Z[fZCutsNumber]==0){
		fVessel_Z[fZCutsNumber]=-fVessel_Z[0];
		BxLog(warning) << "The boundary at pi of the vessel's TF1 was zero." << endlog;
	}

	if (fVessel_Z[0]==0){
		fVessel_Z[0]=-fVessel_Z[fZCutsNumber];
		BxLog(warning) << "The boundary at zero of the vessel's TF1 was zero." << endlog;
	}

	fVessel_Rout[0]=0;
	fVessel_Rin[0]=0;
	fVessel_Rout[fZCutsNumber]=0;
	fVessel_Rin[fZCutsNumber]=0;

	double tmp1=0, tmp2=0;
	for (int i=0; i<=fZCutsNumber;i++){
		if (abs(fVessel_Z[i])>tmp1)
			tmp1=abs(fVessel_Z[i]);
		if (abs(fZoneI_Z[i])>tmp2)
			tmp2=abs(fZoneI_Z[i]);
	}

	fZCutsNumber++;
*/	
	BxLog(routine) << "Vessel Shape Loaded" << endlog;

	return true;
}

std::vector<G4double> BxReadParameters::ConvertTo4DoubleVector(const char* st) {
	G4double v1;
	std::istringstream is(st);
	std::vector<G4double> tmp1, tmp2;
	while (is){
		is >> v1 ;
		tmp1.push_back(v1);
	}

	for (unsigned int i=0; i<tmp1.size()-1; i++)
		tmp2.push_back(tmp1[i]);

	return tmp2;
}
