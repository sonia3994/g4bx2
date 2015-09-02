//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#include "BxMaterial.hh"
#include "BxLogger.hh"
#include "BxOutputVertex.hh"
#include "BxPropertyCollection.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "HistoManager.hh"
using namespace std;

BxMaterial* BxMaterial:: me=0;

BxMaterial* BxMaterial::Get() {
if (!me)
 me = new BxMaterial();
 return me;
}


BxMaterial::BxMaterial() {
  DefineMaterials();
  DefineProperties();
  SetMaterialIndexes();
}

BxMaterial::~BxMaterial(){}


void BxMaterial::DefineMaterials() {

  G4String 	name,
		symbol;
  G4double	a,      	//a=mass of a mole;
		density;  
  G4int		z;    		//z=mean number of protons;	     
  G4int 	ncomponents,	     
		natoms; 					     
  G4double 	fractionmass;


  G4Element* H  = new G4Element(name="Hydrogen", symbol="H" , z= 1, a=1.00794*g/mole);
 
  G4Element* Be  = new G4Element(name="Beryllium",  symbol="Be",  z=4,  a=9.012182*g/mole);
  
  G4Element* C  = new G4Element(name="Carbon",  symbol="C",  z=6,  a=12.0107*g/mole);
   
  G4Element* N  = new G4Element(name="Nitrogen", symbol="N",  z=7,  a=14.00674*g/mole);

  G4Element* O  = new G4Element(name="Oxygen",  symbol="O",  z=8,  a=15.9994*g/mole);
 
  G4Element* Fe = new G4Element(name="Iron",   symbol="Fe", z=26, a=55.845*g/mole);

  G4Element* Si = new G4Element(name="Silicium",symbol="Si", z=14, a=28.0855*g/mole);

  G4Element* K  = new G4Element(name="Potassium",symbol="K", z=19, a=39.0983*g/mole);

  G4Element* Cs = new G4Element(name="Cesium",symbol="Cs", z=55, a=132.90545*g/mole);

  //G4Element* Mg = new G4Element(name="Magnesium",symbol="Mg",z=12, a=24.3050*g/mole);
 
  //G4Element* Al = new G4Element(name="Aluminum", symbol="Al", z=13, a=26.981538*g/mole);
  
  //G4Element* Ca = new G4Element(name="Calcium"  ,symbol="Ca", z=20, a=40.078*g/mole);

  //G4Element* Na = new G4Element(name="Sodium"   ,symbol="Na", z=11, a=22.989770*g/mole); 
 // G4Element* Cr  = new G4Element("Chromium",symbol="Cr",z = 24,a=51.9961*g/mole);
  
  //G4Element* Mn = new G4Element("Manganese",symbol="Mn",z = 25,a=54.938049*g/mole);  
  
  G4Element* Ni = new G4Element("Nickel",symbol="Ni",z = 28,a=58.6934*g/mole);

  G4Element* Cu = new G4Element("Copper",symbol="Cu",z = 29,a=63.546*g/mole);
  
  //G4Element* P = new G4Element(name="Phosphorus", symbol="P", z=15, a=30.973761*g/mole );
  
  //G4Element* S = new G4Element(name="Sulfur", symbol="S", z=16,a=32.066*g/mole );
  
  //G4Element* Ti = new G4Element(name="Titanium", symbol="Ti", z=22, a=47.867*g/mole);

  G4Element* Mo = new G4Element(name="Molybdenum", symbol="Mo", z=42, a=95.94*g/mole);

  G4Element* Sb =new G4Element(name="Antimonium", symbol="Sb", z=51, a=121.760*g/mole);
  
  G4Element* Am =new G4Element(name="Americium", symbol="Am", z=92, a=235.3*g/mole); //SM: Warning!! real Am has z=95 and a=243g/mole. 
  										    // Geant4 complains with z>92, so a is scaled to preserve the same Z/A

// masses taken from http://www.chemicalelements.com/show/mass.html


//
// define a material from elements.   case: chemical molecule
//

  density = 8.96*g/cm3; //from http://www.chemicalelements.com/elements/cu.html
  fCopper = new G4Material(name="Copper", density, ncomponents=1);
  fCopper->AddElement(Cu, natoms=1);
  
  density = 1.06*g/cm3; //from http://www.lxrcxc.chemchina.com/rcfsten/cpyfw/ppysb/webinfo/2012/06/1338782867275246.htm
  fPPO = new G4Material(name="PPO", density, ncomponents=4);
  fPPO->AddElement(C, natoms=15);
  fPPO->AddElement(H, natoms=11);
  fPPO->AddElement(N, natoms=1);
  fPPO->AddElement(O, natoms=1);

  density = 1.19*g/cm3; //from http://en.wikipedia.org/wiki/Dimethyl_phthalate
  fDMP = new G4Material(name="DMP", density, ncomponents=3);
  fDMP->AddElement(C, natoms=10);
  fDMP->AddElement(H, natoms=10);
  fDMP->AddElement(O, natoms=4);

  density = 0.882*g/cm3; //from http://www.chemicalbook.com/ChemicalProductProperty_EN_CB0854702.htm
  fPC = new G4Material(name="PC", density, ncomponents=2);
  fPC->AddElement(C, natoms=9);
  fPC->AddElement(H, natoms=12);

  density = 1.00*g/cm3;
  fWater = new G4Material(name="Water", density, ncomponents=2);
  fWater->AddElement(H,natoms=2);
  fWater->AddElement(O,natoms=1);

  fAluminum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

  density = 2.2*g/cm3; //should be corrected
  fGlass = new G4Material(name="Glass", density, ncomponents=2);
  fGlass->AddElement(Si, natoms=1);
  fGlass->AddElement(O, natoms=2);

    // define a material from elements.   case: mixture by fractional mass


fSteel= G4NistManager::Instance()->FindOrBuildMaterial("G4_STAINLESS-STEEL");
//material for PMT window
//fGlass=G4NistManager::Instance()->FindOrBuildMaterial("G4_GLASS_PLATE");


//Tyvek
  
  density=0.38*g/cm3; //http://www2.dupont.com/Tyvek
  fTyvek=new G4Material(name="Tyvek", density, ncomponents=2); //http://www.academia.edu/5490993/Tyvek_By_Dupont
  fTyvek->AddElement(H, natoms=6);
  fTyvek->AddElement(C, natoms=2);


//check on nylon from J. Benziger et al., Nucl. Instr. & Meth. A 582509 (2007) 
  density = 1.13*g/cm3; //http://www.matweb.com,  material Nylon 6 (DAM)
  fNylon = new G4Material(name="Nylon", density, ncomponents=4);//corrected with formula C6 H11 N O
  fNylon->AddElement(H, fractionmass=0.10);
  fNylon->AddElement(C, fractionmass=0.64);
  fNylon->AddElement(O, fractionmass=0.14);
  fNylon->AddElement(N, fractionmass=0.12);

  /*
  fNylonL = new G4Material(name="NylonL", density, ncomponents=3);
  fNylonL->AddElement(H, fractionmass=0.08);
  fNylonL->AddElement(C, fractionmass=0.60);
  fNylonL->AddElement(O, fractionmass=0.32);
*/
  density = 1.290*mg/cm3;
  fAir = new G4Material(name="Air", density, ncomponents=2);
  fAir->AddElement(N, fractionmass=0.7);
  fAir->AddElement(O, fractionmass=0.3);

 density = 1.29e-20*g/cm3;
 fVacuum = new G4Material(name="Vacuum", density, ncomponents=2);
 fVacuum->AddElement(N, fractionmass=0.7);
 fVacuum->AddElement(O, fractionmass=0.3);

  density = 4.000*g/cm3;  //should be corrected, but not important
  fBialkali = new G4Material(name="Bialkali", density, ncomponents=3);
  fBialkali->AddElement(K, natoms=1);
  fBialkali->AddElement(Cs, natoms=1);
  fBialkali->AddElement(Sb, natoms=1);
  
 density = 8.7*g/cm3; //http://mumetal.co.uk/?p=100
 fMuMetal = new G4Material(name="MuMetal", density, ncomponents=3);
 fMuMetal->AddElement(Ni, fractionmass=0.8);
 fMuMetal->AddElement(Mo, fractionmass=0.045);
 fMuMetal->AddElement(Fe, fractionmass=0.155);

 /*
  density = 4.000*g/cm3;  //should be corrected, but not important
  fPaint = new G4Material(name="Paint", density, ncomponents=2);
  fPaint->AddElement(K, fractionmass=0.5);
  fPaint->AddElement(Cs, fractionmass=0.5);
*/
    // define a material from others materials (mixture of mixtures)
  density = 0.8774*g/cm3; //Aldo and Corrado measurement for T=15 degrees Celcium
  fScintillator = new G4Material(name="Scintillator", density, ncomponents=2);
  fScintillator->AddMaterial(fPC, fractionmass=99.83*perCent);
  fScintillator->AddMaterial(fPPO, fractionmass=0.17*perCent);

    // define a material from others materials (mixture of mixtures)
  density = 0.8776*g/cm3; //Density of the buffer reported by Frank Calaprice is by 0.2 kg/m^3 higher than the scintillator 
                          //2.5 g/l DMP - October 2009 
  fDMPbuffer = new G4Material(name="DMPbuffer", density, ncomponents=2);
  fDMPbuffer->AddMaterial(fPC, fractionmass=99.77*perCent);
  fDMPbuffer->AddMaterial(fDMP, fractionmass=0.23*perCent);

  density= 2.200*g/cm3;// http://heraeus-quarzglas.com/media/webmedia_local/downloads/broschren_mo/dataandproperties_optics_fusedsilica.pdf
  fQuartz = new G4Material (name="Quartz", density, ncomponents=2);
  fQuartz->AddElement(Si, 1);
  fQuartz->AddElement(O , 2);

  //AmBe - density to be checked (currently it's beryllium density)
  density = 1.848*g/cm3;
  fAmBe = new G4Material(name="AmBe", density, ncomponents=3);
  fAmBe->AddElement(O, natoms=2);
  fAmBe->AddElement(Am, natoms=1); 
  fAmBe->AddElement(Be, natoms=500); 
  
  density = 1.41*g/cm3; //http://www.indat.it/PDF/TekDelrin.pdf
  fDerlin = new G4Material(name="Derlin", density, ncomponents=3);
  fDerlin->AddElement(H, natoms=2);
  fDerlin->AddElement(C, natoms=1);
  fDerlin->AddElement(O, natoms=1);  

  fPb=G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
  //fPb =  new G4Material("Pb", z=82, a=207.2*g/mole, density=11.34*g/cm3);
 
 /* 
   //Rock: hep-ex/0312050v2 Wulandari et al.
  density = 2.71*g/cm3; //
  fRock = new G4Material(name="Rock", density, ncomponents=7);
  fRock->AddElement(C,  fractionmass=0.1188);
  fRock->AddElement(O,  fractionmass=0.4892); // <-- here I added 1.01% to obtain a fractionmass = 1
  fRock->AddElement(Mg, fractionmass=0.0558);
  fRock->AddElement(Al, fractionmass=0.0103);
  fRock->AddElement(Si, fractionmass=0.0127);
  fRock->AddElement(K,  fractionmass=0.0103);
  fRock->AddElement(Ca, fractionmass=0.3029);

  //Concrete: hep-ex/0312050v2 Wulandari et al.
  density = 2.4*g/cm3; //
  fConcrete = new G4Material(name="Concrete", density, ncomponents=13);
  fConcrete->AddElement(H,  fractionmass=0.0089);
  fConcrete->AddElement(C,  fractionmass=0.0799);
  fConcrete->AddElement(O,  fractionmass=0.4971);// <-- here I add 1.28% to obtain a fractionmass = 1
  fConcrete->AddElement(Na, fractionmass=0.0006);
  fConcrete->AddElement(Mg, fractionmass=0.0085);
  fConcrete->AddElement(Al, fractionmass=0.0009);
  fConcrete->AddElement(Si, fractionmass=0.0386);
  fConcrete->AddElement(P,  fractionmass=0.0004);
  fConcrete->AddElement(S,  fractionmass=0.0016);
  fConcrete->AddElement(K,  fractionmass=0.0054);
  fConcrete->AddElement(Ca, fractionmass=0.3534);// <-- here I add 1.28% to obtain a fractionmass = 1
  fConcrete->AddElement(Ti, fractionmass=0.0004);
  fConcrete->AddElement(Fe, fractionmass=0.0043);
*/
 

}


void BxMaterial::DefineProperties() {
    // Generate & Add Material Properties Table
  G4int	        theNum;
  G4double	*pointer1;
  G4double	*pointer2;

    //  Properties for Cherenkov Light in Water
    //    This part describes the water in inner detector.
    //The data are taken from SuperK detector water
 
  const G4int NUMENTRIES = 25;

  G4double PPCKOV[NUMENTRIES] =  
            {   
  1.24108*eV, 1.82330*eV, 1.87855*eV, 1.93725*eV, 1.99975*eV, 
  2.06640*eV, 2.13766*eV, 2.21400*eV, 2.29600*eV, 2.38431*eV,  
  2.47968*eV, 2.58301*eV, 2.69531*eV, 2.81782*eV, 2.95201*eV,  
  3.09961*eV, 3.26274*eV, 3.44401*eV, 3.64660*eV, 3.87451*eV,  
  4.13281*eV, 4.42801*eV, 4.76862*eV, 5.16601*eV, 10.9721*eV  	    
		};
//TO BE CONTROLLED 
G4double RINDEX1[NUMENTRIES] =
            { 1.3400, 1.3402, 1.3402, 1.3413, 1.3424, 
	      1.3435, 1.3440, 1.3450, 1.3460, 1.3470,
	      1.3480, 1.3492, 1.3505, 1.3518, 1.3530,
	      1.3540, 1.3550, 1.3560, 1.3572, 1.3585,
	      1.3595, 1.3608, 1.3610, 1.3610, 1.3610
	      };

  G4double ABSORPTION1[NUMENTRIES] =
          {
	   344.8*cm,     666.2*cm, 1250.0*cm,  2000.4*cm,  3334.6*cm,
	   4545.9*cm,   5555.2*cm, 8334.1*cm, 10000.8*cm, 12500.0*cm, 
	   12500.6*cm, 11111.4*cm, 8345.5*cm,  5561.9*cm,  3571.2*cm,
	   2526.2*cm,   2127.6*cm, 1754.2*cm,  1000.2*cm,   500.9*cm,
	   400.5*cm,     357.7*cm,  322.7*cm,   285.3*cm,   200.0*cm
	   };

  // This is for Air - refraction index=1.0 
  G4double RINDEX2[NUMENTRIES] =
            { 
	      1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00
	    };

    //This part describes the water in Outer detector
    //The data are taken from SuperK detector water 



  const G4int NUMENTRIESair = 2;
//TO BE CHANGED AFTER CHECK 
  G4double ABSORPTIONair[NUMENTRIESair] = {100000*cm,100000*cm};
  G4double PPCKOVair[NUMENTRIESair]={1.24108*eV,10.9721*eV};


  BxLog(trace) << "Starting to load material properties" << endlog; 

  //Properties for Vacuum 
  //MaterialIndex == 1
  G4MaterialPropertiesTable *myMPT1 = new G4MaterialPropertiesTable();
  fVacuum->SetMaterialPropertiesTable(myMPT1);
/*
  fRock->SetMaterialPropertiesTable(myMPT1);
  fConcrete->SetMaterialPropertiesTable(myMPT1);

  //Properties for falseWater without Cherenkov
  //MaterialIndex == 2
  //G4MaterialPropertiesTable *myMPT2 = new G4MaterialPropertiesTable();
  //fFalseWater->SetMaterialPropertiesTable(myMPT2);
*/
  //Properties for Water Cherenkov Inner Detector
  //MaterialIndex == 3  
  G4MaterialPropertiesTable *myMPT3 = new G4MaterialPropertiesTable();
  myMPT3->AddProperty("RINDEX", PPCKOV, RINDEX2, NUMENTRIES);
  myMPT3->AddProperty("ABSLENGTH",PPCKOV, ABSORPTION1, NUMENTRIES);
  fWater->SetMaterialPropertiesTable(myMPT3);
  
  //Properties for Nylon
  //MaterialIndex == 5
  G4MaterialPropertiesTable *myMPT5 = new G4MaterialPropertiesTable();

  theNum = BxPropertyCollection::Get()->GetNylonAbsorptionLength().size();
  pointer1 = BxPropertyCollection::Get()->GetNylonAbsorptionLength();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyNylonAb();
  myMPT5->AddProperty("ABSLENGTH",pointer2, pointer1, theNum);
  //HistoManager::Get()->CreateGraph(10,pointer2,pointer1,theNum);

  theNum = BxPropertyCollection::Get()->GetNylonRefractionIndex().size();
  pointer1 = BxPropertyCollection::Get()->GetNylonRefractionIndex();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyNylonRef();
  myMPT5->AddProperty("RINDEX",pointer2, pointer1, theNum);
  //HistoManager::Get()->CreateGraph(11,pointer2,pointer1,theNum);

  fNylon->SetMaterialPropertiesTable(myMPT5);

  //Properties for Glass
  //MaterialIndex == 6
  G4MaterialPropertiesTable *myMPT6 = new G4MaterialPropertiesTable();

  theNum = BxPropertyCollection::Get()->GetGlassAbsorptionLength().size();
  pointer1 = BxPropertyCollection::Get()->GetGlassAbsorptionLength();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyGlassAb();
  //HistoManager::Get()->CreateGraph(10,pointer2,pointer1,theNum);
  myMPT6->AddProperty("ABSLENGTH",pointer2, pointer1, theNum);

  theNum = BxPropertyCollection::Get()->GetGlassRefractionIndex().size();
  pointer1 = BxPropertyCollection::Get()->GetGlassRefractionIndex();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyGlassRef();
  myMPT6->AddProperty("RINDEX",pointer2, pointer1, theNum);
  //HistoManager::Get()->CreateGraph(11,pointer2,pointer1,theNum);
  fGlass->SetMaterialPropertiesTable(myMPT6);

  //Properties for Bialkali
  //MaterialIndex == 7
  G4MaterialPropertiesTable *myMPT7 = new G4MaterialPropertiesTable();
  theNum = BxPropertyCollection::Get()->GetGlassRefractionIndex().size();
  pointer1 = BxPropertyCollection::Get()->GetGlassRefractionIndex();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyGlassRef();
  myMPT7->AddProperty("RINDEX",pointer2, pointer1, theNum);
  fBialkali->SetMaterialPropertiesTable(myMPT7);

  //Properties for Steel
  //MaterialIndex == 8  
  G4MaterialPropertiesTable *myMPT8 = new G4MaterialPropertiesTable();
  fSteel->SetMaterialPropertiesTable(myMPT8);

  //Properties for MuMetal
  //MaterialIndex == 9  
  G4MaterialPropertiesTable *myMPT8a = new G4MaterialPropertiesTable();
  fMuMetal->SetMaterialPropertiesTable(myMPT8a);
  
  //Properties for Tyvek
  //MaterialIndex == 10  
  G4MaterialPropertiesTable *myMPT8b = new G4MaterialPropertiesTable();
  fTyvek->SetMaterialPropertiesTable(myMPT8b);
  
  //Properties for Air
  //MaterialIndex == 11
  G4MaterialPropertiesTable *myMPT9 = new G4MaterialPropertiesTable();
  myMPT9->AddProperty("RINDEX", PPCKOV, RINDEX2, NUMENTRIES);
  myMPT9->AddProperty("ABSLENGTH", PPCKOVair, ABSORPTIONair, NUMENTRIESair);
  fAir->SetMaterialPropertiesTable(myMPT9);

/*
  //Properties for Paint
  //MaterialIndex == 10
  G4MaterialPropertiesTable *myMPT10 = new G4MaterialPropertiesTable();
  fPaint->SetMaterialPropertiesTable(myMPT10);

  //Properties for NylonL
  //MaterialIndex == 11
  G4MaterialPropertiesTable *myMPT11 = new G4MaterialPropertiesTable();
  fNylonL->SetMaterialPropertiesTable(myMPT11);
*/


  //Properties for PC
  //MaterialIndex == 12

  G4MaterialPropertiesTable *myMPT12 = new G4MaterialPropertiesTable();
  
  theNum = BxPropertyCollection::Get()->GetPCAbsorptionLength().size();
  pointer1 = BxPropertyCollection::Get()->GetPCAbsorptionLength();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyPCAb();
  //HistoManager::Get()->CreateGraph(0,pointer2,pointer1,theNum);
  myMPT12->AddProperty("ABSLENGTH",pointer2, pointer1, theNum);

  theNum = BxPropertyCollection::Get()->GetPCRefractionIndex().size();
  pointer1 = BxPropertyCollection::Get()->GetPCRefractionIndex();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyPCRef();
  //HistoManager::Get()->CreateGraph(1,pointer2,pointer1,theNum);
  myMPT12->AddProperty("RINDEX",pointer2, pointer1, theNum);

  theNum = BxPropertyCollection::Get()->GetPCEmissionSpectra().size();
  pointer1 = BxPropertyCollection::Get()->GetPCEmissionSpectra();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyPCEm();
  //HistoManager::Get()->CreateGraph(2,pointer2,pointer1,theNum);
  myMPT12->AddProperty("EMISSION",pointer2, pointer1, theNum);

  theNum = BxPropertyCollection::Get()->GetPCReemissionProb().size();
  pointer1 = BxPropertyCollection::Get()->GetPCReemissionProb();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyPCRP();
  //HistoManager::Get()->CreateGraph(3,pointer2,pointer1,theNum);
  myMPT12->AddProperty("REEMISPROB",pointer2, pointer1, theNum);

  theNum = BxPropertyCollection::Get()->GetPCScatteringLength().size();
  pointer1 = BxPropertyCollection::Get()->GetPCScatteringLength();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyPCAb();
  //HistoManager::Get()->CreateGraph(4,pointer2,pointer1,theNum);
  myMPT12->AddProperty("RAYLEIGH",pointer2, pointer1, theNum);

  fPC->SetMaterialPropertiesTable(myMPT12);

  //Properties for PPO
  //MaterialIndex == 13

  G4MaterialPropertiesTable *myMPT13 = new G4MaterialPropertiesTable();

  theNum = BxPropertyCollection::Get()->GetPPOAbsorptionLength().size();
  pointer1 = BxPropertyCollection::Get()->GetPPOAbsorptionLength();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyPPOAb();
  //HistoManager::Get()->CreateGraph(5,pointer2,pointer1,theNum);
  myMPT13->AddProperty("ABSLENGTH",pointer2, pointer1, theNum);
for(int i=0;i<theNum;i++)

  theNum = BxPropertyCollection::Get()->GetPPOEmissionSpectra().size();
  pointer1 = BxPropertyCollection::Get()->GetPPOEmissionSpectra();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyPPOEm();
  //HistoManager::Get()->CreateGraph(6,pointer2,pointer1,theNum);
  myMPT13->AddProperty("EMISSION",pointer2, pointer1, theNum);

  theNum = BxPropertyCollection::Get()->GetPPOReemissionProb().size();
  pointer1 = BxPropertyCollection::Get()->GetPPOReemissionProb();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyPPORP();
  //HistoManager::Get()->CreateGraph(7,pointer2,pointer1,theNum);
  myMPT13->AddProperty("REEMISPROB",pointer2, pointer1, theNum);


  fPPO->SetMaterialPropertiesTable(myMPT13);
  
//Properties for Scintillator
  //MaterialIndex == 14

  G4MaterialPropertiesTable *myMPT14 = new G4MaterialPropertiesTable();

  theNum = BxPropertyCollection::Get()->GetScintRefractionIndex().size();
  pointer1 = BxPropertyCollection::Get()->GetScintRefractionIndex();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyScintRef();
  //HistoManager::Get()->CreateGraph(8,pointer2,pointer1,theNum);
  myMPT14->AddProperty("RINDEX",pointer2, pointer1, theNum);

  theNum = BxPropertyCollection::Get()->GetPCScatteringLength().size();
  pointer1 = BxPropertyCollection::Get()->GetPCScatteringLength();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyPCAb();
  //HistoManager::Get()->CreateGraph(9,pointer2,pointer1,theNum);
  myMPT14->AddProperty("RAYLEIGH",pointer2, pointer1, theNum);

  fScintillator->SetMaterialPropertiesTable(myMPT14);

  //Properties for DMPbuffer
  //MaterialIndex == 15

  G4MaterialPropertiesTable *myMPT15 = new G4MaterialPropertiesTable();

  theNum = BxPropertyCollection::Get()->GetScintRefractionIndex().size();
  pointer1 = BxPropertyCollection::Get()->GetScintRefractionIndex();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyScintRef();
  myMPT15->AddProperty("RINDEX",pointer2, pointer1, theNum);
 
  theNum = BxPropertyCollection::Get()->GetPCScatteringLength().size();
  pointer1 = BxPropertyCollection::Get()->GetPCScatteringLength();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyPCAb();
  myMPT15->AddProperty("RAYLEIGH",pointer2, pointer1, theNum);

  fDMPbuffer->SetMaterialPropertiesTable(myMPT15);
  
  //Properties for DMP
  //MaterialIndex == 16  
  
  G4MaterialPropertiesTable *myMPT16 = new G4MaterialPropertiesTable();

  theNum = BxPropertyCollection::Get()->GetDMPAbsorptionLength().size();
  pointer1 = BxPropertyCollection::Get()->GetDMPAbsorptionLength();
  pointer2 = BxPropertyCollection::Get()->GetPhotonEnergyDMPAb();
  myMPT16->AddProperty("ABSLENGTH",pointer2, pointer1, theNum);

  fDMP->SetMaterialPropertiesTable(myMPT16);

//Properties of Quartz

  const G4int NUMENTRIESQ = 3;
 //these numbers have to be checked! 
  G4double quartz_PP[NUMENTRIESQ]   = { 1.24108*eV, 2.55*eV,  10.9721*eV  }; // lambda range 4 ri
  G4double quartz_RIND[NUMENTRIESQ] = { 1.4565 , 1.4632,   1.52   };     // ref index
  G4double quartz_ABSL[NUMENTRIESQ] = { 10.0*cm, 10.0*cm,  10.0*cm };// atten length
  
  G4MaterialPropertiesTable *quartz_mt = new G4MaterialPropertiesTable();
  quartz_mt->AddProperty("RINDEX", quartz_PP, quartz_RIND, NUMENTRIESQ);
  quartz_mt->AddProperty("ABSLENGTH", quartz_PP, quartz_ABSL, NUMENTRIESQ);
  fQuartz->SetMaterialPropertiesTable(quartz_mt);

  G4MaterialPropertiesTable *AmBe_Container = new G4MaterialPropertiesTable(); 
  fAmBe->SetMaterialPropertiesTable(AmBe_Container);

  G4MaterialPropertiesTable *Pb_Container = new G4MaterialPropertiesTable(); 
  fPb->SetMaterialPropertiesTable(Pb_Container);

  G4MaterialPropertiesTable *Derlin_Container = new G4MaterialPropertiesTable(); 
  fDerlin->SetMaterialPropertiesTable(Derlin_Container);

  G4MaterialPropertiesTable *aluminum = new G4MaterialPropertiesTable(); 
  fAluminum->SetMaterialPropertiesTable(aluminum);

  G4MaterialPropertiesTable *copper = new G4MaterialPropertiesTable(); 
  fCopper->SetMaterialPropertiesTable(copper);


 BxPropertyCollection::Get()->BuildPhysicsTable();  
  BxLog(trace) << "Material properties loaded" << endlog;

}


void BxMaterial::SetMaterialIndexes() { 
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  for(int i=0;i< (int) G4Material::GetNumberOfMaterials(); ++i ) {
    G4Material* aMaterial = (*theMaterialTable)[i];
    if(aMaterial->GetName() == "Scintillator") BxOutputVertex::Get()->SetScintillatorIndex((int)aMaterial->GetIndex()) ;
    if(aMaterial->GetName() == "DMPbuffer")    BxOutputVertex::Get()->SetDMPBufferIndex((int)aMaterial->GetIndex()) ;
    
  }

}
