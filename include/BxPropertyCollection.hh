//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef	BXPROPERTYCOLLECTION_H
#define BXPROPERTYCOLLECTION_H
#include "globals.hh"
#include "templates.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4PhysicalConstants.hh"
#include <fstream>
#include "BxSizedArray.hh"
class	G4PhysicsOrderedFreeVector;

/**
 * This class reads the parameters about optical properties and scintillation from BxReadParameters and stores them
*/
class BxPropertyCollection {
  friend class BxDetectorConstruction;

  public:

	virtual ~BxPropertyCollection();
	bool	ReadProperties();

	 sized_array<G4double>&    GetPCAbsorptionLength() {return 
	   					   fPCAbsorptionLength; }
	 sized_array<G4double>&    GetPCScatteringLength() {return 
	   					   fPCScatteringLength; }
	 sized_array<G4double>&    GetPPOAbsorptionLength() {return
	   					   fPPOAbsorptionLength; }
	 sized_array<G4double>&    GetWaterAbsorptionLength() {return
	   					   fWaterAbsorptionLength;}
	 sized_array<G4double>&    GetNylonAbsorptionLength() {return
	   					   fNylonAbsorptionLength;}
	 sized_array<G4double>&    GetGlassAbsorptionLength() {return
	   					   fGlassAbsorptionLength;}
	 sized_array<G4double>&    GetDMPAbsorptionLength() {return
	   					   fDMPAbsorptionLength; }



	 sized_array<G4double>&    GetPCEmissionSpectra() {return
	   					   fPCEmissionSpectra; }
	 sized_array<G4double>&    GetPPOEmissionSpectra() {return
	   					   fPPOEmissionSpectra; }


	 sized_array<G4double>&    GetPCRefractionIndex() {return
	   					   fPCRefractionIndex; }
	 sized_array<G4double>&    GetScintRefractionIndex() {return
	   					   fScintRefractionIndex; }
	 sized_array<G4double>&    GetWaterRefractionIndex() {return
	   					   fWaterRefractionIndex; }
	 sized_array<G4double>&    GetNylonRefractionIndex() {return
	   					   fNylonRefractionIndex; }
	 sized_array<G4double>&    GetGlassRefractionIndex() {return
	   					   fGlassRefractionIndex; }


	 sized_array<G4double>&    GetPCReemissionProb() {return
	   				   fPCReemissionProbability; }
	 sized_array<G4double>&    GetPPOReemissionProb() {return
	   				   fPPOReemissionProbability; }
	 sized_array<G4double>&    GetPMTQuantumEfficiency() {return    fPMTQuantumEfficiency; }

	 sized_array<G4double>&    GetPhotonEnergyPCAb() {return
	   				   fPhotonEnergyPCAb; }
	 sized_array<G4double>&    GetPhotonEnergyPPOAb() {return
	   				   fPhotonEnergyPPOAb; }
	 sized_array<G4double>&    GetPhotonEnergyWaterAb() {return
	   				   fPhotonEnergyWaterAb; }
	 sized_array<G4double>&    GetPhotonEnergyNylonAb() {return
	   				   fPhotonEnergyNylonAb; }
	 sized_array<G4double>&    GetPhotonEnergyGlassAb() {return
	   				   fPhotonEnergyGlassAb; }
	 sized_array<G4double>&    GetPhotonEnergyDMPAb() {return
	   				   fPhotonEnergyPPOAb; }

	 sized_array<G4double>&    GetPhotonEnergyPCEm() {return
	   				   fPhotonEnergyPCEm; }
	 sized_array<G4double>&    GetPhotonEnergyPPOEm() {return
	   				   fPhotonEnergyPPOEm; }


	 sized_array<G4double>&    GetPhotonEnergyPCRef() {return
	   				   fPhotonEnergyPCRef; }
	 sized_array<G4double>&    GetPhotonEnergyScintRef() {return
	   				   fPhotonEnergyScintRef; }
	 sized_array<G4double>&    GetPhotonEnergyWaterRef() {return
	   				   fPhotonEnergyWaterRef; }
	 sized_array<G4double>&    GetPhotonEnergyNylonRef() {return
	   				   fPhotonEnergyNylonRef; }
	 sized_array<G4double>&    GetPhotonEnergyGlassRef() {return
	   				   fPhotonEnergyGlassRef; }

	 sized_array<G4double>&    GetPhotonEnergyPCRP() {return
	   				   fPhotonEnergyPCReemProb; }
	 sized_array<G4double>&    GetPhotonEnergyPPORP() {return
	   				   fPhotonEnergyPPOReemProb; }
	 sized_array<G4double>&    GetPhotonEnergyQE() {return fPhotonEnergyQE; }
	
	 sized_array<G4double>&    GetPhotonEnergyCopperStrutRef() {return fPhotonEnergyCopperStrutRef; }
	 
	 sized_array<G4double>&    GetPhotonEnergyCopperStrutEff() {return fPhotonEnergyCopperStrutEff; }
	 
	 sized_array<G4double>&    GetPhotonEnergyTubeSteelRef() {return fPhotonEnergyTubeSteelRef; }
	 
	 sized_array<G4double>&    GetPhotonEnergyTubeSteelEff() {return fPhotonEnergyTubeSteelEff; }
	 
	 sized_array<G4double>&    GetPhotonEnergyNylonReflectivity() {return fPhotonEnergyNylonReflectivity; }
	 
	 sized_array<G4double>&    GetPhotonEnergyNylonEff() {return fPhotonEnergyNylonEff; }
	 
	 sized_array<G4double>&    GetSSSReflectivity() {return
	   				   fSSSReflectivity;}
	 
	 
	 sized_array<G4double>&    GetPMTReflectivity() {return
	   				   fPMTReflectivity;}
	 sized_array<G4double>&    GetPMTShieldReflectivity() {return
	   				   fPMTShieldReflectivity;}
	 sized_array<G4double>&    GetPMTHousingReflectivity() {return
	   				   fPMTHousingReflectivity;}
	 sized_array<G4double>&    GetPMTringReflectivity() {return
	   				   fPMTringReflectivity;}
	 sized_array<G4double>&    GetPMTringspecLOBE() {return
	   				   fPMTringspecLOBE;}
	 sized_array<G4double>&    GetPMTringspecSPIKE() {return
	   				   fPMTringspecSPIKE;}
	 sized_array<G4double>&    GetPMTringbackSCATT() {return
	   				   fPMTringbackSCATT;}

	 sized_array<G4double>&    GetExternalLGReflectivity() {return
	   				   fExternalLGReflectivity;}
	 sized_array<G4double>&    GetInternalLGReflectivity() {return
	   				   fInternalLGReflectivity;}

	 sized_array<G4double>&    GetDerlinReflectivity() {return
	   				   fDerlinReflectivity;}
	 sized_array<G4double>&    GetDerlinEfficiency() {return
	   				   fDerlinEfficiency;}


	 sized_array<G4double>&    GetCopperStrutReflectivity() {return
	   				   fCopperStrutReflectivity;}
	 sized_array<G4double>&    GetCopperStrutEfficiency() {return
	   				   fCopperStrutEfficiency;}
	 
	 sized_array<G4double>&    GetTubeSteelReflectivity() {return
	   				   fTubeSteelReflectivity;}
	 sized_array<G4double>&    GetTubeSteelEfficiency() {return
	   				   fTubeSteelEfficiency;}
	 
	 sized_array<G4double>&    GetNylonReflectivity() {return
	   				   fNylonReflectivity;}
	 sized_array<G4double>&    GetNylonEfficiency() {return
	   				   fNylonEfficiency;}
	 
	 sized_array<G4double>&    GetExternalLGspecLOBE() {return
	   				   fExternalLGspecLOBE;}
	 sized_array<G4double>&    GetExternalLGspecSPIKE() {return
	   				   fExternalLGspecSPIKE;}
	 sized_array<G4double>&    GetExternalLGbackSCATT() {return
	   				   fExternalLGbackSCATT;}
	 sized_array<G4double>&    GetInternalLGspecLOBE() {return
	   				   fInternalLGspecLOBE;}
	 sized_array<G4double>&    GetInternalLGspecSPIKE() {return
	   				   fInternalLGspecSPIKE;}
	 sized_array<G4double>&    GetInternalLGbackSCATT() {return
	   				   fInternalLGbackSCATT;}

	 sized_array<G4double>&    GetTyvekspecLOBE() {return
	   				   fTyvekspecLOBE;}
	 sized_array<G4double>&    GetTyvekspecSPIKE() {return
	   				   fTyvekspecSPIKE;}
	 sized_array<G4double>&    GetTyvekbackSCATT() {return
	   				   fTyvekbackSCATT;}

	 sized_array<G4double>&    GetTyvekEfficiency() {return
	   				   fTyvekEfficiency;}

	 sized_array<G4double>&    GetTyvekReflectivity() {return
	   				   fTyvekReflectivity;}

	 sized_array<G4double>&    GetSSSspecLOBE() {return
	   				   fSSSspecLOBE;}
	 sized_array<G4double>&    GetSSSspecSPIKE() {return
	   				   fSSSspecSPIKE;}
	 sized_array<G4double>&    GetSSSbackSCATT() {return
	   				   fSSSbackSCATT;}

	 void SetSSSReflectivity(G4double);
	 void SetSSSspecLOBE (G4double);
	 void SetSSSspecSPIKE (G4double);
	 void SetSSSbackSCATT (G4double);
	 
	 void SetNylonReflectivity (G4double);
	 void SetPMTringReflectivity (G4double);
	 void SetPMTringspecSPIKE (G4double);
	 void SetConcentratorExternalReflectivity (G4double);
	 void SetConcentratorExternalspecSPIKE (G4double);
	 void SetConcentratorInternalReflectivity (G4double);
	 void SetConcentratorInternalspecSPIKE (G4double);
	 void SetCathodeReflectivity (G4double);
	 void SetShieldReflectivity (G4double);
	 
	 sized_array<G4double>&    GetSSSEfficiency() {return
	   				   fSSSEfficiency;}
	 sized_array<G4double>&    GetPMTEfficiency() {return
	   				   fPMTEfficiency;}
	 sized_array<G4double>&    GetPMTShieldEfficiency() {return
	   				   fPMTShieldEfficiency;}
	 sized_array<G4double>&    GetPMTHousingEfficiency() {return
	   				   fPMTHousingEfficiency;}
	 sized_array<G4double>&    GetPMTringEfficiency() {return
	   				   fPMTringEfficiency;}
	 sized_array<G4double>&    GetExternalLGEfficiency() {return
	   				   fExternalLGEfficiency;}
	 sized_array<G4double>&    GetInternalLGEfficiency() {return
	   				   fInternalLGEfficiency;}

	 sized_array<G4double>&    GetPhotonEnergySSSRef() {return
	   			   fPhotonEnergySSSRef; }
	 sized_array<G4double>&    GetPhotonEnergyPMTRef() {return
	   			   fPhotonEnergyPMTRef; }
	 sized_array<G4double>&    GetPhotonEnergyPMTShieldRef() {return
	   			   fPhotonEnergyPMTShieldRef; }
	 sized_array<G4double>&    GetPhotonEnergyPMTHousingRef() {return
	   			   fPhotonEnergyPMTHousingRef; }
	 sized_array<G4double>&    GetPhotonEnergyPMTringRef() {return
	   			   fPhotonEnergyPMTringRef; }
	 sized_array<G4double>&    GetPhotonEnergyExternalLGRef() {return
	   			   fPhotonEnergyExternalLGRef; }
	 sized_array<G4double>&    GetPhotonEnergyInternalLGRef() {return
	   			   fPhotonEnergyInternalLGRef; }

	 sized_array<G4double>&    GetPhotonEnergySSSspecLOBE() {return
	   			   fPhotonEnergySSSspecLOBE; }
	 sized_array<G4double>&    GetPhotonEnergySSSspecSPIKE() {return
	   			   fPhotonEnergySSSspecSPIKE; }
	 sized_array<G4double>&    GetPhotonEnergySSSbackSCATT() {return
	   			   fPhotonEnergySSSbackSCATT; }

	 sized_array<G4double>&    GetPhotonEnergyPMTringspecSPIKE() {return
	   			   fPhotonEnergyPMTringspecSPIKE; }
	 sized_array<G4double>&    GetPhotonEnergyPMTringspecLOBE() {return
	   			   fPhotonEnergyPMTringspecLOBE; }
	 sized_array<G4double>&    GetPhotonEnergyPMTringbackSCATT() {return
	   			   fPhotonEnergyPMTringbackSCATT; }

	 sized_array<G4double>&    GetPhotonEnergySSSEff() {return
	   			   fPhotonEnergySSSEff; }
	 sized_array<G4double>&    GetPhotonEnergyPMTEff() {return
	   			   fPhotonEnergyPMTEff; }
	 sized_array<G4double>&    GetPhotonEnergyPMTShieldEff() {return
	   			   fPhotonEnergyPMTShieldEff; }
	 sized_array<G4double>&    GetPhotonEnergyPMTHousingEff() {return
	   			   fPhotonEnergyPMTHousingEff; }
	 sized_array<G4double>&    GetPhotonEnergyPMTringEff() {return
	   			   fPhotonEnergyPMTringEff; }
	 sized_array<G4double>&    GetPhotonEnergyExternalLGEff() {return
	   			   fPhotonEnergyExternalLGEff; }
	 sized_array<G4double>&    GetPhotonEnergyInternalLGEff() {return
	   			   fPhotonEnergyInternalLGEff; }

	 sized_array<G4double>&    GetPhotonEnergyDerlinEff() {return
			           fPhotonEnergyDerlinEff; }
	 sized_array<G4double>&    GetPhotonEnergyDerlinRef() {return
	   			   fPhotonEnergyDerlinRef; }

	
	sized_array<G4double>&    GetPhotonEnergyTyvek() {return
	   			   fPhotonEnergyTyvek; }

	
	sized_array<G4int>&          GetChannelNumber() {return
				        fChannelNumber; }

        sized_array<G4int>&          GetChannelNumberOD() {return fChannelNumberOD; }


       G4double*      GetChannelRelQE() {return  fChannelRelQE; }

///It returns the Physics Table
	G4PhysicsTable*	 	GetPhysicsTable() const {return fPhysicsTable; }


///It builds the table in which the properties of the materials are stored
	void	BuildPhysicsTable();

///Getter of the singleton	
    static BxPropertyCollection* Get();

    void ReadC14CorrectionFactor();
    G4double GetC14CorrectionFactor(G4int element) {return fC14ShapeFactorCorrection[element];}

  private:
///singleton
    static BxPropertyCollection * me;
    BxPropertyCollection();
	
    sized_array<G4double> fPCAbsorptionLength;
	sized_array<G4double> fPPOAbsorptionLength;
	sized_array<G4double> fWaterAbsorptionLength;
	sized_array<G4double> fNylonAbsorptionLength;
	sized_array<G4double> fGlassAbsorptionLength;
	sized_array<G4double> fDMPAbsorptionLength;
	sized_array<G4double> fPCScatteringLength;

	sized_array<G4double> fPCEmissionSpectra;
	sized_array<G4double> fPPOEmissionSpectra;

  	sized_array<G4double> fPCRefractionIndex;
  	sized_array<G4double> fScintRefractionIndex;
  	sized_array<G4double> fWaterRefractionIndex;
  	sized_array<G4double> fNylonRefractionIndex;
  	sized_array<G4double> fGlassRefractionIndex;

	sized_array<G4double> fPCReemissionProbability;
	sized_array<G4double> fPPOReemissionProbability;
  	sized_array<G4double> fPMTQuantumEfficiency;

  	sized_array<G4double> fPhotonEnergyPCAb;
  	sized_array<G4double> fPhotonEnergyPPOAb;
  	sized_array<G4double> fPhotonEnergyWaterAb;
  	sized_array<G4double> fPhotonEnergyNylonAb;
  	sized_array<G4double> fPhotonEnergyGlassAb;
  	sized_array<G4double> fPhotonEnergyDMPAb;

  	sized_array<G4double> fPhotonEnergyPCEm;
  	sized_array<G4double> fPhotonEnergyPPOEm;

  	sized_array<G4double> fPhotonEnergyPCRef;
  	sized_array<G4double> fPhotonEnergyScintRef;
  	sized_array<G4double> fPhotonEnergyWaterRef;
  	sized_array<G4double> fPhotonEnergyNylonRef;
  	sized_array<G4double> fPhotonEnergyGlassRef;

  	sized_array<G4double> fPhotonEnergyPCReemProb;
  	sized_array<G4double> fPhotonEnergyPPOReemProb;
  	sized_array<G4double> fPhotonEnergyQE;
  	sized_array<G4double> fPhotonEnergyTyvek;

  	sized_array<G4double> fPhotonEnergyCopperStrutRef;
  	sized_array<G4double> fPhotonEnergyCopperStrutEff;
  	sized_array<G4double> fPhotonEnergyTubeSteelRef;
  	sized_array<G4double> fPhotonEnergyTubeSteelEff;
  	sized_array<G4double> fPhotonEnergyNylonReflectivity;
  	sized_array<G4double> fPhotonEnergyNylonEff;
    	
	sized_array<G4double> fSSSReflectivity;
    	sized_array<G4double> fPMTReflectivity;
    	sized_array<G4double> fPMTShieldReflectivity;
   	sized_array<G4double> fPMTHousingReflectivity;
   	sized_array<G4double> fPMTringReflectivity;

   	sized_array<G4double> fPMTringspecSPIKE;
   	sized_array<G4double> fPMTringspecLOBE;
   	sized_array<G4double> fPMTringbackSCATT;
    	
	sized_array<G4double> fTyvekReflectivity;
    	sized_array<G4double> fTyvekEfficiency;
	sized_array<G4double> fTyvekspecLOBE;
	sized_array<G4double> fTyvekspecSPIKE;
	sized_array<G4double> fTyvekbackSCATT;

	sized_array<G4double> fCopperStrutReflectivity;
	sized_array<G4double> fCopperStrutEfficiency;
	sized_array<G4double> fTubeSteelReflectivity;
	sized_array<G4double> fTubeSteelEfficiency;
	sized_array<G4double> fNylonReflectivity;
	sized_array<G4double> fNylonEfficiency;

	sized_array<G4double> fExternalLGspecLOBE;
	sized_array<G4double> fExternalLGspecSPIKE;
	sized_array<G4double> fExternalLGbackSCATT;
	sized_array<G4double> fInternalLGspecLOBE;
	sized_array<G4double> fInternalLGspecSPIKE;
	sized_array<G4double> fInternalLGbackSCATT;
	
	sized_array<G4double> fExternalLGReflectivity;
	sized_array<G4double> fInternalLGReflectivity;

	sized_array<G4double> fSSSspecLOBE;
	sized_array<G4double> fSSSspecSPIKE;
	sized_array<G4double> fSSSbackSCATT;

	sized_array<G4double> fSSSEfficiency;
    	sized_array<G4double> fPMTEfficiency;
    	sized_array<G4double> fPMTShieldEfficiency;
   	sized_array<G4double> fPMTHousingEfficiency;
   	sized_array<G4double> fPMTringEfficiency;

    	sized_array<G4double> fExternalLGEfficiency;
    	sized_array<G4double> fInternalLGEfficiency;

    	sized_array<G4double> fDerlinReflectivity;
    	sized_array<G4double> fDerlinEfficiency;

    	sized_array<G4double> fPhotonEnergySSSRef;
    	sized_array<G4double> fPhotonEnergyPMTRef;
    	sized_array<G4double> fPhotonEnergyPMTShieldRef;
    	sized_array<G4double> fPhotonEnergyPMTHousingRef;
    	sized_array<G4double> fPhotonEnergyPMTringRef;

    	sized_array<G4double> fPhotonEnergyPMTringspecSPIKE;
    	sized_array<G4double> fPhotonEnergyPMTringspecLOBE;
    	sized_array<G4double> fPhotonEnergyPMTringbackSCATT;
    	
	sized_array<G4double> fPhotonEnergyExternalLGRef;
    	sized_array<G4double> fPhotonEnergyInternalLGRef;

	sized_array<G4double> fPhotonEnergySSSspecLOBE;
	sized_array<G4double> fPhotonEnergySSSspecSPIKE;
	sized_array<G4double> fPhotonEnergySSSbackSCATT;

	sized_array<G4double> fPhotonEnergySSSEff;
    	sized_array<G4double> fPhotonEnergyPMTEff;
    	sized_array<G4double> fPhotonEnergyPMTShieldEff;
  	sized_array<G4double> fPhotonEnergyPMTHousingEff;
  	sized_array<G4double> fPhotonEnergyPMTringEff;

    	sized_array<G4double> fPhotonEnergyExternalLGEff;
    	sized_array<G4double> fPhotonEnergyInternalLGEff;

    	sized_array<G4double> fPhotonEnergyDerlinRef;
   	sized_array<G4double> fPhotonEnergyDerlinEff;

        sized_array<G4int>   fChannelNumber;
        sized_array<G4int>   fChannelNumberOD;
        sized_array<G4double>   fChannelRelQE;
	

	G4PhysicsTable*	 fPhysicsTable;

	G4double  thePhotonWaveLength;
	G4double  thePhotonEnergy;
	G4double  thePCAbsorptionLength;	
	G4double  thePCScatteringLength;	
  	G4double  thePPOAbsorptionLength;
	G4double  theDMPAbsorptionLength;
	G4double  theNylonAbsorptionLength;
        G4double  theScintRefractionIndex;
        G4double  theGlassRefractionIndex;
    
	G4double fC14ShapeFactorCorrection [401];
};


#endif
