//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014

#ifndef BXDATACOLLECTION_H
#define BXDATACOLLECTION_H

#include "G4ThreeVector.hh"
#include <vector>

class PhotonData  {
   private:

	G4int 		fID;		   //Photon ID (number)
  	//G4double	fStartWavelength;  //Photon initial Wave Length
  	//G4double	fEndWavelength;    //Absorpted Photon Wave Length
	G4double	fTrackLength;      //Total Track Length
	G4double	fMeanStepLength;   //Mean Value of Step Length
	G4int		fNumberOfSteps;    //Number of Steps
	G4int		fAbsZoneStatus;    //At what Zone it was absorpted
	G4int		fAbsMatStatus;     //At what Material it was absorpted
  	G4int		fPMTNumber;	   //At which PMT it was absorpted
  	G4int		fPMT_OD_Number;	   //At which PMT in OUTER Detector it was absorpted
	G4double	fFlightTime;	   //Time of flight only
//	G4double	fPMTSpreadTime;	   //PMT Jitter time
//	G4double	fScintillatorTime; //Time of Emission
//	G4double	fPPOReemissionTime;//Time of PPO ReEmission
//	G4double	fPCReemissionTime; //Time of PC ReEmission
  //Addition to control the reemission processes
  //      G4int           fNumPPOReemissions;//Number of PPO ReEmissions
    //    G4int           fNumPCReemissions; //Number of PC ReEmissions
      //  G4int           fNumPCRelays;      //Number of PC RelayScatterings
  //Addition to control the reflection processes
        G4int           fNbOfPhReflections;//Number of Photon Reflections

  
  public:
	PhotonData();
        ~PhotonData() { }
	
  	void 	SetID(G4int id)		 	       {fID = id; }
 // 	void	SetStartWavelength(G4double swl)       {fStartWavelength = swl;}
//  	void	SetEndWavelength(G4double ewl)         {fEndWavelength = ewl; }
	void	SetTrackLength(G4double tl)  	       {fTrackLength = tl; }
	void	SetMeanStepLength(G4double msl)        {fMeanStepLength = msl;}
	void	SetNumberOfSteps(G4int nos)  	       {fNumberOfSteps = nos; }
	void	SetAbsZoneStatus(G4int isd)	       {fAbsZoneStatus = isd; }
	void	SetAbsMatStatus(G4int isd)	       {fAbsMatStatus = isd; }
  	void	SetPMTNumber(G4int pmtnum)  	       {fPMTNumber = pmtnum; }
  	void	SetPMT_OD_Number(G4int pmtODnum)       {fPMT_OD_Number = pmtODnum; }
	void	SetFlightTime(G4double time) 	       {fFlightTime = time; }
//	void	SetPMTSpreadTime(G4double stime)       {fPMTSpreadTime = stime; }
//	void	SetScintillatorTime(G4double sctime)   {fScintillatorTime = sctime;}
//	void	SetPPOReemissionTime(G4double ppotime) {fPPOReemissionTime = ppotime;}
//	void	SetPCReemissionTime(G4double pctime)   {fPCReemissionTime = pctime;}
  //Addition to control the reemission processes
  //      void    SetNumPPOReemissions(G4int NumPPO)     {fNumPPOReemissions = NumPPO;}    
    //    void    SetNumPCReemissions(G4int NumPCreem)   {fNumPCReemissions = NumPCreem;}
      //  void    SetNumPCRelays(G4int NumPCrel)         {fNumPCRelays = NumPCrel;}
  //Addition to control the reflection processes
        void    SetNbOfPhReflections(G4int NPhRefl)     {fNbOfPhReflections = NPhRefl;}  
  
  	G4int 		GetID()			const {return fID; }
  //	G4double	GetStartWavelength() 	const {return fStartWavelength;}
 // 	G4double	GetEndWavelength()  	const {return fEndWavelength; }
	G4double	GetTrackLength()  	const {return fTrackLength; }
	G4double	GetMeanStepLength()  	const {return fMeanStepLength;}
	G4int		GetNumberOfSteps()  	const {return fNumberOfSteps; }
	G4int		GetAbsZoneStatus() 	const {return fAbsZoneStatus;}
	G4int		GetAbsMatStatus() 	const {return fAbsMatStatus;}
  	G4int		GetPMTNumber()  	const {return fPMTNumber; }
  	G4int		GetPMT_OD_Number()  	const {return fPMT_OD_Number; }	
	G4double	GetFlightTime()  	const {return fFlightTime; }
//	G4double	GetPMTSpreadTime()  	const {return fPMTSpreadTime; }
//	G4double	GetScintillatorTime()  	const {return fScintillatorTime;} 
  //      G4double        GetPPOReemissionTime()  const {return fPPOReemissionTime;}
    //    G4double        GetPCReemissionTime()   const {return fPCReemissionTime;}
  //Addition to control the reemission processes
      //  G4int           GetNumPPOReemissions()  const {return fNumPPOReemissions;}
        //G4int           GetNumPCReemissions()  	const {return fNumPCReemissions;}
        //G4int           GetNumPCRelays()  	const {return fNumPCRelays;}
  //Addition to control the reflection processes
        G4int           GetNbOfPhReflections()  const {return fNbOfPhReflections;}
};


class BxDataCollection  {

	private:
		///constructor
		BxDataCollection();

		static BxDataCollection* me;

		//For standart Event at one point
		G4int		fEID;			//Primary particle ID (number)
		G4int 		fEParticle;		//Primary particle Code (Type)

		G4ThreeVector	fEPosition;     	//Primary particle position
		G4ThreeVector	fEDirection;     	//Primary particle direction
		//	G4double 	fShootDist;             //Distance to the Detector Center from Muon Trajectory

		G4double 	fEEnergy;		//Primary particle Energy
		
	   G4int		fNumberOfPhotons;	//Number of generated Photons

	G4int		fNumberOfHitPMT;	//Total number of hit PMT
						//in INNER Detector
	G4int		fNumberOfHitPMTmu;	//Total number of hit PMT
						//in OUTER Detector
	G4int 		fVerboseFlag;

	G4double        fNeutronTime;           //Neutron time untill capture
	G4int		fNeutronCaptureFlag;    //Neutron capture flag
	
        G4double        fDTEMPOASS;             //Taken from Geneb
	std::vector<PhotonData> 	fPhData;	//Collection of Photons

   public:
	virtual ~BxDataCollection();
	static BxDataCollection* Get();
	void Clear();

//***************************Data-access**************************

//For standart Event at one point
	void	SetEID(G4int i)			{fEID = i; }
	void	SetNumberOfPhotons(G4int num)	{fNumberOfPhotons = num; }

	void	SetNumberOfHitPMT(G4int num)	{fNumberOfHitPMT = num; }
	void	SetNumberOfHitPMTmu(G4int num)	{fNumberOfHitPMTmu = num; }
	
	void 	SetPosition (G4ThreeVector p )	{fEPosition = p; }
	void 	SetDirection (G4ThreeVector d )	{fEDirection = d; }	

//	void 	SetShootDist ( G4double dist )	{fShootDist = dist; }	
	void 	SetParticle (G4int i )		{fEParticle = i; }
	void 	SetEnergy ( G4double e )	{fEEnergy = e; }
 	void 	SetVerbose (G4int v )		{fVerboseFlag = v; }
	void 	SetNeutronTime (G4double tn )	{fNeutronTime = tn; }
	void 	SetNeutronCaptureFlag (G4int ncf )	{fNeutronCaptureFlag = ncf; }
  	void 	SetfDTEMPOASS(G4double dtemp)	{fDTEMPOASS = dtemp;}

	//For standart Event at one point
	G4int 		GetEID() 		const {return fEID; }
	G4int		GetNumberOfPhotons()	const {return fNumberOfPhotons; }
	G4int		GetNumberOfHitPMT()	const {return fNumberOfHitPMT; }

	G4int		GetNumberOfHitPMTmu()	const {return fNumberOfHitPMTmu; }
	
	G4ThreeVector 	GetPosition ()		const {return fEPosition; }
	G4ThreeVector 	GetDirection ()		const {return fEDirection; }
	
//	G4double 	GetShootDist ()		const {return fShootDist; }	
	G4int 		GetParticle ()		const {return fEParticle; }
	G4double 	GetEnergy ()		const {return fEEnergy; }
     	G4int 		GetVerbose ()		const {return fVerboseFlag; }
	G4double 	GetNeutronTime ()	const {return fNeutronTime; }
	G4int   	GetNeutronCaptureFlag () 	const {return fNeutronCaptureFlag; }	
  	G4double 	GetfDTEMPOASS()		const {return fDTEMPOASS;}
   	std::vector<PhotonData>&	GetPhotonData()       {return fPhData;}

};

#endif
