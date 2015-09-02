//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: BxStackingRDMChain.cc,v 1.2 2015/01/16 10:25:57 acaminata Exp $
// GEANT4 tag $Name:  $
//

#include "BxReadParameters.hh"
#include "BxOutputVertex.hh"
//#include "BxGenebReader.hh"
#include "BxDataCollection.hh"
#include "BxStackingRDMChainMessenger.hh"
#include "BxStackingRDMChain.hh"
#include "BxLogger.hh"
#include "G4StackManager.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4ParticleTable.hh"
#include "G4NeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4OpticalPhoton.hh" 
#include "G4GenericIon.hh" 
#include "G4Ions.hh" 
#include "globals.hh"
#include "CLHEP/Random/RandExponential.h"
#include "G4SystemOfUnits.hh"
using namespace std;

BxStackingRDMChain::BxStackingRDMChain() {
  changedParticle = new G4ParticleChange;
   IsShort = false ;
   BxReadParameters::Get()->SetRDMDecay(true);
   BxReadParameters::Get()->SetRDMChain(true);
   IsReclassify = false ;
   fRadius = BxReadParameters::Get()->GetVesselExternalRadius();//425.*cm ;
   fMaxLifeTime = 10*24*3600*s; //10 days
   fVthick = BxReadParameters::Get()->GetVesselThickness();//125/1000000.*m;
   VSloaded = false;
   ffRhocut = false;
   ffZcut = false;
   frhoCut = 10*m;
   fzCut = 10*m;
   fMessenger = new BxStackingRDMChainMessenger(this);
   fNSequence = 0;
   isFirst=true;
   BxLog(routine) << "RDMChain Stacking Methode Active" << endlog ;   
}


BxStackingRDMChain::~BxStackingRDMChain(){;}

G4ClassificationOfNewTrack BxStackingRDMChain::BxClassifyNewTrack (const G4Track* aTrack) {
//Track photons
  if (aTrack->GetDefinition()->GetPDGEncoding() == 50) return fUrgent; 

if((aTrack->GetDefinition()->GetParticleName()!="opticalphoton" && aTrack->GetParentID() != 0) || (aTrack->GetDefinition()->GetAtomicNumber () > 2 )){
   (const_cast<G4Track *>(aTrack))->SetGlobalTime( 0. );
   (const_cast<G4Track *>(aTrack))->SetLocalTime( 0. );
   (const_cast<G4Track *>(aTrack))->SetProperTime( 0. );
}
  if(abs(aTrack->GetParentID()) == 1 && aTrack->GetProperTime() == 0.0 && aTrack->GetDefinition()->GetPDGEncoding() != 50) {
	  BxLog(development) 
	  <<   aTrack->GetDefinition()->GetParticleName() << " "
	  <<  " gtime: " << aTrack->GetGlobalTime() << " "
	  <<  " ptime: " << aTrack->GetProperTime() << " "
	  <<  " ltime: " << aTrack->GetLocalTime()/ns  << " "
	  <<  " steps: " << aTrack->GetCurrentStepNumber() << " "
	  <<  " E: " <<aTrack->GetKineticEnergy()/eV<< " "
	  <<  " ID: " <<aTrack->GetTrackID() << " "
	  <<  " PID: " <<aTrack->GetParentID() << " "
	  << endlog ;
   }
  G4ClassificationOfNewTrack status = fUrgent;
  if(ffRhocut){
    zz = aTrack->GetPosition().z();  
    rr = aTrack->GetPosition().mag();
    if(sqrt(rr*rr - zz*zz) > frhoCut) return fKill;
    
  }
  
  if(ffZcut){
    zz = aTrack->GetPosition().z();  
    if(zz > fzCut) return fKill;
  }

//TODO  
  if(fSupDef){
G4cout<<"CUT TO BE IMPLEMENTED "<<G4endl;
/*
    if(aTrack->GetParentID()==0){
      rr = aTrack->GetPosition().mag();
      zz = aTrack->GetPosition().z(); 
      if (VSloaded == false){ 
         BxReadParameters::Get()->ReadVesselGeometry();
         RZDim = G4int((BxReadParameters::Get()->GetZDVessel()).size()); 
	 VSloaded = true;
      }   
        
      G4double Rvessel = 0;
      G4double Rhovessel = 0; 
      for(G4int i = 0; i < RZDim; i++) { //il file parte da z massimo e scende
          if(zz < (BxReadParameters::Get()->GetZDVessel())[i] && zz > (BxReadParameters::Get()->GetZDVessel())[i+1]){
          G4double dr = (BxReadParameters::Get()->GetRDVessel())[i+1]-(BxReadParameters::Get()->GetRDVessel())[i];
          G4double dz = (BxReadParameters::Get()->GetZDVessel())[i+1]-(BxReadParameters::Get()->GetZDVessel())[i];
          Rhovessel = (BxReadParameters::Get()->GetRDVessel())[i] + dr/dz*(zz-(BxReadParameters::Get()->GetZDVessel())[i]);
          Rvessel = pow(Rhovessel*Rhovessel + zz*zz,0.5);
	  break;
        } 
      }
      if(rr > Rvessel-(1/1000000.)*m || rr < Rvessel-fVthick ) { 
        return fKill ;
      }
      
    } 
 */ }
  
 //TODO 
  if(fBufDef){
//G4cout<<"CUT TO BE IMPLEMENTED "<<G4endl;
/*
      if(aTrack->GetParentID() == 0){
      rr = aTrack->GetPosition().mag();
      zz = aTrack->GetPosition().z(); 
      G4double Rmax = 4.7*m;
      if (VSloaded == false){ 
         BxReadParameters::Get()->ReadVesselGeometry();
         RZDim = G4int((BxReadParameters::Get()->GetZDVessel()).size()); 
	 VSloaded = true;
      }  
      if(rr < Rmax){
        G4double Rvessel = 0;
        G4double Rhovessel = 0; 
        for(G4int i = 0; i < RZDim; i++) { //il file parte da z massimo e scende
          if(zz < (BxReadParameters::Get()->GetZDVessel())[i] && zz > (BxReadParameters::Get()->GetZDVessel())[i+1]){
            G4double dr = (BxReadParameters::Get()->GetRDVessel())[i+1]-(BxReadParameters::Get()->GetRDVessel())[i];
            G4double dz = (BxReadParameters::Get()->GetZDVessel())[i+1]-(BxReadParameters::Get()->GetZDVessel())[i];
            Rhovessel = (BxReadParameters::Get()->GetRDVessel())[i] + dr/dz*(zz-(BxReadParameters::Get()->GetZDVessel())[i]);
            Rvessel = pow(Rhovessel*Rhovessel + zz*zz,0.5);
	     
            break;
          }
          
        }
        if(rr < Rvessel) { return fKill ;}
      }
    }
 */ }
  
  G4Ions *ion = (G4Ions*)aTrack->GetDefinition();  

  if( aTrack->GetDefinition()->GetAtomicNumber() > 2) {
	  //Kill if the particle is stable and not excitated
	  if(aTrack->GetDefinition()->GetPDGStable() && ion->GetExcitationEnergy () == 0) {
		  BxLog(debugging)<<" Stable isotope NAME "<<aTrack->GetDefinition()->GetParticleName()<<endlog;
		  isFirst=true;
		  return fKill;
	  }

	  //Kill if the particle lives too long
	  if(aTrack->GetDefinition()->GetPDGLifeTime()>fMaxLifeTime){
		  isFirst=true;
		  return fKill;
	  }


	  if(!isDaughter && aTrack->GetParentID() == -1){ // Postponed event: now processed
		  fNSequence++;
		  isDaughter = true ;     

		  if(aTrack->GetDefinition()->GetPDGLifeTime() < 2*ms) IsShort = true ;
		  else IsShort = false ;

		  BxReadParameters::Get()->SetRealPDGMeanLife(aTrack->GetDefinition()->GetPDGLifeTime());
		  BxLog(debugging)<<" Postponed isotope NAME "<<aTrack->GetDefinition()->GetParticleName()<<endlog;
		  BxLog(debugging)<<" Postponed isotope CODE "<<aTrack->GetDefinition()->GetPDGEncoding()<<endlog;
		  BxLog(debugging)<<" Postponed isotope realPDGMeanLife "<<aTrack->GetDefinition()->GetPDGLifeTime()<<endlog;

		  G4double coinctime;
if(BxReadParameters::Get()->GetRealPDGMeanLife()>0 ){
G4double decayTime= CLHEP::RandExponential::shoot(BxReadParameters::Get()->GetRealPDGMeanLife());
coinctime=decayTime+BxReadParameters::Get()->GetPreAbsTime();
}else{
coinctime= BxReadParameters::Get()->GetPreAbsTime();
BxLog(debugging)<<"GetRealPDGMeanLife()<0"<<endlog;
}
		  BxDataCollection::Get()->SetfDTEMPOASS(coinctime);
		  BxReadParameters::Get()->SetPreAbsTime(BxDataCollection::Get()->GetfDTEMPOASS());

      BxOutputVertex::Get()->SetPDG(aTrack->GetDefinition()->GetPDGEncoding());
      BxOutputVertex::Get()->SetTime(coinctime);
      BxOutputVertex::Get()->SetNSequence(fNSequence);
      BxOutputVertex::Get()->SetPosition(aTrack->GetPosition ()/cm );
      BxOutputVertex::Get()->SetDirection(aTrack->GetMomentumDirection () );
 
      BxOutputVertex::Get()->SetIsotopeCoinc(fNSequence) ;     
      //aTrack->GetDefinition()->SetPDGLifeTime(0*s);
      //theDataCollection->SetPosition(aTrack->GetPosition());
      //theDataCollection->SetDirection(aTrack->GetMomentumDirection () );
      return fUrgent ;
    }

    if(isDaughter && aTrack->GetParentID() != 0 ) {// Postpone the daughter of a postponed event
      isSecondDaughter = true ;
      if(aTrack->GetDefinition()->GetPDGLifeTime() < fMaxLifeTime) {
	      if(aTrack->GetDefinition()->GetPDGLifeTime() < 0) 
		      return fUrgent;        
	      return fPostpone;//Prepare next isotope in chain 
      }  else {
	      isFirst=true;
	      return fKill ; //Stop the chain 	
      }
    }

    if(isSecondDaughter){ 
	    return fKill;} 

/*
    if(isSecondDaughter  || isDaughter   )       
      if(aTrack->GetParentID() == 0) return fKill ; //Kill the head chain isotope if processing a daughter 
  */                                                  // g4bx shoots it at every event, even if a postponed event is shooted

    if(aTrack->GetParentID() == 0)     {  //Process the head chain isotope 
	    if(isFirst==true){ 
		    fNSequence = 0;

		    if(aTrack->GetDefinition()->GetPDGLifeTime() < 2*ms) IsShort = true ;
		    else IsShort = false ;
      
		    BxReadParameters::Get()->SetRealPDGMeanLife(aTrack->GetDefinition()->GetPDGLifeTime());
		    BxDataCollection::Get()->SetfDTEMPOASS(0);
		    BxReadParameters::Get()->SetPreAbsTime(BxDataCollection::Get()->GetfDTEMPOASS());
		    BxLog(debugging)<<"Head isotope realPDGMeanLife "<<aTrack->GetDefinition()->GetPDGLifeTime()<<endlog;
		  BxLog(debugging)<<" Head isotope CODE "<<aTrack->GetDefinition()->GetPDGEncoding()<<endlog;
		  BxLog(debugging)<<" Head isotope NAME "<<aTrack->GetDefinition()->GetParticleName()<<endlog;

		    aTrack->GetDefinition()->SetPDGLifeTime(0*s);

		    BxOutputVertex::Get()->SetPosition(aTrack->GetPosition ()/cm );
		    BxOutputVertex::Get()->SetDirection(aTrack->GetMomentumDirection () );
		    BxOutputVertex::Get()->SetPDG(aTrack->GetDefinition()->GetPDGEncoding());
		    BxOutputVertex::Get()->SetTime(0);
		    BxOutputVertex::Get()->SetNSequence(fNSequence);  
		    BxOutputVertex::Get()->SetIsotopeCoinc(0) ;
		    isFirst=false;
		    return fUrgent ;
	    }
	    else{
		    return fKill;
	    }
    }  else if (aTrack->GetParentID() > 0)  {

      if(aTrack->GetDefinition()->GetPDGLifeTime() < fMaxLifeTime)  {
	return fPostpone;	//Prepare next isotope in chain 				        
      } else  {
	      isFirst=true;	
	      return fKill ;         //Stop the chain 	
      }

    }
  }
  
  if (aTrack->GetDefinition() == G4NeutrinoE::NeutrinoE())	       return fKill; 
  if (aTrack->GetDefinition() == G4NeutrinoMu::NeutrinoMu())	       return fKill; 
  if (aTrack->GetDefinition() == G4NeutrinoTau::NeutrinoTau())         return fKill; 
  if (aTrack->GetDefinition() == G4AntiNeutrinoE::AntiNeutrinoE())     return fKill; 
  if (aTrack->GetDefinition() == G4AntiNeutrinoMu::AntiNeutrinoMu())   return fKill; 
  if (aTrack->GetDefinition() == G4AntiNeutrinoTau::AntiNeutrinoTau()) return fKill; 

  if(aTrack->GetCurrentStepNumber() == 0 
		  &&  aTrack->GetDefinition()->GetPDGEncoding() != 50 
		  &&  abs(aTrack->GetDefinition()->GetPDGEncoding()) != 12 
		  // && abs(aTrack->GetParentID()) == 1
		  && aTrack->GetProperTime () == 0) {

	  BxLog(trace)  <<   aTrack->GetDefinition()->GetParticleName() << " "
		  <<  " gtime: " << aTrack->GetGlobalTime() << " "
		  <<  " ltime: " << aTrack->GetLocalTime()/ns  << " "
		  <<  " steps: " << aTrack->GetCurrentStepNumber() << " "
		  <<  " E: " <<aTrack->GetKineticEnergy()/eV<< " "
		  <<  " ID: " <<aTrack->GetTrackID() << " "
		  <<  " PID: " <<aTrack->GetParentID() << " "
		  <<  " Proper time: " <<aTrack->GetProperTime ()
		  << endlog ;

	  BxOutputVertex::Get()->SetDId(fDaughters);		
	  BxOutputVertex::Get()->SetDPDG(aTrack->GetDefinition()->GetPDGEncoding());	     
	  BxOutputVertex::Get()->SetDTime(0.0);	     
	  BxOutputVertex::Get()->SetDEnergy(aTrack->GetKineticEnergy ()/MeV);
	  BxOutputVertex::Get()->SetDPosition(aTrack->GetPosition ()/cm );
	  BxOutputVertex::Get()->SetDDirection(aTrack->GetMomentumDirection());
	  BxOutputVertex::Get()->SetDaughters();
	  fDaughters++;

  }

return status;
}

void BxStackingRDMChain::BxNewStage() { 
	return ;
}

void BxStackingRDMChain::BxPrepareNewEvent() {  
  fCounter         = 0;
  IsValid          = false ;
  stage            = 0;
  isDaughter       = false; 
  isSecondDaughter = false;
  fDaughters       = 0;
  fAlpha           = 0;
  fBeta            = 0;
  fGamma           = 0;
}

/*
 * $Log: BxStackingRDMChain.cc,v $
 * Revision 1.2  2015/01/16 10:25:57  acaminata
 * Stacking revised
 *
 * Revision 1.1  2014/10/20 14:53:09  marcocci
 * g4bx2 added to repository
 *
 * Revision 1.5  2014/10/10 06:21:57  acaminata
 * *** empty log message ***
 *
 * Revision 1.4  2014/09/11 10:08:43  acaminata
 * *** empty log message ***
 *
 * Revision 1.3  2014/09/09 09:12:02  acaminata
 * *** empty log message ***
 *
 * Revision 1.2  2014/06/03 13:30:23  acaminata
 * *** empty log message ***
 *
 * Revision 1.1  2014/05/03 10:33:31  marcocci
 * CMakeLists.txt
 *
 * Revision 1.13  2011/01/20 17:07:06  buizza
 * added some minor modification for z and rho cut
 *
 * Revision 1.10  2008-10-23 12:24:05  dfranco
 * fixed nsequence in postponed RDM chain events
 *
 * Revision 1.9  2007-11-07 14:10:12  dfranco
 * removed neutrinos and antineutrinos from output binary file for beta- and beta+ decays
 *
 * Revision 1.8  2007-11-05 09:56:21  dfranco
 * icoinc in the rdm chain decay is set
 *
 * Revision 1.7  2007-04-27 15:05:31  dfranco
 * RDM generators adapted to the new output
 *
 * Revision 1.6  2007-04-27 14:48:23  dfranco
 * Filling output data for each generator in order to complete the new output
 * format
 *
 * Revision 1.5  2007-04-26 08:51:25  dfranco
 * Added the structure for the new ouput format. It is not yet active!!
 *
 * Revision 1.4  2007-03-30 15:58:21  dfranco
 * Added the cvs log at the end of file
 *
 */

