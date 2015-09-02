//---------------------------------------------------------------------------//
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//                                                                           //
//                                                                           //
//                         MaGe Simulation                                   //
//                                                                           //
//      This code implementation is the intellectual property of the         //
//      MAJORANA and Gerda Collaborations. It is based on Geant4, an         //
//      intellectual property of the RD44 GEANT4 collaboration.              //
//                                                                           //
//                        *********************                              //
//                                                                           //
//    Neither the authors of this software system, nor their employing       //
//    institutes, nor the agencies providing financial support for this      //
//    work  make  any representation or  warranty, express or implied,       //
//    regarding this software system or assume any liability for its use.    //
//    By copying, distributing or modifying the Program (or any work based   //
//    on on the Program) you indicate your acceptance of this statement,     //
//    and all its terms.                                                     //
//                                                                           //
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//---------------------------------------------------------------------------//
/**
*
*---------------------------------------------------------------------------//
*
* DESCRIPTION:
*
*/
// Begin description of class here
/**
*
*
*
*/
// End class description
//
/**
* SPECIAL NOTES:
*
* This is part of the G4 rdm demo package. Original contact is P. Truscott (DERA)
*/
//
// --------------------------------------------------------------------------//
/**
* CONTACT: ydc@lbl.gov
* FIRST SUBMISSION: Wed Mar 10 14:42:07 PST 2004
*
* REVISION:Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#ifndef BxGeneratorRDMNucleus_h
#define BxGeneratorRDMNucleus_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4MGGeneratorRDMNucleus.hh
//
// Version:             0.b.3
// Date:                29/02/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "iostream"
////////////////////////////////////////////////////////////////////////////////
//
///It contains the information identifying an isotope
class BxGeneratorRDMNucleus
{
  // class description
  // The G4BxGeneratorRDMNucleus class is used to contain information identifying an
  // isotope (a,z,e)
  //
  // class description - end
public: // with description
  BxGeneratorRDMNucleus ();
///ctor  
BxGeneratorRDMNucleus (G4int a, G4int z, G4double e);
  //    Constructor defining new isotope with A,Z.E
  ///dtor
  ~BxGeneratorRDMNucleus();
  //  Destructor
  
private:
  ///mass number
  G4int a;
  ///atomic number
  G4int z;
  ///Energy
  G4double e;

  //
  //
  // INLINE DECLARATIONS/DEFINITIONS:
  //
public: // with description
  ///    Returns the value of a
  inline  G4int GetA () const {return a;}
  ///    Returns the value of z
  inline  G4int GetZ () const {return z;}
  ///    Returns the value of e
  inline  G4double GetE () const {return e;}

  //
  //
  // DECLARATIONS OF FRIENDS OF THE CLASS.
  //
  friend std::ostream &operator << (std::ostream &s, const BxGeneratorRDMNucleus &q);

};
////////////////////////////////////////////////////////////////////////////////
#endif



