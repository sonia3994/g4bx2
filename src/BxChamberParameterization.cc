//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// Created by D. Franco
// Revised by A. Caminata and S. Marcocci, Sept. 2014
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "BxChamberParameterization.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BxChamberParameterization::BxChamberParameterization(  
        G4int    NoChambers, 
        G4double startX,          //  X of center of first 
        G4double spacingX,        //  X spacing of centers
        G4double widthChamber, 
        G4double lengthInitial, 
        G4double lengthFinal )
{
   fNoChambers =  NoChambers; 
   fStartX     =  startX; 
   fHalfWidth  =  widthChamber*0.5;
   fSpacing    =  spacingX;
   fHalfLengthFirst = 0.5 * lengthInitial; 
   if( NoChambers > 0 ){
      fHalfLengthIncr =  0.5 * (lengthFinal-lengthInitial)/NoChambers;
      if (spacingX < widthChamber) {
       //G4Exception("BxChamberParameterization construction: Width>Spacing");
      }
   }
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BxChamberParameterization::~BxChamberParameterization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BxChamberParameterization::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  G4double      Xposition= fStartX + (copyNo+1) * fSpacing;
  G4ThreeVector origin(Xposition, 0, 0);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BxChamberParameterization::ComputeChamberDimensions
(G4Box& trackerChamber, const G4int copyNo, const G4VPhysicalVolume*) const
{
  G4double  halfLength= fHalfLengthFirst + copyNo * fHalfLengthIncr;
  trackerChamber.SetZHalfLength(halfLength);
  trackerChamber.SetYHalfLength(halfLength);
  trackerChamber.SetXHalfLength(fHalfWidth);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
