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
//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
//
//  A parameterisation that describes a series of boxes along X
//    The boxes have equal width, & their lengths are a linear equation.
//    They are spaced an equal distance apart, starting from given location.
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef BXCHAMBERPARAMETERIZATIONH_H
#define BXCHAMBERPARAMETERIZATIONH_H 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
class G4Box;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class BxChamberParameterization : public G4VPVParameterisation
{
/**This class provides the parameterisation for the tracking boxes of Opera Detector.
 *
 */
  public:
 ///Constructor 
    BxChamberParameterization(G4int    NoChambers, 
                                 G4double startX, 
                                 G4double spacing,
                                 G4double widthChamber, 
                                 G4double lengthInitial,
                                 G4double lengthFinal );

    	 ~BxChamberParameterization();
   ///It returns the traslation which gives the position of a single copy
    virtual void ComputeTransformation (const G4int copyNo,
                                G4VPhysicalVolume* physVol) const ;
    ///It attributes the spatial dimensions to the total ensemble of chambers
    virtual void ComputeChamberDimensions (G4Box & trackerLayer, const G4int copyNo,
                            const G4VPhysicalVolume* physVol) const;

  private:

    G4int    fNoChambers;   
    G4double fStartX;
    G4double fHalfWidth;    //  The half-width of each tracker chamber
    G4double fSpacing;      //  The distance between the chambers' center
    G4double fHalfLengthFirst;  //  The first half-length 
    G4double fHalfLengthIncr;   //  The Increment for the half-length 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
