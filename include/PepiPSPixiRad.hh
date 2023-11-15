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
/////////////////////////////////////////////////////
// Code by 
// - Akinori Kimura, Tsukasa ASO, 14-11-2005   
// - Luca Brombal, INFN - Trieste - 20.06.2020
// - Camilo Sevilla, EAFIT - Medellin - 14.11.2023
// - Cristian Tibambre, UniAndes - Bogota D.C. - 14.11.2023
/////////////////////////////////////////////////////
//
/// \file PepiPSPixiRad.hh
/// \brief Definition of the PepiPSPixiRad class

#ifndef PepiPSPixiRad_h
#define PepiPSPixiRad_h 1

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"

/// PixiRad sensitive detector class

namespace PEPI2
{

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This is a primitive scorer class for scoring energy deposit.
//
///////////////////////////////////////////////////////////////////////////////

class PepiPSPixiRad : public G4VPrimitiveScorer
{
 
 public: // with description
     /*PepiPSPixiRad(G4String name, G4double threshold, 
          G4int nPixelsX, G4int nPixelsY, 
          G4double pixelSizeX, G4double pixelSizeY, G4double pixelSizeZ,
          G4int depth=0); // default unit*/
     PepiPSPixiRad(G4String name, G4double threshold, 
          G4int nPixelsX, G4int nPixelsY, 
          G4double pixelSizeX, G4double pixelSizeY, G4double pixelSizeZ,
          const G4String& unit, G4int depth=0);
     virtual ~PepiPSPixiRad();

 protected: // with description
     virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);

 public: 
     virtual void Initialize(G4HCofThisEvent*);
     virtual void EndOfEvent(G4HCofThisEvent*);
     virtual void clear();
     virtual void DrawAll();
     virtual void PrintAll();

     virtual void SetUnit(const G4String& unit);

  private:
     G4int GetIdx(G4double, G4double);

     G4int HCID;
     G4THitsMap<G4double>* EvtMap;

     G4double fThreshold;
     G4int fnPixelsX;
     G4int fnPixelsY;
     G4double fPixelSizeX;
     G4double fPixelSizeY;
     G4double fPixelSizeZ;
     G4ThreeVector fHitpos;
     G4int fIndex;
     G4double feDep;
      
};

}

#endif // PepiPSPixiRad_h
