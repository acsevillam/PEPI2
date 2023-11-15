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
// - Luca Brombal, INFN - Trieste - 20.06.2020
// - Camilo Sevilla, EAFIT - Medellin - 14.11.2023
// - Cristian Tibambre, UniAndes - Bogota D.C. - 14.11.2023
/////////////////////////////////////////////////////
//
/// \file RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "G4THitsMap.hh"
#include "globals.hh"

class G4Run;

/// Run action class

namespace PEPI2
{

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    ~RunAction() override = default;

    G4Run* GenerateRun() override;
    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    // - Get all HitsMaps
    std::vector<G4THitsMap<G4double>*> GetHCofThisRun() const {return fRunMaps;};
    G4double GetObjectDose(G4String);

  private:
    G4long  fSeed;
    G4int   fLuxury;
    G4int   fSaveRndm;

    std::vector<G4THitsMap<G4double>*> fRunMaps;

    G4double fObjDose;
    G4double fObjMGDose;

    std::vector<G4String> fSDName; 

};

}

#endif // RunAction_h

