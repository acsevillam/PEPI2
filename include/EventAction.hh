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
/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

/// Event action class

namespace PEPI2
{

class RunAction;
class PepiDetectorConstruction;

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction* runAction,PepiDetectorConstruction* detector);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    void SumDose(G4double dose) {fObjDose += dose;};
    void SumMGDose(G4double dose) {fObjMGDose += dose;};

    void ResetDose() {fObjDose = 0; fObjMGDose = 0;}; 

    void SetDoseLimit(G4double dose) {fDoseLimit = dose;};
    void SetCurrentLimit(G4double current) {fCurrentLimit = current;};
    void SetAvgPixLimit(G4double avgPix) {fDetCountsLimit = avgPix;};
    void SetAvgPixWindow(G4double window) {fWindow = window;};

    G4double GetDoseLimit() const {return fDoseLimit;};
    G4double GetCurrentLimit() const {return fCurrentLimit;};
    G4double GetAvgPixLimit() const {return fDetCountsLimit;};
    G4int    GetAvgPixWindow() const {return fWindow;};

    std::vector<G4double>* GetKinEnergyVector() const {return fKinEnergyVector;};

  private:
    RunAction* fRunAction = nullptr;

    G4double fObjDose;
    G4double fObjMGDose;

  	G4double fDoseLimit;
    G4double fCurrentLimit;
    G4double fDetCountsLimit;
    G4int    fWindow;

    std::vector<G4double>* fKinEnergyVector;

    PepiDetectorConstruction* fDetector;

};

}

#endif


