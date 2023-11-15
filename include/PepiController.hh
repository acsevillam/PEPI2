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
/// \file PepiController.hh
/// \brief Definition of the PepiController class

#ifndef PepiController_h
#define PepiController_h 1

#include <vector>
#include <string>
#include <map>

#include "G4THitsMap.hh"

#include "PepiControllerMessenger.hh"

class G4Run;
class G4Timer;

/// PEPI controller class

namespace PEPI2
{

class PepiControllerMessenger;
class PepiDetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PepiController
{
  public:

    PepiController(PepiDetectorConstruction*);
    ~PepiController();

    // - Start a Simulation
    void Simulate(G4int);
    // - Load Rotation Angles from an external file
    void LoadConfig(G4String);
    // - Set the Output File Name
    void SetOutputFileName(G4String);
    // - Abort run when dose is sufficient
    void SetMaxDose(G4double);
    // - Abort run when current is sufficient
    void SetMaxCurrent(G4double);
    // - Abort run when average flat field pixel count is sufficient
    void SetMaxPixCount(G4double);

  private:

    // - Rotate and translate the Object
    void EIMovements(G4int);
    // - End of Projection results
    void RunOutput();
    // - End of Simulation results
    void SimulationOutput();
    // - Save data to disk
    void SaveData(G4THitsMap<G4double>*);
    // - Start a Projection
    void BeamOn(G4int);

    // - Private data members
    G4bool  fLoaded;

    std::vector<G4double> fSx;
    std::vector<G4double> fDx;
    std::vector<G4double> fThetaAngles;
    G4String              fSubSim;
    G4int                 fSize;
    G4double              fThetaMin;
    G4double              fThetaMax;
    G4int		  fRunID;
    G4String fDoseType;

    G4double fDose;

    G4double fMaxDose;
    G4double fMaxCurrent;
    G4double fMaxPixCount;

    G4int    fNx;
    G4int    fNy;

    G4double fSumTimeReal;
    G4double fSumTimeUser;

    G4String fBaseName;
    G4String fTxtName;
    G4String fRawName;

    PepiDetectorConstruction* fDetector;
    PepiControllerMessenger*  fMessenger;
    G4Timer*                   fTimer;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif // PepiController_h
