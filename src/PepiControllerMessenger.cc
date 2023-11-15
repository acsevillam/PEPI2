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
/// \file PepiControllerMessenger.cc
/// \brief Implementation of the PepiControllerMessenger class

#include "PepiControllerMessenger.hh"
#include "G4SystemOfUnits.hh"

namespace PEPI2
{

PepiControllerMessenger::PepiControllerMessenger(PepiController* controller)
{
  fController = controller;

  fPepiDirectory = new G4UIdirectory( "/Pepi/");
  fPepiDirectory->SetGuidance("UI commands specific to this example.");

  fControllerDirectory = new G4UIdirectory("/Pepi/cont/");
  fControllerDirectory->SetGuidance("Controller commands");

  beamOnCmd = new G4UIcmdWithAnInteger( "/Pepi/cont/beamOn", this);
  beamOnCmd->SetGuidance("Run the Pepi program with this many events");
  beamOnCmd->SetParameterName("beamOn",false);
  beamOnCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  loadConfigCmd = new G4UIcmdWithAString( "/Pepi/cont/loadConfig", this);
  loadConfigCmd->SetGuidance("Read the list of rotation angles and translation distances from a file");
  loadConfigCmd->SetParameterName("loadConfig",false);
  loadConfigCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  setBaseNameCmd = new G4UIcmdWithAString( "/Pepi/cont/setBaseName", this);
  setBaseNameCmd->SetGuidance("Define the BASE output file name");
  setBaseNameCmd->SetGuidance("For example '../analysis/Smc_WaterPhantom'. Do not provide an extension!");
  setBaseNameCmd->SetParameterName("setBaseName",false);
  setBaseNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  setMaxDoseCmd = new G4UIcmdWithADoubleAndUnit("/Pepi/cont/setMaxDose",this);
  setMaxDoseCmd->SetGuidance("Impose a restriction on the object dose");
  setMaxDoseCmd->SetParameterName("setMaxDose",false);
  setMaxDoseCmd->SetUnitCategory("Dose");
  setMaxDoseCmd->AvailableForStates(G4State_Idle);

  setMaxCurrentCmd = new G4UIcmdWithAnInteger("/Pepi/cont/setMaxCurrent",this);
  setMaxCurrentCmd->SetGuidance("Impose a restriction on the object dose");
  setMaxCurrentCmd->SetParameterName("setMaxCurrent",false);
  setMaxCurrentCmd->AvailableForStates(G4State_Idle);

  setMaxPixCountCmd = new G4UIcmdWithAnInteger("/Pepi/cont/setMaxPixCount",this);
  setMaxPixCountCmd->SetGuidance("Impose a restriction on the object dose");
  setMaxPixCountCmd->SetParameterName("setMaxPixCount",false);
  setMaxPixCountCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PepiControllerMessenger::~PepiControllerMessenger()
{
    delete fPepiDirectory;
    delete fControllerDirectory;
    delete beamOnCmd;
    delete loadConfigCmd;
    delete setBaseNameCmd;
    delete setMaxDoseCmd;
    delete setMaxCurrentCmd;
    delete setMaxPixCountCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiControllerMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if (command == beamOnCmd)
    {
        fController->Simulate(beamOnCmd->GetNewIntValue(newValue));
    }
    else if (command == loadConfigCmd)
    {
        fController->LoadConfig(newValue);
    }
    else if (command == setBaseNameCmd)
    {
        fController->SetOutputFileName(newValue);
    }
    else if (command == setMaxDoseCmd)
    {
        fController->SetMaxDose(setMaxDoseCmd->GetNewDoubleValue(newValue));
    }
    else if (command == setMaxCurrentCmd)
    {
        fController->SetMaxCurrent(setMaxCurrentCmd->GetNewIntValue(newValue));
    }
    else if (command == setMaxPixCountCmd)
    {
        fController->SetMaxPixCount(setMaxPixCountCmd->GetNewIntValue(newValue));
    }
}

}