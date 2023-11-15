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
/// \file PepiPhysicsList.cc
/// \brief Implementation of the PepiPhysicsList class

#include "PepiPhysicsList.hh"
#include "PepiPhysicsListMessenger.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysics_option4.hh"

#include "PepiPhysicsXrayRefraction.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

namespace PEPI2
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PepiPhysicsList::PepiPhysicsList() 
: G4VModularPhysicsList()
{
  fMessenger = new PepiPhysicsListMessenger(this);

  SetVerboseLevel(-1);

  // - Build the Default physics...
  // Electromagnetic Physics
  fEmPhysicsList = new G4EmLivermorePhysics();
  // and X-Ray Refraction Physics
  fXrPhysicsList = new PepiPhysicsXrayRefraction();

  // - Electromagnetic physics options
   G4EmParameters::Instance()->SetBuildCSDARange(true);
   defaultCutValue = 1*mm;
   SetDefaultCutValue(defaultCutValue);  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PepiPhysicsList::~PepiPhysicsList()
{
	delete fEmPhysicsList;
	delete fXrPhysicsList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiPhysicsList::ConstructParticle()
{
  fEmPhysicsList->ConstructParticle();
  fXrPhysicsList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ProcessManager.hh"

void PepiPhysicsList::ConstructProcess()
{
  // - Transportation
  AddTransportation();
  
  // - Electromagnetic physics list
  fEmPhysicsList->ConstructProcess();
  
  // - X-rays Refraction physics list
  fXrPhysicsList->ConstructProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiPhysicsList::SetCuts()
{
  //G4VUserPhysicsList::SetCuts();

  if (verboseLevel > 0)
  {
    G4cout << "PepiPhysicsList::SetCuts: default cut length : "
           << G4BestUnit(defaultCutValue, "Length") << G4endl;
  }

  // - These values are used as the default production thresholds
  // for the world volume
  SetCutsWithDefault();

  // - Production thresholds for detector region
  G4Region* aRegion = G4RegionStore::GetInstance()->GetRegion("PixiRad");
  G4ProductionCuts* cuts = new G4ProductionCuts;
  cuts->SetProductionCut(5*um, "gamma");
  cuts->SetProductionCut(1*mm, "e-");
  cuts->SetProductionCut(1*mm, "e+");
  aRegion->SetProductionCuts(cuts);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiPhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel > -1) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == fEmName) return;

  if (name == "default")
  {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLivermorePhysics(); 
  } 
  else if (name == "emstandard_opt4")
  {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option4(); 
  } 
  else if (name == "emlivermore")
  {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLivermorePhysics();
  } 
  else if (name == "empenelope")
  {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmPenelopePhysics();
  }
  else if (name == "emstandard_opt0")
  {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics();
  }
  else
  {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined" << G4endl;
  }

  // G4EmParameters::Instance()->SetBuildCSDARange(true);
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
