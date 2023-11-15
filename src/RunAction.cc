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
/////////////////////////////////////////////////////
// Code by 
// - Luca Brombal, INFN - Trieste - 20.06.2020
// - Camilo Sevilla, EAFIT - Medellin - 14.11.2023
// - Cristian Tibambre, UniAndes - Bogota D.C. - 14.11.2023
/////////////////////////////////////////////////////
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class
//

#include "RunAction.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "PepiRun.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4AnalysisManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

namespace PEPI2
{

RunAction::RunAction()
{
  // add new units for dose
  //
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;
  const G4double picogray  = 1.e-12*gray;

  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);

  // Register accumulable to the accumulable manager
  //G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  //accumulableManager->RegisterAccumulable(fEdep);
  //accumulableManager->RegisterAccumulable(fEdep2);

  // Create analysis manager
  // The choice of the output format is done via the specified
  // file extension.
  auto analysisManager = G4AnalysisManager::Instance();

  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //

  // Creating histograms
  analysisManager->CreateH1("Eabs","Edep in absorber", 100, 0., 800*MeV);
  analysisManager->CreateH1("Egap","Edep in gap", 100, 0., 100*MeV);
  analysisManager->CreateH1("Labs","trackL in absorber", 100, 0., 1*m);
  analysisManager->CreateH1("Lgap","trackL in gap", 100, 0., 50*cm);

  // Creating ntuple
  //
  analysisManager->CreateNtuple("B4", "Edep and TrackL");
  analysisManager->CreateNtupleDColumn("Eabs");
  analysisManager->CreateNtupleDColumn("Egap");
  analysisManager->CreateNtupleDColumn("Labs");
  analysisManager->CreateNtupleDColumn("Lgap");
  analysisManager->FinishNtuple();

  fSDName.push_back(G4String("PixiRadSD"));   
  fSDName.push_back(G4String("IonChamberSD"));

}

G4Run* RunAction::GenerateRun()
{
  // - Generate new RUN object, which is specially
  // dedicated for MultiFunctionalDetector scheme.
  // Detailed description can be found in PepiRun.hh/cc.
  return new PepiRun(fSDName);
}

void RunAction::BeginOfRunAction(const G4Run* aRun)
{

  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  //G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  //accumulableManager->Reset();

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "PEPI.root";
  // Other supported output types:
  // G4String fileName = "PEPI.csv";
  // G4String fileName = "PEPI.hdf5";
  // G4String fileName = "PEPI.xml";
  analysisManager->OpenFile(fileName);
  G4cout << "Using " << analysisManager->GetType() << G4endl;
  
  if (!IsMaster()) //it is a slave, do nothing else
  {
    G4cout << ">>> Run " << std::setw(4) << aRun->GetRunID() << " starts on slave." << G4endl;
    return;
  }
 
  // - Each time a NEW RUN starts, reseed the random numbers engine
  // with a seed taken from the time in seconds from the Unix Epoch...
  fSeed = time(NULL);
  G4Random::setTheSeed(fSeed,fLuxury);
  // G4Random::setTheSeed(123456789,fLuxury);
  // G4Random::showEngineStatus();
      
  // - Save Rndm status
  if (fSaveRndm > 0) G4Random::saveEngineStatus("beginOfRun.rndm");
 
  //EventAction* eventAction = (EventAction*)(G4RunManager::GetRunManager()->GetUserEventAction());
  //eventAction->ResetDose();

  G4cout << ">>> Frame " << std::setw(4) << aRun->GetRunID() + 1 << G4endl;
  
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  auto analysisManager = G4AnalysisManager::Instance();

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
  
  // Merge accumulables
	G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	accumulableManager->Merge();

  // Merge HC
  PepiRun* pepiRun = (PepiRun*) aRun;
  pepiRun->Merge(aRun);
  
  if (!IsMaster()){
    G4cout << ">>> Run " << std::setw(4) << aRun->GetRunID() << " ends on slave." << G4endl;
    return;
  }

  // - Get accumulated quantities for this RUN.
  fRunMaps = pepiRun->GetHCofThisRun();

}

}
