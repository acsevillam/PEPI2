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
/// \file PepiController.cc
/// \brief Implementation of the PepiController class

#include "PepiController.hh"
#include "PepiRun.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "PepiDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Timer.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <stdio.h>
#include <ios>
#include <iomanip>
#include <string>

namespace PEPI2
{

PepiController::PepiController(PepiDetectorConstruction* Det)
{
    fDetector = Det;

    fBaseName = "image";
    //fTxtName  = fBaseName + ".txt";
    fRawName  = fBaseName + ".raw";
    fSubSim = "";

    fThetaMin = 0.*deg;
    fThetaMax = 0.*deg;

    // - Add a few "Dose" Units. By defaut Geant4 only
    // defines 1 Gray. Addining some more
    const G4double milliGray = 1.E-3  * gray;
    const G4double microGray = 1.E-6  * gray;
    const G4double nanoGray  = 1.E-9  * gray;  
    const G4double picoGray  = 1.E-12 * gray;

    new G4UnitDefinition("milligray", "mGy", "Dose", milliGray);
    new G4UnitDefinition("microgray", "uGy", "Dose", microGray);
    new G4UnitDefinition("nanogray" , "nGy", "Dose", nanoGray);
    new G4UnitDefinition("picogray" , "pGy", "Dose", picoGray);

    fDoseType = "MGD";

    fDose = 0;

    fMaxDose = 0;
    fMaxCurrent = 0;
    fMaxPixCount = 0;

    // - Add some bigger "Time" Units. By defaut Geant4 only
    // defines up to seconds. Our simulations will run for longer
    const G4double minute = 60      * s;
    const G4double hour   = 3600    * s;
    const G4double day    = 24*3600 * s;  

    new G4UnitDefinition("minutes", "min"  , "Time", minute);
    new G4UnitDefinition("hours"  , "hours", "Time", hour);
    new G4UnitDefinition("days"   , "days" , "Time", day);

    // - Set a timer to measure elapsed time. The G4Timer class
    // uses <sys/times.h> and <unistd.h>
    fTimer = new G4Timer;
    
    fMessenger = new PepiControllerMessenger(this);
    //fRunID;
}

PepiController::~PepiController()
{
  delete fTimer;
  delete fMessenger;
  fSx.clear();
  fDx.clear();
  fThetaAngles.clear();
}

void PepiController::Simulate(G4int nEvents)
{
  fSumTimeReal = 0;
  fSumTimeUser = 0;

  if (fSize > 0)
  {
    for (G4int i = 0; i < fSize; i++)
    {
      EIMovements(i);
      BeamOn(nEvents);
      RunOutput();
    }
  }
  else 
  {
    G4cout << "Missing or empty projections file." << G4endl
           << "reverting to single projection"     << G4endl;

    BeamOn(nEvents);
    RunOutput();
  }

  SimulationOutput();
}

void PepiController::EIMovements(G4int position)
{
  // Rotate the object
  fDetector->SetEIMovements(fSx[position],fDx[position],fThetaAngles[position]);
}

void PepiController::BeamOn(G4int nEvents)
{
  fTimer->Start();
  G4RunManager::GetRunManager()->BeamOn(nEvents);
  fTimer->Stop();
}

void PepiController::RunOutput()
{
  // - Get RunID() and NumberOfEvent() of the current RUN
  PepiRun* pepiRun = (PepiRun*)(G4RunManager::GetRunManager()->GetCurrentRun());

  G4int runID = pepiRun->GetRunID();
  G4int nEvents = pepiRun->GetNumberOfEvent();
  fRunID=runID;

  // - Get the Dose deposited in the Object from UserRunAction and accumulate
  RunAction* runAction = (RunAction*)(G4RunManager::GetRunManager()->GetUserRunAction());
  
  /*G4double runDose = pepiRunAction->GetObjectDose(fDoseType);
  fDose += runDose; 
*/

   //G4THitsMap<G4double>* detHits;
  G4double* current;

  std::vector<G4THitsMap<G4double>*> runMaps = runAction->GetHCofThisRun();

  // - Number of HitsMap in this RUN.
  G4int nMaps = runMaps.size();
  // Loop over HitMaps
  for (G4int i = 0; i < nMaps; i++)
  {
    G4THitsMap<G4double>* runMap = runMaps[i];
    if (runMap)
    {
      if (runMap->GetName() == "TotalSurfCurrent")
      {
        SaveData(runMap);

      }
      else if (runMap->GetName() == "Threshold1")
      {
        SaveData(runMap);
      }
      else if (runMap->GetName() == "Threshold2")
      {
        SaveData(runMap);
      }
      else if (runMap->GetName() == "CurrentIoC")
      {
        current = (*runMap)[0] ? (*runMap)[0] : new G4double(0);
      }
    }
  }

  /*// - Get the Primary Generator Action and extract the particle energy
  PepiPrimaryGeneratorAction* pepiPrimGenAction 
  = (PepiPrimaryGeneratorAction*)(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

  G4double pEnergy = (pepiPrimGenAction->GetParticleGun())->GetParticleEnergy();

  std::ofstream txtFile;
  G4String txtName = "MGDt_G0_E" + G4UIcommand::ConvertToString(pEnergy/keV) + "keV.txt";
  txtFile.open(txtName, std::ofstream::app);
  txtFile << runDose/(1.E-6*gray) << "\n";
  txtFile.close();*/

  G4double runTimeUser = fTimer->GetUserElapsed();

  fSumTimeUser += runTimeUser;

  G4double timeLeft = (fSize - (runID+1)) * runTimeUser;

  G4cout                                                                                        << G4endl
    << "|---### End of Frame " << std::setw(4) << runID+1 << " of " << fSize << " ###---"  << G4endl
    << "|"                                                                                      << G4endl
    << "| The sample was at position: "   << (fDetector->GetObjectTranslation()+fDetector->GetObjectDithering())/um  << " um" << G4endl
    << "| The sample mask was at position: "   << (fDetector->GetObjectTranslation())/um  << " um" << G4endl
    << "| The sample was at angle: "   << (fDetector->GetObjectRotation())/deg  << " deg" << G4endl
    << "| The Run consisted of: "   << nEvents << " Events"                                     << G4endl
    << "| User elapsed time is: "   << G4BestUnit(fSumTimeUser * 1.E+9, "Time")                 << G4endl
    << "| User time remaining is: " << G4BestUnit(timeLeft     * 1.E+9, "Time")                 << G4endl
    << "|"                                                                                      << G4endl    
//    << "| The Dose in the Object for this run is: "   << G4BestUnit(runDose, "Dose")            << G4endl
//    << "| The Total expected dose in the Object is: " << G4BestUnit(runDose*fSize, "Dose")      << G4endl 
//    << "|"                                                                                      << G4endl
    << "| The photon current per unit area (mm2) on the Object is: " << *current                << G4endl
    << "|"                                                                                      << G4endl    
//    << "| Saving data in: '" << fTxtName << "'"                                                 << G4endl
    << "|     Saving data in: '" << fRawName << "'"                                                 << G4endl
    << "|"                                                                                      << G4endl
    << "|---### End of Frame " << std::setw(4) << runID+1 << " of " << fSize << " ###---"  << G4endl
    << G4endl;
}

void PepiController::SimulationOutput()
{
  // - This ROUTINE prints to screen a summary of the
  // simulation
  G4cout                                                                             << G4endl
    << "|---------------### End of Simulation ###---------------"                    << G4endl
    << "|"                                                                           << G4endl
    << "| The Simulation consisted of: " << std::setw(4) << fSize << " Projections"  << G4endl
    << "| User time elapsed is: " << G4BestUnit(fSumTimeUser * 1.E+9, "Time")        << G4endl
    << "|"                                                                           << G4endl    
//    << "| The Final Dose in the Object is: " << G4BestUnit(fDose, "Dose")            << G4endl
//    << "|"                                                                           << G4endl
//    << "| Data Saved in: '" << fTxtName << "'"                                       << G4endl
//    << "|     Data Saved in: '" << fRawName << "'"                                       << G4endl
//    << "|"                                                                           << G4endl
    << "| The Frames are " << fNx << " x " << fNy << " pixels"                  << G4endl
    << "|"                                                                           << G4endl
    << "|---------------### End of Simulation ###---------------"                    << G4endl
    << G4endl;
}

void PepiController::SaveData(G4THitsMap<G4double>* runMap)
{
   // - Get the hits collection on the Detector Pixels
  // PepiRunAction* pepiRunAction = (PepiRunAction*)(G4RunManager::GetRunManager()->GetUserRunAction());

  // G4THitsMap<G4double>* detHits = pepiRunAction->GetPixelCount();

  fDetector->GetNumberOfPixelsInDetector(fNx,fNy);

  G4int pixel = 0;
  std::vector<G4float> buffer;

  std::ofstream txtFile;
  //txtFile.open(fTxtName, std::ofstream::app);
  
  SetOutputFileName(fBaseName);
  G4cout<<"\n\n\n\n\n" << fRawName << G4endl;
  FILE *rawFile = fopen(fRawName, "ab");

  for (G4int i = 0; i < fNy; i++)
  {
    for (G4int j = 0; j < fNx; j++)
    {
      G4double* nCounts = (*runMap)[pixel] ? (*runMap)[pixel] : new G4double(0);
      //txtFile  << G4float(*nCounts) << "\n";
      buffer.push_back(G4float(*nCounts));

      pixel++;
    }
    //txtFile << "\n";
  }

  //txtFile.close();

  // - Set position at the end of the file. Should NOT be needed if using
  // "a" (append) option while opening. See above
  // fseek(rawFile, 0, SEEK_END);

  if (fwrite(&buffer[0], sizeof(G4float), buffer.size(), rawFile) != buffer.size())
  {
    G4cout << "Error writing float image file!" << G4endl;
  }

  fclose(rawFile);
}

void PepiController::SetOutputFileName(G4String baseName)
{
  // - This ROUTINE build the output filename using a few
  // simulation parameters as well as a basename provided by the user
  if (!fLoaded)
  {
    G4cout << "The User should run the UI command '/Pepi/cont/loadAngles <basename>'" << G4endl
           << "before attemption to issue this command. Using default filename" << G4endl;
    return;
  }


  // - Get the Detector Construction and extract the number of x and y pixels
  fDetector->GetNumberOfPixelsInDetector(fNx,fNy);
  G4double distance = fDetector->GetObjectDetDistance();
  //G4double pixelsize = fDetector->GetBeamSizeY();
  G4double aperture = fDetector->GetM2Aperture();
  G4double pitch = fDetector->GetM2Pitch();
  G4double step = fDetector->GetObjectTranslation();
  G4double dith = fDetector->GetObjectDithering();
  //G4double ang = fDetector->GetObjectRotation();
  G4String type = fDetector->GetDetType();
  G4double th1 = fDetector->GetThreshold1();
  G4double th2 = fDetector->GetThreshold2();
  
  std::string th1_str, th2_str;
  std::ostringstream convert1, convert2;
  convert1 << std::setfill('0') << std::setw(5) << th1/eV;
  th1_str = convert1.str();
  convert2 << std::setfill('0') << std::setw(5) << th2/eV;
  th2_str = convert2.str();
  
  
  fBaseName = baseName;

  G4String sep = "_";

  //G4double maxDose = fMaxDose/(1.E-6 * gray);
  G4String fileName = fBaseName                                                     + sep +
                      "d"  + G4UIcommand::ConvertToString(distance/m)  + "m"        + sep +
                      "a"  + G4UIcommand::ConvertToString(aperture/um)  + "um"      + sep +
                      "p"  + G4UIcommand::ConvertToString(pitch/um)  + "um"         + sep +
                      "step"  + G4UIcommand::ConvertToString(step/um)               + sep +
                      "dith"  + G4UIcommand::ConvertToString(dith/um)           ;
if(type=="1COL" || type=="2COL")
{
           fileName = fileName                                        + sep +            
                      "th1_"  +  th1_str + "eV";          		       
}
if(type=="2COL")
{
           fileName = fileName                                        + sep +	      
                      "th2_"  + th2_str  + "eV";		       
}	       
           fileName = fileName + sep + G4UIcommand::ConvertToString(fRunID+1);
           		       
  fRawName = fileName + ".raw";
}

void PepiController::LoadConfig(G4String fileName)
{
  // - NEW ROUTINE (USE THIS ONE!!!!!!!) that reads a configuration file containig
  // the translation, dithering and rotation angles for all the frames
  // this version is way more robust than the previous one 
  G4cout << G4endl << "Reading configuration file..." << G4endl;

  G4String delimiter = "_";

  if (fileName.find(delimiter) != G4String::npos)
  {
    fSubSim = fileName.substr(fileName.find(delimiter)+1, 2); 
  }

  if (fLoaded) G4cout << "Another file was already loaded. Changing now..." << G4endl; 

  std::ifstream file(fileName);
  G4String line;

  while (getline(file, line))
  {
    // if (line.find("conf") == G4String::npos) continue;

    std::istringstream ss(line);

    G4double sx;
    G4double dx;
    G4double theta;

    ss >> sx >> dx >> theta;

    sx    *= um;
    dx    *= um;
    theta *= deg;

    G4cout << "Step is: " << sx/um << " um and Dither is is: " << dx/um << " um and angle is: " << theta/deg<< " deg" << G4endl;
    fSx.push_back(sx);
    fDx.push_back(dx);
    fThetaAngles.push_back(theta);

    if (!ss) continue;
  }

  file.close();

  fSize = static_cast<G4int>(fThetaAngles.size());

  // - Get min & max angles for this simulation
  if (fSize > 0)
  {
    fThetaMin = *min_element(fThetaAngles.begin(), fThetaAngles.end());
    fThetaMax = *max_element(fThetaAngles.begin(), fThetaAngles.end());
  }

  fLoaded = true;

  G4cout << "... done " << G4endl;
}

void PepiController::SetMaxDose(G4double maxDose)
{
  fMaxDose = maxDose;

  EventAction* eventAction = (EventAction*)(G4RunManager::GetRunManager()->GetUserEventAction());
  eventAction->SetDoseLimit(fMaxDose);
}

void PepiController::SetMaxCurrent(G4double maxCurrent)
{
  fMaxCurrent = maxCurrent;

  EventAction* eventAction = (EventAction*)(G4RunManager::GetRunManager()->GetUserEventAction());
  eventAction->SetCurrentLimit(fMaxCurrent);
}

void PepiController::SetMaxPixCount(G4double maxPixCount)
{
  fMaxPixCount = maxPixCount;

  EventAction* eventAction = (EventAction*)(G4RunManager::GetRunManager()->GetUserEventAction());
  eventAction->SetAvgPixLimit(fMaxPixCount);
}

}
