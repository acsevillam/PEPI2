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
/// \file PepiDetectorMessenger.cc
/// \brief Implementation of the PepiDetectorMessenger class

#include "PepiDetectorMessenger.hh"
#include "PepiDetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4Material.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace PEPI2
{

PepiDetectorMessenger::PepiDetectorMessenger(PepiDetectorConstruction* Det)
: G4UImessenger(),
  fDetectorConstruction(Det),
  fPepiDirectory(0),
  fDetDirectory(0),
  fObjectMaterialCmd(0),
  fObjectDetDistanceCmd(0),
  fMaskThicknessCmd(0),
  fM2PitchCmd(0),
  fM2ApertureCmd(0),
  fSourcePosZCmd(0),
  fRotationAngleCmd(0),
  fEISteppingCmd(0),
  fObjectRadiusCmd(0),
  //fBeamHeightCmd(0),
  fAcquisitionTypeCmd(0),
  fDetTypeCmd(0),
  fSetBidimensionalCmd(0),
  fThreshold1Cmd(0),
  fThreshold2Cmd(0),
  fCheckOverlapsCmd(0)
{
  fPepiDirectory = new G4UIdirectory("/Pepi/");
  fPepiDirectory->SetGuidance("UI commands specific to this example.");

  fDetDirectory = new G4UIdirectory("/Pepi/det/");
  fDetDirectory->SetGuidance("Detector construction control");

  G4String matList;
  const G4MaterialTable* matTbl = G4Material::GetMaterialTable();
  for(size_t i=0;i<G4Material::GetNumberOfMaterials();i++)
  {
    matList += (*matTbl)[i]->GetName();
    matList += " ";
    G4cout << "\nMatlist "<< (*matTbl)[i]->GetName() << G4endl;
  }

  fObjectMaterialCmd = new G4UIcmdWithAString("/Pepi/det/setObjMat", this);
  fObjectMaterialCmd->SetGuidance("Select the material of the object");
  fObjectMaterialCmd->SetParameterName("ObjMat",false);
  fObjectMaterialCmd->SetCandidates(matList);
  fObjectMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fObjectDetDistanceCmd = new G4UIcmdWithADoubleAndUnit("/Pepi/det/setObjectDetDistance",this);
  fObjectDetDistanceCmd->SetGuidance("Select the distance between the object and the detector. Only use it before geometry initialization");
  fObjectDetDistanceCmd->SetParameterName("objectDetDistance",false);
  fObjectDetDistanceCmd->SetUnitCategory("Length");
  fObjectDetDistanceCmd->AvailableForStates(G4State_PreInit);
// new for PEPI
  fSourcePosZCmd = new G4UIcmdWithADoubleAndUnit("/Pepi/det/setSourcePosZ",this);
  fSourcePosZCmd->SetGuidance("Select the size of the source (round). Only use it before geometry initialization");
  fSourcePosZCmd->SetParameterName("sourcePosZ",false);
  fSourcePosZCmd->SetUnitCategory("Length");
  fSourcePosZCmd->AvailableForStates(G4State_PreInit);

  fMaskThicknessCmd = new G4UIcmdWithADoubleAndUnit("/Pepi/det/setMaskThickness",this);
  fMaskThicknessCmd->SetGuidance("Select the thickness of the mask. Only use it before geometry initialization");
  fMaskThicknessCmd->SetParameterName("maskThickness",false);
  fMaskThicknessCmd->SetUnitCategory("Length");
  fMaskThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fM2PitchCmd = new G4UIcmdWithADoubleAndUnit("/Pepi/det/setM2Pitch",this);
  fM2PitchCmd->SetGuidance("Select the pitch of the detector mask. Only use it before geometry initialization");
  fM2PitchCmd->SetParameterName("m2Pitch",false);
  fM2PitchCmd->SetUnitCategory("Length");
  fM2PitchCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fM2ApertureCmd = new G4UIcmdWithADoubleAndUnit("/Pepi/det/setM2Aperture",this);
  fM2ApertureCmd->SetGuidance("Select the aperture of the detector mask. Only use it before geometry initialization");
  fM2ApertureCmd->SetParameterName("m2Aperture",false);
  fM2ApertureCmd->SetUnitCategory("Length");
  fM2ApertureCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fSrcObjDistanceCmd = new G4UIcmdWithADoubleAndUnit("/Pepi/det/setSrcObjDistance",this);
  fSrcObjDistanceCmd->SetGuidance("Select the distance between the source and the object. Only use it before geometry initialization");
  fSrcObjDistanceCmd->SetParameterName("srcObjDistance",false);
  fSrcObjDistanceCmd->SetUnitCategory("Length");
  fSrcObjDistanceCmd->AvailableForStates(G4State_PreInit);

  fEISteppingCmd = new G4UIcmdWithADoubleAndUnit("/Pepi/det/setEIStepping",this);
  fEISteppingCmd->SetGuidance("Define the stepping of the object mask. Default: 0 um");
  fEISteppingCmd->SetGuidance("The conversion is automatic");
  fEISteppingCmd->SetParameterName("setEIStepping",false);
  fEISteppingCmd->SetUnitCategory("Length");
  fEISteppingCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fEIDitheringCmd = new G4UIcmdWithADoubleAndUnit("/Pepi/det/setEIDithering",this);
  fEIDitheringCmd->SetGuidance("Define the setEIDithering of the object. Default: 0 um");
  fEIDitheringCmd->SetGuidance("The conversion is automatic");
  fEIDitheringCmd->SetParameterName("setEIDithering",false);
  fEIDitheringCmd->SetUnitCategory("Length");
  fEIDitheringCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
//--------------------------------------------
  fRotationAngleCmd = new G4UIcmdWithADoubleAndUnit("/Pepi/det/setRotAngle",this);
  fRotationAngleCmd->SetGuidance("Define the rotation angle for the object. Default: 0 deg");
  fRotationAngleCmd->SetGuidance("You can use both deg and rad units. The conversion is automatic");
  fRotationAngleCmd->SetParameterName("rotAngle",false);
  fRotationAngleCmd->SetUnitCategory("Angle");
  fRotationAngleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fSetBidimensionalCmd = new G4UIcmdWithABool("/Pepi/det/setBidimensional",this);
  fSetBidimensionalCmd->SetGuidance("Activate a bidimensional acquisition. Default: OFF(false)");
  fSetBidimensionalCmd->SetParameterName("setBidimensional",false);
  fSetBidimensionalCmd->AvailableForStates(G4State_PreInit);
// new for PEPI
  fThreshold1Cmd = new G4UIcmdWithADoubleAndUnit("/Pepi/det/setThreshold1",this);
  fThreshold1Cmd->SetGuidance("Set detector's energy threshold 1. Default: 3 keV");
  fThreshold1Cmd->SetParameterName("setThreshold1",false);
  fThreshold1Cmd->SetUnitCategory("Energy");
  fThreshold1Cmd->AvailableForStates(G4State_PreInit);

  fThreshold2Cmd = new G4UIcmdWithADoubleAndUnit("/Pepi/det/setThreshold2",this);
  fThreshold2Cmd->SetGuidance("Set detector's energy threshold 2. Default: 3 keV");
  fThreshold2Cmd->SetParameterName("setThreshold2",false);
  fThreshold2Cmd->SetUnitCategory("Energy");
  fThreshold2Cmd->AvailableForStates(G4State_PreInit);

  fAcquisitionTypeCmd = new G4UIcmdWithAString("/Pepi/det/setAcquisitionType",this);
  fAcquisitionTypeCmd->SetGuidance("Set acquisition type. Available states are 'doublemask', 'singlemask', 'conventional' (i.e., no mask). Default: doublemask");
  fAcquisitionTypeCmd->SetParameterName("setAcquisitionType",false);
  fAcquisitionTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDetTypeCmd = new G4UIcmdWithAString("/Pepi/det/setDetType",this);
  fDetTypeCmd->SetGuidance("Set detector type. '0COL' for ideal photon counting, '1COL' for single threshold, '2COL' for single threshold.");
  fDetTypeCmd->SetGuidance(" Default: '0COL'");
  fDetTypeCmd->SetParameterName("setDetType",false);
  fDetTypeCmd->AvailableForStates(G4State_PreInit);

  fCheckOverlapsCmd = new G4UIcmdWithABool("/Pepi/det/setCheckOverlaps",this);
  fCheckOverlapsCmd->SetGuidance("Enable checking of overlaps between geometries. Default: OFF(false)");
  fCheckOverlapsCmd->SetGuidance("It will greatly slow down the simulation");
  fCheckOverlapsCmd->SetParameterName("setCheckOverlaps",false);
  fCheckOverlapsCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PepiDetectorMessenger::~PepiDetectorMessenger()
{
  delete fPepiDirectory;
  delete fDetDirectory;
  delete fObjectMaterialCmd;
  delete fObjectDetDistanceCmd;
  delete fSrcObjDistanceCmd;
  delete fSourcePosZCmd;
  delete fMaskThicknessCmd;
  delete fM2PitchCmd;
  delete fM2ApertureCmd;
  delete fEISteppingCmd;
  delete fEIDitheringCmd;
  delete fRotationAngleCmd;
  delete fSetBidimensionalCmd;
  delete fThreshold1Cmd;
  delete fThreshold2Cmd;
  delete fAcquisitionTypeCmd;
  delete fDetTypeCmd;
  delete fCheckOverlapsCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if (command == fObjectMaterialCmd)
  {
    fDetectorConstruction->SetObjectMaterial(newValue);
  }
  else if (command == fObjectDetDistanceCmd)
  {
    fDetectorConstruction->SetObjectDetDistance(fObjectDetDistanceCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fSrcObjDistanceCmd)
  {
    fDetectorConstruction->SetSrcObjDistance(fSrcObjDistanceCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fSourcePosZCmd)
  {
    fDetectorConstruction->SetSourcePosZ(fSourcePosZCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fMaskThicknessCmd)
  {
    fDetectorConstruction->SetMaskThickness(fMaskThicknessCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fM2PitchCmd)
  {
    fDetectorConstruction->SetM2Pitch(fM2PitchCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fM2ApertureCmd)
  {
    fDetectorConstruction->SetM2Aperture(fM2ApertureCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fEISteppingCmd)
  { 
    fDetectorConstruction->SetEIMovements(fEISteppingCmd->GetNewDoubleValue(newValue),fDetectorConstruction->GetObjectDithering(),fDetectorConstruction->GetObjectRotation());
  }
  else if (command == fEIDitheringCmd)
  { 
    fDetectorConstruction->SetEIMovements(fDetectorConstruction->GetObjectTranslation(), fEIDitheringCmd->GetNewDoubleValue(newValue),fDetectorConstruction->GetObjectRotation());
  }
  else if (command == fRotationAngleCmd)
  { 
    fDetectorConstruction->SetEIMovements(fDetectorConstruction->GetObjectTranslation(), fDetectorConstruction->GetObjectDithering(),fRotationAngleCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fSetBidimensionalCmd)
  {
    fDetectorConstruction->SetBidimensionalAcquisition(fSetBidimensionalCmd->GetNewBoolValue(newValue));
  }
  else if (command == fThreshold1Cmd)
  {
    fDetectorConstruction->SetThreshold(fThreshold1Cmd->GetNewDoubleValue(newValue), fDetectorConstruction->GetThreshold2());
  }
  else if (command == fThreshold2Cmd)
  {
    fDetectorConstruction->SetThreshold(fDetectorConstruction->GetThreshold1(), fThreshold2Cmd->GetNewDoubleValue(newValue));
  }
  else if (command == fAcquisitionTypeCmd)
  {
    fDetectorConstruction->SetAcquisitionType(newValue);
  }
  else if (command == fDetTypeCmd)
  {
    fDetectorConstruction->SetDetType(newValue);
  }
  else if (command == fCheckOverlapsCmd)
  {
    fDetectorConstruction->SetCheckOverlaps(fCheckOverlapsCmd->GetNewBoolValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String PepiDetectorMessenger::GetCurrentValue(G4UIcommand* command)
{
  G4String answer;

  if(command == fRotationAngleCmd)
  { 
    answer = fRotationAngleCmd->ConvertToString(fDetectorConstruction->GetObjectRotation());
  }
  else if(command == fEISteppingCmd)
  {
    answer = fEISteppingCmd->ConvertToString(fDetectorConstruction->GetObjectTranslation());  
  }
  else if(command == fEIDitheringCmd)
  {
    answer = fEIDitheringCmd->ConvertToString(fDetectorConstruction->GetObjectDithering());  
  }
  else if(command == fSetBidimensionalCmd)
  {
    answer = fSetBidimensionalCmd->ConvertToString(fDetectorConstruction->GetBidimensionalAcquisition());  
  }
  else if(command == fObjectDetDistanceCmd)
  {
    answer = fObjectDetDistanceCmd->ConvertToString(fDetectorConstruction->GetObjectDetDistance());  
  }
  else if(command == fSrcObjDistanceCmd)
  {
    answer = fSrcObjDistanceCmd->ConvertToString(fDetectorConstruction->GetSrcObjDistance());  
  }
  else if(command == fSourcePosZCmd)
  {
    answer = fSourcePosZCmd->ConvertToString(fDetectorConstruction->GetSourcePosZ());  
  }
  else if(command == fMaskThicknessCmd)
  {
    answer = fMaskThicknessCmd->ConvertToString(fDetectorConstruction->GetMaskThickness());  
  }
  else if(command == fM2PitchCmd)
  {
    answer = fM2PitchCmd->ConvertToString(fDetectorConstruction->GetM2Pitch());  
  }
  else if(command == fM2ApertureCmd)
  {
    answer = fM2ApertureCmd->ConvertToString(fDetectorConstruction->GetM2Aperture());  
  }
  else if(command == fObjectMaterialCmd)
  {
    answer = fDetectorConstruction->GetObjectMaterial();  
  }
  else if(command == fAcquisitionTypeCmd)
  {
    answer = fDetectorConstruction->GetAcquisitionType();  
  }
  else if(command == fDetTypeCmd)
  {
    answer = fDetectorConstruction->GetDetType();  
  }
  else if(command == fCheckOverlapsCmd)
  {
    answer = fCheckOverlapsCmd->ConvertToString(fDetectorConstruction->GetCheckOverlaps());
  }
  else if(command == fThreshold1Cmd)
  {
    answer = fThreshold1Cmd->ConvertToString(fDetectorConstruction->GetThreshold1());
  }
  else if(command == fThreshold2Cmd)
  {
    answer = fThreshold2Cmd->ConvertToString(fDetectorConstruction->GetThreshold2());
  }

  return answer;
}

}
