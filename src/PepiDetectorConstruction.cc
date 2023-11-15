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
/// \file PepiDetectorConstruction.cc
/// \brief Implementation of the PepiDetectorConstruction class

#include "PepiDetectorConstruction.hh"
#include "PepiDetectorMessenger.hh"
//#include "PepiImageQualityPhantomParam.hh"

#include "G4RunManager.hh"

#include "G4NistManager.hh"
#include "G4MaterialPropertiesTable.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSFlatSurfaceCurrent.hh"
#include "G4PSCylinderSurfaceCurrent.hh"
#include "PepiPSPixiRad.hh"
//#include "PepiPSIoC.hh"
#include "G4SDParticleFilter.hh"
#include "G4SDParticleWithEnergyFilter.hh"

#include "G4VisAttributes.hh"
#include "G4UIcommand.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <tuple>

namespace PEPI2
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//G4ThreadLocal G4bool PepiDetectorConstruction::fConstructedSDandField = false;

PepiDetectorConstruction::PepiDetectorConstruction()
: G4VUserDetectorConstruction(),
  fConstructed(false),
  fConstructedSDandField(false),
  fBidimensional(false),
  fWorldMaterial(0),
  fDetectorMaterial(0),
  fMaskMaterial(0),
  fObjectMaterial(0),
  fObject2Material(0),
  fObject3Material(0),
  fSubMaterial(0),
  fSphereMaterial(0),
  fWorldSolid(0),
  fPixiRadSolid(0),
  fPixelSolid(0),
  fM2subSolid(0),
  fM1subSolid(0),
  fM2Solid(0),
  fM1Solid(0),
  fEnvelopeM2Solid(0),
  fEnvelopeM1Solid(0),
  fObjectSolid(0),
  fObject2Solid(0),
  fObject3Solid(0),
  fSphereSolid(0),
  fWorldLogical(0),
  fPixiRadLogical(0),
  fPixelLogical(0),
  fM2Logical(0),
  fM1Logical(0),
  fM1subLogical(0),
  fM2subLogical(0),  
  fEnvelopeM2Logical(0),
  fEnvelopeM1Logical(0),
  fObjectLogical(0),
  fObject2Logical(0),
  fObject3Logical(0),
  fSphereLogical(0),
  fWorldPhysical(0),
  fPixiRadPhysical(0),
  fPixelPhysical(0),
  fM2Physical(0),
  fM1Physical(0),
  fM1subPhysical(0),
  fM2subPhysical(0),  
  fEnvelopeM2Physical(0),
  fEnvelopeM1Physical(0),
  fObjectPhysical(0),
  fObject2Physical(0),
  fObject3Physical(0),
  fSpherePhysical(0),
  fRotAngle(0*deg),
  fTrans(0*um),
  fDith(0*um),
  fObjectDetDistance(10*cm),
  fSrcObjDistance(0.7*m),
  fWorldSizeX(0),
  fWorldSizeY(0),
  fWorldSizeZ(0),
  fOffset(50*cm),
  fnPixelsX(0),
  fnPixelsY(0),
  fPixiRadSizeX(0),
  fPixiRadSizeY(0),
  fPixiRadSizeZ(0),
  fMaskThickness(300*um),
  fM2Aperture(15*um),
  fM2Pitch(62*um),
  fSubThickness(525*um),
  fSourcePosZ(-85*cm),
  fThreshold1(3*keV),
  fThreshold2(3*keV),
  fAcquisitionType("doublemask"),
  fDetType("0COL"),
  fCheckOverlaps(false),
  fMessenger(0)
{
  // - All geometrical parameters depend on the object size
  // The object is defined as a full cylinder 
  // Inside the object there will be details of different 
  // materials 

  fWorldSizeX = 1*m;
  fWorldSizeY = 1*m;
  fWorldSizeZ = 2.3*m;

  fObjSizeR = 0.1*cm;
  fObjSizeY = 2*cm;  

  fSkinThickness = 1.45*mm;

  fDetailSizeR = 2.5*cm;

  fPixelSizeX =  62*um;
  fPixelSizeY =  62*um;
  fPixelSizeZ =  650*um;

  fMessenger = new PepiDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PepiDetectorConstruction::~PepiDetectorConstruction()
{ 
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* PepiDetectorConstruction::Construct()
{ 
    fConstructed = true;
    // - Define the Materials
    DefineMaterials();
    // - Define the Volumes
    DefineVolumes();

  return fWorldPhysical;
}

void PepiDetectorConstruction::ConstructSDandField()
{

    fConstructedSDandField = true;
    // - Define the Sensitive Detectors
    DefineDetectors();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiDetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nist = G4NistManager::Instance();

  // ========================================
  //              MATERIALS
  // ========================================

  G4Material* CdTe          = nist->FindOrBuildMaterial("G4_CADMIUM_TELLURIDE");
  G4Material* Air           = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* PlexiGlass    = nist->FindOrBuildMaterial("G4_PLEXIGLASS");
  G4Material* Water         = nist->FindOrBuildMaterial("G4_WATER");
  G4Material* Nylon  	    = nist->FindOrBuildMaterial("G4_NYLON-6-6");
  G4Material* PolyCarbonate = nist->FindOrBuildMaterial("G4_POLYCARBONATE");
  G4Material* Silicon       = nist->FindOrBuildMaterial("G4_Si");
  G4Material* Graphite      = nist->FindOrBuildMaterial("G4_GRAPHITE");

  // - Useful materials for TEST purposes...
  //G4Material* Vacuum1 = new G4Material("VACUUM1", 1., 1.01*g/mole, universe_mean_density); 
  //G4Material* Lead = nist->FindOrBuildMaterial("G4_Pb");
  G4Material* Gold = nist->FindOrBuildMaterial("G4_Au");

  // ========================================
  //         REFRACTION COEFFICIENTS
  // ========================================

  
  // Load delta coefficients and the respective energy interval
  // Values taken from: http://ts-imaging.science.unimelb.edu.au/Services/Simple/ICUtilXdata.aspx
  std::vector<double> GraphiteDelta = LoadDelta("data/Graphite_delta.txt");
  std::vector<double> SiliconDelta = LoadDelta("data/Silicon_delta.txt");
  std::vector<double> PlexiGlassDelta = LoadDelta("data/PMMA_delta.txt");
  std::vector<double> NylonDelta = LoadDelta("data/Nylon_delta.txt");
  std::vector<double> PolyCarbonateDelta = LoadDelta("data/Polycharbonate_delta.txt");
  std::vector<double> WaterDelta = LoadDelta("data/Water_delta.txt");
  std::vector<double> energies = LoadDelta("data/Energy.txt");
  	
  
  
  G4int NumEntries = static_cast<int>(energies.size());
  std::vector<double> GraphiteRindex(NumEntries);
  std::vector<double> SiliconRindex(NumEntries);
  std::vector<double> PlexiGlassRindex(NumEntries);
  std::vector<double> NylonRindex(NumEntries);
  std::vector<double> PolyCarbonateRindex(NumEntries);
  std::vector<double> WaterRindex(NumEntries);
  std::vector<double> CdTeRindex(NumEntries);
  std::vector<double> AirRindex(NumEntries);
  std::vector<double> GoldRindex(NumEntries);
  
  for (G4int i = 0; i < NumEntries; ++i)
  {
    PlexiGlassRindex[i]     = 1 - PlexiGlassDelta[i];
    GraphiteRindex[i]       = 1 - GraphiteDelta[i];
    SiliconRindex[i]        = 1 - SiliconDelta[i];

    WaterRindex[i]          = 1 - WaterDelta[i];
    NylonRindex[i]	    = 1 - NylonDelta[i];
    PolyCarbonateRindex[i]  = 1 - PolyCarbonateDelta[i];
    energies[i] = energies[i]*keV;
    CdTeRindex[i]           = 1;
    AirRindex[i]            = 1;
    GoldRindex[i]           = 1; 
 }

  G4MaterialPropertiesTable* GraphiteMatPropTbl = new G4MaterialPropertiesTable();
  GraphiteMatPropTbl->AddProperty("RINDEX",energies.data(),GraphiteRindex.data(),NumEntries);
  Graphite->SetMaterialPropertiesTable(GraphiteMatPropTbl);  

  G4MaterialPropertiesTable* SiliconMatPropTbl = new G4MaterialPropertiesTable();
  SiliconMatPropTbl->AddProperty("RINDEX",energies.data(),SiliconRindex.data(),NumEntries);
  Silicon->SetMaterialPropertiesTable(SiliconMatPropTbl);  

  G4MaterialPropertiesTable* WaterMatPropTbl = new G4MaterialPropertiesTable();
  WaterMatPropTbl->AddProperty("RINDEX",energies.data(),WaterRindex.data(),NumEntries);
  Water->SetMaterialPropertiesTable(WaterMatPropTbl);    

  G4MaterialPropertiesTable* PlexiGlassMatPropTbl = new G4MaterialPropertiesTable();
  PlexiGlassMatPropTbl->AddProperty("RINDEX",energies.data(),PlexiGlassRindex.data(),NumEntries);
  PlexiGlass->SetMaterialPropertiesTable(PlexiGlassMatPropTbl);

  G4MaterialPropertiesTable* NylonMatPropTbl = new G4MaterialPropertiesTable();
  NylonMatPropTbl->AddProperty("RINDEX",energies.data(),NylonRindex.data(),NumEntries);
  Nylon->SetMaterialPropertiesTable(NylonMatPropTbl);

  G4MaterialPropertiesTable* PolyCarbonateMatPropTbl = new G4MaterialPropertiesTable();
  PolyCarbonateMatPropTbl->AddProperty("RINDEX",energies.data(),PolyCarbonateRindex.data(),NumEntries);
  PolyCarbonate->SetMaterialPropertiesTable(PolyCarbonateMatPropTbl);

  G4MaterialPropertiesTable* CdTeMatPropTbl = new G4MaterialPropertiesTable();
  CdTeMatPropTbl->AddProperty("RINDEX",energies.data(),CdTeRindex.data(),NumEntries);
  CdTe->SetMaterialPropertiesTable(CdTeMatPropTbl);

  G4MaterialPropertiesTable* AirMatPropTbl = new G4MaterialPropertiesTable();
  AirMatPropTbl->AddProperty("RINDEX",energies.data(),AirRindex.data(),NumEntries);
  Air->SetMaterialPropertiesTable(AirMatPropTbl);

  G4MaterialPropertiesTable* GoldMatPropTbl = new G4MaterialPropertiesTable();
  GoldMatPropTbl->AddProperty("RINDEX",energies.data(),GoldRindex.data(),NumEntries);
  Gold->SetMaterialPropertiesTable(GoldMatPropTbl); 
  
  // ========================================
  //           DEFAULT MATERIALS
  // ========================================

  fWorldMaterial    = Air;// Vacuum1;
  fIonCMaterial     = Air;// Vacuum1;
  fDetectorMaterial = CdTe;
  fMaskMaterial     = Gold;
  fObjectMaterial   = PlexiGlass;
  fObject2Material  = PlexiGlass;
  fObject3Material  = PlexiGlass;
  fSubMaterial	    = Silicon;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiDetectorConstruction::DefineVolumes()
{

  // - PIXIRAD Parameters

  fnPixelsX = 512;

  // - The y-dimension of PixiRad changes if we want
  // a bidimensional acquisition or not
  if (fBidimensional)
  {
    fnPixelsY = 402;
  }
  else
  {
    fnPixelsY = 1;
  }

  G4int nx, ny;
  GetNumberOfPixelsInDetector(nx,ny);
  G4cout<<"sensor0: "<<nx<<"x"<<ny<<G4endl;

  fPixiRadSizeX = fPixelSizeX * fnPixelsX;
  fPixiRadSizeY = fPixelSizeY * fnPixelsY;
  fPixiRadSizeZ = fPixelSizeZ;

  // ========================================
  //                  WORLD
  // ========================================
  
  // - Build the WORLD as an unrotated Box in (0,0,0)
  fWorldSolid = new G4Box("World",                         //its name
                          fWorldSizeX/2,                   //its size
                          fWorldSizeY/2,
                          fWorldSizeZ/2);          
   
  fWorldLogical = new G4LogicalVolume(fWorldSolid,         //its solid
                                      fWorldMaterial,      //its material
                                      "World");            //its name
                       
  fWorldPhysical =  new G4PVPlacement(0,                   //no rotation
                                      G4ThreeVector(),     //at (0,0,0)
                                      fWorldLogical,       //its logical volume
                                      "World",             //its name
                                      0,                   //its mother  volume
                                      false,               //no boolean operation
                                      0,                   //copy number
                                      fCheckOverlaps);     //checking overlaps 


  // ========================================
  //                DETECTOR
  // ========================================

  // - Build the DETECTOR PIXEL UNIT as a Box
  fPixelSolid = new G4Box("Pixel",                          //its name
                          fPixelSizeX/2,                    //its size
                          fPixelSizeY/2,
                          fPixelSizeZ/2);
                     
  fPixelLogical = new G4LogicalVolume(fPixelSolid,          //its solid
                                      fDetectorMaterial,    //its material
                                      "PixelLV");           //its name

  // Build the Detector Envelope 
  fPixiRadSolid =  new G4Box("PixiRad",                    //its name                 
                             fPixiRadSizeX/2,              //its size
                             fPixiRadSizeY/2,
                             fPixiRadSizeZ/2);        
      
  fPixiRadLogical =  new G4LogicalVolume(fPixiRadSolid,    //its solid
                                         fWorldMaterial,   //its material
                                         "PixiRad");       //its name
                    
  // - Place the physical copies of the pixel in a x-y matrix
  // The full detector is built from the top-left corner in
  // [*][*][*][*][*][*][*][*][*][*][*][*] #1 row (fnPixelsX long)
  // [*][*][*][*]........................ #2 row (fnPixelsX long)
  // ....................................
  // .................................... # fnPixelsX * fnPixels Y
  G4int copy_no=0;  

  for (G4int iY = 0; iY < fnPixelsY ; iY++)
  {
    for (G4int iX = 0; iX < fnPixelsX ; iX++)
    {
      G4double x = + iX*fPixelSizeX - fPixiRadSizeX/2;// + fPixelSizeX/2;
      G4double y = - iY*fPixelSizeY + fPixiRadSizeY/2;// - fPixelSizeY/2;
      
      G4ThreeVector position = G4ThreeVector(x, y, 0);
      G4String  name = "Pixel_" + G4UIcommand::ConvertToString(copy_no);

      fPixelPhysical =  new G4PVPlacement(0,                           //its rotation
                                          position,                    //its position
                                          fPixelLogical,               //its logical volume
                                          name,                        //its name
                                          fPixiRadLogical,             //its mother volume
                                          false,                       //no boolean operation
                                          copy_no,                     //copy number
                                          fCheckOverlaps);             //checking overlaps 

      copy_no++;                              
    }
  }

  // - Place the Detector Envelope in the World
  G4ThreeVector positionPixirad = G4ThreeVector(0, 0, fSourcePosZ+fSrcObjDistance+fObjectDetDistance);
  fPixiRadPhysical = new G4PVPlacement(0,                                                  //its rotation
                                       positionPixirad,					   //its position
                                       fPixiRadLogical,                                    //its logical volume
                                       "PixiRad",                                          //its name
                                       fWorldLogical,                                      //its mother volume
                                       false,                                              //no boolean operation
                                       0,                                                  //copy number
                                       fCheckOverlaps);                                    //checking overlaps


  G4Region* aRegion = new G4Region("PixiRad");
  fPixiRadLogical->SetRegion(aRegion);
  aRegion->AddRootLogicalVolume(fPixiRadLogical);

  // ========================================
  //                DETECTOR MASK M2 and SUBSTRATE
  // ========================================

  if (fAcquisitionType=="doublemask")
  {
  G4double mag_M2 = (fSrcObjDistance+fObjectDetDistance)/(fSrcObjDistance+fObjectDetDistance-(fMaskThickness/2+fPixiRadSizeZ/2)); // magnification of the mask M2
  G4ThreeVector M2Position = positionPixirad-G4ThreeVector(0*cm,0,(fMaskThickness+fPixiRadSizeZ)/2+1*nm);
  // - Build the MASK APERTURE UNIT as a Box
  CreateMask("M2", mag_M2,fM2Pitch, fM2Aperture, M2Position, fMaskThickness, fMaskMaterial, fEnvelopeM2Logical, fEnvelopeM2Physical);
  CreateSubstrate("M2sub", mag_M2, M2Position-G4ThreeVector(0,0,fSubThickness/2+fMaskThickness/2), fSubThickness, fSubMaterial, fM2subLogical, fM2subPhysical);
  }
  // ========================================
  //                 Objects
  // ========================================

  // - Build the Object 1 as a Trapezoid
  fObjectSolid = new G4Trd("Trap",                         //its name
                            2.1/2*mm,                     //its half y1
                            0.5/2*mm,                     //its half y2
                            0.5*fObjSizeY,                 //its half x1
                            0.5*fObjSizeY,               //its half x2
                            fObjSizeR);                //its half height

  fObjectLogical = new G4LogicalVolume(fObjectSolid,       //its solid
                                       fObjectMaterial,    //its material
                                       "TrapLV");          //its name
  
  G4ThreeVector objectPosition = G4ThreeVector(1.5*mm, 0, fSourcePosZ+fSrcObjDistance);
  G4RotationMatrix* rotMat =  new G4RotationMatrix();
  //rotMat->rotateZ(90*deg);
  rotMat->rotateZ(fRotAngle);
  /*
  fObjectPhysical = new G4PVPlacement(rotMat,              //its rotation
                                      objectPosition,       //its translation
                                      fObjectLogical,      //its logical volume
                                      "Trap",              //its name
                                      fWorldLogical,       //its mother volume
                                      false,               //no boolean operation
                                      0,                   //copy number
                                      fCheckOverlaps);     //checking overlaps
  */
 // - Build the Object 2 as a Trapezoid
  fObject2Solid = new G4Trd("Trap2",                         //its name
                            0.82/2*mm,                             //its half x1
                            0.5/2*mm,                     //its half x2
                            0.5*fObjSizeY,                     //its half y1
                            0.5*fObjSizeY,                     //its half y2
                           fObjSizeR);                //its half height

  fObject2Logical = new G4LogicalVolume(fObject2Solid,       //its solid
                                       fObjectMaterial,
//                                       fObject2Material,		    //its material
                                       "Trap2LV");          //its name
  
  G4ThreeVector objectPosition2 = G4ThreeVector(0.*mm, 0, fSourcePosZ+fSrcObjDistance);
  /*
  fObject2Physical = new G4PVPlacement(rotMat,              //its rotation
                                      objectPosition2,       //its translation
                                      fObject2Logical,      //its logical volume
                                      "Trap2",              //its name
                                      fWorldLogical,       //its mother volume
                                      false,               //no boolean operation
                                      0,                   //copy number
                                      fCheckOverlaps);     //checking overlaps
  */
 // - Build the Object 3 as a Trapezoid
  fObject3Solid = new G4Trd("Trap3",                         //its name
                            0.66/2*mm,                             //its half x1
                            0.5/2*mm,                     //its half x2
                            0.5*fObjSizeY,                     //its half y1
                            0.5*fObjSizeY,                     //its half y2
                            fObjSizeR);                //its half height
                            
  fObject3Logical = new G4LogicalVolume(fObject3Solid,       //its solid
  					fObjectMaterial,     //its material
                                       "Trap3LV");          //its name
  
  G4ThreeVector objectPosition3 = G4ThreeVector(-1.5*mm, 0, fSourcePosZ+fSrcObjDistance);
  /*
  fObject3Physical = new G4PVPlacement(rotMat,              //its rotation
                                      objectPosition3,       //its translation
                                      fObject3Logical,      //its logical volume
                                      "Trap3",              //its name
                                      fWorldLogical,       //its mother volume
                                      false,               //no boolean operation
                                      0,                   //copy number
                                      fCheckOverlaps);     //checking overlaps
  */

  // ========================================
  //       SAMPLE MASK M1 and substrate
  // ========================================
  if (fAcquisitionType=="doublemask"|| fAcquisitionType=="singlemask")
  {
  G4double mag_M1 = (fSrcObjDistance+fObjectDetDistance)/(fSrcObjDistance-(fMaskThickness/2+fObjSizeR)); // magnification of the mask M1
  G4ThreeVector M1Position = objectPosition-G4ThreeVector(0,0,(fMaskThickness+2*fObjSizeR)/2)+G4ThreeVector(0*cm,0,0);
  
  std::tie(fEnvelopeM1Logical,fEnvelopeM1Physical) = CreateMask("M1", mag_M1,fM2Pitch, fM2Aperture, M1Position, fMaskThickness, fMaskMaterial, fEnvelopeM1Logical, fEnvelopeM1Physical);
  
  std::tie(fM1subLogical,fM1subPhysical) = CreateSubstrate("M1sub", mag_M1, M1Position-G4ThreeVector(0,0,fSubThickness/2+fMaskThickness/2), fSubThickness, fSubMaterial, fM1subLogical,fM1subPhysical);
  }

  // ========================================
  //               ION CHAMBER
  // ========================================

  // - Build the ION CHAMBER as an unrotated Box
  fIonCSolid = new G4Box("IonChamber",                     //its name
                         fPixiRadSizeX/2,                            //its size
                         fPixiRadSizeY/2,
                         2*mm);          
   
  fIonCLogical = new G4LogicalVolume(fIonCSolid,          //its solid
                                     fIonCMaterial,       //its material
                                     "IonChamberLV");     //its name
                       
  G4ThreeVector IOCposition = positionPixirad - G4ThreeVector(0, 0, 10*mm);
  
  fIonCPhysical =  new G4PVPlacement(0,                   //no rotation
                                     IOCposition,            //at position
                                     fIonCLogical,        //its logical volume
                                     "IonChamber",        //its name
                                     fWorldLogical,       //its mother  volume
                                     false,               //no boolean operation
                                     0,                   //copy number
                                     fCheckOverlaps);     //checking overlaps 
  

  
  // ========================================
  //              VISUALIZATION
  // ========================================

  G4VisAttributes* worldVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  fWorldLogical->SetVisAttributes(worldVisAtt);  
  // logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());  

  G4VisAttributes* objectVisAtt = new G4VisAttributes(G4Colour(0.6,0.8,1.0));
//  objectVisAtt->SetForceWireframe(true);
  objectVisAtt->SetForceSolid(true);
  // objectVisAtt->SetForceAuxEdgeVisible(true);     
  fObjectLogical->SetVisAttributes(objectVisAtt);

  G4VisAttributes* objectVisAtt2 = new G4VisAttributes(G4Colour(0.6,0.8,1.0));
//  objectVisAtt2->SetForceWireframe(true);
  objectVisAtt2->SetForceSolid(true);
  // objectVisAtt->SetForceAuxEdgeVisible(true);     
  fObject2Logical->SetVisAttributes(objectVisAtt2);

  G4VisAttributes* objectVisAtt3 = new G4VisAttributes(G4Colour(0.6,0.8,1.0));
//  objectVisAtt3->SetForceWireframe(true);
  objectVisAtt3->SetForceSolid(true);
  // objectVisAtt->SetForceAuxEdgeVisible(true);     
  fObject3Logical->SetVisAttributes(objectVisAtt3);

 
  G4VisAttributes* ionCVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  ionCVisAtt->SetVisibility(false);
  ionCVisAtt->SetForceWireframe(true);
  fIonCLogical->SetVisAttributes(ionCVisAtt);  
 
  G4VisAttributes* pixelVisAtt = new G4VisAttributes(G4Colour(0.0,0,1.0));
  pixelVisAtt->SetForceWireframe(true);
  //pixelVisAtt->SetForceSolid(true);
  pixelVisAtt->SetVisibility(false);
  // pixelVisAtt->SetForceAuxEdgeVisible(true);
  fPixelLogical->SetVisAttributes(pixelVisAtt);

  G4VisAttributes* pixiradVisAtt = new G4VisAttributes(G4Colour(0.0,0,1.0));
  pixiradVisAtt->SetForceWireframe(true);
  pixiradVisAtt->SetVisibility(true);
  pixiradVisAtt->SetForceSolid(true);
  // pixelVisAtt->SetForceAuxEdgeVisible(true);
  fPixiRadLogical->SetVisAttributes(pixiradVisAtt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiDetectorConstruction::DefineDetectors()
{
  // ========================================
  //                 SCORING
  // ========================================

  G4SDManager* pSDMan = G4SDManager::GetSDMpointer();
  pSDMan->SetVerboseLevel(1);
  
  // - Create 2 filters:
  // Simple particle filter --> Only gamma are detected
  // G4SDParticleFilter* gammaFilter = new G4SDParticleFilter("gammaFilter", "gamma");

  // - Particle filter with energy thresholds --> Only gamma with energy
  // between eKMin and eKMax are detected
  G4double eKMin = 1*keV;
  G4double eKMax = 150*keV;

  G4SDParticleWithEnergyFilter* gammaEKinFilter = new G4SDParticleWithEnergyFilter("gammaEKinFilter",eKMin,eKMax);
  gammaEKinFilter->add("gamma");

  // ========================================
  // - Define a Multifunctional detector 
  // ----> Ideal photon counter 
  // ========================================

  G4MultiFunctionalDetector* pixiRadSD = new G4MultiFunctionalDetector("PixiRadSD");
  pSDMan->AddNewDetector(pixiRadSD);
  SetSensitiveDetector("PixelLV",pixiRadSD);

  // - Ideal photon counter  scores the number of gammas that hit its -Z surface from outside
  // Surface is defined at the -Z surface.
  // Direction                  -Z   +Z
  //   0  IN || OUT            ->|<-  |
  //   1  IN                   ->|    |
  //   2  OUT                    |<-  |
  G4PSFlatSurfaceCurrent* sTotSurfCurrent = new G4PSFlatSurfaceCurrent("TotalSurfCurrent", 1);
  sTotSurfCurrent->SetFilter(gammaEKinFilter);
  sTotSurfCurrent->DivideByArea(false);
  sTotSurfCurrent->Weighted(false);


  pixiRadSD->RegisterPrimitive(sTotSurfCurrent);


  // -----------------------------------------------
  // - Photon counter (Pixirad) with realistic energy response and up to 2 thresholds per pixel
  if(fDetType=="1COL"|| fDetType=="2COL")
  {
      PepiPSPixiRad* sPixiRad = new PepiPSPixiRad("Threshold1", fThreshold1, 
                                                fnPixelsX, fnPixelsY,
                                                fPixelSizeX, fPixelSizeY, fPixelSizeZ, "keV");

      pixiRadSD->RegisterPrimitive(sPixiRad);
      
      if(fDetType=="2COL")
      {
      PepiPSPixiRad* sPixiRad2 = new PepiPSPixiRad("Threshold2", fThreshold2, 
                                                fnPixelsX, fnPixelsY,
                                                fPixelSizeX, fPixelSizeY, fPixelSizeZ, "keV");

      pixiRadSD->RegisterPrimitive(sPixiRad2);
      
      }
  }
  
  // ========================================  
  // - Define a Multifunctional detector 
  // ----> ION CHAMBER to check if there is flux in front of the detector 
  // ========================================

  G4MultiFunctionalDetector* ionChamberSD = new G4MultiFunctionalDetector("IonChamberSD");
  pSDMan->AddNewDetector(ionChamberSD);
  SetSensitiveDetector("IonChamberLV",ionChamberSD);

  // - Ion Chamber scores the number of photons that enter its surface
  // Surface is defined at the -Z surface.
  // Direction                  -Z   +Z
  //   0  IN || OUT            ->|<-  |
  //   1  IN                   ->|    |
  //   2  OUT                    |<-  |
  G4PSFlatSurfaceCurrent* sCurrentIoC = new G4PSFlatSurfaceCurrent("CurrentIoC", 1, "permm2");
  sCurrentIoC->SetFilter(gammaEKinFilter);
  
  ionChamberSD->RegisterPrimitive(sCurrentIoC);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiDetectorConstruction::SetObjectDetDistance(G4double objectDetDistance)
{
  if(!fConstructed)
  {
    if(objectDetDistance < 10*cm || objectDetDistance > 2*m)
    {
      G4cerr << objectDetDistance << " is out of bounds (Must be > 0.1 m AND < 2 m. - Command is ignored." << G4endl;
    }
    else
    {
      fObjectDetDistance = objectDetDistance;
    }
  }
  else
  {
    G4cerr << "Cannot change geometry after initialization. - Command is ignored." << G4endl;
    return;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// NEW for PEPI
void PepiDetectorConstruction::SetSourcePosZ(G4double sourcePosZ)
{
  if(!fConstructed)
  {
    if(sourcePosZ < -2.3/2*m)
    {
      G4cerr << sourcePosZ << "Source is probably outside the World volume - Command is ignored." << G4endl;
    }
    else
    {
      fSourcePosZ = sourcePosZ;
    }
  }
  else
  {
    G4cerr << "Cannot change geometry after initialization. - Command is ignored." << G4endl;
    return;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// NEW for PEPI
void PepiDetectorConstruction::SetSrcObjDistance(G4double srcObjDistance)
{
  if(!fConstructed)
  {
    if(srcObjDistance < 0*cm || srcObjDistance > 2.3*m)
    {
      G4cerr << srcObjDistance << "Source object distance is out of bounds (Must be > 0 m AND < World size. - Command is ignored." << G4endl;
    }
    else
    {
      fSrcObjDistance = srcObjDistance;
    }
  }
  else
  {
    G4cerr << "Cannot change geometry after initialization. - Command is ignored." << G4endl;
    return;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiDetectorConstruction::SetMaskThickness(G4double maskThickness)
{
  if(!fConstructed)
  {
    if(maskThickness < 0*um || maskThickness > 1000*um)
    {
      G4cerr << maskThickness << "Mask thickness is out of bounds (Must be > 0 um AND < 1000 um (default 300 um).- Command is ignored." << G4endl;
    }
    else
    {
      fMaskThickness = maskThickness;
    }
    G4cout <<"The mask thickness is "<< fMaskThickness/um <<"um"<< G4endl;
  }
  else
  {
    G4cerr << "Cannot change geometry after initialization. - Command is ignored." << G4endl;
    return;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiDetectorConstruction::SetM2Pitch(G4double m2Pitch)
{
  if(!fConstructed)
  {
    if(m2Pitch < 0*um || m2Pitch > 1000*um)
    {
      G4cerr << m2Pitch << "Mask pitch is out of bounds (Must be > 0 um AND < 1000 um (default 62 um).- Command is ignored." << G4endl;
    }
    else
    {
      fM2Pitch = m2Pitch;
    }
    G4cout <<"The detector mask pitch is "<< fM2Pitch/um <<"um"<< G4endl;
  }
  else
  {
    G4cerr << "Cannot change geometry after initialization. - Command is ignored." << G4endl;
    return;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiDetectorConstruction::SetM2Aperture(G4double m2Aperture)
{
  if(!fConstructed)
  {
    if(m2Aperture < 0*um || m2Aperture > fM2Pitch)
    {
      G4cerr << m2Aperture << "Mask aperture is out of bounds (Must be > 0 um AND < pitch (default 15 um).- Command is ignored." << G4endl;
    }
    else
    {
      fM2Aperture = m2Aperture;
    }
    G4cout <<"The detector mask aperture is "<< fM2Aperture/um <<"um"<< G4endl;
  }
  else
  {
    G4cerr << "Cannot change geometry after initialization. - Command is ignored." << G4endl;
    return;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiDetectorConstruction::SetThreshold(G4double threshold1, G4double threshold2)
{  
  if(!fConstructed)
  {
    fThreshold1=threshold1;
    fThreshold2=threshold2;
  }
  else
  {
    G4cerr << "Cannot change threshold after inizialization. - Command is ignored." << G4endl;
    return;
  }  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiDetectorConstruction::SetBidimensionalAcquisition(G4bool bidimensional)
{
  if(fBidimensional) return;

  fBidimensional = bidimensional;
  //G4cout<< bidimensional << G4endl;
  if(!fConstructed) return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiDetectorConstruction::SetEIMovements(G4double trans, G4double dith, G4double rotAngle)
{
  fTrans    = trans;
  fDith     = dith;
  fRotAngle = rotAngle;

  if(!fConstructed) return;
// Stepping and Dithering of the sample
 
 // object 1
  G4ThreeVector position1 = G4ThreeVector(1.5*mm+fTrans+fDith, 0, fSourcePosZ+fSrcObjDistance);
  Move("Trap", fObjectLogical, fObjectPhysical, position1);

 // object 2
  G4ThreeVector position2 = G4ThreeVector(0.*mm+fTrans+fDith, 0, fSourcePosZ+fSrcObjDistance);
  Move("Trap2", fObject2Logical, fObject2Physical, position2);

 // object 3
  G4ThreeVector position3 = G4ThreeVector(-1.5*mm+fTrans+fDith, 0, fSourcePosZ+fSrcObjDistance);
  Move("Trap3", fObject3Logical, fObject3Physical, position3);

  G4cout<<"Sample translated to " << (fTrans+fDith)/um << " um" <<G4endl;
  
  if(fAcquisitionType=="singlemask" || fAcquisitionType=="doublemask")
  {
 // Translation Sample Mask M1 and substrate
  G4double rel_mag = fSrcObjDistance/(fSrcObjDistance-(fMaskThickness+2*fObjSizeR)/2);
  G4ThreeVector position = G4ThreeVector(fTrans/rel_mag, 0, fSourcePosZ+fSrcObjDistance-(fMaskThickness+2*fObjSizeR)/2);
  Move("EnvelopeM1", fEnvelopeM1Logical, fEnvelopeM1Physical, position);
  Move("M1sub", fM1subLogical, fM1subPhysical, position-G4ThreeVector(0,0,fSubThickness/2+fMaskThickness/2));
 // print  
  G4cout<<"Sample Mask translated to " << fTrans/rel_mag/um << " um" <<G4endl;
  }

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiDetectorConstruction::SetObjectMaterial(G4String materialChoice)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);

  if(pttoMaterial){
    fObjectMaterial = pttoMaterial;
    if(fConstructed) fObjectLogical->SetMaterial(fObjectMaterial);

    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
  else G4cerr << materialChoice << " is not defined. - Command is ignored." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiDetectorConstruction::SetAcquisitionType(G4String acquisitionType)
{
  fAcquisitionType = acquisitionType;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiDetectorConstruction::SetDetType(G4String detType)
{
  fDetType = detType;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiDetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<double> PepiDetectorConstruction::LoadDelta(G4String FileName) 
{
G4cout<<"reading delta values from "<< FileName << G4endl;

  std::ifstream myData(FileName, std::ios::binary);
  std::vector<G4double> deltas;
  if(!myData.is_open())//file not open
    {
        G4cout<< FileName << "file does not exist in the specified path" << G4endl;
	return deltas;
    }
  double num = 0.0;
  while (myData >> num){
      deltas.push_back(num);
  }
  return deltas;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::tuple<G4LogicalVolume*,G4VPhysicalVolume*> PepiDetectorConstruction::CreateMask(G4String name, G4double mag, G4double pitch, G4double aperture, G4ThreeVector position, G4double thickness, G4Material* material, G4LogicalVolume* EnvelopeLogical, G4VPhysicalVolume* EnvelopePhysical)
{

  // - Build the MASK APERTURE UNIT as a Box
  G4Box* MSolid =     new G4Box(name,                               //its name
                          ((pitch-aperture)/mag)/2,   //its size
                          1.1*fPixiRadSizeY/2,
                          thickness/2);

  G4String lvname = name+"LV";                   
  G4LogicalVolume* MLogical = new G4LogicalVolume(MSolid,       //its solid
                                   material,  //its material
                                   lvname);        //its name

  // Build the MASK ENVELOPE 
  G4String envname = "Envelope"+name;                     
  G4Box* EnvelopeSolid =  new G4Box(envname,                  //its name                 
                        (1.1*fPixiRadSizeX/mag)/2,              //its size
                        1.1*fPixiRadSizeY/2,
                        thickness/2);        
  G4String lvenvname = "Envelope"+name+"LV";                         
  EnvelopeLogical =  new G4LogicalVolume(EnvelopeSolid,    //its solid
                                            fWorldMaterial,      //its material
                                            lvenvname);     //its name

  // - Place the physical copies of the mask aperture unit
   G4int copy_no=0;  

    for (G4int iX = 0; iX < int(fnPixelsX*fPixelSizeX/pitch)+2; iX++)
    {
          
      G4double x = + iX*pitch/mag - ((fPixiRadSizeX/mag)/2 + (pitch)/mag + ((pitch)/mag)/2);
      G4double y = 0;
 // G4cout << "position \n" << x <<G4endl;
      
      G4ThreeVector px_position = G4ThreeVector(x, y, 0);
      G4String  name1 = name + "_" + G4UIcommand::ConvertToString(copy_no);

      /*G4VPhysicalVolume* MPhysical =*/  new G4PVPlacement(0,                           //its rotation
                                       px_position,                    //its position
                                       MLogical,                  //its logical volume
                                       name1,                        //its name
                                       EnvelopeLogical,          //its mother volume
                                       false,                       //no boolean operation
                                       copy_no,                     //copy number
                                       fCheckOverlaps);             //checking overlaps 
      copy_no++;                              
    }

  // - Place the Sample Mask Envelope in the World
  EnvelopePhysical = new G4PVPlacement(0,                                                  //its rotation
                                          position,   				      //its position
                                          EnvelopeLogical,                                    //its logical volume
                                          envname,                                       //its name
                                          fWorldLogical,                                      //its mother volume
                                          false,                                              //no boolean operation
                                          0,                                                  //copy number
                                          fCheckOverlaps);                                    //checking overlaps


  G4Region* aRegion = new G4Region(name);
  EnvelopeLogical->SetRegion(aRegion);
  aRegion->AddRootLogicalVolume(EnvelopeLogical);
  
  

  G4VisAttributes* MVisAtt = new G4VisAttributes(G4Colour(0.8,0.6,0.));
//  M1VisAtt->SetForceWireframe(false);
  MVisAtt->SetForceSolid(true);
  MVisAtt->SetVisibility(true);
  // pixelVisAtt->SetForceAuxEdgeVisible(true);
  MLogical->SetVisAttributes(MVisAtt);

  G4VisAttributes* envelopeMVisAtt = new G4VisAttributes(G4Colour(1.0,0.,0.));
//  envelopeM1VisAtt->SetForceWireframe(false);
  envelopeMVisAtt->SetForceSolid(true);
  envelopeMVisAtt->SetVisibility(false);
  envelopeMVisAtt->SetForceAuxEdgeVisible(false);
  EnvelopeLogical->SetVisAttributes(envelopeMVisAtt);
  
  return std::make_tuple(EnvelopeLogical, EnvelopePhysical);  
  
}//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::tuple<G4LogicalVolume*,G4VPhysicalVolume*> PepiDetectorConstruction::CreateSubstrate(G4String name, G4double mag, G4ThreeVector position, G4double thickness, G4Material* material, G4LogicalVolume* SubLogical, G4VPhysicalVolume* SubPhysical)
{

  // - Build the substrate as a box
  G4Box* SubSolid = new G4Box(name,                               //its name
                          (1.2*fPixiRadSizeX/mag)/2,   //its size
                          (1.2*fPixiRadSizeY)/2,
                          thickness/2);
  G4String lvname = name + "LV";
     
  SubLogical = new G4LogicalVolume(SubSolid,       //its solid
                                       material,    //its material
                                       lvname);          //its name
 

  SubPhysical = new G4PVPlacement(0,              //its rotation
                                      position,       //its translation
                                      SubLogical,      //its logical volume
                                      name,              //its name
                                      fWorldLogical,       //its mother volume
                                      false,               //no boolean operation
                                      0,                   //copy number
                                      fCheckOverlaps);     //checking overlaps

  G4VisAttributes* SubVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
//  M1subVisAtt->SetForceWireframe(false);
  SubVisAtt->SetForceSolid(true);
  SubVisAtt->SetVisibility(true);
  SubVisAtt->SetForceAuxEdgeVisible(false);
  SubLogical->SetVisAttributes(SubVisAtt);
  
 return std::make_tuple(SubLogical, SubPhysical);
}  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiDetectorConstruction::Move(G4String name, G4LogicalVolume* Logical, G4VPhysicalVolume* Physical, G4ThreeVector position)
{
// for lateral translation of masks and sample

  Logical->RemoveDaughter(Physical);
  delete Physical;
  Physical = new G4PVPlacement(0,                                                  //its rotation
                                          position,   				      //its position
                                          Logical,                                    //its logical volume
                                          name,                                       //its name
                                          fWorldLogical,                                      //its mother volume
                                          false,                                              //no boolean operation
                                          0,                                                  //copy number
                                          fCheckOverlaps);                                    //checking overlaps                                
} 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}