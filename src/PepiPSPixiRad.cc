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
/// \file PepiPSPixiRad.cc
/// \brief Implementation of the PepiPSPixiRad class
//

#include "PepiPSPixiRad.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4SteppingManager.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4RunManager.hh"

#include <stdio.h>
#include <ios>
#include <iomanip>
#include <string>

namespace PEPI2
{

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This is a primitive scorer class for scoring energy deposit.
// 
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 2010-07-22   Introduce Unit specification.
// 2021-06-09   Modified by Luca Brombal and Nicola Poles, University of Trieste
// 
///////////////////////////////////////////////////////////////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PepiPSPixiRad::PepiPSPixiRad(G4String name, G4double threshold, 
                               G4int nPixelsX, G4int nPixelsY,
                               G4double pixelSizeX, G4double pixelSizeY, G4double pixelSizeZ,
                               const G4String& unit, G4int depth)
: G4VPrimitiveScorer(name,depth),
  HCID(-1),
  EvtMap(0),
  fThreshold(threshold),
  fnPixelsX(nPixelsX),
  fnPixelsY(nPixelsY),
  fPixelSizeX(pixelSizeX),
  fPixelSizeY(pixelSizeY),
  fPixelSizeZ(pixelSizeZ),
  fHitpos(0),
  fIndex(0),
  feDep(0)
{
  SetUnit(unit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PepiPSPixiRad::~PepiPSPixiRad()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool PepiPSPixiRad::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double eDep = aStep->GetTotalEnergyDeposit();

  if (eDep == 0) return FALSE;

  eDep *= aStep->GetPreStepPoint()->GetWeight(); // (Particle Weight)

  G4int index = GetIndex(aStep);

  G4Track* aTrack = aStep->GetTrack();

  // G4int PDGencoding = aTrack->GetDefinition()->GetPDGEncoding();
   G4int parentID    = aTrack->GetParentID();
  //G4int trackID     = aTrack->GetTrackID();
  // G4double kinetic = aStep->GetPreStepPoint()->GetKineticEnergy();
   G4ThreeVector stpos = aStep->GetPreStepPoint()->GetPosition();

  // - Get the track's position in global coordinates
  G4ThreeVector vtpos = aTrack->GetPosition();

  // - Convert the track position to local coordinates
  G4TouchableHandle theTouchable = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4ThreeVector ltpos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(vtpos);

    EvtMap->add(index, eDep);
    
    // store the position of the hit where most of the energy is released
    if (parentID==0) {feDep=0;}
    if (feDep<eDep)
    {
    feDep = eDep;
    fHitpos = G4ThreeVector(ltpos.x(), ltpos.y(), 0);//hit position
    fIndex = index;
    }

    //store the local position of the first hit
    /*if (parentID==0)
    {
    fHitpos = G4ThreeVector(ltpos.x(), ltpos.y(), 0);//first hit position
    fIndex = index;
        //G4cout << "\n\n first hit position" << fHitpos << "\n\n" << G4endl;
        //G4cout << "\n\n first hit pixel index " << fIndex << "\n\n" << G4endl;
    }*/
    
    
  // - Debugging Output
   /* G4cout << "paricle " << PDGencoding << G4endl
    << "parentID " << parentID << G4endl
    << "pre step pos " << stpos/um << G4endl
    << "kinetic energy " << kinetic/keV << " keV" << G4endl
    << "trackID " << trackID << G4endl
    << "process " << aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl
    << "post step pos in pixel" << ltpos/um << G4endl
    << "deposited energy " << eDep/keV << "keV" << G4endl
    << "3D displacement " <<  vtpos/um - stpos/um << "um" << G4endl
    << "lateral displacement " <<  std::sqrt((vtpos.x()-stpos.x())*(vtpos.x()-stpos.x())+(vtpos.y()-stpos.y())*(vtpos.y()-stpos.y()))/um << "um" << G4endl
    << "pixel index " << index << G4endl;*/
    //return TRUE;
  //}
/*


*/
/*
  G4cout  << "\tSTEP BEGIN\n"
    << "eDep " << eDep/keV << " keV" << G4endl
    << "nCouples " << nCouples << G4endl
    << "paricle " << PDGencoding << G4endl
    << "process " << aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl
    << "length " << G4BestUnit(aStep->GetStepLength(), "Length") << G4endl
    << "vtx " << vtpos.x()/um << " um" << G4endl
    << "vty " << vtpos.y()/um << " um" << G4endl
    << "vtz " << vtpos.z()/mm << " mm" << G4endl

    << "stx " << stpos.x()/um << " um" << G4endl
    << "sty " << stpos.y()/um << " um" << G4endl
    << "stz " << stpos.z()/mm << " mm" << G4endl

    << "ltx " << ltpos.x()/um << " um" << G4endl
    << "lty " << ltpos.y()/um << " um" << G4endl
    << "ltz " << ltpos.z()/um << " um" << G4endl
    << "depth " << depth/um << " um" << G4endl
    << "sigma " << sigma/um << " um" << G4endl
*/

  return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int PepiPSPixiRad::GetIdx(G4double pos, G4double size)
{
  size /= 2;

  if (G4int(pos/size) == 1 || G4int(pos/size) == 2)
  {
    return 1;
  }
  else if(G4int(pos/size) > 2)
  {
    return 2; 
  }
  else if (G4int(pos/size) == -1 || G4int(pos/size) == -2)
  {
    return -1;
  }
  else if (G4int(pos/size) < -2)
  {
    return -2;
  }
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiPSPixiRad::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(),
				    "low_th");			    
  if(HCID < 0) {HCID = GetCollectionID(0);}
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiPSPixiRad::EndOfEvent(G4HCofThisEvent*)
{ 
  // - Debugging Output
  /*G4cout << "\tlateral position biggest hit: " << fHitpos/um << " um" << G4endl;      
  G4cout << "\tpixel index biggest hit: " << fIndex << G4endl;      
  */
  G4int idx2, idx3, idx4;
  // computes the cluster of 4 pixels involved in the charge summing  
  if(fHitpos.x() < 0 && fHitpos.y() < 0)
  {    
      idx2=fIndex -1;
      idx3=fIndex -1 + fnPixelsX;
      idx4=fIndex + fnPixelsX;
 }
  else if(fHitpos.x() < 0 && fHitpos.y() > 0)
  {
      idx2=fIndex -1;
      idx3=fIndex -1 - fnPixelsX;
      idx4=fIndex - fnPixelsX;
  }
  else if(fHitpos.x() > 0 && fHitpos.y() < 0)
  {
      idx2=fIndex +1;
      idx3=fIndex +1 + fnPixelsX;
      idx4=fIndex + fnPixelsX;
  }
  else
  {
      idx2=fIndex +1;
      idx3=fIndex +1 - fnPixelsX;
      idx4=fIndex - fnPixelsX;
  }  
           
  std::map<G4int,G4double*>::iterator itr;
  G4double edep2=0;

  
  // sums the energy of the cluster
  for (itr = EvtMap->GetMap()->begin(); itr != EvtMap->GetMap()->end(); itr++)
  {  
    if ((itr->first) == fIndex || (itr->first) == idx2 || (itr->first) == idx3 || (itr->first) == idx4)
    {
       edep2 += *(itr->second)/keV;
    }    
  }
  
 // G4cout << "\tPixel: " << index2 << G4endl;
 // G4cout << "\taccumulated " << edep2 << G4endl;

  
  for (itr = EvtMap->GetMap()->begin(); itr != EvtMap->GetMap()->end(); itr++)
  {
    G4double newval = 0;
    // energy blurring calibration from Di Trapani et al. 2020
    G4double radThreshold = CLHEP::RandGauss::shoot(fThreshold, (2.635+0.02918*(*(itr->second)/keV))/2.35*keV); 
    G4int i=0;
    // if the pixel index belong to the cluster it assigns all the signal to the first pixel to be hit (fIndex)
    if ((itr->first) == fIndex || (itr->first) == idx2 || (itr->first) == idx3 || (itr->first) == idx4)
    {
       if(edep2 >= radThreshold/keV && i<1)
       {
          newval = 1;
          //G4cout << "1st if \tpixel no.: " << fIndex << " absorbed energy: " << radThreshold << "keV" << G4endl;   
          EvtMap->set(fIndex,newval);
          if(fIndex != itr->first)
          { 
            newval = 0.0;
            EvtMap->set(itr->first,newval);
          }
          //G4cout << "\tpixel no.: " << fIndex << " Events counted: " << *(itr->second) << G4endl;
          i++;
       }
       else
       {  
          newval = 0.0;
          EvtMap->set(itr->first,newval);
       }
    }
    else
    {    
        if (*(itr->second) >= radThreshold)
        {
          newval = 1;
  //        G4cout << "2nd if \tpixel no.: " << itr->first << " absorbed energy: " << *(itr->second)/keV << "keV" << G4endl;
  //    txtFile << itr->first << "\t" << *(itr->second)/keV << "\n";
          EvtMap->set(itr->first,newval);

          //G4cout << "\tpixel no.: " << itr->first << " Events counted: " << *(itr->second) << G4endl;
        }
        else
        {   
            newval = 0.0;
	    EvtMap->set(itr->first,newval);
        }
    }
    //G4cout << "out \tpixel no.: " << itr->first << " Events counted: " << *(itr->second) << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiPSPixiRad::clear()
{
  EvtMap->clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiPSPixiRad::DrawAll()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiPSPixiRad::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  energy deposit: " 
	   << *(itr->second)/GetUnitValue()
	   << " [" << GetUnit()<<"]"
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiPSPixiRad::SetUnit(const G4String& unit)
{
	CheckAndSetUnit(unit,"Energy");
}

}