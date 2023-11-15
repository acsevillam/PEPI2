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
/// \file PepiRun.cc
/// \brief Implementation of the PepiRun class
//

//=====================================================================
//
//  (Description)
//    PepiRun Class is for accumulating scored quantities which is 
//  scored using G4MutiFunctionalDetector and G4VPrimitiveScorer.
//  Accumulation is done using G4THitsMap object.
//
//    The constructor PepiRun(const std::vector<G4String> mfdName)
//  needs a vector filled with MultiFunctionalDetector names which
//  was assigned at instantiation of MultiFunctionalDetector(MFD).
//  Then PepiRun constructor automatically scans primitive scorers
//  in the MFD, and obtains collectionIDs of all collections associated
//  to those primitive scorers. Futhermore, the G4THitsMap objects 
//  for accumulating during a RUN are automatically created too.
//  (*) Collection Name is same as primitive scorer name.
// 
//    The resultant information is kept inside PepiRun objects as
//  data members.
//  std::vector<G4String> fCollName;            // Collection Name,
//  std::vector<G4int> fCollID;                 // Collection ID,
//  std::vector<G4THitsMap<G4double>*> fRunMap; // HitsMap for RUN.
//
//  The resualtant HitsMap objects are obtain using access method,
//  GetHitsMap(..).
//
//=====================================================================

#include "PepiRun.hh"
#include "G4SDManager.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//  Constructor. 
//   (The vector of MultiFunctionalDetector name has to be given.)

namespace PEPI2
{

PepiRun::PepiRun(const std::vector<G4String> mfdName) : G4Run()
{
  G4SDManager* pSDman = G4SDManager::GetSDMpointer();
  //G4SDManager* pSDman = G4SDManager::GetSDMpointerIfExist();

  //=================================================
  //  Initalize RunMaps for accumulation.
  //  Get CollectionIDs for HitCollections.
  //=================================================
  G4int nMfd = mfdName.size();

  // Loop for all MultiFunctionalDetectors
  for (G4int idet = 0; idet < nMfd ; idet++)
  {  
    
    G4String detName = mfdName[idet];
    
    // Seek and Obtain MFD objects from SDmanager.
    G4MultiFunctionalDetector* mfd = (G4MultiFunctionalDetector*)(pSDman->FindSensitiveDetector(detName));
    if (mfd)
    {
      // Loop over the registered primitive scorers.
      for (G4int icol = 0; icol < mfd->GetNumberOfPrimitives(); icol++)
      {
        // Get Primitive Scorer object.
        G4VPrimitiveScorer* scorer=mfd->GetPrimitive(icol);

        // collection name and collectionID for HitsCollection,
        // where type of HitsCollection is G4THitsMap in case of primitive 
        // scorer.
        // The collection name is given by <MFD name>/<Primitive Scorer name>.
        G4String collectionName     = scorer->GetName();
        G4String fullCollectionName = detName+"/"+collectionName;
        G4int    collectionID       = pSDman->GetCollectionID(fullCollectionName);
        
        if (collectionID >= 0)
        {
           G4cout << "++ "<<fullCollectionName<< " id " << collectionID << G4endl;

          // Store obtained HitsCollection information into data members.
          // And, creates new G4THitsMap for accumulating quantities during RUN.
          fCollName.push_back(fullCollectionName);
          fCollID.push_back(collectionID);
          fRunMap.push_back(new G4THitsMap<G4double>(detName,collectionName));
        }
        else
        {
           G4cout << "** collection " << fullCollectionName << " not found. " << G4endl;
        }
      }
    }
  }

}

// Destructor
// Clear all data members.
PepiRun::~PepiRun()
{
  // Clear HitsMap for RUN
  G4int nMap = fRunMap.size();

  for (G4int i = 0; i < nMap; i++)
  {
    if(fRunMap[i]) fRunMap[i]->clear();
  }
  fCollName.clear();
  fCollID.clear();
  fRunMap.clear();
}

//  RecordEvent is called at end of event.
//  For scoring purpose, the resultant quantity in a event,
//  is accumulated during a Run.
void PepiRun::RecordEvent(const G4Event* aEvent)
{
  numberOfEvent++;  // This is an original line.

  //=============================
  // HitsCollection of This Event
  //============================
  G4HCofThisEvent* pHCE = aEvent->GetHCofThisEvent();

  if (!pHCE) return;

  //=======================================================
  // Sum up HitsMap of this Event  into HitsMap of this RUN
  //=======================================================
  G4int nCol = fCollID.size();
  //G4cout<<"\n\n\n\n\n" <<nCol <<G4endl;
  for ( G4int i = 0; i < nCol ; i++ ){  // Loop over HitsCollection
    G4THitsMap<G4double>* evtMap=0;
    if (fCollID[i] >= 0)
    {           // Collection is attached to pHCE
      evtMap = (G4THitsMap<G4double>*)(pHCE->GetHC(fCollID[i]));
    } else
    {
      G4cout <<" Error evtMap Not Found "<< i << G4endl;
    }
    if (evtMap)
    {
      //=== Sum up HitsMap of this event to HitsMap of RUN.===
      *fRunMap[i] += *evtMap;
      //======================================================
    }
  }
  G4Run::RecordEvent(aEvent);

}

void PepiRun::Merge(const G4Run * aRun)
{
  const PepiRun* localRun = static_cast<const PepiRun *>(aRun);
  
  //=======================================================
  // Merge HitsMap of working threads
  //=======================================================
  G4int nCol = localRun->fCollID.size();

  for (G4int i = 0; i < nCol ; i++)
  {  
    // Loop over HitsCollection
    if (localRun->fCollID[i] >= 0)
    {
      *fRunMap[i] += *localRun->fRunMap[i];
    }
  }

  G4Run::Merge(aRun);
}

// Access method for HitsMap of the RUN
// By MultiFunctionalDetector name and Collection Name.
G4THitsMap<G4double>* PepiRun::GetHitsMap(const G4String& detName, const G4String& colName)
{
  G4String fullName = detName+"/"+colName;
      //G4cout<<"\n\n\n\n\n\n\n\n\n\n\n" <<detName<< colName << fullName << G4endl;
  return GetHitsMap(fullName);
}

// Access HitsMap.
// By full description of collection name, that is
// <MultiFunctional Detector Name>/<Primitive Scorer Name>
G4THitsMap<G4double>* PepiRun::GetHitsMap(const G4String& fullName)
{
  G4int nCol = fCollName.size();

  for (G4int i = 0; i < nCol; i++)
  {
    if (fCollName[i] == fullName)
    {
        return fRunMap[i];
    }
  }
  return NULL;
}

// Dump All HitsMap of this RUN. (for debuging and monitoring of quantity).
// This method calls G4THisMap::PrintAll() for individual HitsMap.
void PepiRun::DumpAllScorers()
{
  // - Number of HitsMap in this RUN.
  G4int n = GetNumberOfHitsMap();
  // - GetHitsMap and dump values.
  for ( G4int i = 0; i < n ; i++ ){
    G4THitsMap<G4double>* runMap =GetHitsMap(i);
    if ( runMap ) {
      G4cout << " PrimitiveScorer RUN " 
             << runMap->GetSDname() <<","<< runMap->GetName() << G4endl;
      G4cout << " Number of entries " << runMap->entries() << G4endl;
      std::map<G4int,G4double*>::iterator itr = runMap->GetMap()->begin();
      for(; itr != runMap->GetMap()->end(); itr++) {
        G4cout << "  copy no.: " << itr->first
               << "  Run Value : " << *(itr->second) 
               << G4endl;
      }
    }
  }
}

}

