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
// Original code by 
// - Zhentian Wang - 12.06.2009
// If you use this code plese cite
// Wang, Zhentian, et al. "Implement X-ray refraction effect in Geant4 for phase contrast imaging." IEEE, 2009.
//
// Modified by
// - Luca Brombal, INFN - Trieste - 20.06.2020
// - Camilo Sevilla, EAFIT - Medellin - 14.11.2023
// - Cristian Tibambre, UniAndes - Bogota D.C. - 14.11.2023
/////////////////////////////////////////////////////
//
/// \file PepiXrayRefraction.cc
/// \brief Implementation of the PepiXrayRefraction class
//


#include "G4ios.hh"
#include "G4Gamma.hh"
#include "PepiXrayRefraction.hh"
#include "G4GeometryTolerance.hh"
#include "G4SystemOfUnits.hh"

#include "G4PhysicalConstants.hh"
#include "G4OpProcessSubType.hh"

#include "G4OpBoundaryProcess.hh"

#include "G4VSensitiveDetector.hh"
#include "G4ParallelWorldProcess.hh"

namespace PEPI2
{

PepiXrayRefraction::PepiXrayRefraction(const G4String& processName,
		G4ProcessType type) :
	G4VDiscreteProcess(processName, type) {

  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
//  G4cout << "Tolerance " << kCarTolerance << G4endl;
}

PepiXrayRefraction::~PepiXrayRefraction() {
}

G4VParticleChange*PepiXrayRefraction::PostStepDoIt(const G4Track& aTrack,
		const G4Step& aStep) {

	G4cout.precision(3);

	aParticleChange.Initialize(aTrack);

	G4StepPoint* pPreStepPoint = aStep.GetPreStepPoint();
	G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();


 	if (pPostStepPoint->GetStepStatus() != fGeomBoundary) {
 		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
 	}

	if (aTrack.GetStepLength()<=kCarTolerance/2) {
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	Material1 = pPreStepPoint -> GetMaterial();
	Material2 = pPostStepPoint -> GetMaterial();

	const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

	thePhotonMomentum = aParticle->GetTotalMomentum();
	//G4cout<<"Photon energy "<< thePhotonMomentum<<G4endl;
	OldMomentum = aParticle->GetMomentumDirection();

	G4MaterialPropertiesTable* aMaterialPropertiesTable = 0;
	G4MaterialPropertyVector* Rindex = 0;

 	aMaterialPropertiesTable = Material1->GetMaterialPropertiesTable();
 	if (aMaterialPropertiesTable) {
		Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
	} else {
		// If no RINDEX is contained in the MaterialPropertiesTable Rindex = 1
	        Rindex1=1;
		//aParticleChange.ProposeTrackStatus(fStopAndKill);
		//return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	if (Rindex) {
	  G4bool bIsOutOfRange;
	  Rindex1 = Rindex->GetValue(thePhotonMomentum, bIsOutOfRange);
	} else {
		// if the photon energy is out of the range for which Rindex is specified Rindex = 1
	        Rindex1=1;
		//aParticleChange.ProposeTrackStatus(fStopAndKill);
		//return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	Rindex = NULL;
	if (Material1 == Material2) {
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}
	aMaterialPropertiesTable =Material2->GetMaterialPropertiesTable();
	if (aMaterialPropertiesTable)
		Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");

	if (Rindex) {
	  G4bool bIsOutOfRange;
	  Rindex2 = Rindex->GetValue(thePhotonMomentum, bIsOutOfRange);
	} else {
		// if the photon energy is out of the range for which Rindex is specified Rindex = 1
		Rindex2=1;
		//aParticleChange.ProposeTrackStatus(fStopAndKill);
		//return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	G4ThreeVector theGlobalPoint = pPostStepPoint->GetPosition();

	G4Navigator
			* theNavigator =G4TransportationManager::GetTransportationManager()->
			GetNavigatorForTracking();

	G4ThreeVector theLocalPoint = theNavigator->
	GetGlobalToLocalTransform().
	TransformPoint(theGlobalPoint);

	G4ThreeVector theLocalNormal;

	G4bool valid;
	theLocalNormal = theNavigator->GetLocalExitNormal(&valid);

	if (valid) {
		theLocalNormal = -theLocalNormal;
	} else {
		G4cerr << " Err in G4XrayRefraction/PostStepDoIt(). "<< G4endl;
	}
	
	theGlobalNormal = theNavigator->GetLocalToGlobalTransform().
	TransformAxis(theLocalNormal);
	if (OldMomentum * theGlobalNormal > 0.0) {

		theGlobalNormal = -theGlobalNormal;
	}
	
    Snell();//Snell's law
	NewMomentum = NewMomentum.unit();
	
	aParticleChange.ProposeMomentumDirection(NewMomentum);
	
	return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

void PepiXrayRefraction::Snell() {

	G4bool Swap = false;
	G4bool Through = false;


		if (Through)
		{
			Swap = !Swap;
			Through = false;
			theGlobalNormal = -theGlobalNormal;
			G4Swap(Material1, Material2);
			G4Swap(&Rindex1, &Rindex2);
		}
		theFacetNormal = theGlobalNormal;
		G4double PdotN = OldMomentum * theFacetNormal;
		cost1 = -PdotN;
		if (std::abs(cost1) < 1.0-kCarTolerance) {
			sint1 = std::sqrt(1.-cost1*cost1);
			sint2 = sint1*Rindex1/Rindex2; 

		} else {
			sint1 = 0.0;
			sint2 = 0.0;
		}
		// TOTAL REFLECTION
		if (sint2 >= 1.0) {
			if (Swap) Swap = !Swap;
			PdotN = OldMomentum * theFacetNormal;
			NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;

		// REFRACTION
		} else if (sint2 < 1.0) {
			if (cost1 > 0.0) {
				cost2 = std::sqrt(1.-sint2*sint2);
			} else {
				cost2 = -std::sqrt(1.-sint2*sint2);
			}
			Through = true;
			if (sint1 > 0.0) {
				G4double alpha = cost1 - cost2*(Rindex2/Rindex1);
				NewMomentum = OldMomentum + alpha*theFacetNormal;
				NewMomentum = NewMomentum.unit();
		 		/*G4cout << "OldMomentum " << OldMomentum << G4endl;
				 G4cout << "NewMomentum " << NewMomentum << G4endl;
				 printf ("Rindex1: %4.12f \n", Rindex1);
				 printf ("Rindex2: %4.12f \n", Rindex2);
				 printf ("delta1: %4.12f \n", 1-Rindex1);
				 printf ("delta2: %4.12f \n", 1-Rindex2);
				 printf ("sint1: %4.12f \n", sint1);
				 printf ("sint2: %4.12f \n", sint2);
				 G4cout << "Material1 " << Material1 << G4endl;
				 G4cout << "Material2 " << Material2 << G4endl;*/
				
//				PdotN = -cost2;
			} else { 
				NewMomentum = OldMomentum;

			}	
		}
		OldMomentum = NewMomentum.unit();

}

G4double PepiXrayRefraction::GetMeanFreePath(const G4Track&, G4double,
		G4ForceCondition* condition) {
	*condition = Forced;

	return DBL_MAX;
}

}
