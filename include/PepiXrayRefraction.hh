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
/// \file PepiXrayRefraction.hh
/// \brief Definition of the PepiXrayRefraction class

#ifndef PepiXrayRefraction_h
#define PepiXrayRefraction_h 1

#include "globals.hh"
#include "templates.hh"
#include "geomdefs.hh"
#include "Randomize.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4TransportationManager.hh"
#include "G4Gamma.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4OpticalPhoton.hh"

/// PEPI x-ray refraction class

namespace PEPI2
{

class PepiXrayRefraction : public G4VDiscreteProcess 
{

private:

public: 

    PepiXrayRefraction(const G4String& processName = "XrayRefraction",
                                     G4ProcessType type = fOptical);

	~PepiXrayRefraction();

public: 

    G4bool IsApplicable(const G4ParticleDefinition& aParticleType);

	G4double GetMeanFreePath(const G4Track& ,
				 G4double ,
				 G4ForceCondition* condition);

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
				       const G4Step&  aStep);

private:

	void Snell();

	void G4Swap(G4double* a, G4double* b) const;

	void G4Swap(G4Material* a, G4Material* b) const;

	void G4VectorSwap(G4ThreeVector* vec1, G4ThreeVector* vec2) const;

private:

	G4Material* Material1;
	G4Material* Material2;

	G4double Rindex1;
	G4double Rindex2;

	G4double cost1, cost2, sint1, sint2;
	
	G4double thePhotonMomentum;
	G4ThreeVector OldMomentum;
	G4ThreeVector NewMomentum;

	G4ThreeVector theGlobalNormal;
	G4ThreeVector theFacetNormal;

        G4double kCarTolerance;
};

inline
G4bool PepiXrayRefraction::IsApplicable(const G4ParticleDefinition& 
					               aParticleType)
{
   return  ( &aParticleType == G4Gamma::Gamma() );
}

inline
void PepiXrayRefraction::G4Swap(G4double* a, G4double* b) const
{


  G4double temp;

  temp = *a;
  *a = *b;
  *b = temp;
}

inline
void PepiXrayRefraction::G4Swap(G4Material* a, G4Material* b) const
{

   G4Material* temp = a;

   a = b;
   b = temp;
}

inline
void PepiXrayRefraction::G4VectorSwap(G4ThreeVector* vec1,
				       G4ThreeVector* vec2) const
{

  G4ThreeVector temp;

  temp = *vec1;
  *vec1 = *vec2;
  *vec2 = temp;
}

}

#endif // PepiXrayRefraction_h
