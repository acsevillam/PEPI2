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
/// \file PepiPhysicsXrayRefraction.cc
/// \brief Implementation of the PepiPhysicsXrayRefraction class
//

#include "PepiPhysicsXrayRefraction.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "PepiXrayRefraction.hh"
#include "G4Gamma.hh"

#include "G4PhysicsConstructorFactory.hh"

// #include "G4PhysicsListHelper.hh"
// #include "G4BuilderType.hh"

namespace PEPI2
{

G4_DECLARE_PHYSCONSTR_FACTORY(PepiPhysicsXrayRefraction);
// ---------------------------------------------------------------------------


PepiPhysicsXrayRefraction::PepiPhysicsXrayRefraction()
  : G4VPhysicsConstructor("G4XrayRefraction")
{ }

// ---------------------------------------------------------------------------
PepiPhysicsXrayRefraction::~PepiPhysicsXrayRefraction()
{ }

// ---------------------------------------------------------------------------
void PepiPhysicsXrayRefraction::ConstructParticle()
{
 G4Gamma::Gamma();
}

void PepiPhysicsXrayRefraction::ConstructProcess()
{
  // - Add XrayRefraction for gammas
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleTable::G4PTblDicIterator* particleIterator = particleTable->GetIterator();

  particleIterator->reset();

  while ((*particleIterator)())
  {
    G4ParticleDefinition* particle = particleIterator->value();

    if (particle == G4Gamma::Gamma()) 
    {
      G4ProcessManager* pManager = particle->GetProcessManager();
      pManager -> AddDiscreteProcess(new PepiXrayRefraction);
    }   
  }
}

}