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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"

#include "G4VProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "CreateTree.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction(const G4int& scint, const G4int& cher)
: G4UserStackingAction(),
fScintillationCounter(0), fCerenkovCounter(0),
propagateScintillation(scint),propagateCerenkov(cher)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  G4ParticleDefinition * particleType = aTrack->GetDefinition();
  if(particleType == G4OpticalPhoton::OpticalPhotonDefinition())
  { // particle is optical photon
    G4int parentID = aTrack->GetParentID();
    if(parentID != 0)
    { // particle is secondary
      G4String processName = aTrack->GetCreatorProcess()->GetProcessName();
      if(processName == "Scintillation")
      {
        fScintillationCounter++;
        if(!propagateScintillation)
        {
          return fKill;
        }
      }
      if(processName == "Cerenkov")
      {
        fCerenkovCounter++;
        if(!propagateCerenkov)
        {
          return fKill;
        }
      }
      G4LogicalVolume * thePostVolume  = aTrack->GetVolume()->GetLogicalVolume();
      G4String thePostPVname = thePostVolume->GetName();
      if((processName == "Cerenkov") && (thePostPVname == "photoDet_log"))
      {
        return fKill;
      }
    }
  }
  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StackingAction::NewStage()
{
  if(CreateTree::Instance()->verbose > 0) G4cout << "Number of Scintillation photons produced in this event : " << fScintillationCounter << G4endl;
  if(CreateTree::Instance()->verbose > 0) G4cout << "Number of Cerenkov photons produced in this event : "<< fCerenkovCounter << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StackingAction::PrepareNewEvent()
{
  fScintillationCounter = 0;
  fCerenkovCounter = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
