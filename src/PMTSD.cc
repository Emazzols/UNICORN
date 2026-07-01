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
//
#include "PMTSD.hh"

#include "PMTHit.hh"
// #include "UserTrackInformation.hh"

#include "G4OpticalPhoton.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PMTSD::PMTSD(G4String name)
  : G4VSensitiveDetector(name)
{
  collectionName.insert("pmtHitCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PMTSD::Initialize(G4HCofThisEvent* HCE)
{
  fPMTHitCollection =
    new PMTHitsCollection(SensitiveDetectorName, collectionName[0]);

  if(fHCID < 0)
  {
    fHCID =
      G4SDManager::GetSDMpointer()->GetCollectionID(fPMTHitCollection);
  }
  HCE->AddHitsCollection(fHCID, fPMTHitCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool PMTSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  if(!aStep)
    return false;
  G4Track* theTrack = aStep->GetTrack();
  theTrack->SetTrackStatus(fStopAndKill);

  // Need to know if this is an optical photon
  if(theTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
  {
    return false;
  }

  // Find out information regarding the hit
  G4StepPoint* thePostPoint = aStep->GetPostStepPoint();

  G4LogicalVolume* volume
    = aStep->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
  // UserTrackInformation* trackInformation =
  //   (UserTrackInformation*) theTrack->GetUserInformation();

  G4TouchableHistory* theTouchable =
    (G4TouchableHistory*) (thePostPoint->GetTouchable());

  G4ThreeVector photonArrive = thePostPoint->GetPosition();
  //G4ThreeVector photonExit   = trackInformation->GetExitPosition();


  G4double arrivalTime       = theTrack->GetGlobalTime();
  G4double energy            = theTrack->GetTotalEnergy();

  G4double edepStep          = aStep->GetTotalEnergyDeposit();
  // G4cout << "total energy deposited in step: " << G4BestUnit(edepStep, "Energy") << G4endl;


  // Convert the global coordinate for arriving photons into
  // the local coordinate of the detector
  photonArrive =
    theTouchable->GetHistory()->GetTopTransform().TransformPoint(photonArrive);

  // Creating the hit and add it to the collection
  fPMTHitCollection->insert(
    new PMTHit(arrivalTime, energy, edepStep, photonArrive));//, photonExit));

  // if (edepStep==0.) return true;
  // auto hit = new PMTHit();
  //
  // hit->AddEdep(edepStep);

  return true;
}

G4bool PMTSD::ProcessHits_boundary(const G4Step* aStep, G4TouchableHistory*)
{
  // need to know if this is an optical photon
  if(aStep->GetTrack()->GetDefinition() !=
     G4OpticalPhoton::OpticalPhotonDefinition())
    return false;

  G4int pmtNumber =
    aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(1);
  G4VPhysicalVolume* physVol =
    aStep->GetPostStepPoint()->GetTouchable()->GetVolume(1);
  // G4cout << "phys volume is " << physVol << G4endl;

  size_t n       = fPMTHitCollection->entries();
  PMTHit* hit = nullptr;
  for(size_t i = 0; i < n; ++i)
  {
    if((*fPMTHitCollection)[i]->GetPMTNumber() == pmtNumber)
    {
      hit = (*fPMTHitCollection)[i];
      break;
    }
  }

  if(hit == nullptr)
  {                         // this pmt wasn't previously hit in this event
    hit = new PMTHit();  // so create new hit
    hit->SetPMTNumber(pmtNumber);
    hit->SetPMTPhysVol(physVol);
    fPMTHitCollection->insert(hit);
    // hit->SetPMTPos((*fPMTPositionsX)[pmtNumber], (*fPMTPositionsY)[pmtNumber],
    //                (*fPMTPositionsZ)[pmtNumber]);
  }

  hit->IncPhotonCount();  // increment hit for the selected pmt


  return true;
}

void PMTSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel > 0 ) {
     G4int nofHits = fPMTHitCollection->entries();
     G4cout << G4endl
            << "-------->Hits Collection: in this event there are " << nofHits
            << " hits in the photon detector: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fPMTHitCollection)[i]->Print();
  }
}
