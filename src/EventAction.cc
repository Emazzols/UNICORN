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
// $Id: EventAction.cc 75604 2013-11-04 13:17:26Z gcosmo $
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"

#include "RunAction.hh"
#include "PMTHit.hh"
#include "ScintHit.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4SDManager.hh"
#include "G4AnalysisManager.hh"

#include "Randomize.hh"
#include <iomanip>
#include "CreateTree.hh"

#include <fstream>
#include <iostream>

// const G4int number_of_detectors;
using std::array;

namespace {

// Utility function which finds a hit collection with the given Id
// and print warnings if not found
G4VHitsCollection* GetHC(const G4Event* event, G4int collId) {
  auto hce = event->GetHCofThisEvent();
  if (!hce) {
      G4ExceptionDescription msg;
      msg << "No hits collection of this event found." << G4endl;
      G4Exception("EventAction::EndOfEventAction()",
                  "Code001", JustWarning, msg);
      return nullptr;
  }

  auto hc = hce->GetHC(collId);
  if ( ! hc) {
    G4ExceptionDescription msg;
    msg << "Hits collection " << collId << " of this event not found." << G4endl;
    G4Exception("EventAction::EndOfEventAction()",
                "Code001", JustWarning, msg);
  }
  return hc;
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction,int logEscape)
: G4UserEventAction(),
fRunAction(runAction),
fLogEscape(logEscape)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* event)
{

  //print event number
  G4int eventID = event->GetEventID();
  // if(CreateTree::Instance()->verbose > 0)
  // {
  //   G4cout << "//***********************************************//" << G4endl;
  //   G4cout << "//              Event Number " << eventID		<< G4endl;
  //   G4cout << "//***********************************************//" << G4endl;
  // }
  //
  // if(CreateTree::Instance()->verbose == 0)
  // {
  //   for (G4int i = 1; i < 100; i++)
  //   if(eventID == i*1000){
  //     G4cout << "//***********************************************//" << G4endl;
  //     G4cout << "//              Reached processing of " << i*1000	<<  " events!" << G4endl;
  //     G4cout << "//***********************************************//" << G4endl;
  //   }
  // }
  auto sdManager = G4SDManager::GetSDMpointer();

  totalNumOfPhotons = 0;
  totalEnergyDeposited = 0;
  fHitCount = 0;

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if(fScintCollID < 0)
    fScintCollID = SDman->GetCollectionID("scintHitCollection");
  if(fPMTCollID < 0)
    fPMTCollID = SDman->GetCollectionID("pmtHitCollection");


  // switch (number_of_detectors)
  // {
  //   case 1: {
  //     array<G4String, oneDet> hHCNameSingle;
  //     hHCNameSingle = {"photonDet0/PhotonDetHitCollection"};
  //     fPhotSingleHCID[0]   = sdManager->GetCollectionID(hHCNameSingle[0]);
  //     break;
  //   }
  //   case 2: {
  //     array<G4String, twoDet> hHCNameDouble;
  //     hHCNameDouble = {{"photonDet0/PhotonDetHitCollection", "photonDet1/PhotonDetHitCollection"}};
  //     fPhotDoubleHCID[0]   = sdManager->GetCollectionID(hHCNameDouble[0]);
  //     fPhotDoubleHCID[1]   = sdManager->GetCollectionID(hHCNameDouble[1]);
  //     break;
  //   }
  //   case 6: {
  //     array<G4String, sixDet> hHCName;
  //     hHCName = {{"photonDet0/PhotonDetHitCollection", "photonDet1/PhotonDetHitCollection",
  //                 "photonDet2/PhotonDetHitCollection", "photonDet3/PhotonDetHitCollection",
  //                 "photonDet4/PhotonDetHitCollection", "photonDet5/PhotonDetHitCollection"}};
  //     for (G4int iDet = 0; iDet < sixDet; ++iDet) {
  //       // hit collections IDs
  //       fPhotHCID[iDet]   = sdManager->GetCollectionID(hHCName[iDet]);
  //     }
  //     break;
  //   }
  //   default:
  //     G4cout << "Select number of detectors as 1, 2 or 6! " << G4endl;
  // }
  //
  // for (G4int iDet = 0; iDet < number_of_detectors; ++iDet) {
  //   // hit collections IDs
  //   fPhotHCID[iDet]   = sdManager->GetCollectionID(hHCName[iDet]);
  // }



  Int_t run = CreateTree::Instance()->Run;
  long int seed = CreateTree::Instance()->Seed;
  CreateTree::Instance()->Clear();


  CreateTree::Instance()->Seed = seed;
  CreateTree::Instance()->Run = run;
  CreateTree::Instance()->Event = event->GetEventID();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  // Accumulate statistics
  auto analysisManager = G4AnalysisManager::Instance();

  CreateTree::Instance()->FoundGammaVertexMomentum = false;
  CreateTree::Instance()->primaryFirstEntranceFound = false;
  CreateTree::Instance()->SourceX = (Float_t) event->GetPrimaryVertex()->GetX0();
  CreateTree::Instance()->SourceY = (Float_t) event->GetPrimaryVertex()->GetY0();
  CreateTree::Instance()->SourceZ = (Float_t) event->GetPrimaryVertex()->GetZ0();

  if(CreateTree::Instance()->totalEnergyDepositedInCrystals > 0)
  {
    CreateTree::Instance()->primaryInteracted = 1;
  }

  // G4cout << "energy deposited in crystal = " << totalEnergyDeposited/CLHEP::MeV << " MeV." << G4endl;

  if(fLogEscape) // always save
  {
    CreateTree::Instance()->Fill();
    if(CreateTree::Instance()->verbose > 0) G4cout << "Total Energy deposited in this event = "<< CreateTree::Instance()->totalEnergyDeposited << " MeV" << G4endl;
  }
  else
  {
    if(CreateTree::Instance()->totalEnergyDeposited > 0)
    {
      CreateTree::Instance()->Fill();
      if(CreateTree::Instance()->verbose > 0) G4cout << "Total Energy deposited in this event = "<< CreateTree::Instance()->totalEnergyDeposited << " MeV" << G4endl;
    }
    else
    {
      if(CreateTree::Instance()->verbose > 0) G4cout << "No Energy deposited in this event" << G4endl;
    }
  }

  ScintHitsCollection* scintHC = nullptr;
  PMTHitsCollection* pmtHC     = nullptr;
  G4HCofThisEvent* hitsCE      = event->GetHCofThisEvent();

  if(hitsCE)
  {
    if(fScintCollID >= 0)
    {
      scintHC = (ScintHitsCollection*) (hitsCE->GetHC(fScintCollID));
    }
    if(fPMTCollID >= 0)
    {
      pmtHC = (PMTHitsCollection*) (hitsCE->GetHC(fPMTCollID));
    }
  }
  // Hits in scintillator
  if(scintHC)
  {
    size_t n_hit = scintHC->entries();
    G4ThreeVector eWeightPos(0.);
    G4double edep;
    G4double edepMax = 0;

    for(size_t i = 0; i < n_hit; ++i)
    {  // gather info on hits in scintillator
      edep = (*scintHC)[i]->GetEdep();
      totalEnergyDeposited += edep;
    }
    // G4cout << "\tTotal energy deposition in scintillator : " << totalEnergyDeposited / MeV << " (MeV)" << G4endl;
  }

  if(pmtHC)
  {
    size_t pmts = pmtHC->entries();
    for(size_t i = 0; i < pmts; ++i)
    {
      fHitCount += (*pmtHC)[i]->GetPhotonCount();
    }
    // G4cout << "\tTotal number of hits : " << fHitCount << G4endl;
    CreateTree::Instance()->NumOptPhotonsReachedPMT = fHitCount;
    CreateTree::Instance()->TotalNumOptPhotonsReachedPMT = fHitCount;
    CreateTree::Instance()->photons->Fill();
    CreateTree::Instance()->photonsTotal->Fill();

  }

  // analysisManager->FillH1(0, n_hit);
  // CreateTree::Instance()->NumOptPhotonsReachedPMT = n_hit;
  // CreateTree::Instance()->TotalNumOptPhotonsReachedPMT = n_hit;
  // CreateTree::Instance()->photons->Fill();
  // CreateTree::Instance()->photonsTotal->Fill();
  // G4cout << "Photodetector has " << n_hit  << " hits." << G4endl;
  // switch (number_of_detectors)
  // {
  //   case 1: {
  //     auto hc = GetHC(event, fPhotSingleHCID[0]);
  //     if ( ! hc ) return;
  //     auto n_hit = hc->GetSize();
  //
  //     // for(int i = 0; i < n_hit; i++)
  //     // {
  //     //   auto hit = static_cast<PhotonDetHit*>(hc->GetHit(i));
  //     //   hit->Print();
  //     // }
  //
  //     // analysisManager->FillH1(0, n_hit);
  //     CreateTree::Instance()->NumOptPhotonsReachedPMT = n_hit;
  //     CreateTree::Instance()->TotalNumOptPhotonsReachedPMT = n_hit;
  //     CreateTree::Instance()->photons->Fill();
  //     CreateTree::Instance()->photonsTotal->Fill();
  //     // G4cout << "Photodetector has " << n_hit  << " hits." << G4endl;
  //     break;
  //   }
  //   case 2: {
  //     for (G4int iDet = 0; iDet < twoDet; ++iDet)
  //     {
  //       auto hc = GetHC(event, fPhotDoubleHCID[iDet]);
  //       if ( ! hc ) return;
  //       auto n_hit = hc->GetSize();
  //
  //       // for(int i = 0; i < n_hit; i++)
  //       // {
  //       //   auto hit = static_cast<PhotonDetHit*>(hc->GetHit(i));
  //       //   hit->Print();
  //       // }
  //
  //       // analysisManager->FillH1(0, n_hit);
  //       CreateTree::Instance()->NumOptPhotonsReachedPMT = n_hit;
  //       CreateTree::Instance()->photons->Fill();
  //       totalNumOfPhotons += n_hit;
  //       // G4cout << "Photodetector " << iDet << " has " << n_hit  << " hits." << G4endl;
  //     }
  //     CreateTree::Instance()->TotalNumOptPhotonsReachedPMT = totalNumOfPhotons;
  //     CreateTree::Instance()->photonsTotal->Fill();
  //     break;
  //   }
  //   case 6: {
  //     for (G4int iDet = 0; iDet < sixDet; ++iDet)
  //     {
  //       auto hc = GetHC(event, fPhotHCID[iDet]);
  //       if ( ! hc ) return;
  //       auto n_hit = hc->GetSize();
  //
  //       // for(int i = 0; i < n_hit; i++)
  //       // {
  //       //   auto hit = static_cast<PhotonDetHit*>(hc->GetHit(i));
  //       //   hit->Print();
  //       // }
  //
  //       // analysisManager->FillH1(0, n_hit);
  //       CreateTree::Instance()->NumOptPhotonsReachedPMT = n_hit;
  //       CreateTree::Instance()->photons->Fill();
  //       totalNumOfPhotons += n_hit;
  //       // G4cout << "Photodetector " << iDet << " has " << n_hit  << " hits." << G4endl;
  //     }
  //     CreateTree::Instance()->TotalNumOptPhotonsReachedPMT = totalNumOfPhotons;
  //     CreateTree::Instance()->photonsTotal->Fill();
  //     break;
  //   }
  //   default:
  //     G4cout << "Select number of detectors as 1, 2 or 6! " << G4endl;
  // }
  // for (G4int iDet = 0; iDet < number_of_detectors; ++iDet) {
  //   auto hc = GetHC(event, fPhotHCID[iDet]);
  //   if ( ! hc ) return;
  //   auto n_hit = hc->GetSize();
  //   G4cout << "Photodetector " << iDet << " has " << n_hit  << " hits." << G4endl;

    // auto n_hit = hc->GetSize();
    // analysisManager->FillH1(2, n_hit);
    // for (unsigned int i = 0; i < n_hit; ++i) {
    //   auto hit = static_cast<PhotonDetHit*>(hc->GetHit(i));
    //   analysisManager->FillH1(0, hit->GetEnergy());
    //   analysisManager->FillH1(1, hit->GetArrivalTime());
    // }
    // }


}
