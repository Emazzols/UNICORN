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
/// \file optical/wls/include/WLSPMTSD.hh
/// \brief Definition of the WLSPMTSD class
//
//
#include "PMTHit.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


G4ThreadLocal G4Allocator<PMTHit>* PMTHitAllocator = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PMTHit::PMTHit()
{
  fArrivalTime = 0.;
  fEnergy = 0.;
  fEdep = 0.;
  fPosArrive   = G4ThreeVector(0., 0., 0.);
  //fPosExit     = G4ThreeVector(0., 0., 0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PMTHit::PMTHit(G4double pTime, G4double pEnergy, G4double pEdep, G4ThreeVector pArrive)//, G4ThreeVector pExit)
{
  //fPosExit     = pExit;
  fPosArrive   = pArrive;
  fArrivalTime = pTime;
  fEnergy      = pEnergy;
  fEdep        = pEdep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PMTHit::~PMTHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PMTHit::PMTHit(const PMTHit& right)
  : G4VHit()
{
  fPmtNumber = right.fPmtNumber;
  fPhysVol   = right.fPhysVol;
  fPhotons   = right.fPhotons;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const PMTHit& PMTHit::operator=(const PMTHit& right)
{
  //fPosExit     = right.fPosExit;
  fPmtNumber = right.fPmtNumber;
  fPhysVol   = right.fPhysVol;
  fPhotons   = right.fPhotons;

  fPosArrive   = right.fPosArrive;
  fArrivalTime = right.fArrivalTime;
  fEnergy      = right.fEnergy;
  fEdep        = right.fEdep;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool PMTHit::operator==(const PMTHit& right) const
{
  return fArrivalTime == right.fArrivalTime && fEnergy == right.fEnergy && fEdep == right.fEdep && fPmtNumber == right.fPmtNumber;
}

void PMTHit::Print()
{
  G4cout
    << "Arrival time: "
    << std::setw(7) << G4BestUnit(fArrivalTime, "Time")
    << "Arrival position: ("
    << std::setw(7) << G4BestUnit(fPosArrive.x(), "Length") << ", "
    << std::setw(7) << G4BestUnit(fPosArrive.y(), "Length") << "); "
    // << "Exit position: ("
    // << std::setw(7) << G4BestUnit(fPosExit.x(), "Length") << ", "
    // << std::setw(7) << G4BestUnit(fPosExit.y(), "Length") << "); "
    << "Particle energy: "
    << std::setw(7) << G4BestUnit(fEnergy, "Energy")
    << "Energy deposited: "
    << std::setw(7) << G4BestUnit(fEdep, "Energy")
    // << fEdep/keV << " (keV) "
    << G4endl;
}
