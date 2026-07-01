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
/// \file optical/wls/include/WLSPMTHit.hh
/// \brief Definition of the WLSPMTHit class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef PMTHit_h
#define PMTHit_h 1

#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "G4VPhysicalVolume.hh"

class PMTHit : public G4VHit
{
 public:
  PMTHit();
  PMTHit(G4double pTime, G4double pEnergy, G4double pEdep, G4ThreeVector pArrive);//, G4ThreeVector pExit);
  ~PMTHit();

  PMTHit(const PMTHit& right);
  const PMTHit& operator=(const PMTHit& right);

  G4bool operator==(const PMTHit& right) const;

  inline void* operator new(size_t);
  inline void operator delete(void* aHit);

  inline void SetArrivalPos(G4ThreeVector xyz) { fPosArrive = xyz; }
  inline G4ThreeVector GetArrivalPos() { return fPosArrive; }

  inline void SetExitPos(G4ThreeVector xyz) { fPosExit = xyz; }
  inline G4ThreeVector GetExitPos() { return fPosExit; }

  inline void SetArrivalTime(G4double t) { fArrivalTime = t; }
  inline G4double GetArrivalTime() { return fArrivalTime; }

  inline void SetEnergy(G4double en) { fEnergy = en; }
  inline G4double GetEnergy() { return fEnergy; }

  inline void SetEdep(G4double de) { fEdep = de; }
  inline void AddEdep(G4double de) { fEdep += de; }
  inline G4double GetEdep() const { return fEdep; }

  inline void IncPhotonCount() { ++fPhotons; }
  inline G4int GetPhotonCount() { return fPhotons; }

  inline void SetPMTNumber(G4int n) { fPmtNumber = n; }
  inline G4int GetPMTNumber() { return fPmtNumber; }

  inline void SetPMTPhysVol(G4VPhysicalVolume* physVol)
  {
    fPhysVol = physVol;
  }
  inline G4VPhysicalVolume* GetPMTPhysVol() { return fPhysVol; }

  inline void SetPMTPos(G4double x, G4double y, G4double z)
  {
    fPos = G4ThreeVector(x, y, z);
  }

  inline G4ThreeVector GetPMTPos() { return fPos; }

  void Print();

 private:
  // the arrival time of the photon
  G4double fArrivalTime;
  // where the photon hit the detector (detector's coordinate)
  G4ThreeVector fPosArrive;
  // // where the photon exited the fiber (world's coordinate)
  G4ThreeVector fPosExit;
  // energy of photon
  G4double fEnergy;
  G4double fEdep;
  G4int fPmtNumber = -1;
  G4int fPhotons = 0;
  G4ThreeVector fPos;
  G4VPhysicalVolume* fPhysVol = nullptr;
};

//--------------------------------------------------
// Type Definitions
//--------------------------------------------------

using PMTHitsCollection = G4THitsCollection<PMTHit>;

extern G4ThreadLocal G4Allocator<PMTHit>* PMTHitAllocator;

//--------------------------------------------------
// Operator Overloads
//--------------------------------------------------

inline void* PMTHit::operator new(size_t)
{
  if(!PMTHitAllocator)
    PMTHitAllocator = new G4Allocator<PMTHit>;
  return (void*) PMTHitAllocator->MallocSingle();
}

inline void PMTHit::operator delete(void* aHit)
{
  PMTHitAllocator->FreeSingle((PMTHit*) aHit);
}

#endif
