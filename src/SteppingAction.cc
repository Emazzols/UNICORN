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
// $Id: SteppingAction.cc 71007 2013-06-09 16:14:59Z maire $
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "PMTSD.hh"
// #include <sstream>
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
// #include "B4Analysis.hh"
#include "CreateTree.hh"
// #include "struct.hh"

#include "ConfigFile.h"
#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"

#include <iostream>
#include <fstream>
#include <filesystem>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


SteppingAction::SteppingAction(  ConfigFile& config,MyMaterials* fMaterials, MySurfaces* fSurfaces)
: G4UserSteppingAction()
//fEventAction(eventAction)
{
  Materials = fMaterials;
  Surfaces = fSurfaces;
  fScintillationCounter = 0;
  fCerenkovCounter      = 0;
  fPhotonsReachedPMTCounter = 0;
  fEventNumber = -1;

  // quantumEff = config.read<double>("quantumEff");
  logEscape = config.read<int>("log_escaping_particles",0);
  particle_sim = config.read<int>("particle");
  // take materials list save in output struct
  // for(unsigned int i = 0; i < CreateTree::Instance()->materials.size(); i++)
  // {
  //   materialsList.push_back(CreateTree::Instance()->materials[i].name);
  // }

  filename = config.read<std::string>("output","out0");
  filename = std::filesystem::path(filename).stem().string();
  ExtraInfoFile = config.read<int>("ExtraInfoFile",0);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ ; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // G4cout << "Step begin " << G4endl;

  G4int eventNumber = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  if (eventNumber != fEventNumber)
  {
    if(CreateTree::Instance()->verbose > 0)
    {
      G4cout << "Number of Scintillation Photons in previous event: " << fScintillationCounter << G4endl;
      G4cout << "Number of Cerenkov Photons in previous event: " << fCerenkovCounter << G4endl;
      G4cout << "Number of Photons Reached PMT in previous event: " << fPhotonsReachedPMTCounter << G4endl;
    }

    fEventNumber = eventNumber;
    fScintillationCounter = 0;
    fCerenkovCounter = 0;
    fTotalEnDep = 0;
  }



  G4Track* track = step->GetTrack();

  G4String ParticleName = track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
  // G4int PDGEncoding     = track->GetDynamicParticle()->GetParticleDefinition()->GetPDGEncoding();
  G4int TrackID         = track->GetTrackID();
  G4int ParentID        = track->GetParentID();
  G4String materialName = step->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetName();

  G4StepPoint * thePrePoint  = step->GetPreStepPoint();
  G4StepPoint * thePostPoint = step->GetPostStepPoint();

  G4VPhysicalVolume* thePrePV  = step->GetPreStepPoint()->GetPhysicalVolume();
  G4VPhysicalVolume* thePostPV = step->GetPostStepPoint()->GetPhysicalVolume();

  G4String PreMaterialName;
  G4String PostMaterialName;
  G4String PreVolumeName;
  G4String PostVolumeName;
  G4String ProcessDefinedName = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

  if(thePrePV)
  {
    PreMaterialName = thePrePV->GetLogicalVolume()->GetMaterial()->GetName();
    PreVolumeName   = thePrePV->GetName();
  }
  if(thePostPV)
  {
    PostMaterialName = thePostPV->GetLogicalVolume()->GetMaterial()->GetName();
    PostVolumeName   = thePostPV->GetName();
  }

  G4OpBoundaryProcessStatus boundaryStatus=Undefined;
  //getting information from the particles
  //for opticalphoton, we want to know where they stopped
  //find the boundary process only once if it's an opticalphoton
  if (fBoundary == nullptr && ParticleName == "opticalphoton")
  {
    G4ProcessManager* pm = track->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    for(G4int i=0;i<nprocesses;i++)
    {
      if(nullptr != (*pv)[i] && (*pv)[i]->GetProcessName() == "OpBoundary")
      {
        fBoundary = dynamic_cast<G4OpBoundaryProcess*>((*pv)[i]);
        break;
      }
    }
  }

  // if ( (particle_sim != 3 && track->GetParentID() == 0) || (particle_sim == 3) )
  // {
  auto secondaries = step->GetSecondaryInCurrentStep();
  if (nullptr != secondaries && !secondaries->empty())
  {
    for(G4int i=0; i<secondaries->size(); ++i)
    {
      if (secondaries->at(i)->GetParentID()>0)
      {
        if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
        {
          if (secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Scintillation")
          {
            CreateTree::Instance()->NumOptPhotons++;
            fScintillationCounter++;
          }
          else if (secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Cerenkov")
          {
            CreateTree::Instance()->NumCherenkovPhotons++;
            fCerenkovCounter++;
          }
        }
      }
    }
  }
  // }

  if(ParticleName == "opticalphoton")
  {
    if(nullptr != fBoundary)
      boundaryStatus = fBoundary->GetStatus();

    if(thePostPoint->GetStepStatus() == fGeomBoundary)
    {
      switch(boundaryStatus)
      {
        case Absorption:
        case Detection:  // Note, this assumes that the volume causing detection
                         // is the photocathode because it is the only one with
                         // non-zero efficiency
        {
          // Trigger sensitive detector manually since photon is
          // absorbed but status was Detection
          G4SDManager* SDman = G4SDManager::GetSDMpointer();
          G4String sdName    = "/Det/pmtSD";
          PMTSD* pmtSD    = (PMTSD*) SDman->FindSensitiveDetector(sdName);
          if(pmtSD)
            pmtSD->ProcessHits_boundary(step, nullptr);
          // trackInformation->AddTrackStatusFlag(hitPMT);
          break;
        }
        case FresnelReflection:
        case TotalInternalReflection:
        case LambertianReflection:
        case LobeReflection:
        case SpikeReflection:
        case BackScattering:
        default:
          break;
      }
    }

    // bool OnDetectorFlag = false;

    // G4cout << PostVolumeName << G4endl;

    // if(   (PostVolumeName.contains("TILECAL")) && (boundaryStatus == FresnelRefraction)) //whatever enters a detector volume
    // {
    //   // G4cout << "Step should kill the photon " << G4endl;
    //
    //   //kill the photon - it will be recorded later
    //   OnDetectorFlag = true;
    //   track->SetTrackStatus(fStopAndKill);
    // }
    // G4cout << "PreVolumeName = " << PreVolumeName << G4endl;
    if(PreVolumeName == "opticalCoupling_phys" && PostVolumeName == "airGap_phys")
    {
      // G4cout << "Found track to be killed" << G4endl;
      track->SetTrackStatus(fStopAndKill);
    }

    if(PreVolumeName == "airGap_phys" && PostVolumeName == "opticalCoupling_phys")
    {
      // G4cout << "Found track to be killed" << G4endl;
      track->SetTrackStatus(fStopAndKill);
    }

    if(PreVolumeName == "pmtWindow_phys" && PostVolumeName == "airGap_phys")
    {
      // G4cout << "Found track to be killed" << G4endl;
      track->SetTrackStatus(fStopAndKill);
    }

    // if(PreVolumeName == "composite_phys" && PostVolumeName == "airGap_phys")
    // {
    //   // G4cout << "Found track to be killed" << G4endl;
    //   track->SetTrackStatus(fStopAndKill);
    // }

    // if(PreVolumeName == "airGap_phys" && PostVolumeName == "opticalCoupling_phys")
    // {
    //   G4cout << "Found track to be killed" << G4endl;
    //   track->SetTrackStatus(fStopAndKill);
    // }



    // and record
    // if((track->GetTrackStatus() == fStopAndKill) && OnDetectorFlag == true) // || (track->GetTrackStatus() == fKillTrackAndSecondaries)
    // {
    //   G4ThreeVector OnDetectorPosition = track->GetStep()->GetPostStepPoint()->GetPosition(); //get the position vector
    //   G4ThreeVector PreOnDetectorMomentum = track->GetStep()->GetPreStepPoint()->GetMomentumDirection(); //get the momentum, unit vector
    //   G4ThreeVector PostOnDetectorMomentum = track->GetStep()->GetPostStepPoint()->GetMomentumDirection(); //get the momentum, unit vector
    //   G4double globalTime = track->GetGlobalTime();
    //   CreateTree::Instance()->PositionX = OnDetectorPosition.getX();
    //   CreateTree::Instance()->PositionY = OnDetectorPosition.getY();
    //   CreateTree::Instance()->PositionZ = OnDetectorPosition.getZ();
    //   CreateTree::Instance()->PreMomentumX  = PreOnDetectorMomentum.getX();
    //   CreateTree::Instance()->PreMomentumY  = PreOnDetectorMomentum.getY();
    //   CreateTree::Instance()->PreMomentumZ  = PreOnDetectorMomentum.getZ();
    //   CreateTree::Instance()->globalTime = globalTime;
    //   CreateTree::Instance()->NumOptPhotonsReachedPMT++;//Is this the correct way?
    //   fPhotonsReachedPMTCounter++;
    //   if(CreateTree::Instance()->verbose > 0)
    //   {
    //     G4cout << "An optical photon reached the PMT, current count = " << fPhotonsReachedPMTCounter << G4endl;
    //   }
    //
    //   //fPhotonsReachedPMTCounter++;
    //   CreateTree::Instance()->photons->Fill();
    // }

  }
  else //there's only gammas and electrons. for them, we want to know where they left energy (opticalphoton don't really leave any)
  {
  //   if(track->GetTrackID() == 1) // primary particle (gamma)
  //   // Vojtech: But we dont have primary gamma, rather primary electron?
  //   {
  //     if(!CreateTree::Instance()->FoundGammaVertexMomentum)
  //     {
  //       CreateTree::Instance()->SourceMomentumX = (Float_t) track->GetVertexMomentumDirection().getX();
  //       CreateTree::Instance()->SourceMomentumY = (Float_t) track->GetVertexMomentumDirection().getY();
  //       CreateTree::Instance()->SourceMomentumZ = (Float_t) track->GetVertexMomentumDirection().getZ();
  //
  //       // std::ofstream myfile;
  //       // myfile.open ("SourceMomentum.txt",std::ios::app);
  //       // myfile << track->GetVertexMomentumDirection().getX() << " "
  //       //        << track->GetVertexMomentumDirection().getY() << " "
  //       //        << track->GetVertexMomentumDirection().getZ() << std::endl;
  //       // myfile.close();
  //       CreateTree::Instance()->FoundGammaVertexMomentum = true;
  //     }
  //   }

    // no reason to do it only for primary! search for the first particle that enters

    // if(  !CreateTree::Instance()->primaryFirstEntranceFound) // not yet found
    // {
    //   if(PreVolumeName == "World") // coming from exp hall...
    //   {
    //     if(
    //        PostVolumeName == "Module" ||
    //        // PostVolumeName.contains("layer") ||
    //        PostVolumeName == "composite_phys" ||
    //        PostVolumeName == "airGap_phys" ||
    //        // PostVolumeName.contains("spreader") ||
    //        // PostVolumeName.contains("epoxy") ||
    //        PostVolumeName == "photoDet_phys"
    //        // PostVolumeName.contains("esr") ||
    //        // PostVolumeName.contains("air_to")
    //       ) //.. but not staying in exp hall!
    //     {
    //   // if(PreVolumeName.contains("crystal")) // entering now in a crystal
    //   // {
    //       G4ThreeVector particlePosition = thePostPoint->GetPosition();
    //       G4ThreeVector particleMomentum = thePostPoint->GetMomentumDirection();
    //       G4double      particleEnergy   = track->GetTotalEnergy ();
    //
    //       CreateTree::Instance()->primaryPositionOnFirstCrystalX = particlePosition.getX();
    //       CreateTree::Instance()->primaryPositionOnFirstCrystalY = particlePosition.getY();
    //       CreateTree::Instance()->primaryPositionOnFirstCrystalZ = particlePosition.getZ();
    //       CreateTree::Instance()->primaryMomentumOnFirstCrystalX = particleMomentum.getX();
    //       CreateTree::Instance()->primaryMomentumOnFirstCrystalY = particleMomentum.getY();
    //       CreateTree::Instance()->primaryMomentumOnFirstCrystalZ = particleMomentum.getZ();
    //       CreateTree::Instance()->primaryEnergyOnFirstCrystal    = particleEnergy;
    //
    //       CreateTree::Instance()->primaryFirstEntranceFound = true;
    //   // }
    //     }
    //   }
    // }


    G4double edep; 
    edep = step->GetTotalEnergyDeposit()/CLHEP::MeV; //check if there was an energy deposition

    if(edep != 0) //if there was en dep, save the data
    {
        //add total energy deposited
        CreateTree::Instance()->totalEnergyDeposited += edep;
        // fTotalEnDep += edep;
        // totalEnergyDeposited += totalEnergyDeposited;
        // G4cout << "energy deposited = " << fTotalEnDep << " MeV." << G4endl;

        if(ExtraInfoFile == 1)
        {
            const G4TouchableHandle& touchable = step->GetPreStepPoint()->GetTouchableHandle();

            G4Track* currentTrack = step->GetTrack();

            G4ThreeVector pos = thePostPoint->GetPosition();

            G4double stepLength = step->GetStepLength()/CLHEP::mm;

            G4double trackLength = step->GetTrack()->GetTrackLength()/CLHEP::mm;

            auto tree = CreateTree::Instance();

            tree->Edep_vec.push_back(edep);

            tree->X_vec.push_back(pos.x()/CLHEP::mm);
            tree->Y_vec.push_back(pos.y()/CLHEP::mm);
            tree->Z_vec.push_back(pos.z()/CLHEP::mm);

            tree->StepLength_vec.push_back(stepLength);
            tree->TrackLength_vec.push_back(trackLength);
        }
        G4int IsInCase = 0;

        if(PreVolumeName.contains("Capsula"))
        {
          //add total energy deposited
          CreateTree::Instance()->totalEnergyDepositedInExtraVolume += edep;
        }
        if(PreVolumeName.contains("sample"))
        {
          //add total energy deposited
          CreateTree::Instance()->totalEnergyDepositedInSample += edep;
        }
        if(PreVolumeName.contains("case"))
        {
          //add total energy deposited
          CreateTree::Instance()->totalEnergyDepositedInCase += edep;
        }
        if(PreVolumeName.contains("composite_phys0"))
        {
          //add total energy deposited
          CreateTree::Instance()->totalEnergyDepositedInFirstCrystal += edep;
        }
        if(PreVolumeName.contains("composite_phys"))
        {
          //add total energy deposited
          CreateTree::Instance()->totalEnergyDepositedInCrystals += edep;
          // G4cout << "energy deposited in crystal = " << edep << " MeV." << G4endl;
        }
    //
    }
    // if user requires, log the particles (but not optical photons) when they escape from the module, even if they never deposited energy
    // if(logEscape)
    // {
    //   if(PostVolumeName == "World") // going to exp hall...
    //   {
    //     if(PreVolumeName != "World") //.. but not already in exp hall!
    //     {
    //       // log and kill
    //       G4ThreeVector positionVector = track->GetStep()->GetPreStepPoint()->GetPosition(); //get the position vector
    //       G4ThreeVector particleMomentum = thePrePoint->GetMomentumDirection();
    //       G4double globalTime = track->GetGlobalTime();
    //
    //       std::stringstream nameOfTheParticle;
    //       nameOfTheParticle << ParticleName;
    //
    //
    //     }
    //   }
  }




}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
