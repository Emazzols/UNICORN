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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
// #include "PrimaryGeneratorMessenger.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ChargedGeantino.hh"
#include "G4IonTable.hh"
#include "CreateTree.hh"

#include "G4GeneralParticleSource.hh" 
#include "G4SPSEneDistribution.hh"

#include "G4DecayPhysics.hh" 
#include "G4RadioactiveDecayPhysics.hh" 
#include "G4RadioactiveDecay.hh" 
#include "G4NucleusLimits.hh" 
#include "G4ProcessTable.hh" 
#include "G4UAtomicDeexcitation.hh" 
#include "G4DeexPrecoParameters.hh" 
#include "G4NuclearLevelData.hh"

#include <iostream>
#include <fstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(ConfigFile& config)
: G4VUserPrimaryGeneratorAction(),
fParticleGun(0),
fGPS(0)
//    ,fConfig(config)
{
  G4int n_particle = 1;
  source_from_file = config.read<G4int>("source_from_file",0);
  energyFileName = config.read<G4String>("energy_file", "PulseTube_10kV_20mA.txt");
  if(source_from_file==0){
    fGPS = new G4GeneralParticleSource();
    fGPS->SetNumberOfParticles(n_particle);
    G4SPSEneDistribution *energyDist = fGPS->GetCurrentSource()->GetEneDist();
    energyDist->SetEnergyDisType("Arb");
    energyDist->ArbEnergyHistoFile(energyFileName);
    energyDist->ArbInterpolate("Lin");
  } else {
    fParticleGun = new G4ParticleGun(n_particle);
  }
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  // G4int n_particle = 1;
  // fParticleGun = new G4ParticleGun(n_particle);
  // G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();


  //--------------------------------------------------------//
  // GAMMA SOURCE                                           //
  //--------------------------------------------------------//
  //source position and direction
  sourceT5 = config.read<bool>("sourceT5",false);
  if (!sourceT5)
  {
    DBDrandom                  = config.read<double>("DBDrandomize");
    sourcex                    = config.read<double>("sourcex");
    sourcey                    = config.read<double>("sourcey");
    sourcez                    = config.read<double>("sourcez");
    direction_x                = config.read<double>("direction_x");
    direction_y                = config.read<double>("direction_y");
    direction_z                = config.read<double>("direction_z");
  }
  else // set to 0, calc later
  {
    sourcex                    = 0;
    sourcey                    = 0;
    sourcez                    = 0;
    direction_x                = 0;
    direction_y                = 0;
    direction_z                = 0;
  }
  crystal_dim_x              = (config.read<double>("composite_dim_x",50.0)) * mm;
  crystal_dim_y              = (config.read<double>("composite_dim_y",50.0)) * mm;
  crystal_dim_z              = (config.read<double>("composite_dim_z", 6.0)) * mm;
  case_dx                    = (config.read<double>("case_dim_x",50.0)) * mm; 
  case_dy                    = (config.read<double>("case_dim_y",50.0)) * mm; 
  case_dz                    = (config.read<double>("case_dim_z",50.0)) * mm; 
  case_geometry              = config.read<int>("case_geometry", 0); 
  source_dim                 = (config.read<double>("source_dim",50.0)) * mm; 
  multiple_scint             = config.read<int>("multiple_scint", 0); 
  sample_dx                  = (config.read<double>("sample_dim_x", 0)) * mm;
  layers                     = config.read<int>("layers",1);
  gap_between_layers         = (config.read<double>("gap_between_layers",2.0)) * mm;
  air_gap_cry_esr_thickness  = (config.read<double>("air_gap_cry_esr_thickness",0.1)) * mm;
  esr_thickness              = (config.read<double>("esr_thickness",0.1)) * mm;
  layer_dim_z                = crystal_dim_z + 2.0*air_gap_cry_esr_thickness + 2.0*esr_thickness;

  Z                          = config.read<int>("Z",1); 
  A                          = config.read<int>("A",1); 
  full_chain                 = config.read<int>("full_chain", 0);

  module_dx = CreateTree::Instance()->module_dx;
  module_dy = CreateTree::Instance()->module_dy;
  module_dz = CreateTree::Instance()->module_dz;

  x_d = -module_dx/2.0;
  x_u = +module_dx/2.0;
  y_d = -module_dy/2.0;
  y_u = +module_dy/2.0;
  z_d = -module_dz/2.0;
  z_u = +module_dz/2.0;

  opticalCoupling_thickness = config.read<double>("opticalCoupling_thickness", 0.1) * mm; 
  pmtWindowThickness = config.read<double>("pmtWindowThickness", 0.1) * mm; 

  composite_center_x = 0.0;
  composite_center_y =-(opticalCoupling_thickness + pmtWindowThickness)/2.0;
  composite_center_z = 0.0;


  // layer_z = new G4double[layers];
  //
  // for(int iLayer = 0; iLayer < layers ; iLayer++)
  // {
  //   layer_z[iLayer] = (layer_dim_z + gap_between_layers) * (iLayer) - ((layer_dim_z + gap_between_layers)*(layers -1)/2.0);
  // }
  // // take first layer seen by the gamma
  // if(sourcez < 0)
  // {
  //   // it's the first
  //   // first_layer_z = layer_z[0];
  //   first_entrance_z = layer_z[0] - (crystal_dim_z/2.0);
  // }
  // else // cannot be 0, it would be in the middle of a crystal
  // {
  //   // it's the last
  //   // first_layer_z = layer_z[layers-1];
  //   first_entrance_z = layer_z[layers-1] + (crystal_dim_z/2.0);
  // }
  // std::cout << layer_z[0] << " "
  //           << layer_z[layers-1] << " "
  //           << first_entrance_z << " "
  //           << std::endl;




  //energy of gammas
  if (source_from_file == 0) {
    energy = fGPS->GetParticleEnergy();
  }
  else {
    energy = config.read<double>("energy");
  }
  //energy = config.read<double>("energy");
  particle_type = config.read<double>("particle");
  position = config.read<double>("position");

  G4ParticleDefinition* particle;
  if(particle_type == 0)
  {
    particle = particleTable->FindParticle("gamma");
  }
  else if(particle_type == 1)
  {
    particle = particleTable->FindParticle("e-");
  }
  else if(particle_type == 2)                                   
  {            
    particle = particleTable->FindParticle("alpha"); 
  }
  else if(particle_type == 3)                                   
  {                              
    particle = particleTable->FindParticle("geantino");  
  } 

  if (source_from_file==0){
    fGPS->GetCurrentSource()->SetParticleDefinition(particle);
    fGPS->GetCurrentSource()->SetParticleTime(0.0*ns);
    fGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0.,0.,0.));
    fGPS->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(direction_x, direction_y, direction_z));
  }
  else{
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleTime(0.0*ns);
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    fParticleGun->SetParticleEnergy(energy * keV);
  }
  // fParticleGun->SetParticleDefinition(particle);
  // fParticleGun->SetParticleTime(0.0*ns);

  // // default
  // fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  // fParticleGun->SetParticleEnergy(energy*keV);

  std::cout << "end of primary generator ctor" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGPS;
  //   delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  if(source_from_file != 0 && fParticleGun->GetParticleDefinition() == G4Geantino::Geantino())
  {
    G4double IonCharge = 0.*eplus;
    G4double ExcitEnergy = 0.*keV;

    G4ParticleDefinition *ion = G4IonTable::GetIonTable()->GetIon(Z, A, ExcitEnergy);
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(IonCharge);
    fParticleGun->SetParticleTime(0.0*ns);

    if(full_chain==0){
      auto radioactiveDecay = static_cast<G4RadioactiveDecay*>(G4ProcessTable::GetProcessTable()->FindProcess("RadioactiveDecay", "GenericIon"));
      if (radioactiveDecay) {
        G4NucleusLimits nucleusLimits(A, A, Z, Z);
        radioactiveDecay->SetNucleusLimits(nucleusLimits);
        G4DeexPrecoParameters* deex = G4NuclearLevelData::GetInstance()->GetParameters();
        deex->SetCorrelatedGamma(true);
        deex->SetStoreAllLevels(true);
        // deex->SetInternalConversionFlag(true);
        // deex->SetIsomerProduction(true);
      } else {
        G4cerr << "RadioactiveDecay process not found for GenericIon!" << G4endl;
      }
    }
  }

  //now direction
  G4ThreeVector directionVector;
  // std::cout << "sourceT5 = " << sourceT5 << std::endl;
  if(sourceT5)
  {
    // start from center of module, i.e. center of world
    G4double angleLimit;
    double acosMin,vx,vy,vz,randomNum,theta,phi,xc,yc,zc;
    double xp,yp,zp;
    angleLimit = CLHEP::pi;
    acosMin = cos(angleLimit);
    randomNum =  G4UniformRand()*(1.0 - acosMin)  + (acosMin);
    theta = acos(randomNum);
    phi = G4UniformRand() * 2.0 * CLHEP::pi;
    vx = sin(theta)*sin(phi);
    vy = sin(theta)*cos(phi);
    vz = cos(theta);

    // std::cout << vx << " " << vy << " " << vz << std::endl;

    // project on the sphere of radius = 2*module diagonal
    double diagonal = sqrt(pow(module_dx,2)+pow(module_dy,2)+pow(module_dz,2));
    double radius = 2.0*diagonal;

    xp = radius * vx;
    yp = radius * vy;
    zp = radius * vz;
    // std::cout << xp << " " << yp << " " << zp << std::endl;

    // find intersection
    double intersection[3];
    CVec3 B1,B2,L1,L2,Hit;
    // B1 and B2 are the two "extremes" of the module box
    B1.x = -module_dx/2.0;
    B1.y = -module_dy/2.0;
    B1.z = -module_dz/2.0;
    B2.x = +module_dx/2.0;
    B2.y = +module_dy/2.0;
    B2.z = +module_dz/2.0;
    // L1 and L2 are the two extremes of the segment
    L1.x = 0.0;
    L1.y = 0.0;
    L1.z = 0.0;
    L2.x = xp;
    L2.y = yp;
    L2.z = zp;
    // find hit
    int foundHit = CheckLineBox( B1, B2, L2, L1, Hit);
    intersection[0] = Hit.x;
    intersection[1] = Hit.y;
    intersection[2] = Hit.z;
    // std::cout << "found " << foundHit << " "<< intersection[0] << " " << intersection[1] << " " << intersection[2] << std::endl;

    CVec3 normalToSurf;
    // std::cout << "hit " << Hit.x << " " << Hit.y << " " << Hit.z << " "  << std::endl;
    // find normal to intersection surface
    // need tolerance!!!
    float epsilonx = module_dx/1000000.0;
    float epsilony = module_dy/1000000.0;
    float epsilonz = module_dz/1000000.0;
    if((Hit.x > (-module_dx/2.0 - epsilonx)) && (Hit.x < (-module_dx/2.0 + epsilonx)) ) // -x plane
    {
      normalToSurf.x = -module_dx/2.0;
      normalToSurf.y = 0;
      normalToSurf.z = 0;
    }
    if((Hit.x > (module_dx/2.0 - epsilonx)) && (Hit.x < (module_dx/2.0 + epsilonx)) ) // +x plane
    {
      normalToSurf.x = +module_dx/2.0;
      normalToSurf.y = 0;
      normalToSurf.z = 0;
    }
    if((Hit.y > (-module_dy/2.0 - epsilony)) && (Hit.y < (-module_dy/2.0 + epsilony)) ) // -y plane
    {
      normalToSurf.x = 0;
      normalToSurf.y = -module_dy/2.0;
      normalToSurf.z = 0;
    }
    if((Hit.y > (module_dy/2.0 - epsilony)) && (Hit.y < (module_dy/2.0 + epsilony)) ) // +y plane
    {
      normalToSurf.x = 0;
      normalToSurf.y = +module_dy/2.0;
      normalToSurf.z = 0;
    }
    if((Hit.z > (-module_dz/2.0 - epsilonz)) && (Hit.z < (-module_dz/2.0 + epsilonz)) ) // -z plane
    {
      normalToSurf.x = 0;
      normalToSurf.y = 0;
      normalToSurf.z = -module_dz/2.0;
    }
    if((Hit.z > (module_dz/2.0 - epsilonz)) && (Hit.z < (module_dz/2.0 + epsilonz)) ) // +z plane
    {
      normalToSurf.x = 0;
      normalToSurf.y = 0;
      normalToSurf.z = +module_dz/2.0;
    }
    // std::cout << "norm " << normalToSurf.x << " " << normalToSurf.y << " " << normalToSurf.z << " "  << std::endl;

    float dot  = 1.0;
    while(dot > 0)
    {
      // generate random direction
      // produce a random direction from this point
      // use the same variables
      angleLimit = CLHEP::pi;
      acosMin = cos(angleLimit);
      randomNum =  G4UniformRand()*(1.0 - acosMin)  + (acosMin);
      theta = acos(randomNum);
      phi = G4UniformRand() * 2.0 * CLHEP::pi;
      vx = sin(theta)*sin(phi);
      vy = sin(theta)*cos(phi);
      vz = cos(theta);


      CVec3 versor;
      versor.x = vx;
      versor.y = vy;
      versor.z = vz;

      // dot product normal with versor
      dot = versor.Dot(normalToSurf);
      // std::cout << "in dot " << dot << " "<< vx << " " << vy << " " << vz << " "<< theta << " " << phi << std::endl;

      // if dot product between direction and normal is positive, direction points away from module, so produce another direction
    }

    // std::cout << "dot " << dot << " "<< vx << " " << vy << " " << vz << std::endl;

    // set source position and direction
    sourcex = intersection[0];
    sourcey = intersection[1];
    sourcez = intersection[2];
    directionVector = G4ThreeVector(vx,vy,vz);




  }
  else
  {
    if( (direction_x == 0) && (direction_y == 0) && (direction_z == 0) ) // random towards first crystal
    {
      G4double halfDiagonal,distance,angleLimit;
      bool notCrossing = true;
      double acosMin,vx,vy,vz,randomNum,theta,phi,xc,yc,zc;
      while(notCrossing)
      {
        // if(sourcex == 0 && sourcey == 0)// if sourcex  and sourcey == 0, faster way
        // {
          // generate theta
          // find the limit for acos
        //   halfDiagonal = sqrt(pow(crystal_dim_x,2.0) + pow( crystal_dim_y ,2.0) ) / 2.0;
        //   distance = fabs(sourcez - first_entrance_z);
        //   angleLimit = atan(halfDiagonal / distance);
        // }
        // else // otherwise slower
        // {
          // generate theta on all pi, then it will be cut by check of direction and intersection
          // this will be (maybe) a bit slower but will keep being valid also if the source is not in x=0 y=0
          angleLimit = CLHEP::pi;
        // }

        // generate on 4pi, then check direction and crossing
        // angleLimit = CLHEP::pi;
        acosMin = cos(angleLimit);

        //so acos will have to be generated uniformely between acosMin and +1
        randomNum =  G4UniformRand()*(1.0 - acosMin)  + (acosMin);
        theta = acos(randomNum);
        phi = G4UniformRand() * 2.0 * CLHEP::pi;
        vx = sin(theta)*sin(phi);
        vy = sin(theta)*cos(phi);
        vz = cos(theta);

        // G4cout << module_dx/2.0 << " "
        //        << module_dy/2.0 << " "
        //        << module_dz/2.0 << " "
        //        // << vx << " "
        //        // << vy << " "
        //        // << vz << " "
        //        << G4endl;
        //
        // G4cout << sourcex << " "
        //        << sourcey << " "
        //        << sourcez << " "
        //        << vx << " "
        //        << vy << " "
        //        << vz << " "
        //        << G4endl;

        //if source negative, vz needs to be positive, and vice versa
        // if((sourcex < 0) && (vx < 0) )
        // {
        //   continue;
        // }
        // if( (sourcex > 0) && (vx > 0) )
        // {
        //   continue;
        // }
        // //if source negative, vz needs to be positive, and vice versa
        // if((sourcey < 0) && (vy < 0) )
        // {
        //   continue;
        // }
        // if( (sourcey > 0) && (vy > 0) )
        // {
        //   continue;
        // }
        //if source negative, vz needs to be positive, and vice versa
        if((sourcez < 0) && (vz < 0) )
        {
          continue;
        }
        if( (sourcez > 0) && (vz > 0) )
        {
          continue;
        }
        // otherwise...
        // find if it crosses any module face
        // z up
        // G4cout << xc << " "
               // << yc << " "
               // << zc << " "
               // << vx << " "
               // << vy << " "
               // << vz << " "
               // << G4endl;

        xc = sourcex + vx * ((z_u-sourcez)/vz);
        yc = sourcey + vy * ((z_u-sourcez)/vz);
        // G4cout << "zup "
        // << "x "<< xc << " "
        // << "y "<< yc << " "
        // << "z "<< zc << " "
        // << G4endl;
        if( (xc >= -module_dx/2.0) && (xc <= module_dx/2.0) && (yc >= -module_dy/2.0) && (yc <= module_dy/2.0) )
        {
          notCrossing = false;
          break;
        }

        // z down
        xc = sourcex + vx * ((z_d-sourcez)/vz);
        yc = sourcey + vy * ((z_d-sourcez)/vz);

        // G4cout << "zdo "
        // << "x "<< xc << " "
        // << "y "<< yc << " "
        // << "z "<< zc << " "
        // << G4endl;
        if( (xc >= -module_dx/2.0) && (xc <= module_dx/2.0) && (yc >= -module_dy/2.0) && (yc <= module_dy/2.0) )
        {
          notCrossing = false;
          break;
        }
        // y up
        xc = sourcex + vx * ((y_u-sourcey)/vy);
        zc = sourcez + vz * ((y_u-sourcey)/vy);
        // G4cout << "yup "
        // << "x "<< xc << " "
        // << "y "<< yc << " "
        // << "z "<< zc << " "
        // << G4endl;
        if( (zc >= -module_dz/2.0) && (zc <= module_dz/2.0) && (xc >= -module_dx/2.0) && (xc <= module_dx/2.0) )
        {
          notCrossing = false;
          break;
        }

        // y down
        xc = sourcex + vx * ((y_d-sourcey)/vy);
        zc = sourcez + vz * ((y_d-sourcey)/vy);
        // G4cout << "ydo "
        // << "x "<< xc << " "
        // << "y "<< yc << " "
        // << "z "<< zc << " "
        // << G4endl;
        if( (zc >= -module_dz/2.0) && (zc <= module_dz/2.0) && (xc >= -module_dx/2.0) && (xc <= module_dx/2.0) )
        {
          notCrossing = false;
          break;
        }

        // x up
        yc = sourcey + vy * ((x_u-sourcex)/vx);
        zc = sourcez + vz * ((x_u-sourcex)/vx);
        // G4cout << "xup "
        // << "x "<< xc << " "
        // << "y "<< yc << " "
        // << "z "<< zc << " "
        // << G4endl;
        if( (zc >= -module_dz/2.0) && (zc <= module_dz/2.0) && (yc >= -module_dy/2.0) && (yc <= module_dy/2.0) )
        {
          notCrossing = false;
          break;
        }

        // x down
        yc = sourcey + vy * ((x_d-sourcex)/vx);
        zc = sourcez + vz * ((x_d-sourcex)/vx);
        // G4cout << "xdo "
        // << "x "<< xc << " "
        // << "y "<< yc << " "
        // << "z "<< zc << " "
        // << G4endl;
        if( (zc >= -module_dz/2.0) && (zc <= module_dz/2.0) && (yc >= -module_dy/2.0) && (yc <= module_dy/2.0) )
        {
          notCrossing = false;
          break;
        }

        // // otherwise, check if the line crosses the first crystal surf
        // if( (xc >= -crystal_dim_x/2.0) && (xc <= crystal_dim_x/2.0) && (yc >= -crystal_dim_y/2.0) && (yc <= crystal_dim_y/2.0) )
        // {
        //   notCrossing = false;
        // }
      }
      directionVector = G4ThreeVector(vx,vy,vz);


    }
    else // user choice
    {
      directionVector = G4ThreeVector(direction_x,direction_y,direction_z);
    }
  }

  // std::cout << sourcex << " " << sourcey << " " <<sourcez << " ";
  // std::cout << directionVector.x()  << " " << directionVector.y() << " " << directionVector.z()<< std::endl;
  if(DBDrandom)
  {
    if (position == 0){
      sourcex = composite_center_x + (G4UniformRand()- 0.5) * source_dim;
      sourcey = composite_center_y + (G4UniformRand()- 0.5) * source_dim;
      sourcez = composite_center_z + (G4UniformRand()- 0.5) * source_dim;
      //sourcez = case_dx/2-source_dim*sqrt(3)/2+(G4UniformRand()- 0.5) * source_dim;
      direction_x =  (((G4UniformRand()*2))-1);
      direction_y =  (((G4UniformRand()*2))-1);
      direction_z =  (((G4UniformRand()*2))-1);
      directionVector = G4ThreeVector(direction_x,direction_y,direction_z);
    }
    else if (position == 1){
      if(case_geometry == 0){
        sourcex = composite_center_x + (G4UniformRand()- 0.5) * case_dx;
        sourcey = composite_center_y + (G4UniformRand()- 0.5) * case_dy;
        sourcez = composite_center_z + (G4UniformRand()- 0.5) * case_dz;
      }
      else if(case_geometry == 1){

        radius = case_dx / 2.0;
        r = radius * std::pow(G4UniformRand(), 1.0 / 3.0);  
        theta_ang = std::acos(2.0 * G4UniformRand() - 1.0);
        phi_ang = G4UniformRand() * 2.0 * CLHEP::pi;

        // Coordinate cartesiane del punto sorgente
        sourcex = composite_center_x + r * std::sin(theta_ang) * std::cos(phi_ang);
        sourcey = composite_center_y + r * std::sin(theta_ang) * std::sin(phi_ang);
        sourcez = composite_center_z + r * std::cos(theta_ang);
      }
      direction_x =  (((G4UniformRand()*2))-1);
      direction_y =  (((G4UniformRand()*2))-1);
      direction_z =  (((G4UniformRand()*2))-1);
      directionVector = G4ThreeVector(direction_x,direction_y,direction_z);
    }
    else if (position == 2){
      if(multiple_scint == 0){
        if(case_geometry == 0){
          sourcex = composite_center_x + G4UniformRand() * case_dx - case_dx / 2.0;
          sourcey = composite_center_y + G4UniformRand() * case_dy - case_dy / 2.0;
          //sourcez = -(2*case_dz);
        }
        else if(case_geometry == 1){
          radius = case_dx / 2.0;  
          r = radius * std::cbrt(G4UniformRand());
          phi_ang = G4UniformRand() * 2.0 * CLHEP::pi;
          sourcex = composite_center_x + r * std::cos(phi_ang);
          sourcey = composite_center_y + r * std::sin(phi_ang);
          sourcez = composite_center_z - (2.0 * case_dx);
        }
      }
      else if(multiple_scint == 1){
        sourcex = composite_center_x + G4UniformRand() * sample_dx - sample_dx / 2.0;
        sourcey = composite_center_y + G4UniformRand() * sample_dx - sample_dx / 2.0;
        //sourcez = -(2*sample_dx);        
      }

      direction_x = 0;
      direction_y = 0;
      direction_z = 1;

      directionVector = G4ThreeVector(direction_x, direction_y, direction_z);
    }
  }

  if (source_from_file==0){
    fGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(sourcex,sourcey,sourcez));
    fGPS->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(directionVector);
    //G4double generatedEnergy = fGPS->GetParticleEnergy();
    //G4cout << "Generated energy: " << generatedEnergy << " MeV" << G4endl;
    fGPS->GeneratePrimaryVertex(anEvent);
  }
  else{
    fParticleGun->SetParticlePosition(G4ThreeVector(sourcex,sourcey,sourcez));
    fParticleGun->SetParticleMomentumDirection(directionVector);
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
  // fParticleGun->SetParticlePosition(G4ThreeVector(sourcex,sourcey,sourcez));
  // fParticleGun->SetParticleMomentumDirection(directionVector);
  // fParticleGun->GeneratePrimaryVertex(anEvent);
}



int PrimaryGeneratorAction::GetIntersection( float fDst1, float fDst2, CVec3 P1, CVec3 P2, CVec3 &Hit) {
  if ( (fDst1 * fDst2) >= 0.0f) return 0;
  if ( fDst1 == fDst2) return 0;
  Hit = P1 + (P2-P1) * ( -fDst1/(fDst2-fDst1) );
  return 1;
}

int PrimaryGeneratorAction::InBox( CVec3 Hit, CVec3 B1, CVec3 B2, const int Axis) {
  if ( Axis==1 && Hit.z > B1.z && Hit.z < B2.z && Hit.y > B1.y && Hit.y < B2.y) return 1;
  if ( Axis==2 && Hit.z > B1.z && Hit.z < B2.z && Hit.x > B1.x && Hit.x < B2.x) return 1;
  if ( Axis==3 && Hit.x > B1.x && Hit.x < B2.x && Hit.y > B1.y && Hit.y < B2.y) return 1;
  return 0;
}

// returns true if line (L1, L2) intersects with the box (B1, B2)
// returns intersection point in Hit
int PrimaryGeneratorAction::CheckLineBox( CVec3 B1, CVec3 B2, CVec3 L1, CVec3 L2, CVec3 &Hit)
{
  if (L2.x < B1.x && L1.x < B1.x) return false;
  if (L2.x > B2.x && L1.x > B2.x) return false;
  if (L2.y < B1.y && L1.y < B1.y) return false;
  if (L2.y > B2.y && L1.y > B2.y) return false;
  if (L2.z < B1.z && L1.z < B1.z) return false;
  if (L2.z > B2.z && L1.z > B2.z) return false;
  if (L1.x > B1.x && L1.x < B2.x && L1.y > B1.y && L1.y < B2.y && L1.z > B1.z && L1.z < B2.z){
    Hit = L1;
    return true;
  }
  if ( (GetIntersection( L1.x-B1.x, L2.x-B1.x, L1, L2, Hit) && InBox( Hit, B1, B2, 1 ))
  || (GetIntersection( L1.y-B1.y, L2.y-B1.y, L1, L2, Hit) && InBox( Hit, B1, B2, 2 ))
  || (GetIntersection( L1.z-B1.z, L2.z-B1.z, L1, L2, Hit) && InBox( Hit, B1, B2, 3 ))
  || (GetIntersection( L1.x-B2.x, L2.x-B2.x, L1, L2, Hit) && InBox( Hit, B1, B2, 1 ))
  || (GetIntersection( L1.y-B2.y, L2.y-B2.y, L1, L2, Hit) && InBox( Hit, B1, B2, 2 ))
  || (GetIntersection( L1.z-B2.z, L2.z-B2.z, L1, L2, Hit) && InBox( Hit, B1, B2, 3 )))
  return true;

  return false;
}
