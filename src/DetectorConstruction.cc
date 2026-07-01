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


#include "DetectorConstruction.hh"
#include "G4SDManager.hh"
#include "PMTSD.hh"
#include "ScintSD.hh"

#include "G4UserLimits.hh"
#include "Randomize.hh"
#include "G4Exception.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "CreateTree.hh"
#include "ConfigFile.h"
#include "G4NistManager.hh"
#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(ConfigFile& config,MyMaterials* Materials, MySurfaces* Surfaces)
: G4VUserDetectorConstruction(),
fConfig(config),
fMaterials(Materials),
fSurfaces(Surfaces)
{
  //Set default physical dimensions
  //all dimensions in G4 are semi dimensions (so smart...)
  //so dimensions are defined here, but will have to be divided by 2 when creating the objects
  fCheckOverlaps = true;

  manual_step   = config.read<int>("manual_step",0); 
  maxStep = (config.read<double>("max_step_size",1.0)) * mm;

  sourcex                    = config.read<double>("sourcex");
  sourcey                    = config.read<double>("sourcey");
  sourcez                    = config.read<double>("sourcez");

  world_material        = config.read<int>("world_material");

  composite_material        = config.read<int>("composite_material");
  opticalCoupling_material  = config.read<int>("opticalCoupling_material");
  G4cout << "Gap to SiPM material selected as: " << opticalCoupling_material << G4endl;

  number_of_detectors    = config.read<int>("number_of_detectors");
  photoDetSize = (config.read<double>("photoDetSize")) * mm; // active area of the PMT only
  pmtWindowThickness = (config.read<double>("pmtWindowThickness")) * mm;
  photoCathThickness = (config.read<double>("photoCathThickness")) * mm;
  G4cout << "Number of detectors selected as: " << number_of_detectors << G4endl;

  // Case
  case_material = config.read<int>("case_material", 0); 
  case_geometry = config.read<int>("case_geometry", 0);
  case_dx = (config.read<double>("case_dim_x",50.0)) * mm; 
  case_dy = (config.read<double>("case_dim_y",50.0)) * mm; 
  case_dz = (config.read<double>("case_dim_z",50.0)) * mm; 

  // Structure
  structure = config.read<int>("structure", 0); 

  multiple_scint = config.read<int>("multiple_scint", 0); 
  sample_dx = (config.read<double>("sample_dim_x", 0)) * mm; 
  numCoordsCase = config.read<int>("scint_num", 1); 
  sample_material       = config.read<int>("sample_material",0);

  numCoords = config.read<int>("composite_num", 1); 

  numCoords_x = config.read<int>("composite_num_x", 1); 
  numCoords_y = config.read<int>("composite_num_y", 1); 
  numCoords_z = config.read<int>("composite_num_z", 1);

  composite_dx = (config.read<double>("composite_dim_x")) * mm;
  composite_dy = (config.read<double>("composite_dim_y")) * mm;
  composite_dz = (config.read<double>("composite_dim_z")) * mm;

  composite_pos = (config.read<double>("composite_pos",0));

  capsula = config.read<int>("capsula",0); 
  capsula_dx = (config.read<double>("capsula_dim_x",50.0)) * mm; 
  capsula_dy = (config.read<double>("capsula_dim_y",50.0)) * mm; 
  capsula_dz = (config.read<double>("capsula_dim_z", 5)) * mm; 

  capsulax                    = config.read<double>("capsulax"); 
  capsulay                    = config.read<double>("capsulay"); 
  capsulaz                    = config.read<double>("capsulaz");

  // GAPS dimensions
  opticalCoupling_thickness    = (config.read<double>("opticalCoupling_thickness")) * mm;
  // epoxy_thickness          = (config.read<double>("epoxy_thickness")) * mm;
  airGapThickness = (config.read<double>("airGapThickness")) * mm;
  reflectorThickness = (config.read<double>("reflectorThickness")) * mm;
  // total_internal_reflection = config.read<int>("total_internal_reflection");

  opticalCoupling_refr_index    = config.read<double>("opticalCoupling_refr_index");

  // fastrisetime  = config.read<double>("fastrisetime");
  // fastdecaytime = config.read<double>("fastdecaytime");
  // fastratio     = config.read<double>("fastratio");
  // slowrisetime  = config.read<double>("slowrisetime");
  // slowdecaytime = config.read<double>("slowdecaytime");

  // Material properties of the Nanocomposite

  host_fractMass = config.read<double>("host_fractMass");
  NP_fractMass = config.read<double>("NP_fractMass");
  lightyield = config.read<double>("lightyield");
  //G4cout << "Scintillator Light Yield [Ph/MeV]: " << lightyield << G4endl;
  risetime  = config.read<double>("risetime");
  decaytime = config.read<double>("decaytime");
  rayleigh_length = config.read<double>("rayleigh_length");
  mie_length = config.read<double>("mie_length");

  if((host_fractMass + NP_fractMass) != 100 )
  {
    G4cerr << " ERROR: (host_fractMass + NP_fractMass) MUST BE = 100 !!!!" << G4endl;
    G4cerr << " Aborting..." << G4endl;
    exit(0);
  }

  G4cout << G4endl;
  G4cout << "******************************************************************** " << G4endl;
  G4cout << "*                          Nanocomposite properties                * " << G4endl;
  G4cout << "* host fractional mass         =                              " << host_fractMass << " % *"<< G4endl;
  G4cout << "* nanoparticle fractional mass =                              " << NP_fractMass << " % *"<< G4endl;
  G4cout << "* nanoparticle light yield     =                       " << lightyield << " ph/MeV *"<< G4endl;
  G4cout << "* nanoparticle rise time       =                          " << risetime << " ns *"<< G4endl;
  G4cout << "* nanoparticle decay time      =                            " << decaytime << " ns *"<< G4endl;
  G4cout << "* rayleigh scattering length for 500 nm photons =            " << rayleigh_length << " mm *"<< G4endl;
  G4cout << "******************************************************************** " << G4endl;
  G4cout << G4endl;



  // energy resolution
  resolutionScale = config.read<double>("resolutionScale");

  esr_type = config.read<int>("esr_type",0);
  simplified_esr_transmittance = config.read<double>("simplified_esr_transmittance");
  simplified_esr_reflectivity  = config.read<double>("simplified_esr_reflectivity");
  simplified_esr_efficiency    = config.read<double>("simplified_esr_efficiency");

  if((simplified_esr_transmittance + simplified_esr_reflectivity) > 1.0)
  {
    G4cerr << " ERROR: (simplified_esr_transmittance + simplified_esr_reflectivity) MUST BE <= 1 !!!!" << G4endl;
    G4cerr << " Aborting..." << G4endl;
    exit(0);
  }

  esrTransmittance = config.read<double>("esrTransmittance",-1);
  if(esrTransmittance == -1)
  {
    G4cout << "Using esr transmittance values from measurement" << G4endl;
  }
  else
  {
    G4cout << "Esr transmittance fixed to " << esrTransmittance << " for all wavelenghts" << G4endl;
  }

  // depolished     = config.read<int>   ("depolished"    ,0);
  // sigmaalpha     = config.read<double>("sigmaalpha"    ,0.05);
  // specularlobe   = config.read<double>("specularlobe"  ,1);
  // specularspike  = config.read<double>("specularspike" ,0);
  // backscattering = config.read<double>("backscattering",0);

  // G4cout << G4endl;
  // G4cout << ">>> Surface state:" << G4endl;
  // if(depolished == 0 ) G4cout << "polished" << G4endl;
  // if(depolished == 1 ) G4cout << "depolished" << G4endl;
  // G4cout << "sigmaalpha     = " << sigmaalpha     << G4endl;
  // G4cout << "specularlobe   = " << specularlobe   << G4endl;
  // G4cout << "specularspike  = " << specularspike  << G4endl;
  // G4cout << "backscattering = " << backscattering << G4endl;
  // G4cout << G4endl;

  // sides
  //          0
  //        ______
  // y     |      |
  // ^  3  |      | 1
  // |     |______|
  // |        2
  // -----> x
  //

  // in the detector phy log etc we will build the PMT line starting from one edge of the side and then iterativly moving to the right
  // we have to define the starting position and the step dimension as a function of the chPerSide

  if(multiple_scint==0){
    module_dx = case_dx;
    module_dy = case_dy;
    module_dz = case_dz;
  }
  else if(multiple_scint==1){
    module_dx = sample_dx;
    module_dy = sample_dx;
    module_dz = sample_dx;
  }

  PMTtotalGap = module_dx-photoDetSize;
  PMTgapBetwCenters = PMTtotalGap+photoDetSize; // an entire gap is being kept between PMT and half a gap between PMT and the corner of the tile
  PMTbuildStartPoint = module_dx/2 - photoDetSize/2 - ((PMTgapBetwCenters-photoDetSize)/2);

  if(PMTtotalGap<0){std::cout << "####  Warning! PMTs overlap, decrese number of chPerSide of the dimension of the PMT.  ####" << std::endl;};
  switch(number_of_detectors)
  {
    case 1:
      module_x = 0;
      module_y = -(opticalCoupling_thickness + pmtWindowThickness)/2;
      module_z = 0;

      opticalCoupling_dx.push_back(photoDetSize);
      opticalCoupling_dy.push_back(opticalCoupling_thickness);
      opticalCoupling_dz.push_back(photoDetSize);

      opticalCoupling_x.push_back(0);
      opticalCoupling_y.push_back(module_dy/2 - (pmtWindowThickness)/2);
      opticalCoupling_z.push_back(0);

      pmtWindow_x.push_back(-PMTbuildStartPoint);
      pmtWindow_y.push_back(module_dy/2.0 + opticalCoupling_thickness/2); // + epoxy_dy[0]
      pmtWindow_z.push_back(-PMTbuildStartPoint);

      photoCath_x.push_back(-PMTbuildStartPoint);
      photoCath_y.push_back(pmtWindowThickness/2.0-photoCathThickness/2.0); // + epoxy_dy[0]
      photoCath_z.push_back(-PMTbuildStartPoint);

      airGap_dx = module_dx + airGapThickness;
      airGap_dy = module_dy + opticalCoupling_thickness + pmtWindowThickness;
      airGap_dz = module_dz + airGapThickness;

      reflector_dx = module_dx + airGapThickness + reflectorThickness;
      reflector_dy = module_dy + opticalCoupling_thickness  + pmtWindowThickness+ reflectorThickness;
      reflector_dz = module_dz + airGapThickness + reflectorThickness;
      break;
    case 2:
      opticalCoupling_dx.push_back(photoDetSize);
      opticalCoupling_dx.push_back(photoDetSize);

      opticalCoupling_dy.push_back(opticalCoupling_thickness);
      opticalCoupling_dy.push_back(opticalCoupling_thickness);

      opticalCoupling_dz.push_back(photoDetSize);
      opticalCoupling_dz.push_back(photoDetSize);

      opticalCoupling_x.push_back(0);
      opticalCoupling_x.push_back(0);

      opticalCoupling_y.push_back(module_dy/2.0 + opticalCoupling_thickness/2.0);
      opticalCoupling_y.push_back(-(module_dy/2.0  + opticalCoupling_dy[1]/2.0));

      opticalCoupling_z.push_back(0);
      opticalCoupling_z.push_back(0);

      pmtWindow_x.push_back(-PMTbuildStartPoint);
      pmtWindow_x.push_back(-PMTbuildStartPoint);

      pmtWindow_y.push_back(module_dy/2.0  + opticalCoupling_thickness  + pmtWindowThickness/2.0); // + epoxy_dy[0]
      pmtWindow_y.push_back(-(module_dy/2.0 + opticalCoupling_thickness + pmtWindowThickness/2.0)); // + epoxy_dy[2]

      pmtWindow_z.push_back(-PMTbuildStartPoint);
      pmtWindow_z.push_back(-PMTbuildStartPoint);

      photoCath_x.push_back(-PMTbuildStartPoint);
      photoCath_x.push_back(-PMTbuildStartPoint);

      photoCath_y.push_back(pmtWindowThickness/2.0-photoCathThickness/2.0); // + epoxy_dy[0]
      photoCath_y.push_back(-(pmtWindowThickness/2.0-photoCathThickness/2.0)); // + epoxy_dy[2]

      photoCath_z.push_back(-PMTbuildStartPoint);
      photoCath_z.push_back(-PMTbuildStartPoint);

      airGap_dx = module_dx + airGapThickness;
      airGap_dy = module_dy + 2*opticalCoupling_thickness + 2*pmtWindowThickness;
      airGap_dz = module_dz + airGapThickness;

      reflector_dx = module_dx + airGapThickness + reflectorThickness;
      reflector_dy = module_dy + 2*opticalCoupling_thickness + 2*pmtWindowThickness + 2*reflectorThickness;
      reflector_dz = module_dz + airGapThickness + reflectorThickness;
      break;
    case 6:
      opticalCoupling_dx.push_back(photoDetSize);
      opticalCoupling_dx.push_back(opticalCoupling_thickness);
      opticalCoupling_dx.push_back(photoDetSize);
      opticalCoupling_dx.push_back(opticalCoupling_thickness);
      opticalCoupling_dx.push_back(photoDetSize);
      opticalCoupling_dx.push_back(photoDetSize);

      opticalCoupling_dy.push_back(opticalCoupling_thickness);
      opticalCoupling_dy.push_back(photoDetSize);
      opticalCoupling_dy.push_back(opticalCoupling_thickness);
      opticalCoupling_dy.push_back(photoDetSize);
      opticalCoupling_dy.push_back(photoDetSize);
      opticalCoupling_dy.push_back(photoDetSize);

      opticalCoupling_dz.push_back(photoDetSize);
      opticalCoupling_dz.push_back(photoDetSize);
      opticalCoupling_dz.push_back(photoDetSize);
      opticalCoupling_dz.push_back(photoDetSize);
      opticalCoupling_dz.push_back(opticalCoupling_thickness);
      opticalCoupling_dz.push_back(opticalCoupling_thickness);

      opticalCoupling_x.push_back(0);
      opticalCoupling_x.push_back(module_dx/2.0 + opticalCoupling_dx[1]/2.0);
      opticalCoupling_x.push_back(0);
      opticalCoupling_x.push_back(-(module_dx/2.0 + opticalCoupling_dx[3]/2.0));
      opticalCoupling_x.push_back(0);
      opticalCoupling_x.push_back(0);

      opticalCoupling_y.push_back(module_dy/2.0 + opticalCoupling_dy[0]/2.0);
      opticalCoupling_y.push_back(0);
      opticalCoupling_y.push_back(-(module_dy/2.0 + opticalCoupling_dy[2]/2.0));
      opticalCoupling_y.push_back(0);
      opticalCoupling_y.push_back(0);
      opticalCoupling_y.push_back(0);

      opticalCoupling_z.push_back(0);
      opticalCoupling_z.push_back(0);
      opticalCoupling_z.push_back(0);
      opticalCoupling_z.push_back(0);
      opticalCoupling_z.push_back(-(module_dz/2.0 + opticalCoupling_dz[4]/2.0));
      opticalCoupling_z.push_back(module_dz/2.0 + opticalCoupling_dz[5]/2.0);

      pmtWindow_x.push_back(-PMTbuildStartPoint);
      pmtWindow_x.push_back(module_dx/2.0  + opticalCoupling_thickness + pmtWindowThickness/2.0);
      pmtWindow_x.push_back(-PMTbuildStartPoint);
      pmtWindow_x.push_back(-(module_dx/2.0  + opticalCoupling_thickness  + pmtWindowThickness/2.0));
      pmtWindow_x.push_back(-PMTbuildStartPoint);
      pmtWindow_x.push_back(-PMTbuildStartPoint);

      pmtWindow_y.push_back(module_dy/2.0  + opticalCoupling_thickness  + pmtWindowThickness/2.0);
      pmtWindow_y.push_back(-PMTbuildStartPoint);
      pmtWindow_y.push_back(-(module_dy/2.0 + opticalCoupling_thickness + pmtWindowThickness/2.0));
      pmtWindow_y.push_back(-PMTbuildStartPoint);
      pmtWindow_y.push_back(-PMTbuildStartPoint);
      pmtWindow_y.push_back(-PMTbuildStartPoint);

      pmtWindow_z.push_back(-PMTbuildStartPoint);
      pmtWindow_z.push_back(-PMTbuildStartPoint);
      pmtWindow_z.push_back(-PMTbuildStartPoint);
      pmtWindow_z.push_back(-PMTbuildStartPoint);
      pmtWindow_z.push_back(module_dz/2.0  + opticalCoupling_thickness + pmtWindowThickness/2.0);
      pmtWindow_z.push_back(-(module_dz/2.0  + opticalCoupling_thickness + pmtWindowThickness/2.0));

      photoCath_x.push_back(-PMTbuildStartPoint);
      photoCath_x.push_back(pmtWindowThickness/2.0-photoCathThickness/2.0);
      photoCath_x.push_back(-PMTbuildStartPoint);
      photoCath_x.push_back(-(pmtWindowThickness/2.0-photoCathThickness/2.0));
      photoCath_x.push_back(-PMTbuildStartPoint);
      photoCath_x.push_back(-PMTbuildStartPoint);

      photoCath_y.push_back(pmtWindowThickness/2.0-photoCathThickness/2.0);
      photoCath_y.push_back(-PMTbuildStartPoint);
      photoCath_y.push_back(-(pmtWindowThickness/2.0-photoCathThickness/2.0));
      photoCath_y.push_back(-PMTbuildStartPoint);
      photoCath_y.push_back(-PMTbuildStartPoint);
      photoCath_y.push_back(-PMTbuildStartPoint);

      photoCath_z.push_back(-PMTbuildStartPoint);
      photoCath_z.push_back(-PMTbuildStartPoint);
      photoCath_z.push_back(-PMTbuildStartPoint);
      photoCath_z.push_back(-PMTbuildStartPoint);
      photoCath_z.push_back(pmtWindowThickness/2.0-photoCathThickness/2.0);
      photoCath_z.push_back(-(pmtWindowThickness/2.0-photoCathThickness/2.0));

      airGap_dx = module_dx + 2*opticalCoupling_thickness + 2*pmtWindowThickness;
      airGap_dy = module_dy + 2*opticalCoupling_thickness + 2*pmtWindowThickness;
      airGap_dz = module_dz + 2*opticalCoupling_thickness + 2*pmtWindowThickness;

      reflector_dx = module_dx + 2*opticalCoupling_thickness + 2*pmtWindowThickness + 2*reflectorThickness;
      reflector_dy = module_dy + 2*opticalCoupling_thickness + 2*pmtWindowThickness + 2*reflectorThickness;
      reflector_dz = module_dz + 2*opticalCoupling_thickness + 2*pmtWindowThickness + 2*reflectorThickness;
      break;
  }


  // module_dx = composite_dx;
  // module_dy = composite_dy;
  // module_dz = composite_dz;

  // save module dimensions
  CreateTree::Instance()->module_dx = module_dx;
  CreateTree::Instance()->module_dy = module_dy;
  CreateTree::Instance()->module_dz = module_dz;

  // make the world larger than N times the module, or large enough to contain the source...

  fExpHall_dx = std::max(fabs(module_dx * 10),fabs(sourcex*4));
  fExpHall_dy = std::max(fabs(module_dy * 10),fabs(sourcey*4));
  fExpHall_dz = std::max(fabs(module_dz * 10),fabs(sourcez*4));



  //-------------------------------------------------------------------//
  //                                                                   //
  //                             MATERIALS                             //
  //                                                                   //
  //-------------------------------------------------------------------//
  if(composite_material == 0) // LYSO
  {
    CompositeMaterial = Materials->LYSO(lightyield, resolutionScale);
  }
  else if(composite_material == 1) // CdSe
  {
    // CompositeMaterial = Materials->CuDopedCdSNCs(lightyield,
    //                                              risetime,
    //                                              decaytime,
    //                                              rayleigh_length,
    //                                              mie_length);
    CompositeMaterial = Materials->CdSeCdSinPMMA(host_fractMass,
                                                 NP_fractMass,
                                                 lightyield,
                                                 risetime,
                                                 decaytime,
                                                 rayleigh_length,
                                                 mie_length);
    //
    // G4cout << "DotInPlateCdSeS chosen as scintillator material" << G4endl;

  }
  else if(composite_material == 2) // BC422
  {
    CompositeMaterial = Materials->BC422(lightyield,resolutionScale);
  }
  else
  {
    G4cout << "ERROR: No scintillator material chosen???" << G4endl;
    // return 1;
  }

  if(case_material == 0){ 
    CaseMaterial = Materials->Plastic();
  }
  else if(case_material == 1){
    CaseMaterial = Materials->FusedSilica();
  }


  if(sample_material == 0){ 
    SampleMaterial = Materials->Plastic();
  }
  else if(sample_material == 1){
    SampleMaterial = Materials -> FusedSilica();
  }

  if(opticalCoupling_material == 0) // air
  {
    OpticalCouplingMaterial = Materials->AirThinLayer();
  }
  else if(opticalCoupling_material == 1) // grease
  {
    OpticalCouplingMaterial =  Materials->OpticalGrease(opticalCoupling_refr_index);
    G4cout << "OpticalCouplingMaterial chosen as: OpticalGrease" << G4endl;
  }
  else
  {
    G4cout << "ERROR: No OpticalCouplingMaterial material chosen???" << G4endl;
    // return 1;
  }

  if(world_material == 0) 
  {
    WorldMaterial = Materials->Vacuum();
  }
  else if(world_material == 1)
  {
    WorldMaterial = Materials->Air();
  }
  PhotoCathodeMaterial = Materials->Aluminum();
  PmtWindowMaterial = Materials->FusedSilica();
  AirGapMaterial   = Materials->AirThinLayer();
  ReflectorMaterial = Materials->Plastic();


  //-------------------------------------------------------------------//
  //                                                                   //
  //                             REFLECTOR                             //
  //                                                                   //
  //-------------------------------------------------------------------//
  // choose ESR

  if(esr_type == 0) // Standard ESR
  {
    ESR_reflector = Surfaces->ESR(esrTransmittance);
    G4cout << "Reflector type chosen as esr_type 0" << G4endl;
  }
  else if(esr_type == 1)
  {
    ESR_reflector = Surfaces->Simplified_ESR(simplified_esr_reflectivity,simplified_esr_efficiency);
    G4cout << "Reflector type chosen as esr_type 1:" << G4endl;
    G4cout << "simplified_esr_reflectivity = " << simplified_esr_reflectivity << G4endl;
    G4cout << "simplified_esr_efficiency = " << simplified_esr_efficiency << G4endl;


  }
  else
  {
    G4cout << "ERROR: No ESR type chosen???" << G4endl;
    // return 1;
  }


  // if(total_internal_reflection == 0)
  // {
  //   G4cout << "Using ESR on both composite and light guides" << G4endl;
  // }
  // else
  // {
  //   G4cout << "ESR only on light guides" << G4endl;
  // }

  //-------------------------------------------------------------------//
  //                                                                   //
  //                             SURFACES                              //
  //                                                                   //
  //-------------------------------------------------------------------//
  // if(depolished == 1)
  // {
  //   DepolishedSurface = Surfaces->DepolishedSurface(specularlobe,specularspike,backscattering);
  // }



  //-------------------------------------------------------------------//
  //                                                                   //
  //                             VISUALIZATION                         //
  //                                                                   //
  //-------------------------------------------------------------------//
  worldVis            = new G4VisAttributes(G4Colour(0.0,0.0,0.0));                 // white
  compositeVis        = new G4VisAttributes(G4Colour(0.0,0.0,1.0));                 // cyan
  opticalCouplingVis  = new G4VisAttributes(G4Colour(0.0,1.0,0.0)); opticalCouplingVis->SetVisibility(true);              // green
  pmtWindowVis         = new G4VisAttributes(G4Colour(0.0,0.0,0.0));                 // red
  photoCathVis        = new G4VisAttributes(G4Colour(1.0,0.0,1.0));  photoCathVis->SetVisibility(true);               // red
  airGapVis           = new G4VisAttributes(G4Colour(0.2,0.1,0.6)); airGapVis->SetVisibility(true);
  reflectorVis        = new G4VisAttributes(G4Colour(0.6,0.1,0.4)); reflectorVis->SetVisibility(true);
  capsulaVis         = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  sampleVis          = new G4VisAttributes(G4Colour(0.0,1.0,1.0));                 // cyan
  caseVis          = new G4VisAttributes(G4Colour(0.0,0.0,0.0));

  materialsList = Materials->GetMaterialsList();
  std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<" << std::endl;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::vector<G4ThreeVector>& CompositeParameterisation::GetPositions() const {
    return positions;
}

CompositeParameterisation::CompositeParameterisation(G4double case_dim_x, G4double case_dim_y, G4double case_dim_z, G4int num_composite_x, G4int num_composite_y, G4int num_composite_z, G4double composite_dim) : case_dim_x(case_dim_x), case_dim_y(case_dim_y), case_dim_z(case_dim_z), composite_dim(composite_dim){

  G4double initial_offset_x = -case_dim_x/2 + composite_dim;
  G4double initial_offset_y = -case_dim_y/2 + composite_dim;
  G4double initial_offset_z = -case_dim_z/2 + composite_dim;

  G4double spacing_x = (num_composite_x > 1) ? (case_dim_x - 2*composite_dim) / (num_composite_x - 1) : 0;
  G4double spacing_y = (num_composite_y > 1) ? (case_dim_y - 2*composite_dim) / (num_composite_y - 1) : 0;
  G4double spacing_z = (num_composite_z > 1) ? (case_dim_z - 2*composite_dim) / (num_composite_z - 1) : 0;

  for (int i = 0; i < num_composite_x; ++i) {
      for (int j = 0; j < num_composite_y; ++j) {
          for (int k = 0; k < num_composite_z; ++k) {
              G4double x = initial_offset_x + i * spacing_x;
              G4double y = initial_offset_y + j * spacing_y;
              G4double z = initial_offset_z + k * spacing_z;

              positions.push_back(G4ThreeVector(x, y, z));                
          }
      }
  }
}


void CompositeParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const {
    // Posiziona ogni cristallo nella posizione generata
    G4ThreeVector position = positions[copyNo];
    physVol->SetTranslation(position);
}

void CompositeParameterisation::ComputeDimensions(G4Box& box, const G4int copyNo, const G4VPhysicalVolume* physVol) const {
    // Imposta le dimensioni del cristallo
    box.SetXHalfLength(composite_dim / 2.0);
    box.SetYHalfLength(composite_dim / 2.0);
    box.SetZHalfLength(composite_dim / 2.0);
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{

  if(manual_step == 1){
    maxStepLimit = new G4UserLimits(maxStep);
  }

  // id of detectors
  G4int detector_id = 0;

  // The experimental Hall
  expHall_box  = new G4Box("World",fExpHall_dx/2.0,fExpHall_dy/2.0,fExpHall_dz/2.0);
  expHall_log  = new G4LogicalVolume(expHall_box,WorldMaterial,"World",0,0,0);
  expHall_phys = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0,fCheckOverlaps);
  expHall_log->SetVisAttributes(worldVis);   // set vis color

  reflector_box = new G4Box("reflector_box", reflector_dx/2.0, reflector_dy/2.0, reflector_dz/2.0);
  reflector_log = new G4LogicalVolume(reflector_box,AirGapMaterial,"reflector_log",0,0,0);
  reflector_phys = new G4PVPlacement(0,G4ThreeVector(reflector_x,reflector_y,reflector_z),reflector_log,"reflector_phys",expHall_log,false,0,0);
  reflector_log->SetVisAttributes(reflectorVis);

  airGap_box = new G4Box("airGap_box", airGap_dx/2.0, airGap_dy/2.0, airGap_dz/2.0);
  airGap_log = new G4LogicalVolume(airGap_box,AirGapMaterial,"airGap_log",0,0,0);
  airGap_phys = new G4PVPlacement(0,G4ThreeVector(airGap_x,airGap_y,airGap_z),airGap_log,"airGap_phys",reflector_log,false,0,0);
  airGap_log->SetVisAttributes(airGapVis);

  // composite_box = new G4Box("composite_box",composite_dx/2.0,composite_dy/2.0,composite_dz/2.0);
  // composite_log = new G4LogicalVolume(composite_box,CompositeMaterial,"composite_log",0,0,0);
  // composite_phys = new G4PVPlacement(0,G4ThreeVector(composite_x,composite_y,composite_z),composite_log,"composite_phys",airGap_log,false,0,fCheckOverlaps);
  // composite_log->SetVisAttributes(compositeVis);   // set vis color

  if (structure == 0){
    if (case_geometry == 0){
      composite_box  = new G4Box("composite_box",composite_dx/2.0,composite_dy/2.0,composite_dz/2.0);
    }
    else if(case_geometry == 1){ 
      composite_box  = new G4Sphere("composite_box", 0, composite_dx/2.0, 0.0 * deg, 360.0 * deg, 0.0 * deg, 180.0 * deg);
    } 
    composite_log  = new G4LogicalVolume(composite_box,CompositeMaterial,"composite_log",0,0,0);
    composite_phys = new G4PVPlacement(0,G4ThreeVector(module_x,module_y,module_z),composite_log,"composite_phys",airGap_log,false,0,fCheckOverlaps);
    composite_log -> SetVisAttributes(compositeVis);   // set vis color
    // compositeVis -> SetForceSolid(true);
    if(manual_step == 1){ 
      composite_log->SetUserLimits(maxStepLimit); 
    } 
  }

  else if (structure == 1){
    if (composite_pos == 0){
      if(case_geometry == 0){
        case_volume = new G4Box("case_volume", case_dx/2.0, case_dy/2.0, case_dz/2.0);
      }
      else if(case_geometry == 1){
        case_volume = new G4Sphere("case_volume", 0, case_dx/2.0, 0.0 * deg, 360.0 * deg, 0.0 * deg, 180.0 * deg);
      }
      case_log = new G4LogicalVolume(case_volume, CaseMaterial, "case_log");
      case_log -> SetVisAttributes(caseVis);   // set vis color
      if(multiple_scint == 0){// If there is an external volume that the case must be placed in, you must define the physical volume of the case after you have defined the external volume
        case_phys = new G4PVPlacement(0, G4ThreeVector(module_x,module_y,module_z), case_log, "case_phys", airGap_log, false, 0, fCheckOverlaps);
      }
      if(manual_step == 1){
        case_log->SetUserLimits(maxStepLimit); 
      }

      if (numCoords == 1){
        composite_box  = new G4Box("composite_box",composite_dx/2.0,composite_dy/2.0,composite_dz/2.0);
        composite_log  = new G4LogicalVolume(composite_box,CompositeMaterial,"composite_log");
        composite_phys = new G4PVPlacement(0,G4ThreeVector(),composite_log,"composite_phys",case_log,false,0,fCheckOverlaps);
        compositeVis -> SetForceSolid(true);               
        composite_log -> SetVisAttributes(compositeVis);   // set vis color
        if(manual_step == 1){
          composite_log->SetUserLimits(maxStepLimit);
        }
      }
      else {
          CompositeParameterisation* param = new CompositeParameterisation(case_dx, case_dy, case_dz, numCoords_x, numCoords_y, numCoords_z, composite_dx);

          const std::vector<G4ThreeVector>& positions = param->GetPositions();
          std::vector<std::string> composite_names(numCoords);

          composite_box = new G4Box("composite_box", composite_dx/2.0, composite_dx/2.0, composite_dx/2.0);
          composite_log = new G4LogicalVolume(composite_box, CompositeMaterial, "composite_log");
          compositeVis->SetForceSolid(true);
          composite_log->SetVisAttributes(compositeVis);
          if(manual_step == 1){
            composite_log->SetUserLimits(maxStepLimit);
          }
          
          auto start = std::chrono::high_resolution_clock::now();
          if(case_geometry==0){
            for (G4int copyNo = 0; copyNo < positions.size(); ++copyNo) {
              std::stringstream composite_stream;
              composite_stream << "composite_phys" << copyNo;
              composite_names[copyNo] = composite_stream.str();
              composite_phys = new G4PVPlacement(0, positions[copyNo], composite_log, composite_names[copyNo].c_str(), case_log, false, copyNo, fCheckOverlaps==false);
            }
            //crystal_phys = new G4PVParameterised(crystal_stream.str().c_str(), crystal_log, case_log, kUndefined, numCoords, param, fCheckOverlaps==false);
          }
          
          //G4cout << "Number of nanocrystals: " << crystal_number << G4endl;
          auto end = std::chrono::high_resolution_clock::now();
          std::chrono::duration<double> elapsed = end - start;
          G4cout << "Time taken to create G4PVParameterised: " << elapsed.count() << " seconds" << G4endl;
      }
    }

    if (composite_pos == 1){
      if(case_geometry == 0){
        case_volume = new G4Box("case_volume", case_dx/2.0, case_dy/2.0, case_dz/2.0);
      }
      else if(case_geometry == 1){
        case_volume = new G4Sphere("case_volume", 0, case_dx/2.0, 0.0 * deg, 360.0 * deg, 0.0 * deg, 180.0 * deg);
      }
      case_log = new G4LogicalVolume(case_volume, CaseMaterial, "case_log");
      case_log -> SetVisAttributes(caseVis);
      if(multiple_scint == 0){ 
        case_phys = new G4PVPlacement(0, G4ThreeVector(module_x,module_y,module_z), case_log, "case_phys", airGap_log, false, 0, fCheckOverlaps);
      }
      if(manual_step == 1){
        case_log->SetUserLimits(maxStepLimit); 
      }

      int copyNoFirstComposite = 0;

      std::stringstream composite_stream;
      composite_stream << "composite_phys" << copyNoFirstComposite;
      composite_box = new G4Box("composite_box", composite_dx / 2.0, composite_dx / 2.0, composite_dx / 2.0);
      composite_log = new G4LogicalVolume(composite_box, CompositeMaterial, "composite_log", 0, 0, 0);
      composite_phys = new G4PVPlacement(0, G4ThreeVector(), composite_log, composite_stream.str().c_str(), case_log, false, 0, fCheckOverlaps);
      //composite_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,case_dx/2-composite_dx*sqrt(3)/2), composite_log, "composite_phys", case_log, false, 0, fCheckOverlaps);
      compositeVis->SetForceSolid(true); 
      composite_log->SetVisAttributes(compositeVis);
      if(manual_step == 1){
        composite_log->SetUserLimits(maxStepLimit);
      }

      int i = 0;
      int copyNo = 1;
      while (i < numCoords-1) {
          G4double coord_x = G4UniformRand() * case_dx - case_dx / 2.0;
          G4double coord_y = G4UniformRand() * case_dx - case_dx / 2.0;
          G4double coord_z = G4UniformRand() * case_dx - case_dx / 2.0;

          G4ThreeVector currentPosition(coord_x, coord_y, coord_z);
          G4double distance_from_center = currentPosition.mag();

          G4VisAttributes* visAttributes;
          std::stringstream composite_stream;
          composite_stream << "composite_phys" << copyNo;

          composite_box = new G4Box("composite_box", composite_dx / 2.0, composite_dx / 2.0, composite_dx / 2.0);
          composite_log = new G4LogicalVolume(composite_box, CompositeMaterial, "composite_log", 0, 0, 0);
          compositeVis->SetForceSolid(true); 
          composite_log->SetVisAttributes(compositeVis);
          composite_phys = new G4PVPlacement(0, currentPosition, composite_log, composite_stream.str().c_str(), case_log, false, 0, fCheckOverlaps==false);

          if(manual_step == 1){
            composite_log->SetUserLimits(maxStepLimit);
          }

          bool checkOverlap = composite_phys->CheckOverlaps();

          if (checkOverlap) {
              delete composite_phys;
              delete composite_log;
              continue;
          } else {
              ++i;
              ++copyNo;
              G4cout << "Number of nanocomposites: " << i+1 << G4endl;
          }
          //composite_log->SetVisAttributes(compositeVis); 
      }
    }
  }

  if(multiple_scint == 1){

    sample_volume = new G4Box("sample_volume", sample_dx/2.0, sample_dx/2.0, sample_dx/2.0);
    sample_log = new G4LogicalVolume(sample_volume, SampleMaterial, "sample_log");
    sample_phys = new G4PVPlacement(0, G4ThreeVector(module_x,module_y,module_z), sample_log, "sample_phys", airGap_log, false, 0, fCheckOverlaps);
    sample_log->SetVisAttributes(sampleVis);
    // sampleVis->SetForceSolid(false);
    if(manual_step == 1){
      sample_log->SetUserLimits(maxStepLimit); 
    }

    int copyNoFirstCase = 0;

    std::stringstream case_stream;
    case_stream << "case_phys" << copyNoFirstCase;
    case_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), case_log, case_stream.str().c_str(), sample_log, false, 0, fCheckOverlaps);

    int i = 0;
    int copyNoCase = 1;
    while (i < numCoordsCase-1) {
      G4double coord_x = G4UniformRand() * sample_dx - sample_dx / 2.0;
      G4double coord_y = G4UniformRand() * sample_dx - sample_dx / 2.0;
      G4double coord_z = G4UniformRand() * sample_dx - sample_dx / 2.0;

      std::stringstream case_stream;
      case_stream << "case_phys" << copyNoCase;

      G4ThreeVector currentPosition(coord_x, coord_y, coord_z);
      case_phys = new G4PVPlacement(0, currentPosition, case_log, case_stream.str().c_str(), sample_log, false, 0, fCheckOverlaps==false);

      bool checkOverlap = case_phys->CheckOverlaps();

      if (checkOverlap) {
          delete case_phys;
          continue;
      } else {
          ++i;
          ++copyNoCase;
          G4cout << "Number of cases: " << i+1 << G4endl;
      }
    }

  }


  if(capsula==1){
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* CapsulaMaterial = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
    G4Box*             capsula_box  = new G4Box("Capsula",capsula_dx/2.0,capsula_dy/2.0, capsula_dz/2.0);
    G4LogicalVolume*   capsula_log  = new G4LogicalVolume(capsula_box,CapsulaMaterial,"Capsula",0,0,0);
    G4VPhysicalVolume* capsula_phys = new G4PVPlacement(0,G4ThreeVector(capsulax,capsulay,capsulaz),capsula_log,"Capsula",expHall_log,false,0,fCheckOverlaps);
    // capsulaVis->SetForceSolid(true);
    capsula_log->SetVisAttributes(capsulaVis);   // set vis color
    if(manual_step == 1){ 
      capsula_log->SetUserLimits(maxStepLimit); 
    }
  }



  for(G4int iSide = 0; iSide < number_of_detectors; iSide++)
  {
    if(opticalCoupling_thickness > 0.)
    {
      opticalCoupling_box  = new G4Box("opticalCoupling_box",opticalCoupling_dx[iSide]/2.0,opticalCoupling_dy[iSide]/2.0,opticalCoupling_dz[iSide]/2.0);
      opticalCoupling_log  = new G4LogicalVolume(opticalCoupling_box,OpticalCouplingMaterial,"opticalCoupling_log");
      opticalCoupling_phys = new G4PVPlacement(0,G4ThreeVector(opticalCoupling_x[iSide],opticalCoupling_y[iSide],opticalCoupling_z[iSide]),opticalCoupling_log,"opticalCoupling_phys",airGap_log,false,iSide,fCheckOverlaps);

      // opCompositeToGapSurface[iSide] = new G4OpticalSurface("CompositeToGapSurface");
      // opCompositeToGapSurface[iSide]->SetType(dielectric_dielectric);
      // opCompositeToGapSurface[iSide]->SetFinish(polished);
      // opCompositeToGapSurface[iSide]->SetModel(unified);
      // new G4LogicalBorderSurface("CompositeToWorldSurface",composite_phys,opticalCoupling_phys[iSide],opCompositeToGapSurface[iSide]);

    }

    switch(number_of_detectors)
    {
      case 1:
        pmtWindow_box = new G4Box("pmtWindow_box",photoDetSize/2.0,pmtWindowThickness/2.0,photoDetSize/2.0);
        photoCath_box = new G4Box("photoCath_box",photoDetSize/2.0,photoCathThickness/2.0,photoDetSize/2.0);
        break;
      case 2:
        pmtWindow_box = new G4Box("pmtWindow_box",photoDetSize/2.0,pmtWindowThickness/2.0,photoDetSize/2.0);
        photoCath_box = new G4Box("photoCath_box",photoDetSize/2.0,photoCathThickness/2.0,photoDetSize/2.0);
        break;
      case 6:
        if ((iSide == 0) || (iSide == 2))
        {
          pmtWindow_box = new G4Box("pmtWindow_box",photoDetSize/2.0,pmtWindowThickness/2.0,photoDetSize/2.0);
          photoCath_box = new G4Box("photoCath_box",photoDetSize/2.0,photoCathThickness/2.0,photoDetSize/2.0);
        } else if ((iSide == 1) || (iSide == 3))
        {
          pmtWindow_box = new G4Box("pmtWindow_box",pmtWindowThickness/2.0,photoDetSize/2.0,photoDetSize/2.0);
          photoCath_box = new G4Box("photoCath_box",photoCathThickness/2.0,photoDetSize/2.0,photoDetSize/2.0);
        } else if ((iSide == 4) || (iSide == 5))
        {
          pmtWindow_box = new G4Box("pmtWindow_box",photoDetSize/2.0,photoDetSize/2.0,pmtWindowThickness/2.0);
          photoCath_box = new G4Box("photoCath_box",photoDetSize/2.0,photoDetSize/2.0,photoCathThickness/2.0);
        }

    }

    opticalCoupling_log->SetVisAttributes(opticalCouplingVis);

    pmtWindow_log = new G4LogicalVolume(pmtWindow_box,PmtWindowMaterial,"pmtWindow_log");
    pmtWindow_phys = new G4PVPlacement(0,G4ThreeVector(pmtWindow_x[iSide],pmtWindow_y[iSide],pmtWindow_z[iSide]),pmtWindow_log,"pmtWindow_phys",airGap_log,false,iSide,fCheckOverlaps);
    pmtWindow_log->SetVisAttributes(pmtWindowVis);

    photoCath_log = new G4LogicalVolume(photoCath_box,PhotoCathodeMaterial,"photoCath_log");
    photoCath_phys = new G4PVPlacement(0,G4ThreeVector(photoCath_x[iSide],photoCath_y[iSide],photoCath_z[iSide]),photoCath_log,"photoCath_phys",pmtWindow_log,false,iSide,fCheckOverlaps);
    photoCath_log->SetVisAttributes(photoCathVis);

       // set vis color
    // pmtWindow_log[1]->SetVisAttributes(pmtWindowVis);   // set vis color
    // std::stringstream Surfname;
    // Surfname << "Composite_to_gap_surface" << iLayer << "_" << iSide;

    opReflectorSurface = new G4OpticalSurface("reflectorSurface");
    opReflectorSurface->SetType(dielectric_metal);
    opReflectorSurface->SetFinish(polished);
    opReflectorSurface->SetModel(glisur);
    opReflectorSurface->SetMaterialPropertiesTable(ESR_reflector);
    new G4LogicalBorderSurface("reflectorSurface",airGap_phys,reflector_phys,opReflectorSurface);


    std::vector<G4double> ephoton_airPmt = { 2.034 * eV, 4.136 * eV };
    std::vector<G4double> reflectivity_airPmt = { 0.0, 0.0 };
    std::vector<G4double> efficiency_airPmt   = { 1.0, 1.0 };

    opAirPmtSurface = new G4OpticalSurface("airPmtSurface");
    opAirPmtSurface->SetType(dielectric_metal);
    opAirPmtSurface->SetFinish(polished);
    opAirPmtSurface->SetModel(glisur);
    auto opAirPmtSurfaceMPT = new G4MaterialPropertiesTable();
    opAirPmtSurfaceMPT->AddProperty("REFLECTIVITY", ephoton_airPmt, reflectivity_airPmt);
    opAirPmtSurfaceMPT->AddProperty("EFFICIENCY", ephoton_airPmt, efficiency_airPmt);
    new G4LogicalBorderSurface("airPmtSurface",airGap_phys,pmtWindow_phys,opAirPmtSurface);
    new G4LogicalBorderSurface("airPmtSurface",airGap_phys,photoCath_phys,opAirPmtSurface);

    // opCompositeToPhotDetSurface = new G4OpticalSurface("opCompositeToPhotDetSurface");
    // opCompositeToPhotDetSurface->SetType(dielectric_dielectric);
    // opCompositeToPhotDetSurface->SetFinish(polished);
    // opCompositeToPhotDetSurface->SetModel(unified);
    // new G4LogicalBorderSurface("CompositeToPhotDetSurface",composite_phys,pmtWindow_phys[iSide],opCompositeToPhotDetSurface);

    std::vector<G4double> ephoton = { 1.4 * eV, 5 * eV };
    std::vector<G4double> photoCath_EFF     = { 1., 1. };
    std::vector<G4double> photoCath_ReR     = { 1.92, 1.92 };
    std::vector<G4double> photoCath_ImR     = { 1.69, 1.69 };
    auto photoCath_mt = new G4MaterialPropertiesTable();
    photoCath_mt->AddProperty("EFFICIENCY", ephoton, photoCath_EFF);
    photoCath_mt->AddProperty("REALRINDEX", ephoton, photoCath_ReR);
    photoCath_mt->AddProperty("IMAGINARYRINDEX", ephoton, photoCath_ImR);
    auto photoCath_opsurf = new G4OpticalSurface(
      "photoCath_opsurf", glisur, polished, dielectric_metal);
    photoCath_opsurf->SetMaterialPropertiesTable(photoCath_mt);

    //**Create logical skin surfaces
    new G4LogicalSkinSurface("photoCath_surf", photoCath_log, photoCath_opsurf);

    // save detector volumes in output file
    CreateTree::Instance()->Detector_boxSide = photoDetSize;
    CreateTree::Instance()->detectorID = detector_id;
    CreateTree::Instance()->detectorX  = pmtWindow_x[iSide];
    CreateTree::Instance()->detectorY  = pmtWindow_y[iSide];
    CreateTree::Instance()->detectorZ  = pmtWindow_z[iSide];
    CreateTree::Instance()->detectorDX = photoDetSize;
    CreateTree::Instance()->detectorDY = photoDetSize;
    CreateTree::Instance()->detectorDZ = photoDetSize;
    CreateTree::Instance()->detectors->Fill();


    detector_id++;

    }


  // finally, add all used materials to the material list struct in output

  for(G4int iMat = 0; iMat <  materialsList.size(); iMat++)
  {
    // material_t tempMaterial;
    // tempMaterial.number = iMat;
    // tempMaterial.name = materialsList[iMat];

    CreateTree::Instance()->materialsList.push_back(materialsList[iMat]);
  }

  std::cout << ">>>>>>>>>> Geometry built." << std::endl;
  // always return the world
  return expHall_phys;



}

void DetectorConstruction::ConstructSDandField()
{
  PMTSD* pmt = pmtDet.Get();
  if(!pmt)
  {
    // Created here so it exists as pmts are being placed
    G4cout << "Construction /Det/pmtSD" << G4endl;
    auto pmt_SD = new PMTSD("/Det/pmtSD");
    pmtDet.Put(pmt_SD);
  }
  G4SDManager::GetSDMpointer()->AddNewDetector(pmtDet.Get());
  SetSensitiveDetector(photoCath_log, pmtDet.Get());

  ScintSD* scint = scintDet.Get();
  if(!scint)
  {
    G4cout << "Construction /Det/scintSD" << G4endl;
    auto scint_SD = new ScintSD("/Det/scintSD");
    scintDet.Put(scint_SD);
  }
  G4SDManager::GetSDMpointer()->AddNewDetector(scintDet.Get());
  SetSensitiveDetector(composite_log, scintDet.Get());


}
