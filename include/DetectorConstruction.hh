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

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "MyMaterials.hh"
#include "MySurfaces.hh"
#include "G4VPVParameterisation.hh"
#include "ConfigFile.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PMTSD;
class ScintSD;

class CompositeParameterisation : public G4VPVParameterisation {
public:
    CompositeParameterisation(G4double case_dim_x, G4double case_dim_y, G4double case_dim_z, int composite_num_x, int composite_num_y, int composite_num_z, G4double composite_dim);
    virtual ~CompositeParameterisation() {}

    virtual void ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const override;
    virtual void ComputeDimensions(G4Box& box, const G4int copyNo, const G4VPhysicalVolume* physVol) const override;

    const std::vector<G4ThreeVector>& GetPositions() const;

private:
    G4double case_dim_x;
    G4double case_dim_y;
    G4double case_dim_z;
    std::vector<G4ThreeVector> positions;
    G4double composite_dim;
    G4double initial_offset_x;
    G4double initial_offset_y;
    G4double initial_offset_z;
    G4double spacing_x;
    G4double spacing_y;
    G4double spacing_z;
    G4double x;
    G4double y;
    G4double z;
};


class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction(ConfigFile& config,MyMaterials* Materials, MySurfaces* Surfaces);
  virtual ~DetectorConstruction();

public:
  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

  G4LogicalVolume* GetLogPhotoCath() { return photoCath_log; }
  G4LogicalVolume* GetLogScint() { return composite_log; }

  std::vector<G4ThreeVector> GetPmtPositions() { return fPmtPositions; }

private:

  ConfigFile& fConfig;
  MyMaterials* fMaterials;
  MySurfaces* fSurfaces;

  G4int manual_step; 
  G4double maxStep;
  G4UserLimits* maxStepLimit;

  G4bool fCheckOverlaps;
  G4double sourcex;
  G4double sourcey;
  G4double sourcez;

  G4double capsulax; 
  G4double capsulay; 
  G4double capsulaz;

  G4double PMTbuildStartPoint;
  G4double PMTgapBetwCenters;
  G4double PMTtotalGap;

  G4int number_of_detectors;
  G4double photoDetSize;
  G4double pmtWindowThickness;
  G4double photoCathThickness;

  G4double resolutionScale;

  // G4double fastrisetime ;
  // G4double fastdecaytime;
  // G4double fastratio    ;
  // G4double slowrisetime ;
  // G4double slowdecaytime;

  G4double host_fractMass;
  G4double NP_fractMass;
  G4double lightyield;
  G4double risetime;
  G4double decaytime;
  G4double rayleigh_length;
  G4double mie_length;


  G4int esr_type ;
  G4double simplified_esr_transmittance ;
  G4double simplified_esr_reflectivity  ;
  G4double simplified_esr_efficiency  ;
  G4double esrTransmittance;


  // MyMaterials* Materials;
  // MySurfaces* Surfaces;

  G4Material* WorldMaterial;
  G4Material* CompositeMaterial;
  G4Material* OpticalCouplingMaterial;
  G4Material* AirGapMaterial;
  G4Material* ReflectorMaterial;
  G4Material* PhotoCathodeMaterial;
  G4Material* PmtWindowMaterial;
  G4Material *CapsulaMaterial;
  G4Material *CaseMaterial;
  G4Material *SampleMaterial;


  G4int world_material;
  G4int composite_material;
  G4int opticalCoupling_material;
  G4int case_material; 
  G4int sample_material;

  std::vector<std::string> materialsList;

  G4MaterialPropertiesTable* ESR_reflector;

  G4double fExpHall_dx;
  G4double fExpHall_dy;
  G4double fExpHall_dz;

  G4double module_dx;
  G4double module_dy;
  G4double module_dz;

  G4double composite_dx;
  G4double composite_dy;
  G4double composite_dz;

  G4double case_dx; 
  G4double case_dy;
  G4double case_dz;

  G4double module_x;
  G4double module_y;
  G4double module_z;

  // G4double composite_x;
  // G4double composite_y;
  // G4double composite_z;

  G4double sample_dx; 

  G4int structure; 

  G4int multiple_scint; 

  G4int case_geometry; 

  G4int numCoords; 
  G4int numCoords_x; 
  G4int numCoords_y; 
  G4int numCoords_z; 
  G4int numCoordsCase; 

  G4double coord_x; 
  G4double coord_y; 
  G4double coord_z; 

  G4double composite_pos;

  G4int capsula; 
  G4double capsula_dx; 
  G4double capsula_dy; 
  G4double capsula_dz;

  G4double airGap_dx;
  G4double airGap_dy;
  G4double airGap_dz;

  G4double airGap_x;
  G4double airGap_y;
  G4double airGap_z;

  G4double reflector_dx;
  G4double reflector_dy;
  G4double reflector_dz;

  G4double reflector_x;
  G4double reflector_y;
  G4double reflector_z;

  std::vector<G4double> opticalCoupling_dx;
  std::vector<G4double> opticalCoupling_dy;
  std::vector<G4double> opticalCoupling_dz;

  std::vector<G4double> opticalCoupling_x;
  std::vector<G4double> opticalCoupling_y;
  std::vector<G4double> opticalCoupling_z;

  std::vector<G4double> pmtWindow_x;
  std::vector<G4double> pmtWindow_y;
  std::vector<G4double> pmtWindow_z;

  std::vector<G4double> photoCath_x;
  std::vector<G4double> photoCath_y;
  std::vector<G4double> photoCath_z;

  // G4int total_internal_reflection;


  G4double airGapThickness;
  G4double reflectorThickness;

  G4double opticalCoupling_thickness;

  G4double opticalCoupling_refr_index;
  // G4double gap_to_composite_refr_index;
  // G4double epoxy_refr_index;


  // G4int    depolished;
  // G4double sigmaalpha;
  // G4double specularlobe;
  // G4double specularspike;
  // G4double backscattering;

  G4VisAttributes* worldVis;
  G4VisAttributes* moduleVis;
  G4VisAttributes* compositeVis;
  G4VisAttributes* opticalCouplingVis;
  G4VisAttributes* pmtWindowVis;
  G4VisAttributes* photoCathVis;
  G4VisAttributes* airGapVis;
  G4VisAttributes* reflectorVis;
  G4VisAttributes* capsulaVis;
  G4VisAttributes* sampleVis;
  G4VisAttributes* caseVis;

  G4Box*             opticalCoupling_box;
  G4LogicalVolume*   opticalCoupling_log;
  G4VPhysicalVolume* opticalCoupling_phys;

  G4Box*             expHall_box;
  G4LogicalVolume*   expHall_log;
  G4VPhysicalVolume* expHall_phys;

  G4Box*             pmtWindow_box;
  G4LogicalVolume*   pmtWindow_log;
  G4VPhysicalVolume* pmtWindow_phys;

  G4Box*             photoCath_box;
  G4LogicalVolume*   photoCath_log;
  G4VPhysicalVolume* photoCath_phys;

  G4VSolid*          composite_box;
  G4LogicalVolume*   composite_log;
  G4VPhysicalVolume* composite_phys;

  G4VSolid*           case_volume ;
  G4LogicalVolume*    case_log ; 
  G4VPhysicalVolume*  case_phys ; 

  G4VSolid*           sample_volume ;
  G4LogicalVolume*    sample_log ; 
  G4VPhysicalVolume*  sample_phys ; 

  G4Box*             airGap_box;
  G4LogicalVolume*   airGap_log;
  G4VPhysicalVolume* airGap_phys;

  G4Box*             reflector_box;
  G4LogicalVolume*   reflector_log;
  G4VPhysicalVolume* reflector_phys;

  // G4OpticalSurface*  opCompositeToGapSurface;
  G4OpticalSurface*  opReflectorSurface;
  G4OpticalSurface*  opAirPmtSurface;
  // G4OpticalSurface*  opAirPhotoCathSurface;
  // G4OpticalSurface*  opCompositeToPhotDetSurface;
  // G4LogicalSkinSurface* compositeSurface;

  // G4MaterialPropertiesTable* DepolishedSurface;

  G4Cache<PMTSD*> pmtDet;
  G4Cache<ScintSD*> scintDet;

  // Sensitive Detectors positions
  std::vector<G4ThreeVector> fPmtPositions;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*DetectorConstruction_h*/
