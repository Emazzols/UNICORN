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

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "globals.hh"
#include "ConfigFile.h"
#include "G4Geantino.hh"
#include "G4IonTable.hh"

class G4ParticleGun;
class G4Event;
class PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// class to define 3d vectors and the operations between them
class CVec3
{
public:
  // Data
  float x, y, z;
  // Ctors
  CVec3( float InX, float InY, float InZ ) : x( InX ), y( InY ), z( InZ )
  {
  }
  CVec3( ) : x(0), y(0), z(0)
  {
  }
  // Operator Overloads
  inline bool operator== (const CVec3& V2) const
  {
    return (x == V2.x && y == V2.y && z == V2.z);
  }

  inline CVec3 operator+ (const CVec3& V2) const
  {
    return CVec3( x + V2.x,  y + V2.y,  z + V2.z);
  }
  inline CVec3 operator- (const CVec3& V2) const
  {
    return CVec3( x - V2.x,  y - V2.y,  z - V2.z);
  }
  inline CVec3 operator- ( ) const
  {
    return CVec3(-x, -y, -z);
  }

  inline CVec3 operator/ (float S ) const
  {
    float fInv = 1.0f / S;
    return CVec3 (x * fInv , y * fInv, z * fInv);
  }
  inline CVec3 operator/ (const CVec3& V2) const
  {
    return CVec3 (x / V2.x,  y / V2.y,  z / V2.z);
  }
  inline CVec3 operator* (const CVec3& V2) const
  {
    return CVec3 (x * V2.x,  y * V2.y,  z * V2.z);
  }
  inline CVec3 operator* (float S) const
  {
    return CVec3 (x * S,  y * S,  z * S);
  }

  inline void operator+= ( const CVec3& V2 )
  {
    x += V2.x;
    y += V2.y;
    z += V2.z;
  }
  inline void operator-= ( const CVec3& V2 )
  {
    x -= V2.x;
    y -= V2.y;
    z -= V2.z;
  }

  inline float operator[] ( int i )
  {
    if ( i == 0 ) return x;
    else if ( i == 1 ) return y;
    else return z;
  }

  // Functions
  inline float Dot( const CVec3 &V1 ) const
  {
    return V1.x*x + V1.y*y + V1.z*z;
  }

  inline CVec3 CrossProduct( const CVec3 &V2 ) const
  {
    return CVec3(
      y * V2.z  -  z * V2.y,
      z * V2.x  -  x * V2.z,
      x * V2.y  -  y * V2.x   );
  }

  // Return vector rotated by the 3x3 portion of matrix m
  // (provided because it's used by bbox.cpp in article 21)
  CVec3 RotByMatrix( const float m[16] ) const
  {
    return CVec3(
    x*m[0] + y*m[4] + z*m[8],
    x*m[1] + y*m[5] + z*m[9],
    x*m[2] + y*m[6] + z*m[10] );
  }

  // These require math.h for the sqrtf function
  float Magnitude( ) const
  {
    return sqrtf( x*x + y*y + z*z );
  }

  float Distance( const CVec3 &V1 ) const
  {
    return ( *this - V1 ).Magnitude();
  }

  inline void Normalize()
  {
    float fMag = ( x*x + y*y + z*z );
    if (fMag == 0) {return;}

    float fMult = 1.0f/sqrtf(fMag);
    x *= fMult;
    y *= fMult;
    z *= fMult;
    return;
  }
};
// end of class CVect3




class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(ConfigFile& config);
  virtual ~PrimaryGeneratorAction();

public:
  virtual void GeneratePrimaries(G4Event*);
  int GetIntersection( float fDst1, float fDst2, CVec3 P1, CVec3 P2, CVec3 &Hit);
  int InBox( CVec3 Hit, CVec3 B1, CVec3 B2, const int Axis);
  int CheckLineBox( CVec3 B1, CVec3 B2, CVec3 L1, CVec3 L2, CVec3 &Hit);

  //     void SetOptPhotonPolar();
  //     void SetOptPhotonPolar(G4double);

private:
  G4ParticleGun* fParticleGun;
  // PrimaryGeneratorMessenger* fGunMessenger;
  G4GeneralParticleSource* fGPS; 
  G4int source_from_file;
  G4String energyFileName;
  bool sourceT5;
  G4double DBDrandom;
  G4double sourcex;
  G4double sourcey;
  G4double sourcez;

  G4double composite_center_x;
  G4double composite_center_y;
  G4double composite_center_z;

  G4double opticalCoupling_thickness;
  G4double pmtWindowThickness;

  G4int case_geometry; 

  G4double case_dx; 
  G4double case_dy; 
  G4double case_dz; 

  G4double sample_dx; 
  G4int multiple_scint; 

  G4double source_dim; 

  G4double radius; 
  G4double r; 
  G4double theta_ang; 
  G4double phi_ang;

  G4double direction_x;
  G4double direction_y;
  G4double direction_z;

  G4double crystal_dim_x;
  G4double crystal_dim_y;
  G4double crystal_dim_z;

  int layers;
  G4double gap_between_layers       ;
  G4double air_gap_cry_esr_thickness;
  G4double esr_thickness            ;
  G4double layer_dim_z              ;

  // G4double *layer_z;
  // G4double first_layer_z ;
  // G4double first_entrance_z;

  // G4double crystalz;
  // G4double EventSourcex;
  // G4double EventSourcey;
  // G4double EventSourcez;
  // G4double EventTargetx;
  // G4double EventTargety;
  // G4double EventTargetz;
  // G4double greaseBack;
  // G4double glassBack;
  // G4double airBack;

  // G4double distance;
  G4double energy;
  G4double position;
  G4double ExcitEnergy; 
  G4double IonCharge; 
  G4int Z; 
  G4int A; 
  G4int full_chain;
  G4int particle_type;
  // G4double direction;
  // G4double userIrradiationDiagonal;

  G4double module_dx;
  G4double module_dy;
  G4double module_dz;

  G4double x_d;
  G4double x_u;
  G4double y_d;
  G4double y_u;
  G4double z_d;
  G4double z_u;

  // G4double sphereRadius;

  // G4double theta;
  // G4double phi ;
  // G4double esrThickness;
  // G4double crystalx;
  // G4double crystaly;
  // G4int ncrystalx;
  // G4int ncrystaly;
  // G4bool backgroudSimulation;






};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*PrimaryGeneratorAction_h*/
