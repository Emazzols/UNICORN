#ifndef MySurfaces_hh
#define MySurfaces_hh

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4MaterialPropertiesTable.hh"
// #include "globals.hh"
// #include "ConfigFile.h"


class MySurfaces
{
private:

public:
  MySurfaces();
  ~MySurfaces();
  G4MaterialPropertiesTable* ESR(G4double esrTransmittance);
  G4MaterialPropertiesTable * DepolishedSurface(G4double specularlobe,G4double specularspike,G4double backscattering);
  G4MaterialPropertiesTable * TeflonSurface(G4double teflonReflectivity);
  G4MaterialPropertiesTable * Simplified_ESR(G4double simplified_esr_transmittance,G4double simplified_esr_reflectivity);
};

#endif
