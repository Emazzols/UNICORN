#include "MySurfaces.hh"
#include "G4NistManager.hh"

using namespace CLHEP;

MySurfaces::MySurfaces()
{
  // default constructor
}

MySurfaces::~MySurfaces()
{}


G4MaterialPropertiesTable * MySurfaces::ESR(G4double esrTransmittance)
{
  //ESR
  std::vector<G4double> vikuiti_energy_im = { 1.25*eV, 1.31*eV, 1.34*eV, 1.38*eV, 1.42*eV, 1.46*eV, 1.51*eV, 1.60*eV, 1.66*eV, 1.75*eV, 1.80*eV, 1.90*eV, 1.95*eV, 2.05*eV, 2.10*eV, 2.20*eV, 2.26*eV, 2.35*eV, 2.45*eV, 2.51*eV, 2.71*eV, 2.81*eV, 3.01*eV, 3.21*eV, 3.41*eV, 3.55*eV, 3.71*eV, 3.91*eV, 4.14*eV };
  //   G4double vikuiti_RIndex_im[NUMvikuiti] =
  //   {
  //     9.49/1.47, 8.88/1.47, 8.49/1.47, 8.30/1.47, 8.18/1.47, 8.22/1.47, 8.31/1.47,
  //     8.60/1.47, 8.62/1.47, 8.39/1.47, 8.21/1.47, 7.82/1.47, 7.65/1.47, 7.31/1.47,
  //     7.15/1.47, 6.85/1.47, 6.69/1.47, 6.42/1.47, 6.15/1.47, 6.03/1.47, 5.58/1.47,
  //     5.38/1.47, 5.02/1.47, 4.71/1.47, 4.43/1.47, 4.24/1.47, 4.06/1.47, 3.84/1.47, 3.61/1.47
  //
  //   };
  // we specify the imaginary index WITHOUT dividing by 1.47 because in this matrix
  // esr is in air contact with everything, so n1=1
  std::vector<G4double> vikuiti_RIndex_im = { 9.49, 8.88, 8.49, 8.30, 8.18, 8.22, 8.31, 8.60, 8.62, 8.39, 8.21, 7.82, 7.65, 7.31, 7.15, 6.85, 6.69, 6.42, 6.15, 6.03, 5.58, 5.38, 5.02, 4.71, 4.43, 4.24, 4.06, 3.84, 3.61 };
  std::vector<G4double> vikuiti_energy_re = { 1.24*eV, 1.27*eV, 1.31*eV, 1.34*eV, 1.38*eV, 1.42*eV, 1.46*eV, 1.51*eV, 1.55*eV, 1.60*eV, 1.66*eV, 1.71*eV, 1.77*eV, 1.84*eV, 1.91*eV, 1.99*eV, 2.07*eV, 2.16*eV, 2.26*eV, 2.37*eV, 2.48*eV, 2.62*eV, 2.76*eV, 2.92*eV, 3.11*eV, 3.31*eV, 3.55*eV, 3.82*eV, 4.14*eV };

  std::vector<G4double> vikuiti_RIndex_re = { 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.1, 0.6, 1.2, 0.4, 0.23, 0.17, 0.26, 1.3, 0.24, 0.11, 0.19, 0.24, 0.25, 0.87, 0.595, 0.18, 0.178, 0.6, 8.0, 8.0, 8.0, 8.0 };

  std::vector<G4double> pp = {2.038*eV, 4.144*eV};

  G4MaterialPropertiesTable *ESR_surf = new G4MaterialPropertiesTable();
  std::vector<G4double> FixedTransmittance_surf = {esrTransmittance, esrTransmittance};

  std::vector<G4double> esrTransmittance_energy = { 2.7008152174*eV, 2.7304945055*eV, 2.7608333333*eV, 2.7918539326*eV,	2.8235795455*eV, 2.8560344828*eV, 2.889244186*eV , 2.9232352941*eV,2.9580357143*eV, 2.9936746988*eV, 3.0301829268*eV, 3.0675925926*eV, 3.1059375*eV   , 3.1452531646*eV, 3.1855769231*eV, 3.2269480519*eV, 3.2694078947*eV};
  std::vector<G4double> esrTransmittance_values = {0.0024,0.0025,0.002 ,0.0024,0.0015,0.0011,0.0011,0.0011,0.0009,0.0007,0.001 ,0.0016,0.005 ,0.0164,0.0303,0.0365,0.0389};


  //modfication for double surface. since the Transmittance measured is for a real material, with 2 surfaces, we need to scale the measurement to have a correct result when we use the two of them
  //this is anyway a big question mark..
  // for(int esr_n = 0; esr_n < numEsrTrans ; esr_n++)
  // {
  //   esrTransmittance_values[esr_n] = pow(esrTransmittance_values[esr_n],0.6666666666666666666);
  // }

  ESR_surf->AddProperty ("REALRINDEX",        vikuiti_energy_re,  vikuiti_RIndex_re);
  ESR_surf->AddProperty ("IMAGINARYRINDEX",   vikuiti_energy_im,  vikuiti_RIndex_im);
  if(esrTransmittance == -1)
  {
    ESR_surf->AddProperty ("TRANSMITTANCE", esrTransmittance_energy , esrTransmittance_values);
  }
  else
  {
    ESR_surf->AddProperty ("TRANSMITTANCE", pp , FixedTransmittance_surf);
  }

  return ESR_surf;

}

G4MaterialPropertiesTable * MySurfaces::Simplified_ESR(G4double simplified_esr_reflectivity,G4double simplified_esr_efficiency)
{
  G4MaterialPropertiesTable *ESR_surf = new G4MaterialPropertiesTable();

  std::vector<G4double> pp = {0.5*eV, 8*eV};

  // std::vector<G4double> FixedTransmittance = {simplified_esr_transmittance, simplified_esr_transmittance};
  std::vector<G4double> FixedReflectivity  = {simplified_esr_reflectivity , simplified_esr_reflectivity };
  std::vector<G4double> FixedEfficiency  = {simplified_esr_efficiency , simplified_esr_efficiency };


  // ESR_surf->AddProperty ("TRANSMITTANCE", pp , FixedTransmittance);
  ESR_surf->AddProperty ("REFLECTIVITY", pp , FixedReflectivity);
  ESR_surf->AddProperty("EFFICIENCY", pp, FixedEfficiency);
  return ESR_surf;
}


G4MaterialPropertiesTable * MySurfaces::DepolishedSurface(G4double specularlobe,
                                                          G4double specularspike,
                                                          G4double backscattering
                                                        )
{
  G4MaterialPropertiesTable *crystalDepolished_surf = new G4MaterialPropertiesTable();
  std::vector<G4double> depoEnergy = {0.5*eV, 8*eV};
  std::vector<G4double> specularlobeVector = {specularlobe, specularlobe};
  std::vector<G4double> specularspikeVector = {specularspike, specularspike};
  std::vector<G4double> backscatterVector = {backscattering,backscattering};
  crystalDepolished_surf->AddProperty("SPECULARLOBECONSTANT",depoEnergy,specularlobeVector);
  crystalDepolished_surf->AddProperty("SPECULARSPIKECONSTANT",depoEnergy,specularspikeVector);
  crystalDepolished_surf->AddProperty("BACKSCATTERCONSTANT",depoEnergy,backscatterVector);
  return crystalDepolished_surf;
}


G4MaterialPropertiesTable * MySurfaces::TeflonSurface(G4double teflonReflectivity)
{
  G4MaterialPropertiesTable *teflonSurfPropertiesTable = new G4MaterialPropertiesTable();
  std::vector<G4double> teflonEnergy = {0.5*eV, 8*eV};
  std::vector<G4double> teflonlobeVector = {0., 0.};
  std::vector<G4double> teflonSpecularspikeVector = {0., 0.};
  std::vector<G4double> teflonBackscatterVector = {0.,0.};
  std::vector<G4double> teflonReflectivityVector = {teflonReflectivity,teflonReflectivity};
  teflonSurfPropertiesTable->AddProperty("SPECULARLOBECONSTANT", teflonEnergy,teflonlobeVector);
  teflonSurfPropertiesTable->AddProperty("SPECULARSPIKECONSTANT", teflonEnergy,teflonSpecularspikeVector);
  teflonSurfPropertiesTable->AddProperty("BACKSCATTERCONSTANT", teflonEnergy,teflonBackscatterVector);
  teflonSurfPropertiesTable->AddProperty("REFLECTIVITY", teflonEnergy, teflonReflectivityVector);

  return teflonSurfPropertiesTable;
}





// //------------------------------------//
// //            DEPOLISHING             //
// //------------------------------------//
//
// // -------- Lateral faces ------------//
//
// //when crystal is depolished, we can define a surface between crystal and air thin layer
// G4MaterialPropertiesTable *crystalDepolished_surf = new G4MaterialPropertiesTable();
// const G4int Ndepo = 2;
// G4double depoEnergy[Ndepo] = {2.038*eV, 4.144*eV};
// G4double specularlobeVector[Ndepo] = {latspecularlobe, latspecularlobe};
// G4double specularspikeVector[Ndepo] = {latspecularspike, latspecularspike};
// G4double backscatterVector[Ndepo] = {latbackscattering,latbackscattering};
// crystalDepolished_surf->AddProperty("SPECULARLOBECONSTANT",depoEnergy,specularlobeVector,Ndepo);
// crystalDepolished_surf->AddProperty("SPECULARSPIKECONSTANT",depoEnergy,specularspikeVector,Ndepo);
// crystalDepolished_surf->AddProperty("BACKSCATTERCONSTANT",depoEnergy,backscatterVector,Ndepo);
// if(latsurfaceroughness != 0)
//   crystalDepolished_surf->AddConstProperty("SURFACEROUGHNESS", latsurfaceroughness);
//
// // ---- Front and Back faces --------//
//
// G4MaterialPropertiesTable *crystalReal_surf = new G4MaterialPropertiesTable();
// G4double specularlobeRealVector[Ndepo] = {realspecularlobe, realspecularlobe};
// G4double specularspikeRealVector[Ndepo] = {realspecularspike, realspecularspike};
// G4double backscatterRealVector[Ndepo] = {realbackscattering,realbackscattering};
// crystalReal_surf->AddProperty("SPECULARLOBECONSTANT",depoEnergy,specularlobeRealVector,Ndepo);
// crystalReal_surf->AddProperty("SPECULARSPIKECONSTANT",depoEnergy,specularspikeRealVector,Ndepo);
// crystalReal_surf->AddProperty("BACKSCATTERCONSTANT",depoEnergy,backscatterRealVector,Ndepo);
// if(realsurfaceroughness != 0)
//   crystalReal_surf->AddConstProperty("SURFACEROUGHNESS", realsurfaceroughness);
