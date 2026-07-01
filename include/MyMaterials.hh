#ifndef MyMaterials_hh
#define MyMaterials_hh

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4NistManager.hh"
// #include "struct.hh"
// #include "ConfigFile.h"


class MyMaterials
{
private:
  G4double a, z;
  G4Element *N  ;
  G4Element *H  ;
  G4Element *C  ;
  G4Element* O  ;
  G4Element *Lu ;
  G4Element *Si ;
  G4Element* Bi ;
  G4Element* Ge ;
  G4Element* B  ;
  G4Element* Na ;
  G4Element* K  ;
  G4Element* Al ;
  G4Element *Cd ;
  G4Element *Se ;
  G4Element *S  ;
  G4Element* F  ;
  G4Element* Te ;
  G4Element* Ga ;
  G4Element* In ;
  G4NistManager* fNistMan = nullptr;

public:
  MyMaterials();
  ~MyMaterials();


  G4Material* Vacuum();
  G4Material* Air();
  G4Material* AirThinLayer();
  G4Material* LYSO(G4double lightyield, G4double resolutionScale);
  G4Material* BGO(G4double lightyield, G4double resolutionScale);
  G4Material* OpticalGrease(G4double refr_index);
  G4Material* FusedSilica();
  G4Material* Silicium();
  G4Material* Glass();
  G4Material* Aluminum();
  // G4Material* SilicioMPPC();
  G4Material* BorosilicateGlass();
  G4Material* Epoxy(G4double refr_index);
  G4Material* Teflon();
  G4Material *BC422(G4double lightyield,G4double resolutionScale);
  G4Material *CdSeCdSinPMMA(G4double host_fractMass,
                            G4double NP_fractMass,
                            G4double lightyield,
                            G4double risetime,
                            G4double decaytime,
                            G4double rayleigh_length,
                            G4double mie_length);
  G4Material *CuDopedCdSNCs(G4double lightyield,
                            G4double risetime,
                            G4double decaytime,
                            G4double rayleigh_length,
                            G4double mie_length);
  G4Material* Plastic();
  G4Material* GalliumNitride();
  G4Material* InGaNitride();
  G4Material* AluminumOxide();
  G4Material* fPMMA = nullptr;
  std::vector<std::string> GetMaterialsList(){return materials;};

  std::vector<std::string> materials;


  void AddMaterialToList(std::string name);

};

#endif
