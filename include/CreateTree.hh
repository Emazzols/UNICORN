#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "struct.hh"
#include "globals.hh"

class CreateTree
{
private:

  TTree*              ftree;
  TString             fname;

  Bool_t              HITS;
  Bool_t              ABSORPTIONS;

  static const Int_t  MaxNum = 2000000;
  static const Int_t  MaxNumPro = 100;

  int nCrystalsX, nCrystalsY, nDetectorsX, nDetectorsY;

public:

  CreateTree(TString name,int x, int y, int z, int k);
  ~CreateTree();


  TTree*              GetTree() const { return ftree; };
  TString             GetName() const { return fname;};
  Int_t               Fill() { return this->GetTree()->Fill(); };
  Bool_t              Write();
  void                Clear();
  static CreateTree*  Instance() { return fInstance; };
  static CreateTree*  fInstance;

  Bool_t              Hits() const { return this->HITS; };
  Bool_t              Absorptions() const { return this->ABSORPTIONS; };

  //Short_t             tempOriginCrystal;

  double module_dx;
  double module_dy;
  double module_dz;

  // Global part
  long int            Seed;
  Int_t               Run;
  Int_t               Event;
  Float_t             SourceX;
  Float_t             SourceY;
  Float_t             SourceZ;
  Float_t             SourceMomentumX;
  Float_t             SourceMomentumY;
  Float_t             SourceMomentumZ;
  Float_t             totalEnergyDeposited;
  Float_t             totalEnergyDepositedInExtraVolume;
  Float_t             totalEnergyDepositedInSample;
  Float_t             totalEnergyDepositedInCase;
  Float_t             totalEnergyDepositedInFirstCrystal;
  Float_t             totalEnergyDepositedInCrystals;
  Int_t               NumOptPhotons;
  Int_t               NumCherenkovPhotons;
  std::vector<Float_t> Edep_vec;
  std::vector<Float_t> X_vec;
  std::vector<Float_t> Y_vec;
  std::vector<Float_t> Z_vec;
  std::vector<Float_t> StepLength_vec;
  std::vector<Float_t> TrackLength_vec;
  int verbose;

  //Hits on detectors
  // Short_t*            DetectorHit;

  // Optical Photons part

  //energy deposition part


  // materials list part
  // std::vector<material_t> materials;
  std::vector<std::string> materialsList;
  // bool exitFound;
  // int exitFace;
  // bool FoundGammaFirstEntrance;
  // bool FoundGammaFirstDeposition;
  bool FoundGammaVertexMomentum;
  bool primaryFirstEntranceFound ;
  Float_t primaryPositionOnFirstCrystalX;
  Float_t primaryPositionOnFirstCrystalY;
  Float_t primaryPositionOnFirstCrystalZ;
  Float_t primaryMomentumOnFirstCrystalX;
  Float_t primaryMomentumOnFirstCrystalY;
  Float_t primaryMomentumOnFirstCrystalZ;
  Float_t primaryEnergyOnFirstCrystal   ;
  Int_t   primaryInteracted   ;

  TTree*  detectors;
  Float_t Detector_boxSide;
  Int_t   chPerSide;
  Int_t   detectorID;
  Float_t detectorX ;
  Float_t detectorY ;
  Float_t detectorZ ;
  Float_t detectorDX;
  Float_t detectorDY;
  Float_t detectorDZ;

  TTree* photons;
  TTree* photonsTotal;
  //Double_t pde;
  //Int_t primaries_per_event;
  //Int_t events_per_run;
  //Float_t vertX         ;
  //Float_t vertY         ;
  //Float_t vertZ         ;
  //Float_t vertMomentumX  ;
  //Float_t vertMomentumY  ;
  //Float_t vertMomentumZ  ;
  Float_t PositionX     ;
  Float_t PositionY     ;
  Float_t PositionZ     ;
  Float_t PreMomentumX  ;
  Float_t PreMomentumY  ;
  Float_t PreMomentumZ  ;
  //Float_t PostMomentumX ;// No need to do it, photons are killed -> PostMomentum is automatically 0
  //Float_t PostMomentumY ;
  //Float_t PostMomentumZ ;
  Int_t NumOptPhotonsReachedPMT ;
  Int_t TotalNumOptPhotonsReachedPMT ;
  Float_t globalTime    ;
  //Float_t localTime    ;
  //Float_t PhotonEnergy  ;





  //exit face legend:
  // 0 - no exit was found among the next possibilities
  // 1 - exit from front face (face coupled to MPPC)
  // 2 - exit from back face (face opposite to MPPC)


  // std::vector<float>*  CryEnergyDeposited;
  // std::vector<float>** pCryEnergyDeposited;
  //
  // std::vector<float>*  PosXEnDep;
  // std::vector<float>** pPosXEnDep;
  // std::vector<float>*  PosYEnDep;
  // std::vector<float>** pPosYEnDep;
  // std::vector<float>*  PosZEnDep;
  // std::vector<float>** pPosZEnDep;

  // std::vector<float> PositionX;
  // std::vector<float> *pPositionX;
  // std::vector<float> PositionY;
  // std::vector<float> *pPositionY;
  // std::vector<float> PositionZ;
  // std::vector<float> *pPositionZ;
  //
  // std::vector<float> PreMomentumX;
  // std::vector<float> *pPreMomentumX;
  // std::vector<float> PreMomentumY;
  // std::vector<float> *pPreMomentumY;
  // std::vector<float> PreMomentumZ;
  // std::vector<float> *pPreMomentumZ;
  //
  // std::vector<float> PostMomentumX;
  // std::vector<float> *pPostMomentumX;
  // std::vector<float> PostMomentumY;
  // std::vector<float> *pPostMomentumY;
  // std::vector<float> PostMomentumZ;
  // std::vector<float> *pPostMomentumZ;
  //
  // std::vector<float> GlobalTime;
  // std::vector<float> *pGlobalTime;
  //
  // std::vector<int> PhotonType;
  // std::vector<int> *pPhotonType;
  //
  // std::vector<float> PhotonEnergy;
  // std::vector<float> *pPhotonEnergy;
  //
  // std::vector<int> OriginCrystalI;
  // std::vector<int> *pOriginCrystalI;
  // std::vector<int> OriginCrystalJ;
  // std::vector<int> *pOriginCrystalJ;


};
