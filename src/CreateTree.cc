#include "CreateTree.hh"
#include <sstream>
#include "TROOT.h"
// #include "struct.hh"

CreateTree* CreateTree::fInstance = NULL;

using namespace std;

CreateTree::CreateTree(TString name, int x, int y, int z, int k)
: nCrystalsX(x),
nCrystalsY(y),
nDetectorsX(z),
nDetectorsY(k)
{

  photons = new TTree("photons", "photons");
  //photons->Branch("run"           , &this->Run           , "run/I");
  //photons->Branch("event"         , &this->Event         , "event/I");
  //photons->Branch("vertX"         , &this->vertX         , "vertX/F");
  //photons->Branch("vertY"         , &this->vertY         , "vertY/F");
  //photons->Branch("vertZ"         , &this->vertZ         , "vertZ/F");
  photons->Branch("PreMomentumX"  , &this->PreMomentumX  , "PreMomentumX/F");
  photons->Branch("PreMomentumY"  , &this->PreMomentumY  , "PreMomentumY/F");
  photons->Branch("PreMomentumZ"  , &this->PreMomentumZ  , "PreMomentumZ/F");
  //photons->Branch("PostMomentumX" , &this->PostMomentumX , "PostMomentumX/F");// No need to do it, photons are killed -> PostMomentum is automatically 0
  //photons->Branch("PostMomentumY" , &this->PostMomentumY , "PostMomentumY/F");
  //photons->Branch("PostMomentumZ" , &this->PostMomentumZ , "PostMomentumZ/F");
  photons->Branch("PositionX"     , &this->PositionX     , "PositionX/F");
  photons->Branch("PositionY"     , &this->PositionY     , "PositionY/F");
  photons->Branch("PositionZ"     , &this->PositionZ     , "PositionZ/F");
  photons->Branch("globalTime"    , &this->globalTime    , "globalTime/F");
  photons->Branch("NumOptPhotonsReachedPMT",&this->NumOptPhotonsReachedPMT,"NumOptPhotonsReachedPMT/I");

  photonsTotal = new TTree("photonsTotal", "photonsTotal");
  photonsTotal->Branch("TotalNumOptPhotonsReachedPMT",&this->TotalNumOptPhotonsReachedPMT,"TotalNumOptPhotonsReachedPMT/I");

  gROOT->ProcessLine("#include <vector>"); //this is needed otherwise ROOT will complain about not knowing what a std::vector is...

  if(fInstance)
  {
    return;
  }

  this->HITS=true;
  this->ABSORPTIONS=true;

  this->fInstance = this;
  this->fname = name;
  this->ftree = new TTree("tree","name");


  this->GetTree()->Branch("Seed",&this->Seed,"Seed/L");
  this->GetTree()->Branch("Run",&this->Run,"Run/I");
  this->GetTree()->Branch("Event",&this->Event,"Event/I");
  this->GetTree()->Branch("SourceX",&this->SourceX,"SourceX/F");
  this->GetTree()->Branch("SourceY",&this->SourceY,"SourceY/F");
  this->GetTree()->Branch("SourceZ",&this->SourceZ,"SourceZ/F");
  this->GetTree()->Branch("SourceMomentumX",&this->SourceMomentumX,"SourceMomentumX/F");
  this->GetTree()->Branch("SourceMomentumY",&this->SourceMomentumY,"SourceMomentumY/F");
  this->GetTree()->Branch("SourceMomentumZ",&this->SourceMomentumZ,"SourceMomentumZ/F");
  this->GetTree()->Branch("primaryPositionOnFirstCrystalX",&this->primaryPositionOnFirstCrystalX,"primaryPositionOnFirstCrystalX/F");
  this->GetTree()->Branch("primaryPositionOnFirstCrystalY",&this->primaryPositionOnFirstCrystalY,"primaryPositionOnFirstCrystalY/F");
  this->GetTree()->Branch("primaryPositionOnFirstCrystalZ",&this->primaryPositionOnFirstCrystalZ,"primaryPositionOnFirstCrystalZ/F");
  this->GetTree()->Branch("primaryMomentumOnFirstCrystalX",&this->primaryMomentumOnFirstCrystalX,"primaryMomentumOnFirstCrystalX/F");
  this->GetTree()->Branch("primaryMomentumOnFirstCrystalY",&this->primaryMomentumOnFirstCrystalY,"primaryMomentumOnFirstCrystalY/F");
  this->GetTree()->Branch("primaryMomentumOnFirstCrystalZ",&this->primaryMomentumOnFirstCrystalZ,"primaryMomentumOnFirstCrystalZ/F");
  this->GetTree()->Branch("primaryEnergyOnFirstCrystal",&this->primaryEnergyOnFirstCrystal,"primaryEnergyOnFirstCrystal/F");
  this->GetTree()->Branch("primaryInteracted",&this->primaryInteracted,"primaryInteracted/I");
  this->GetTree()->Branch("Edep", &this->Edep_vec);
  this->GetTree()->Branch("X", &this->X_vec);
  this->GetTree()->Branch("Y", &this->Y_vec);
  this->GetTree()->Branch("Z", &this->Z_vec);
  this->GetTree()->Branch("StepLength", &this->StepLength_vec);
  this->GetTree()->Branch("TrackLength", &this->TrackLength_vec);
  this->GetTree()->Branch("totalEnergyDeposited",&this->totalEnergyDeposited,"totalEnergyDeposited/F");
  this->GetTree()->Branch("totalEnergyDepositedInExtraVolume",&this->totalEnergyDepositedInExtraVolume,"totalEnergyDepositedInExtraVolume/F");
  this->GetTree()->Branch("totalEnergyDepositedInSample",&this->totalEnergyDepositedInSample,"totalEnergyDepositedInSample/F");
  this->GetTree()->Branch("totalEnergyDepositedInCase",&this->totalEnergyDepositedInCase,"totalEnergyDepositedInCase/F");
  this->GetTree()->Branch("totalEnergyDepositedInFirstCrystal",&this->totalEnergyDepositedInFirstCrystal,"totalEnergyDepositedInFirstCrystal/F");
  this->GetTree()->Branch("totalEnergyDepositedInCrystals",&this->totalEnergyDepositedInCrystals,"totalEnergyDepositedInCrystals/F");
  this->GetTree()->Branch("NumOptPhotons",&this->NumOptPhotons,"NumOptPhotons/I");
  this->GetTree()->Branch("NumCherenkovPhotons",&this->NumCherenkovPhotons,"NumCherenkovPhotons/I");

  // this->GetTree()->Branch("materials","std::vector<material_t>",&materials);

  FoundGammaVertexMomentum = false;

  detectors = new TTree("detectors","detectors");
  detectors->Branch("side",&this->Detector_boxSide ,"side/F");
  detectors->Branch("chPerSide",&this->chPerSide ,"chPerSide/I");
  detectors->Branch("id",&this->detectorID ,"id/I");
  detectors->Branch("x" ,&this->detectorX  ,"x/F" );
  detectors->Branch("y" ,&this->detectorY  ,"y/F" );
  detectors->Branch("z" ,&this->detectorZ  ,"z/F" );
  detectors->Branch("dx",&this->detectorDX ,"dx/F");
  detectors->Branch("dy",&this->detectorDY ,"dy/F");
  detectors->Branch("dz",&this->detectorDZ ,"dz/F");


  this->Clear();
}

CreateTree::~CreateTree()
{

}

Bool_t CreateTree::Write()
{
  TString filename = this->GetName();
  filename+=".root";
  TFile* file = new TFile(filename,"RECREATE");
  this->GetTree()->Write();
  file->Close();
  return true;
}

void CreateTree::Clear()
{
  Run=0;
  Event=0;
  SourceX=0;
  SourceY=0;
  SourceZ=0;
  SourceMomentumX=0;
  SourceMomentumY=0;
  SourceMomentumZ=0;
  primaryPositionOnFirstCrystalX=0;
  primaryPositionOnFirstCrystalY=0;
  primaryPositionOnFirstCrystalZ=0;
  primaryMomentumOnFirstCrystalX=0;
  primaryMomentumOnFirstCrystalY=0;
  primaryMomentumOnFirstCrystalZ=0;
  primaryEnergyOnFirstCrystal   =0;
  primaryInteracted             =0;
  primaryFirstEntranceFound     = false;
  NumOptPhotons=0;
  NumOptPhotonsReachedPMT=0;
  TotalNumOptPhotonsReachedPMT=0;
  NumCherenkovPhotons=0;
  totalEnergyDeposited=0;
  totalEnergyDepositedInExtraVolume=0;
  totalEnergyDepositedInSample=0;
  totalEnergyDepositedInCase=0;
  totalEnergyDepositedInFirstCrystal=0;
  totalEnergyDepositedInCrystals=0;
  Edep_vec.clear();
  X_vec.clear();
  Y_vec.clear();
  Z_vec.clear();
  StepLength_vec.clear();
  TrackLength_vec.clear();
}
