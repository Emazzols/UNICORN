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
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//#ifdef G4MULTITHREADED
#include "G4RunManagerFactory.hh"
//#else
//#include "G4RunManager.hh"
//#include "G4MTRunManager.hh"
//#endif

#include "G4UImanager.hh"
#include "G4VModularPhysicsList.hh"
#include "G4OpticalPhysics.hh"
//#include "FTFP_BERT.hh"
#include "FTFP_BERT_HP.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4Cerenkov.hh"

#include "G4StepLimiterPhysics.hh"

//#include "G4DecayPhysics.hh"
//#include "G4RadioactiveDecayPhysics.hh"



#include "MyMaterials.hh"
#include "MySurfaces.hh"

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

// #ifdef G4VIS_USE
#include "G4VisExecutive.hh"
// #endif

// #ifdef G4UI_USE
#include "G4UIExecutive.hh"
// #endif

// #include "ConfigFile.hh"
#include "ConfigFile.h"
#include "CreateTree.hh"

#include <TNamed.h>
#include <TString.h>
#include <TROOT.h>
#include <time.h>

#include "TRandom1.h"


//#include <sys/syscall.h>    // SYS_sysctl
//#include <sys/sysctl.h>     // CTL_KERN, KERN_HOSTNAME
//#include <unistd.h>         // syscall


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// namespace {
//   void PrintUsage() {
//     G4cerr << " Usage: " << G4endl;
//     G4cerr << " dbd [-m macro ] [-u UIsession] [-t nThreads] [-r seed] "
//     << G4endl;
//     G4cerr << "   note: -t option is available only for multi-threaded mode."
//     << G4endl;
//   }
// }

G4int CreateSeed();
/*int sysctl(int *name, u_int namelen, void *oldp, size_t *oldlenp, void *newp, size_t newlen)
{
  return syscall(SYS_sysctl, name, namelen, oldp, oldlenp, newp, newlen);
}

int gethostname(char *buf, size_t buflen)
{
  int name[] = { CTL_KERN, KERN_HOSTNAME };
  size_t namelen = 2;

  return sysctl(name, namelen, buf, &buflen, NULL, 0);
}*/


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{

  if(argc<2) {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " dbd <config_file> [-m macro ] [-r seed] [-o output]" << G4endl;
    G4cerr << "   note: -m -r -o options are optional" << G4endl;
    exit(0);
  }

  // //HACK to use the dictionary easily
  // std::string executable = argv[0];
  // std::string execPath = executable.substr(0,executable.find_last_of("/")+1); //very platform dependent..
  // std::string command = ".L " + execPath + "structDictionary.C+";
  // G4cout << "command " << command << G4endl;
  // gROOT->ProcessLine(command.c_str());
  //----------------------------------------------------------//
  //  Parse the config file                                   //
  //----------------------------------------------------------//

  std::cout<<"\n"<<std::endl;
  std::cout<<"###########################################################"<<std::endl;
  std::cout<<"#                                                         #"<<std::endl;
  std::cout<<"#           New Clear PEM matrix simulation               #"<<std::endl;
  std::cout<<"#                                                         #"<<std::endl;
  std::cout<<"###########################################################"<<std::endl;
  std::cout<<"\n\n"<<std::endl;
  std::cout<<"=====>   C O N F I G U R A T I O N   <====\n"<<std::endl;

  // GET INFO ON PC AND FILES
  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);
  std::string HostNameString(hostname);
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::string PWDstring(cwd);
  //PUT THE ENTIRE CONFIG FILE IN A STRINGSTREAM, TO SAVE IT LATER IN THE OUTPUT FILE
  std::string ConfigFileName ;
  ConfigFileName = argv[1];
  std::stringstream streamConfigFile;
  std::string strConfig;
  std::ifstream inConfig;
  inConfig.open (ConfigFileName.c_str(), std::ifstream::in);
  while (std::getline(inConfig, strConfig))
  {
    streamConfigFile << strConfig << std::endl;
  }
  inConfig.close();

  //config file
  G4cout << "Configuration file: '" << argv[1] << "'"<< G4endl;
  ConfigFile config(argv[1]);

  //output root file
  std::string filename;
  config.readInto(filename,"output");
  //G4cout << "Writing data to file '" << filename << ".root '..."<< G4endl;

  //run macro
  G4String macro;
  config.readInto(macro,"macro");

  //seed for random generator
  G4long myseed = config.read<long int>("seed");

  //type of session
  G4String session;
  #ifdef G4MULTITHREADED
  G4int nThreads = 1;
  #endif

  //if passed by command line, override the filename, macro and myseed (the others are not important for us)
  for (G4int i = 0; i < argc ; i++)
  {
    if(G4String(argv[i]) == "-m") macro   = argv[i+1];
    else if(G4String(argv[i]) == "-r") myseed  = G4UIcommand::ConvertToInt(argv[i+1]);
    else if(G4String(argv[i]) == "-o") filename = argv[i+1];
  }

  if(myseed==-1)
  {
    G4cout<<"Creating random seed..."<<G4endl;
    myseed=CreateSeed();
  }

  std::stringstream streamSeed;
  streamSeed << myseed;

  //PUT THE ENTIRE MACRO FILE IN A STRINGSTREAM, TO SAVE IT LATER IN THE OUTPUT FILE
  std::string MacroFileName ;
  MacroFileName = macro;
  std::stringstream streamMacroFile;
  std::string strMacro;
  std::ifstream inMacro;
  inMacro.open (MacroFileName.c_str(), std::ifstream::in);
  while (std::getline(inMacro, strMacro))
  {
    streamMacroFile << strMacro << std::endl;
  }
  inMacro.close();

  // also save seed


  G4cout << "Reading macro '" << macro << "' ..." << G4endl;
  G4cout << "Random seed : "  << myseed <<G4endl;
  G4cout << "Output file '"   << filename << ".root '..." << G4endl;




  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  auto* runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
#ifdef G4MULTITHREADED
  if ( nThreads > 0 ) {
    runManager->SetNumberOfThreads(nThreads);
  }
#endif

  // Seed the random number generator manually
  G4Random::setTheSeed(myseed);


  //create output ttree
  //   CreateTree* mytree = new CreateTree("StandardTree",ncrystalx,ncrystaly,nmppcx,nmppcy,true,true);

  CreateTree* mytree = new CreateTree("StandardTree",1,1,1,1);
  //((CreateTree*)mytree)->SetModuleElements(ncrystalx,ncrystaly,nmppcx,nmppcy);
  //mytree->SetModuleElements(ncrystalx,ncrystaly,nmppcx,nmppcy);

  //save the seed in the ttree
  CreateTree::Instance()->Seed = myseed;

  // verbose
  int verbose = config.read<int> ("verbose",0);
  CreateTree::Instance()->verbose = verbose;

  // initialize the materials and surfaces classes
  MyMaterials *Materials = new MyMaterials();
  MySurfaces  *Surfaces  = new MySurfaces();
  int scint = config.read<int> ("propagateScintillation",1);
  int cher  = config.read<int> ("propagateCerenkov",1);

  //create the detector construction
  DetectorConstruction* detector = new DetectorConstruction(config,Materials,Surfaces);


  // Set mandatory initialization classes
  //
  // Detector construction
  runManager-> SetUserInitialization(detector);
  // Physics list
  G4VModularPhysicsList* physicsList = new FTFP_BERT_HP();//option with HP automatically activates radioactive decay

  // EM physics
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());

  // Need to register optical physics as well
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  physicsList->RegisterPhysics(opticalPhysics);

  int manual_step = config.read<int>("manual_step",0); 
  if(manual_step == 1){ 
    G4StepLimiterPhysics* steplimiterPhysics = new G4StepLimiterPhysics(); 
    physicsList->RegisterPhysics(steplimiterPhysics); 
    } 

  int manual_cuts = config.read<int>("manual_cuts",0); 
  G4double cuts = config.read<double>("cuts",0.025);
  std::stringstream streamCut;
  streamCut << cuts;
  if(manual_cuts == 1){
    physicsList->SetDefaultCutValue(cuts); 
    physicsList->SetCuts();
  }





//  G4OpticalParameters::Instance()->SetScintillationExcitationRatio(0.0);

  G4OpticalParameters::Instance()->SetCerenkovMaxPhotonsPerStep(1000);
  G4OpticalParameters::Instance()->SetCerenkovMaxBetaChange(10.0);

  G4OpticalParameters::Instance()->SetCerenkovTrackSecondariesFirst(true);
  G4OpticalParameters::Instance()->SetScintTrackSecondariesFirst(true);
  G4OpticalParameters::Instance()->SetScintFiniteRiseTime(true);
  G4OpticalParameters::Instance()->SetScintVerboseLevel(0);
  //G4OpticalParameters::Instance()->SetCerenkovVerboseLevel(2);

  runManager->SetUserInitialization(physicsList);


  // User action initialization
  runManager->SetUserInitialization(new ActionInitialization(config,Materials,Surfaces,scint,cher));
  // Initialize G4 kernel
  //
  runManager->Initialize(); 

  G4ProductionCutsTable::GetProductionCutsTable()->DumpCouples(); 

  // #ifdef G4VIS_USE
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
  // #endif


  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( macro.size() ) {
    // Batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else // Define UI session for interactive mode
  {
    // #ifdef G4UI_USE
    G4UIExecutive* ui = 0;
    ui = new G4UIExecutive(argc,argv,session);
    // #ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute vis.mac");
    // #else
    // UImanager->ApplyCommand("/control/execute g4matrix.in");
    // #endif
    // if (ui->IsGUI())
    //   UImanager->ApplyCommand("/control/execute gui.mac");
    ui->SessionStart();
    delete ui;
    // #endif
  }


  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !




  // #ifdef G4VIS_USE
  // delete visManager;

  // #endif
  delete visManager;
  delete runManager;

  //output root file
  std::stringstream outfileName;
  outfileName << filename << ".root";

  TFile* outfile = new TFile(TString(outfileName.str().c_str()),"RECREATE");
  outfile->cd();
  G4cout<<"Writing tree to file "<< outfileName.str() <<" ..."<<G4endl;
  mytree->GetTree()->Write();
  std::vector<std::string> materialsList = mytree->materialsList;
  gDirectory->WriteObject(&materialsList, "materials");
  mytree->detectors->Write();
  mytree->photons->Write();
  mytree->photonsTotal->Write();
  TDirectory *configDir = outfile->mkdir("Configuration");
  configDir->cd();
  TNamed SeedNameD("Seed",streamSeed.str().c_str());
  TNamed CutNameD("Cuts",streamCut.str().c_str());
  TNamed HostNameD("Hostname",HostNameString.c_str());
  TNamed PWDNameD("PWD",PWDstring.c_str());
  TNamed ConfigNameD("ConfigFile",streamConfigFile.str().c_str());
  TNamed GpsNameD("GpsFile",streamMacroFile.str().c_str());
  SeedNameD.Write();
  CutNameD.Write();
  HostNameD.Write();
  PWDNameD.Write();
  ConfigNameD.Write();
  GpsNameD.Write();
  //outfile->Write();
  outfile->Close();


  return 0;
}


G4int CreateSeed()
{
  TRandom1 rangen;
  long int ss = time(0);
  std::cout<<"Time : "<<ss<<std::endl;
  ss+=getpid();
  std::cout<<"PID  : "<<getpid()<<std::endl;
  //
  FILE * fp = fopen ("/proc/uptime", "r");
  int uptime,upsecs;
  if (fp != NULL)
  {
    char buf[BUFSIZ];
    int res;
    char *b = fgets (buf, BUFSIZ, fp);
    if (b == buf)
    {
      /* The following sscanf must use the C locale.  */
      setlocale (LC_NUMERIC, "C");
      res = sscanf (buf, "%i", &upsecs);
      setlocale (LC_NUMERIC, "");
      if (res == 1)
        uptime = (time_t) upsecs;
    }

    fclose (fp);
  }

  std::cout<<"Uptime: "<<upsecs<<std::endl;
  ss+=upsecs;
  //
  std::cout<<"Seed for srand: "<<ss<<std::endl;
  srand(ss);
  rangen.SetSeed(rand());
  long int seed = round(1000000*rangen.Uniform());
  return seed;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
