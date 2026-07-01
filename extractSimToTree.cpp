// compile with
// g++ -o ../build/extractSimToTree extractSimToTree.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore

// small program to extrac parameters and run info from the translated output of a simulation
// run extractSimToTree without arguments for usage info

#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3I.h"
#include "TString.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TTree.h"
#include "TF2.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TTreeFormula.h"
#include "TMath.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TError.h"
#include "TEllipse.h"
#include "TFormula.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TMultiGraph.h"
#include "TCutG.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TGraphDelaunay.h"
#include "TVector.h"

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <getopt.h>

void usage()
{
  std::cout << "\t\t" << "[-i | --input]   <input file> " << std::endl
            << "\t\t" << std::endl;
}

int main (int argc, char** argv)
{
  if(argc < 2)
  {
    std::cout << argv[0] << std::endl;
    usage();
    return 1;
  }

  std::string inputFileName;

  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
			{ NULL, 0, 0, 0 }
	};

  while(1) {
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "i:o:", longOptions, &optionIndex);
		if (c == -1) {
			break;
		}
		if (c == 'i'){
			inputFileName = (char *)optarg;
    }
		else if (c == 0 && optionIndex == 0){
      inputFileName = (char *)optarg;
    }
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}

  TFile *file = TFile::Open(inputFileName.c_str());
  file->cd("Configuration");

  std::cout << std::endl;
  std::cout << ">>>>>>>> Hostname" << std::endl;
  std::cout << std::endl;
  std::cout << ((TNamed*) gDirectory->Get("Hostname"))->GetTitle() << std::endl;
  std::cout << std::endl;
  std::cout << ">>>>>>>> Working Directory" << std::endl;
  std::cout << std::endl;
  std::cout << ((TNamed*) gDirectory->Get("PWD"))->GetTitle() << std::endl;
  std::cout << std::endl;
  std::cout << ">>>>>>>> Command line" << std::endl;
  std::cout << std::endl;
  std::cout << ((TNamed*) gDirectory->Get("CommandLine"))->GetTitle() << std::endl;
  std::cout << std::endl;
  std::cout << ">>>>>>>> Parameters" << std::endl;
  std::cout << std::endl;
  std::cout << ((TNamed*) gDirectory->Get("Parameters"))->GetTitle() << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  return 0;
}
