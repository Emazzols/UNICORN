#include <dirent.h>
// list files in directory
// taken from
// http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
void read_directory(const std::string& name, std::vector<std::string> &v)
{
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}


void simpleData(std::string inputFilePrefix)
{
  std::string outFileName = "summary_" + inputFilePrefix + ".root";
  TFile* fOut = new TFile(outFileName.c_str(),"RECREATE");

  std::vector<Double_t> chNumber;
  std::vector<Double_t> PeakAvg;
  std::vector<Double_t> EnRes;

  std::vector<std::string> v;
  read_directory("./", v);
  // extract files with correct prefix
  std::vector<std::string> listInputFiles;
  for(unsigned int i = 0 ; i < v.size() ; i++)
  {
    if(!v[i].compare(0,inputFilePrefix.size(),inputFilePrefix))
    {
      listInputFiles.push_back("./" + v[i]);
    }
  }

  for(unsigned int i = 0 ; i < listInputFiles.size(); i++)
  {
    std::cout << "Reading file " << listInputFiles[i] << std::endl;
    TFile *_file0 = TFile::Open(listInputFiles[i].c_str());
    TTree *tree = (TTree*) _file0->Get("adc");
    // find number of channels
    TObjArray *leavescopy = tree->GetListOfLeaves();
    int nLeaves = leavescopy->GetEntries();
    std::vector<std::string> leavesName;
    // fill a vector with the leaves names
    for(int i = 0 ; i < nLeaves ; i++)
    {
      leavesName.push_back(leavescopy->At(i)->GetName());
    }
    // count the entries that start with "ch"
    int numOfCh = 0;
    std::string det_prefix("ch");
    for(int i = 0 ; i < nLeaves ; i++)
    {
      if (!leavesName[i].compare(0, det_prefix.size(), det_prefix))
      numOfCh++;
    }

    // sum channels
    std::stringstream sumString;
    sumString << "(ch0";
    for(int j = 1; j < numOfCh ; j++)
    {
      sumString << "+ch" << j;
    }
    sumString << ")";
    std::cout << sumString.str() << std::endl;

    // photopeak
    std::stringstream sname;
    sname << numOfCh << "ch spectrum";
    TH1F* spectrum = new TH1F(sname.str().c_str(),sname.str().c_str(),250,0,20000);
    sname.str("");
    sname << sumString.str() << " >> "<< numOfCh << "ch spectrum";
    tree->Draw(sname.str().c_str());

    TSpectrum *sPeak;
    sPeak = new TSpectrum(5);
    Int_t PeaksN = sPeak->Search(spectrum,1,"",0.5);
    Double_t *Peaks  = sPeak->GetPositionX();
    Double_t *PeaksY = sPeak->GetPositionY();
    float maxPeak = 0.0;
    int peakID = 0;
    for (int peakCounter = 0 ; peakCounter < PeaksN ; peakCounter++ )
    {
      if(Peaks[peakCounter] > maxPeak) // rightmost peak
      {
        maxPeak = Peaks[peakCounter];
        peakID = peakCounter;
      }
    }
    // assume en res aroun 10% fwhm
    TF1 *gaussPeak = new TF1("gaussPeak","gaus", Peaks[peakID] - 0.10*Peaks[peakID], Peaks[peakID] + 0.10*Peaks[peakID]);
    gaussPeak->SetParameter(0,PeaksY[peakID]); // heigth
    gaussPeak->SetParameter(1,Peaks[peakID]); //mean
    gaussPeak->SetParameter(2,(0.10*Peaks[peakID])/2.355); // sigma
    spectrum->Fit(gaussPeak,"Q","", Peaks[peakID] - 0.10*Peaks[peakID], Peaks[peakID] + 0.10*Peaks[peakID]);

    chNumber.push_back(numOfCh);
    PeakAvg.push_back(gaussPeak->GetParameter(1));
    EnRes.push_back((gaussPeak->GetParameter(2)*2.355)/gaussPeak->GetParameter(1));
    // save spectrum
    fOut->cd();
    spectrum->Write();



  }

  TGraph *gPeakAvg = new TGraph(chNumber.size(),&chNumber[0],&PeakAvg[0]);
  gPeakAvg->SetName("Peak position");
  gPeakAvg->GetXaxis()->SetTitle("Number of SiPM per side");
  gPeakAvg->GetYaxis()->SetTitle("Peak Position [N of Photons Detected]");
  TGraph *gEnRes= new TGraph(chNumber.size(),&chNumber[0],&EnRes[0]);
  gEnRes->SetName("Energy Resolution");
  gEnRes->GetXaxis()->SetTitle("Number of SiPM per side");
  gEnRes->GetYaxis()->SetTitle("En Res FWHM [fraction]");

  gPeakAvg->Write();
  gEnRes->Write();
  fOut->Close();
}
