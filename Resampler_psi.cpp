// receive a given filestem, creates psi manifest files, reads raw .dat data, and exports resampled .h5 file (and individual h5 samples) and save central value as a .dat file for plotting

// minimal compilation:
// clang++ ../Resampler_psi.cpp -o resamplerpsi -I$HOME/phd/Grid/paboyle/include -I/usr/local/include -I$HOME/phd/latan/include  -L$HOME/phd/Grid/paboyle/lib -lGrid -lHadrons -fopenmp -L/usr/local/opt/llvm/lib -L/usr/local/lib -lhdf5_cpp -lhdf5 -L$HOME/phd/latan/lib -lLatAnalyze -lcrypto -L/usr/local/opt/openssl/lib


// individual usage example:
// ./resamplerpsi ../../data/C1/weak --nvec=1

// combine with for loop:
// for i in 1 2 4 6 8 10 15 20 25 30 35 40 45 50 55 60 70 80 90 100 110 ; do ./resamplerpsi ../../data/C1/weak --nvec=${i} ; done

#include <Grid/Grid.h>
#include <Hadrons/Global.hpp>
#include <Hadrons/Modules/MContraction/Meson.hpp>
#include <LatAnalyze/Core/OptParser.hpp>
#include <LatAnalyze/Statistics/Dataset.hpp>
#include <LatAnalyze/Io/Io.hpp>

using namespace std;
using namespace Latan;
using namespace Grid;
typedef Hadrons::MContraction::Meson::Result Result;

#define DEF_NSAMPLE "10000"

int main(int argc, char *argv[])
{
  // parse command line ////////////////////////////////////////////////////////
  OptParser              opt;
  bool                   parsed;
  string                 manFilename , inFilestem , outIndividualFilename;;
  Latan::Index           binSize, nSample;
  random_device          rd;
  SeedType               seed = rd();
  
  opt.addOption("n", "nvec"     , OptParser::OptType::value,   false,
              "number of Lap eigenvectors" , "");
  opt.addOption("b", "bin"       , OptParser::OptType::value,   true,
                "bin size", "1");
  
  parsed = opt.parse(argc, argv);
  binSize     = opt.optionValue<Latan::Index>("b");
  inFilestem  = opt.getArgs()[0];
  int  nvec   = opt.optionValue<Latan::Index>("n");

  // generating manifest of psi files
  std::string command = "cd " + inFilestem + "; find `pwd` | grep psi_nvec" + std::to_string(nvec) + "\\\\..*\\\\.dat | sort > psi_nvec" + std::to_string(nvec) + ".manifest";
  std::system(command.c_str());
  manFilename = inFilestem + "/psi_nvec" + std::to_string(nvec) + ".manifest";
  std::string outResampleFilename = manFilename;

  std::cout << "-- loading data... " << manFilename << std::endl;
  std::vector<std::string> inFilename = readManifest(manFilename);             //gets a list of the file paths contained in the manifest file and writes this to inFilename

  //load data /////////////////////////////////////////////////////////////////
  std:vector<double> data;
  std::string stemp;
  std::ifstream f; 

  manFilename.erase(0,manFilename.find("nvec") + 4);
  manFilename.erase(manFilename.find(".") , manFilename.size());

  unsigned int nFile = inFilename.size();                              //gets the number of files listed in that manifest file

  int Nr;
  Dataset<DMat>  dataset(nFile);

  for (unsigned int i = 0; i < nFile; ++i)                                    //loop trough the list of configuration files
  {
    std::cout << "Reading " << inFilename[i] << std::endl;
    std::cout << '\r' << Latan::ProgressBar(i + 1, nFile);
    std::cout << std::endl;

    f.open(inFilename[i]);
    while(f)
    {
      std::getline(f , stemp);
      if(stemp[0] != '#')
      {                                                                           // # is reserved for metadata lines
        stemp.erase(0, 12);
        std::cout << stemp << std::endl;    //DUMP
        data.push_back(std::stod(stemp));
      }
    }
    f.close();
    f.clear();

    Nr = data.size();

    //converting to a Latan type
    dataset[i].resize(Nr,1);
    for(int j=0 ; j<Nr ; j++){
      dataset[i](j) = data[j];
    }
    // converting my data to .h5 for each sample
    outIndividualFilename = inFilename[i];
    outIndividualFilename.insert(outIndividualFilename.find("/psi_nvec") , "/h5");  // save it to h5/ subfolder
    outIndividualFilename.erase(outIndividualFilename.find(".dat"),outIndividualFilename.size());
    outIndividualFilename += ".h5";

    std::cout << "Saving trajectory " << std::to_string(i)  << " to " << outIndividualFilename << std::endl;
    Io::save<DMat>(dataset[i] , outIndividualFilename);
    data.clear();
  }
  std::cout << "Saving to " << outIndividualFilename << std::endl;
  std::cout << std::endl;

  //resampling it
  nSample = std::stoi(DEF_NSAMPLE);
  cout << "-- resampling (" << nSample << " samples)..." << endl;
  DMatSample out;
  out = dataset.bootstrapMean(nSample, seed);
  DMat central = out.mean();
  DMat bootErr = out.variance().cwiseSqrt();
  cout << "Mean = " << endl << central << endl;
  cout << "Std deviation = " << endl << bootErr << endl;
  // save central+-err to .dat files for plotting
  const int IDENT = 16;
  std::string outMeanFilepath = inFilestem;
  outMeanFilepath.replace(outMeanFilepath.find("data") , 4 , "figures");
  outMeanFilepath += "/psi_nvec" + std::to_string(nvec);
  Latan::mkdir(outMeanFilepath);
  outMeanFilepath += "/points_3.dat";
  std::cout << "Writing central value to " << outMeanFilepath << std::endl;
  std::ofstream fileMeanOut(outMeanFilepath);
  for(int i=0 ; i<Nr ; i++){
      fileMeanOut << std::left        << std::setw(IDENT) << std::to_string(i)
                  << std::scientific  << std::setw(IDENT) << std::to_string(central(i))
                  << std::scientific  << std::setw(IDENT) << std::to_string(bootErr(i)) << std::endl;
  }
  fileMeanOut.close();
  // saving resamplers
  outResampleFilename.erase( outResampleFilename.find_last_of(".") , outResampleFilename.size() ) ;
  outResampleFilename.replace(outResampleFilename.find("data") , 4 , "samples");
  outResampleFilename += ".h5";
  Io::save<DMatSample>(out, outResampleFilename);
  cout << "Sample saved to '" << outResampleFilename << "'" << endl;


  return EXIT_SUCCESS;
}