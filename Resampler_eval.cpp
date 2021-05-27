// reads raw eval data, resample (save resamples and individual samples), exports central value to .dat file

// minimal compilation :
//clang++ ../Resampler_eval.cc -o resamplereval -I$HOME/phd/Grid/paboyle/include -I/usr/local/include -I$HOME/phd/latan/include  -L$HOME/phd/Grid/paboyle/lib -lGrid -lHadrons -fopenmp -L/usr/local/opt/llvm/lib -L/usr/local/lib -lhdf5_cpp -lhdf5 -L$HOME/phd/latan/lib -lLatAnalyze -lcrypto -L/usr/local/opt/openssl/lib -lz -llime 

//usage example (assuming my folder tree, reading .h5 files)
// ./resamplereval ../../eval/evecs_C0_strong.evals 

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>

#include <Hadrons/Modules/MDistil/Distil.hpp>

#include <LatAnalyze/Statistics/Dataset.hpp>
#include <LatAnalyze/Io/Io.hpp>

using namespace std;

int main(int argc, char *argv[]){
    if(argc <= 1)
    {
        if(argv[0])
            std::cout << "Usage: " << argv[0] << " <radius stem>" << '\n';
        else
            std::cout << "Usage: <program name> <radius stem> " << '\n';
        return 0;
    }

    std::string inFilestem = argv[1];
    std::cout << "Filestem is " << inFilestem << std::endl;

    //################ inserting eigenvalues ##################


    std::string ensembleLabel = inFilestem.substr(inFilestem.find("evecs_")+6,2);
    std::cout << "Detected ensemble: " << ensembleLabel << std::endl;

    // C0 is standard
    int Nt = 96;
    int Nvec = 500;
    unsigned int TRAJECTORY_START = 1100;
    unsigned int TRAJECTORY_END = 1480;
    unsigned int TRAJECTORY_STEP = 20;
    if(ensembleLabel=="C1"){ //in case is C1
        Nt = 64;
        Nvec = 120;
        TRAJECTORY_START = 3140;
        TRAJECTORY_END = 3900;
        TRAJECTORY_STEP = 40;
    }

    Grid::Hadrons::MDistil::TimesliceEvals Evals(Nt,Nvec);
    int nFile = 1 + (double)(TRAJECTORY_END - TRAJECTORY_START)/TRAJECTORY_STEP;
    int nData = nFile * Nt;                 // using timeslice data as if they were uncorrelated (will bin the timslices at the end)
    Latan::Dataset<Latan::DMat>  dataset(nData);

    int iSample = 0;
    for(unsigned int s=TRAJECTORY_START; s<=TRAJECTORY_END; s+=TRAJECTORY_STEP){
        std::string filename = inFilestem + "." + std::to_string(s);
        Evals.read( filename.c_str() );
        std::cout << "Reading evals from trajectory " + std::to_string(s) + "..." << std::endl;
        for(int t=0 ; t<Nt ; t++){  // dataset[s*Nt + t](i) {conf1_t1, conf1_t2, ... , conf1_tN, conf2_t1 , ... , conf2_tN, ...} -> bin it 
            dataset[iSample].resize(Nvec, 1);       //effectively a vector
            for(int i =0 ;i<Nvec ; i++){
                dataset[iSample](i) = Evals.tensor(t,i);
            }
            iSample++;
        }
    }
    dataset.bin(Nt);    //binning time direction, remains nFile samples; write them out
    std::string outSampleFilepath = inFilestem;
    outSampleFilepath.insert(6 , "data/");
    outSampleFilepath.insert(outSampleFilepath.find("eval/")+5 , "h5/");
    outSampleFilepath.erase(outSampleFilepath.find_last_of(".") , outSampleFilepath.size());
    outSampleFilepath.replace(outSampleFilepath.find("evecs"), 5 , "evals");
    int i_conf=0;
    std::cout << "Saving individual config sample to " << outSampleFilepath + ".####.h5" << endl;
    for(unsigned int s=TRAJECTORY_START; s<=TRAJECTORY_END; s+=TRAJECTORY_STEP){
        Latan::Io::save<Latan::DMat>(dataset[i_conf], outSampleFilepath + "." + std::to_string(s) + ".h5");
        i_conf++;
    }

    //bootstraping and writing out resamples
    Latan::Index nSample = 10000;
    std::random_device  rd;
    Latan::SeedType     seed = rd();
    Latan::DMatSample   bootstrapEval = dataset.bootstrapMean(nSample , seed);
    std::string outResampleFilepath = inFilestem;
    outResampleFilepath.insert(6 , "samples/");
    outResampleFilepath.erase(outResampleFilepath.find_last_of(".") , outResampleFilepath.size());
    outResampleFilepath.replace(outResampleFilepath.find("evecs"), 5 , "evals");
    outResampleFilepath += ".h5";
    std::cout << "Writing resamples to " << outResampleFilepath << std::endl;
    Latan::Io::save<Latan::DMatSample>(bootstrapEval, outResampleFilepath);

    //wcentral value+bootErr and writing it out
    std::vector<Latan::DMat> resultEval(2);
    resultEval[0] = bootstrapEval.mean();
    resultEval[1] = bootstrapEval.variance().cwiseSqrt();
    std::string outMeanFilepath = inFilestem;
    outMeanFilepath.insert(6 , "data/");
    outMeanFilepath.erase(outMeanFilepath.find_last_of(".") , outMeanFilepath.size());
    outMeanFilepath.replace(outMeanFilepath.find("evecs"), 5 , "evals");
    outMeanFilepath += ".dat";

    const int IDENT = 16;
    std::cout << "Writing mean+-err to " << outMeanFilepath << std::endl;
    std::ofstream fileMeanOut(outMeanFilepath);
    for(int i=0 ; i<Nvec ; i++){
        fileMeanOut << std::left        << std::setw(IDENT) << std::to_string(i+1) 
                << std::scientific  << std::setw(IDENT) << resultEval[0](i) 
                << std::scientific  << std::setw(IDENT) << resultEval[1](i) << std::endl;
    }
    fileMeanOut << "### Evals from " << inFilestem;
    fileMeanOut.close();

    return 0;
}