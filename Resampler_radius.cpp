// 1) Receives a filestem, reads the psi data per config per Nvec
// 2) Reads evals per config in .h5 somewhere else
// 3) Computes radius per config and write the individual samples and resamples out in .h5
// 4) Puts radius x eval in a XYSample and fit it to a power-law; save fit out

// minimal compilation:
// clang++ ../Resampler_radius.cpp -o radiusresampler -I$HOME/phd/Grid/paboyle/include -I/usr/local/include -I$HOME/phd/latan/include  -L$HOME/phd/Grid/paboyle/lib -lGrid -lHadrons -fopenmp -L/usr/local/opt/llvm/lib -L/usr/local/lib -lhdf5_cpp -lhdf5 -L$HOME/phd/latan/lib -lLatAnalyze -lcrypto -L/usr/local/opt/openssl/lib

// usage example (suppose certain folder structure): ./radiusresampler C0 weak 

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>

#include <Hadrons/Modules/MDistil/Distil.hpp>

#include <LatAnalyze/Statistics/Dataset.hpp>
#include <LatAnalyze/Io/Io.hpp>
#include <LatAnalyze/Statistics/MatSample.hpp>
#include <LatAnalyze/Statistics/XYSampleData.hpp>
#include <LatAnalyze/Numerical/GslMinimizer.hpp>
#include <LatAnalyze/Numerical/MinuitMinimizer.hpp>

using namespace std;

double trapezoidalRule(Latan::DMat yData, int x_low , int x_high){
    // assumming equaly spaced intervals
    double dx = 1.;

    double result;
    result = (yData(x_high) + yData(x_low))/2.0;
    for(int i = x_low+1 ; i<x_high ; i++){
        result += yData(i);
    }
    result *= dx;

    return result;
}

double ratio_fine(double area_low , double area_total , double x_low , double y_low , double x_radius , double y_radius){
    double area_trapezoid = (y_radius + y_low) * (x_radius - x_low)/2.0;
    return (area_low + area_trapezoid)/area_total;
}

int main(int argc, char *argv[]){

    if(argc <= 2){
        if(argv[0])
            std::cout << "Usage: " << argv[0] << " <ensemble> <stout>" << '\n';
        else
            std::cout << "Usage: <program name> <psi stem> " << '\n';
        return 0;
    }
    std::string ensembleLabel = argv[1];
    std::string stoutStrength = argv[2];

    std::string inFilestem = "../../data/" + ensembleLabel + "/" + stoutStrength;
    std::cout << "Input ensemble: " << ensembleLabel << std::endl;
    std::cout << "Input stout strength: " << stoutStrength << std::endl;
    std::cout << "Raw psi data filestem: " << inFilestem << std::endl;
    std::cout << "Raw eval data filestem: " << "../../data/eval/h5" << std::endl;

    // C0 is the standard
    int Nt = 96;
    int Nvec = 500;
    unsigned int TRAJECTORY_START = 1100;
    unsigned int TRAJECTORY_END = 1480;
    unsigned int TRAJECTORY_STEP = 20;
    // std::vector<int> nvec_list = {1,4,8,12,16,24,32,48,64,88,112,128,152,200,260,310,380,450,490};
    std::vector<int> nvec_list = {1,4,8,12,16,24,32,48,64,128,200,260,310,380,450,500};
    if(ensembleLabel=="C1"){ //in case is C1
        Nt = 64;
        Nvec = 120;
        TRAJECTORY_START = 3140;
        TRAJECTORY_END = 3900;
        TRAJECTORY_STEP = 40;
        nvec_list = {1,2,4,6,8,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,110};
    }
    std::vector<int> nconf_list;
    for(unsigned int s=TRAJECTORY_START; s<=TRAJECTORY_END; s+=TRAJECTORY_STEP){
        nconf_list.push_back(s);
    }

    int nConfigs = nconf_list.size();
    
    //getting psi from .h5
    int nNvec = nvec_list.size();
    Latan::Dataset<Latan::Dataset<Latan::DMat>> psi(nNvec);
    std::string inPsiFilename;
    for(int i_nvec=0 ; i_nvec<nNvec ; i_nvec++){
        psi[i_nvec].resize(nConfigs);
        std::cout << std::endl<< "Loading Nvec=" << nvec_list[i_nvec] << " data" << std::endl;
        for(int i_conf=0 ; i_conf<nConfigs ; i_conf++){
            std::cout << '\r' << Latan::ProgressBar(i_conf + 1, nConfigs);
            inPsiFilename = inFilestem + "/h5/psi_nvec" + std::to_string(nvec_list[i_nvec]) + "." + std::to_string(nconf_list[i_conf]) + ".h5";
            psi[i_nvec][i_conf] = Latan::Io::load<Latan::DMat>(inPsiFilename);
            // cout << psi[i_nvec][i_conf] << endl; //DUMP
        }
    }

    //load binned eval from .h5
    std::string inEvalFilename;
    Latan::Dataset<Latan::DMat> evalWhole(nConfigs);
    Latan::Dataset<Latan::DMat> eval(nConfigs);
    std::cout << std::endl << "Loading eval data" << std::endl;
    for(int i_conf=0 ; i_conf<nConfigs ; i_conf++){
        std::cout << '\r' << Latan::ProgressBar(i_conf + 1, nConfigs);
        inEvalFilename = "../../data/eval/h5/evals_" + ensembleLabel + "_" + stoutStrength
                            + "." + std::to_string(nconf_list[i_conf]) + ".h5";
        evalWhole[i_conf] = Latan::Io::load<Latan::DMat>(inEvalFilename);
        eval[i_conf].resize(nNvec,1);
        int i_nvec = 0;
        for(int i=0 ; i<evalWhole[i_conf].size() ; i++){
            if(std::find(nvec_list.begin() , nvec_list.end() , i+1) != nvec_list.end()){
                eval[i_conf](i_nvec) = evalWhole[i_conf](i);
                i_nvec++;
            }
        }
        // cout << "########" << endl;
        // cout << eval[i_conf] << endl; //DUMP
        // cin.get();
    }
    std::cout << std::endl;

    //compute radius by integration per Nvec per config, cover 34.1% of the area from the centre
    Latan::Dataset<Latan::DMat> radius(nConfigs);
    double area_total, area_tmp , area_ratio , tol = 1e-8, y_high, y_low, y_radius, x_radius,xa,xb,ya,yb,ratio_a,ratio_radius;
    int x_low, x_high;
    for(int i_conf=0 ; i_conf<nConfigs ; i_conf++){
        radius[i_conf].resize(nNvec,1);
        for(int i_nvec=0 ; i_nvec<nNvec ; i_nvec++){
            x_high = psi[i_nvec][i_conf].size()-1;
            area_total = 2*trapezoidalRule(psi[i_nvec][i_conf], 0, x_high);
            x_high = 0;
            area_ratio = .0;
            while(area_ratio < .341){   //coarse radius
                x_high++;
                area_tmp = trapezoidalRule(psi[i_nvec][i_conf], 0, x_high);
                area_ratio = area_tmp/area_total;
            }
            x_low = x_high-1;           //fine radius is between these two
            // cout << " " << x_high << " " << area_tmp/area_total << endl;
            // fine radius
            area_tmp = trapezoidalRule(psi[i_nvec][i_conf], 0, x_low);
            xa = x_low;
            xb = x_high;
            ya = y_low = psi[i_nvec][i_conf](x_low);
            yb = y_high = psi[i_nvec][i_conf](x_high);
            x_radius = (xa + xb)/2.0;
            y_radius = (ya + yb)/2.0;
            do{
                ratio_a = ratio_fine(area_tmp , area_total , x_low , y_low , xa , ya)-0.341;
                ratio_radius = ratio_fine(area_tmp , area_total , x_low , y_low , x_radius , y_radius)-0.341;
                // cout << xa << " " << ratio_a << " "  << xb << " "  << x_radius << " "  << ratio_radius << " "  << endl;
                if( ratio_radius*ratio_a > 0 ){
                    xa = x_radius;
                    ya = y_radius;
                }
                else{
                    xb = x_radius;
                    yb = y_radius;
                }
                x_radius = (xa + xb)/2.0;
                y_radius = (ya + yb)/2.0;
                // cout << xb-xa << " "  << std::scientific << x_radius << endl;
            }while(fabs(xb-xa)>tol);
            // cout << x_radius << endl;
            radius[i_conf](i_nvec) = x_radius;
        }
        // cout << "######" << endl;
        // cout << setprecision(10) << radius[i_conf] << endl;
        // cin.get();
    }

    //resample eval and radius
    int nSample = 3000;
    random_device       rd;
    Latan::SeedType     seed = rd();
    Latan::DMatSample evalResample(nSample) , radiusResample(nSample);
    evalResample = eval.bootstrapMean(nSample, seed);
    radiusResample = radius.bootstrapMean(nSample, seed);
    // TODO save resamples .h5
    
    // TODO save central in .dat (radius x eval x nvec)
    std::string outCentralPath = inFilestem + "/radius/radius.dat";
    outCentralPath.replace(outCentralPath.find("data"),4,"figures");
    cout << "Writing to " << outCentralPath << endl;
    std::ofstream outCentralFile(outCentralPath);
    const int IDENT=16;
    std::vector<Latan::DMat> central(2);
    std::vector<Latan::DMat> bootErr(2);
    central[0] = eval.mean();   central[1] = radius.mean();
    bootErr[0] = eval.variance().cwiseSqrt();   bootErr[1] = radius.variance().cwiseSqrt();
    for(int i_nvec=0 ; i_nvec<nNvec ; i_nvec++){
        outCentralFile  << std::left    << std::setw(IDENT) << nvec_list[i_nvec]      << std::scientific
                        << std::left    << std::setw(IDENT) << central[0](i_nvec)   << std::scientific 
                        << std::left    << std::setw(IDENT) << bootErr[0](i_nvec)   << std::scientific
                        << std::left    << std::setw(IDENT) << central[1](i_nvec)   << std::scientific
                        << std::left    << std::setw(IDENT) << bootErr[1](i_nvec)   << std::endl;
    }

    // //power-law fit
    // Latan::XYSampleData radius_eval(nSample);
    // radius_eval.addXDim(nNvec);
    // radius_eval.addYDim();
    // radius_eval.setUnidimData(evalResample , radiusResample);

    // Latan::DoubleModel           power_law([](const double *x, const double *p)
    //                         {return p[0]*pow(x[0]-p[2],p[1]);}, 1, 3);
    // Latan::DVec                 init(3);
    //  init(0)=1; init(1)=-0.5 ; init(2)=0.5;
    // Latan::SampleFitResult      fit;
    // Latan::MinuitMinimizer      locMin;
    
    // // fit
    // cout << "-- fit..." << endl;

    // radius_eval.assumeYYCorrelated(false, 0, 0);
    // fit = radius_eval.fit(locMin, init, power_law);
    // fit.print();
    // cin.get();



    return 0;
}