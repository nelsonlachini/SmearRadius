/*
minimal compilation : clang++ ../Smear_radius_logistic.cc -o logisticsmearradius -I$HOME/phd/Grid/paboyle/include -I/usr/local//include -I/Users/s2000761/phd/include -llime -L/usr/local/lib -lGrid -L/Users/s2000761/phd/Grid/paboyle/lib -lHadrons -fopenmp -L/usr/local/opt/llvm/lib -lcrypto -L/usr/local/opt/openssl/lib -lz
*/

#include <Grid/Grid.h>
// #include <Hadrons/Global.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Modules/MDistil/Distil.hpp>

#include <Grid/qcd/action/gauge/GaugeImplTypes.h>

using namespace std;
using namespace Grid;

double cutoff(const double xi , double lambda, const double xi_cutoff, const double xi_i){
    lambda /= xi_cutoff;
    double p = 80; // 15;
    // lambda*(Nv - Ncutoff) should be ~ 10 so that the contributions from the tail are negligible
    // return 1.0/(1.0 + exp(lambda * (xi - xi_cutoff)));
    // return 1.0/(1.0 + exp(lambda * (xi - xi_cutoff))) * exp(-0.5*pow((xi - xi_cutoff)/xi_cutoff,2));
    // if(xi <= xi_cutoff)
    //     return 1.2 * exp(lambda*(xi/xi_cutoff-1));
    // else
    //     return 0;
    return (0.2 * exp(-lambda * (xi/xi_i - 1)) + 1.75 * exp(lambda * (xi/xi_cutoff - 1))) / (1 + exp(p * lambda * (xi/xi_cutoff - 1)));
}

int main (int argc, char ** argv)
{
    bool doMock = false;

    // core
    Grid_init(&argc,&argv);
    Coordinate latt_size  = GridDefaultLatt();
    Coordinate simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    GridCartesian   grid(latt_size,simd_layout,mpi_layout);
    const int vol = grid.gSites();
    
    // temporaries
    using LapEvecs = Grid::Hadrons::EigenPack<LatticeColourVector>;
    LapEvecs eigen4d(1 , &grid);
    LapEvecs eigenTemp(1 , &grid);
    LatticeColourMatrix outerPsi_r(&grid);
    

    // file handling
    const int IDENT = 16;
    std::string stout_strength = "weak";
    std::string ensembleLabel = "C0";
    std::string eigenFilepath, outFilepath;
    std::string eigenFilestem = "/tessfs1/work/dp008/dp008/dc-lach1/lapevecs/evecs_"
                                 + ensembleLabel + "_" + stout_strength;
    std::ofstream f;
    std::ifstream evalFile;


    // parameters
    unsigned int TRAJECTORY_START;
    unsigned int TRAJECTORY_END;
    unsigned int TRAJECTORY_STEP;
    std::vector<int> nvec_list;
    if(ensembleLabel == "C0"){
        nvec_list = {1,4,8,12,16,24,32,48,64,88,112,128,152,200,260,310,380,450,490};
        TRAJECTORY_START = 1100;
        TRAJECTORY_END = 1481;
        TRAJECTORY_STEP = 20;
    }
    else{
        nvec_list = {1,2,4,6,8,10,12,15,17,20,23,25,30,35,40,45,50,55,60,70,80,90,100,110};
        unsigned int TRAJECTORY_START = 3140;
        unsigned int TRAJECTORY_END = 3901;
        unsigned int TRAJECTORY_STEP = 40;
    }
    if(doMock){
        TRAJECTORY_START = 2000;
        TRAJECTORY_END = 2000;
        TRAJECTORY_STEP = 1;
        eigenFilestem = "../../evec/eigen_test";
        nvec_list = {1,5};
    }

    std::vector<RealD> Psi_r(latt_size[Xdir]/2 + 1 , 0.0);    //final Psi

    int Nvec_total;
    double lambda = 2;
    int delta_hard_cutoff = 10;
    if(doMock){
        delta_hard_cutoff = 0;
        lambda = 0;
    }

    const int Nvec_max = nvec_list.back() + delta_hard_cutoff;

    for(unsigned int i_conf=TRAJECTORY_START ; i_conf<=TRAJECTORY_END ; i_conf+=TRAJECTORY_STEP){
        std::cout << GridLogMessage << "Starting measurement on trajectory " << i_conf << "..." << std::endl;

        // reseting Psi_r

        //reading eigenvectors
        eigenFilepath = eigenFilestem + "." + std::to_string(i_conf);
        std::cout << GridLogMessage << "Reading Nvec=" << Nvec_max << " eigenvectors from " << eigenFilepath << std::endl;
        eigen4d.resize(Nvec_max , &grid);
        eigen4d.read(eigenFilepath , false);

        //reading central eigenvalue
        std::vector<double> evalNvec;
        std::vector<double> eval;
        std::vector<double> err_eval;
        if(!doMock){

        std::string inEvalstem = "../../data/eval/evals_";
        inEvalstem += ensembleLabel + "_" + stout_strength + ".dat";

        evalFile.open(inEvalstem);
        std::cout << "Reading eigenvalues from " << inEvalstem << "..." << std::endl;
        std::string line;
        while(evalFile){
            std::getline(evalFile , line);
            if(line[0] != '#'){       // # is reserved for metadata lines
                evalNvec.push_back(std::stoi(line.substr(0 , IDENT) ));
                eval.push_back(std::stod(line.substr(IDENT , 2*IDENT)));
                err_eval.push_back(std::stod(line.substr(2*IDENT , 3*IDENT)));
                std::cout << GridLogMessage 
                    << evalNvec.back() << " " << eval.back() << " " << err_eval.back() << endl;
            }
        }
        evalFile.clear();
        evalFile.close();
        std::cout << "#############" << std::endl;
        }

        

        std::cout << GridLogMessage << "Lattice time extension = " << latt_size[Tdir] << std::endl;
        for(int Nvec_smooth : nvec_list){
            Nvec_total = Nvec_smooth + delta_hard_cutoff;

            //reweighting eigenvectors according to cutoff()
            eigenTemp.resize(Nvec_total, &grid);
            if(!doMock){
                if(i_conf==TRAJECTORY_START){   //save only first ocurrence
                    outFilepath = "../../data/filter/" + ensembleLabel + "/" + stout_strength + "/cutoff";
                    outFilepath += "/eval" + std::to_string(Nvec_smooth) + ".dat";
                    f.open(outFilepath);
                }
                for(int i=0 ; i<Nvec_total ; i++){
                    eigenTemp.evec[i] = eigen4d.evec[i] * cutoff(eval[i],lambda,eval[Nvec_smooth-1],eval[0]);
                    cout << eval[i] << " " << cutoff(eval[i],lambda,eval[Nvec_smooth-1],eval[0]) << endl;
                    if(i_conf==TRAJECTORY_START)
                        f << eval[i] << " " << cutoff(eval[i],lambda,eval[Nvec_smooth-1],eval[0]) << endl;
                }
                if(i_conf==TRAJECTORY_START){
                    cout << "Saved cutoff profile to " << outFilepath << endl;
                    f.close();
                }
            }
            else
            for(int i=0 ; i<Nvec_total ; i++){
                eigenTemp.evec[i] = eigen4d.evec[i];
            }


            //measuring psi
            // std::cout << GridLogMessage << "Measuring spatial distribution of distil operator with logistic cuoff, Nvec_smooth=" <<  Nvec_smooth << "..." << std::endl;
            for(int r_shift=0 ; r_shift<latt_size[Xdir]/2+1 ; r_shift++){   //assuming spatial lattice is symmetric and pbc
                Psi_r[r_shift] = 0.0;
                for(int spatial_dir=0 ; spatial_dir<3 ; spatial_dir++){     //assuming pbc
                    outerPsi_r = Zero();
                    for(int i=0 ; i<Nvec_total ; i++){
                        outerPsi_r += outerProduct(eigenTemp.evec[i] , Cshift(eigenTemp.evec[i] , spatial_dir , r_shift));
                        // std::cout<<GridLogMessage<< "Summing out Psi (trajectory " << i_conf << ", Nvec=" << i+1 << ", r_shift=" << r_shift << ", spatial_dir=" << spatial_dir << "):" << std::endl;
                    }
                    std::cout<<GridLogMessage<< "Summing out Psi (trajectory " << i_conf << ", Nvec=" << Nvec_total << ", r_shift=" << r_shift << ", spatial_dir=" << spatial_dir << "):" << std::endl;
                    Psi_r[r_shift] +=  toReal(
                                        TensorRemove(
                                        sum(sqrt(
                                        trace(outerPsi_r * adj(outerPsi_r)
                                            )))))/(3.0*vol*(Nvec_total)*sqrt(Nc));
                }
            }

            std::cout<<GridLogMessage<< "Final result dump (trajectory " << i_conf << "):" << std::endl;
            
            for(int r_shift=0 ; r_shift<latt_size[Xdir]/2+1 ; r_shift++){
                std::cout << GridLogMessage << left << setw(IDENT) << r_shift << setw(IDENT) << Psi_r[r_shift] << std::endl;
            }

            // writing out
            outFilepath = "../../data/filter/" + ensembleLabel + "/" + stout_strength;
            if(doMock)
                outFilepath = "../../data/test/";
            outFilepath += "/psi_nvec" + std::to_string(Nvec_smooth) + "." + std::to_string(i_conf) + ".dat";
            f.open(outFilepath);
            std::cout << GridLogMessage << "Writing to " << outFilepath << std::endl;
            for(int r_shift=0 ; r_shift<=latt_size[Xdir]/2 ; r_shift++){
                f << left << setw(IDENT) << r_shift << setw(IDENT) << std::scientific << Psi_r[r_shift]/Psi_r[0] << std::endl;
            }
            f << "#################################################################" << std::endl;
            f << "# Spatial distribution of distillation operator (unnormalised)#" << std::endl;
            f << "#################################################################" << std::endl;
            for(int r_shift=0 ; r_shift<=latt_size[Xdir]/2 ; r_shift++){
                f << "### " << left << setw(IDENT) << r_shift << setw(IDENT) << std::scientific << Psi_r[r_shift] << std::endl;
            }
            f << "### Eigenvectors from " + eigenFilestem << std::endl;
            f << "### Lattice = " << latt_size[Xdir] << "^3 x " << latt_size[Tdir] << std::endl;
            f << "### Logistic cutoff lambda =" << lambda << std::endl;
            f << "### Nvec_smooth =" << Nvec_smooth ;
            f.close();
            f.clear();
        
            std::cout << GridLogMessage << "Finished" << std::endl;
            std::cout << GridLogMessage << "################" << std::endl;
        }
    }

    std::cout << GridLogMessage << "Process finished, finalizing Grid..." << std::endl;

    Grid_finalize();
}



 