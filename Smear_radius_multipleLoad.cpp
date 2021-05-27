/*
minimal compilation : clang++ ../Smear_radius.cc -o smearradius -I$HOME/phd/Grid/paboyle/include -I/usr/local//include -I/Users/s2000761/phd/include -llime -L/usr/local/lib -lGrid -L/Users/s2000761/phd/Grid/paboyle/lib -lHadrons -fopenmp -L/usr/local/opt/llvm/lib -lcrypto -L/usr/local/opt/openssl/lib -lz
*/

#include <Grid/Grid.h>
#include <Hadrons/Global.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Modules/MDistil/Distil.hpp>

using namespace std;
using namespace Grid;
using namespace Hadrons;

int main (int argc, char ** argv)
{
    bool doMock = true;

    // core
    Grid_init(&argc,&argv);
    Coordinate latt_size  = GridDefaultLatt();
    Coordinate simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    GridCartesian   grid(latt_size,simd_layout,mpi_layout);
    const int vol = grid.gSites();
    
    /////////////////////// temporaries/eigen
    std::unique_ptr<LatticeColourVector> ioBuf{nullptr};
    ScidacReader binReader;
    Hadrons::PackRecord record;
    LatticeColourVector evec(&grid);
    Real eval;
    int i_nvec;
    std::vector<LatticeColourMatrix> outerPsi_r(latt_size[Xdir]/2 + 1 , &grid);


    // file handling
    std::string stout_strength = "weak";
    std::string eigenFilepath, outFilepath;
    std::string eigenFilestem = "/tessfs1/work/dp008/dp008/dc-lach1/lapevecs/evecs_C0_" + stout_strength;
    // std::string eigenFilestem = "/tessfs1/work/dp008/dp008/dc-lach1/lapevecs/evecs_C1_" + stout_strength;
    std::ofstream f;

    // parameters
    // std::vector<int> nvec_list = {1,4,8,12,16,24,32,48,64,88,112,128,152,200,260,310,380,450,490};
    std::vector<int> nvec_list = {1,2,4,6,8,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,110,120};
    std::vector<std::vector<RealD>> Psi_r(latt_size[Xdir]/2 + 1 , std::vector<RealD>(nvec_list.size(), 0.0));    //final Psi
    // unsigned int TRAJECTORY_START = 1100;
    // unsigned int TRAJECTORY_END = 1481;
    // unsigned int TRAJECTORY_STEP = 20;
    unsigned int TRAJECTORY_START = 3140;
    unsigned int TRAJECTORY_END = 3901;
    unsigned int TRAJECTORY_STEP = 40;
    if(doMock){
        TRAJECTORY_START = 2000;
        TRAJECTORY_END = 2001;
        TRAJECTORY_STEP = 1;
        eigenFilestem = "../../evec/eigen_test";
        nvec_list = {1,5};
    }


    for(unsigned int i_conf=TRAJECTORY_START ; i_conf<=TRAJECTORY_END ; i_conf+=TRAJECTORY_STEP){
        eigenFilepath = eigenFilestem + "." + std::to_string(i_conf) + ".bin";

        std::cout << GridLogMessage << "Starting measurement on trajectory " << i_conf << "..." << std::endl;
        std::cout << GridLogMessage << "Lattice time extension = " << latt_size[Tdir] << std::endl;
        std::cout << GridLogMessage << "Measuring spatial distribution of distil operator, Nvec=" << nvec_list << "..." << std::endl;
        
        for(int r_shift=0 ; r_shift<latt_size[Xdir]/2 + 1; r_shift++){
            for(int i=0 ; i < nvec_list.size() ; i++){
                Psi_r[r_shift][i] = 0.0;
            }
        }

        //execution
        for(int spatial_dir=0 ; spatial_dir<3 ; spatial_dir++){
            std::cout << GridLogMessage << "Summing spatial direction " << spatial_dir <<  std::endl;
            //initialising grids to zero
            for(int r_shift=0 ; r_shift<latt_size[Xdir]/2 + 1 ; r_shift++){
                outerPsi_r[r_shift] = Zero();
            }
            std::cout << GridLogMessage << "Temp initialised "<< std::endl;

            eigenFilepath = eigenFilestem + "." + std::to_string(i_conf) + ".bin";  // need to open file twice to account for 'bug' on reader
            // binReader.open(eigenFilepath); binReader.close();
            binReader.open(eigenFilepath);
            i_nvec = 0;
            for(int i=0 ; i < nvec_list.back() ; i++){
                Hadrons::EigenPackIo::readElement(evec, eval, i, binReader, ioBuf.get());
                std::cout << GridLogMessage << "Adding eigenvector " << i+1 << "/" << nvec_list.back() << std::endl;
                for(int r_shift=0 ; r_shift<latt_size[Xdir]/2 + 1; r_shift++){   
                    outerPsi_r[r_shift] += outerProduct(evec , Cshift(evec , spatial_dir , r_shift));
                }
                if(std::find(nvec_list.begin(), nvec_list.end() , i+1) != nvec_list.end()){
                    for(int r_shift=0 ; r_shift<latt_size[Xdir]/2 + 1; r_shift++){
                        Psi_r[r_shift][i_nvec] += toReal(
                                                    TensorRemove(
                                                    sum(
                                                    sqrt(
                                                    trace(
                                                        outerPsi_r[r_shift] *
                                                        adj(outerPsi_r[r_shift])
                                                    ))))) / (3.0*vol*(i+1)*sqrt(Nc)) ;
                    }
                    i_nvec++;
                }
            } //Nvec loop
            binReader.close();
            
        } //spatial direction loop

        // result dump
        std::cout<<GridLogMessage<< "###### DUMP ########" << std::endl;
        for(i_nvec=0; i_nvec < nvec_list.size() ; i_nvec++){
            std::cout<<GridLogMessage<< "Nvec=" << nvec_list[i_nvec] << std::endl;
            for(int r_shift=0 ; r_shift<latt_size[Xdir]/2 + 1; r_shift++){
                std::cout<<GridLogMessage<< "(r , Psi_r) = (" << r_shift << " , " << std::scientific << Psi_r[r_shift][i_nvec] << ")" << std::endl;
            }
            
        }

        for(i_nvec=0 ; i_nvec<nvec_list.size() ; i_nvec++){
            // output files
            // outFilepath = "../../data/C0/" + stout_strength;
            outFilepath = "../../data/C1/" + stout_strength;
            if(doMock)
                outFilepath = "../../data/test/";
            outFilepath += "/psi_nvec" + std::to_string(nvec_list[i_nvec]) + "." + std::to_string(i_conf) + ".dat";
            f.open(outFilepath);
            std::cout << GridLogMessage << "Writing to " << outFilepath << std::endl;
            const int IDENT = 16;
            for(int r_shift=0 ; r_shift<=latt_size[Xdir]/2 ; r_shift++){
                f << left << setw(IDENT) << r_shift << setw(IDENT) << std::scientific << Psi_r[r_shift][i_nvec]/Psi_r[0][i_nvec] << std::endl;
            }
            f << "#################################################################" << std::endl;
            f << "# Spatial distribution of distillation operator (unnormalised)#" << std::endl;
            f << "#################################################################" << std::endl;
            for(int r_shift=0 ; r_shift<=latt_size[Xdir]/2 ; r_shift++){
                f << "### " << left << setw(IDENT) << r_shift << setw(IDENT) << std::scientific << Psi_r[r_shift][i_nvec] << std::endl;
            }
            f << "### Eigenvectors from " + eigenFilestem << std::endl;
            f << "### Lattice = " << latt_size[Xdir] << "^3 x " << latt_size[Tdir] << std::endl;
            f << "### Nvec =" << nvec_list[i_nvec] ;
            f.close();
            f.clear();
        
        }
        std::cout << GridLogMessage << "###################################" << std::endl;
        std::cout << GridLogMessage << "Finished trajectory " << i_conf << std::endl;
        std::cout << GridLogMessage << "###################################" << std::endl;
    } // trajectory loop

    std::cout << GridLogMessage << "Process finished, finalizing Grid..." << std::endl;

    Grid_finalize();
}



 