/*
minimal compilation : clang++ ../Smear_radius_inefficient.cc -o inefficientsmearradius -I$HOME/phd/Grid/paboyle/include -I/usr/local//include -I/Users/s2000761/phd/include -llime -L/usr/local/lib -lGrid -L/Users/s2000761/phd/Grid/paboyle/lib -lHadrons -fopenmp -L/usr/local/opt/llvm/lib -lcrypto -L/usr/local/opt/openssl/lib -lz
*/

#include <Grid/Grid.h>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Modules/MDistil/Distil.hpp>

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
    // bool doMock = true;

    // core
    Grid_init(&argc,&argv);
    Coordinate latt_size  = GridDefaultLatt();
    Coordinate simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    GridCartesian   grid(latt_size,simd_layout,mpi_layout);
    std::cout << GridLogMessage << "Lattice time extension = " << latt_size[Tdir] << std::endl;
    const int vol = grid.gSites();
    

    // temporaries
    using LapEvecs = Grid::Hadrons::EigenPack<LatticeColourVector>;
    LapEvecs eigen4d(1 , &grid);
    LatticeColourMatrix outerPsi_r(&grid);
    
    // file handling
    std::string eigenFilepath;
    std::string eigenFilestem = "/tessfs1/work/dp008/dp008/dc-lach1/lapevec_new/result/eigen_rho0.200000_n12.1100.bin";
    // std::string eigenFilestem = "/tessfs1/work/dp008/dp008/dc-lach1/lapevecs/evecs_C1_" + stout_strength;
    // std::string eigenFilestem = "../../eigenmock";
    std::ofstream f;


    // parameters
    std::vector<int> nvec_list = {1,4,8,12,16,24,32,48,64,88,112,128,152,200,260,310,380,450,490};
    // std::vector<int> nvec_list = {1,2,4,6,8,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,110,120};
    // std::vector<int> nvec_list = {1};
    // unsigned int TRAJECTORY_START = 1100;
    // unsigned int TRAJECTORY_END = 1481;
    // unsigned int TRAJECTORY_STEP = 20;
    // unsigned int TRAJECTORY_START = 3140;
    // unsigned int TRAJECTORY_END = 3901;
    // unsigned int TRAJECTORY_STEP = 40;
    unsigned int TRAJECTORY_START = 1100;
    unsigned int TRAJECTORY_END = 1460;
    unsigned int TRAJECTORY_STEP = 40;
    // if(doMock){
    //     TRAJECTORY_START = 2000;
    //     TRAJECTORY_END = 2000;
    //     TRAJECTORY_STEP = 1;
    //     eigenFilestem = "../../evec/eigen_test";
    //     nvec_list = {1,5};
    // }
    
    std::vector<std::vector<RealD>> Psi_r(latt_size[Xdir]/2 + 1 , std::vector<RealD>(nvec_list.size(), 0.0));    //final Psi
    int i_nvec;
    
    for(unsigned int i_conf=TRAJECTORY_START ; i_conf<TRAJECTORY_END ; i_conf+=TRAJECTORY_STEP){
        std::cout << GridLogMessage << "Starting measurement on trajectory " << i_conf << "..." << std::endl;

        // reseting Psi_r
        for(int r_shift=0 ; r_shift<latt_size[Xdir]/2 + 1; r_shift++){
            for(int i=0 ; i < nvec_list.size() ; i++){
                Psi_r[r_shift][i] = 0.0;
            }
        }

        //reading eigenvectors
        eigenFilepath = eigenFilestem + "." + std::to_string(i_conf);
        std::cout << GridLogMessage << "Reading Nvec=" << nvec_list.back() << " eigenvectors from " << eigenFilepath << std::endl;
        eigen4d.resize(nvec_list.back() , &grid);
        eigen4d.read(eigenFilepath , false);

        //measuring psi
        std::cout << GridLogMessage << "Measuring spatial distribution of distil operator, Nvec =" << nvec_list << "..." << std::endl;
        for(int r_shift=0 ; r_shift<=latt_size[Xdir]/2 ; r_shift++){   //assuming spatial lattice is symmetric and pbc
            for(int spatial_dir=0 ; spatial_dir<3 ; spatial_dir++){     //assuming pbc
                outerPsi_r = Zero();
                i_nvec = 0;
                for(int i=0 ; i<nvec_list.back() ; i++){
                    outerPsi_r += outerProduct(eigen4d.evec[i] , Cshift(eigen4d.evec[i] , spatial_dir , r_shift));  // v_x X v_x+r
                    if(std::find(nvec_list.begin(), nvec_list.end() , i+1) != nvec_list.end()){
                        Psi_r[r_shift][i_nvec] += toReal(TensorRemove(sum(sqrt(trace(outerPsi_r * adj(outerPsi_r))))))/(3.0*vol*(i+1)*sqrt(Nc));
                        cout << i+1 << endl;
                        cout << toReal(TensorRemove(sum(sqrt(trace(outerPsi_r * adj(outerPsi_r))))))/(3.0*vol*(i+1)*sqrt(Nc)) << endl;
                        std::cout<<GridLogMessage<< "Summing out Psi (trajectory " << i_conf << ", Nvec=" << i+1 << ", r_shift=" << r_shift << ", spatial_dir=" << spatial_dir << "):" << std::endl;
                        // cin.get();
                        i_nvec++;
                    }
                }
            }
        }

        // result dump
        std::cout<<GridLogMessage<< "###### DUMP ######## (trajectory " << i_conf << "):" << std::endl;
        for(i_nvec=0; i_nvec < nvec_list.size() ; i_nvec++){
            std::cout<<GridLogMessage<< "Nvec=" << nvec_list[i_nvec] << std::endl;
            for(int r_shift=0 ; r_shift<latt_size[Xdir]/2 + 1; r_shift++){
                std::cout<<GridLogMessage<< "(r , Psi_r) = (" << r_shift << " , " << Psi_r[r_shift][i_nvec] << ")" << std::endl;
            }
        }

        // writing out
        for(i_nvec=0 ; i_nvec<nvec_list.size() ; i_nvec++){
            // output files
            // outFilepath = "../../data/C0/" + stout_strength;
            // string outFilepath = "../../data/C1/" + stout_strength;
            string outFilepath = "./";
            // if(doMock)
            //     outFilepath = "../../data/test/";
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
    
        // Psi_r.clear();
        std::cout << GridLogMessage << "Finished" << std::endl;
        std::cout << GridLogMessage << "################" << std::endl;
    }

    std::cout << GridLogMessage << "Process finished, finalizing Grid..." << std::endl;

    Grid_finalize();
}



 