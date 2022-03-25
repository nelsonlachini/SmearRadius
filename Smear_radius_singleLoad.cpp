/*
minimal compilation : clang++ ../Smear_radius_inefficient.cc -o inefficientsmearradius -I$HOME/phd/Grid/paboyle/include -I/usr/local//include -I/Users/s2000761/phd/include -llime -L/usr/local/lib -lGrid -L/Users/s2000761/phd/Grid/paboyle/lib -lHadrons -fopenmp -L/usr/local/opt/llvm/lib -lcrypto -L/usr/local/opt/openssl/lib -lz
*/

#include <Grid/Grid.h>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/NamedTensor.hpp>

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
    // bool doMock = true;

    //if (argc != 4)
    //{
        std::cerr << "usage: " << argv[0] << " <evec stem> <output stem> <max nvec> <config list>";
        //std::cerr << std::endl;
        
        //return EXIT_FAILURE;
    //}
    string eigenFilestem = argv[1];
    string outStem = argv[2];
    int nVec_max = stoi(argv[3]);
    vector<int> config_list = strToVec<int>(argv[4]);
    const int step = 20;

    // core
    Grid_init(&argc,&argv);
    Coordinate latt_size  = GridDefaultLatt();
    Coordinate simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    GridCartesian   grid(latt_size,simd_layout,mpi_layout);
    std::cout << GridLogMessage << "Lattice time extension = " << latt_size[Tdir] << std::endl;
    std::cout << GridLogMessage << "Lattice spatial extension = " << latt_size[Xdir] << std::endl; 
    const int vol = grid.gSites();
    
    // temporaries
    using LapEvecs = Grid::Hadrons::EigenPack<LatticeColourVector>;
    LapEvecs eigen4d(1 , &grid);
    LatticeColourMatrix outerPsi_r(&grid);
    
    // file handling
    std::string eigenFilepath;
    //std::string eigenFilestem = "/tessfs1/work/dp008/dp008/dc-lach1/lapevec_new/result/eigen_rho0.200000_n12"; //C0
    //std::string eigenFilestem = "t96l32_b3.80_kls0.138963_csw1.95524/evec_t96l32_b3.80_kls0.138963_csw1.95524";
    std::ofstream f;
    //std::string outStem = "t96l32/";

    // parameters
    std::vector<int> nvec_list;// = {1,4,8,12,16,24,32,48,64,88,112,128,152,200,260,310,380,450,500};
    for(int i=step;i<nVec_max;i+=step) nvec_list.push_back(i);
    // std::vector<int> nvec_list = {1,3};
    // unsigned int TRAJECTORY_START = 1100;
    // unsigned int TRAJECTORY_END = 1481;
    // unsigned int TRAJECTORY_STEP = 20;
    // unsigned int TRAJECTORY_START = 3140;
    // unsigned int TRAJECTORY_END = 3901;
    // unsigned int TRAJECTORY_STEP = 40;
    unsigned int TRAJECTORY_START = 80;
    unsigned int TRAJECTORY_END = 110;
    unsigned int TRAJECTORY_STEP = 10;
    
    std::vector<std::vector<RealD>> Psi_r(latt_size[Xdir]/2 + 1 , std::vector<RealD>(nvec_list.size(), 0.0));    //final Psi
    int i_nvec;
    
    //for(unsigned int i_conf=TRAJECTORY_START ; i_conf<TRAJECTORY_END ; i_conf+=TRAJECTORY_STEP){
    for(auto i_conf : config_list){  
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
                    outerPsi_r += outerProduct( eigen4d.evec[i] , adj(Cshift(eigen4d.evec[i] , spatial_dir , r_shift)) ) ;  // v_x X v_x+r
                    // cout << outerPsi_r << endl;
                    // cin.get();
                    if(std::find(nvec_list.begin(), nvec_list.end() , i+1) != nvec_list.end()){
                        Psi_r[r_shift][i_nvec] += toReal(TensorRemove( sum(sqrt(trace(outerPsi_r * adj(outerPsi_r)))) ))/(3.0*vol*(i+1)*sqrt(Nc));
                        // cout << sum(trace( outerPsi_r * adj(outerPsi_r) )) << endl;
                        // cin.get();
                        // cout << i+1 << endl;
                        // cout << toReal(TensorRemove(sum(sqrt(trace(outerPsi_r * adj(outerPsi_r))))))/(3.0*vol*(i+1)*sqrt(Nc)) << endl;
                        std::cout<<GridLogMessage<< "Summing out Psi (trajectory " << i_conf << ", Nvec=" << i+1 << ", r_shift=" << r_shift << ", spatial_dir=" << spatial_dir << "):" << std::endl;
                        // cin.get();s
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
            string outFilepath = outStem;
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



 
