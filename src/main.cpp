#include<chrono>
#include<filesystem>

#include <omp.h>

#include "util.h"
#include "cmd_actions.h"
#include "pt4.h"

std::chrono::time_point<std::chrono::high_resolution_clock> tik,tok;

void print_usage(std::ostream& err){
    err << "usage: pt4 [options] [.pt4 file]" << std::endl;
    err << "options:" << std::endl;
    err << "-v              generate volumes from .pt4" << std::endl;
    err << "-p              generate projections from .pt4" << std::endl;
    err << "-i              ignore out of date projection date from .pt4" << std::endl;
}

struct arg_switches_t{
    bool do_pt4_volumes;
    bool do_pt4_projections;
};

int main(int argc, char* argv[]) {
    out.log(INF) << "================ pt4 ===============" << std::endl;

#if _OPENMP
    omp_set_max_active_levels(8);
#endif
    std::string pt4_name;
    arg_switches_t arg_switches = {false,false};

    if (!(argc > 1)) {
        out.log(ERR) << "fatal error: no input files" << std::endl;
        print_usage(out.log(ERR));
        exit(1);
    }
    for(size_t argi = 1; argi < argc; ++argi){
        std::string arg = argv[argi];
        if(arg == "-v"){arg_switches.do_pt4_volumes = true; continue;}
        if(arg == "-p"){arg_switches.do_pt4_projections = true; continue;}
        if(arg[0] == '-'){
            out.log(ERR) << "invalid argument: " << arg << std::endl;
            print_usage(out.log(ERR));
            exit(1);
        }
        if(arg[0] != '-'){pt4_name = argv[argi]; continue;}//if multiple files are called then they get replaced
        else{
            out.log(ERR) << "fatal error: no input files" << std::endl;
            exit(1);
        }
    }
    if(     !arg_switches.do_pt4_volumes 
        &&  !arg_switches.do_pt4_projections 
    ){
        out.log(ERR) << "fatal error: no function flags" << std::endl;
        print_usage(out.log(ERR));
        exit(1);
    }

    std::string pt4_fp = pt4_name + ".pt4";//input file
    std::string ncf_dir = pt4_name + '/';//output directory

    //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //

    if(arg_switches.do_pt4_volumes || arg_switches.do_pt4_projections){
        pt4 pt4_0 = load_pt4(pt4_fp);

        if(!std::filesystem::exists(pt4_fp)){
            out.log(ERR) << "fatal: file does not exists: " << pt4_fp << std::endl;
            exit(2);
        }else{
            //remove old files
            std::filesystem::remove_all(ncf_dir);
            //create output dir if we find and preprocess the .pt4 file
            std::filesystem::create_directories(ncf_dir);
            std::filesystem::create_directories(ncf_dir+"ztime/");
            std::filesystem::create_directories(ncf_dir+"vol/");
            std::filesystem::create_directories(ncf_dir+"proj/");
        }
        out.log(INF) << "SIMULATING PHANTOM..." << std::endl;
        tik = std::chrono::high_resolution_clock::now();

        if(arg_switches.do_pt4_volumes){
            out.log(INF) << "GENERATING VOLUMES..." << std::endl;
            write_volumes(pt4_name, pt4_0, true);
        }
        if(arg_switches.do_pt4_projections){
            out.log(INF) << "GENERATING PROJECTIONS..." << std::endl;
            write_projections(pt4_name, pt4_0);
        }
        
        tok = std::chrono::high_resolution_clock::now();
        out.log(INF) << "SIMULATION TIME: " << std::chrono::duration_cast<std::chrono::milliseconds>(tok - tik).count() << " ms" << std::endl;
    }

    return 0;
}