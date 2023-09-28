#include "scan_spec.h"

int uninitalised_proj_type(){
        out.log(ERR) << "fatal: proj_type not initialised" << std::endl;
        return 2;
}