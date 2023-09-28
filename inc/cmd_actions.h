#ifndef CMD_ACTIONS_H
#define CMD_ACTIONS_H

#include <string>
#include <optional>

#include "util.h"
#include "pt4.h"
#include "ndvec.h"

pt4 load_pt4(std::string& pt4_fp);
int write_volumes(const std::string& odir, const pt4 pt4_0, bool zproject);
int write_projections(const std::string& pt4_name, const pt4 pt4_0);

#endif