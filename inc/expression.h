#ifndef EXPRESSION_H
#define EXPRESSION_H

#include<string>
#include<memory>

#include "expression_cint.h"

struct expression_variables{
    double t;
    double dt;
    double s;
    double x;
    double y;
    double z;
};

namespace exprtk{
    template <typename T> class symbol_table; 
    template <typename T> class expression; 
    template <typename T> class parser; 
}
typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double>   expression_t;
typedef exprtk::parser<double>       parser_t;

class expression_tp{
    std::shared_ptr<expression_t> etp;
    public:
        double value();
        
        expression_tp(void) = default;
        expression_tp(std::string& expr_str, expression_variables* ev);
};

double evaluate_string(std::string& expr_str, expression_variables& ev);

#endif