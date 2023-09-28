
#include <numbers>


#include "exprtk.hpp"
#pragma message("INCLUDE OF 'exprtk.hpp'. COMPILE TIMES WILL BE LONG")

#include "expression.h"
#include "util.h"

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double> expression_t;
typedef exprtk::parser<double>       parser_t;

double expression_tp::value(){return etp->value();}

expression_tp::expression_tp(std::string& expr_str, expression_variables* ev){
    etp = std::make_unique<expression_t>();
    //compile expression
    symbol_table_t symbol_table;
    symbol_table.add_variable("t",ev->t);    //global time
    symbol_table.add_variable("dt",ev->dt);  //local time
    symbol_table.add_variable("s",ev->s);   //texture value
    symbol_table.add_variable("x",ev->x);   //local x
    symbol_table.add_variable("y",ev->y);   //local y
    symbol_table.add_variable("z",ev->z);   //local z

    double pi_workaround = std::numbers::pi;
    double e_workaround = std::numbers::e;
    double ln2_workaround = std::numbers::ln2;
    double sqrt2_workaround = std::numbers::sqrt2;
    symbol_table.add_variable("pi",pi_workaround);
    symbol_table.add_variable("e",e_workaround);
    symbol_table.add_variable("ln2",ln2_workaround);
    symbol_table.add_variable("sqrt2",sqrt2_workaround);
    
    etp->register_symbol_table(symbol_table);

    parser_t parser;

    if (!parser.compile(expr_str,*etp)){
        out.log(ERR) << "COMPILE ERROR:~~~~~~~~~~~~~" << std::endl << expr_str << std::endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    }
}

double evaluate_string(std::string& expr_str, expression_variables& ev){

    symbol_table_t symbol_table;
    symbol_table.add_variable("t",ev.t);    //global time
    symbol_table.add_variable("dt",ev.dt);  //local time

    double pi_workaround = std::numbers::pi;
    double e_workaround = std::numbers::e;
    double ln2_workaround = std::numbers::ln2;
    double sqrt2_workaround = std::numbers::sqrt2;
    symbol_table.add_variable("pi",pi_workaround);
    symbol_table.add_variable("e",e_workaround);
    symbol_table.add_variable("ln2",ln2_workaround);
    symbol_table.add_variable("sqrt2",sqrt2_workaround);
     
    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;

    if (!parser.compile(expr_str,expression)){
        out.log(ERR) << "COMPILE ERROR:~~~~~~~~~~~~~" << std::endl << expr_str << std::endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        return NAN;
    }

    double result = expression.value();
    return result;
}


double evaluate_conststring(const char* str){

    symbol_table_t symbol_table;

    double pi_workaround = std::numbers::pi;
    double e_workaround = std::numbers::e;
    double ln2_workaround = std::numbers::ln2;
    double sqrt2_workaround = std::numbers::sqrt2;
    symbol_table.add_variable("pi",pi_workaround);
    symbol_table.add_variable("e",e_workaround);
    symbol_table.add_variable("ln2",ln2_workaround);
    symbol_table.add_variable("sqrt2",sqrt2_workaround);
     
    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;

    if (!parser.compile(std::string(str),expression)){
        out.log(ERR) << "COMPILE ERROR:~~~~~~~~~~~~~" << std::endl << str << std::endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        return NAN;
    }

    double result = expression.value();
    return result;
}