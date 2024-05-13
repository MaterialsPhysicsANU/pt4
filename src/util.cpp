
#include "util.h"

#include "simplecpp.h"

#include <ostream>
#include<tuple>
#include<vector>
#include <numbers>

double pi = std::numbers::pi_v<double>;
logger out = logger(INF, std::cout);
auto nullstream = std::ostream(nullptr);

std::string operator/(const std::string& a, const std::string& b){
    const std::string sep = "/";
    return a + sep + b;
}

logger::logger(print_mode, std::ostream& output){
  printer = &output;
}

std::ostream& logger::log(print_mode mode){
  (*printer).clear(mode >= level
  ? std::ios_base::goodbit
  : std::ios_base::badbit);
  return *printer;
}

int preprocess(const char* fpi, std::FILE* fo){
    simplecpp::DUI dui;
    bool quiet = false;
    bool error_only = false;
    // Perform preprocessing
    simplecpp::OutputList outputList;
    std::vector<std::string> files;
    simplecpp::TokenList *rawtokens;

    auto infile = std::ifstream(fpi);
    if (!infile.is_open()) {
        std::cout << "error: could not open file '" << fpi << "'" << std::endl;
        return 1;
    }
    rawtokens = new simplecpp::TokenList(infile, files,fpi,&outputList);
    

    rawtokens->removeComments();

    std::map<std::string, simplecpp::TokenList*> included = simplecpp::load(*rawtokens, files, dui, &outputList);
    for (std::pair<std::string, simplecpp::TokenList *> i : included)
        i.second->removeComments();
    simplecpp::TokenList outputTokens(files);
    simplecpp::preprocess(outputTokens, *rawtokens, files, included, dui, &outputList);
    delete rawtokens;
    rawtokens = nullptr;

    bool has_error = false;
    for (const simplecpp::Output &output : outputList) {
        print_mode level = WRN;
        std::string output_type_str;
        switch (output.type) {
        case simplecpp::Output::ERROR:
            level = ERR;
            output_type_str = "#error: ";
            break;
        case simplecpp::Output::WARNING:
            level = WRN;
            output_type_str = "#warning: ";
            break;
        case simplecpp::Output::MISSING_HEADER:
            level = ERR;
            output_type_str = "missing header: ";
            break;
        case simplecpp::Output::INCLUDE_NESTED_TOO_DEEPLY:
            level = ERR;
            output_type_str = "include nested too deeply: ";
            break;
        case simplecpp::Output::SYNTAX_ERROR:
            level = ERR;
            output_type_str = "syntax error: ";
            break;
        case simplecpp::Output::PORTABILITY_BACKSLASH:
            level = WRN;
            output_type_str = "portability: ";
            break;
        case simplecpp::Output::UNHANDLED_CHAR_ERROR:
            level = ERR;
            output_type_str = "unhandled char error: ";
            break;
        case simplecpp::Output::EXPLICIT_INCLUDE_NOT_FOUND:
            level = ERR;
            output_type_str = "explicit include not found: ";
            break;
        }
        if(level == ERR){has_error = true;}

        out.log(level) << output.location.file() << ':' << output.location.line << ": ";
        out.log(level) << output_type_str;
        out.log(level) << output.msg << std::endl;
    }
    if(has_error){return 1;}

    std::string output_string = outputTokens.stringify();
    for (auto it = output_string.begin(); it != output_string.end(); ++it) {
        if ((*it == '\\' || *it == '.') && std::next(it) != output_string.end() && *(std::next(it)) == ' ') {
            it = output_string.erase(std::next(it));
        }
    }

    fputs(output_string.c_str(), fo);
    rewind(fo);

    simplecpp::cleanup(included);
    return 0;
}