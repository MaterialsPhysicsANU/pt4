#ifndef PT4_H
#define PT4_H

#include<vector>
#include<array>
#include<string>
typedef std::string std__string;

#include "scan_spec.h"

#include "pt4_cint.h"
#include "expression.h"


int uninitalised_fill_type();

class domain_list{
    public:
        std::vector<domain> domains;

        void push_back(domain);
        double t_length() const;
};

class pt4{
    public:
        std::vector<lazy_primative> primitives;
        double vol_t_step;
        size_t size[3];//xyz
        scan_spec_t scan_spec;
        version file_ver;
        static inline const version code_ver = {0,5,4};
        

        void push_back(lazy_primative);
        double t_length() const;
        pt4(void) = default;
        pt4(std::FILE*);
};

struct location{
    std::array<double,3> pos;
	std::array<double,3> sma;
	std::array<double,3> axis;
	double angle;

    std::array<double,3> inv_affine(const std::array<double,3>& p) const;

    location(void);//construts the identity location
    location(const lazy_location& p, expression_variables ev);
};

std::ostream& operator<<(std::ostream& os, const lazy_location& obj);

struct attenuation{
	blend_type blend;
    fill_type fill;
    location loc;
	string_handle atnf;
    expression_variables* ev;
	expression_tp expression;

    float texture_at(std::array<double,3>);
    float atn_at(std::array<double,3>);

    attenuation(void) = default;
    attenuation(const attenuation &);                // Copy Constructor
    attenuation(attenuation &&);                     // Move Constructor
    attenuation& operator=(const attenuation &);     // Copy Assignment
    attenuation& operator=(attenuation &&);          // Move Assignment
    ~attenuation(void);
    attenuation(const lazy_attenuation& p,const expression_variables ev);
};

#endif