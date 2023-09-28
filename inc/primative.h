#ifndef PRIMATIVE_H
#define PRIMATIVE_H

#include "pt4.h"
#include <memory>

struct primative{
    location loc;
	attenuation atn;

    virtual int inside(const std::array<double,3>&) const =0;

    virtual ~primative() = default;
};

std::unique_ptr<primative> primative_from_lazy(const lazy_primative& p, double t);

struct elipsoid : public primative{
    int inside(const std::array<double,3>&) const override;
};

struct cylinder : public primative{
    int inside(const std::array<double,3>&) const override;
};

struct cubeoid : public primative{
    int inside(const std::array<double,3>&) const override;
};



#endif