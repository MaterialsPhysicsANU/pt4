#include "primative.h"
#include "vec.h"

#include <cmath>

std::tuple<size_t, expression_variables> current_domain(const lazy_primative& p, double t){
    //search for current domain
    size_t domain_index = 0; 
    double len_time = p.domains_h->domains[domain_index].len;
    double end_time = len_time;
    while(t > end_time){
        domain_index++;
        len_time = p.domains_h->domains[domain_index].len;
        end_time += len_time;
    }
    double dt = t - end_time + len_time;
    expression_variables ev = {t,dt};
    return std::make_tuple(domain_index, ev);
}


std::unique_ptr<primative> primative_from_lazy(const lazy_primative& p, double t){
    auto [domain_index, ev] = current_domain(p,t);
    lazy_primitive_base_parameters l_prim = p.domains_h->domains[domain_index].val.l_prim;

    std::unique_ptr<primative> prim_p;
    switch (p.type)
    {
    case (ELIPSOID):
        prim_p = std::make_unique<elipsoid>();
        break;
    case (CYLINDER):
        prim_p = std::make_unique<cylinder>();
        break; 
    case (CUBEOID):
        prim_p = std::make_unique<cubeoid>();
        break;
    }
    prim_p->loc = location(l_prim.loc,ev);
    prim_p->atn = attenuation(l_prim.atn, ev);
    return prim_p;
}

int elipsoid::inside(const std::array<double,3>& p) const{
    return (mag_sq(p) <= 1.0);
}

int cylinder::inside(const std::array<double,3>& p) const{
    auto q = abs(p);
    return (q[0]*q[0]+q[1]*q[1] < 1.0) * (std::abs(q[2]) <= 1.0);
}

int cubeoid::inside(const std::array<double,3>& p) const{
    auto q = abs(p);
    return ((*std::max_element(q.begin(), q.end())) <= 1.0);
}
