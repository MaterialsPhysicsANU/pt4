#include "pt4.h"

#include "util.h"
#include "quaternion.h"

#include <numbers>
#include <cstring>

#include "quick_rng.h"
#include "expression.h"


int uninitalised_fill_type(){
        out.log(ERR) << "fatal: fill_type not initialised" << std::endl;
        return 2;
}

std::array<double,3> location::inv_affine(const std::array<double,3>& p) const{
    std::array<double,3> r; 
    //unmove
    r = p - pos;
    //unrotate
    rotate(r,axis,-angle);
    //unscale
    r = r/sma;
    return r;
}

location::location(void){
    pos = {0.0,0.0,0.0};
    sma = {1.0,1.0,1.0};
	axis = {0.0,0.0,1.0};
	angle = 0.0;
}

location::location(const lazy_location& l, expression_variables ev){
    pos = {evaluate_string(*l.pos[0],ev),evaluate_string(*l.pos[1],ev),evaluate_string(*l.pos[2],ev)};
    sma = {evaluate_string(*l.sma[0],ev),evaluate_string(*l.sma[1],ev),evaluate_string(*l.sma[2],ev)};
	axis = {evaluate_string(*l.axis[0],ev),evaluate_string(*l.axis[1],ev),evaluate_string(*l.axis[2],ev)};
	angle = evaluate_string(*l.angle,ev);
}

std::ostream& operator<<(std::ostream& os, const lazy_location& obj) {
    // Define how the object should be represented in the output stream
    os << std::endl;
    os << "init: " << obj.init << std::endl;
    os << "pos: {" << *(obj.pos[0]) << "," << *(obj.pos[1]) << "," << *(obj.pos[2]) << ",}" << std::endl;
    os << "sma: {" << *(obj.sma[0]) << "," << *(obj.sma[1]) << "," << *(obj.sma[2]) << ",}" << std::endl;
    os << "axis: {" << *(obj.axis[0]) << "," << *(obj.axis[1]) << "," << *(obj.axis[2]) << ",}" << std::endl;
    os << "angle: " << *(obj.angle) << std::endl;
    return os;
}

attenuation::attenuation(const attenuation& atn){
    blend = atn.blend;
	atnf  = atn.atnf;
    fill = atn.fill;
    loc = atn.loc;
    ev = new expression_variables(*(atn.ev));
    expression = atn.expression;
}
attenuation::attenuation(attenuation && atn){
    blend = atn.blend;
	atnf  = atn.atnf;
    fill = atn.fill;
    loc = atn.loc;
    std::swap(ev,atn.ev);
    expression = atn.expression;
}
attenuation& attenuation::operator=(const attenuation & atn){
    blend = atn.blend;
	atnf  = atn.atnf;
    fill = atn.fill;
    loc = atn.loc;
    ev = new expression_variables(*(atn.ev));
    expression = atn.expression;
    return *this;
}
attenuation& attenuation::operator=(attenuation && atn){
    blend = atn.blend;
	atnf  = atn.atnf;
    fill = atn.fill;
    loc = atn.loc;
    std::swap(ev,atn.ev);
    expression = atn.expression;
    return *this;
}

attenuation::~attenuation(void){
    delete ev;
}
attenuation::attenuation(const lazy_attenuation& a, const expression_variables ev0){
    blend = a.blend;
	fill  = a.fill;                
    ev = new expression_variables(ev0);
    loc   = location(a.loc,*ev);
	atnf  = a.atnf;
    expression = expression_tp(*(a.atnf), ev);
}

float attenuation::texture_at(std::array<double,3> pos){
    
    switch (fill) {
    case INIT_FILL:
        exit(uninitalised_fill_type());
        break;
    case FILL_SOLID:
        return 1;
        break;
    case FILL_PLANES:
        return 2*(std::abs(fmod(pos[2], 1.0)) > 0.5) - 1;
        break;
    case FILL_LINES:
        return 2*(std::abs(fmod(pos[2], 1.0)) > 0.5)*(std::abs(fmod(pos[1], 1.0)) > 0.5) - 1;
        break;
    case FILL_POINTS:
        return 2*(std::abs(fmod(pos[2], 1.0)) > 0.5)*(std::abs(fmod(pos[1], 1.0)) > 0.5)*(std::abs(fmod(pos[0], 1.0)) > 0.5) - 1;
        break;
    case FILL_NOISE:
        std::array<int,3> index = floor<int,double,3>(pos);
        return texture_noise_volume[index];
        break;
    }
    return NAN;
}

float attenuation::atn_at(std::array<double,3> pos){
    std::array<double,3> tex_pos = loc.inv_affine(pos);
    ev->x = tex_pos[0];
    ev->y = tex_pos[1];
    ev->z = tex_pos[2];
    ev->s = texture_at(tex_pos);

    float result = expression.value();
    return result;
}

extern "C"
{
    string_handle create_string(char* s) {
        string_handle result = new std::string();
        for(; *s != '\0'; s++){
        if (*s != '\\') {
            *result += *s;
            }
        }
        return result;
    }
    string_handle create_quoted_string(char* s) {
        std::string stds = std::string(s);
        std::string s_unquote = stds.substr(1, stds.size() - 2);
        string_handle result = new std::string(s_unquote);
        return result;
    }
    void push_back_quoted_string(string_handle p, char* s) {
        std::string stds = std::string(s);
        std::string s_unquote = stds.substr(1, stds.size() - 2);
        p->append('\n' + s_unquote);
    }
    void    free_string(string_handle p) { delete p; }
    const char*    c_str(string_handle p) {return p->c_str(); }
}

extern "C"
{
    domain_list_handle create_domain_list() { return new domain_list(); }
    void    free_domain_list(domain_list_handle p) { delete p; }
    void    push_back_domain_list(domain_list_handle p, domain d) { return p->push_back(d); }
}

void domain_list::push_back(domain d){
    domains.push_back(d);
}
double domain_list::t_length() const{
    double t_len = 0.0;
    for(domain d : domains){
        t_len += d.len;
    }
    return t_len;
}

extern "C"
{
    pt4_handle create_pt4() { return new pt4(); }
    void    free_pt4(pt4_handle p) { delete p; }
    void    push_back_pt4(pt4_handle p, lazy_primative f) { return p->push_back(f); }

    void	set_version_pt4(pt4_handle p, version v){p->file_ver = v;}
    void	assert_version_pt4(pt4_handle p){
                                                version f = p->file_ver;
                                                version c = p->code_ver;
                                                if(f.major != c.major || f.minor != c.minor || f.patch != c.patch){
                                                    out.log(ERR) << "fatal .pt4 verions mismatch!" << std::endl;
                                                    out.log(ERR) << "    file version: " << f.major << '.' << f.minor << '.' << f.patch << std::endl;
                                                    out.log(ERR) << "    code version: " << c.major << '.' << c.minor << '.' << c.patch << std::endl;
                                                    exit(2);
                                                }
                                            }

    void	set_size_x_pt4(pt4_handle p, size_t n){p->size[0] = n;}
    void	set_size_y_pt4(pt4_handle p, size_t n){p->size[1] = n;}
    void	set_size_z_pt4(pt4_handle p, size_t n){p->size[2] = n;}

    void	set_projection_integrand_pt4            (pt4_handle p, proj_type v){p->scan_spec.proj_int = v;}
    void	set_noise_seed_pt4                      (pt4_handle p, int v){p->scan_spec.seed = v;}
    void	set_photon_flux_pt4                     (pt4_handle p, double v){p->scan_spec.I0 = v;}
    void	set_noise_quantisation_pt4               (pt4_handle p, int v){p->scan_spec.noise_quantisation = v;}
    void	set_noise_poisson_pt4                   (pt4_handle p, int v){p->scan_spec.noise_poisson = v;}
    void	set_noise_gaussian_pt4                  (pt4_handle p, double v){p->scan_spec.noise_gaussian = v;}
    void	set_unit_time_per_volume_pt4           (pt4_handle p, double v){p->vol_t_step = v;}
    void	set_revolutions_per_unit_time_pt4       (pt4_handle p, double v){p->scan_spec.revolutions_per_unit_time = v;}
    void	set_projections_per_revolution_pt4      (pt4_handle p, double v){p->scan_spec.projections_per_revolution = v;}
    void	set_projection_supersampling_ratio_pt4  (pt4_handle p, double v){p->scan_spec.projection_supersampling_ratio = v;}
}

void pt4::push_back(lazy_primative p){
       primitives.push_back(p);
}

double pt4::t_length() const{
    double max_length = 0.0;
    for(lazy_primative p : primitives){
        max_length = std::max(max_length, p.domains_h->t_length());
    }
    return max_length;
}


void fill_string_handle(string_handle& dst, string_handle& src){
    if(!dst){
        dst = src;
    }
}

void fill_string_handle_triplet(string_handle_triplet& dst, string_handle_triplet& src){
    if(!dst[0] || !dst[1] || !dst[2]){
        memcpy(dst,src,sizeof(string_handle_triplet));
    }
}

void fill_lazy_location(lazy_location& dst, lazy_location& src){
    if(!dst.init){
        fill_string_handle_triplet(dst.pos,src.pos);
        fill_string_handle_triplet(dst.sma,src.sma);
        fill_string_handle_triplet(dst.axis,src.axis);
        fill_string_handle(dst.angle,src.angle);
        dst.init = 1;
    }
}

void fill_blend_type(blend_type& dst, blend_type& src){
    if(dst == INIT_BLEND){
        dst = src;
    }
}

void fill_fill_type(fill_type& dst, fill_type& src){
    if(dst == INIT_FILL){
        dst = src;
    }
}

void fill_lazy_attenuation(lazy_attenuation& dst, lazy_attenuation& src){
    if(!dst.init){
        fill_blend_type(dst.blend,src.blend);
        fill_fill_type(dst.fill,src.fill);
        fill_string_handle(dst.atnf,src.atnf);
        fill_lazy_location(dst.loc,src.loc);
        dst.init = 1;
    }
}

extern "C" int yyparse(void);
extern "C" FILE* yyin;
extern "C" pt4_handle parser_pt4_handle;

pt4::pt4(std::FILE* f){

    yyin = f;
	if ( yyin == NULL ){
        char ebuf[1024];
        #ifdef _WIN32
	strerror_s(ebuf, sizeof(ebuf), errno);
	#else
	strerror_r(errno, ebuf, sizeof(ebuf));
	#endif
        out.log(ERR) << ebuf << std::endl;
        exit(errno);
    }
    parser_pt4_handle = this;
    if(yyparse()){exit(1);};

    assert_version_pt4(this);

    //reverse reversed lists while parsing
    for(lazy_primative& p : primitives){
        std::reverse(p.domains_h->domains.begin(), p.domains_h->domains.end());
    }
    //fill in any empty locations
    for (int j = 0; j < primitives.size(); j++) {
        lazy_primative& p = primitives[j]; 
        std::vector<domain>& domains = p.domains_h->domains;

        if(!domains[0].val.l_prim.init){
            printf("ERROR: first domain of primative %d missing parameters\n", j); exit(2);
        }
        for(size_t i = 1; i < domains.size(); i++){
            lazy_primitive_base_parameters& l_prim_init = domains[i].val.l_prim;
            lazy_primitive_base_parameters& l_prim_prev = domains[i-1].val.l_prim;

            fill_lazy_location(l_prim_init.loc,l_prim_prev.loc);
            fill_lazy_attenuation(l_prim_init.atn,l_prim_prev.atn);
            l_prim_init.init = 1;
        }
    }
}
