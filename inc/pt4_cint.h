#ifndef PT4_CINT_H
#define PT4_CINT_H

#include "string_handle.h"
#include "scan_spec.h"


typedef string_handle string_handle_triplet[3];

typedef enum{
	ELIPSOID = 0,
	CYLINDER,
	CUBEOID,
} primative_type;

typedef enum{
	INIT_BLEND = 0,
	REPLACE,
	ADD,
	MULTIPLY,
	MASK,
} blend_type;

typedef enum{
	INIT_FILL = 0,
	FILL_SOLID,     //solid fill
	FILL_PLANES,    //alternating black white planes
	FILL_LINES,     //checkerboard of lines
	FILL_POINTS,    //3D checkerboard
	FILL_NOISE,     //gausian noise
} fill_type;

typedef struct{
	int init;
	string_handle_triplet pos;
	string_handle_triplet sma;
	string_handle_triplet axis;
	string_handle angle;
} lazy_location;

typedef struct{
	int init;
	blend_type blend;
	fill_type fill;
	string_handle atnf;
	lazy_location loc;
} lazy_attenuation;

typedef struct{
	int init;
	lazy_location loc;
	lazy_attenuation atn;
} lazy_primitive_base_parameters;

typedef union{
	int temp;
	// lazy_elipsoid_parameters loc;
	// lazy_cubeoid_parameters atn;
} lazy_primitive_derived_parameters;

typedef struct {
	lazy_primitive_base_parameters l_prim;
	lazy_primitive_derived_parameters l_dir;
} lazy_primitive_parameters;

typedef struct{
	double len;
	lazy_primitive_parameters val;
} domain;

#ifdef __cplusplus
extern "C" {
#endif

typedef struct domain_list * domain_list_handle;
domain_list_handle create_domain_list();
void    free_domain_list(domain_list_handle);
void    push_back_domain_list(domain_list_handle, domain);

#ifdef __cplusplus
}
#endif

typedef struct {
	primative_type type;
	domain_list_handle domains_h;
} lazy_primative;

typedef struct {
	int major;
	int minor;
	int patch;
} version;

#ifdef __cplusplus
extern "C" {
#endif

typedef struct pt4 * pt4_handle;
pt4_handle create_pt4();
void    free_pt4(pt4_handle);
void    push_back_pt4(pt4_handle, lazy_primative);
void	assert_version_pt4(pt4_handle);

void	set_version_pt4(pt4_handle, version);

void	set_size_x_pt4(pt4_handle, size_t);
void	set_size_y_pt4(pt4_handle, size_t);
void	set_size_z_pt4(pt4_handle, size_t);
void	set_unit_time_per_volume_pt4(pt4_handle, double);

void	set_projection_integrand_pt4(pt4_handle, proj_type);
void	set_noise_seed_pt4(pt4_handle, int);
void	set_photon_flux_pt4(pt4_handle, double);
void	set_noise_quantisation_pt4(pt4_handle, int);
void	set_noise_poisson_pt4(pt4_handle, int);
void	set_noise_gaussian_pt4(pt4_handle, double);
void	set_unit_time_per_volume_pt4(pt4_handle, double);
void	set_revolutions_per_unit_time_pt4(pt4_handle, double);
void	set_projections_per_revolution_pt4(pt4_handle, double);
void	set_projection_supersampling_ratio_pt4(pt4_handle, double);

#ifdef __cplusplus
}
#endif


#endif