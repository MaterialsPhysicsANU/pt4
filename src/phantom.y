%{
#include <stdio.h>     /* C declarations used in actions */
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
int isatty(int);
int _get_osfhandle(int);

#ifdef _MSC_VER
	#include <windows.h>
#endif

#include "pt4_cint.h"
#include "expression_cint.h"

#include "y.tab.h"

void yyerror (char *s);
int yylex();
extern int yylineno;
extern FILE* yyin;

const string_handle_triplet string_handle_triplet_NULL = {NULL,NULL,NULL};

int yydebug = 1;

extern void swex(int);
extern void swpr(int);

extern enum yytokentype current_primative_tk;


pt4_handle parser_pt4_handle;
%}

%define parse.error verbose

%union {
		int num; 
		double dbl; 
		char* cstr;
		version ver;
		proj_type e_proj;
		blend_type e_blend;
		fill_type  e_fill;
		string_handle str;
		string_handle_triplet str_trp;
		lazy_primitive_base_parameters base_prmt;
		lazy_primitive_derived_parameters  dir_prmt;
		lazy_primitive_parameters prim_prmt;
		lazy_primative prim;
		lazy_location loc;
		lazy_attenuation atn;
		domain dom;
		domain_list_handle dom_lst_hand;
		}         /* Yacc definitions */


%start line

%token version_tk

%token size_x_tk
%token size_y_tk
%token size_z_tk

%token projection_integrand_tk
%token noise_seed_tk
%token photon_flux_tk
%token noise_quanisation_tk
%token noise_poisson_tk
%token noise_gaussian_tk
%token unit_time_per_volume_tk
%token revolutions_per_unit_time_tk
%token projections_per_revolution_tk
%token projection_supersampling_ratio_tk

%token elipsoid_tk
%token cylinder_tk
%token cubeoid_tk
%token loc_tk
%token atn_tk
%token blend_tk
%token fill_tk
%token atnf_tk
%token pos_tk
%token sma_tk
%token axis_tk
%token angle_tk
%token attenuation_tk
%token len_tk
%token name
%token domain_tk
%token val_tk

%token <num> num_int
%token <dbl> num_flt
%token <ver> num_ver
%token <cstr> str_expr
%token <cstr> str_expr_quoted
%token <e_proj> enum_proj
%token <e_blend> enum_blend
%token <e_fill> enum_fill

%type <num> plus_minus 
%type <dbl> sign_flt 
%type <dbl> number 
%type <str> expression
%type <str> eq_expression
%type <prim_prmt> eq_val_members
%type <base_prmt> val_members_base
%type <dir_prmt> val_members_derived
%type <str_trp> expression_trp 
%type <str_trp> eq_expression_trp

%type <str_trp> pos_option
%type <str_trp> sma_option
%type <str_trp> axis_option
%type <str> angle_option

%type <e_blend> blend_option
%type <e_fill> fill_option
%type <str> atnf_option

%type <loc> loc_members
%type <loc> loc_option
%type <atn> atn_members
%type <atn> atn_option



%type <str> str_expr_quoted_elems

%type <dom> domain_members
%type <dom> domain

%type <dom_lst_hand> domain_list_elems
%type <dom_lst_hand> domain_list
%type <prim> primitive


%%


/* descriptions of expected inputs     corresponding actions (in C) */

line    : primitive ';'			{push_back_pt4(parser_pt4_handle,$1);}
		| metadata ';'			{;}
		| line primitive ';'	{push_back_pt4(parser_pt4_handle,$2);}
		| line metadata ';'		{;}
        ;

metadata : version_tk '=' num_ver {set_version_pt4(parser_pt4_handle, $3); assert_version_pt4(parser_pt4_handle);}
		 | size_x_tk eq_expression {set_size_x_pt4(parser_pt4_handle, evaluate_conststring(c_str($2)));}
		 | size_y_tk eq_expression {set_size_y_pt4(parser_pt4_handle, evaluate_conststring(c_str($2)));}
		 | size_z_tk eq_expression {set_size_z_pt4(parser_pt4_handle, evaluate_conststring(c_str($2)));}
		 | projection_integrand_tk '=' enum_proj {set_projection_integrand_pt4(parser_pt4_handle, $3);}
		 | noise_seed_tk eq_expression {set_noise_seed_pt4(parser_pt4_handle, evaluate_conststring(c_str($2)));}
		 | photon_flux_tk eq_expression {set_photon_flux_pt4(parser_pt4_handle, evaluate_conststring(c_str($2)));}
		 | noise_quanisation_tk '=' num_int {set_noise_quanisation_pt4(parser_pt4_handle, $3);}
		 | noise_poisson_tk '=' num_int {set_noise_poisson_pt4(parser_pt4_handle, $3);}
		 | noise_gaussian_tk eq_expression {set_noise_gaussian_pt4(parser_pt4_handle, evaluate_conststring(c_str($2)));}
		 | unit_time_per_volume_tk eq_expression {set_unit_time_per_volume_pt4(parser_pt4_handle, evaluate_conststring(c_str($2)));}
		 | revolutions_per_unit_time_tk eq_expression {set_revolutions_per_unit_time_pt4(parser_pt4_handle, evaluate_conststring(c_str($2)));}
		 | projections_per_revolution_tk eq_expression {set_projections_per_revolution_pt4(parser_pt4_handle, evaluate_conststring(c_str($2)));}
		 | projection_supersampling_ratio_tk eq_expression {set_projection_supersampling_ratio_pt4(parser_pt4_handle, evaluate_conststring(c_str($2)));}
		 ;


primitive : elipsoid_tk {current_primative_tk = elipsoid_tk;} name '=' domain_list	{$$ = (lazy_primative){.type =  ELIPSOID, .domains_h = $5,};}
		  | cylinder_tk  {current_primative_tk = cylinder_tk;} name '=' domain_list	{$$ = (lazy_primative){.type =  CYLINDER, .domains_h = $5,};}
		  | cubeoid_tk  {current_primative_tk = cubeoid_tk;} name '=' domain_list	{$$ = (lazy_primative){.type =  CUBEOID, .domains_h = $5,};}
		  ;

domain_list : '{' domain_list_elems '}'		{$$ = $2;}
			;

domain_list_elems	: domain	comma_terminator	{$$ = create_domain_list(); push_back_domain_list($$,$1); }
					| domain ',' domain_list_elems	{$$ = $3; push_back_domain_list($3,$1);} 
					;			

domain : domain_tk name '=' domain_members			{$$ = $4;}
	   ;

domain_members : '{'
					len_tk eq_expression ','
					{swpr(1);} val_tk eq_val_members comma_terminator
				 '}'
				{$$ = (domain){
						.len = evaluate_conststring(c_str($3)),
						.val = $7,
				};}
				;

/* this should be a union type eventually */
eq_val_members : elipsoid_tk {swpr(0);} val_members_base {$$ = (lazy_primitive_parameters){.l_prim = $3};}
			   | cylinder_tk {swpr(0);} val_members_base  {$$ = (lazy_primitive_parameters){.l_prim = $3};}
			   | cubeoid_tk {swpr(0);} val_members_base  {$$ = (lazy_primitive_parameters){.l_prim = $3};}
			   ;


loc_members : '{' 	pos_option 
					sma_option
					axis_option
					angle_option 
 
					{$$ = (lazy_location){
						.init = $2[0] && $2[1] && $2[2] && $3[0] && $3[1] && $3[2] && $4[0] && $4[1] && $4[2] && $5,
						.pos =  {$2[0],$2[1],$2[2]},
						.sma =  {$3[0],$3[1],$3[2]},
						.axis = {$4[0],$4[1],$4[2]},
						.angle = $5,
					};}
				;

pos_option: pos_tk eq_expression_trp option_end	{memcpy($$,$2,sizeof(string_handle_triplet));}
		  | 									{memcpy($$,string_handle_triplet_NULL,sizeof(string_handle_triplet));}
		  ;

sma_option: sma_tk eq_expression_trp option_end	{memcpy($$,$2,sizeof(string_handle_triplet));}
		  | 									{memcpy($$,string_handle_triplet_NULL,sizeof(string_handle_triplet));}
		  ;

axis_option: axis_tk eq_expression_trp option_end	{memcpy($$,$2,sizeof(string_handle_triplet));}
		  | 										{memcpy($$,string_handle_triplet_NULL,sizeof(string_handle_triplet));}
		  ;

angle_option: angle_tk eq_expression option_end	{$$ = $2;}
		  | 									{$$ = NULL;}
		  ;

atn_members : '{' 	blend_option 
					fill_option
					atnf_option 
					loc_option 
					{$$ = (lazy_attenuation){
						.init = ($2 != INIT_BLEND) && ($3 != INIT_FILL) && ($4 != NULL) && ($5.init),
						.blend =  $2,
						.fill =  $3,
						.atnf = $4,
						.loc = $5,
					};}
				;

blend_option: blend_tk '=' enum_blend option_end	{$$ = $3;}
		  | 										{$$ = INIT_BLEND;}
		  ;

fill_option: fill_tk '=' enum_fill option_end		{$$ = $3;}
		  | 										{$$ = INIT_FILL;}
		  ;

atnf_option: atnf_tk eq_expression option_end		{$$ = $2;}
		  | 										{$$ = NULL;}
		  ;


val_members_base : '{'  loc_option 
					    atn_option 
 
					{$$ = (lazy_primitive_base_parameters){
						.init = $2.init && $3.init,
						.loc =  $2,
						.atn =  $3,
					};}
				;

loc_option: loc_tk '=' loc_members option_end	{$$ = $3;}
		  | 									{$$ = (lazy_location){.init = 0};}
		  ;

atn_option: atn_tk '=' atn_members option_end	{$$ = $3;}
		  | 									{$$ = (lazy_attenuation){.init = 0};}
		  ;

eq_expression_trp : {swex(1);} '='  expression_trp	{swex(0); memcpy($$,$3,sizeof(string_handle_triplet));}
			;

eq_expression		: {swex(1);} '='  expression	{swex(0); $$ = $3;}
			;

expression_trp : '{' expression ',' expression ',' expression comma_terminator '}' {$$[0] = $2; $$[1] = $4; $$[2] = $6;}
		 ;

expression : str_expr					{$$ = create_string($1);}
		   | str_expr_quoted_elems		{$$ = $1;}
	 	   ;

str_expr_quoted_elems	: str_expr_quoted	{$$ = create_quoted_string($1); }
						| str_expr_quoted_elems str_expr_quoted  {$$ = $1; push_back_quoted_string($1,$2);} 
						;			

sign_flt : number								{$$ = $1;}
		 | plus_minus number					{$$ = $1 * $2;}
		 ; 

number	 : num_flt								{$$ = $1;}
		 | num_int								{$$ = $1;}
		 ;

plus_minus : '+' 			{$$ = 1;}
		   | '-'			{$$ = -1;}
		   ;

option_end		 : ',' 			{;}
				 | ',' '}'		{;}
				 ; 

comma_terminator : ',' 		{;}
				 |			{;}
				 ;

%%                     /* C code */

void yyerror (char *s) {
		const char* errorfp = "e.pt4";
		FILE* pt4f = yyin;
		FILE* of = fopen(errorfp, "w");
		rewind(pt4f);
		char ch = fgetc(pt4f);
		while (ch != EOF)
		{
			fputc(ch, of);
			ch = fgetc(pt4f);
		}
		fclose(of);

		fprintf (stderr, "%s(%d): %s\n", errorfp,yylineno, s); 
}

