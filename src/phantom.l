%{
#include "pt4_cint.h"
#include "y.tab.h"
int isatty(int);

void yyerror (char *s);
int yylex();
enum yytokentype current_primative_tk;
version atover(const char*);
%}


%option nounistd
%option yylineno
%option stack

%s PRIM
%s EXPR
%%

\/\/.*				                     {}
#.*				                  	     {}
\"([^\\\"]|\\.)*\"                       {yylval.cstr = yytext; return str_expr_quoted;}
<EXPR>[-+*\/^()._a-zA-Z0-9]([-+*\/^()._a-zA-Z0-9 \t\r\n]|\\,)*               {yylval.cstr = yytext; return str_expr;}
<PRIM>"="								 {return current_primative_tk;}
"version"                                {return version_tk;}
"size_x"                                 {return size_x_tk;}
"size_y"                                 {return size_y_tk;}
"size_z"                                 {return size_z_tk;}
"projection_integrand"                   {return projection_integrand_tk;}
"INTENSITY"							 	 {yylval.e_proj = INTENSITY; return enum_proj;}
"ATTENUATION"							 {yylval.e_proj = ATTENUATION; return enum_proj;}
"noise_seed"                             {return noise_seed_tk;}
"photon_flux"                            {return photon_flux_tk;}
"noise_quantisation"                      {return noise_quantisation_tk;}
"noise_poisson"                          {return noise_poisson_tk;}
"noise_gaussian"                         {return noise_gaussian_tk;}
"unit_time_per_volume"                  {return unit_time_per_volume_tk;}
"revolutions_per_unit_time"              {return revolutions_per_unit_time_tk;}
"projections_per_revolution"             {return projections_per_revolution_tk;}
"projection_supersampling_ratio"         {return projection_supersampling_ratio_tk;}
"ellipsoid"				                 {return ellipsoid_tk;}
"cylinder"				                 {return cylinder_tk;}
"cubeoid"				                 {return cubeoid_tk;}
".loc"									 {return loc_tk;}
".atn"									 {return atn_tk;}
".blend"								 {return blend_tk;}
"REPLACE"								 {yylval.e_blend = REPLACE; return enum_blend;}
"ADD"								 	 {yylval.e_blend = ADD; return enum_blend;}
"MULTIPLY"								 {yylval.e_blend = MULTIPLY; return enum_blend;}
"MASK"								 	 {yylval.e_blend = MASK; return enum_blend;}
".fill"									 {return fill_tk;}
"FILL_SOLID"							 {yylval.e_fill = FILL_SOLID; return enum_fill;}
"FILL_PLANES"							 {yylval.e_fill = FILL_PLANES; return enum_fill;}
"FILL_LINES"							 {yylval.e_fill = FILL_LINES; return enum_fill;}
"FILL_POINTS"							 {yylval.e_fill = FILL_POINTS; return enum_fill;}
"FILL_NOISE"							 {yylval.e_fill = FILL_NOISE; return enum_fill;}
".atnf"									 {return atnf_tk;}
".pos"				                     {return pos_tk;}
".scale"				                     {return scale_tk;}
".axis"				                     {return axis_tk;}
".angle"				                 {return angle_tk;}
"domain"				                 {return domain_tk;}
".len"      				             {return len_tk;}
".val"     				                 {return val_tk;}
[a-zA-Z_$][a-zA-Z_0-9]*           	     {/*printf(" %s ", yytext)*/;return name;}
[0-9]+\.[0-9]+\.[0-9]+   				 {yylval.ver = atover(yytext); return num_ver;}
[0-9]+                                   {yylval.num = atoi(yytext); return num_int;}
[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?   {yylval.dbl = atof(yytext); return num_flt;}
[ \t\n]                                  ;
[-+=;^*{},]           	                 {return yytext[0];}
.                                        {ECHO; yyerror ("unexpected character");}

%%
int yywrap (void) {return 1;}

void swex(int on){
	if(on){
		BEGIN(EXPR);
	}else{
		BEGIN(INITIAL);
	}
}


void swpr(int on){
	if(on){
		BEGIN(PRIM);
	}else{
		BEGIN(INITIAL);
	}
}

version atover(const char* a){
	version v;
	char* minor;
	char* patch;

	minor = strchr(a, '.')+1;
	patch = strchr(minor, '.')+1;

	*(minor-1) = '\0';
	*(patch-1) = '\0';

	v.major = atoi(a);
	v.minor = atoi(minor);
	v.patch = atoi(patch);

	return v;
}
