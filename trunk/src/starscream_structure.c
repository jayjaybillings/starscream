/*-----------------------------------------------------------------------------
/
/ Filename: starscream_structure.c
/ Author: Jay Billings
/ Author's email: jayjaybillings@gmail.com
/ Description: Routines for describing the structure of a galaxy,
/              specifically the mass, density, scale lengths, etc. are
/              defined in this file. Starscream creates galaxies. 
/
/	       Starscream uses the GNU Scientific Library (GSL). You can
/	       download the GSL source code from:
/
/		http://www.gnu.org/software/gsl
/
/	       or replace it with another math library.
/
/ Copyright Information:
/
/ Copyright (c) 2008, 2009       Jay Jay Billings
/
/ This program is free software; you can redistribute it and/or modify
/ it under the terms of the GNU General Public License as published by
/ the Free Software Foundation; either version 2 of the License, or
/ (at your option) any later version.
/
/ This program is distributed in the hope that it will be useful,
/ but WITHOUT ANY WARRANTY; without even the implied warranty of
/ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/ GNU General Public License for more details.
/
/ You should have received a copy of the GNU General Public License
/ along with this program; if not, write to the Free Software
/ Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
/
/ The license is also available at:
/
/		http://www.gnu.org/copyleft/gpl.html .
/
/ Date: 2009/06/07
/
*///---------------------------------------------------------------------------

#include "starscream.h"

// Distribution function for an exponential disk
double f_exp_disk(double radius, double mu) {

     return radius * exp(-radius/mu) / (mu*mu);
}

// Distribution function for a Hernquist dark matter halo profile
double f_hernquist_halo(double radius, double mu) {

     return mu * radius / (2.0*(radius + mu)*(radius + mu)*(radius + mu));
}

// The density of the disk.
double density_disk_func(galaxy *gal, double radius, double z) {

     double h, density;

     h = gal->disk_scale_length;
     density = gal->disk_mass/(4.0*pi*0.2*h*h*h *
               pow(kpc,3.0))*exp(-radius/h) /
               cosh(z/(0.2*h));

     return density;
}

// The density of the disk.
double density_disk_cyl_func(galaxy *gal, double radius, double z) {

     double h, density;

     h = gal->disk_scale_length*kpc;
     density = gal->disk_mass/(4.0*pi*0.2*h*h*h)*exp(-radius/h) /
               cosh(z/(0.2*h));

     return density;
}

// The density of the halo.
double density_halo_func(galaxy *gal, double radius) {

     double a, density;

     a = gal->halo_a_value;
     density = gal->halo_mass/(2.0*pi) *
               a/(radius*pow((radius + a), 3.0) *
               pow(kpc,3.0));

     return density;
}

// The density of the halo in cylindrical coordinates.
double density_halo_cyl_func(galaxy *gal, double radius, double z) {

     double a, b, density;

     a = gal->halo_a_value*kpc;
     b = sqrt(radius*radius+z*z);
     density = gal->halo_mass/(2.0*pi) *
               //a/(b*pow((b*kpc + a*kpc), 3.0));
               a/(b*pow((b + a), 3.0));

     return density;
}

// This function determines the scale length of the disk to within 15% or so using
// the fitting formula provided in MMW97.
double disk_scale_length_func(galaxy *gal) {

     int i;
     double base, power, c, f_conc, f, scale, old_scale, v_c;

     c = gal->halo_concentration;
     f_conc = 2.0/3.0+ pow((c/21.5),0.7);
     base = (gal->j_d*gal->lambda)/(0.1*gal->m_d);
     power = (-0.06+2.71*gal->m_d+0.0047*gal->m_d/(gal->j_d*gal->lambda));
     f = pow(base,power)*(1.0-3.0*gal->m_d+5.2*gal->m_d*gal->m_d) *
         (1.0-0.019*c+0.00025*c*c+0.52/c);
     gal->disk_scale_length = (1.0/sqrt(2.0))*(gal->j_d/gal->m_d)*gal->lambda*
                             gal->r200*(f/sqrt(f_conc));

     // Now try to iterate for the correct value.
     for (i = 0; i < 50; ++i) { 
         scale = sqrt(1.0/2.0)*(gal->j_d/gal->m_d)*gal->lambda*gal->r200*
                 2.0*(gal->v200/j_d_func(gal));
                 sqrt(1.0/f_c_func(gal->halo_concentration));
         gal->disk_scale_length = scale;
     } 

     return scale;
}

// Mass distribution for a halo with a Hernquist profile
double halo_mass_func(galaxy *gal, double radius) {

     return gal->halo_mass*radius*radius /
            ((radius+gal->halo_a_value)*(radius+gal->halo_a_value));
}

// Mass distribution for an exponential disk
double disk_mass_func(galaxy *gal, double radius) {
       
    double total_mass, scale;

    scale = gal->disk_scale_length;
    total_mass = gal->disk_mass;
 
    return total_mass*(1.0-(1.0+radius/scale)*exp(-radius/scale));
}

// The next two functions are for calculating the angular momentum
// of the disk. The first is the integral of the disk and the 
// second is the integrand. Please note that this returns the
// specific angular momentum, J_d/M_d, not J_d. The units are cm/s.
double j_d_func(galaxy *gal) {

    int status;
    double result, error, j_d;
    gsl_integration_workspace *w
      = gsl_integration_workspace_alloc(1000); 
    gsl_function F;

    F.function = &dj_d_func;
    F.params = gal;
    gsl_integration_qags(&F,0.0,100.0*gal->disk_scale_length
                         ,0.001,1.0e-7,1000,w,&result,&error);

    gsl_integration_workspace_free(w); 

    return result; 
}

double dj_d_func(double radius, void *params) {

    int status;
    galaxy *gal = (galaxy *) params;
    double I_0, K_0, I_1, K_1, halo_mass, disk_mass, y, a;
    double total_disk_mass, v_c, scale;
    gsl_sf_result result;

    scale = gal->disk_scale_length;
    total_disk_mass = gal->disk_mass;
    a = gal->halo_a_value;
    halo_mass = gal->halo_mass;

    y = radius/(2.0*scale); 
    status = gsl_sf_bessel_I0_e(y,&result);
    I_0 = result.val;
    status = gsl_sf_bessel_K0_e(y,&result);
    K_0 = result.val;
    status = gsl_sf_bessel_I1_e(y,&result);
    I_1 = result.val;
    status = gsl_sf_bessel_K1_e(y,&result);
    K_1 = result.val;
    disk_mass = 2.0*total_disk_mass*y*y*(I_0*K_0 - I_1*K_1); 
    halo_mass = halo_mass*radius*radius/((radius+a)*(radius+a));
    v_c = sqrt(G*(halo_mass/(radius*kpc) + disk_mass/(scale*kpc)));

    return v_c*(radius/scale)*(radius/scale)*exp(-radius/scale);
}

// This function returns the energy fraction function of MMW for the
// angular momentum.
double f_c_func(double c) {

    double a, upper, lower;

    upper = c*(1.0 - 1.0/pow(1.0+c,2.0) - 2.0*log(1.0+c)/(1.0+c));
    a = log(1.0+c) - c/(1.0+c);
    lower = (2.0*pow(a,2.0));

    return upper/lower;
}

// This is the integral g_c for determining the halo azimuthal circular
// velocity fraction f_s.
double g_c_func(double c) {

    int status;
    double result, error;
    gsl_integration_workspace *w
      = gsl_integration_workspace_alloc(1000); 
    gsl_function F;

    F.function = &dg_c_func;
    gsl_integration_qags(&F,0.0,c,0.001,1.0e-7,1000,w,&result,&error);

    gsl_integration_workspace_free(w); 
    
    return result;
}

// This is the integrand for the integral g_c above.
double dg_c_func(double x, void *params) {

    double a, b;

    a = sqrt(log(1.0 + x) - x/(1.0 + x));
    b = pow(x,1.5)/pow(1.0+x,2.0); 

    return a*b;
}

// This is the azimuthal circular velocity fraction, f_s. It is needed
// to determine the azimuthal streaming velocity.
double f_s_func(double c, double lambda) {

    double f_s, g_c, f_c;

    g_c = g_c_func(c);
    f_c = f_c_func(c);
    f_s = 1.5*lambda*sqrt(2.0*c/f_c)*pow(log(1.0+c)-c/(1.0+c),1.5)/g_c;

    return f_s;
}

