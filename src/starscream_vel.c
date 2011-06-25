/*-----------------------------------------------------------------------------
/
/ Filename: starscream_vel.c
/ Author: Jay Billings
/ Author's email: jayjaybillings@gmail.com
/ Description: Starscream requires a lot of functions to describe the
/              velocity structure of a galaxy. These functions are located
/              here. Starscream creates galaxies. 
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

// This function calculates the z-axis velocity moment at a given radius
// for the halo. The scheme developed be Springel, White, Di Matteo, and 
// Hernquist is used.
//
// See the references for more details, but be warned that this routine and
// its children are dark creatures from the depths of the Abyss. Honestly,
// no kidding, I saw it eat a baby.
double v2a_z_halo_func(galaxy *gal) {

    int status;
    double integral, error;
    gsl_integration_workspace *w
       = gsl_integration_workspace_alloc(10000); 
    gsl_function F;

    F.function = &dv2a_z_halo_func;
    F.params = gal;
    gsl_integration_qag(&F,sqrt(pow(gal->storage1[2]*kpc,2)),1.0E4*kpc,
                        epsabs,epsrel,10000,key,w,&integral,&error);

    gsl_integration_workspace_free(w); 
    return integral/density_halo_cyl_func(gal,gal->storage1[0]*kpc,
           gal->storage1[2]*kpc);
}

// This is the integrand for the previous function. It is setup to work with
// the GSL_qags structures.
static double dv2a_z_halo_func(double z, void *params) {

    double integrand, radius;
    galaxy *gal = (galaxy *) params;

    radius = gal->storage1[0]*kpc;
    integrand = density_halo_cyl_func(gal,radius,z)*
                vrz_force_func(gal,radius,z);

    return integrand;
}

// This function calculates PART of the phi axis velocity moment. In
// particular, it calculates the derivative of the radial velocity 
// moment, which is equal to the z-axis velocity moment.
//
// This function uses the galaxy storage variables a lot! If you have
// editted the code to ill effect, check that the storage variables 
// properly reassigned after this function.
//
// The derivative is calculated in the x,y plane using the
// 5-point stencil method.
double v2a_theta_halo_func(galaxy *gal, double radius) {

    double x, y, z, h, v2a_z1, v2a_z2, v2a_z3, v2a_z4, v2a_theta;
    double new_radius;

    // Set the derivative stepsize.
    h = 1.0;

    // Store x,y,z to replace later.
    z = gal->storage1[2];

    // Calculate the radius
    new_radius = radius+2.0*h;
    gal->storage1[0] = new_radius; 
    // Find the velocity moment
    v2a_z1 = get_dispersion(gal->r200,new_radius,z)*
             density_halo_cyl_func(gal,new_radius*kpc,z*kpc);

    new_radius = radius+h;
    gal->storage1[0] = new_radius; 
    v2a_z2 = get_dispersion(gal->r200,new_radius,z)*
             density_halo_cyl_func(gal,new_radius*kpc,z*kpc);
    
    new_radius = radius-h;
    gal->storage1[0] = new_radius; 
    v2a_z3 = get_dispersion(gal->r200,new_radius,z)*
             density_halo_cyl_func(gal,new_radius*kpc,z*kpc);
    
    new_radius = radius-2.0*h;
    gal->storage1[0] = new_radius; 
    v2a_z4 = get_dispersion(gal->r200,new_radius,z)*
             density_halo_cyl_func(gal,new_radius*kpc,z*kpc);

    // Calculate the derivative 
    v2a_theta = (-v2a_z1 + 8.0*v2a_z2 - 8.0*v2a_z3 + v2a_z4)/(12.0*h*kpc);

    // Restore the storage variables
    gal->storage1[0] = radius;
 
    return v2a_theta*radius*kpc/density_halo_cyl_func(gal,radius*kpc,z*kpc);  
}

// This function calculates the z-axis velocity moment at a given radius
// for the disk. The scheme developed by Springel, White, Di Matteo, and 
// Hernquist is used.
//
// See the references for more details, but be warned that this routine and
// its children are dark creatures from the depths of the Abyss. Honestly,
// no kidding, I saw it eat a baby.
double v2a_z_disk_func(galaxy *gal) {

    int status;
    double integral, error;
    gsl_integration_workspace *w
       = gsl_integration_workspace_alloc(10000); 
    gsl_function F;

    F.function = &dv2a_z_disk_func;
    F.params = gal;
    gsl_integration_qag(&F,sqrt(pow(gal->storage1[2]*kpc,2.0)),1.0E4*kpc,
                        epsabs,epsrel,10000,key,w,&integral,&error);

    gsl_integration_workspace_free(w); 
    return integral/density_disk_cyl_func(gal,gal->storage1[0]*kpc,
           gal->storage1[2]*kpc);
}

// This is the integrand for the previous function. It is setup to work with
// the GSL_qags structures.
static double dv2a_z_disk_func(double z, void *params) {

    double integrand, radius;
    galaxy *gal = (galaxy *) params;

    radius = gal->storage1[0]*kpc;
    integrand = density_disk_cyl_func(gal,radius,z)*
                vrz_force_func(gal,radius,z);

    return integrand;
}

// This function calculates the first velocity moment in the azimuthal 
// direction for the disk.
double v2a_theta_disk_func(galaxy *gal, double radius, double v2a_z,
                           double v_c) {

    double gamma_sqrd, kappa_sqrd, h, x, y, z, dforcedr;
    double df1, df2, df3, df4, dens, Q;
    double disk_force, abserr, force;
    gsl_function F;

    // Setup the gsl function for finding the derivative
    F.function = &disk_potential_wrapper_func;
    F.params = gal;
    gsl_deriv_central(&F,radius*kpc,1.0E-8,&disk_force,&abserr);

    // Set the force
    force = halo_force_func(gal,radius) + disk_force/kpc;

    // Save the storage variables
    x = gal->storage1[0]*cos(gal->storage1[1]);
    y = gal->storage1[0]*sin(gal->storage1[1]);
    z = gal->storage1[2];

    // Calculate force and force derivative
    h = 1.0;
    dforcedr = (pow(v_c_func(gal,radius+h),2.0)/((radius+h)*kpc) -
               pow(v_c_func(gal,radius),2.0)/(radius*kpc))/(h*kpc);

    kappa_sqrd = 3.0*force/(radius*kpc) + dforcedr;
    gamma_sqrd = 4.0*force/(kappa_sqrd*radius*kpc);

    dens = gal->disk_mass/(2.0*pi*pow(gal->disk_scale_length*kpc,2.0))*
           exp(-radius/gal->disk_scale_length);
    Q = sqrt(v2a_z*sqrt(kappa_sqrd*kappa_sqrd))/(3.36*G*dens);
    Q_min = (Q < Q_min && Q > 0.0) ? Q : Q_min;

    return v2a_z/gamma_sqrd;
}

// A wrapper for the disk potential function.
double disk_potential_wrapper_func(double radius, void *params) {

    double theta, ans, x, y;
    galaxy *gal = (galaxy *) params;

    x = radius*cos(gal->storage1[1])/kpc;
    y = radius*sin(gal->storage1[1])/kpc;
    ans = disk_potential_func(gal,x,y,gal->storage1[2]);   

    return ans;
}

// This function calculates the force in the radial direction for the
// previous function. It uses the 2D five-point stencil.
static double vtheta_force_func(galaxy *gal, double radius) {

    double force, x, y, z, h, theta, p1, p2, p3, p4, new_radius;
    double force2;

    z = gal->storage1[0];
    theta = atan(gal->storage1[2]/gal->storage1[1]);
    h = 1.0/kpc;

    // Calculate the potential at the different radii.
    new_radius = radius + 2.0*h;
    x = new_radius*cos(theta);
    y = new_radius*sin(theta);
    p1 = -disk_potential_func(gal,x,y,z);
    new_radius = radius + h;
    x = new_radius*cos(theta);
    y = new_radius*sin(theta);
    p2 = 8.0*disk_potential_func(gal,x,y,z);
    new_radius = radius - h;
    x = new_radius*cos(theta);
    y = new_radius*sin(theta);
    p3 = -8.0*disk_potential_func(gal,x,y,z);
    new_radius = radius - 2.0*h;
    x = new_radius*cos(theta);
    y = new_radius*sin(theta);
    p4 = disk_potential_func(gal,x,y,z);

    // Calculate the force. Use the five point rule for the disk.
    force = halo_force_func(gal,radius) + 
           (p1+p2+p3+p4)/(12.0*h*kpc);

    return force;
}

// The circular velocity function. This function returns the velocity
// in units of cm/s.
double v_c_func(galaxy *gal, double radius) {
      
    int status; 
    double I_0, K_0, I_1, K_1, halo_mass, disk_mass, y;
    gsl_sf_result result;

    y = radius/(2.0*gal->disk_scale_length); 
    status = gsl_sf_bessel_I0_e(y,&result);
    I_0 = result.val;
    status = gsl_sf_bessel_K0_e(y,&result);
    K_0 = result.val;
    status = gsl_sf_bessel_I1_e(y,&result);
    I_1 = result.val;
    status = gsl_sf_bessel_K1_e(y,&result);
    K_1 = result.val;
    disk_mass = 2.0*gal->disk_mass*y*y*(I_0*K_0 - I_1*K_1); 
    halo_mass = halo_mass_func(gal,radius);

    return sqrt(G*(halo_mass/(radius*kpc) + disk_mass/(gal->disk_scale_length*kpc)));

}

// The circular velocity function for the disk
double v_c_disk_func(galaxy *gal, double radius) {
      
    int status; 
    double I_0, K_0, I_1, K_1, halo_mass, disk_mass, y;
    gsl_sf_result result;

    y = radius/(2.0*gal->disk_scale_length); 
    status = gsl_sf_bessel_I0_e(y,&result);
    I_0 = result.val;
    status = gsl_sf_bessel_K0_e(y,&result);
    K_0 = result.val;
    status = gsl_sf_bessel_I1_e(y,&result);
    I_1 = result.val;
    status = gsl_sf_bessel_K1_e(y,&result);
    K_1 = result.val;
    disk_mass = 2.0*gal->disk_mass*y*y*(I_0*K_0 - I_1*K_1); 

    return sqrt(G*disk_mass/(gal->disk_scale_length*kpc));

}
  
// The circular velocity function for the halo
double v_c_halo_func(galaxy *gal, double radius) {
      
    int status; 
    double I_0, K_0, I_1, K_1, halo_mass, disk_mass, y;
    gsl_sf_result result;

    halo_mass = halo_mass_func(gal,radius);

    return sqrt(G*(halo_mass/(radius*kpc)));

}

// This function calculates the force on a test particle due to the disk
// at some point z. This is used to calculate the velocity dispersion.
//
// This function simply performs a five point differentiation around 
// x,y,z in the z direction with a small stepsize.
//
// This function is specially formatted for the velocity dispersion
// routines.
static double vrz_disk_force_func(galaxy *gal, double z) {

    int i;
    double force, h, x, y;

    z = z/kpc;
    x = gal->storage1[0]*cos(gal->storage1[1]);
    y = gal->storage1[0]*sin(gal->storage1[1]);
    h = 1.0;
    force = (-disk_potential_func(gal,x,y,z+2.0*h) +
            8.0*disk_potential_func(gal,x,y,z+1.0*h) -
            8.0*disk_potential_func(gal,x,y,z-1.0*h) +
            disk_potential_func(gal,x,y,z-2.0*h))/(12.0*h*kpc);

    //printf("%le %le disk force\n",force);

    return force;

}

// This function returns the cumulative force along the z axis for the
// halo and at the point r,z for the disk.
//
// This function is specially formatted for the velocity dispersion
// routines.
static double vrz_force_func(galaxy *gal, double r_p, double z) {

    return halo_force_cyl_func(gal,r_p,z) + vrz_disk_force_func(gal,z);

}
