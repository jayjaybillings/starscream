/*-----------------------------------------------------------------------------
/
/ Filename: starscream_pf.c
/ Author: Jay Billings
/ Author's email: jayjaybillings@gmail.com
/ Description: The routines in the file calculate the potentials and force
/              fields on a galaxy. Starscream creates galaxies.
/
/              There are several rountines in this file for calculating
/              forces and potentials using direct summation. They are not
/              directly used by any higher routines but are left here as
/              tools for calculating potentials independently. 
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

// This function calculates the potential due to the disk at a point x,y,z
// by interpolating between grid points on the particle mesh. The 
// interpolation routine uses the CIC kernel, which oddly enough is just
// a bilinear interpolation scheme...
//
// If the point lies off of the particle mesh, it approximates the potential
// as a function of 1/r.
double disk_potential_func(galaxy *gal, double x, double y, double z) {

    int i = 0, node_x, node_y, node_z;
    double pot1, pot2, pot3, pot4, pot5, pot6, pot7, pot8;
    double a, r_p, r_max, dx, dy, dz, tx, ty, tz, potential;

    // Scale the coordinates
    r_p = sqrt(x*x + y*y + z*z);
    x = x/gal->space[0] + (double) (Ng/2);
    y = y/gal->space[0] + (double) (Ng/2);
    z = z/gal->space[0] + (double) (Ng/2);

    // Determine the parent node.
    node_x = (int) x;
    node_y = (int) y;
    node_z = (int) z;
    r_max = gal->space[0]*(double) (Ng/2 - Ng/4 + 2);

    // Consider points off the grid or continue
    if (r_p > r_max) {
       if (node_x < Ng/4 + 2) node_x = Ng/4 + 2;
       if (node_x > 3*Ng/4 - 2) node_x = 3*Ng/4 - 2;
       if (node_y < Ng/4 + 2) node_y = Ng/4 + 2;
       if (node_y > 3*Ng/4 - 2) node_y = 3*Ng/4 - 2;
       if (node_z < Ng/4 + 2) node_z = Ng/4 + 2;
       if (node_z > 3*Ng/4 + 2) node_z = 3*Ng/4 - 2;
       a = gal->potential[node_x][node_y][node_z]*r_max;
       potential = a/r_p;
    } else {
       // Check to see if (x,y,z) is a grid point.
       if (x == (double) node_y && y == (double) node_y &&
          z == (double) node_z) {
          // If (x,y,z) is a grid point, return its potential.
          potential = gal->potential[node_x][node_y][node_z];
       } else { 
          // If (x,y,z) is not a grid point, use the CIC
          // interpolation function to calculate the potential.
          pot1 = gal->potential[node_x][node_y][node_z];
          pot2 = gal->potential[node_x+1][node_y][node_z];
          pot3 = gal->potential[node_x][node_y+1][node_z];
          pot4 = gal->potential[node_x][node_y][node_z+1];
          pot5 = gal->potential[node_x][node_y+1][node_z+1];
          pot6 = gal->potential[node_x+1][node_y+1][node_z];
          pot7 = gal->potential[node_x+1][node_y][node_z+1];
          pot8 = gal->potential[node_x+1][node_y+1][node_z+1];
          // CIC fractions
          dx = 1.0 - (x - (double) node_x); 
          dy = 1.0 - (y - (double) node_y); 
          dz = 1.0 - (z - (double) node_z);
          tx = 1.0 - dx;
          ty = 1.0 - dy;
          tz = 1.0 - dz; 
          // Return the interpolated potential.
          potential = dx*dy*dz*pot1 + tx*dy*dz*pot2 +
                      dx*ty*dz*pot3 + dx*dy*tz*pot4 +
                      dx*ty*tz*pot5 + tx*ty*dz*pot6 +
                      tx*dy*tz*pot7 + tx*ty*tz*pot8;
       }
    }

    return potential;
}

// This function calculates the potential due to the disk at a point
// x,y,z using direct summation.
double disk_potential_nbody_func(galaxy *gal, double x, double y, double z) {

    int i;
    double potential, x_p, y_p, z_p, r_p;

    potential = 0.0;
    for (i = 0; i < gal->num_part[0]; ++i) {
        x_p = pow(gal->x[i]-x,2.0);
        y_p = pow(gal->y[i]-y,2.0);
        z_p = pow(gal->z[i]-z,2.0);
        r_p = sqrt(x_p + y_p + z_p);
        potential = potential + G*gal->mass[i]/(r_p*kpc);
    } 

    return -potential;
}

// This is an analytic function for the potential of a halo with mass
// distributed according to a Hernquist profile.
double halo_potential_func(galaxy *gal, double x, double y, double z) {

    double potential, a, r_p;

    a = gal->halo_a_value;
    r_p = sqrt(x*x+y*y+z*z);
    potential = G*gal->halo_mass/((r_p+a)*kpc);

    return -potential;
}

// This function calculates the potential due to the disk at a point x,y,z
// using direct summation.
double halo_potential_nbody_func(galaxy *gal, double x, double y, double z) {
   
    int i; 
    double potential, x_p, y_p, z_p, r_p;

    potential = 0.0;
    for (i = gal->num_part[0]; i < gal->num_part[2]; ++i) {
        x_p = pow(gal->x[i]-x,2.0);
        y_p = pow(gal->y[i]-y,2.0);
        z_p = pow(gal->z[i]-z,2.0);
        r_p = sqrt(x_p + y_p + z_p);
        potential = potential + G*gal->mass[i]/(r_p*kpc);
    } 

    return -potential;
}

// This function returns the cumulative potential along the z axis
// due to the disk and the halo.
double potential_func(galaxy *gal, double x, double y, double z) {

    return halo_potential_func(gal,x,y,z) + disk_potential_func(gal,x,y,z);
}

// This function calculates the force on a test particle due to the disk
// at some point using direct summation. This is used to calculate the 
// velocity dispersion.
double disk_force_nbody_func(galaxy *gal, double x, double y, double z) {

    int i;
    double force, r_p, r_pp;

    force = 0.0;
    r_p = sqrt(x*x + y*y + z*z);
    for (i = 0; i < gal->num_part[0]; ++i) {
        r_pp = sqrt(gal->x[i]*gal->x[i] + gal->y[i]*gal->y[i] +
               gal->z[i]*gal->z[i]);
        force = force + G*(gal->mass[i]*r_p*kpc) / 
                    pow(((r_pp*r_pp + r_p*r_p)*kpc*kpc),3.0/2.0);
    } 

    return force;

}

// This function calculates the force on a test particle due to the halo
// along the r axis. This is used to calculate the velocity dispersion.
double halo_force_func(galaxy *gal, double r_p) {

    double force, mass, a, b;

    a = gal->halo_a_value;
    mass = gal->halo_mass;
    force = G*mass /
          pow((r_p*kpc + a*kpc),2.0);
    
    return force;
}

// This function calculates the force on a test particle due to the halo
// in the rz plane in cylindrical coordinates. This is used to calculate 
// the velocity dispersion.
double halo_force_cyl_func(galaxy *gal, double r_p, double z) {

    double force, mass, a, b;

    a = gal->halo_a_value;
    b = sqrt(r_p*r_p + z*z);
    mass = gal->halo_mass;
    force = G*mass /
          //pow((b*kpc + a*kpc),2.0);
          pow((b + a*kpc),2.0);
    
    return force;
}

