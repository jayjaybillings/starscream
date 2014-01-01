/*-----------------------------------------------------------------------------
/
/ Filename: starscream_init.c
/ Author: Jay Billings
/ Author's email: jayjaybillings@gmail.com
/ Description: These are the initialization routines for Starscream. 
/              They include routines to set up the memory structure and 
/              allocate particle positions and velocities. Starscream 
/              creates galaxies. 
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

// Get the configuration header
//#include "../config.h"
#include "starscreamConfig.h"

// Include the starscream header
#include "starscream.h"

// Create a galaxy
int create_galaxy(galaxy *gal, int parts[2], double m_d, double j_d, double lambda, 
                 double c, double v200, int Ngrid, double space, int info) {

    int i,j,num_part;
    long seed;

    /* A list of universal constants. Most constants were calculated on an AMD Turion
    64x2 TL-50 with gfortran v.4.1.2, a Texas Instruments TI-83 Plus, or a Texas 
    Instruments TI-36X Solar. Some of the constants were looked up in the literature.
    
    These values are global, but assigned values here.                               */
    pi = 3.1415926535897932;   // pi taken from "100000 digits of pi," 
 			       // http://www.geom.uiuc.edu/~huberty/math5337/groupe/digits.html
    G = 6.67428E-8;            // G = 6.67428E-8 cm^3 g^-1 s^-2
    H = 2.300952983428601E-18; // 7.1E6/(1.0E3*kpc) s^-1 or 71.0 km s^-1 Mpc^-1
    unit_mass = 1.989E43;      // 1.989E43 = 1E10 Solar masses
    kpc = 3.085678E21;         // 1 kpc = 3.085678E21 cm
    crit_dens = 3.0*H*H/(8.0*pi*G); // The background critical density, (overdensity)
    Ng = 2*Ngrid;

    // Create the random number generator environment.
    if (random_number_set != 1) {
       seed = time(NULL);
       gsl_rng_env_setup();
       T = gsl_rng_default;
       r = gsl_rng_alloc(T);
       gsl_rng_set(r,seed);
       random_number_set = 1;
    }

    // Reset Q_min to some arbitrary value so it can be calculated later.
    Q_min = 10.0;

    // Set the parameters.
    gal->num_part[0] = parts[0];
    gal->num_part[1] = parts[1];
    gal->m_d = m_d;
    gal->j_d = j_d;
    gal->lambda = lambda;
    gal->halo_concentration = c;
    gal->v200 = v200;
    gal->space[0] = space;
    gal->num_part[2] = gal->num_part[0] + gal->num_part[1];
    num_part = gal->num_part[2];

    // Set a bunch of constants that define the galaxy's structure.
    gal->total_mass = (gal->v200*gal->v200*gal->v200)/(10.0*G*H);
    gal->disk_mass = gal->m_d*gal->total_mass;
    gal->halo_mass = (1.0 - gal->m_d)*gal->total_mass;
    gal->r200 = gal->v200/(10.0*H*kpc);
    gal->halo_scale_length = gal->r200/gal->halo_concentration;
    gal->halo_a_value = gal->halo_scale_length*sqrt(2.0*(log(1.0+c) - c/(1.0+c)));
    // Set the disk scale length, but check to see if it is a copy first.
    if (a_copy == 0) disk_scale_length_func(gal);

    // Print some information to the screen.
    if (info != 0) {
       fprintf(stderr,"Creating a galaxy with a total of %d particles.\n",
              gal->num_part[2]);
       fprintf(stderr,"\tNumber of disk particles: %d\n",gal->num_part[0]);
       fprintf(stderr,"\tNumber of halo particles: %d\n",gal->num_part[1]);
       fprintf(stderr,"The total mass of this galaxy is %le grams.\n",
              gal->total_mass);
       fprintf(stderr,"\tDisk mass: %le g\n\tHalo mass: %le g\n",gal->disk_mass,
              gal->halo_mass);
       fprintf(stderr,"The virial radius of this galaxy is %lf kpc.\n",
              gal->r200);
       fprintf(stderr,"\tDisk scale length: %lf kpc\n",gal->disk_scale_length);
       fprintf(stderr,"\tHalo scale length: %lf kpc\n",gal->halo_scale_length);
       fprintf(stderr,"\tHalo 'a' parameter = %lf kpc\n",gal->halo_a_value);
       fprintf(stderr,"Potential grid has dimensions x=%d y=%d z=%d,",
              Ngrid,Ngrid,Ngrid);
       fprintf(stderr,"(%d points)\n",Ngrid*Ngrid*Ngrid);
       fprintf(stderr,"\tDimensions are padded to x=%d y=%d z=%d, (%d points)",
             Ng,Ng,Ng,Ng*Ng*Ng);
       fprintf(stderr,"\n");
    }     

    // Allocate particle id numbers array.
    if (!(gal->id=calloc(num_part,sizeof(int)))) {
       fprintf(stderr,"Unable to allocate particle ID numbers.\n");
       return -1;
    }
    // Allocate x coordinates for all the particles.
    if (!(gal->x=calloc(num_part,sizeof(double)))) {
       fprintf(stderr,"Unable to allocate particle x coordinates.\n");
       return -1;
    }
    // Allocate y coordinates for all the particles.
    if (!(gal->y=calloc(num_part,sizeof(double)))) {
       fprintf(stderr,"Unable to allocate particle y coordinates.\n");
       return -1;
    }
    // Allocate z coordinates for all the particles.
    if (!(gal->z=calloc(num_part,sizeof(double)))) {
       fprintf(stderr,"Unable to allocate particle z coordinates.\n");
       return -1;
    }
    // Allocate x velocities for all the particles.
    if (!(gal->vel_x=calloc(num_part,sizeof(double)))) {
       fprintf(stderr,"Unable to allocate particle x coordinates.\n");
       return -1;
    }
    // Allocate y velocities for all the particles.
    if (!(gal->vel_y=calloc(num_part,sizeof(double)))) {
       fprintf(stderr,"Unable to allocate particle y coordinates.\n");
       return -1;
    }
    // Allocate z velocities for all the particles.
    if (!(gal->vel_z=calloc(num_part,sizeof(double)))) {
       fprintf(stderr,"Unable to allocate particle z coordinates.\n");
       return -1;
    }
    // Allocate masses for all the particles.
    if (!(gal->mass=calloc(num_part,sizeof(double)))) {
       fprintf(stderr,"Unable to allocate particle masses.\n");
       return -1;
    }
    // Turn on and allocate the potential grid, start with x-axis
    gal->potential_defined = 1;
    if (!(gal->potential=calloc(Ng,sizeof(double *)))) {
       fprintf(stderr,"Unable to create potential x axis.\n");
       return -1;
    }
    for (i = 0; i < Ng; ++i) {
        // y-axis
        if (!(gal->potential[i] = calloc(Ng,sizeof(double *)))) {
           fprintf(stderr,"Unable to create potential y axis.\n");
           return -1;
        }
        // z-axis
        for (j = 0; j < Ng; ++j) {
            if (!(gal->potential[i][j] = calloc(Ng,sizeof(double)))) {
               fprintf(stderr,"Unable to create potential z axis.\n");
               return -1;
            }
        }
    }

    // Set the masses and ids for all the particles.
    for (i = 0; i < gal->num_part[0]; ++i) {
        gal->mass[i] = gal->disk_mass/gal->num_part[0];
        gal->id[i] = i;
    }
    for (i = gal->num_part[0]; i < gal->num_part[2]; ++i) {
        gal->mass[i] = gal->halo_mass/gal->num_part[1];
        gal->id[i] = i;
    }

    return 0;
}

// Use this to allocate the galaxy storage variable on the fly.
void allocate_galaxy_storage_variable(galaxy *gal, int size) {

     // Allocate the storage array.
     if (!(gal->storage1=calloc(size,sizeof(double)))) {
        fprintf(stderr,"Unable to allocate particle masses.\n");
        return;
     }

    return;
}

// Use this to allocate the galaxy storage variable on the fly.
void allocate_dispersion_grid() {

     int i;

     // Allocate the dispersion grid.
     if (!(v2rz=calloc(Ng/2,sizeof(double *)))) {
        fprintf(stderr,"Unable to allocate dispersion grid.\n");
        return;
     }
     for (i = 0; i < Ng/2; ++i) {
         if (!(v2rz[i]=calloc(Ng/2,sizeof(double)))) {
            fprintf(stderr,"Unable to allocate dispersion grid.\n");
         }
     }

    return;
}

// Use this to destroy the dispersion grid.
void destroy_dispersion_grid() {

     int i;

     // Deallocate the dispersion grid.
     for (i = 0; i < Ng/2; ++i) {
         free(v2rz[i]);
     }
     free(v2rz);

     return;
}

// Set the positions of the particles in the galaxy to be consistent with the 
// various particle distributions. This function uses the Metropolis algorithm
// to allocate particle positions from the distribution functions of each 
// galaxy component. 
int set_galaxy_coords(galaxy *gal) {

    int i;
    double mu, w_x, w_y, prob, *radius;
    double theta, phi, a, b, c;

    // Allocate some space for storing the previous radii.
    if (!(radius = calloc(gal->num_part[2],sizeof(double)))) {
       fprintf(stderr,"Unable to allocate radial position vector in ");
       fprintf(stderr,"set_galaxy_coords.\n");
       return -1;    
    } 

    // Use the Metropolis algorithm to place the disk particles.
    
    // r and theta
    mu = gal->disk_scale_length;
    w_x = f_exp_disk(mu,mu);
    radius[0] = mu;
    theta = 2.0*pi*gsl_rng_uniform_pos(r);
    gal->x[0] = radius[0]*cos(theta);
    gal->y[0] = radius[0]*sin(theta);
    for (i = 1; i < gal->num_part[0]; ++i) {
        a = 10.0*mu*gsl_rng_uniform_pos(r);
        w_y = f_exp_disk(a,mu);
        prob = w_y/w_x;
        if (prob >= 1.0) {
           radius[i] = a;
        } else {
           b = gsl_rng_uniform_pos(r);
           if (b <= prob) {
              radius[i] = a;
           } else {
              radius[i] = radius[i-1];
           }
        }
        theta = 2.0*pi*gsl_rng_uniform_pos(r);
        gal->x[i] = radius[i]*cos(theta);
        gal->y[i] = radius[i]*sin(theta);
    }

    // z
    w_x = 1.0/(0.4*mu*cosh(0.0)*cosh(0.0));
    gal->z[0] = 0.0;
    for (i = 1; i < gal->num_part[0]; ++i) {
        a = -0.2*mu + 0.4*mu*gsl_rng_uniform_pos(r);
        w_y = 1.0/(0.4*mu*cosh(a/(0.2*mu))*cosh(a/(0.2*mu)));
        prob = w_y/w_x;
        if (prob >= 1.0) {
           gal->z[i] = a;
       } else {
           b = gsl_rng_uniform_pos(r);
           if (b <= prob) {
             gal->z[i] = a;
           } else {
             gal->z[i] = gal->z[i-1];
           }
       }
    }

    // Use the Metropolis algorithm to place the halo particles. Note that the 
    // radius changed from a cylindrical radius to a spherical radius!
    c = gal->halo_concentration;
    mu = gal->halo_a_value;
    w_x = f_hernquist_halo(mu,mu);
    radius[gal->num_part[0]] = mu;
    theta = 2.0*pi*gsl_rng_uniform_pos(r);
    phi = acos(-1.0 + 2.0*gsl_rng_uniform_pos(r));
    gal->x[gal->num_part[0]] = radius[gal->num_part[0]]*cos(theta)*sin(phi);
    gal->y[gal->num_part[0]] = radius[gal->num_part[0]]*sin(theta)*sin(phi);
    gal->z[gal->num_part[0]] = radius[gal->num_part[0]]*cos(phi);
    for (i = gal->num_part[0] + 1; i < gal->num_part[2]; ++i) {
        a = 10.0*mu*gsl_rng_uniform_pos(r);
        w_y = f_hernquist_halo(a,mu);
        prob = w_y/w_x;
        if (prob >= 1.0) {
           radius[i] = a;
        } else {
           b = gsl_rng_uniform_pos(r);
           if (b <= prob) {
              radius[i] = a;
           } else {
              radius[i] = radius[i-1];
           }
        }
        theta = 2.0*pi*gsl_rng_uniform_pos(r);
        phi = acos(-1.0 + 2.0*gsl_rng_uniform_pos(r));
        gal->x[i] = radius[i]*cos(theta)*sin(phi);
        gal->y[i] = radius[i]*sin(theta)*sin(phi);
        gal->z[i] = radius[i]*cos(phi);
    }

    free(radius);
    return 0;
}

// Set the velocities of the particles in the galaxy. Users can choose from
// Keplerian orbits or pick velocities from a distribution.
int set_galaxy_velocity(galaxy *gal, int vel_type) {

    int i, status;
    double alpha, v_c, v_theta, v_phi, y, halo_mass, disk_mass, radius, theta;
    double phi, z, f_s, v_r, v_z, sigma_theta;
    double I_0, K_0, I_1, K_1, vel_rad_disp, v2a_r, v2a_z, v2a_theta, va_theta;
    gsl_sf_result result;

    // Allocate and utilize the storage variable of the galaxy type for passing
    // information to the dispersion function. It's a hack, I know...
    allocate_galaxy_storage_variable(gal,3);

    // Check to see if the user wants to use Keplerian motion.
    if (vel_type == 0) {
       // Disk particle velocities.
       for (i = 0; i < gal->num_part[0]; ++i) {
           radius = sqrt(gal->x[i]*gal->x[i] + gal->y[i]*gal->y[i]);
           theta = atan2(gal->y[i],gal->x[i]) ;
           v_c = v_c_func(gal,radius);
           //Make sure to divide by 1.0E5 to put the velocities in km/s.
           gal->vel_x[i] = (-v_c*sin(theta))/1.0E5;
           gal->vel_y[i] = (v_c*cos(theta))/1.0E5;
           gal->vel_z[i] = 0.0;
       }
       // Halo particle velocities.
       for (i = gal->num_part[0]; i < gal->num_part[2]; ++i) {
           radius = sqrt(gal->x[i]*gal->x[i] + gal->y[i]*gal->y[i] + gal->z[i]*gal->z[i]);
           theta = atan2(gal->y[i],gal->x[i]) ;
           phi = acos(gal->y[i]/radius);
           v_c = v_c_func(gal,radius);
           alpha = 2.0*pi*gsl_rng_uniform_pos(r);
           v_theta = v_c*cos(alpha);
           v_phi = v_c*sin(alpha);
           //Make sure to divide by 1.0E5 to put the velocities in km/s.
           gal->vel_x[i] = (v_theta*sin(-theta) + v_phi*cos(theta)*cos(phi))/1.0E5;
           gal->vel_y[i] = (v_theta*cos(theta) + v_phi*sin(theta)*cos(phi))/1.0E5;
           gal->vel_z[i] = (v_phi*sin(-phi))/1.0E5;
       }
    // Check to see if the user wants to use the Springel-Hernquist description
    // of velocity dispersion to calculate the velocity structure.
    } else {
       // Allocate dispersion grid
       allocate_dispersion_grid();
       // Disk particle velocities.
       set_dispersion_grid(gal,gal->disk_scale_length*10.0,v2a_z_disk_func);
       for (i = 0; i < gal->num_part[0]; ++i) {
           radius = sqrt(gal->x[i]*gal->x[i] + gal->y[i]*gal->y[i]);
           theta = atan2(gal->y[i],gal->x[i]);
           z = gal->z[i];
           v_c = v_c_func(gal,radius);
           gal->storage1[0] = radius;
           gal->storage1[1] = theta;
           gal->storage1[2] = z;
           v2a_z = get_dispersion(gal->disk_scale_length*10.0,radius,z);
           v2a_r = v2a_z;
           v2a_theta = v2a_theta_disk_func(gal,radius,v2a_z,v_c);
           v_r = gsl_ran_gaussian_ziggurat(r,sqrt(v2a_z));
           v_theta = gsl_ran_gaussian_ziggurat(r,sqrt(v2a_theta)) + v_c;
           v_z = gsl_ran_gaussian_ziggurat(r,sqrt(v2a_z));
           //Make sure to divide by 1.0E5 to put the velocities in km/s.
           gal->vel_x[i] = (v_r*cos(theta)-v_theta*sin(theta))/1.0E5;
           gal->vel_y[i] = (v_r*sin(theta)+v_theta*cos(theta))/1.0E5;
           gal->vel_z[i] = v_z/1.0E5;
       }
       // Halo particle velocities.
       f_s = f_s_func(gal->halo_concentration,gal->lambda);
       set_dispersion_grid(gal,gal->r200,v2a_z_halo_func);
       for (i = gal->num_part[0]; i < gal->num_part[2]; ++i) {
           radius = sqrt(gal->x[i]*gal->x[i] + gal->y[i]*gal->y[i]);
           theta = atan2(gal->y[i],gal->x[i]) ;
           z = gal->z[i];
           v_c = v_c_func(gal,radius);
           gal->storage1[0] = radius;
           gal->storage1[1] = theta;
           gal->storage1[2] = z;
           v2a_z = get_dispersion(gal->r200,radius,z);
           v2a_z = sqrt(pow(v2a_z,2.0));
           v2a_r = sqrt(pow(v2a_z,2.0));
           v2a_theta = v2a_r + v2a_theta_halo_func(gal,radius) + v_c*v_c; 
           va_theta = f_s*v_c;
           sigma_theta = sqrt(sqrt(pow(v2a_theta - va_theta*va_theta,2.0)));
           v_r = gsl_ran_gaussian_ziggurat(r,sqrt(v2a_r));
           v_theta = gsl_ran_gaussian_ziggurat(r,sigma_theta) + va_theta;
           v_z = gsl_ran_gaussian_ziggurat(r,sqrt(v2a_z));
           //Make sure to divide by 1.0E5 to put the velocities in km/s.
           gal->vel_x[i] = (v_r*cos(theta)-v_theta*sin(theta))/1.0E5;
           gal->vel_y[i] = (v_r*sin(theta)+v_theta*cos(theta))/1.0E5;
           gal->vel_z[i] = v_z/1.0E5;
       }
       destroy_dispersion_grid();
    }

    free(gal->storage1);
    return 0;
}

// This function calculates and stores the velocity dispersion in the rz
// plane on a grid. These values can then be used to interpolate the 
// velocity dispersion instead of calculating it about 5*N times for N
// particles.
static int set_dispersion_grid(galaxy *gal, double radius, 
                               double disp_func(galaxy *)) {

    int i, j;
    double h;

    h = 2.0*(3.0*radius)/Ng;
    // Skip r = 0 and break the grid in half to avoid the singularity at 0.
    for (i = 1; i < Ng/2; ++i) {
        // Top half
        for (j  = Ng/4; j < Ng/2; ++j) {
            gal->storage1[0] = ((double) i)*h;
            gal->storage1[1] = 0.0;
            gal->storage1[2] = ((double) (j - Ng/4))*2.0*h;
            v2rz[i][j] = disp_func(gal);
        }
        // Bottom half
        for (j  = 0; j < Ng/4; ++j) {
            gal->storage1[0] = ((double) i)*h;
            gal->storage1[1] = 0.0;
            gal->storage1[2] = ((double) j - Ng/4)*2.0*h;
            v2rz[i][j] = disp_func(gal);
        }
    }

    return 0;
}

// This routine calculates the dispersion at a point by using bilinear
// interpolation on the dispersion grid.
double get_dispersion(double cutoff, double radius, double z) {

    int x1, x2, y1, y2;
    double v2a, h, r1, r2, z1, z2, denom;

    // Set stepsize
    h = 2.0*(3.0*cutoff)/Ng;

    // Check for bad radius
    if (radius < h) {
       return 0.0;
    }

    // Set constants
    r1 = (double) ((int) radius);
    r2 = r1 + h;
    z1 = (double) ((int) z);
    z2 = z1 + 2.0*h;
    denom = (r2 - r1)*(z2 - z1);
    x1 = (int) (r1/h);
    x2 = (int) (r2/h);
    y1 = (int) (z1/(2.0*h)) + Ng/4;
    y2 = (int) (z2/(2.0*h)) + Ng/4;

    v2a = (v2rz[x1][y1]/denom)*(r2-radius)*(z2-z) +
          (v2rz[x2][y1]/denom)*(radius-r1)*(z2-z) +
          (v2rz[x1][y2]/denom)*(r2-radius)*(z-z1) +
          (v2rz[x2][y2]/denom)*(radius-r1)*(z-z1);

    return v2a;
}

// This function calculates the potential due to a galactic disk using Cloud-
// In-Cell mass assignment with vacuum (isolated) boundary conditions on a 
// Cartesian grid. The method was adapted from the discussion found in Hockney
// and Eastwood, "Computer Simulation Using Particles," 1981.
//
// This function should not be called to find the potential of a single point
// because it does not return such a value! Use the <component>_potential_func() 
// routines for that!
//
// The function requires a pointer to a galaxy as one argument and an optional
// dump argument. If a user sets the dump option, the disk potential in the x,y,
// z=0, plane is dumped to standard out. If the user catches the output using a
// pipe to a file, it can be plotted simply in gnuplot by typing at the prompt,
//
// 	splot '<piped_file>' w lines
//
// where <piped_file> is the name of the file containing the output and the 
// apostrophes are required.
//
// Users may notice that this function has space_x, space_y, and space_z
// defined instead of a single space variable. I leave these here because
// an earlier version had adaptive grid sizes and some users may want to try
// playing around with them now. On the flip side, there's a lot of other 
// stuff from that early version that I took out. Good luck...
//
int set_galaxy_potential(galaxy *gal, int dump) {

    int i, j, k, l, node_x, node_y, node_z;
    double space_x, space_y, space_z, dx, dy, dz, tx, ty, tz, n;
    double x[gal->num_part[0]], y[gal->num_part[0]], z[gal->num_part[0]];
    size_t p_x[gal->num_part[0]], p_y[gal->num_part[0]], p_z[gal->num_part[0]];
    fftw_plan fft_green, fft_potential;

    // Setup fftw threads
    #if USE_FFTW_THREADS == 1 
    printf("Setting up threaded FFT call..."); 
    fftw_init_threads();
    fftw_plan_with_nthreads(2);
    printf(" Done!\n"); 
    #endif

    // Allocate grid storage variables
    /*if (!(green = calloc(pow(Ng,3),sizeof(double)))) {
       fprintf(stderr,"Unable to allocate space for Green's function.\n");
       return -1;
    }
    if (!(potential = calloc(pow(Ng,3),sizeof(double)))) {
       fprintf(stderr,"Unable to allocate space for potential buffer.\n");
       return -1;
    }*/
    if (!(green = fftw_malloc(Ng*Ng*Ng*sizeof(double)))) {
       fprintf(stderr,"Unable to allocate space for Green's function.\n");
       return -1;
    }
    if (!(potential = fftw_malloc(Ng*Ng*Ng*sizeof(double)))) {
       fprintf(stderr,"Unable to allocate space for potential buffer.\n");
       return -1;
    }

    if (!(green_grid=calloc(Ng,sizeof(double *)))) {
       fprintf(stderr,"Unable to create Green's function x axis.\n");
       return -1;
    }
    for (i = 0; i < Ng; ++i) {
        // y-axis
        if (!(green_grid[i] = calloc(Ng,sizeof(double *)))) {
           fprintf(stderr,"Unable to create Green's function y axis.\n");
           return -1;
        }
        // z-axis
        for (j = 0; j < Ng; ++j) {
            if (!(green_grid[i][j] = calloc(Ng,sizeof(double)))) {
               fprintf(stderr,"Unable to create Green's function z axis.\n");
               return -1;
            }
        }
    }

    // Allocate the fftw complex output value and the fftw dft plan.
    fft_potential = fftw_plan_r2r_3d(Ng,Ng,Ng,potential,potential,FFTW_REDFT00,
                                    FFTW_REDFT00,FFTW_REDFT00,FFTW_ESTIMATE);
    fft_green = fftw_plan_r2r_3d(Ng,Ng,Ng,green,green,FFTW_REDFT00,
                                FFTW_REDFT00,FFTW_REDFT00,FFTW_ESTIMATE);
    // Normalization constant
    n = (2.0*(double) (Ng) - 2.0) *
        (2.0*(double) (Ng) - 2.0) *
        (2.0*(double) (Ng) - 2.0);

    // Check for bad grids
    if (Ng <= 0) {
       fprintf(stderr,"Grid dimensions must be greater than zero!");
       return -1;
    }

    // Sort the position arrays and figure out the spacing between grid points.
    // Subtract 2 for a.) the C offset and b.) the CIC offset. Finally, store
    // the values in the galaxy for later use.
    gsl_sort_index(p_x,gal->x,1,gal->num_part[0]);
    gsl_sort_index(p_y,gal->y,1,gal->num_part[0]);
    gsl_sort_index(p_z,gal->z,1,gal->num_part[0]);
    space_x = gal->space[0];
    space_y = space_x;
    space_z = space_x;

    // Make the coordinate information unitless.
    for (i = 0; i < gal->num_part[0]; ++i) {
        x[i] = gal->x[i]/(space_x) + (double) (Ng/4);
        y[i] = gal->y[i]/(space_y) + (double) (Ng/4);
        z[i] = gal->z[i]/(space_z) + (double) (Ng/4);
    }

    // Calculate the density using the CIC routine. The positions are shifted 
    // such that the particles are in the +x,+y,+z octant. space_* is added 
    // to take care of the vacuum boundary conditions. The density values are
    // stored in the potential structure for now.
    for (i = 0; i < gal->num_part[0]; ++i) {
        // Figure out which node owns the particle
        node_x = (int) x[i];
        node_y = (int) y[i];
        node_z = (int) z[i];
        // Set the CIC size fractions
        dx = 1.0 - (x[i] - (double) node_x);
        dy = 1.0 - (y[i] - (double) node_y);
        dz = 1.0 - (z[i] - (double) node_z);
        tx = 1.0 - dx;
        ty = 1.0 - dy;
        tz = 1.0 - dz;
        // Calculate the CIC densities
        gal->potential[node_x][node_y][node_z] += gal->mass[i]*(dx*dy*dz) /
                                                  (unit_mass*space_x*space_y*
                                                  space_z);
        gal->potential[node_x+1][node_y][node_z] += gal->mass[i]*(tx*dy*dz) /
                                                    (unit_mass*space_x*
                                                    space_y*space_z);
        gal->potential[node_x][node_y+1][node_z] += gal->mass[i]*(dx*ty*dz) /
                                                    (unit_mass*space_x*
                                                    space_y*space_z);
        gal->potential[node_x][node_y][node_z+1] += gal->mass[i]*(dx*dy*tz) /
                                                    (unit_mass*space_x*
                                                    space_y*space_z);
        gal->potential[node_x][node_y+1][node_z+1] += gal->mass[i]*(dx*ty*tz) /
                                                      (unit_mass*space_x*
                                                    space_y*space_z);
        gal->potential[node_x+1][node_y+1][node_z] += gal->mass[i]*(tx*ty*dz) /
                                                      (unit_mass*space_x*
                                                    space_y*space_z);
        gal->potential[node_x+1][node_y][node_z+1] += gal->mass[i]*(tx*dy*tz) /
                                                      (unit_mass*space_x*
                                                    space_y*space_z);
        gal->potential[node_x+1][node_y+1][node_z+1] += gal->mass[i]*(tx*ty*tz) /
                                                        (unit_mass*space_x*
                                                    space_y*space_z);
        // NGP if anyone needs it
        //gal->potential[node_x][node_y][node_z] += gal->mass[i]/(unit_mass *
        //                                          (space_x*space_y*space_z));
    }

    // Define Green's function. dx, dy, dz, tx, ty, and tz are reused because the author is 
    // lazy -- these are the grid points as measured from the center of Green's function 
    // and the local value of the truncation function. The density is also packed into a 
    // buffer here.
    //
    // Green's function is defined eight times here, once per octant, to take care of the
    // isolated (vacuum) boundary conditions. See Hockney and Eastwood, 1980, ch. 6 for a
    // discussion. The octants start in the lower left at (p,q) = (0,0) and progress
    // counter clockwise.
    for (i = 0; i < Ng/2; ++i) {
        for (j = 0; j < Ng/2; ++j) {
            for (k = 0; k < Ng/2; ++k) {
                dx = sqrt(pow((double) (i-Ng/4),2.0));
                dy = sqrt(pow((double) (j-Ng/4),2.0));
                dz = sqrt(pow((double) (k-Ng/4),2.0));
                tx = dx/((double) (Ng/2));
                ty = dy/((double) (Ng/2));
                tz = dz/((double) (Ng/2));
                // Octant 1
                green_grid[i][j][k] = truncation_func(tx,ty,tz)
                           /(4.0*pi*sqrt(dx*dx + dy*dy + dz*dz));
                // Octant 2
                green_grid[Ng-1-i][j][k] = green_grid[i][j][k];
                // Octant 3
                green_grid[Ng-1-i][Ng-1-j][k] = green_grid[i][j][k];
                // Octant 4
                green_grid[i][Ng-1-j][k] = green_grid[i][j][k];
                // Octant 5
                green_grid[i][j][Ng-1-k] = green_grid[i][j][k];
                // Octant 6
                green_grid[Ng-1-i][j][Ng-1-k] = green_grid[i][j][k];
                // Octant 7
                green_grid[Ng-1-i][Ng-1-j][Ng-1-k] = green_grid[i][j][k];
                // Octant 8
                green_grid[i][Ng-1-j][Ng-1-k] = green_grid[i][j][k];
            }
        }
    }
    // Set singularities in Green's function to 1.0
    green_grid[Ng/4][Ng/4][Ng/4] = 1.0;
    green_grid[Ng-1-Ng/4][Ng/4][Ng/4] = 1.0;
    green_grid[Ng-1-Ng/4][Ng-1-Ng/4][Ng/4] = 1.0;
    green_grid[Ng/4][Ng-1-Ng/4][Ng/4] = 1.0;
    green_grid[Ng/4][Ng/4][Ng-1-Ng/4] = 1.0;
    green_grid[Ng-1-Ng/4][Ng/4][Ng-1-Ng/4] = 1.0;
    green_grid[Ng-1-Ng/4][Ng-1-Ng/4][Ng-1-Ng/4] = 1.0;
    green_grid[Ng/4][Ng-1-Ng/4][Ng-1-Ng/4] = 1.0;
 
    // Pack Green's function and the density into 1D arrays
    l = 0;
    for (i = 0; i < Ng; ++i) {
        for (j = 0; j < Ng; ++j) {
            for (k = 0; k < Ng; ++k) {
                green[l] = green_grid[i][j][k];
                potential[l] = gal->potential[i][j][k];
                ++l;
            }
        }   
    }
    
    // Perform the fourier transforms. Density first, Green's function second.
    fftw_execute(fft_potential);
    fftw_execute(fft_green);

    // Multiply the density by Green's function to find the k-space potential and 
    // invert for the real potenital. Second, normalize the system and, finally,
    // put the potential information into the grid.
    for (i = 0; i < (Ng)*(Ng)*(Ng); ++i) {
        // Convolve the potential
        potential[i] = green[i]*potential[i];
    }
    // Inversion
    fftw_execute(fft_potential);
    // Normalization
    for (i = 0; i < (Ng)*(Ng)*(Ng); ++i) {
        potential[i] = potential[i]/n;
    }
    l = 0;
    for (i = 0; i < Ng; ++i) {
        for (j = 0; j < Ng; ++j) {
            for (k = 0; k < Ng; ++k) {
                // Fix the grid info
                gal->potential[i][j][k] = -4.0*pi*potential[l]*G*unit_mass/
					  kpc;
                ++l;
            }
        }
    }
     
    // Print the potential in the xy-plane for z = 0 if the option is set.
    fprintf(stderr,"Grid cell spacings (kpc): x = %f y = %f z = %f\n",
           space_x,space_y,space_z);
    if (dump != 0) {
       for (i = Ng/4; i < 3*Ng/4; ++i) {
           for (j = Ng/4; j < 3*Ng/4; ++j) {
               printf("%lf %lf %le \n",(double) (i-Ng/2)*space_x,(double) 
                     (j-Ng/2)*space_y,gal->potential[i][j][Ng/2]);
           }
           printf("\n");
       }
    }

    // Free fftw plan.
    fftw_destroy_plan(fft_potential);
    fftw_destroy_plan(fft_green);

    // Kill the storage arrays since they are no longer needed.
    fftw_free(green); fftw_free(potential);
    for (i = 0; i < Ng; ++i) {
        for (j = 0; j < Ng; ++j) {
            free(green_grid[i][j]);
            }
        free(green_grid[i]); 
    }
    free(green_grid);

    #ifdef OPENMP_CFLAGS
    fftw_cleanup_threads();
    #endif

    return 0;
}

// The truncation function defined by Hockney and Eastwood for
// discretizing Green's function.
static double truncation_func(double x, double y, double z) {

    double tx, ty, tz;

    // Make sure to take the absolute value
    x = sqrt(x*x);
    y = sqrt(y*y);
    z = sqrt(z*z);
    // Truncation in the x direction
    tx = 1.0;
    if (x > 0.5) {
       tx = 0.0;
    } else if (x == 0.5) {
       tx = 0.5;
    }
    // Truncation in the y direction
    ty = 1.0;
    if (y > 0.5) {
       ty = 0.0;
    } else if (y == 0.5) {
       ty = 0.5;
    }
    // Truncation in the z direction
    tz = 1.0;
    if (z > 0.5) {
       tz = 0.0;
    } else if (z == 0.5) {
       tz = 0.5;
    }

    return tx*ty*tz;
}

// Destroy a galaxy
void destroy_galaxy(galaxy *gal, int info) {

     int i,j;

     if (info != 0) {
        fprintf(stderr,"Destroying galaxy...");
     }

     // Deallocate the potential grid to be really nice to the memory.
     for (i = 0; i < Ng; ++i) {
         for (j = 0; j < Ng; ++j) {
             free(gal->potential[i][j]);
         }
         free(gal->potential[i]);
     }
     free(gal->potential);

     // Deallocate all the small parts of the galaxy to be nice to the memory.
     free(gal->id); free(gal->x); free(gal->y); free(gal->z); free(gal->mass);
     free(gal->vel_x); free(gal->vel_y); free(gal->vel_z);

     if (info != 0) {
        fprintf(stderr," All memory deallocated.\n");
     }
     
     return;
}

// Destroy a system of galaxies.
void destroy_galaxy_system(galaxy *gal, int info) {

     if (info != 0) {
        fprintf(stderr,"Destroying galaxy system...");
     }

     free(gal->id); free(gal->x); free(gal->y); free(gal->z); free(gal->mass);
     free(gal->vel_x); free(gal->vel_y); free(gal->vel_z);

     if (info != 0) {
        fprintf(stderr," All memory deallocated.\n");
     }
     
     return;
}

// In order to perform a galaxy collision, it is necessary to combine 
// two galaxies into one set of orbiting galaxies. The orbits should 
// be defined by a user-defined orbit or the set_orbit_parabolic()
// function and the following function should be called to combine 
// the galaxies into one object.
// 
// It is necessary to allocate the different parts of the memory here
// because this galactic system does not have quantities like halo
// or disk scale lengths. These are quantities in the parent galaxy,
// not the galaxy-orbit-galaxy combination.
int create_galaxy_system(galaxy *gal_1, galaxy *gal_2, galaxy *gal_3) {

     int i, a, num_part;

     // Create the galaxy system first!
     gal_3->num_part[0] = gal_1->num_part[0] + gal_2->num_part[0];
     gal_3->num_part[1] = gal_1->num_part[1] + gal_2->num_part[1];
     gal_3->num_part[2] = gal_3->num_part[0] + gal_3->num_part[1];
     num_part = gal_3->num_part[2];
     
     // Allocate particle id numbers array.
     if (!(gal_3->id=calloc(num_part,sizeof(int)))) {
        printf("Unable to allocate particle ID numbers.\n");
        return -1;
     }
     // Allocate x coordinates for all the particles.
     if (!(gal_3->x=calloc(num_part,sizeof(double)))) {
        printf("Unable to allocate particle x coordinates.\n");
        return -1;
     }
     // Allocate y coordinates for all the particles.
     if (!(gal_3->y=calloc(num_part,sizeof(double)))) {
        printf("Unable to allocate particle y coordinates.\n");
        return -1;
     }
     // Allocate z coordinates for all the particles.
     if (!(gal_3->z=calloc(num_part,sizeof(double)))) {
        printf("Unable to allocate particle z coordinates.\n");
        return -1;
     }
     // Allocate x velocities for all the particles.
     if (!(gal_3->vel_x=calloc(num_part,sizeof(double)))) {
        printf("Unable to allocate particle x coordinates.\n");
        return -1;
     }
     // Allocate y velocities for all the particles.
     if (!(gal_3->vel_y=calloc(num_part,sizeof(double)))) {
        printf("Unable to allocate particle y coordinates.\n");
        return -1;
     }
     // Allocate z velocities for all the particles.
     if (!(gal_3->vel_z=calloc(num_part,sizeof(double)))) {
        printf("Unable to allocate particle z coordinates.\n");
        return -1;
     }
     // Allocate masses for all the particles.
     if (!(gal_3->mass=calloc(num_part,sizeof(double)))) {
        printf("Unable to allocate particle masses.\n");
        return -1;
     }
     // Turn off the galaxy potential.
     gal_3->potential_defined = 0;

     // Copy all of the galaxy information.
     // Disk 1
     a = 0;
     for (i = 0; i < gal_1->num_part[0]; ++i) {
         gal_3->mass[a] = gal_1->mass[i];
         gal_3->id[a] = i;
         gal_3->x[a] = gal_1->x[i];
         gal_3->y[a] = gal_1->y[i];
         gal_3->z[a] = gal_1->z[i];
         gal_3->vel_x[a] = gal_1->vel_x[i];
         gal_3->vel_y[a] = gal_1->vel_y[i];
         gal_3->vel_z[a] = gal_1->vel_z[i];
         ++a;
     }
     // Disk 2
     for (i = 0; i < gal_2->num_part[0]; ++i) {
         gal_3->mass[a] = gal_2->mass[i];
         gal_3->id[a] = a;
         gal_3->x[a] = gal_2->x[i];
         gal_3->y[a] = gal_2->y[i];
         gal_3->z[a] = gal_2->z[i];
         gal_3->vel_x[a] = gal_2->vel_x[i];
         gal_3->vel_y[a] = gal_2->vel_y[i];
         gal_3->vel_z[a] = gal_2->vel_z[i];
         ++a;
     }
     // Halo 1
     for (i = gal_1->num_part[0]; i < gal_1->num_part[2]; ++i) {
         gal_3->mass[a] = gal_1->mass[i];
         gal_3->id[a] = a;
         gal_3->x[a] = gal_1->x[i];
         gal_3->y[a] = gal_1->y[i];
         gal_3->z[a] = gal_1->z[i];
         gal_3->vel_x[a] = gal_1->vel_x[i];
         gal_3->vel_y[a] = gal_1->vel_y[i];
         gal_3->vel_z[a] = gal_1->vel_z[i];
         ++a;
     }
     // Halo 2
     for (i = gal_2->num_part[0]; i < gal_2->num_part[2]; ++i) {
         gal_3->mass[a] = gal_2->mass[i];
         gal_3->id[a] = a;
         gal_3->x[a] = gal_2->x[i];
         gal_3->y[a] = gal_2->y[i];
         gal_3->z[a] = gal_2->z[i];
         gal_3->vel_x[a] = gal_2->vel_x[i];
         gal_3->vel_y[a] = gal_2->vel_y[i];
         gal_3->vel_z[a] = gal_2->vel_z[i];
         ++a;
     }

    return 0;
}

