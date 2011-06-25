/*-----------------------------------------------------------------------------
/
/ Filename: starscream.h
/ Author: Jay Billings
/ Author's email: jayjaybillings@gmail.com
/ Description: This is the header file for Starscream. Starscream creates
/              galaxies.
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

/* Header files to include							      */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <fftw3.h>

/* This is a type definition of a galaxy. The thought here is to create galaxies
as objects and, hopefully, make the code extremely clean. It is essentially a
a collection of arrays.
        lambda: The spin paramater.
        m_d: The mass fraction of the disk relative to the total mass.
        j_d: The angular momentum fraction of the disk relative to the total.
        total_mass: The total mass of the galaxy.
        disk_mass: The mass of the disk.
        disk_scale_length: The scale length of the disk.
        halo_mass: The mass of the halo.
        halo_scale_length: The scale length of the halo.
        halo_concentration: The concentration of the halo. 
        halo_a_value: The value that relates a Hernquist profile to a NFW profile
                      through the halo scale length and halo concentration.
        v200: The halo's virial velocity.
        r200: The halo's virial radius.
	x: The x positions of the particle.
        y: The y positions of the particle.
        z: The z positions of the particle.
        vel_x, vel_y, vel_z: The velocity components of the particles.
        mass: The masses of the particle.					
        potential: The potential of the disk evaluated at the grid points. The
		   size of this array should be gx*gy*gz.	      
        storage1: A special variable for storing important info. Could be anything!   
	id: The particle id number. 
	num_part: An integer array containing the number of particles in the disk, [0],
		  and in the halo, [1], and the total number of particles (optional). */
typedef struct {
    double lambda;
    double m_d;
    double j_d;
    double total_mass;
    double disk_mass;
    double disk_scale_length;
    double halo_mass;
    double halo_scale_length;
    double halo_concentration;
    double halo_a_value;
    double v200;
    double r200;
    double *x;
    double *y;
    double *z;
    double *vel_x;
    double *vel_y;
    double *vel_z;
    double *mass;
    double ***potential;
    double space[3];
    double *storage1;
    int *id;
    int num_part[3];
    int potential_defined;
} galaxy;

//Gadget2-style header for Gadget2 snapshots.
struct io_header_1 {
   int npart[6];
   double mass[6];
   double time;
   double redshift;
   int flag_sfr;
   int flag_feedback;
   int npartTotal[6];
   int flag_cooling;
   int num_files;
   double BoxSize;
   double Omega0;
   double OmegaLambda;
   double HubbleParam;
   char fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8];
} header1;

//Gadget2-style particle data structure.
struct particle_data {
   float Pos[3];
   float Vel[3];
   float Mass;
   int Type;

   float Rho, U, Temp, Ne;
} *P;

//Particle IDs.
int *Id;

//Time and Redshift variables. Both will probably be 0.
double Time, Redshift;

//Good numbers to have on file: Total number of particles and total number of gas
//particles.
int NumPart, Ngas;

// Global variables for the random number environment.
const gsl_rng_type *T;
gsl_rng *r;
int random_number_set;

// A tolerance parameter for the GSL QAG integration.
#define epsabs 1.0E-10
#define epsrel 1.0E-7
#define key 6

// A variable to check if the galaxy is being copied.
int a_copy;

// A list of universal constants. They are set in create_galaxy(). 
double pi, G, H, unit_mass, kpc, crit_dens;

// The particle-mesh grid size and the Green's function and potential 
// storage buffers. Global to keep from hitting the stack limit for large grid 
// size.
int Ng; 
double *green, ***green_grid, *potential;

// A grid for the velocity dispersion.
double **v2rz;

// A global variable for storing the minimum local stability parameter, the
// Toomre Q value.
double Q_min;

/*----- These are the function prototypes for Starscream. 		-----*/

// Initialization and destruction functions
int create_galaxy(galaxy *, int [], double, double, double, double, 
                 double, int, double, int);
void allocate_galaxy_storage_variable(galaxy *, int);
void allocate_dispersion_grid();
void destroy_dispersion_grid();
int set_galaxy_coords(galaxy *);
int set_galaxy_velocity(galaxy *, int);
static int set_dispersion_grid(galaxy *, double, double disp_func(galaxy *));
double get_dispersion(double, double, double);
int set_galaxy_potential(galaxy *, int);
static double truncation_func(double, double, double);
void destroy_galaxy(galaxy *, int);
void destroy_galaxy_system(galaxy *, int);
int create_galaxy_system(galaxy *, galaxy *, galaxy *); 

// Structure functions
double f_exp_disk(double,double); 
double f_hernquist_halo(double,double);
double density_disk_func(galaxy *, double, double);
double density_disk_cyl_func(galaxy *, double, double);
double density_halo_func(galaxy *, double);
double density_halo_cyl_func(galaxy *, double, double);
double disk_scale_length_func(galaxy *);
double halo_mass_func(galaxy *,double);
double disk_mass_func(galaxy *,double);
double j_d_func(galaxy *);
static double dj_d_func(double, void *);
double f_c_func(double);
double g_c_func(double);
static double dg_c_func(double, void *);
double f_s_func(double, double);

// Velocity functions
double v_c_func(galaxy *,double);
double v_c_disk_func(galaxy *,double);
double v_c_halo_func(galaxy *,double);
double v2a_z_disk_func(galaxy *);
double v2a_z_halo_func(galaxy *);
static double dv2a_z_halo_func(double, void *);
static double dv2a_z_disk_func(double, void *);
double v2a_theta_disk_func(galaxy *, double, double, double);
double v2a_theta_halo_func(galaxy *, double);
static double vtheta_force_func(galaxy *, double);
static double vrz_disk_force_func(galaxy *, double);
static double vrz_force_func(galaxy *, double, double);
double disk_potential_wrapper_func(double, void *);

// Potential and force functions
double disk_potential_func(galaxy *, double, double, double);
double disk_potential_nbody_func(galaxy *, double, double, double);
double halo_potential_func(galaxy *, double, double, double);
double halo_potential_nbody_func(galaxy *, double, double, double);
double potential_func(galaxy *, double, double, double);
double disk_force_nbody_func(galaxy *, double, double, double);
double halo_force_func(galaxy *, double);
double halo_force_cyl_func(galaxy *, double, double);
double halo_force_dp_func(galaxy *, double);
double force_func(galaxy *, double);

// Input, output, and manipulation functions
void write_galaxy_position(galaxy *,int);
void write_galaxy_position_disk(galaxy *);
void write_galaxy_position_halo(galaxy *);
void write_galaxy_velocity(galaxy *,int);
void write_galaxy_velocity_disk(galaxy *);
void write_galaxy_velocity_halo(galaxy *);
void write_galaxy_potential(galaxy *);
void write_galaxy_rotation_curve(galaxy *);
int write_gadget_ics(galaxy *, char *);
void copy_galaxy(galaxy *, galaxy *, int);
void set_orbit_parabolic(galaxy *, galaxy *, double, double);
int load_snapshot(char *, int);
int allocate_memory();
int reordering();
int unload_snapshot();
