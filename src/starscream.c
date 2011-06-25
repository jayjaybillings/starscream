/*-----------------------------------------------------------------------------
/
/ Filename: starscream.c
/ Author: Jay Billings
/ Author's email: jayjaybillings@gmail.com
/ Description: Starscream creates galaxies. 
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
/ Date: 2008/04/15
/
*///---------------------------------------------------------------------------

/* The following program is an example of how to use Starscream to create initial
   conditions for a galaxy.                                                           */

// Be sure to include the starscream header file, wherever it is located.
#include "starscream.h"

// The All-Powerful Main!
int main (int argc, char **argv) {

     int i, j, k, parts[2], N;
     double m_d, j_d, lambda, c, v200, space;
     // Galaxy pointers
     galaxy *galaxy_1,*galaxy_2,*galaxy_3, galaxy_4;

     // Set some important variables: The number of disk particles and
     // the number of halo particles, the disk mass fraction, the disk
     // angular momentum fraction, the spin parameter, the halo 
     // concentration, the virial velocity, the potential grid size,
     // and the grid spacing.
     parts[0] = 10000;
     parts[1] = 20000;
     m_d = 0.025;
     j_d = m_d;
     lambda = 0.050;
     c = 15.0;
     v200 = 1.6E7;
     N = 128;
     space = 2.0;

     // Since this example is using pointers to galaxies, allocate them.
     // Alternatively, you could declare a non-pointer galaxy and pass its
     // memory address.
     if (!(galaxy_1=calloc(1,sizeof(galaxy)))) {
        fprintf(stderr,"Unable to allocate galaxy. Aborting.\n");
        return 0;
     }
     if (!(galaxy_2=calloc(1,sizeof(galaxy)))) {
        fprintf(stderr,"Unable to allocate galaxy. Aborting.\n");
        return 0;
     }
     if (!(galaxy_3=calloc(1,sizeof(galaxy)))) {
        fprintf(stderr,"Unable to allocate galaxy. Aborting.\n");
        return 0;
     }


     // Provide values for the free parameters of the galaxy. If the galaxy can
     // not be created, return an error.
     if ((i = create_galaxy(galaxy_1,parts,m_d,j_d,lambda,c,v200,N,space,1)) 
        != 0) {
        fprintf(stderr,"Unable to create galaxy. Aborting.\n");
        exit(0);
     }

     // Set up the particles positions, disk potential, and particle velocities of 
     // the particles in the galaxy. The option on set_galaxy_velocity tells the 
     // function to use dispersion.
     if ((i = set_galaxy_coords(galaxy_1)) != 0) {
        fprintf(stderr,"Unable to set coordinates. Aborting.\n");
        return 0;
     }
     set_galaxy_potential(galaxy_1,0);
     if ((i = set_galaxy_velocity(galaxy_1,1)) != 0) {
        printf("Unable to set velocities. Aborting.\n");
        return 0;
     }
     fprintf(stderr,"Q_min = %f\n",Q_min);
     
     // There are a set of functions for writing physical information about the galaxy.
     //write_galaxy_position_disk(galaxy_1);
     //write_galaxy_velocity(galaxy_1,0);
     //write_galaxy_potential(galaxy_1);
     //write_galaxy_rotation_curve(galaxy_1);
   
     /*for (i = -galaxy_1->r200; i < galaxy_1->r200; i=i+4) {
         for (j = -galaxy_1->r200; j < galaxy_1->r200; j=j+4) {
             for (k = -galaxy_1->r200; k < galaxy_1->r200; k=k+4) {
               fprintf(stdout,"%f %f %f %le\n",
                      (float) i, (float) j, (float) k,
                      potential_func(galaxy_1,(double) i,(double) j, 0.0));
             }
          }
      }*/
 
     // Load a Gadget2 snapshot and read it into a galaxy, halo first and then the disk.
     /*load_snapshot("snapshot_080",1);
     i = 1; 
     for (j = galaxy_1->num_part[0]; j < galaxy_1->num_part[2]; ++j) {
         galaxy_1->x[j] = P[i].Pos[0];
         galaxy_1->y[j] = P[i].Pos[1];
         galaxy_1->z[j] = P[i].Pos[2];
         galaxy_1->vel_x[j] = P[i].Vel[0];
         galaxy_1->vel_y[j] = P[i].Vel[1];
         galaxy_1->vel_z[j] = P[i].Vel[2];
         galaxy_1->mass[j] = P[i].Mass*unit_mass;
         ++i;
     }
     for (j = 0; j < galaxy_1->num_part[0]; ++j) {
         galaxy_1->x[j] = P[i].Pos[0];
         galaxy_1->y[j] = P[i].Pos[1];
         galaxy_1->z[j] = P[i].Pos[2];
         galaxy_1->vel_x[j] = P[i].Vel[0];
         galaxy_1->vel_y[j] = P[i].Vel[1];
         galaxy_1->vel_z[j] = P[i].Vel[2];
         galaxy_1->mass[j] = P[i].Mass*unit_mass;
         ++i;
     }
     // It is important to always unload the snapshot as soon as you use it.
     unload_snapshot();*/
     
     // Create a second galaxy by copying the first.
     copy_galaxy(galaxy_1,galaxy_2,0);

     // This function sets two galaxies on a parabolic orbit.
     set_orbit_parabolic(galaxy_1,galaxy_2,galaxy_1->r200,3.5);

     // Now that the galaxies have been put on an orbit with each other,
     // store them in the same galaxy object.
     create_galaxy_system(galaxy_1,galaxy_2,galaxy_3);

     // Write initial conditions for the massively parallel N-body code Gadget2.
     // If you do not use Gadget, do not use this function -- the output is in binary!
     //write_gadget_ics(galaxy_1,"initial_conditions.dat");
     //write_gadget_ics(galaxy_2,"initial_conditions.dat");
     write_gadget_ics(galaxy_3,"initial_conditions.dat");

     // Destroy a galaxy. If the galaxy can not be destroyed, return an error. This
     // function will SEGFAULT if the arrays in the galaxy can not be freed.
     destroy_galaxy(galaxy_1,0);
     destroy_galaxy(galaxy_2,0);
     destroy_galaxy_system(galaxy_3,0);

     // Destroy the random number environment.
     gsl_rng_free(r);

     // Print an exit statement and call it quits.
     fprintf(stderr,"\nDone.\n");
     return 0;
}

