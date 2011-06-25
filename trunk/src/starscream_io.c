/*-----------------------------------------------------------------------------
/
/ Filename: starscream_io.c
/ Author: Jay Billings
/ Author's email: jayjaybillings@gmail.com
/ Description: These are the input, output, and manipulation routines for
/              Starscream. "Manipulation" of a galaxy includes and action
/              that operates on the galaxy as a whole, such as placing the
/              entire galaxy on an orbit or adding it into a structure 
/              with another galaxy. Starscream creates galaxies. 
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

// The next three functions write the position coordinates of a galaxy to the 
// screen in xyz format.
void write_galaxy_position(galaxy *gal,int info) {

     int i;
     
     if (info != 0) {
        fprintf(stderr,"Writing galaxy positions to standard out...\n");
        printf("#Disk: \n");
     }
     
     //Print disk positions
     if (gal->num_part[0] > 0) {
        for (i = 0; i < gal->num_part[0]; ++i) {
            printf("%le %le %le \n",gal->x[i],gal->y[i],gal->z[i]);
        }
     }

     //Print halo positions
     if (info != 0) {
        printf("#Halo: \n");
     }
     if (gal->num_part[1] != 0) {    
        for (i = gal->num_part[0]; i < gal->num_part[2] ; ++i) {
            printf("%le %le %le \n",gal->x[i],gal->y[i],gal->z[i]);
        }
     }
     return;
}

void write_galaxy_position_disk(galaxy *gal) {

     int i;

     fprintf(stderr,"Writing disk positions to standard out...\n");

     printf("#Disk: \n");
     //Print disk positions
     if (gal->num_part[0] > 0) {
        for (i = 0; i < gal->num_part[0]; ++i) {
            printf("%le %le %le \n",gal->x[i],gal->y[i],gal->z[i]);
        }
     }

     return;
}

void write_galaxy_position_halo(galaxy *gal) {

     int i;

     fprintf(stderr,"Writing halo positions to standard out...\n");

     //Print halo positions
     printf("#Halo: \n");
     if (gal->num_part[1] != 0) {    
        for (i = gal->num_part[0]; i < gal->num_part[2] ; ++i) {
            printf("%le %le %le \n",gal->x[i],gal->y[i],gal->z[i]);
        }
     }

     return;
}

// The next three functions write the velocity coordinates of a galaxy to 
// the screen in xyz format.
void write_galaxy_velocity(galaxy *gal, int info) {

     int i;

     //Print disk positions
     if (info != 0) {
        fprintf(stderr,"Writing galaxy velocities to standard out...\n");
        printf("#Disk: \n");
     }
     if (gal->num_part[0] > 0) {
        for (i = 0; i < gal->num_part[0]; ++i) {
            printf("%le %le %le \n",gal->vel_x[i],gal->vel_y[i],gal->vel_z[i]);
        }
     }

     //Print halo positions
     if (info != 0) {
        printf("#Halo: \n");
     }
     if (gal->num_part[1] != 0) {    
        for (i = gal->num_part[0]; i < gal->num_part[2] ; ++i) {
            printf("%lf %lf %lf \n",gal->vel_x[i],gal->vel_y[i],gal->vel_z[i]);
        }
     }

     return;
}

void write_galaxy_velocity_disk(galaxy *gal) {

     int i;

     fprintf(stderr,"Writing disk velocities to standard out...\n");
     printf("#Disk: \n");
     if (gal->num_part[0] > 0) {
        for (i = 0; i < gal->num_part[0]; ++i) {
            printf("%le %le %le \n",gal->vel_x[i],gal->vel_y[i],gal->vel_z[i]);
        }
     }

     return;
}

void write_galaxy_velocity_halo(galaxy *gal) {

     int i;

     fprintf(stderr,"Writing halo velocities to standard out...\n");
     printf("#Halo: \n");
     if (gal->num_part[1] != 0) {    
        for (i = gal->num_part[0]; i < gal->num_part[2] ; ++i) {
            printf("%lf %lf %lf \n",gal->vel_x[i],gal->vel_y[i],gal->vel_z[i]);
        }
     }

     return;
}

// Function to write initial conditions to file in default Gadget2 format.
// The code was originally an input routine, read_snapshot.c, provided by
// Volker Springel with the Gadget2 source code. It has been hacked into a 
// write routine.
int write_gadget_ics(galaxy *gal, char *fname) {

    FILE *fp1, *fp2;
    int i, j, k, dummy, ntot_withmasses;
    int t,n,off,pc,pc_new,pc_sph;
    int files = 1;
    int Ids[gal->num_part[2]];
    char buf[200];
    #define SKIP2 fwrite(&dummy, sizeof(dummy), 1, fp1);

    NumPart = gal->num_part[2];
    // Set everything to zero and overwrite it later if needed.
    header1.npart[0] = 0;
    header1.npart[1] = 0;
    header1.npart[2] = 0;
    header1.npart[3] = 0;
    header1.npart[4] = 0;
    header1.npart[5] = 0;
    header1.npartTotal[0] = 0;
    header1.npartTotal[1] = 0;
    header1.npartTotal[2] = 0;
    header1.npartTotal[3] = 0;
    header1.npartTotal[4] = 0;
    header1.npartTotal[5] = 0;
    header1.mass[0] = 0.0;
    header1.mass[1] = 0.0;
    header1.mass[2] = 0.0;
    header1.mass[3] = 0.0;
    header1.mass[4] = 0.0;
    header1.mass[5] = 0.0;
    
    // Set the header values to some defaults.
    header1.npart[1] = gal->num_part[1];
    header1.npart[2] = gal->num_part[0];
    header1.npartTotal[1] = gal->num_part[1];
    header1.npartTotal[2] = gal->num_part[0];
    header1.time = 0.0;
    header1.redshift = 0.0;
    header1.flag_sfr = 0.0;
    header1.flag_feedback = 0.0;
    header1.flag_cooling = 0.0;
    header1.num_files = 1; 
    header1.BoxSize = 0.0;
    header1.Omega0 = 0.0;
    header1.OmegaLambda = 0.0;
    header1.HubbleParam = 1.0;

    if (!(P=malloc(NumPart*sizeof(struct particle_data)))) {
       fprintf(stderr,"Unable to create particle data structure in memory.");
       exit(0);
    }
    P--;

    // Transfer the particle data from the Starscream galaxy data structure to the Gadget
    // position_data structure.
    j = 1;
    while(j < gal->num_part[2]) {
        for (i = gal->num_part[0]; i < gal->num_part[2]; ++i) {
            P[j].Pos[0] = gal->x[i];
            P[j].Pos[1] = gal->y[i];
            P[j].Pos[2] = gal->z[i];
            P[j].Vel[0] = gal->vel_x[i];
            P[j].Vel[1] = gal->vel_y[i];
            P[j].Vel[2] = gal->vel_z[i];
            P[j].Mass = gal->mass[i]/unit_mass;
            P[j].Type = 1;
            ++j;
        }
        for (i = 0; i < gal->num_part[0]; ++i) {
            P[j].Pos[0] = gal->x[i];
            P[j].Pos[1] = gal->y[i];
            P[j].Pos[2] = gal->z[i];
            P[j].Vel[0] = gal->vel_x[i];
            P[j].Vel[1] = gal->vel_y[i];
            P[j].Vel[2] = gal->vel_z[i];
            P[j].Mass = gal->mass[i]/unit_mass;
            P[j].Type = 2;
            ++j;
        }
    }

    fprintf(stderr,"Writing initial conditions... \n");
   
    for(i=0, pc=1; i<files; i++, pc=pc_new)
      {
        if(files>1)
	  sprintf(buf,"%s.%d",fname,i);
        else
	  sprintf(buf,"%s",fname);

      if(!(fp1=fopen(buf,"w")))
	{
	  fprintf(stderr,"can't open file `%s`\n",buf);
	  exit(0);
	}

      fflush(stdout);

      dummy = sizeof(header1);
      fwrite(&dummy, sizeof(dummy), 1, fp1);
      fwrite(&header1, sizeof(header1), 1, fp1);
      fwrite(&dummy, sizeof(dummy), 1, fp1);

      for(k=0, ntot_withmasses=0; k<6; k++)
	{
	  if(header1.mass[k]==0)
	    ntot_withmasses+= header1.npart[k];
	}

      dummy = 3*sizeof(float)*NumPart;
      SKIP2;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fwrite(&P[pc_new].Pos[0], sizeof(float), 3, fp1);
	      pc_new++;
	    }
	}
      SKIP2;

      SKIP2;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fwrite(&P[pc_new].Vel[0], sizeof(float), 3, fp1);
	      pc_new++;
	    }
	}
      SKIP2;
    

      dummy = sizeof(int)*NumPart;
      SKIP2;

      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
              Ids[pc_new] = pc_new;
	      fwrite(&Ids[pc_new], sizeof(int), 1, fp1);
	      pc_new++;
	    }
	}
      SKIP2;


      if(ntot_withmasses>0) {
        dummy = sizeof(float)*NumPart;
	SKIP2;
      }
      for(k=0, pc_new=pc; k<6; k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      P[pc_new].Type=k;

	      if(header1.mass[k]==0)
		fwrite(&P[pc_new].Mass, sizeof(float), 1, fp1);
	      else
		P[pc_new].Mass= header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses>0) {
	SKIP2;
      }

      if(header1.npart[0]>0)
	{
          dummy = sizeof(float)*header1.npart[0];
	  SKIP2;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fwrite(&P[pc_sph].U, sizeof(float), 1, fp1);
	      pc_sph++;
	    }
	  SKIP2;

/*	  SKIP2;
          printf("Writing particle density data with buffer size = %d \n",dummy);
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fwrite(&P[pc_sph].Rho, sizeof(float), 1, fp1);
	      pc_sph++;
	    }
	  SKIP2;

	  if(header1.flag_cooling)
	    {
	      SKIP2;
              printf("Writing particle smoothing data with buffer size = %d \n",dummy);
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		  fwrite(&P[pc_sph].Ne, sizeof(float), 1, fp1);
		  pc_sph++;
		}
	      SKIP2;
	    }*/
	}

    }
    
    fprintf(stderr,"Initial Conditions Written. \n");
    P++; free(P);
    fclose(fp1);
    return 0;
}

// This function calculates and writes the rotation curve for a particular galaxy.
// The circular velocity is calculated out to the virial radius.
//
// Unlike some of the other output functions, this function writes the output to
// a set of files -- rcurve.dat, rcurve_disk.dat, rcurve_halo.dat -- These files
// can be plotted in gnuplot with the command
//
//    plot 'rcurve.dat','rcurve_disk.dat','rcurve_halo.dat'
//
// to produce a cumulative rotation curve and its parts on the same graph.
void write_galaxy_rotation_curve(galaxy *gal) {

    int i;
    double radius, theta, v_c;
    char filename[200];
    FILE *fp1;

    // Total rotation curve
    sprintf(filename, "rcurve.dat");
    fp1 = fopen(filename, "w");
    if (fp1 == NULL) {
       printf("Could not open rcurve.dat for writing the rotation curve. Aborting.\n");
    }
    for (i = 1; i < (int) gal->r200; ++i) {
        radius = (double) i;
        v_c = v_c_func(gal,radius)/1.0E5;
        // Write the radius and circular velocity to file in kpc and km/s respectively.
        fprintf(fp1,"%lf %lf\n",radius,v_c);
    }
    fclose(fp1);

    // Disk rotation curve
    sprintf(filename, "rcurve_disk.dat");
    fp1 = fopen(filename, "w");
    if (fp1 == NULL) {
       printf("Could not open rcurve_disk.dat for writing the rotation curve.");
       printf("Aborting.\n");
    }
    for (i = 1; i < (int) gal->r200; ++i) {
        radius = (double) i;
        v_c = v_c_disk_func(gal,radius)/1.0E5;
        // Write the radius and circular velocity to file in kpc and km/s respectively.
        fprintf(fp1,"%lf %lf\n",radius,v_c);
    }
    fclose(fp1);

    // Halo rotation curve
    sprintf(filename, "rcurve_halo.dat");
    fp1 = fopen(filename, "w");
    if (fp1 == NULL) {
       fprintf(stderr,"Could not open rcurve_halo.dat for writing the");
       fprintf(stderr,"rotation curve. Aborting.\n");
    }
    for (i = 0; i < (int) gal->r200; ++i) {
        radius = (double) i;
        v_c = v_c_halo_func(gal,radius)/1.0E5;
        // Write the radius and circular velocity to file in kpc and km/s respectively.
        fprintf(fp1,"%lf %lf\n",radius,v_c);
    }
    fclose(fp1);

    return;
}

// This function writes the potential for a particular galaxy. The function 
// calculates the potentials of the disk and halo and their sum, (the total
// gravitational potential), and outputs in an nxy format to a file. This
// data can be plotted with xmgrace using a command like
//
//    xmgrace -nxy potential.dat
//
// or with gnuplot using splot and a script like
//
//    splot 'potential.dat' using 1:2:3, \
//         'potential.dat' using 1:2:4, \
//         'potential.dat' using 1:2:5
//
// to produce a detailed picture of the potential structure of the system.
// The user may want to consider adding "w lines" to that script since 
// there are a lot of points in the data file.
//
// Please note that the potential of the disk is calculated only in the x,y
// plane. The potentials are calculated to the +/- the virial radius.
void write_galaxy_potential(galaxy *gal) {
   
    int i,j;
    int r200 = (int) gal->r200;
    char filename[200];
    FILE *fp1;

    if (gal->num_part[2] > 0) {
       sprintf(filename, "potential.dat");
       fp1 = fopen(filename, "w");
       if (fp1 == NULL) {
          fprintf(stderr,"Could not open potential.dat for writing the");
          fprintf(stderr," potential. Skipping it...\n");
          return;
       }
       for (i = -r200; i < r200; i=i+4) {
           for (j = -r200; j < r200; j=j+4) {
               fprintf(fp1,"%d %d %le",i,j,disk_potential_func(gal,(double) i,
                      (double) j,0.0));
               fprintf(fp1," %le",halo_potential_func(gal,(double) i,
                      (double) j,0.0));
               fprintf(fp1," %le\n",potential_func(gal,(double) i,(double) j,
                      0.0));
           }
       }
       fclose(fp1);
    } else {
       fprintf(stderr,"There are no particles in the galaxy!\n");
    }

    return;
}

// This function copies one galaxy to another.
void copy_galaxy(galaxy *gal_1, galaxy *gal_2, int info) {

     int i, parts[2];
     double c;

     // Create the galaxy first!
     parts[0] = gal_1->num_part[0];
     parts[1] = gal_1->num_part[1];
     c = gal_1->halo_concentration;
     a_copy = 1;
     create_galaxy(gal_2,parts,gal_1->m_d,gal_1->j_d,gal_1->lambda,c,
                  gal_1->v200,Ng/2,gal_1->space[0],info);
     gal_2->disk_scale_length = gal_1->disk_scale_length;
     a_copy = 0;

     // Copy all the coordinate information.
     for (i = 0; i < gal_1->num_part[2]; ++i) {
         gal_2->x[i] = gal_1->x[i];
         gal_2->y[i] = gal_1->y[i];
         gal_2->z[i] = gal_1->z[i];
         gal_2->vel_x[i] = gal_1->vel_x[i];
         gal_2->vel_y[i] = gal_1->vel_y[i];
         gal_2->vel_z[i] = gal_1->vel_z[i];
         gal_2->mass[i] = gal_1->mass[i];
     }

     return;
}

// This function places two galaxies on a parabolic orbit in the x,y 
// plane. For a detailed discussion of parabolic orbits, see Elements 
// of Astromechanics, Peter van de Kamp, 1962.
void set_orbit_parabolic(galaxy *gal_1, galaxy *gal_2, double radius, double min_r) {
    
     int i;
     double nu, x_o, y_o, z_o, mu;
    
     nu = acos((min_r/radius) - 1.0);
     x_o = radius*cos(nu);
     y_o = radius*sin(nu);
     z_o = 0.0;
     // Mu is the standard gravitational parameter, calculated here in 
     // terms of the reduced mass.
     mu = G*(gal_1->total_mass*gal_2->total_mass) /
          (gal_1->total_mass + gal_2->total_mass);

     // Galaxy 1, displaced by adding x_o and y_o.
     for (i = 0; i < gal_1->num_part[2]; ++i) {
         gal_1->x[i] += x_o;
         gal_1->y[i] += y_o;
         gal_1->vel_x[i] = -gal_1->vel_x[i] + sqrt(mu/(2.0*min_r*kpc))*
                           sin(nu)/1.0E5;
         gal_1->vel_y[i] = -gal_1->vel_y[i] - sqrt(mu/(2.0*min_r*kpc))*
                           (1.0+cos(nu))/1.0E5;
     }
     // Galaxy 2, displaced by subtracting x_o and y_o.
     for (i = 0; i < gal_2->num_part[2]; ++i) {
         gal_2->x[i] -= x_o;
         gal_2->y[i] -= y_o;
         gal_2->vel_x[i] = -gal_2->vel_x[i] - sqrt(mu/(2.0*min_r*kpc))*
                           sin(nu)/1.0E5;
         gal_2->vel_y[i] = -gal_2->vel_y[i] + sqrt(mu/(2.0*min_r*kpc))*
                           (1.0+cos(nu))/1.0E5;
     }

     return;
}

/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 *
 * Author: Volker Springel, Max Planck Institute
 */
int load_snapshot(char *fname, int files)
{
  FILE *fd;
  char   buf[200];
  int    i,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc,pc_new,pc_sph;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
      if(files>1)
	sprintf(buf,"%s.%d",fname,i);
      else
	sprintf(buf,"%s",fname);

      if(!(fd=fopen(buf,"r")))
	{
	  printf("can't open file `%s`\n",buf);
	  exit(0);
	}

      printf("reading `%s' ...\n",buf); fflush(stdout);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      if(files==1)
	{
	  for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	    NumPart+= header1.npart[k];
	  Ngas= header1.npart[0];
	}
      else
	{
	  for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	    NumPart+= header1.npartTotal[k];
	  Ngas= header1.npartTotal[0];
	}

      for(k=0, ntot_withmasses=0; k<5; k++)
	{
	  if(header1.mass[k]==0)
	    ntot_withmasses+= header1.npart[k];
	}

      if(i==0)
	allocate_memory();

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;
    

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&Id[pc_new], sizeof(int), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;


      if(ntot_withmasses>0)
	SKIP;
      for(k=0, pc_new=pc; k<6; k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      P[pc_new].Type=k;

	      if(header1.mass[k]==0)
		fread(&P[pc_new].Mass, sizeof(float), 1, fd);
	      else
		P[pc_new].Mass= header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses>0)
	SKIP;
      

      if(header1.npart[0]>0)
	{
	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].U, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		  fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		P[pc_sph].Ne= 1.0;
		pc_sph++;
	      }
	}

      fclose(fd);
    }


  Time= header1.time;
  Redshift= header1.time;
}

/* this routine allocates the memory for the 
 * particle data.
 *
 * Author: Volker Springel, Max Planck Institute
 */
int allocate_memory(void)
{
  printf("allocating memory...\n");

  if(!(P=malloc(NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  P--;   /* start with offset 1 */

  
  if(!(Id=malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  Id--;   /* start with offset 1 */

  printf("allocating memory...done\n");
}

/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 *
 * Author: Volker Springel, Max Planck Institute
 */
int reordering(void)
{
  int i,j;
  int idsource, idsave, dest;
  struct particle_data psave, psource;


  printf("reordering....\n");

  for(i=1; i<=NumPart; i++)
    {
      if(Id[i] != i)
	{
	  psource= P[i];
	  idsource=Id[i];
	  dest=Id[i];

	  do
	    {
	      psave= P[dest];
	      idsave=Id[dest];

	      P[dest]= psource;
	      Id[dest]= idsource;
	      
	      if(dest == i) 
		break;

	      psource= psave;
	      idsource=idsave;

	      dest=idsource;
	    }
	  while(1);
	}
    }

  printf("done.\n");

  Id++;   
  free(Id);

  printf("space for particle ID freed\n");
}

// This function frees the memory used when the load_snapshot function is 
// called. In general, users will not know that P has to be incremented
// because it was decremented and an attempt to free it will result in
// a segmentation fault. This prevents the problem.
int unload_snapshot() {

    P++; free(P);

    return 0;
}

// Load a Gadget2 snapshot into a galaxy object.
int load_gadget2_galaxy(char *filename,galaxy *galaxy_1) {

    int i,j,parts[2];

    load_snapshot(filename,1);
    parts[0] = header1.npart[2];
    parts[1] = header1.npart[1];
 
    // Create the galaxy, but pass false information for the irrelevant parts
    // at the moment.
    create_galaxy(galaxy_1,parts,0.025,0.025,0.05,15.0,1.6E7,64,1.0,0);

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
    unload_snapshot();

    return 0;
}

// This function prints the release date of the current version of Starscream.
void write_starscream_version() {

    printf("Starscream version 20090607.\n");

    return;
}
