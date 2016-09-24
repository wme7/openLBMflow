// openLBMflow v1.0.1 Copyright (C) 2013 LBMflow
// Open Source Lattice Boltzmann Solver
// www.lbmflow.com
// open@lbmflow.com

// LICENSE
// The openLBMflow code is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// DISCLAIMER OF WARRANTY
// The code is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// This software may contain errors that could cause failures or loss of data,
// and may be incomplete or contain inaccuracies.  You expressly acknowledge and agree
// that use of the openLBMflow software is at your sole risk.
// The openLBMflow software is provided 'AS IS' and without warranty of any kind.

// To comile use: gcc -O2 -o openLBMflow openLBMflow.c -lm

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <string.h>
#include <dirent.h>


//include file with initial parameters
#include "openLBMFlow_conf.c"


//define lattice constants D3Q19 1/3 1/18 1/36, D2Q9 4/9 1/9 1/36
#ifdef Lattice3D
double D=3, Q=19, w1 = 18.0, w2 = 36;
double wi[19]= {1.0/3.0, 1.0/18.0, 1.0/36.0, 1.0/18.0, 1.0/36.0, 1.0/18.0, 1.0/36.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
#else
int dz;
double D=2, Q=9, w1 = 9.0, w2 = 36;
double wi[9]= {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0};
#endif
double ex[19]= {0.0, 1.0, 1.0, 0.0,-1.0,-1.0,-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0,-1.0, 1.0};
double ey[19]= {0.0, 0.0, 1.0, 1.0, 1.0, 0.0,-1.0,-1.0,-1.0, 1.0, 0.0,-1.0,-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
double ez[19]= {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0};

//define global variables
int ***solid;
double ***rho, ***ux, ***uy, ***uz, ***p, ***phi;
double ****fp, ****fn, ****tmp_f;
FILE *errorfile;
int x, y, z, k, t, bot, top, lef, rig, fro, bac, error;
double feq[19], colision_operator[19], tmp_fn[19];
double tmp, tmp_rho, tmp_phi, tmp_ux, tmp_uy, tmp_uz, vel1, vel2, vel3, uxyz2, tmp_fp, ux2, uy2, uz2, uxy, uxz, uyz, uxy2, uxz2, uyz2;
double grad_phi_x, grad_phi_y, grad_phi_z, body_force_x, body_force_y, body_force_z;
double Speed, step_now, time_now, T, rh, press, mass;

//Function prototypes
int*** init_mem_int3d(int nx, int ny, int nz);
double*** init_mem_float3d(int nx, int ny, int nz);
double**** init_mem_float4d(int nx, int ny, int nz, int nu);

//define function prototypes
void initialize_memory();
void initialize_boundary(int boundary_bot, int boundary_top, int boundary_lef, int boundary_rig, int boundary_fro, int boundary_bac, double rho_boundary);
void initialize_density(double density);
void initialize_distrFunc();
void finalise();
void initialize_droplet(int dx, int dy, int dz, int dr, int drop);
void update();
void streaming(int x, int y, int z);
void boundary(int x, int y, int z);
void streaming_boundary(int x, int y, int z);
void equlibrium();
void outputSave();
void massConservation();
void writeVTK(int t, int nx, int ny, int nz, double ***rho, int write_rho, double ***pre, int write_pre, double ***ux, double ***uy, double ***uz, int write_vel, char *directory, char *filename);
void write_collection_pvd(int t, int nx, int ny, int nz, double ***array, char *directory, char *filename);



int main(int argc, char **argv)
{
    printf("openLBMflow v1.0.0 (c) 2010 www.lbmflow.com\n");
    initialize_memory();
    initialize_boundary(boundary_bot, boundary_top, boundary_lef, boundary_rig, boundary_fro, boundary_bac, rho_boundary);
    initialize_density( rhol );
    if (d1r>0) initialize_droplet(d1x, d1y, d1z, d1r, drop1);  //droplet 1
    if (d2r>0) initialize_droplet(d2x, d2y, d2z, d2r, drop2);  //droplet 2
    //droplet 3 initialize_droplet(dx, dy, dz, dr, drop);
    initialize_distrFunc();

    ///main iteration loop
    for (t=0; t<=time_total; t++)
    {
        if (t%time_save==0)
        {
            outputSave(); //save output to VTK image file
        }
        update(); //calculate distribution function for next time step
        if (error != 0)
        {
            // exit with error
            return 1;
        }
    }
    finalise(); //free allocated memory
    return 0; //exit program
}




void initialize_memory()
{
    mass = 0;
    error = 0;
#ifndef Lattice3D
    nz=1;
    dz = 1;
#endif
    /// Allocate RAM memory
    solid = init_mem_int3d(nx,ny,nz);
    fp    = init_mem_float4d(nx,ny,nz,Q);
    fn    = init_mem_float4d(nx,ny,nz,Q);
    rho   = init_mem_float3d(nx,ny,nz);
    ux    = init_mem_float3d(nx,ny,nz);
    uy    = init_mem_float3d(nx,ny,nz);
#ifdef Lattice3D
    uz    = init_mem_float3d(nx,ny,nz);
#endif

#ifdef MultiPhase
    phi = init_mem_float3d(nx,ny,nz);
#else
    rhol = 1.0;
    rhoh = 1.0;
    G = 0.0;
#endif

    write_collection_pvd(t, nx, ny, nz, &rho[0], "output", "openLBMflow");
}

void initialize_boundary(int boundary_bot, int boundary_top, int boundary_lef, int boundary_rig, int boundary_fro, int boundary_bac, double rho_boundary)
{
    ///Initialize type of nodes
    for (x=0; x<nx; x++)
    {
        for (y=0; y<ny; y++)
        {
            for (z=0; z<nz; z++)
            {
                solid[x][y][z] = 0;
            }
        }
    }

    ///Define Bounce Back Boundary
    for (x=0; x<nx; x++)
    {
        for (z=0; z<nz; z++)
        {
            if (boundary_bot>0)
            {
                y=0;      // node on bottom boundary
                solid[x][y][z] = 1;
                rho[x][y][z] = rho_boundary*(rhoh-rhol)+rhol;
                if (bot_wall_speed!=0)
                {
                    ux[x][y][z]=bot_wall_speed;
                    uy[x][y][z]=0;
                    uz[x][y][z]=0;
                }
            }
            if (boundary_top>0)
            {
                y=ny-1;    // node on top boundary
                solid[x][y][z] = 2;
                rho[x][y][z] = rho_boundary*(rhoh-rhol)+rhol;
                if (top_wall_speed!=0)
                {
                    ux[x][y][z]=top_wall_speed;
                    uy[x][y][z]=0;
                    uz[x][y][z]=0;
                }
            }
        }
    }

    for (y=0; y<ny; y++)
    {
        for (z=0; z<nz; z++)
        {
            if (boundary_lef>0)
            {
                x=0;      // node on left boundary
                solid[x][y][z] = 3;
                rho[x][y][z] = rho_boundary*(rhoh-rhol)+rhol;
            }
            if (boundary_rig>0)
            {
                x=nx-1;    // node on right boundary
                solid[x][y][z] = 4;
                rho[x][y][z] = rho_boundary*(rhoh-rhol)+rhol;
            }
        }
    }

    for (x=0; x<nx; x++)
    {
        for (y=0; y<ny; y++)
        {
            if (boundary_fro>0)
            {
                z=0;       // node on front boundary
                solid[x][y][z] = 5;
                rho[x][y][z] = rho_boundary*(rhoh-rhol)+rhol;
            }
            if (boundary_bac>0)
            {
                z=nz-1;    // node on back boundary
                solid[x][y][z] = 6;
                rho[x][y][z] = rho_boundary*(rhoh-rhol)+rhol;
            }
        }
    }
}

void initialize_density(double density)
{
    for (x=0; x<nx; x++)
    {
        for (y=0; y<ny; y++)
        {
            for (z=0; z<nz; z++)
            {
                if (solid[x][y][z] == 0)
                {
                    rho[x][y][z] = density;
                }
            }
        }
    }
}


void initialize_distrFunc()
{
    /// body_force vector
    body_force_x =  body_force*sin(body_force_dir/(180.0/M_PI));
    body_force_y = -body_force*cos(body_force_dir/(180.0/M_PI));
    body_force_z = 0;

    tmp_ux = 0.0;
    tmp_uy = 0.0;
    tmp_uz = 0.0;

    for (x=0; x<nx; x++)
    {
        for (y=0; y<ny; y++)
        {
            for (z=0; z<nz; z++)
            {
                if (solid[x][y][z] == 0)
                {
                    ux[x][y][z] = tmp_ux;
                    uy[x][y][z] = tmp_uy;
#ifdef Lattice3D
                    uz[x][y][z] = tmp_uz;
#endif
                    for (k=0; k<Q; k++)
                    {
#ifdef Lattice3D
                        vel1 = ex[k]*ux[x][y][z] + ey[k]*uy[x][y][z] + ez[k]*uz[x][y][z];
                        vel3 = ux[x][y][z]*ux[x][y][z] + uy[x][y][z]*uy[x][y][z] +  uz[x][y][z]*uz[x][y][z];
#else
                        vel1 = ex[k]*ux[x][y][z] + ey[k]*uy[x][y][z];
                        vel3 = ux[x][y][z]*ux[x][y][z] + uy[x][y][z]*uy[x][y][z];
#endif
                        vel2 = vel1 * vel1;
                        feq[k] = wi[k]*rho[x][y][z]*(1.0 + 3.0*vel1 + 4.5*vel2 - 1.5*vel3);
                        fp[x][y][z][k] = feq[k];
                        fn[x][y][z][k] = feq[k];
                    } //end of k loop
                } //end of solid
            }//end of z loop
        } //end of y loop
    } //end of x loop
}


void initialize_droplet(int dx, int dy, int dz, int dr, int drop)
{
    for (x=0; x<nx; x++)
    {
        for (y=0; y<ny; y++)
        {
            for (z=0; z<nz; z++)
            {
                if (solid[x][y][z] == 0)
                {
                    double dx2, dy2, dz2, radius;
                    dx2 = ((double)(x - dx))*((double)(x - dx));
                    dy2 = ((double)(y - dy))*((double)(y - dy));
#ifdef Lattice3D
                    dz2 = ((double)(z - dz))*((double)(z - dz));
#else
                    dz2 = 0.0;
#endif
                    radius  = sqrt(dx2 + dy2 + dz2);
                    tmp = 0.5 * ( (rhoh + rhol) - drop * (rhoh - rhol) * tanh ((radius - dr) / ifaceW * 2.0) );
                    if (tmp>rho[x][y][z]) rho[x][y][z] = tmp;
                }
            }
        }
    }
}

void update()
{
    //Calculate rho, u,v
    for (x=0; x<nx; x++)
    {
        for (y=0; y<ny; y++)
        {
            for (z=0; z<nz; z++)
            {
                if (solid[x][y][z]==0)
                {
                    //calculate rho and ux, uy
#ifdef Lattice3D
                    tmp_rho =  fp[x][y][z][0]+fp[x][y][z][1]+fp[x][y][z][2]+fp[x][y][z][3]+fp[x][y][z][4]+fp[x][y][z][5]+fp[x][y][z][6]+fp[x][y][z][7]+fp[x][y][z][8]+fp[x][y][z][9]+fp[x][y][z][10]+fp[x][y][z][11]+fp[x][y][z][12]+fp[x][y][z][13]+fp[x][y][z][14]+fp[x][y][z][15]+fp[x][y][z][16]+fp[x][y][z][17]+fp[x][y][z][18];
                    tmp_ux  = (fp[x][y][z][1]+fp[x][y][z][2]+fp[x][y][z][8]-fp[x][y][z][4]-fp[x][y][z][5]-fp[x][y][z][6]+fp[x][y][z][15]+fp[x][y][z][18]-fp[x][y][z][16]-fp[x][y][z][17])/tmp_rho;
                    tmp_uy  = (fp[x][y][z][2]+fp[x][y][z][3]+fp[x][y][z][4]-fp[x][y][z][6]-fp[x][y][z][7]-fp[x][y][z][8]+fp[x][y][z][9]+fp[x][y][z][14]-fp[x][y][z][11]-fp[x][y][z][12])/tmp_rho;
                    tmp_uz  = (fp[x][y][z][9]+fp[x][y][z][10]+fp[x][y][z][11]-fp[x][y][z][12]-fp[x][y][z][13]-fp[x][y][z][14]+fp[x][y][z][15]+fp[x][y][z][16]-fp[x][y][z][17]-fp[x][y][z][18])/tmp_rho;
#else
                    tmp_rho =  fp[x][y][z][0]+fp[x][y][z][1]+fp[x][y][z][2]+fp[x][y][z][3]+fp[x][y][z][4]+fp[x][y][z][5]+fp[x][y][z][6]+fp[x][y][z][7]+fp[x][y][z][8];
                    tmp_ux  = (fp[x][y][z][1]+fp[x][y][z][2]+fp[x][y][z][8]-fp[x][y][z][4]-fp[x][y][z][5]-fp[x][y][z][6])/tmp_rho;
                    tmp_uy  = (fp[x][y][z][2]+fp[x][y][z][3]+fp[x][y][z][4]-fp[x][y][z][6]-fp[x][y][z][7]-fp[x][y][z][8])/tmp_rho;
#endif
                    if ( (tmp_rho!=tmp_rho)&&(error == 0) )
                    {
                        error = 1;
                        printf("Exit with error in node x=%d y=%d z=%d at Timestep=%d\n", x, y, z, t);
                        errorfile = fopen("error.log","w");
                        fprintf(errorfile, "Exit with error in node x=%d y=%d z=%d at Timestep=%d\n", x, y, z, t);
                        fclose(errorfile);
                    }
                    //body force
                    tmp_ux += tau*body_force_x;
                    tmp_uy += tau*body_force_y;
#ifdef Lattice3D
                    tmp_uz += tau*body_force_z;
#endif
                    rho[x][y][z] = tmp_rho;
#ifdef MultiPhase
                    ux[x][y][z] = tmp_ux;
                    uy[x][y][z] = tmp_uy;
#ifdef Lattice3D
                    uz[x][y][z] = tmp_uz;
#endif
                }
                //calculate interparticular force in multiphase Schan-Chen model
                phi[x][y][z] = 1.0-exp(-rho[x][y][z]);
            }
        }
    }

    for (x=0; x<nx; x++)
    {
        for (y=0; y<ny; y++)
        {
            for (z=0; z<nz; z++)
            {
                if (solid[x][y][z]==0)
                {
                    bot = (y+ny-1)%ny;
                    top = (y+1)%ny;
                    lef = (x+nx-1)%nx;
                    rig = (x+1)%nx;
#ifdef Lattice3D
                    fro = (z+nz-1)%nz;
                    bac = (z+1)%nz;
#endif
                    tmp_rho = rho[x][y][z];
                    tmp_phi = phi[x][y][z];
                    tmp_ux  = ux[x][y][z];
                    tmp_uy  = uy[x][y][z];
#ifdef Lattice3D
                    tmp_uz  = uz[x][y][z];
#endif
                    // calculate gradient psi-gradient
                    grad_phi_x = (phi[rig][y][z]-phi[lef][y][z])/w1;
                    grad_phi_y = (phi[x][top][z]-phi[x][bot][z])/w1;
                    grad_phi_x+= (phi[rig][top][z]-phi[lef][top][z]+phi[rig][bot][z]-phi[lef][bot][z])/w2;
                    grad_phi_y+= (phi[rig][top][z]+phi[lef][top][z]-phi[lef][bot][z]-phi[rig][bot][z])/w2;
#ifdef Lattice3D
                    // 3d part
                    grad_phi_z = (phi[x][y][bac]-phi[x][y][fro])/w1;
                    grad_phi_z+= (phi[rig][y][bac]+phi[lef][y][bac]-phi[lef][y][fro]-phi[rig][y][fro])/w2;
                    grad_phi_x+= (phi[rig][y][bac]-phi[lef][y][bac]+phi[rig][y][fro]-phi[lef][y][fro])/w2;
                    grad_phi_y+= (phi[x][top][bac]+phi[x][top][fro]-phi[x][bot][bac]-phi[x][bot][fro])/w2;
                    grad_phi_z+= (phi[x][top][bac]+phi[x][bot][bac]-phi[x][bot][fro]-phi[x][top][fro])/w2;
#endif
                    // interparticule potential in equilibrium velocity
                    tmp_ux += tau*(-G*tmp_phi*grad_phi_x)/tmp_rho;
                    tmp_uy += tau*(-G*tmp_phi*grad_phi_y)/tmp_rho;
#ifdef Lattice3D
                    tmp_uz += tau*(-G*tmp_phi*grad_phi_z)/tmp_rho;
#endif

#else //SinglePhase
                    bot = (y+ny-1)%ny;
                    top = (y+1)%ny;
                    lef = (x+nx-1)%nx;
                    rig = (x+1)%nx;
#ifdef Lattice3D
                    fro = (z+nz-1)%nz;
                    bac = (z+1)%nz;
#endif

#endif //end Single/Multi Phase

                    ux[x][y][z] = tmp_ux;
                    uy[x][y][z] = tmp_uy;
#ifdef Lattice3D
                    uz[x][y][z] = tmp_uz;
                    uxyz2 = (tmp_ux)*(tmp_ux) + (tmp_uy)*(tmp_uy) + (tmp_uz)*(tmp_uz);
#else
                    uxyz2 = (tmp_ux)*(tmp_ux) + (tmp_uy)*(tmp_uy);
#endif
                    ux2 = tmp_ux*tmp_ux;
                    uy2 = tmp_uy*tmp_uy;
                    uz2 = tmp_uz*tmp_uz;
                    uxy2 = ux2+uy2;
                    uxz2 = ux2+uz2;
                    uyz2 = uy2+uz2;
                    uxy = 2.0*tmp_ux*tmp_uy;
                    uxz = 2.0*tmp_ux*tmp_uz;
                    uyz = 2.0*tmp_uy*tmp_uz;

                    tmp_fn[0] = fp[x][y][z][0] - (fp[x][y][z][0] - (wi[0]*tmp_rho*(1.0 - 1.5*uxyz2)))/tau;
                    tmp_fn[1] = fp[x][y][z][1] - (fp[x][y][z][1] - (wi[1]*tmp_rho*(1.0 + 3.0*tmp_ux             + 4.5*ux2        - 1.5*uxyz2)))/tau;
                    tmp_fn[2] = fp[x][y][z][2] - (fp[x][y][z][2] - (wi[2]*tmp_rho*(1.0 + 3.0*(+tmp_ux+tmp_uy)   + 4.5*(uxy2+uxy) - 1.5*uxyz2)))/tau;
                    tmp_fn[3] = fp[x][y][z][3] - (fp[x][y][z][3] - (wi[3]*tmp_rho*(1.0 + 3.0*tmp_uy             + 4.5*uy2        - 1.5*uxyz2)))/tau;
                    tmp_fn[4] = fp[x][y][z][4] - (fp[x][y][z][4] - (wi[4]*tmp_rho*(1.0 + 3.0*(-tmp_ux+tmp_uy)   + 4.5*(uxy2-uxy) - 1.5*uxyz2)))/tau;
                    tmp_fn[5] = fp[x][y][z][5] - (fp[x][y][z][5] - (wi[5]*tmp_rho*(1.0 - 3.0*tmp_ux             + 4.5*ux2        - 1.5*uxyz2)))/tau;
                    tmp_fn[6] = fp[x][y][z][6] - (fp[x][y][z][6] - (wi[6]*tmp_rho*(1.0 + 3.0*(-tmp_ux-tmp_uy)   + 4.5*(uxy2+uxy) - 1.5*uxyz2)))/tau;
                    tmp_fn[7] = fp[x][y][z][7] - (fp[x][y][z][7] - (wi[7]*tmp_rho*(1.0 - 3.0*tmp_uy             + 4.5*uy2        - 1.5*uxyz2)))/tau;
                    tmp_fn[8] = fp[x][y][z][8] - (fp[x][y][z][8] - (wi[8]*tmp_rho*(1.0 + 3.0*(+tmp_ux-tmp_uy)   + 4.5*(uxy2-uxy) - 1.5*uxyz2)))/tau;
#ifdef Lattice3D
                    tmp_fn[9]  = fp[x][y][z][9]  - (fp[x][y][z][9]  - (wi[9] *tmp_rho*(1.0 + 3.0*(+tmp_uy+tmp_uz)   + 4.5*(uyz2+uyz) - 1.5*uxyz2)))/tau;
                    tmp_fn[10] = fp[x][y][z][10] - (fp[x][y][z][10] - (wi[10]*tmp_rho*(1.0 + 3.0*(+tmp_uz       )   + 4.5*(uz2     ) - 1.5*uxyz2)))/tau;
                    tmp_fn[11] = fp[x][y][z][11] - (fp[x][y][z][11] - (wi[11]*tmp_rho*(1.0 + 3.0*(-tmp_uy+tmp_uz)   + 4.5*(uyz2-uyz) - 1.5*uxyz2)))/tau;
                    tmp_fn[12] = fp[x][y][z][12] - (fp[x][y][z][12] - (wi[12]*tmp_rho*(1.0 + 3.0*(-tmp_uy-tmp_uz)   + 4.5*(uyz2+uyz) - 1.5*uxyz2)))/tau;
                    tmp_fn[13] = fp[x][y][z][13] - (fp[x][y][z][13] - (wi[13]*tmp_rho*(1.0 + 3.0*(-tmp_uz       )   + 4.5*(uz2     ) - 1.5*uxyz2)))/tau;
                    tmp_fn[14] = fp[x][y][z][14] - (fp[x][y][z][14] - (wi[14]*tmp_rho*(1.0 + 3.0*(+tmp_uy-tmp_uz)   + 4.5*(uyz2-uyz) - 1.5*uxyz2)))/tau;
                    tmp_fn[15] = fp[x][y][z][15] - (fp[x][y][z][15] - (wi[15]*tmp_rho*(1.0 + 3.0*(+tmp_ux+tmp_uz)   + 4.5*(uxz2+uxz) - 1.5*uxyz2)))/tau;
                    tmp_fn[16] = fp[x][y][z][16] - (fp[x][y][z][16] - (wi[16]*tmp_rho*(1.0 + 3.0*(-tmp_ux+tmp_uz)   + 4.5*(uxz2-uxz) - 1.5*uxyz2)))/tau;
                    tmp_fn[17] = fp[x][y][z][17] - (fp[x][y][z][17] - (wi[17]*tmp_rho*(1.0 + 3.0*(-tmp_ux-tmp_uz)   + 4.5*(uxz2+uxz) - 1.5*uxyz2)))/tau;
                    tmp_fn[18] = fp[x][y][z][18] - (fp[x][y][z][18] - (wi[18]*tmp_rho*(1.0 + 3.0*(+tmp_ux-tmp_uz)   + 4.5*(uxz2-uxz) - 1.5*uxyz2)))/tau;
#endif

                    streaming(x, y, z);
                }
            }
        }
    }

    for (x=0; x<nx; x++)
    {
        for (y=0; y<ny; y++)
        {
            for (z=0; z<nz; z++)
            {
                if (solid[x][y][z]!=0)
                {
                    boundary(x, y, z);
                    streaming_boundary(x, y, z);
                }
            }
        }
    }
    //swiching the pointers of distribution functions
    tmp_f = fp;
    fp  = fn;
    fn  = tmp_f;
}

void streaming(int x, int y, int z)
{
    top = (y+1)%ny;
    bot = (y+ny-1)%ny;
    rig = (x+1)%nx;
    lef = (x+nx-1)%nx;

    fn[x][y][z][0]     = tmp_fn[0];
    fn[rig][y][z][1]   = tmp_fn[1];
    fn[rig][top][z][2] = tmp_fn[2];
    fn[x][top][z][3]   = tmp_fn[3];
    fn[lef][top][z][4] = tmp_fn[4];
    fn[lef][y][z][5]   = tmp_fn[5];
    fn[lef][bot][z][6] = tmp_fn[6];
    fn[x][bot][z][7]   = tmp_fn[7];
    fn[rig][bot][z][8] = tmp_fn[8];

#ifdef Lattice3D
    bac = (z+1)%nz;
    fro = (z+nz-1)%nz;

    fn[x][top][bac][9]  = tmp_fn[9];
    fn[x][y][bac][10]   = tmp_fn[10];
    fn[x][bot][bac][11] = tmp_fn[11];
    fn[x][bot][fro][12] = tmp_fn[12];
    fn[x][y][fro][13]   = tmp_fn[13];
    fn[x][top][fro][14] = tmp_fn[14];
    fn[rig][y][bac][15] = tmp_fn[15];
    fn[lef][y][bac][16] = tmp_fn[16];
    fn[lef][y][fro][17] = tmp_fn[17];
    fn[rig][y][fro][18] = tmp_fn[18];
#endif
}

void boundary(int x, int y, int z)
{
    // Half-Way bounce back
    tmp = fn[x][y][z][1];
    fn[x][y][z][1] = fn[x][y][z][5];
    fn[x][y][z][5] = tmp;

    tmp = fn[x][y][z][2];
    fn[x][y][z][2] = fn[x][y][z][6];
    fn[x][y][z][6] = tmp;

    tmp = fn[x][y][z][3];
    fn[x][y][z][3] = fn[x][y][z][7];
    fn[x][y][z][7] = tmp;

    tmp = fn[x][y][z][4];
    fn[x][y][z][4] = fn[x][y][z][8];
    fn[x][y][z][8] = tmp;

    if (top_wall_speed != 0) //mouving top wall
    {
        fn[x][y][z][6] -= top_wall_speed*rho[x][y][z]/6;
        fn[x][y][z][8] += top_wall_speed/rho[x][y][z]/6;
    }
    if (bot_wall_speed != 0) //mouving bottom wall
    {
        fn[x][y][z][2] += bot_wall_speed*rho[x][y][z]/6;
        fn[x][y][z][4] -= bot_wall_speed/rho[x][y][z]/6;
    }

#ifdef Lattice3D
    tmp = fn[x][y][z][9];
    fn[x][y][z][9] = fn[x][y][z][12];
    fn[x][y][z][12] = tmp;

    tmp = fn[x][y][z][10];
    fn[x][y][z][10] = fn[x][y][z][13];
    fn[x][y][z][13] = tmp;

    tmp = fn[x][y][z][14];
    fn[x][y][z][14] = fn[x][y][z][11];
    fn[x][y][z][11] = tmp;

    tmp = fn[x][y][z][15];
    fn[x][y][z][15] = fn[x][y][z][17];
    fn[x][y][z][17] = tmp;

    tmp = fn[x][y][z][18];
    fn[x][y][z][18] = fn[x][y][z][16];
    fn[x][y][z][16] = tmp;
#endif
}

void streaming_boundary(int x, int y, int z)
{
    top = (y+1)%ny;
    bot = (y+ny-1)%ny;
    rig = (x+1)%nx;
    lef = (x+nx-1)%nx;

    fn[x][y][z][0]     = fn[x][y][z][0];
    fn[rig][y][z][1]   = fn[x][y][z][1];
    fn[rig][top][z][2] = fn[x][y][z][2];
    fn[x][top][z][3]   = fn[x][y][z][3];
    fn[lef][top][z][4] = fn[x][y][z][4];
    fn[lef][y][z][5]   = fn[x][y][z][5];
    fn[lef][bot][z][6] = fn[x][y][z][6];
    fn[x][bot][z][7]   = fn[x][y][z][7];
    fn[rig][bot][z][8] = fn[x][y][z][8];

#ifdef Lattice3D
    bac = (z+1)%nz;
    fro = (z+nz-1)%nz;

    fn[x][top][bac][9]  = fn[x][y][z][9];
    fn[x][y][bac][10]   = fn[x][y][z][10];
    fn[x][bot][bac][11] = fn[x][y][z][11];
    fn[x][bot][fro][12] = fn[x][y][z][12];
    fn[x][y][fro][13]   = fn[x][y][z][13];
    fn[x][top][fro][14] = fn[x][y][z][14];
    fn[rig][y][bac][15] = fn[x][y][z][15];
    fn[lef][y][bac][16] = fn[x][y][z][16];
    fn[lef][y][fro][17] = fn[x][y][z][17];
    fn[rig][y][fro][18] = fn[x][y][z][18];
#endif
}


int*** init_mem_int3d(int nx, int ny, int nz)
{
    int ***p; /*declaration of p as: pointer-to-pointer-to-pointer of int */
    int x, y;

    p = malloc( nx * sizeof(*p) ); /*Allocate pointers for the nx */
    if (p != NULL)
    {
        for (x = 0; x < nx; x++)
        {
            p[x] = malloc( ny * sizeof **p );/*Allocate pointers for the ny */
            if (p[x] == NULL)
            {
                printf("Memory allocation failed. Exiting....");
                //return (1);
            }
            else
            {
                for (y = 0; y < ny; y++)
                {
                    p[x][y] = malloc( nz * sizeof ***p ); /*Allocate pointers for the nz */
                    if (p[x][y] == NULL)
                    {
                        printf("Memory allocation failed. Exiting....");
                        //return (1);
                    }
                }
            }
        }
    }
    else
    {
        printf("Memory allocation failed. Exiting....");
        //return (1);
    }
    return p;
}


double*** init_mem_float3d(int nx, int ny, int nz)
{
    double ***p; /*declaration of p as: pointer-to-pointer-to-pointer of int */
    int x, y;

    p = malloc( nx * sizeof(*p) ); /*Allocate pointers for the nx */
    if (p != NULL)
    {
        for (x = 0; x < nx; x++)
        {
            p[x] = malloc( ny * sizeof **p );/*Allocate pointers for the ny */
            if (p[x] == NULL)
            {
                printf("Memory allocation failed. Exiting....");
                //return (1);
            }
            else
            {
                for (y = 0; y < ny; y++)
                {
                    p[x][y] = malloc( nz * sizeof ***p ); /*Allocate pointers for the nz */
                    if (p[x][y] == NULL)
                    {
                        printf("Memory allocation failed. Exiting....");
                        //return (1);
                    }
                }
            }
        }
    }
    else
    {
        printf("Memory allocation failed. Exiting....");
        //return (1);
    }
    return p;
}


double**** init_mem_float4d(int nx, int ny, int nz, int nu)
{
    double ****p; /*declaration of p as: pointer-to-pointer-to-pointer of int */
    int x, y, z;

    p = malloc( nx * sizeof(*p) ); /*Allocate pointers for the nx */
    if (p != NULL)
    {
        for (x = 0; x < nx; x++)
        {
            p[x] = malloc( ny * sizeof **p );/*Allocate pointers for the ny */
            if (p[x] == NULL)
            {
                printf("Memory allocation failed. Exiting....");
                //return (1);
            }
            else
            {
                for (y = 0; y < ny; y++)
                {
                    p[x][y] = malloc( nz * sizeof ***p ); /*Allocate pointers for the nz */
                    if (p[x][y] == NULL)
                    {
                        printf("Memory allocation failed. Exiting....");
                        //return (1);
                    }
                    else
                    {
                        for (z = 0; z < nz; z++)
                        {
                            p[x][y][z] = malloc( nu * sizeof ****p ); /*Allocate pointers for the nz */
                            if (p[x][y][z] == NULL)
                            {
                                printf("Memory allocation failed. Exiting....");
                                //return (1);
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        printf("Memory allocation failed. Exiting....");
        //return (1);
    }
    return p;

}

void outputSave()
{
    writeVTK(t, nx, ny, nz, &rho[0], save_rho, &rho[0], save_pre, &ux[0], &uy[0], &uz[0], save_vel, "output", "openLBMflow");

    //check mass conservation for debug only
    //massConservation();

    //Calculate Mega Lattice Site Update per second MLSU/s
    Speed = (nx*ny*nz)*(t-step_now)/((clock() - time_now)/CLOCKS_PER_SEC)/1000000.0;
    step_now = t;
    time_now = clock();
    if (mass == 0) printf("t=%d\tSpeed=%f MLUP/s\n", t, Speed);
    else printf("t=%d\tSpeed=%f MLUP/s mass=%f\n", t, Speed, mass);
}

void massConservation()
{
    mass = 0.0;
    for (x=0; x<nx; x++)
    {
        for (y=0; y<ny; y++)
        {
            for (z=0; z<nz; z++)
            {
                mass += rho[x][y][z];
            }
        }
    }
}

void writeVTK(int t, int nx, int ny, int nz, double ***rho, int write_rho, double ***pre, int write_pre, double ***ux, double ***uy, double ***uz, int write_vel, char *directory, char *filename)

{
    int x,y,z,dir;
    char dataFileName[255];
    FILE *dataFile;
#ifdef WIN32
    dir = mkdir(directory);
#else
    dir = mkdir(directory,0777);
#endif
    if (dir==0) printf("Error: Can't create output directory!\n");
    sprintf(dataFileName,"%s/%s_%07d.vti",directory,filename,t);
    dataFile = fopen(dataFileName,"w");
    fprintf(dataFile, "<?xml version=\"1.0\"?>\n");
    fprintf(dataFile, "<!-- openLBMflow v1.0.1, www.lbmflow.com -->\n");
    fprintf(dataFile, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(dataFile, "  <ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n",nx-1,ny-1,nz-1);
    fprintf(dataFile, "  <Piece Extent=\"0 %d 0 %d 0 %d\">\n",nx-1,ny-1,nz-1);
    fprintf(dataFile, "    <PointData Scalars=\"scalars\">\n");

    //write density
    if (write_rho != 0)
    {
        fprintf(dataFile, "      <DataArray type=\"Float32\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for (z=0; z<nz; z++)
        {
            for (y=0; y<ny; y++)
            {
                for (x=0; x<nx; x++)
                {
                    fprintf(dataFile,"%.4e ", rho[x][y][z]);
                }
                fprintf(dataFile, "\n");
            }
        }
        fprintf(dataFile, "      </DataArray>\n");
    }

    //write pressure
    if (write_pre != 0)
    {
        fprintf(dataFile, "      <DataArray type=\"Float32\" Name=\"Pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for (z=0; z<nz; z++)
        {
            for (y=0; y<ny; y++)
            {
                for (x=0; x<nx; x++)
                {
#ifdef MultiPhase
                    fprintf(dataFile,"%.4e ", rho[x][y][z]/3.0 + G*phi[x][y][z]*phi[x][y][z]/6.0);
#else
                    fprintf(dataFile,"%.4e ", rho[x][y][z]/3.0);
#endif
                }
                fprintf(dataFile, "\n");
            }
        }
        fprintf(dataFile, "      </DataArray>\n");
    }

    //fprintf(dataFile, "    <PointData Vectors=\"Velocity\">\n");
    //write velocity
    if (write_vel != 0)
    {
        fprintf(dataFile, "      <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for (z=0; z<nz; z++)
        {
            for (y=0; y<ny; y++)
            {
                for (x=0; x<nx; x++)
                {
                    fprintf(dataFile,"%.4e ", ux[x][y][z]);
                    fprintf(dataFile,"%.4e ", uy[x][y][z]);
#ifdef Lattice3D
                    fprintf(dataFile,"%.4e ", uz[x][y][z]);
#else
                    fprintf(dataFile,"%d ", 0);
#endif

                }
                fprintf(dataFile, "\n");
            }
        }
        fprintf(dataFile, "      </DataArray>\n");
    }
    //fprintf(dataFile, "    </PointData>\n");


    fprintf(dataFile, "    </PointData>\n");

    fprintf(dataFile, "    <CellData>\n");
    fprintf(dataFile, "    </CellData>\n");
    fprintf(dataFile, "  </Piece>\n");
    fprintf(dataFile, "  </ImageData>\n");

    fprintf(dataFile, "</VTKFile>\n");
    fclose(dataFile);
}

void write_collection_pvd(int t, int nx, int ny, int nz, double ***array, char *directory, char *filename)
{
    int x,dir;
    char dataFileName[255];
    FILE *dataFile;
#ifdef WIN32
    dir = mkdir(directory);
#else
    dir = mkdir(directory,0777);
#endif
    if (dir==0) printf("Error: Can't create output directory!\n");
    sprintf(dataFileName,"%s/%s.pvd",directory,filename);
    dataFile = fopen(dataFileName,"w");
    fprintf(dataFile, "<?xml version=\"1.0\"?>\n");
    fprintf(dataFile, "<!-- openLBMflow v1.0.1, www.lbmflow.com -->\n");
    fprintf(dataFile, "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(dataFile, "  <Collection>\n");
    for (x=0; x<=time_total; x += time_save)
    {
        fprintf(dataFile, "    <DataSet  timestep=\"%d\" group=\"\" part=\"%d\" file=\"%s_%07d.vti\"/>\n",x,0,filename,x);
    }
    fprintf(dataFile, "  </Collection>\n");
    fprintf(dataFile, "</VTKFile>\n");
    fclose(dataFile);
}

void finalise()
{
    //free allocated memory
    free(solid);
    free(fp);
    free(fn);
    free(rho);
    free(ux);
    free(uy);
#ifdef Lattice3D
    free(uz);
#endif

#ifdef MultiPhase
    free(phi);
#endif
}

