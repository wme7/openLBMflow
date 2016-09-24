// openLBMflow v1.0.1 Copyright (C) 2013 LBMflow
// Open Source Lattice Boltzmann Solver
// www.lbmflow.com
// open@lbmflow.com

//Default configuration is 3D Schan-Chen Multiphase model
#define Lattice3D  // delete this line for 2D code
#define MultiPhase // delete this line for Singlephase code

// initial paramaters
int nx = 60;                //lattice size x
int ny = 100;               //lattice size y
int nz = 60;                //lattice size z (only for 3D code)
float tau = 1.0;           //relaxation time
float rhoh = 2.49;       //high density fluid (rhoh=1 for singlephase)
float rhol = 0.0734;       //low  density fluid (rhoh=1 for singlephase)
float ifaceW = 3.0;        //interface width
float G = -6.0;            //interparticular interaction potential (G=0 for singlephase)
float body_force = 1.0e-4; //gravity force
float body_force_dir = 0;  //gravity direction (0=down 90=right 180=top 270=left)
int time_total = 5500;      //total time step
int time_save = 50;         //save result to VTK image file (*.vti can be open in Paraview)
int save_rho = 1;           //save density  in output file (0=don't save, 1=save)
int save_pre = 1;           //save pressure in output file (0=don't save, 1=save)
int save_vel = 1;           //save velocity in output file (0=don't save, 1=save)
int boundary_bot = 1;       //0=periodic, 1=HBB, set half way bounce back on bottom wall
int boundary_top = 1;       //0=periodic, 1=HBB, set half way bounce back on top wall
int boundary_lef = 0;       //0=periodic, 1=HBB, set half way bounce back on left wall
int boundary_rig = 0;       //0=periodic, 1=HBB, set half way bounce back on right wall
int boundary_fro = 0;       //0=periodic, 1=HBB, set half way bounce back on front wall (only for 3D)
int boundary_bac = 0;       //0=periodic, 1=HBB, set half way bounce back on back wall (only for 3D)
float rho_boundary = 0.2;  //define density of solid walls
float top_wall_speed = 0.0;//speed of top wall for lid driven cavity
float bot_wall_speed = 0.0;//speed of bottom wall
int d1r = 15;               //droplet1 radius
int d1x = 30;               //droplet1 position x (for multiphase model)
int d1y = 20;               //droplet1 position y (for multiphase model)
int d1z = 30;               //droplet1 position z (for multiphase model)
int drop1 = -1.0;            //1=droplet, -1=buble
int d2r = 0;                //droplet2 radius
int d2x = 0;                //droplet2 position x (for multiphase model)
int d2y = 0;                //droplet2 position y (for multiphase model)
int d2z = 0;                //droplet2 position z (for multiphase model)
int drop2 = 1.0;            //1=droplet, -1=buble
///To add more droplet use function initialize_droplet(dx, dy, dz, dr, drop); in main function.
