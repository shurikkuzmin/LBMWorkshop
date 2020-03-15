// This is an example 2D immersed boundary lattice Boltzmann method code.
// It uses the D2Q9 lattice with Guo's forcing term.
// Rigid bottom and top walls are parallel to the x-axis (channel).
// The flow is periodic along the x-axis.
// One initially cylindrical particle is positioned in the flow.
// This particle can be rigid/deformable and stationary/moving.
//
// Last update 22-Aug-2011 by Timm Krüger.
// This code may be changed and distributed freely.
//
// The lattice velocities are defined according to the following scheme:
// index:   0  1  2  3  4  5  6  7  8
// ----------------------------------
// x:       0 +1 -1  0  0 +1 -1 +1 -1
// y:       0  0  0 +1 -1 +1 -1 -1 +1
//
// 8 3 5  ^y
//  \|/   |   x
// 2-0-1   --->
//  /|\
// 6 4 7

/// *********************
/// PREPROCESSOR COMMANDS
/// *********************

#include <vector> // vector containers
#include <cmath> // mathematical library
#include <iostream> // for the use of 'cout'
#include <fstream> // file streams
#include <sstream> // string streams
#include <cstdlib> // standard library
#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code

using namespace std; // permanently use the standard namespace

/// *********************
/// SIMULATION PARAMETERS
/// *********************

// These are the relevant simulation parameters.
// They can be changed by the user.
// If a bottom or top wall shall move in negative x-direction, a negative velocity has to be specified.
// Moving walls and gravity can be switched on simultaneously.

/// Simulation types
// Exactly one of the following options has to be defined.
// RIGID_CYLINDER
// - for simulation of, e.g., Karman vortex street
// - the cylinder is kept in space, it does not move
// DEFORMABLE_CYLINDER
// - for simulation of a deformable cylinder
// - the cylinder moves along with the flow
// DEFORMABLE_RBC
// - for simulation of a deformable red blood cell
// - the cell moves along with the flow

#define RIGID_CYLINDER
//#define DEFORMABLE_CYLINDER
//#define DEFORMABLE_RBC

/// Fluid/lattice properties

#ifdef RIGID_CYLINDER
  const int Nx = 220; // number of lattice nodes along the x-axis (periodic)
  const int Ny = 62; // number of lattice nodes along the y-axis (including two wall nodes)
  const double tau = 0.525; // relaxation time
  const int t_num = 50000; // number of time steps (running from 1 to t_num)
  const int t_disk = 200; // disk write time step (data will be written to the disk every t_disk step)
  const int t_info = 1000; // info time step (screen message will be printed every t_info step)
  const double gravity = 0.000008; // force density due to gravity (in positive x-direction)
  const double wall_vel_bottom = 0; // velocity of the bottom wall (in positive x-direction)
  const double wall_vel_top = 0; // velocity of the top wall (in positive x-direction)
#else
  const int Nx = 30; // number of lattice nodes along the x-axis (periodic)
  const int Ny = 32; // number of lattice nodes along the y-axis (including two wall nodes)
  const double tau = 1; // relaxation time
  const int t_num = 50000; // number of time steps (running from 1 to t_num)
  const int t_disk = 200; // disk write time step (data will be written to the disk every t_disk step)
  const int t_info = 1000; // info time step (screen message will be printed every t_info step)
  const double gravity = 0.0001; // force density due to gravity (in positive x-direction)
  const double wall_vel_bottom = 0; // velocity of the bottom wall (in positive x-direction)
  const double wall_vel_top = 0; // velocity of the top wall (in positive x-direction)
#endif

/// Particle properties

#ifdef RIGID_CYLINDER
  const int particle_num_nodes = 36; // number of surface nodes
  const double particle_radius = 8; // radius
  const double particle_stiffness = 0.03; // stiffness modulus
  const double particle_center_x = 30; // center position (x-component)
  const double particle_center_y = 29; // center position (y-component)
#else
  const int particle_num_nodes = 32; // number of surface nodes
  const double particle_radius = 6; // radius
  const double particle_stiffness = 0.1; // stiffness modulus
  const double particle_bending = 0.01; // bending modulus
  const double particle_center_x = 15; // center position (x-component)
  const double particle_center_y = 8; // center position (y-component)
#endif

/// *****************
/// DECLARE VARIABLES
/// *****************

// The following code should not be modified when it is first used.

const double omega = 1. / tau; // relaxation frequency (inverse of relaxation time)
double ***pop, ***pop_old; // LBM populations (old and new)
double **density; // fluid density
double **velocity_x; // fluid velocity (x-component)
double **velocity_y; // fluid velocity (y-component)
double **force_x; // fluid force (x-component)
double **force_y; // fluid force (y-component)
double force_latt[9]; // lattice force term entering the lattice Boltzmann equation
double pop_eq[9]; // equilibrium populations
const double weight[9] = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.}; // lattice weights

/// ******************
/// PARTICLE STRUCTURE
/// ******************

// The following code handles the object immersed in the flow.
// In the present implementation, only a single object can be put into the flow.

/// Structure for surface nodes
// Each node has a current x- and y-position and a reference x- and y-position.

struct node_struct {

  /// Constructor

  node_struct() {
   x = 0;
   y = 0;
   x_ref = 0;
   y_ref = 0;
   vel_x = 0;
   vel_y = 0;
   force_x = 0;
   force_y = 0;
  }

  /// Elements

  double x; // current x-position
  double y; // current y-position
  double x_ref; // reference x-position
  double y_ref; // reference y-position
  double vel_x; // node velocity (x-component)
  double vel_y; // node velocity (y-component)
  double force_x; // node force (x-component)
  double force_y; // node force (y-component)
};

/// Structure for object (either cylinder or red blood cell)

struct particle_struct {

  /// Constructor

  particle_struct() {
    num_nodes = particle_num_nodes;
    radius = particle_radius;
    stiffness = particle_stiffness;
    center.x = particle_center_x;
    center.y = particle_center_y;
    center.x_ref = particle_center_x;
    center.y_ref = particle_center_y;
    node = new node_struct[num_nodes];

    // The initial shape of the object is set in the following.
    // For a cylinder (rigid or deformable), the nodes define a circle.
    // For a red blood cell, the y-position has to be changed in order to describe a red blood cell.
    // Initially, the current node positions and reference node positions are identical.
    // During the simulation, only the current positions are updated,
    // the reference node positions are fixed.

    for(int n = 0; n < num_nodes; ++n) {
      #if defined RIGID_CYLINDER || defined DEFORMABLE_CYLINDER
        node[n].x = center.x + radius * sin(2. * M_PI * (double) n / num_nodes);
        node[n].x_ref = center.x + radius * sin(2. * M_PI * (double) n / num_nodes);
        node[n].y = center.y + radius * cos(2. * M_PI * (double) n / num_nodes);
        node[n].y_ref = center.y + radius * cos(2. * M_PI * (double) n / num_nodes);
      #endif

      #ifdef DEFORMABLE_RBC
        node[n].x = center.x + radius * sin(2. * M_PI * (double) n / num_nodes);
        node[n].x_ref = center.x + radius * sin(2. * M_PI * (double) n / num_nodes);
        node[n].y = radius * cos(2. * M_PI * (double) n / num_nodes);

        // Parametrization of the red blood cell shape in 2D

        if(node[n].y > 0) {
          node[n].y = center.y + sqrt(1 - SQ((center.x - node[n].x) / radius)) * (0.207 + 2.00 * SQ((center.x - node[n].x) / radius) - 1.12 * SQ(SQ((center.x - node[n].x) / radius))) * radius / 2;
          node[n].y_ref = center.y + sqrt(1 - SQ((center.x - node[n].x) / radius)) * (0.207 + 2.00 * SQ((center.x - node[n].x) / radius) - 1.12 * SQ(SQ((center.x - node[n].x) / radius))) * radius / 2;
        }
        else {
          node[n].y = center.y - sqrt(1 - SQ((center.x - node[n].x) / radius)) * (0.207 + 2.00 * SQ((center.x - node[n].x) / radius) - 1.12 * SQ(SQ((center.x - node[n].x) / radius))) * radius / 2;
          node[n].y_ref = center.y - sqrt(1 - SQ((center.x - node[n].x) / radius)) * (0.207 + 2.00 * SQ((center.x - node[n].x) / radius) - 1.12 * SQ(SQ((center.x - node[n].x) / radius))) * radius / 2;
        }
      #endif
    }
  }

  /// Elements

  int num_nodes; // number of surface nodes
  double radius; // object radius
  double stiffness; // stiffness modulus
  node_struct center; // center node
  node_struct *node; // list of nodes
};

/// *****************
/// DECLARE FUNCTIONS
/// *****************

// The following functions are used in the simulation code.

void initialize(); // allocate memory and initialize variables
void LBM(); // perform LBM operations
void momenta(); // compute fluid density and velocity from the populations
void equilibrium(double, double, double); // compute the equilibrium populations from the fluid density and velocity
void compute_particle_forces(particle_struct); // compute the forces acting on the object nodes
void spread(particle_struct); // spread node forces to fluid lattice
void interpolate(particle_struct); // interpolate node velocities from fluid velocity
void update_particle_position(particle_struct); // update object center position
void write_fluid_vtk(int); // write the fluid state to the disk as VTK file
void write_particle_vtk(int, particle_struct); // write the particle state to the disk as VTK file
void write_data(int, particle_struct); // write data to the disk (drag/lift, center position)

/// *************
/// MAIN FUNCTION
/// *************

// This is the main function, containing the simulation initialization and the simulation loop.

int main() {

  /// ************
  /// PREPARATIONS
  /// ************

  initialize(); // allocate memory and initialize variables
  particle_struct particle; // create immersed object

  /// Compute derived quantities

  const double D = Ny - 2; // inner channel diameter
  const double nu = (tau - 0.5) / 3; // lattice viscosity
  const double umax = gravity / (2 * nu) * SQ(0.5 * D); // expected maximum velocity for Poiseuille flow without immersed object
  const double Re = D * umax / nu; // Reynolds number for Poiseuille flow without immersed object

  /// Report derived parameters

  cout << "simulation parameters" << endl;
  cout << "=====================" << endl;
  cout << "D = " << D << endl;
  cout << "nu = " << nu << endl;
  cout << "umax = " << umax << endl;
  cout << "Re = " << Re << endl;
  cout << endl;

  /// ***************
  /// SIMULATION LOOP
  /// ***************

  // Overview of simulation algorithm:
  // 1) compute the node forces based on the object's deformation
  // 2) spread the node forces to the fluid lattice
  // 3) update the fluid state via LBM
  // 4) interpolate the fluid velocity to the object nodes
  // 5) update node positions and object's center
  // 6) if desired, write data to disk and report status

  cout << "starting simulation" << endl;

  for(int t = 1; t <= t_num; ++t) { // run over all times between 1 and t_num

    compute_particle_forces(particle); // compute particle forces
    spread(particle); // spread forces from the Lagrangian to the Eulerian mesh
    LBM(); // perform collision, propagation, and bounce-back
    interpolate(particle); // interpolate velocity
    update_particle_position(particle); // update particle position

    /// Write fluid and particle to VTK files
    // The data is only written each t_info time step.

    if(t % t_disk == 0) {
      write_fluid_vtk(t);
      write_particle_vtk(t, particle);
      write_data(t, particle);
    }

    /// Report end of time step

    if(t % t_info == 0) {
      cout << "completed time step " << t << " in [1, " << t_num << "]" << endl;
    }
  }

  /// Report successful end of simulation

  cout << "simulation complete" << endl;

  return 0;
} // end of main function

/// ****************************************
/// ALLOCATE MEMORY AND INITIALIZE VARIABLES
/// ****************************************

// The memory for lattice variables (populations, density, velocity, forces) is allocated.
// The variables are initialized.

void initialize() {

  /// Create folders, delete data file
  // Make sure that the VTK folders exist.
  // Old file data.dat is deleted, if existing.

  int ignore; // ignore return value of system calls
  ignore = system("mkdir -p vtk_fluid"); // create folder if not existing
  ignore = system("mkdir -p vtk_particle"); // create folder if not existing
  ignore = system("rm -f data.dat"); // delete file if existing

  /// Allocate memory for the fluid density, velocity, and force

  density = new double*[Nx];
  velocity_x = new double*[Nx];
  velocity_y = new double*[Nx];
  force_x = new double*[Nx];
  force_y = new double*[Nx];

  for(int X = 0; X < Nx; ++X) {
    density[X] = new double[Ny];
    velocity_x[X] = new double[Ny];
    velocity_y[X] = new double[Ny];
    force_x[X] = new double[Ny];
    force_y[X] = new double[Ny];
  }

  /// Initialize the fluid density and velocity
  // Start with unit density and zero velocity.

  for(int X = 0; X < Nx; ++X) {
    for(int Y = 0; Y < Ny; ++Y) {
      density[X][Y] = 1;
      velocity_x[X][Y] = 0;
      velocity_y[X][Y] = 0;
      force_x[X][Y] = 0;
      force_y[X][Y] = 0;
    }
  }

  /// Allocate memory for the populations

  pop = new double**[9];
  pop_old = new double**[9];

  for(int c_i = 0; c_i < 9; ++c_i) {
    pop[c_i] = new double*[Nx];
    pop_old[c_i] = new double*[Nx];

    for(int X = 0; X < Nx; ++X) {
      pop[c_i][X] = new double[Ny];
      pop_old[c_i][X] = new double[Ny];

      for(int Y = 0; Y < Ny; ++Y) {
        pop[c_i][X][Y] = 0;
        pop_old[c_i][X][Y] = 0;
      }
    }
  }

  /// Initialize the populations
  // Use the equilibrium populations corresponding to the initialized fluid density and velocity.

  for(int X = 0; X < Nx; ++X) {
    for(int Y = 0; Y < Ny; ++Y) {
      equilibrium(density[X][Y], velocity_x[X][Y], velocity_y[X][Y]);

      for(int c_i = 0; c_i < 9; ++c_i) {
        pop_old[c_i][X][Y] = pop_eq[c_i];
        pop[c_i][X][Y] = pop_eq[c_i];
      }
    }
  }

  return;
}

/// *******************
/// COMPUTE EQUILIBRIUM
/// *******************

// This function computes the equilibrium populations from the fluid density and velocity.
// It computes the equilibrium only at a specific lattice node: Function has to be called at each lattice node.
// The standard quadratic euilibrium is used.
// reminder: SQ(x) = x * x

void equilibrium(double den, double vel_x, double vel_y) {
  pop_eq[0] = weight[0] * den * (1                                                     - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[1] = weight[1] * den * (1 + 3 * (  vel_x        ) + 4.5 * SQ(  vel_x        ) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[2] = weight[2] * den * (1 + 3 * (- vel_x        ) + 4.5 * SQ(- vel_x        ) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[3] = weight[3] * den * (1 + 3 * (          vel_y) + 4.5 * SQ(          vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[4] = weight[4] * den * (1 + 3 * (        - vel_y) + 4.5 * SQ(        - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[5] = weight[5] * den * (1 + 3 * (  vel_x + vel_y) + 4.5 * SQ(  vel_x + vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[6] = weight[6] * den * (1 + 3 * (- vel_x - vel_y) + 4.5 * SQ(- vel_x - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[7] = weight[7] * den * (1 + 3 * (  vel_x - vel_y) + 4.5 * SQ(  vel_x - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
  pop_eq[8] = weight[8] * den * (1 + 3 * (- vel_x + vel_y) + 4.5 * SQ(- vel_x + vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));

  return;
}

/// **********************
/// PERFORM LBM OPERATIONS
/// **********************

void LBM() {

  /// Swap populations
  // The present code used old and new populations which are swapped at the beginning of each time step.
  // This is sometimes called 'double-buffered' or 'ping-pong' algorithm.
  // This way, the old populations are not overwritten during propagation.
  // The resulting code is easier to write and to debug.
  // The memory requirement for the populations is twice as large.

  double ***swap_temp = pop_old;
  pop_old = pop;
  pop = swap_temp;

  /// Lattice Boltzmann equation
  // The lattice Boltzmann equation is solved in the following.
  // The algorithm includes
  // - computation of the lattice force
  // - combined collision and propagation (faster than first collision and then propagation)

  for(int X = 0; X < Nx; ++X) {
    for(int Y = 1; Y < Ny - 1; ++Y) {

      /// Compute lattice force
      // The following code corresponds to Guo's forcing scheme.
      // Gravity is always along the x-axis.

      force_latt[0] = (1 - 0.5 * omega) * weight[0] * (3 * ((   - velocity_x[X][Y]) * (force_x[X][Y] + gravity) + (   - velocity_y[X][Y]) * force_y[X][Y]));
      force_latt[1] = (1 - 0.5 * omega) * weight[1] * (3 * (( 1 - velocity_x[X][Y]) * (force_x[X][Y] + gravity) + (   - velocity_y[X][Y]) * force_y[X][Y]) + 9 * (velocity_x[X][Y]) * (force_x[X][Y] + gravity));
      force_latt[2] = (1 - 0.5 * omega) * weight[2] * (3 * ((-1 - velocity_x[X][Y]) * (force_x[X][Y] + gravity) + (   - velocity_y[X][Y]) * force_y[X][Y]) + 9 * (velocity_x[X][Y]) * (force_x[X][Y] + gravity));
      force_latt[3] = (1 - 0.5 * omega) * weight[3] * (3 * ((   - velocity_x[X][Y]) * (force_x[X][Y] + gravity) + ( 1 - velocity_y[X][Y]) * force_y[X][Y]) + 9 * (velocity_y[X][Y]) * force_y[X][Y]);
      force_latt[4] = (1 - 0.5 * omega) * weight[4] * (3 * ((   - velocity_x[X][Y]) * (force_x[X][Y] + gravity) + (-1 - velocity_y[X][Y]) * force_y[X][Y]) + 9 * (velocity_y[X][Y]) * force_y[X][Y]);
      force_latt[5] = (1 - 0.5 * omega) * weight[5] * (3 * (( 1 - velocity_x[X][Y]) * (force_x[X][Y] + gravity) + ( 1 - velocity_y[X][Y]) * force_y[X][Y]) + 9 * (velocity_x[X][Y] + velocity_y[X][Y]) * (force_x[X][Y] + gravity + force_y[X][Y]));
      force_latt[6] = (1 - 0.5 * omega) * weight[6] * (3 * ((-1 - velocity_x[X][Y]) * (force_x[X][Y] + gravity) + (-1 - velocity_y[X][Y]) * force_y[X][Y]) + 9 * (velocity_x[X][Y] + velocity_y[X][Y]) * (force_x[X][Y] + gravity + force_y[X][Y]));
      force_latt[7] = (1 - 0.5 * omega) * weight[7] * (3 * (( 1 - velocity_x[X][Y]) * (force_x[X][Y] + gravity) + (-1 - velocity_y[X][Y]) * force_y[X][Y]) + 9 * (velocity_x[X][Y] - velocity_y[X][Y]) * (force_x[X][Y] + gravity - force_y[X][Y]));
      force_latt[8] = (1 - 0.5 * omega) * weight[8] * (3 * ((-1 - velocity_x[X][Y]) * (force_x[X][Y] + gravity) + ( 1 - velocity_y[X][Y]) * force_y[X][Y]) + 9 * (velocity_x[X][Y] - velocity_y[X][Y]) * (force_x[X][Y] + gravity - force_y[X][Y]));

      /// Compute equilibrium
      // The equilibrium populations are computed.

      equilibrium(density[X][Y], velocity_x[X][Y], velocity_y[X][Y]);

      /// Compute new populations
      // This is the lattice Boltzmann equation (combined collision and propagation) including external forcing.
      // Periodicity of the lattice in x-direction is taken into account by the %-operator.

      pop[0][X]                [Y]     = pop_old[0][X][Y] * (1 - omega) + pop_eq[0] * omega + force_latt[ 0];
      pop[1][(X + 1) % Nx]     [Y]     = pop_old[1][X][Y] * (1 - omega) + pop_eq[1] * omega + force_latt[ 1];
      pop[2][(X - 1 + Nx) % Nx][Y]     = pop_old[2][X][Y] * (1 - omega) + pop_eq[2] * omega + force_latt[ 2];
      pop[3][X]                [Y + 1] = pop_old[3][X][Y] * (1 - omega) + pop_eq[3] * omega + force_latt[ 3];
      pop[4][X]                [Y - 1] = pop_old[4][X][Y] * (1 - omega) + pop_eq[4] * omega + force_latt[ 4];
      pop[5][(X + 1) % Nx]     [Y + 1] = pop_old[5][X][Y] * (1 - omega) + pop_eq[5] * omega + force_latt[ 5];
      pop[6][(X - 1 + Nx) % Nx][Y - 1] = pop_old[6][X][Y] * (1 - omega) + pop_eq[6] * omega + force_latt[ 6];
      pop[7][(X + 1) % Nx]     [Y - 1] = pop_old[7][X][Y] * (1 - omega) + pop_eq[7] * omega + force_latt[ 7];
      pop[8][(X - 1 + Nx) % Nx][Y + 1] = pop_old[8][X][Y] * (1 - omega) + pop_eq[8] * omega + force_latt[ 8];
    }
  }

  /// Bounce-back
  // Due to the presence of the rigid walls at y = 0 and y = Ny - 1, the populations have to be bounced back.
  // Ladd's momentum correction term is included for moving walls (wall velocity parallel to x-axis).
  // Periodicity of the lattice in x-direction is taken into account via the %-operator.

  for(int X = 0; X < Nx; ++X) {

    /// Bottom wall (y = 0)

    pop[3][X][1] = pop[4][X]                [0];
    pop[5][X][1] = pop[6][(X - 1 + Nx) % Nx][0] + 6 * weight[6] * density[X][1] * wall_vel_bottom;
    pop[8][X][1] = pop[7][(X + 1) % Nx]     [0] - 6 * weight[7] * density[X][1] * wall_vel_bottom;

    /// Top wall (y = Ny - 1)

    pop[4][X][Ny - 2] = pop[3][X]                [Ny - 1];
    pop[6][X][Ny - 2] = pop[5][(X + 1) % Nx]     [Ny - 1] - 6 * weight[5] * density[X][Ny - 2] * wall_vel_top;
    pop[7][X][Ny - 2] = pop[8][(X - 1 + Nx) % Nx][Ny - 1] + 6 * weight[8] * density[X][Ny - 2] * wall_vel_top;
  }

  /// Compute fluid density and velocity
  // The fluid density and velocity are obtained from the populations.

  momenta();

  return;
}

/// **********************************
/// COMPUTE FLUID DENSITY AND VELOCITY
/// **********************************

// This function computes the fluid density and velocity from the populations.
// The velocity correction due to body force is included (Guo's forcing).

void momenta() {
  for(int X = 0; X < Nx; ++X) {
    for(int Y = 1; Y < Ny - 1; ++Y) {
      density[X][Y] = pop[0][X][Y] + pop[1][X][Y] + pop[2][X][Y] + pop[3][X][Y] + pop[4][X][Y] + pop[5][X][Y] + pop[6][X][Y] + pop[7][X][Y] + pop[8][X][Y];
      velocity_x[X][Y] = (pop[1][X][Y] - pop[2][X][Y] + pop[5][X][Y] - pop[6][X][Y] + pop[7][X][Y] - pop[8][X][Y] + 0.5 * (force_x[X][Y] + gravity)) / density[X][Y];
      velocity_y[X][Y] = (pop[3][X][Y] - pop[4][X][Y] + pop[5][X][Y] - pop[6][X][Y] - pop[7][X][Y] + pop[8][X][Y] + 0.5 * force_y[X][Y]) / density[X][Y];
    }
  }

  return;
}

/// ***********************
/// COMPUTE PARTICLE FORCES
/// ***********************

// The forces acting on the object nodes are computed.
// Depending on the simulation type (rigid/deformable cylinder or red blood cell),
// the force computation is different.

void compute_particle_forces(particle_struct particle) {

  /// Reset forces
  // This way, the force from the previous time step is deleted.
  // This is necessary whenever forces are computed using '+='.

  for(int n = 0; n < particle.num_nodes; ++n) {
    particle.node[n].force_x = 0;
    particle.node[n].force_y = 0;
  }

  /// Compute strain forces
  // The strain forces are proportional to the relative displacement of a node with respect to its two neighbors.
  // This force is invariant under displacements and rotations,
  // i.e., it is only sensitive to deformations.
  // Newton's law is fulfilled: sum of forces is zero.
  // In order to make the force distribution smoother, each node is assigned an equivalent area.
  // This area is the circumference of the cylinder divided by the number of surface nodes.
  // WARNING: This is a simple strain model and not necessarily sufficient for accurate simulations.

  #if defined DEFORMABLE_CYLINDER || defined DEFORMABLE_RBC
    const double area = 2 * M_PI * particle.radius / particle.num_nodes; // area belonging to a node

    for(int n = 0; n < particle.num_nodes; ++n) {
      const double distance = sqrt(SQ(particle.node[n].x - particle.node[(n + 1) % particle.num_nodes].x) + SQ(particle.node[n].y - particle.node[(n + 1) % particle.num_nodes].y)); // current distance between neighboring nodes
      const double distance_ref = sqrt(SQ(particle.node[n].x_ref - particle.node[(n + 1) % particle.num_nodes].x_ref) + SQ(particle.node[n].y_ref - particle.node[(n + 1) % particle.num_nodes].y_ref)); // reference distance between neighboring nodes
      const double f_x = particle.stiffness * (distance - distance_ref) * (particle.node[n].x - particle.node[(n + 1) % particle.num_nodes].x);
      const double f_y = particle.stiffness * (distance - distance_ref) * (particle.node[n].y - particle.node[(n + 1) % particle.num_nodes].y);
      particle.node[n].force_x += -f_x;
      particle.node[n].force_y += -f_y;
      particle.node[(n + 1) % particle.num_nodes].force_x += f_x;
      particle.node[(n + 1) % particle.num_nodes].force_y += f_y;
    }
  #endif

  /// Compute bending forces
  // The bending forces are proportional to the angle deviation (current versus reference angle).
  // Angles are defined by three neighboring points (l = left, m = middle, r = right).
  // It is distinguished between convex (positive angles) and concave (negative angles) shapes.
  // Newton's law is fulfilled: sum of forces is zero.
  // WARNING: This is a simple bending model and not necessarily sufficient for accurate simulations.

  #if defined DEFORMABLE_CYLINDER || defined DEFORMABLE_RBC
    for(int n = 0; n < particle.num_nodes; ++n) {

      /// Get node coordinates for bending
      // Three neighboring nodes are required to compute the bending forces.
      // Both their current and their reference positions are taken.

      const double x_l = particle.node[(n - 1 + particle.num_nodes) % particle.num_nodes].x;
      const double y_l = particle.node[(n - 1 + particle.num_nodes) % particle.num_nodes].y;
      const double x_m = particle.node[n].x;
      const double y_m = particle.node[n].y;
      const double x_r = particle.node[(n + 1) % particle.num_nodes].x;
      const double y_r = particle.node[(n + 1) % particle.num_nodes].y;
      const double x_l_ref = particle.node[(n - 1 + particle.num_nodes) % particle.num_nodes].x_ref;
      const double y_l_ref = particle.node[(n - 1 + particle.num_nodes) % particle.num_nodes].y_ref;
      const double x_m_ref = particle.node[n].x_ref;
      const double y_m_ref = particle.node[n].y_ref;
      const double x_r_ref = particle.node[(n + 1) % particle.num_nodes].x_ref;
      const double y_r_ref = particle.node[(n + 1) % particle.num_nodes].y_ref;

      /// Compute normal vector direction
      // The 'normal' vector of the middle node is defined to be normal to the vector connecting the left and the right nodes (tangential vector).
      // It always points in outward direction.
      // Both the current and the reference normals are computed.

      const double tang_x_ref = x_r_ref - x_l_ref;
      const double tang_y_ref = y_r_ref - y_l_ref;
      double normal_x_ref;
      double normal_y_ref;

      if(abs(tang_x_ref) < abs(tang_y_ref)) {
        normal_x_ref = 1;
        normal_y_ref = -tang_x_ref / tang_y_ref;
      }
      else {
        normal_y_ref = 1;
        normal_x_ref = -tang_y_ref / tang_x_ref;
      }

      const double tang_x = x_r - x_l;
      const double tang_y = y_r - y_l;
      double normal_x;
      double normal_y;

      if(abs(tang_x) < abs(tang_y)) {
        normal_x = 1;
        normal_y = -tang_x / tang_y;
      }
      else {
        normal_y = 1;
        normal_x = -tang_y / tang_x;
      }

      /// Normalize normal vector to unit length and outward direction
      // The normal vectors are defined to have unit length and point in outward direction.
      // In order to check its direction, the cross product of the normal and the tangential vector is computed.

      const double normal_length_ref = sqrt(SQ(normal_x_ref) + SQ(normal_y_ref));
      normal_x_ref /= normal_length_ref;
      normal_y_ref /= normal_length_ref;

      if(normal_x_ref * tang_y_ref - normal_y_ref * tang_x_ref > 0) {
        normal_x_ref *= -1;
        normal_y_ref *= -1;
      }

      const double normal_length = sqrt(SQ(normal_x) + SQ(normal_y));
      normal_x /= normal_length;
      normal_y /= normal_length;

      if(normal_x * tang_y - normal_y * tang_x > 0) {
        normal_x *= -1;
        normal_y *= -1;
      }

      /// Compute bending angles
      // The angles (current and reference) are defined by the three points (left, middle, right).
      // The angle is defined to be zero of all points are on one line.
      // Angles are positive for convex shapes (e.g., circle) and negative else.
      // The angle sign has to be checked explicitly.

      double angle_ref_cos = (x_l_ref - x_m_ref) * (x_m_ref - x_r_ref) + (y_l_ref - y_m_ref) * (y_m_ref - y_r_ref);
      angle_ref_cos /= (sqrt(SQ(x_l_ref - x_m_ref) + SQ(y_l_ref - y_m_ref)) * sqrt(SQ(x_m_ref - x_r_ref) + SQ(y_m_ref - y_r_ref)));
      double angle_ref = acos(angle_ref_cos);

      const double convex_x_ref = (x_l_ref + x_r_ref) / 2 - x_m_ref;
      const double convex_y_ref = (y_l_ref + y_r_ref) / 2 - y_m_ref;

      if(convex_x_ref * normal_x_ref + convex_y_ref * normal_y_ref > 0) {
        angle_ref *= -1;
      }

      double angle_cos = (x_l - x_m) * (x_m - x_r) + (y_l - y_m) * (y_m - y_r);
      angle_cos /= (sqrt(SQ(x_l - x_m) + SQ(y_l - y_m)) * sqrt(SQ(x_m - x_r) + SQ(y_m - y_r)));
      double angle = acos(angle_cos);

      const double convex_x = (x_l + x_r) / 2 - x_m;
      const double convex_y = (y_l + y_r) / 2 - y_m;

      if(convex_x * normal_x + convex_y * normal_y > 0) {
        angle *= -1;
      }

      /// Compute force magnitude
      // The forces are proportional to the angle deviation (current minus reference angle).
      // The bending modulus controls the magnitude of the force.
      // All three nodes defining the angle experience a bending force.
      // The total sum of these forces is zero (Newton's law).

      const double force_mag = particle_bending * (angle - angle_ref);
      const double length_l = abs(tang_x * (x_m - x_l) + tang_y * (y_m - y_l));
      const double length_r = abs(tang_x * (x_m - x_r) + tang_y * (y_m - y_r));

      particle.node[(n - 1 + particle.num_nodes) % particle.num_nodes].force_x += normal_x * force_mag * length_l / (length_l + length_r);
      particle.node[(n - 1 + particle.num_nodes) % particle.num_nodes].force_y += normal_y * force_mag * length_l / (length_l + length_r);
      particle.node[n].force_x += -normal_x * force_mag;
      particle.node[n].force_y += -normal_y * force_mag;
      particle.node[(n + 1) % particle.num_nodes].force_x += normal_x * force_mag * length_r / (length_l + length_r);
      particle.node[(n + 1) % particle.num_nodes].force_y += normal_y * force_mag * length_r / (length_l + length_r);
    }
  #endif

  /// Compute rigid forces
  // Here, the node forces are proportional to the displacement with respect to the reference position.

  #ifdef RIGID_CYLINDER
    const double area = 2 * M_PI * particle.radius / particle.num_nodes; // area belonging to a node

    for(int n = 0; n < particle.num_nodes; ++n) {
      particle.node[n].force_x = -particle.stiffness * (particle.node[n].x - particle.node[n].x_ref) * area;
      particle.node[n].force_y = -particle.stiffness * (particle.node[n].y - particle.node[n].y_ref) * area;
    }
  #endif

  return;
}

/// *************
/// SPREAD FORCES
/// *************

// The node forces are spread to the fluid nodes via IBM.
// The two-point interpolation stencil (bi-linear interpolation) is used in the present code.
// It may be replaced by a higher-order interpolation.

void spread(particle_struct particle) {

  /// Reset forces
  // This is necessary since '+=' is used afterwards.

  for(int X = 0; X < Nx; ++X) {
    for(int Y = 1; Y < Ny - 1; ++Y) {
      force_x[X][Y] = 0;
      force_y[X][Y] = 0;
    }
  }

  /// Spread forces
  // Run over all object nodes.

  for(int n = 0; n < particle.num_nodes; ++n) {

    // Identify the lowest fluid lattice node in interpolation range.
    // 'Lowest' means: its x- and y-values are the smallest.
    // The other fluid nodes in range have coordinates
    // (x_int + 1, y_int), (x_int, y_int + 1), and (x_int + 1, y_int + 1).

    int x_int = (int) (particle.node[n].x - 0.5 + Nx) - Nx;
    int y_int = (int) (particle.node[n].y + 0.5);

    // Run over all neighboring fluid nodes.
    // In the case of the two-point interpolation, it is 2x2 fluid nodes.

    for(int X = x_int; X <= x_int + 1; ++X) {
      for(int Y = y_int; Y <= y_int + 1; ++Y) {

        // Compute distance between object node and fluid lattice node.

        const double dist_x = particle.node[n].x - 0.5 - X;
        const double dist_y = particle.node[n].y + 0.5 - Y;

        // Compute interpolation weights for x- and y-direction based on the distance.

        const double weight_x = 1 - abs(dist_x);
        const double weight_y = 1 - abs(dist_y);

        // Compute lattice force.

        force_x[(X + Nx) % Nx][Y] += (particle.node[n].force_x * weight_x * weight_y);
        force_y[(X + Nx) % Nx][Y] += (particle.node[n].force_y * weight_x * weight_y);
      }
    }
  }

  return;
}

/// **********************
/// INTERPOLATE VELOCITIES
/// **********************

// The node velocities are interpolated from the fluid nodes via IBM.
// The two-point interpolation stencil (bi-linear interpolation) is used in the present code.
// It may be replaced by a higher-order interpolation.

void interpolate(particle_struct particle) {

  // Run over all object nodes.

  for(int n = 0; n < particle.num_nodes; ++n) {

    // Reset node velocity first since '+=' is used.

    particle.node[n].vel_x = 0;
    particle.node[n].vel_y = 0;

    // Identify the lowest fluid lattice node in interpolation range (see spreading).

    int x_int = (int) (particle.node[n].x - 0.5 + Nx) - Nx;
    int y_int = (int) (particle.node[n].y + 0.5);

    // Run over all neighboring fluid nodes.
    // In the case of the two-point interpolation, it is 2x2 fluid nodes.

    for(int X = x_int; X <= x_int + 1; ++X) {
      for(int Y = y_int; Y <= y_int + 1; ++Y) {

        // Compute distance between object node and fluid lattice node.

        const double dist_x = particle.node[n].x - 0.5 - X;
        const double dist_y = particle.node[n].y + 0.5 - Y;

        // Compute interpolation weights for x- and y-direction based on the distance.

        const double weight_x = 1 - abs(dist_x);
        const double weight_y = 1 - abs(dist_y);

        // Compute node velocities.

        particle.node[n].vel_x += (velocity_x[(X + Nx) % Nx][Y] * weight_x * weight_y);
        particle.node[n].vel_y += (velocity_y[(X + Nx) % Nx][Y] * weight_x * weight_y);
      }
    }
  }

  return;
}

/// ************************
/// UPDATE PARTICLE POSITION
/// ************************

// The position of the particle nodes are updated according to their velocity.
// The center position is updated as well.
// The new node position is its old position plus its current velocity (Euler integration).
// The center position is the arithmetic mean of all node positions.
// Periodicity is taken into account:
// If the particle center leaves the system domain (x < 0 or x >= Nx), it reenters from the other side.

void update_particle_position(particle_struct particle) {

  /// Reset center position

  particle.center.x = 0;
  particle.center.y = 0;

  /// Update node and center positions

  for(int n = 0; n < particle.num_nodes; ++n) {
    particle.node[n].x += particle.node[n].vel_x;
    particle.node[n].y += particle.node[n].vel_y;
    particle.center.x += particle.node[n].x / particle.num_nodes;
    particle.center.y += particle.node[n].y / particle.num_nodes;
  }

  /// Check for periodicity along the x-axis

  if(particle.center.x < 0) {
    particle.center.x += Nx;

    for(int n = 0; n < particle.num_nodes; ++n) {
      particle.node[n].x += Nx;
    }
  }
  else if(particle.center.x >= Nx) {
    particle.center.x -= Nx;

    for(int n = 0; n < particle.num_nodes; ++n) {
      particle.node[n].x -= Nx;
    }
  }

  return;
}

/// *****************************
/// WRITE FLUID STATE TO VTK FILE
/// *****************************

// The fluid state is writen to a VTK file at each t_disk step.
// The following data is written:
// - density difference (density - 1)
// - x-component of velocity
// - y-component of velocity
// The following code is designed in such a way that the file can be read by ParaView.

void write_fluid_vtk(int time) {

  /// Create filename

  stringstream output_filename;
  output_filename << "vtk_fluid/fluid_t" << time << ".vtk";
  ofstream output_file;

  /// Open file

  output_file.open(output_filename.str().c_str());

  /// Write VTK header

  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "fluid_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET RECTILINEAR_GRID\n";
  output_file << "DIMENSIONS " << Nx << " " << Ny - 2 << " 1" << "\n";
  output_file << "X_COORDINATES " << Nx << " float\n";

  for(int X = 0; X < Nx; ++X) {
    output_file << X + 0.5 << " ";
  }

  output_file << "\n";
  output_file << "Y_COORDINATES " << Ny - 2 << " float\n";

  for(int Y = 1; Y < Ny - 1; ++Y) {
    output_file << Y - 0.5 << " ";
  }

  output_file << "\n";
  output_file << "Z_COORDINATES " << 1 << " float\n";
  output_file << 0 << "\n";
  output_file << "POINT_DATA " << Nx * (Ny - 2) << "\n";

  /// Write density difference

  output_file << "SCALARS density_difference float 1\n";
  output_file << "LOOKUP_TABLE default\n";

  for(int Y = 1; Y < Ny - 1; ++Y) {
    for(int X = 0; X < Nx; ++X) {
      output_file << density[X][Y] - 1 << "\n";
    }
  }

  /// Write velocity

  output_file << "VECTORS velocity_vector float\n";

  for(int Y = 1; Y < Ny - 1; ++Y) {
    for(int X = 0; X < Nx; ++X) {
      output_file << velocity_x[X][Y] << " " << velocity_y[X][Y] << " 0\n";
    }
  }

  /// Close file

  output_file.close();

  return;
}

/// ********************************
/// WRITE PARTICLE STATE TO VTK FILE
/// ********************************

// The particle state (node positions) is writen to a VTK file at each t_disk step.
// The following code is designed in such a way that the file can be read by ParaView.

void write_particle_vtk(int time, particle_struct particle) {

  /// Create filename

  stringstream output_filename;
  output_filename << "vtk_particle/particle_t" << time << ".vtk";
  ofstream output_file;

  /// Open file

  output_file.open(output_filename.str().c_str());

  /// Write VTK header

  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "particle_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET POLYDATA\n";

  /// Write node positions

  output_file << "POINTS " << particle_num_nodes << " float\n";

  for(int n = 0; n < particle_num_nodes; ++n) {
    output_file << particle.node[n].x << " " << particle.node[n].y << " 0\n";
  }

  /// Write lines between neighboring nodes

  output_file << "LINES " << particle_num_nodes << " " << 3 * particle_num_nodes << "\n";

  for(int n = 0; n < particle_num_nodes; ++n) {
    output_file << "2 " << n << " " << (n + 1) % particle_num_nodes << "\n";
  }

  /// Write vertices

  output_file << "VERTICES 1 " << particle_num_nodes + 1 << "\n";
  output_file << particle_num_nodes << " ";

  for(int n = 0; n < particle_num_nodes; ++n) {
    output_file << n << " ";
  }

  /// Close file

  output_file.close();

  return;
}

/// ************************
/// WRITE DATA TO ASCII FILE
/// ************************

// The following quantities are written to the disk at each t_disk step:
// - drag and lift forces (x- and y-components of the force)
// - object center position (x- and y-components)
// - object center velocity (x- and y-components)
// The data file is readable by gnuplot

void write_data(int time, particle_struct particle) {

  /// Create filename

  string output_filename("data.dat");
  ofstream output_file;

  /// Open file

  output_file.open(output_filename.c_str(), fstream::app);

  /// Compute quantities

  double force_tot_x = 0;
  double force_tot_y = 0;
  double vel_center_x = 0;
  double vel_center_y = 0;

  for(int i = 0; i < particle.num_nodes; ++i) {
    force_tot_x += particle.node[i].force_x;
    force_tot_y += particle.node[i].force_y;
    vel_center_x += particle.node[i].vel_x;
    vel_center_y += particle.node[i].vel_y;
  }

  /// Write data

  output_file << time << " "; // time step
  output_file << force_tot_x << " "; // drag force
  output_file << force_tot_y << " "; // lift force
  output_file << particle.center.x << " "; // center position (x-component)
  output_file << particle.center.y << " "; // center position (y-component)
  output_file << vel_center_x << " "; // center velocity (x-component)
  output_file << vel_center_y << "\n"; // center velocity (y-component)

  /// Close file

  output_file.close();

  return;
}
