#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/timeb.h>

const double pi = 3.14159265358979323846264338328;

long int klok;

FILE * f_trace;
FILE * f_atoms;
char file_error [200] = "Error.dat";

//MPI parameters
int my_rank, comm_size;

//Box parameters
double side_X, side_Y, side_Z, side_XYZ;
double side_hX, side_hY, side_hZ;

//General simulation parameters
double h = 0.05;
double hsq2 = 0.5 * h * h;
double T;//temperature

//Particle numbers
int N_amino = 0; //the number of residues in the protein
int NBP_atom = 0; //the number of backbone atoms in the protein 
int NB_atom = 0; //the total number of backbone atoms in the protein and crowders
int * NS_atom = NULL; //the number of sites per backbone atom or crowder ("NS_atom = m_crwd - 1" for crowders)
int NTP_atom = 0; //the total number of atoms in RNA
int NpTP_atom = 0; //the number of different pairs of atoms in RNA
int NT_atom = 0; //the total number of atoms in RNA and crowders

int ** atom_key = NULL;
char * amino_key = NULL;
int * crowder_key = NULL;
int * part_key = NULL;
int * maxi_key = NULL;

int * INDX = NULL, * JNDX = NULL;

int * residue_mass = NULL;
int ** residue_content = NULL;
int * crowder_mass = NULL;
int ** crowder_content = NULL;

int * RIGID_SET = NULL, * RIGID_END = NULL;

//Individual properties of atoms

const int ns = 6; //the number of different atom types
const int nps = ns * (ns + 1) / 2; //the number of different pairs of atoms
const int npi = 2; //the number of different interactions which require a list
                           //LJ (WCA) -> InInt = 0; Coulomb -> InInt = 1 

const double D [ns] = { 4.2, 5.8, 6.0, 1.5852, 3.8960, 5.3160 }; //phosphate, RIBOSE, BASE, Mg, Cl, K
double PD [nps]; 
double PD_sq [nps];

/////////////////////////////////////
const int na = 9; //the number of different "amber" atom types
const int npa = na * (na + 1) / 2; //the number of different pairs of amber atoms

double MASS [na + 1]; //atom masses
double RLJ [na + 1]; //atom LJ radii
double ELJ [na + 1]; //atom LJ well depths
double QC [na + 1]; //atom charges
double VISC [na + 1];
double K1 [na + 1];
double K2 [na + 1];
double SIGMA_FORCE [na + 1];
double D_LJ [npa + 1]; //LJ pair diameters
double D2_LJ [npa + 1]; //LJ square pair diameters
double D2_LJ_CUTOFF [npa + 1]; //LJ square cutoff distances
double D2_LJ_OVERLAP [npa + 1]; //LJ square overlap distances
double D2_LJ_MINIMUM [npa + 1]; //LJ square minimum distances

double E_LJ [npa + 1]; //LJ energy prefactors
double F_LJ [npa + 1]; //LJ force prefactors
double Q_C [npa + 1]; //Coulomb prefactors

//Crowders
const int n_crwd = 3; //the number of different crowder types
double cM_crwd [ n_crwd ] = { 0.0, 0.0, 0.0 }; //molar concentrations of crowders
double phi_crwd [ n_crwd ] = { 0.0, 0.0, 0.0 }; //volume fractions of crowders
double V_crwd [ n_crwd ]; //crowder volumes
int N_crwd [ n_crwd ]; //the number of crowders of each type
int m_crwd [ n_crwd ]; //the number of sites in a crowder of each type 
double * RX_crwd [ n_crwd ]; //x coordinates of sites in a crowder of each type
double * RY_crwd [ n_crwd ]; //y coordinates of sites in a crowder of each type
double * RZ_crwd [ n_crwd ]; //z coordinates of sites in a crowder of each type
int * part_key_crwd [ n_crwd ]; //site (atom) "part_key" in a crowder of each type
int * maxi_key_crwd [ n_crwd ]; //site (atom) "maxi_key" in a crowder of each type
int RIGID_SET_crwd [ n_crwd ]; //the first site in a crowder rigid sub-unit
int RIGID_END_crwd [ n_crwd ]; //the last site in a crowder rigid sub-unit

//Coordinates, velocities and forces

double * x = NULL;
double * y = NULL; 
double * z = NULL;

double * x_old = NULL;
double * y_old = NULL; 
double * z_old = NULL;

double * vx = NULL;
double * vy = NULL; 
double * vz = NULL;

double * fx = NULL;
double * fy = NULL; 
double * fz = NULL;

double * flx = NULL;
double * fly = NULL; 
double * flz = NULL;

double * fsend = NULL;
double * freceive = NULL;

double maxwell_force [3001];

//Centre of mass motion
double * CMS = NULL;
double * VCMS = NULL;
double * RMASS = NULL;
double * RVISC = NULL;
double * RXYZ = NULL;

//Rotational motion
double * AXES = NULL;
double * W = NULL;
double ** RX = NULL, ** RY = NULL, ** RZ = NULL;
double * IR1 = NULL, * IR2 = NULL, * IR3 = NULL;
double * AV = NULL, * BV = NULL, * CV = NULL;
double * FV = NULL, * GV = NULL, * HV = NULL;

//////////////////////////////////////////////////////////////////////////
/////////////Cells (for list population or pcfs calculation)//////////////

int McX = 0, McY = 0, McZ = 0, McT = 0; //cell numbers

int ** cell_content = NULL, * cell_mass = NULL; //cell content and mass

//////////////////////////////////////////////////////////////////////////
/////////////////////////////////Lists////////////////////////////////////

int * list_content [npi];

int list_mass [npi];

double R1_LIST [npi];

double R2_LIST [npi]; 

double PROGRESS;

double dr1;

//Stacks
double st_D = 0; //adjustment for tertiary stacks
double ss_D = 0; //adjustment for secondary stacks

int s3_N = 0; //number of tertiary stacks
int st_N = 0; //total number of stacks

int * st_i = NULL;
int * st_j = NULL;

double * st_r = NULL;
double * st_theta1 = NULL;
double * st_theta2 = NULL;
double * st_psi = NULL;
double * st_psi1 = NULL;
double * st_psi2 = NULL;
double * st_E = NULL;

double * st_energy = NULL;
int * st_status = NULL;
int * ST_EXCESS = NULL;
int ** ATOM_ST = NULL;

//Hydrogen bonds
double hs_D = 1.0;
int hb_N = 0;
int * hb_dode = NULL;
int * hb_k10 = NULL;
int * hb_k11 = NULL;
int * hb_k12 = NULL;
int * hb_k20 = NULL;
int * hb_k21 = NULL;
int * hb_k22 = NULL;
double * hb_r = NULL;
double * hb_theta1 = NULL;
double * hb_theta2 = NULL;
double * hb_psi = NULL;
double * hb_psi1 = NULL;
double * hb_psi2 = NULL;
double * hb_E = NULL;
double * hb_psi0 = NULL;
double * hb_psi10 = NULL;
double * hb_psi20= NULL;
double * hb_energy = NULL;
int * hb_status = NULL;
int * hb_RES1 = NULL;
int * hb_RES2 = NULL;
char ** hb_code = NULL;
int * HB_NS_atom = NULL; //the number of hydrogen bonding atoms in a "sugar + base" or a phosphate
int HB_NT_atom = 0; //the total number of hydrogen bonding atoms
int ** HB_atom_key = NULL;
int * HB_INDX = NULL, * HB_JNDX = NULL;
int * VALENCE = NULL;
int * HB_EXCESS = NULL;
int * EXCESS = NULL;
int * HB_ATOM_N = NULL;
int ** HB_ATOM = NULL;
int * ATOM_HB_N = NULL;
int ** ATOM_HB = NULL;
int * HB_PAIR_N = NULL;
int ** HB_PAIR = NULL;
int * HB_A = NULL; 
int * hb_K = NULL;

///////////Electrostatics/////////
double lB;

double D_CT_CUTOFF = 50.0;
double D2_CT_CUTOFF = D_CT_CUTOFF * D_CT_CUTOFF;

//Energies
double E2_HB = 0;
double E3_HB = 0;
double E2_STACK = 0;
double E3_STACK = 0;
double ELJ_NON_NATIVE = 0;
double E_YUKAWA = 0;
double E_BOND_B = 0;
double E_BOND_CRWD = 0;
double E_VALENCE_B = 0;
double E_VALENCE_CRWD = 0;
double E_DIHEDRAL_B = 0;
double E_DIHEDRAL_CRWD = 0;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////AUXILIARY ROUTINES//////////////////////////////////////////////////////

int Kdelta ( int i, int j )
{
      if (i == j) return (1);

      else return (0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void generate_unit_vector ( double * n1 )
{
      double a1;

      /////////////////////////////////////////////////
      
      n1 [0] = (double) rand (); 
            
      a1 = (double) rand () / RAND_MAX; 
            
      if (a1 > 0.5) n1 [0] = - n1 [0];

      /////////////////////////////////////////////////
            
      n1 [1] = (double) rand (); 
            
      a1 = (double) rand () / RAND_MAX; 
            
      if (a1 > 0.5) n1 [1] = - n1 [1];

      /////////////////////////////////////////////////
            
      n1 [2] = (double) rand (); 
            
      a1 = (double) rand () / RAND_MAX; 
            
      if (a1 > 0.5) n1 [2] = - n1 [2];

      /////////////////////////////////////////////////
            
      a1 = sqrt ( n1 [0] * n1 [0] + n1 [1] * n1 [1] + n1 [2] * n1 [2] );

      n1 [0] /= a1;

      n1 [1] /= a1;

      n1 [2] /= a1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void generate_perpendicular_vector ( double * n3, double * res )
{
      double p;
      
      
      res [0] = (double) rand () / RAND_MAX; p = (double) rand () / RAND_MAX; if ( p > 0.5 ) res [0] = - res [0];

      res [1] = (double) rand () / RAND_MAX; p = (double) rand () / RAND_MAX; if ( p > 0.5 ) res [1] = - res [1];
            
      res [2] = (double) rand () / RAND_MAX; p = (double) rand () / RAND_MAX; if ( p > 0.5 ) res [2] = - res [2];

      
      p = res [0] * n3 [0] + res [1] * n3 [1] + res [2] * n3 [2];

      res [0] = res [0] - p * n3 [0]; res [1] = res [1] - p * n3 [1]; res [2] = res [2] - p * n3 [2];

      
      p = sqrt (res [0] * res [0] + res [1] * res [1] + res [2] * res [2]);

      res [0] = res [0] / p; res [1] = res [1] / p; res [2] = res [2] / p;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void cross_product ( double * v1, double * v2, double * res )
{
      res [0] = v1 [1] * v2 [2] - v1 [2] * v2 [1];

      res [1] = v1 [2] * v2 [0] - v1 [0] * v2 [2];

      res [2] = v1 [0] * v2 [1] - v1 [1] * v2 [0];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*void maxwell_1 ( double * array, int array_size )
{
      double r1, r2, tmp1, tmp2;

      int i, j;

      if ( (array_size % 2) == 0 ) j = array_size; 
            
      else 
      {
            j = array_size - 1;

gen1: 
            
            r1 = (double) rand () / RAND_MAX; if ( r1 == 0 ) goto gen1;
      
            r2 = (double) rand () / RAND_MAX;
      
            tmp1 = sqrt ( - 2.0 * log (r1) ); tmp2 = 2.0 * pi * r2;

            array [ array_size ] = tmp1 * cos (tmp2); 
      }
      
      for (i = 1; i < j + 1; i+= 2)
      {

gen2: 
      
      r1 = (double) rand () / RAND_MAX; if ( r1 == 0 ) goto gen2;
      
      r2 = (double) rand () / RAND_MAX;
      
      tmp1 = sqrt ( - 2.0 * log (r1) ); tmp2 = 2.0 * pi * r2;

      array [i] = tmp1 * cos (tmp2); array [i + 1] = tmp1 * sin (tmp2);
      
      }
}*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void generate_atom_velocity ( int i, int j )
{
      int k, l;

      double a1, a2, tmp1, tmp2;

      /////////////////////////////////////////////////////////////////
      
      k = atom_key [i][j];

      l = maxi_key [k];

      /////////////////////////////////////////////////////////////////

genvx:
      
      a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvx;
                        
      a2 = (double) rand () / RAND_MAX;
      
      tmp1 = sqrt ( T / MASS [l] ) * sqrt ( -2.0 * log (a1) );
      
      tmp2 = 2.0 * pi * a2;

      vx [k] = tmp1 * cos ( tmp2 );

genvy:
      
      a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvy;
                        
      a2 = (double) rand () / RAND_MAX;
      
      tmp1 = sqrt ( T / MASS [l] ) * sqrt ( -2.0 * log (a1) );
                        
      tmp2 = 2.0 * pi * a2;
                        
      vy [k] = tmp1 * cos ( tmp2 );

genvz:
      
      a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvz;
                        
      a2 = (double) rand () / RAND_MAX;
                        
      tmp1 = sqrt ( T / MASS [l] ) * sqrt ( -2.0 * log (a1) );
                        
      tmp2 = 2.0 * pi * a2;
                        
      vz [k] = tmp1 * cos ( tmp2 );

      /////////////////////////////////////////////////////////////////

      x_old [k] = x [k] - vx [k] * h;
      
      y_old [k] = y [k] - vy [k] * h;

      z_old [k] = z [k] - vz [k] * h;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void leap_frog_move_atom ( int k )
{
      int l;
      
      double x_tmp, y_tmp, z_tmp;

      double vx_tmp, vy_tmp, vz_tmp;

      double vsq;
      
      /////////////////////////////////////////////////////////////////
      
      l = maxi_key [k];

      /////////////////////////////////////////////////////////////////

      x_tmp = x [k]; 
      
      y_tmp = y [k]; 
      
      z_tmp = z [k];

      /////////////////////////////////////////////////////////////////
                  
      vx_tmp = vx [k]; 
      
      vy_tmp = vy [k]; 
      
      vz_tmp = vz [k];

      /////////////////////////////////////////////////////////////////
      
      vx [k] = K1 [l] * vx_tmp + K2 [l] * fx [k];

      vy [k] = K1 [l] * vy_tmp + K2 [l] * fy [k];

      vz [k] = K1 [l] * vz_tmp + K2 [l] * fz [k];

      /////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////

      vsq = MASS [l] * ( vx [k] * vx [k] + vy [k] * vy [k] + vz [k] * vz [k] ) / (300.0 * T);

      if ( vsq > 1.0 ) 
      {
            vsq = 1.0 / sqrt ( vsq );

            vx [k] *= vsq;

            vy [k] *= vsq;

            vz [k] *= vsq;
      }

      /////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////
      
      x [k] = 2.0 * x_tmp - x_old [k] + ( vx [k] - vx_tmp ) * h;

      y [k] = 2.0 * y_tmp - y_old [k] + ( vy [k] - vy_tmp ) * h;

      z [k] = 2.0 * z_tmp - z_old [k] + ( vz [k] - vz_tmp ) * h;

      /////////////////////////////////////////////////////////////////
                  
      x_old [k] = x_tmp; 
      
      y_old [k] = y_tmp; 
      
      z_old [k] = z_tmp;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void half_shift ( double * xx )
{
      if ( xx [0] > side_hX ) xx [0] = xx [0] - side_X; if ( xx [0] < - side_hX) xx [0] = xx [0] + side_X;

      if ( xx [1] > side_hY ) xx [1] = xx [1] - side_Y; if ( xx [1] < - side_hY) xx [1] = xx [1] + side_Y;
      
      if ( xx [2] > side_hZ ) xx [2] = xx [2] - side_Z; if ( xx [2] < - side_hZ) xx [2] = xx [2] + side_Z;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CC_distance ( int i1, int i2 )
{
      double xx [3], r2;
      
      xx [0] = x [i1] - x [i2]; 
                  
      xx [1] = y [i1] - y [i2]; 
                  
      xx [2] = z [i1] - z [i2]; 
      
      r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

      return ( sqrt (r2) );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CC_vector ( int i1, int i2, double * xx )
{
      xx [0] = x [i1] - x [i2]; 
                  
      xx [1] = y [i1] - y [i2]; 
                  
      xx [2] = z [i1] - z [i2]; 
      
      half_shift (xx);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CC_energy_single_pair ( int i1, int i2, int Jn12, double scale_factor, double ( * energy) ( double x2, int Jn12 ), double * observable )
{
      double r2, xx [3];
      
      /////////////////////////////////////////////////////////////////

      xx [0] = x [i1] - x [i2]; 
                  
      xx [1] = y [i1] - y [i2]; 
                  
      xx [2] = z [i1] - z [i2]; 
      
      half_shift (xx);

      /////////////////////////////////////////////////////////////////
      
      r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];
      
      (* observable) += scale_factor * energy ( r2, Jn12 );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CC_force_single_pair ( int i1, int i2, int Jn12, double scale_factor, double ( * force) ( double x2, int Jn12 ) )
{
      double r2, xx [3];
      
      /////////////////////////////////////////////////////////////////

      xx [0] = x [i1] - x [i2]; 
                  
      xx [1] = y [i1] - y [i2]; 
                  
      xx [2] = z [i1] - z [i2]; 
      
      half_shift (xx);

      /////////////////////////////////////////////////////////////////
      
      r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

      r2 = scale_factor * force ( r2, Jn12 );
      
      /////////////////////////////////////////////////////////////////
      
      xx [0] = xx [0] * r2; 
      
      xx [1] = xx [1] * r2; 
      
      xx [2] = xx [2] * r2;

      /////////////////////////////////////////////////////////////////
      
      fx [i1] += xx [0]; fx [i2] -= xx [0]; 
                  
      fy [i1] += xx [1]; fy [i2] -= xx [1]; 
                  
      fz [i1] += xx [2]; fz [i2] -= xx [2];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double WCA_ENERGY ( double x2, int Jn12 )
{
      if ( x2 > D2_LJ_CUTOFF [ Jn12 ] ) return (0);

      else if ( x2 < D2_LJ_MINIMUM [ Jn12 ] ) return ( E_LJ [ Jn12 ] * 1001.0 );
      
      else
      {
            double dr, r1, r2, r4, r6;
            
            r1 = sqrt ( x2 );
            
            dr = r1 + 1.5852 - D_LJ [ Jn12 ];

            r4 = 1.5852 / dr; 
            
            r2 = r4 * r4; r4 = r2 * r2; r6 = r4 * r2; r2 = r6 * r6;
                  
            r4 = E_LJ [ Jn12 ] * ( r2 - 2.0 * r6 + 1.0 );
            
            return ( r4 );
      }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double LJ_ENERGY ( double x2, int Jn12  )
{
      if ( x2 > D2_LJ_CUTOFF [ Jn12 ] ) return (0);

      else if ( x2 < D2_LJ_MINIMUM [ Jn12 ] ) return ( E_LJ [ Jn12 ] * 1000.0 );
      
      else
      {
            double dr, r1, r2, r4, r6;
            
            r1 = sqrt ( x2 );
            
            dr = r1 + 1.5852 - D_LJ [ Jn12 ];

            r4 = 1.5852 / dr; 
            
            r2 = r4 * r4; r4 = r2 * r2; r6 = r4 * r2; r2 = r6 * r6;
                  
            r4 = E_LJ [ Jn12 ] * ( r2 - 2.0 * r6 );
            
            return ( r4 );
      }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double LJ_FORCE ( double x2, int Jn12 )
{
      if ( x2 > D2_LJ_CUTOFF [ Jn12 ] ) return (0);
      
      else
      {
            double dr, r1, r2, r4, r8, r14;

            r1 = sqrt ( x2 );
            
            dr = r1 + 1.5852 - D_LJ [ Jn12 ];

            r4 = 1.5852 / dr; 
            
            r2 = r4 * r4; r4 = r2 * r2; r8 = r4 * r4; r14 = r8 * r4 * r2;
                  
            r4 = F_LJ [ Jn12 ] * ( r14 - r8 ) * dr / r1;
                  
            return ( r4 );
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double COULOMB_TRUNCATED_ENERGY ( double r2, int Jn12 )
{
      if ( r2 > D2_CT_CUTOFF ) return (0);

      else 
      {
            double r1;
            
            r1 = sqrt (r2); 
            
            return ( Q_C [ Jn12 ] * lB / r1 );
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double COULOMB_TRUNCATED_FORCE ( double r2, int Jn12 )
{
      if ( r2 > D2_CT_CUTOFF ) return (0);

      else
      {
            double r1, r3;
            
            r1 = sqrt (r2); 

            r3 = r1 * r2;
            
            return ( Q_C [ Jn12 ] * lB / r3 );
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void reduce_energies ( int i1, int i2, double scale_factor )
{
      int Jn1, Jn2, Jn12;

      /////////////////////////////////////////////////////////////////
      
      Jn1 = maxi_key [i1];

      Jn2 = maxi_key [i2];

      if (Jn1 < Jn2) Jn12 = (Jn1 - 1) * na - (Jn1 - 1) * Jn1 / 2 + Jn2;

      else Jn12 = (Jn2 - 1) * na - (Jn2 - 1) * Jn2 / 2 + Jn1;

      /////////////////////////////////////////////////////////////////
      
      CC_energy_single_pair ( i1, i2, Jn12, scale_factor, WCA_ENERGY, &ELJ_NON_NATIVE );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void reduce_forces ( int i1, int i2, double scale_factor )
{
      int Jn1, Jn2, Jn12;

      /////////////////////////////////////////////////////////////////
      
      Jn1 = maxi_key [i1];

      Jn2 = maxi_key [i2];

      if (Jn1 < Jn2) Jn12 = (Jn1 - 1) * na - (Jn1 - 1) * Jn1 / 2 + Jn2;

      else Jn12 = (Jn2 - 1) * na - (Jn2 - 1) * Jn2 / 2 + Jn1;

      /////////////////////////////////////////////////////////////////
      
      CC_force_single_pair ( i1, i2, Jn12, scale_factor, LJ_FORCE );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double bond_length ( int i1, int i2 )
{
      double v1 [3], r1;

      /////////////////////////////////////////////////////////////////

      v1 [0] = x [i2] - x [i1];

      v1 [1] = y [i2] - y [i1];

      v1 [2] = z [i2] - z [i1];

      /////////////////////////////////////////////////////////////////

      r1 = sqrt ( v1 [0] * v1 [0] + v1 [1] * v1 [1] + v1 [2] * v1 [2] );

      /////////////////////////////////////////////////////////////////
      
      return ( r1 );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void bead_bond_force ( int i1, int i2, double k0, double r0 )
{
      double v1 [3], r1, tmp1;
      
      /////////////////////////////////////////////////////////////////
      
      v1 [0] = x [i2] - x [i1];

      v1 [1] = y [i2] - y [i1];

      v1 [2] = z [i2] - z [i1];

      /////////////////////////////////////////////////////////////////

      r1 = sqrt ( v1 [0] * v1 [0] + v1 [1] * v1 [1] + v1 [2] * v1 [2] );

      tmp1 = 2.0 * k0 * (r0 - r1) / r1;
      
      /////////////////////////////////////////////////////////////////

      r1 = tmp1 * v1 [0]; fx [i1] -= r1; fx [i2] += r1;

      r1 = tmp1 * v1 [1]; fy [i1] -= r1; fy [i2] += r1;
      
      r1 = tmp1 * v1 [2]; fz [i1] -= r1; fz [i2] += r1;

      /////////////////////////////////////////////////////////////////

      reduce_forces ( i1, i2, -1.0 );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void stacking_bond_force ( int i1, int i2, double k0, double r0 )
{
      double v1 [3], r1, tmp1;
      
      /////////////////////////////////////////////////////////////////
      
      v1 [0] = x [i2] - x [i1];

      v1 [1] = y [i2] - y [i1];

      v1 [2] = z [i2] - z [i1];

      /////////////////////////////////////////////////////////////////

      r1 = sqrt ( v1 [0] * v1 [0] + v1 [1] * v1 [1] + v1 [2] * v1 [2] );

      tmp1 = 2.0 * k0 * (r0 - r1) / r1;
      
      /////////////////////////////////////////////////////////////////

      r1 = tmp1 * v1 [0]; fx [i1] -= r1; fx [i2] += r1;

      r1 = tmp1 * v1 [1]; fy [i1] -= r1; fy [i2] += r1;
      
      r1 = tmp1 * v1 [2]; fz [i1] -= r1; fz [i2] += r1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double valence_angle ( int i1, int i2, int i3 )
{
      double v1 [3], v2 [3], v1_norm, v2_norm, v12_norm, kosinus, theta;

      /////////////////////////////////////////////////////////////////

      v1 [0] = x [i1] - x [i2];

      v1 [1] = y [i1] - y [i2];

      v1 [2] = z [i1] - z [i2];

      /////////////////////////////////////////////////////////////////

      v2 [0] = x [i3] - x [i2];

      v2 [1] = y [i3] - y [i2];

      v2 [2] = z [i3] - z [i2];

      /////////////////////////////////////////////////////////////////

      v1_norm = v1 [0] * v1 [0] + v1 [1] * v1 [1] + v1 [2] * v1 [2];

      /////////////////////////////////////////////////////////////////
      
      v2_norm = v2 [0] * v2 [0] + v2 [1] * v2 [1] + v2 [2] * v2 [2];

      /////////////////////////////////////////////////////////////////
      
      v12_norm = sqrt ( v1_norm * v2_norm );

      /////////////////////////////////////////////////////////////////
      
      kosinus = ( v1 [0] * v2 [0] + v1 [1] * v2 [1] + v1 [2] * v2 [2] ) / v12_norm;

      theta = acos ( kosinus );

      /////////////////////////////////////////////////////////////////
      
      return ( theta );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void bead_valence_force ( int i1, int i2, int i3, double k0, double theta0 )
{
      double v1 [3], v2 [3], v1_norm, v2_norm, v12_norm, kosinus, sinus, theta, tmp1, tmp2, tmp3, C1, C2;
      
      int k, l, I [4];

      /////////////////////////////////////////////////////////////////

      I [1] = i1; I [2] = i2; I [3] = i3;
      
      /////////////////////////////////////////////////////////////////
      
      v1 [0] = x [i1] - x [i2];

      v1 [1] = y [i1] - y [i2];

      v1 [2] = z [i1] - z [i2];

      /////////////////////////////////////////////////////////////////

      v2 [0] = x [i3] - x [i2];

      v2 [1] = y [i3] - y [i2];

      v2 [2] = z [i3] - z [i2];

      /////////////////////////////////////////////////////////////////
      
      v1_norm = v1 [0] * v1 [0] + v1 [1] * v1 [1] + v1 [2] * v1 [2];

      v1_norm = 1.0 / v1_norm;

      /////////////////////////////////////////////////////////////////
      
      v2_norm = v2 [0] * v2 [0] + v2 [1] * v2 [1] + v2 [2] * v2 [2];

      v2_norm = 1.0 / v2_norm;

      /////////////////////////////////////////////////////////////////
      
      v12_norm = sqrt ( v1_norm * v2_norm );
      
      /////////////////////////////////////////////////////////////////
      
      kosinus = ( v1 [0] * v2 [0] + v1 [1] * v2 [1] + v1 [2] * v2 [2] ) * v12_norm;

      sinus =  sqrt ( 1.0 - kosinus * kosinus );

      theta = acos ( kosinus );
      
      /////////////////////////////////////////////////////////////////

      v1_norm = kosinus * v1_norm;

      v2_norm = kosinus * v2_norm;
      
      sinus = - 2.0 * k0 * (theta0 - theta) / sinus;
      
      /////////////////////////////////////////////////////////////////
      
      for (k = 0; k <= 2; k++)
      {
            for (l = 1; l <= 3; l++)
            {
                  tmp1 = Kdelta (l, 1);
                  
                  tmp2 = Kdelta (l, 2);

                  tmp3 = Kdelta (l, 3);

                  /////////////////////////////////////////////////////////////////

                  C1 = tmp1 - tmp2;

                  C2 = tmp3 - tmp2;

                  /////////////////////////////////////////////////////////////////
                  
                  tmp1 = C1 * v2 [k] + C2 * v1 [k];

                  tmp2 = C1 * v1 [k];

                  tmp3 = C2 * v2 [k];

                  /////////////////////////////////////////////////////////////////
                  
                  C1 = tmp1 * v12_norm - tmp2 * v1_norm - tmp3 * v2_norm;
                  
                  /////////////////////////////////////////////////////////////////

                  tmp1 = sinus * C1;
                  
                  if (k == 0) fx [ I [l] ] += tmp1;

                  if (k == 1) fy [ I [l] ] += tmp1;

                  if (k == 2) fz [ I [l] ] += tmp1;
            }
      }

      /////////////////////////////////////////////////////////////////
      
      reduce_forces ( i1, i3, -1.0 );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void stacking_valence_force ( int i1, int i2, int i3, double k0, double theta0 )
{
      double v1 [3], v2 [3], v1_norm, v2_norm, v12_norm, kosinus, sinus, theta, tmp1, tmp2, tmp3, C1, C2;
      
      int k, l, I [4];

      /////////////////////////////////////////////////////////////////

      I [1] = i1; I [2] = i2; I [3] = i3;
      
      /////////////////////////////////////////////////////////////////
      
      v1 [0] = x [i1] - x [i2];

      v1 [1] = y [i1] - y [i2];

      v1 [2] = z [i1] - z [i2];

      /////////////////////////////////////////////////////////////////

      v2 [0] = x [i3] - x [i2];

      v2 [1] = y [i3] - y [i2];

      v2 [2] = z [i3] - z [i2];

      /////////////////////////////////////////////////////////////////
      
      v1_norm = v1 [0] * v1 [0] + v1 [1] * v1 [1] + v1 [2] * v1 [2];

      v1_norm = 1.0 / v1_norm;

      /////////////////////////////////////////////////////////////////
      
      v2_norm = v2 [0] * v2 [0] + v2 [1] * v2 [1] + v2 [2] * v2 [2];

      v2_norm = 1.0 / v2_norm;

      /////////////////////////////////////////////////////////////////
      
      v12_norm = sqrt ( v1_norm * v2_norm );
      
      /////////////////////////////////////////////////////////////////
      
      kosinus = ( v1 [0] * v2 [0] + v1 [1] * v2 [1] + v1 [2] * v2 [2] ) * v12_norm;

      /////////////////////////////////////////////////////////////////

      if ( kosinus < -0.99995 )
      {
            sinus = 0;

            theta = pi;
      }

      else
      {
            sinus =  sqrt ( 1.0 - kosinus * kosinus );

            theta = acos ( kosinus );
      }

      /////////////////////////////////////////////////////////////////
      
      v1_norm = kosinus * v1_norm;

      v2_norm = kosinus * v2_norm;

      /////////////////////////////////////////////////////////////////
      
      if ( sinus ) sinus = - 2.0 * k0 * (theta0 - theta) / sinus;
      
      /////////////////////////////////////////////////////////////////
      
      for (k = 0; k <= 2; k++)
      {
            for (l = 1; l <= 3; l++)
            {
                  tmp1 = Kdelta (l, 1);
                  
                  tmp2 = Kdelta (l, 2);

                  tmp3 = Kdelta (l, 3);

                  /////////////////////////////////////////////////////////////////

                  C1 = tmp1 - tmp2;

                  C2 = tmp3 - tmp2;

                  /////////////////////////////////////////////////////////////////
                  
                  tmp1 = C1 * v2 [k] + C2 * v1 [k];

                  tmp2 = C1 * v1 [k];

                  tmp3 = C2 * v2 [k];

                  /////////////////////////////////////////////////////////////////
                  
                  C1 = tmp1 * v12_norm - tmp2 * v1_norm - tmp3 * v2_norm;
                  
                  /////////////////////////////////////////////////////////////////

                  tmp1 = sinus * C1;
                  
                  if (k == 0) fx [ I [l] ] += tmp1;

                  if (k == 1) fy [ I [l] ] += tmp1;

                  if (k == 2) fz [ I [l] ] += tmp1;
            }
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double dihedral_angle ( int i1, int i2, int i3, int i4 )
{
      double v1 [3], v2 [3], v3 [3], m [3], n [3], kosinus, sinus, psi, tmp1;
      
      /////////////////////////////////////////////////////////////////

      v1 [0] = x [i2] - x [i1];

      v1 [1] = y [i2] - y [i1];

      v1 [2] = z [i2] - z [i1];

      /////////////////////////////////////////////////////////////////

      v2 [0] = x [i3] - x [i2];

      v2 [1] = y [i3] - y [i2];

      v2 [2] = z [i3] - z [i2];

      /////////////////////////////////////////////////////////////////

      v3 [0] = x [i4] - x [i3];

      v3 [1] = y [i4] - y [i3];

      v3 [2] = z [i4] - z [i3];

      /////////////////////////////////////////////////////////////////
      
      cross_product ( v1, v2, m );

      cross_product ( v2, v3, n );
      
      /////////////////////////////////////////////////////////////////
      
      tmp1 = sqrt ( v2 [0] * v2 [0] + v2 [1] * v2 [1] + v2 [2] * v2 [2] );

      kosinus = m [0] * n [0] + m [1] * n [1] + m [2] * n [2];

      sinus = ( v1 [0] * n [0] + v1 [1] * n [1] + v1 [2] * n [2] ) * tmp1;

      psi = atan2 ( sinus, kosinus );

      /////////////////////////////////////////////////////////////////
      
      return ( psi );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void stacking_dihedral_force ( int i1, int i2, int i3, int i4, double k0, double psi0 )
{
      double v1 [3], v2 [3], v3 [3], m [3], n [3], m_norm, n_norm, mn_norm, kosinus, sinus, psi, tmp1, tmp2, tmp3, tmp4;
      
      double A11, A22, A33, A12, A23, A13, C1, C2, D1, D2, D3;
      
      int k, l, I [5];

      /////////////////////////////////////////////////////////////////

      I [1] = i1; I [2] = i2; I [3] = i3; I [4] = i4;
      
      /////////////////////////////////////////////////////////////////
      
      v1 [0] = x [i2] - x [i1];

      v1 [1] = y [i2] - y [i1];

      v1 [2] = z [i2] - z [i1];

      /////////////////////////////////////////////////////////////////

      v2 [0] = x [i3] - x [i2];

      v2 [1] = y [i3] - y [i2];

      v2 [2] = z [i3] - z [i2];

      /////////////////////////////////////////////////////////////////

      v3 [0] = x [i4] - x [i3];

      v3 [1] = y [i4] - y [i3];

      v3 [2] = z [i4] - z [i3];

      /////////////////////////////////////////////////////////////////

      cross_product ( v1, v2, m );

      cross_product ( v2, v3, n );
      
      /////////////////////////////////////////////////////////////////

      m_norm = m [0] * m [0] + m [1] * m [1] + m [2] * m [2];
      
      m_norm = 1.0 / m_norm;

      /////////////////////////////////////////////////////////////////
      
      n_norm = n [0] * n [0] + n [1] * n [1] + n [2] * n [2];

      n_norm = 1.0 / n_norm;

      /////////////////////////////////////////////////////////////////
      
      mn_norm = sqrt ( m_norm * n_norm );

      /////////////////////////////////////////////////////////////////
      
      tmp1 = sqrt ( v2 [0] * v2 [0] + v2 [1] * v2 [1] + v2 [2] * v2 [2] );

      kosinus = ( m [0] * n [0] + m [1] * n [1] + m [2] * n [2] ) * mn_norm;

      sinus =  ( v1 [0] * n [0] + v1 [1] * n [1] + v1 [2] * n [2] ) * tmp1 * mn_norm;

      psi = atan2 ( sinus, kosinus );

      /////////////////////////////////////////////////////////////////

      m_norm = kosinus * m_norm;

      n_norm = kosinus * n_norm;

      sinus = - 2.0 * k0 * (psi0 - psi) / sinus;

      /////////////////////////////////////////////////////////////////
      
      for (k = 0; k <= 2; k++)
      {
            A11 = 0; A22 = 0; A33 = 0; A12 = 0; A23 = 0; A13 = 0;

            for (l = 0; l <= 2; l++)
            {
                  tmp1 = 1 - Kdelta (l, k);
                  
                  A11 += tmp1 * v1 [l] * v1 [l];

                  A22 += tmp1 * v2 [l] * v2 [l];

                  A33 += tmp1 * v3 [l] * v3 [l];

                  A12 += tmp1 * v1 [l] * v2 [l];

                  A23 += tmp1 * v2 [l] * v3 [l];

                  A13 += tmp1 * v1 [l] * v3 [l];
            }

            /////////////////////////////////////////////////////////////////
            
            for (l = 1; l <= 4; l++)
            {
                  tmp1 = Kdelta (l, 1);
                  
                  tmp2 = Kdelta (l, 2);

                  tmp3 = Kdelta (l, 3);

                  tmp4 = Kdelta (l, 4);

                  /////////////////////////////////////////////////////////////////

                  C1 = A22 * (tmp3 - tmp4) + A23 * (tmp3 - tmp2);

                  C2 = A12 * (tmp3 - tmp2) + A22 * (tmp1 - tmp2);

                  D1 = A12 * (tmp4 - tmp3) + A23 * (tmp2 - tmp1) + 2 * A13 * (tmp2 - tmp3);

                  D2 = A11 * (tmp2 - tmp3) + A12 * (tmp2 - tmp1);

                  D3 = A33 * (tmp2 - tmp3) + A23 * (tmp4 - tmp3);

                  /////////////////////////////////////////////////////////////////
                  
                  tmp1 = C1 * v1 [k] + D1 * v2 [k] + C2 * v3 [k];

                  tmp2 = C2 * v1 [k] + D2 * v2 [k];

                  tmp3 = C1 * v3 [k] + D3 * v2 [k];

                  /////////////////////////////////////////////////////////////////
                  
                  C1 = tmp1 * mn_norm + tmp2 * m_norm + tmp3 * n_norm;
                  
                  /////////////////////////////////////////////////////////////////

                  tmp1 = sinus * C1;

                  /////////////////////////////////////////////////////////////////

                  if (k == 0) fx [ I [l] ] += tmp1;

                  if (k == 1) fy [ I [l] ] += tmp1;

                  if (k == 2) fz [ I [l] ] += tmp1;
            }
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////LIST ROUTINES/////////////////////////////////////////////////////

void separate_in_cells ( void )
{
      int i, m;
      
      int i1, i2, i3;

      double tmp, lX, lY, lZ;

      /////////////////////////////////////////////////////////////
      
      if ( McT )
      {
            for (m = 1; m < McT + 1; m++) free ( cell_content [m] );

            free ( cell_content ); free ( cell_mass );

            cell_content = NULL; cell_mass = NULL;

            McX = 0; McY = 0; McZ = 0; McT = 0;
      }
      
      /////////////////////////////////////////////////////////////

      McX = (int) ( side_X / R1_LIST [1] ); lX = side_X / McX;
            
      McY = (int) ( side_Y / R1_LIST [1] ); lY = side_Y / McY;
            
      McZ = (int) ( side_Z / R1_LIST [1] ); lZ = side_Z / McZ;
      
      McT = McX * McY * McZ;
      
      /////////////////////////////////////////////////////////////
      
      cell_content = (int **) calloc ( McT + 1, sizeof(int *) );
      
      for (m = 0; m < McT + 1; m++) cell_content [m] = NULL;

      cell_mass = (int *) calloc ( McT + 1, sizeof(int) );
      
      /////////////////////////////////////////////////////////////
      
      for (i = 1; i < NT_atom + 1; i++)
      {
            tmp = x [i];
                  
            if ( tmp < 0 ) tmp += side_X;
            
            else if ( tmp > side_X ) tmp -= side_X;
                  
            i1 = (int) ceil (tmp / lX);

            /////////////////////////////////////////////////////////////
                  
            tmp = y [i]; 
            
            if ( tmp < 0 ) tmp += side_Y;
            
            else if ( tmp > side_Y ) tmp -= side_Y;   
            
            i2 = (int) ceil (tmp / lY);

            /////////////////////////////////////////////////////////////
                  
            tmp = z [i]; 
            
            if ( tmp < 0 ) tmp += side_Z;
            
            else if ( tmp > side_Z ) tmp -= side_Z;   
                  
            i3 = (int) ceil (tmp / lZ);

            /////////////////////////////////////////////////////////////
                  
            m = (i1 - 1) * McY * McZ + (i2 - 1) * McZ + i3;
                  
            cell_mass [m] += 1;
                  
            cell_content [m] = (int *) realloc ( cell_content [m], (cell_mass [m] + 1) * sizeof(int) );

            cell_content [m] [ cell_mass [m] ] = i;
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void populate_lists ( void )
{
      int J1, J2, k;

      int i1, i2;
      
      int j1, j2;
      
      int m1, m2, m3;

      int n1, n2, n3;
      
      double r2, xx [3];
      
      int index [13][3] = { { -1, -1, -1 }, { -1, -1, 0 },  { -1, -1, 1 },  { -1, 0, -1 }, { -1, 0, 0 }, { -1, 0, 1 }, 
      
      { -1, 1, -1 }, { -1, 1, 0 }, { -1, 1, 1 }, { 0, -1, -1 }, { 0, -1, 0 }, { 0, -1, 1 }, { 0, 0, -1 } };  
      
      ///////////////////////////////////////////////////////

      list_mass [0] = 0;

      list_mass [1] = 0;

      ///////////////////////////////////////////////////////

      separate_in_cells ();

      ///////////////////////////////////////////////////////
      
      for (m1 = 1; m1 < McX + 1; m1++)
      {
            for (m2 = 1; m2 < McY + 1; m2++)
            {
                  for (m3 = 1; m3 < McZ + 1; m3++)
                  {
                        J1 = (m1 - 1) * McY * McZ + (m2 - 1) * McZ + m3;
                        
                        for (k = 0; k < 13; k++)
                        {
                              n1 = m1 + index [k][0]; 
            
                              if (n1 > McX) { n1 -= McX; if ( n1 == m1 || n1 == m1 - 1 ) continue; }
            
                              else if (n1 < 1) { n1 += McX; if ( n1 == m1 || n1 == m1 + 1 ) continue; }

                              n2 = m2 + index [k][1]; 
            
                              if (n2 > McY) { n2 -= McY; if ( n2 == m2 || n2 == m2 - 1 ) continue; }
            
                              else if (n2 < 1) { n2 += McY; if ( n2 == m2 || n2 == m2 + 1 ) continue; }
            
                              n3 = m3 + index [k][2]; 
            
                              if (n3 > McZ) { n3 -= McZ; if ( n3 == m3 || n3 == m3 - 1 ) continue; } 
            
                              else if (n3 < 1) { n3 += McZ; if ( n3 == m3 || n3 == m3 + 1 ) continue; }
                              
                              J2 = (n1 - 1) * McY * McZ + (n2 - 1) * McZ + n3;

                              for (j1 = 1; j1 < cell_mass [J1] + 1; j1++)
                              {
                                    i1 = cell_content [J1][j1];

                                    n1 = part_key [i1];
                                    
                                    for (j2 = 1; j2 < cell_mass [J2] + 1; j2++)
                                    {
                                          i2 = cell_content [J2][j2];

                                          n2 = part_key [i2];
                                          
                                          ///////////////////////////////////////////////////////
                                          
                                          CC_vector ( i1, i2, xx );
                        
                                          r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

                                          ///////////////////////////////////////////////////////

                                          if ( r2 <= R2_LIST [0] )
                                          {
                                                list_mass [0] += 2;
                                                
                                                list_content [0][ list_mass [0] - 1 ] = i1; 
                                                
                                                list_content [0][ list_mass [0] ] = i2;
                                          }

                                          ///////////////////////////////////////////////////////

                                          else if ( r2 <= R2_LIST [1] && n1 != 1 && n1 != 2 && n2 != 1 && n2 != 2 )
                                          {
                                                list_mass [1] += 2;
                                                
                                                list_content [1][ list_mass [1] - 1 ] = i1; 
                                                
                                                list_content [1][ list_mass [1] ] = i2;
                                          }
                                    }
                              }
                        }
                              
                        /////////////////////////////////////////
                              
                        for (j1 = 1; j1 < cell_mass [J1] + 1; j1++)
                        {
                              i1 = cell_content [J1][j1];

                              n1 = part_key [i1];
                              
                              for (j2 = j1 + 1; j2 < cell_mass [J1] + 1; j2++)
                              {
                                    i2 = cell_content [J1][j2];

                                    n2 = part_key [i2];

                                    ///////////////////////////////////////////////////////

                                    CC_vector ( i1, i2, xx );
                        
                                    r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

                                    ///////////////////////////////////////////////////////

                                    if ( r2 <= R2_LIST [0] )
                                    {
                                          list_mass [0] += 2;
                                                
                                          list_content [0][ list_mass [0] - 1 ] = i1; 
                                                
                                          list_content [0][ list_mass [0] ] = i2;
                                    }

                                    ///////////////////////////////////////////////////////

                                    else if ( r2 <= R2_LIST [1] && n1 != 1 && n1 != 2 && n2 != 1 && n2 != 2 )
                                    {
                                          list_mass [1] += 2;
                                                
                                          list_content [1][ list_mass [1] - 1 ] = i1; 
                                                
                                          list_content [1][ list_mass [1] ] = i2;
                                    }
                              }
                        }
                  }
            }
      }

      ///////////////////////////////////////////////////////
      
      for (k = 1; k < hb_N + 1; k++) hb_status [k] = -1;
            
      ///////////////////////////////////////////////////////
            
      for (k = 1; k < list_mass [0] + 1; k += 2)
      {
            i1 = list_content [0][k];

            n1 = part_key [i1];

            if ( n1 > 2 ) continue;

            ///////////////////////////////////////////////////////

            i2 = list_content [0][k + 1];

            n2 = part_key [i2];

            if ( n2 > 2 ) continue;

            ///////////////////////////////////////////////////////
                  
            if (i1 < i2) J1 = (i1 - 1) * NTP_atom - (i1 - 1) * i1 / 2 + i2;

            else J1 = (i2 - 1) * NTP_atom - (i2 - 1) * i2 / 2 + i1;
                  
            for (j1 = 1; j1 < HB_PAIR_N [J1] + 1; j1++) hb_status [ HB_PAIR [J1][j1] ] = 0;
      }

      ///////////////////////////////////////////////////////

      for (k = 1; k < HB_NT_atom + 1; k++) HB_EXCESS [k] = HB_ATOM_N [k] - VALENCE [k];

      for (k = 1; k < hb_N + 1; k++)
      {
            if ( hb_status [k] == 0 ) continue;
            
            for (i1 = 1; i1 < ATOM_HB_N [k] + 1; i1++) HB_EXCESS [ ATOM_HB [k][i1] ] -= 1;
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void check_shifts ( void )
{
      int i;
 
      double x_tmp, y_tmp, z_tmp, r2, shift;

      /////////////////////////////////////
      
      shift = 0;
      
      for (i = 1; i < NT_atom + 1; i++)
      {
            x_tmp = x [i] - x_old [i]; 
            
            y_tmp = y [i] - y_old [i]; 
            
            z_tmp = z [i] - z_old [i];

            r2 = x_tmp * x_tmp + y_tmp * y_tmp + z_tmp * z_tmp;

            if ( r2 > shift ) shift = r2;
      }
      
      shift = sqrt ( shift );

      /////////////////////////////////////
      
      PROGRESS += ( shift + shift );

      if ( PROGRESS > dr1 )
      {
            populate_lists ();

            PROGRESS = 0;
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////FORCE SUMMATION ROUTINES////////////////////////////////////////////////

void SR_interactions_lists ( void )
{
      int i, j, k, i1, i2, j1;

      double xx [3], r2, tmp;
      
      ///////////////////////////////////////////////////////
      
      for (k = 1 + 2 * my_rank; k < list_mass [0] + 1; k += 2 * comm_size)
      {
            i = list_content [0][k];

            j = list_content [0][k + 1];

            ///////////////////////////////////////////////////////
                  
            i1 = maxi_key [i];

            i2 = maxi_key [j];

            if (i1 < i2) j1 = (i1 - 1) * na - (i1 - 1) * i1 / 2 + i2;

            else j1 = (i2 - 1) * na - (i2 - 1) * i2 / 2 + i1;

            /////////////////////////////////////////////////////////////////

            xx [0] = x [i] - x [j]; 
                  
            xx [1] = y [i] - y [j]; 
                  
            xx [2] = z [i] - z [j]; 
      
            half_shift (xx);

            r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

            /////////////////////////////////////////////////////////////////

            tmp = LJ_FORCE ( r2, j1 );

            /////////////////////////////////////////////////////////////////

            i1 = part_key [i];

            i2 = part_key [j];

            if ( i1 != 1 && i1 != 2 && i2 != 1 && i2 != 2 ) tmp += COULOMB_TRUNCATED_FORCE ( r2, j1 );
      
            /////////////////////////////////////////////////////////////////
      
            xx [0] = xx [0] * tmp; 
      
            xx [1] = xx [1] * tmp; 
      
            xx [2] = xx [2] * tmp;

            /////////////////////////////////////////////////////////////////
      
            fx [i] += xx [0]; fx [j] -= xx [0]; 
                  
            fy [i] += xx [1]; fy [j] -= xx [1]; 
                  
            fz [i] += xx [2]; fz [j] -= xx [2];
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void LR_interactions_lists ( void )
{
      int i, j, k, i1, i2, j1;

      double xx [3], r2, tmp;
      
      ///////////////////////////////////////////////////////
      
      for (k = 1 + 2 * my_rank; k < list_mass [1] + 1; k += 2 * comm_size)
      {
            i = list_content [1][k];

            j = list_content [1][k + 1];

            ///////////////////////////////////////////////////////
                  
            i1 = maxi_key [i];

            i2 = maxi_key [j];

            if (i1 < i2) j1 = (i1 - 1) * na - (i1 - 1) * i1 / 2 + i2;

            else j1 = (i2 - 1) * na - (i2 - 1) * i2 / 2 + i1;

            /////////////////////////////////////////////////////////////////

            xx [0] = x [i] - x [j]; 
                  
            xx [1] = y [i] - y [j]; 
                  
            xx [2] = z [i] - z [j]; 
      
            half_shift (xx);

            r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

            /////////////////////////////////////////////////////////////////

            tmp = COULOMB_TRUNCATED_FORCE ( r2, j1 );
      
            /////////////////////////////////////////////////////////////////
      
            xx [0] = xx [0] * tmp; 
      
            xx [1] = xx [1] * tmp; 
      
            xx [2] = xx [2] * tmp;

            /////////////////////////////////////////////////////////////////
      
            flx [i] += xx [0]; flx [j] -= xx [0]; 
                  
            fly [i] += xx [1]; fly [j] -= xx [1]; 
                  
            flz [i] += xx [2]; flz [j] -= xx [2];
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void stacking_interactions ( void )
{
      int i, j, k, k_10, k00, k01, k10, k11, k12, k20, k21, k22, k30;
      double r, theta1, theta2, psi, psi1, psi2;
      double psi0, psi10, psi20, r2;
      
      ///////////////////////////////////////////////////////

      for (k = 1; k < st_N + 1; k++) 
      {
            st_energy [k] = 0;
            st_status [k] = 0;
      }

      ///////////////////////////////////////////////////////

      for (k = 1; k < s3_N + 1; k++)
      {
            i = st_i [k];
            j = st_j [k];

            /////////////////////////////////////////////////////////////////
            k10 = atom_key [i][0];
            k11 = atom_key [i][1];
            k12 = atom_key [i + 1][0];
            
            /////////////////////////////////////////////////////////////////
            k20 = atom_key [j][0];
            k21 = atom_key [j][1];
            k22 = atom_key [j + 1][0];
            
            /////////////////////////////////////////////////////////////////
            r = CC_distance ( k11, k21 );

            ///////////////////////////////////////////////////////
            r -= st_r [k];

            if ( r < 10.0 )
            {
                  i = ATOM_ST [k][0];

                  if ( ST_EXCESS [i] == 1 )
                  {
                        k00 = s3_N + (int) ceil (0.5 * i) - 1;
                        st_status [ k00 ] = 1;
                  }

                  ///////////////////////////////////////////////////////

                  i = ATOM_ST [k][1];

                  if ( ST_EXCESS [i] == 1 )
                  {
                        k00 = s3_N + (int) ceil (0.5 * i) - 1;
                        st_status [ k00 ] = 1;
                  }
            }
            
            ///////////////////////////////////////////////////////

            if ( (k % comm_size) == my_rank )
            {
                  theta1 = valence_angle ( k10, k11, k21 );
                  theta2 = valence_angle ( k20, k21, k11 );
                  psi  = dihedral_angle ( k10, k11, k21, k20 );
                  psi1 = dihedral_angle ( k21, k11, k10, k12 );
                  psi2 = dihedral_angle ( k11, k21, k20, k22 );
                        
                  ///////////////////////////////////////////////////////
                  if ( fabs ( psi - st_psi [k] ) <= pi ) psi0 = st_psi [k];
                  else if ( psi - st_psi [k] > pi ) psi0 = st_psi [k] + 2.0 * pi;
                  else psi0 = st_psi [k] - 2.0 * pi;

                  ///////////////////////////////////////////////////////
                  if ( fabs ( psi1 - st_psi1 [k] ) <= pi ) psi10 = st_psi1 [k];
                  else if ( psi1 - st_psi1 [k] > pi ) psi10 = st_psi1 [k] + 2.0 * pi;
                  else psi10 = st_psi1 [k] - 2.0 * pi;

                  ///////////////////////////////////////////////////////
                  if ( fabs ( psi2 - st_psi2 [k] ) <= pi ) psi20 = st_psi2 [k];
                  else if ( psi2 - st_psi2 [k] > pi ) psi20 = st_psi2 [k] + 2.0 * pi;
                  else psi20 = st_psi2 [k] - 2.0 * pi;

                  ///////////////////////////////////////////////////////
                  r2 = 1.0 + 5.00 * r * r;
                  r2 += 1.50 * ( theta1 - st_theta1 [k] ) * ( theta1 - st_theta1 [k] );
                  r2 += 1.50 * ( theta2 - st_theta2 [k] ) * ( theta2 - st_theta2 [k] );
                  r2 += 0.15 * ( psi - psi0 ) * ( psi - psi0 );
                  r2 += 0.15 * ( psi1 - psi10 ) * ( psi1 - psi10 );
                  r2 += 0.15 * ( psi2 - psi20 ) * ( psi2 - psi20 );
            
                  ///////////////////////////////////////////////////////
                  st_energy [k] = st_E [k] / r2;

                  /////////////////////////////////////////////////////////////////
                  r2 = - st_energy [k] / r2;

                  /////////////////////////////////////////////////////////////////
                  stacking_bond_force ( k11, k21, 5.00 * r2, st_r [k] );
                  stacking_valence_force ( k10, k11, k21, 1.50 * r2, st_theta1 [k] );
                  stacking_valence_force ( k20, k21, k11, 1.50 * r2, st_theta2 [k] );
                  stacking_dihedral_force ( k10, k11, k21, k20, 0.15 * r2, psi0 );
                  stacking_dihedral_force ( k21, k11, k10, k12, 0.15 * r2, psi10 );
                  stacking_dihedral_force ( k11, k21, k20, k22, 0.15 * r2, psi20 );

                  ///////////////////////////////////////////////////////
                  E3_STACK += st_energy [k];
            }
      }

      ///////////////////////////////////////////////////////
      for (k = s3_N + 1 + my_rank; k < st_N + 1; k += comm_size)
      {
            if ( st_status [k] ) continue;

            ///////////////////////////////////////////////////////
            i = st_i [k];
                  
            ///////////////////////////////////////////////////////
            k_10 = atom_key [i - 1][0]; 
            k00 = atom_key [i][0]; 
            k01 = atom_key [i][1];
            k10 = atom_key [i + 1][0]; 
            k20 = atom_key [i + 2][0]; 
            k21 = atom_key [i + 2][1]; 
            k30 = atom_key [i + 3][0];
            
            ///////////////////////////////////////////////////////
            r = CC_distance ( k01, k21 );
            psi1 = dihedral_angle ( k_10, k00, k10, k20 );
            psi2 = dihedral_angle ( k30, k20, k10, k00 );
            
            ///////////////////////////////////////////////////////
            if ( fabs ( psi1 + 2.58684 ) <= pi ) psi10 = - 2.58684;
            else if ( psi1 + 2.58684 > pi ) psi10 = - 2.58684 + 2.0 * pi;
            else psi10 = - 2.58684 - 2.0 * pi;

            ///////////////////////////////////////////////////////
            if ( fabs ( psi2 - 3.07135 ) <= pi ) psi20 = 3.07135;
            else if ( psi2 - 3.07135 > pi ) psi20 = 3.07135 + 2.0 * pi;
            else psi20 = 3.07135 - 2.0 * pi;

            ///////////////////////////////////////////////////////
            r2 = 1.0 + 1.40 * ( r - st_r [k] ) * ( r - st_r [k] );
            r2 += 4.00 * ( psi1 - psi10 ) * ( psi1 - psi10 );
            r2 += 4.00 * ( psi2 - psi20 ) * ( psi2 - psi20 );

            ///////////////////////////////////////////////////////
            st_energy [k] = st_E [k] / r2;

            /////////////////////////////////////////////////////////////////
            r2 = - st_energy [k] / r2;

            /////////////////////////////////////////////////////////////////
            stacking_bond_force ( k01, k21, 1.40 * r2, st_r [k] );
            stacking_dihedral_force ( k_10, k00, k10, k20, 4.00 * r2, psi10 );
            stacking_dihedral_force ( k30, k20, k10, k00, 4.00 * r2, psi20 );

            /////////////////////////////////////////////////////////////////
            E2_STACK += st_energy [k];
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hydrogen_bonds_interactions ( void )
{
      int i, j, k;
      
      double r, r2, theta1, theta2, psi, psi1, psi2;

      int nA = 0;

      int nK = 0;

      ///////////////////////////////////////////////////////
      
      for (i = 1; i < HB_NT_atom + 1; i++) EXCESS [i] = HB_EXCESS [i];

      ///////////////////////////////////////////////////////

      for (k = 1; k < hb_N + 1; k++) hb_energy [k] = 0;

      ///////////////////////////////////////////////////////
      
      for (k = 1; k < hb_N + 1; k++)
      {
            if ( hb_status [k] == -1 ) continue;

            ///////////////////////////////////////////////////////
            
            r = CC_distance ( hb_k11 [k], hb_k21 [k] );

            if ( fabs ( r - hb_r [k] ) > 2.0 )
            {
                  hb_status [k] = 1;

                  for (j = 1; j < ATOM_HB_N [k] + 1; j++) EXCESS [ ATOM_HB [k][j] ] -= 1;

                  continue;
            }

            ///////////////////////////////////////////////////////
            
            theta1 = valence_angle ( hb_k10 [k], hb_k11 [k], hb_k21 [k] );
            
            theta2 = valence_angle ( hb_k20 [k], hb_k21 [k], hb_k11 [k] );

            psi = dihedral_angle ( hb_k10 [k], hb_k11 [k], hb_k21 [k], hb_k20 [k] );

            psi1 = dihedral_angle ( hb_k21 [k], hb_k11 [k], hb_k10 [k], hb_k12 [k] );
            
            psi2 = dihedral_angle ( hb_k11 [k], hb_k21 [k], hb_k20 [k], hb_k22 [k] );

            ///////////////////////////////////////////////////////

            if ( fabs ( psi - hb_psi [k] ) <= pi ) hb_psi0 [k] = hb_psi [k];
                  
            else if ( psi - hb_psi [k] > pi ) hb_psi0 [k] = hb_psi [k] + 2.0 * pi;

            else hb_psi0 [k] = hb_psi [k] - 2.0 * pi;

            ///////////////////////////////////////////////////////

            if ( fabs ( psi1 - hb_psi1 [k] ) <= pi ) hb_psi10 [k] = hb_psi1 [k];
                  
            else if ( psi1 - hb_psi1 [k] > pi ) hb_psi10 [k] = hb_psi1 [k] + 2.0 * pi;

            else hb_psi10 [k] = hb_psi1 [k] - 2.0 * pi;

            ///////////////////////////////////////////////////////

            if ( fabs ( psi2 - hb_psi2 [k] ) <= pi ) hb_psi20 [k] = hb_psi2 [k];
                  
            else if ( psi2 - hb_psi2 [k] > pi ) hb_psi20 [k] = hb_psi2 [k] + 2.0 * pi;

            else hb_psi20 [k] = hb_psi2 [k] - 2.0 * pi;

            ///////////////////////////////////////////////////////
            
            r2 = 5.00 * ( r - hb_r [k] ) * ( r - hb_r [k] );

            r2 += 1.50 * ( theta1 - hb_theta1 [k] ) * ( theta1 - hb_theta1 [k] );

            r2 += 1.50 * ( theta2 - hb_theta2 [k] ) * ( theta2 - hb_theta2 [k] );

            r2 += 0.15 * ( psi - hb_psi0 [k] ) * ( psi - hb_psi0 [k] );

            r2 += 0.15 * ( psi1 - hb_psi10 [k] ) * ( psi1 - hb_psi10 [k] );

            r2 += 0.15 * ( psi2 - hb_psi20 [k] ) * ( psi2 - hb_psi20 [k] );
            
            ///////////////////////////////////////////////////////
            
            hb_energy [k] = hb_E [k] * exp ( - r2 );
      }

      ///////////////////////////////////////////////////////
      
      do
      {
            for (i = 1; i < HB_NT_atom + 1; i++)
            {
                  if ( EXCESS [i] > 0 )
                  {
                        nA += 1;
                  
                        HB_A [ nA ] = i;
                  }
            }
            
            if (nA == 0) break;
      
            ///////////////////////////////////////////////////////
            
            r = (double) rand () / RAND_MAX;

            j = (int) ceil ( r * nA ); if ( j == nA + 1 ) j = nA; if ( j == 0 ) j = 1;

            i = HB_A [j]; //generates atom "i"

            ///////////////////////////////////////////////////////

            for (j = 1; j < HB_ATOM_N [i] + 1; j++) //generates a sequence of bonds for atom "i" 
            {
                  k = HB_ATOM [i][j];

                  if ( hb_status [k] == 0 )
                  {
                        nK += 1;
                  
                        hb_K [ nK ] = k;
                  }
            }

            ///////////////////////////////////////////////////////

            for (i = 1; i < nK + 1; i++) //randomly swaps the generated sequence
            {
                  r = (double) rand () / RAND_MAX;

                  j = (int) ceil ( r * nK ); if ( j == nK + 1 ) j = nK; if ( j == 0 ) j = 1;

                  k = hb_K [i];

                  hb_K [i] = hb_K [j];

                  hb_K [j] = k;
            }

            ///////////////////////////////////////////////////////

            k = hb_K [1]; //generates bond "k" to be turned off
            
            for (i = 2; i < nK + 1; i++) 
            {
                  j = hb_K [i];
                  
                  r2 = exp ( ( hb_energy [j] - hb_energy [k] ) / T );

                  r = (double) rand () / RAND_MAX;

                  if ( r < r2 ) k = j;
            }

            ///////////////////////////////////////////////////////
            
            hb_status [k] = 1;

            for (j = 1; j < ATOM_HB_N [k] + 1; j++) EXCESS [ ATOM_HB [k][j] ] -= 1;

            nA = 0;

            nK = 0;
      }
      while (1 > 0);
      
      ///////////////////////////////////////////////////////

      for (k = 1; k < hb_N + 1; k++)
      {
            if ( hb_status [k] == -1 ) continue;
            
            if ( hb_status [k] == 1 ) { hb_status [k] = 0; continue; }

            ///////////////////////////////////////////////////////
            
            if ( (k % comm_size) == my_rank )
            {
                  r2 = - hb_energy [k];

                  if ( hb_dode [k] ) E3_HB -= r2; 

                  else E2_HB -= r2;  
            
                  ///////////////////////////////////////////////////////
                  
                  stacking_bond_force ( hb_k11 [k], hb_k21 [k], 5.00 * r2, hb_r [k] );

                  stacking_valence_force ( hb_k10 [k], hb_k11 [k], hb_k21 [k], 1.50 * r2, hb_theta1 [k] );

                  stacking_valence_force ( hb_k20 [k], hb_k21 [k], hb_k11 [k], 1.50 * r2, hb_theta2 [k] );
            
                  stacking_dihedral_force ( hb_k10 [k], hb_k11 [k], hb_k21 [k], hb_k20 [k], 0.15 * r2, hb_psi0 [k] );

                  stacking_dihedral_force ( hb_k21 [k], hb_k11 [k], hb_k10 [k], hb_k12 [k], 0.15 * r2, hb_psi10 [k] );

                  stacking_dihedral_force ( hb_k11 [k], hb_k21 [k], hb_k20 [k], hb_k22 [k], 0.15 * r2, hb_psi20 [k] );
            }
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void polymer_constraints ( void )
{
      int i, j, k, entry;
      
      /////////////////////////////////////////////////////////////////

      for (entry = 1 + 2 * my_rank; entry < NBP_atom + 1; entry += 2 * comm_size)
      {
            i = atom_key [ entry ][0];

            j = atom_key [ entry ][1];

            /////////////////////////////////////////////////////////////////

            switch ( amino_key [ entry ] )
            {
                  case 'A': 

                        bead_bond_force ( i, j, 10.0, 4.8515 ); //r3
                        
                        //E_BOND_B += bead_bond_energy ( i, j, 10.0, 4.8515 );

                        if ( entry < NBP_atom ) 
                        {
                              k = atom_key [ entry + 1 ][0];

                              bead_valence_force ( j, i, k, 5.0, 1.9259 ); //alpha1
                              
                              //E_VALENCE_B += bead_valence_energy ( j, i, k, 5.0, 1.9259 );
                        }
                        
                        if ( entry > 1 ) 
                        {
                              k = atom_key [ entry - 1 ][0];

                              bead_valence_force ( j, i, k, 5.0, 1.7029 ); //alpha2
                              
                              //E_VALENCE_B += bead_valence_energy ( j, i, k, 5.0, 1.7029 );
                        }

                        break;

                  case 'C': 

                        bead_bond_force ( i, j, 10.0, 4.2738 ); //r3
                        
                        //E_BOND_B += bead_bond_energy ( i, j, 10.0, 4.2738 );

                        if ( entry < NBP_atom ) 
                        {
                              k = atom_key [ entry + 1 ][0];

                              bead_valence_force ( j, i, k, 5.0, 1.9655 ); //alpha1
                              
                              //E_VALENCE_B += bead_valence_energy ( j, i, k, 5.0, 1.9655 );
                        }
                        
                        if ( entry > 1 ) 
                        {
                              k = atom_key [ entry - 1 ][0];

                              bead_valence_force ( j, i, k, 5.0, 1.5803 ); //alpha2
                              
                              //E_VALENCE_B += bead_valence_energy ( j, i, k, 5.0, 1.5803 );
                        }

                        break;

                  case 'G': 

                        bead_bond_force ( i, j, 10.0, 4.9659 ); //r3
                        
                        //E_BOND_B += bead_bond_energy ( i, j, 10.0, 4.9659 );

                        if ( entry < NBP_atom ) 
                        {
                              k = atom_key [ entry + 1 ][0];

                              bead_valence_force ( j, i, k, 5.0, 1.9150 ); //alpha1
                              
                              //E_VALENCE_B += bead_valence_energy ( j, i, k, 5.0, 1.9150 );
                        }
                        
                        if ( entry > 1 ) 
                        {
                              k = atom_key [ entry - 1 ][0];

                              bead_valence_force ( j, i, k, 5.0, 1.7690 ); //alpha2
                              
                              //E_VALENCE_B += bead_valence_energy ( j, i, k, 5.0, 1.7690 );
                        }

                        break;

                  case 'U': 

                        bead_bond_force ( i, j, 10.0, 4.2733 ); //r3
                        
                        //E_BOND_B += bead_bond_energy ( i, j, 10.0, 4.2733 );

                        if ( entry < NBP_atom ) 
                        {
                              k = atom_key [ entry + 1 ][0];

                              bead_valence_force ( j, i, k, 5.0, 1.9663 ); //alpha1
                              
                              //E_VALENCE_B += bead_valence_energy ( j, i, k, 5.0, 1.9663 );
                        }
                        
                        if ( entry > 1 ) 
                        {
                              k = atom_key [ entry - 1 ][0];

                              bead_valence_force ( j, i, k, 5.0, 1.5735 ); //alpha2
                              
                              //E_VALENCE_B += bead_valence_energy ( j, i, k, 5.0, 1.5735 );
                        }

                        break;
                  
                  default: 
                        
                        break;
            }
      }
            
      /////////////////////////////////////////////////////////////////

      for (entry = 2 + 2 * my_rank; entry < NBP_atom + 1; entry += 2 * comm_size)
      {
            i = atom_key [ entry ][0];

            j = atom_key [ entry - 1 ][0];

            bead_bond_force ( i, j, 64.0, 3.8157 ); //r1 
            
            //E_BOND_B += bead_bond_energy ( i, j, 64.0, 3.8157 );

            if ( entry > 2 )
            {
                  k = atom_key [ entry - 2 ][0];
            
                  bead_valence_force ( i, j, k, 20.0, 1.4440 ); //beta2 
                  
                  //E_VALENCE_B += bead_valence_energy ( i, j, k, 20.0, 1.4440 );
            }
      }

      /////////////////////////////////////////////////////////////////

      for (entry = 3 + 2 * my_rank; entry < NBP_atom + 1; entry += 2 * comm_size)
      {
            i = atom_key [ entry ][0];

            j = atom_key [ entry - 1 ][0];

            k = atom_key [ entry - 2 ][0];

            bead_bond_force ( i, j, 23.0, 4.6010 ); //r2 
            
            //E_BOND_B += bead_bond_energy ( i, j, 23.0, 4.6010 );

            bead_valence_force ( i, j, k, 20.0, 1.5256 ); //beta1 
            
            //E_VALENCE_B += bead_valence_energy ( i, j, k, 20.0, 1.5256 );
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void deterministic_forces_MPI ( void )
{
      FILE * f1;
      
      int i, j;
      
      ////////////////////////////////////////////
      
      for (i = 1; i < NT_atom + 1; i++) { fx [i] = 0; fy [i] = 0; fz [i] = 0; }

      ////////////////////////////////////////////

      E2_HB = 0; E3_HB = 0; E2_STACK = 0; E3_STACK = 0;

      ////////////////////////////////////////////

      stacking_interactions ();

      hydrogen_bonds_interactions ();

      polymer_constraints ();

      SR_interactions_lists ();

      ////////////////////////////////////////////
      
      if ( (klok % 10) == 1 || PROGRESS == 0 )
      {
            for (i = 1; i < NT_atom + 1; i++) { flx [i] = 0; fly [i] = 0; flz [i] = 0; }
            
            LR_interactions_lists ();
      }

      ////////////////////////////////////////////

      j = 0;
      
      for (i = 1; i < NT_atom + 1; i++)
      {
            fsend [j] = fx [i] + flx [i]; j++;

            fsend [j] = fy [i] + fly [i]; j++;

            fsend [j] = fz [i] + flz [i]; j++;
      }

      ////////////////////////////////////////////

      //MPI_Allreduce ( fsend, freceive, 3 * NT_atom, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

      ////////////////////////////////////////////

      j = 0;

      for (i = 1; i < NT_atom + 1; i++)
      {
            fx [i] = freceive [j]; j++;

            fy [i] = freceive [j]; j++;

            fz [i] = freceive [j]; j++;
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void full_forces ( void )
{
      int i, j, k;

      double r1, r2, tmp1, tmp2;

      /////////////////////////////////////
      
      for (i = 1; i < 3001; i += 2)
      {

gen2: 
            r1 = (double) rand () / RAND_MAX; if ( r1 == 0 ) goto gen2;
      
            r2 = (double) rand () / RAND_MAX;
      
            tmp1 = sqrt ( - 2.0 * log (r1) ); tmp2 = 2.0 * pi * r2;

            maxwell_force [i] = tmp1 * cos (tmp2); maxwell_force [i + 1] = tmp1 * sin (tmp2);
      }

      /////////////////////////////////////
      
      for (i = 1; i < NT_atom + 1; i++)
      {
            j = maxi_key [i];

            k = i % 1000;
            
            fx [i] = fx [i] + maxwell_force [k + 1] * SIGMA_FORCE [j];
            
            fy [i] = fy [i] + maxwell_force [k + 1001] * SIGMA_FORCE [j];
            
            fz [i] = fz [i] + maxwell_force [k + 2001] * SIGMA_FORCE [j];
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////MOVE ROUTINES///////////////////////////////////////////////////////

void adjust_box ( int kh )
{
      int k;

      double dx, dy, dz;
      
      /////////////////////////////////////
      
      dx = side_hX - x [kh];

      dy = side_hY - y [kh];

      dz = side_hZ - z [kh];

      /////////////////////////////////////
            
      for (k = 1; k < NT_atom + 1; k++)
      {
            x [k] += dx; x_old [k] += dx;
                  
            y [k] += dy; y_old [k] += dy;
                  
            z [k] += dz; z_old [k] += dz;
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void move_rigid_units ( void )
{
      int k;
      
      //////////////////////////////////////////

      for (k = 1; k < NT_atom + 1; k++) leap_frog_move_atom (k);

      /////////////////////////////////////////////////////////////////

      adjust_box ( 500 );

      /////////////////////////////////////////////////////////////////
      
      //Periodic boundary for crowders

      for (k = NTP_atom + 1; k < NT_atom + 1; k++)
      {
            if ( x [k] < 0 ) 
            {
                  x [k] += side_X; 
                        
                  x_old [k] += side_X;
            }
            
            ////////////////////////////////////////////////////

            if ( x [k] > side_X ) 
            {
                  x [k] -= side_X; 
                        
                  x_old [k] -= side_X;
                  
            } 

            ////////////////////////////////////////
            
            if ( y [k] < 0 ) 
            {
                  y [k] += side_Y; 
                        
                  y_old [k] += side_Y;
            }

            ////////////////////////////////////////

            if ( y [k] > side_Y ) 
            {
                  y [k] -= side_Y; 
                        
                  y_old [k] -= side_Y;
            }

            ////////////////////////////////////////
                  
            if ( z [k] < 0 ) 
            {
                  z [k] += side_Z; 
                        
                  z_old [k] += side_Z;
            }
            

            ////////////////////////////////////////
            
            if ( z [k] > side_Z ) 
            {
                  z [k] -= side_Z; 
                        
                  z_old [k] -= side_Z;
            }
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////DIAGNOSTICS ROUTINES//////////////////////////////////////////////////

void velocity_rescale ( double * sc ) 
{
      int k, l, nV;

      double T0;

      ////////////////////////////
      
      T0 = 0; nV = 0;
      
      for (k = 1; k < NT_atom + 1; k++)
      {
            nV += 3;
            
            l = maxi_key [k];
                  
            T0 += MASS [l] * ( vx [k] * vx [k] + vy [k] * vy [k] + vz [k] * vz [k] );
      }
      
      sc [0] = sqrt ( nV * T / T0 );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double check_maxwell ( void )
{
      int i, j, k, l, nV, count;

      double vsq;

      ////////////////////////////////

      nV = 0; count = 0; 

      for (i = 1; i < NB_atom + 1; i++)
      {
            for (j = 0; j < RIGID_SET [i]; j++)
            {
                  nV += 3;
                  
                  k = atom_key [i][j];

                  l = maxi_key [k];
                  
                  ////////////////////////////////

                  vsq = MASS [l] * vx [k] * vx [k];

                  if ( vsq > T ) count += 1;

                  ////////////////////////////////
                  
                  vsq = MASS [l] * vy [k] * vy [k];
                  
                  if ( vsq > T ) count += 1;

                  ////////////////////////////////

                  vsq = MASS [l] * vz [k] * vz [k];
                  
                  if ( vsq > T ) count += 1;
            }

            ////////////////////////////////
            
            for (j = RIGID_END [i] + 1; j < NS_atom [i] + 1; j++)
            {
                  nV += 3;
                  
                  k = atom_key [i][j];

                  l = maxi_key [k];
                  
                  ////////////////////////////////

                  vsq = MASS [l] * vx [k] * vx [k];

                  if ( vsq > T ) count += 1;

                  ////////////////////////////////
                  
                  vsq = MASS [l] * vy [k] * vy [k];
                  
                  if ( vsq > T ) count += 1;

                  ////////////////////////////////

                  vsq = MASS [l] * vz [k] * vz [k];
                  
                  if ( vsq > T ) count += 1;
            }
      }
      
      ////////////////////////////////

      vsq = (double) count * 100 / nV;
      
      return ( vsq );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////OBSERVABLES//////////////////////////////////////////////////////

double radius_of_gyration ( int n1, int n2 )
{
      int i;

      double xx [3], r1 = 0;

      double xcm = 0, ycm = 0, zcm = 0;
      
      int N = 0;

      ////////////////////////////////////////////////////
      
      for (i = n1; i < n2 + 1; i++) { N++; xcm += x [i]; ycm += y [i]; zcm += z [i]; }

      xcm /= N; ycm /= N; zcm /= N;
      
      ////////////////////////////////////////////////////

      for (i = n1; i < n2 + 1; i++)
      {
            xx [0] = x [i] - xcm; 
                  
            xx [1] = y [i] - ycm;
            
            xx [2] = z [i] - zcm;

            r1 += xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];
      }

      r1 = sqrt ( r1 / N );
      
      return ( r1 );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////READ, WRITE & INITIALIZE ROUTINES///////////////////////////////////////////

void read_maxi_key ( char * file_maxi_key )
{
      FILE * f1;

      int Jn1, Jn2, Jn12;
      
      ////////////////////////////////////////////////////

      f1 = fopen ( file_maxi_key, "r" );
      
      for (Jn1 = 1; Jn1 < na + 1; Jn1++) 
      {
            fscanf ( f1, "%d %le %le %le %le\n", &Jn2, MASS + Jn1, RLJ + Jn1, ELJ + Jn1, QC + Jn1 );

            VISC [ Jn1 ] = 0.555 * RLJ [ Jn1 ]; //0.01 x Water viscosity;

            K1 [ Jn1 ] = exp ( - VISC [ Jn1 ] / MASS [ Jn1 ] * h );
            
            K2 [ Jn1 ] = ( 1.0 - K1 [ Jn1 ] ) / VISC [ Jn1 ];
            
            SIGMA_FORCE [ Jn1 ] = sqrt ( 2.0 * VISC [ Jn1 ] * T / h );
      }

      fclose ( f1 );

      ////////////////////////////////////////////////////

      for (Jn1 = 1; Jn1 < na + 1; Jn1++)
      {
            for (Jn2 = Jn1; Jn2 < na + 1; Jn2++)
            {
                  Jn12 = (Jn1 - 1) * na - (Jn1 - 1) * Jn1 / 2 + Jn2;

                  ////////////////////////////////////////////////////

                  if ( Jn1 <= 6 && Jn2 <= 6 ) D_LJ [ Jn12 ] = 3.2;

                  else D_LJ [ Jn12 ] = RLJ [ Jn1 ] + RLJ [ Jn2 ];

                  ////////////////////////////////////////////////////
                  
                  D2_LJ [ Jn12 ] = D_LJ [ Jn12 ] * D_LJ [ Jn12 ];

                  D2_LJ_CUTOFF [ Jn12 ] = D2_LJ [ Jn12 ];

                  ////////////////////////////////////////////////////

                  D2_LJ_OVERLAP [ Jn12 ] = ( D_LJ [ Jn12 ] - 0.2 * 1.5852 ) * ( D_LJ [ Jn12 ] - 0.2 * 1.5852 );

                  D2_LJ_MINIMUM [ Jn12 ] = ( D_LJ [ Jn12 ] - 0.4406141861 * 1.5852 ) * ( D_LJ [ Jn12 ] - 0.4406141861 * 1.5852 );

                  ////////////////////////////////////////////////////

                  E_LJ [ Jn12 ] = sqrt ( ELJ [ Jn1 ] * ELJ [ Jn2 ] );

                  F_LJ [ Jn12 ] = 12.0 * E_LJ [ Jn12 ] / ( 1.5852 * 1.5852 );

                  ////////////////////////////////////////////////////
                  
                  Q_C [ Jn12 ] = QC [ Jn1 ] * QC [ Jn2 ];
            }
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void list_crowder_types ( void )
{
      int i, j;
      
      /////////////////////////////////////////////////

      crowder_mass = (int *) calloc ( n_crwd, sizeof(int) );
      
      crowder_content = (int **) calloc ( n_crwd, sizeof(int *) );

      /////////////////////////////////////////////////////////////////
      
      for (i = 0; i < n_crwd; i++)
      {
            switch (i) 
            {
            
            case 0://Mg

                  m_crwd [i] = 1;
                  
                  V_crwd [i] = m_crwd [i] * 4.0 * pi * pow ( RLJ [7], 3 ) / 3.0;
                  
                  /////////////////////////////////////////////////

                  RIGID_SET_crwd [i] = 0;
      
                  RIGID_END_crwd [i] = -1;

                  /////////////////////////////////////////////////

                  part_key_crwd [i] = NULL;

                  part_key_crwd [i] = (int *) calloc ( m_crwd [i], sizeof(int) );

                  part_key_crwd [i][0] = 3;

                  /////////////////////////////////////////////////

                  maxi_key_crwd [i] = NULL;

                  maxi_key_crwd [i] = (int *) calloc ( m_crwd [i], sizeof(int) );

                  maxi_key_crwd [i][0] = 7;

                  /////////////////////////////////////////////////

                  RX_crwd [i] = NULL;

                  RY_crwd [i] = NULL;

                  RZ_crwd [i] = NULL;

                  RX_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

                  RY_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

                  RZ_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );
                  
                  RX_crwd [i][0] = 0;
                  
                  RY_crwd [i][0] = 0;

                  RZ_crwd [i][0] = 0;

                  /////////////////////////////////////////////////
                  
                  break;

                  /////////////////////////////////////////////////

            case 1://Cl

                  m_crwd [i] = 1;
                  
                  V_crwd [i] = m_crwd [i] * 4.0 * pi * pow ( RLJ [8], 3 ) / 3.0;
                  
                  /////////////////////////////////////////////////

                  RIGID_SET_crwd [i] = 0;
      
                  RIGID_END_crwd [i] = -1;

                  /////////////////////////////////////////////////

                  part_key_crwd [i] = NULL;

                  part_key_crwd [i] = (int *) calloc ( m_crwd [i], sizeof(int) );

                  part_key_crwd [i][0] = 4;

                  /////////////////////////////////////////////////

                  maxi_key_crwd [i] = NULL;

                  maxi_key_crwd [i] = (int *) calloc ( m_crwd [i], sizeof(int) );

                  maxi_key_crwd [i][0] = 8;

                  /////////////////////////////////////////////////

                  RX_crwd [i] = NULL;

                  RY_crwd [i] = NULL;

                  RZ_crwd [i] = NULL;

                  RX_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

                  RY_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

                  RZ_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );
                  
                  RX_crwd [i][0] = 0;
                  
                  RY_crwd [i][0] = 0;

                  RZ_crwd [i][0] = 0;

                  /////////////////////////////////////////////////
                  
                  break;

                  /////////////////////////////////////////////////

            case 2://K

                  m_crwd [i] = 1;
                  
                  V_crwd [i] = m_crwd [i] * 4.0 * pi * pow ( RLJ [9], 3 ) / 3.0;
                  
                  /////////////////////////////////////////////////

                  RIGID_SET_crwd [i] = 0;
      
                  RIGID_END_crwd [i] = -1;

                  /////////////////////////////////////////////////

                  part_key_crwd [i] = NULL;

                  part_key_crwd [i] = (int *) calloc ( m_crwd [i], sizeof(int) );

                  part_key_crwd [i][0] = 5;

                  /////////////////////////////////////////////////

                  maxi_key_crwd [i] = NULL;

                  maxi_key_crwd [i] = (int *) calloc ( m_crwd [i], sizeof(int) );

                  maxi_key_crwd [i][0] = 9;

                  /////////////////////////////////////////////////

                  RX_crwd [i] = NULL;

                  RY_crwd [i] = NULL;

                  RZ_crwd [i] = NULL;

                  RX_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

                  RY_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

                  RZ_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );
                  
                  RX_crwd [i][0] = 0;
                  
                  RY_crwd [i][0] = 0;

                  RZ_crwd [i][0] = 0;

                  /////////////////////////////////////////////////
                  
                  break;

                  /////////////////////////////////////////////////
            
            default: 
                  
                  break;
            }
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_PDB ( char * file_PDB )
{
      FILE * f1;
      
      char str [10], junk [200];

      char chain = 'A';

      int i, j, k, entry;
      
      /////////////////////////////////////////////////////////////////

      residue_mass = (int *) calloc ( 5, sizeof(int) );
      
      residue_content = (int **) calloc ( 5, sizeof(int *) );

      /////////////////////////////////////////////////////////////////
      
      f1 = fopen ( file_PDB, "r" );

      /////////////////////////////////////////////////////////////////

      k = 0;

      entry = 0;
      
      do
      {
            entry ++;
            
            fscanf ( f1, "%s", str );

            switch ( entry ) 
            {
            
            case 3:

                  if ( strcmp ( str, "P" ) == 0 ) 
                  {
                        N_amino += 1;

                        amino_key = (char *) realloc ( amino_key, (N_amino + 1) * sizeof(char) );

                        /////////////////////////////////////////////////////////////////

                        amino_key [ N_amino ] = 'P';
                              
                        residue_mass [0] += 1;

                        residue_content [0] = (int *) realloc ( residue_content [0], (residue_mass [0] + 1) * sizeof(int) ); 
                              
                        residue_content [0][ residue_mass [0] ] = N_amino;
                  }

                  /////////////////////////////////////////////////////////////////
                  
                  else if ( strcmp ( str, "O5'" ) == 0 ) 
                  {
                        N_amino += 1;

                        amino_key = (char *) realloc ( amino_key, (N_amino + 1) * sizeof(char) );

                        /////////////////////////////////////////////////////////////////

                        k = 1;
                  }

                  /////////////////////////////////////////////////////////////////
                  
                  break;

            case 4:
                  
                  if (k)
                  {
                        if ( strcmp ( str, "A" ) == 0 ) 
                        {
                              amino_key [ N_amino ] = 'A';

                              residue_mass [1] += 1;

                              residue_content [1] = (int *) realloc ( residue_content [1], (residue_mass [1] + 1) * sizeof(int) ); 
                              
                              residue_content [1][ residue_mass [1] ] = N_amino;
                        }

                        /////////////////////////////////////////////////////////////////
                  
                        if ( strcmp ( str, "G" ) == 0 ) 
                        {
                              amino_key [ N_amino ] = 'G';

                              residue_mass [2] += 1;

                              residue_content [2] = (int *) realloc ( residue_content [2], (residue_mass [2] + 1) * sizeof(int) ); 
                              
                              residue_content [2][ residue_mass [2] ] = N_amino;
                        }

                        /////////////////////////////////////////////////////////////////
                  
                        if ( strcmp ( str, "C" ) == 0 ) 
                        {
                              amino_key [ N_amino ] = 'C';

                              residue_mass [3] += 1;

                              residue_content [3] = (int *) realloc ( residue_content [3], (residue_mass [3] + 1) * sizeof(int) ); 
                              
                              residue_content [3][ residue_mass [3] ] = N_amino;
                        }

                        /////////////////////////////////////////////////////////////////
                  
                        if ( strcmp ( str, "U" ) == 0 ) 
                        {
                              amino_key [ N_amino ] = 'U';

                              residue_mass [4] += 1;

                              residue_content [4] = (int *) realloc ( residue_content [4], (residue_mass [4] + 1) * sizeof(int) ); 
                              
                              residue_content [4][ residue_mass [4] ] = N_amino;
                        }
                        
                        /////////////////////////////////////////////////////////////////

                        k = 0;
                  }

                  break;

            case 9: 

                  fgets ( junk, 200, f1 );
                  
                  entry = 0; 

                  break;

            default: 
                  
                  break;
            }
      }
      while ( !feof(f1) );

      /////////////////////////////////////////////////////////////////

      NB_atom = N_amino;

      crowder_key = (int *) calloc ( NB_atom + 1, sizeof(int) );
      
      NS_atom = (int *) calloc ( NB_atom + 1, sizeof(int) );

      RIGID_SET = (int *) calloc ( NB_atom + 1, sizeof(int) );

      RIGID_END = (int *) calloc ( NB_atom + 1, sizeof(int) );
      
      atom_key = (int **) calloc ( NB_atom + 1, sizeof(int *) );
      
      /////////////////////////////////////////////////////////////////

      rewind ( f1 );
      
      i = 0;

      j = 0;
      
      k = 0;

      entry = 0;
      
      do
      {
            entry ++;
            
            fscanf ( f1, "%s", str );

            switch ( entry ) 
            {
            
            case 3:

                  if ( strcmp ( str, "P" ) == 0 )
                  {
                        j = 0;

                        /////////////////////////////////////////////////////////////////
                        
                        x [ NT_atom ] /= i;

                        y [ NT_atom ] /= i;

                        z [ NT_atom ] /= i;
                        
                        /////////////////////////////////////////////////////////////////
                        
                        k += 1;

                        NS_atom [k] = 0;
                        
                        RIGID_SET [k] = 0;

                        RIGID_END [k] = -1;

                        atom_key [k] = (int *) calloc ( NS_atom [k] + 1, sizeof(int) );
                        
                        /////////////////////////////////////////////////////////////////
                        
                        NT_atom += 1;
                        
                        part_key = (int *) realloc ( part_key, (NT_atom + 1) * sizeof(int) );

                        maxi_key = (int *) realloc ( maxi_key, (NT_atom + 1) * sizeof(int) );
                        
                        x = (double *) realloc ( x, (NT_atom + 1) * sizeof(double) );
                        
                        y = (double *) realloc ( y, (NT_atom + 1) * sizeof(double) );

                        z = (double *) realloc ( z, (NT_atom + 1) * sizeof(double) );
                        
                        /////////////////////////////////////////////////////////////////

                        atom_key [k][0] = NT_atom;
                              
                        part_key [ NT_atom ] = 0;

                        maxi_key [ NT_atom ] = 1;

                        /////////////////////////////////////////////////////////////////
                        
                        i = 1;
                        
                        x [ NT_atom ] = 0;

                        y [ NT_atom ] = 0;

                        z [ NT_atom ] = 0;
                  }

                  /////////////////////////////////////////////////////////////////

                  else if ( strcmp ( str, "O5'" ) == 0 )
                  {
                        if (i)
                        {
                              x [ NT_atom ] /= i;

                              y [ NT_atom ] /= i;

                              z [ NT_atom ] /= i;
                        }
                        
                        /////////////////////////////////////////////////////////////////
                        
                        k += 1;

                        NS_atom [k] = 1;
                        
                        RIGID_SET [k] = 0;

                        RIGID_END [k] = -1;

                        atom_key [k] = (int *) calloc ( NS_atom [k] + 1, sizeof(int) );

                        /////////////////////////////////////////////////////////////////

                        NT_atom += 1;
                        
                        part_key = (int *) realloc ( part_key, (NT_atom + 1) * sizeof(int) );

                        maxi_key = (int *) realloc ( maxi_key, (NT_atom + 1) * sizeof(int) );
                        
                        x = (double *) realloc ( x, (NT_atom + 1) * sizeof(double) );
                        
                        y = (double *) realloc ( y, (NT_atom + 1) * sizeof(double) );

                        z = (double *) realloc ( z, (NT_atom + 1) * sizeof(double) );
                        
                        /////////////////////////////////////////////////////////////////

                        atom_key [k][0] = NT_atom;

                        part_key [ NT_atom ] = 1;

                        maxi_key [ NT_atom ] = 2;

                        /////////////////////////////////////////////////////////////////

                        i = 1;
                        
                        x [ NT_atom ] = 0;

                        y [ NT_atom ] = 0;

                        z [ NT_atom ] = 0;
                  }
                  
                  /////////////////////////////////////////////////////////////////

                  else if ( strcmp ( str, "N9" ) == 0 && amino_key [k] == 'A' || strcmp ( str, "N9" ) == 0 && amino_key [k] == 'G' || strcmp ( str, "N1" ) == 0 && amino_key [k] == 'C' || strcmp ( str, "N1" ) == 0 && amino_key [k] == 'U' )  
                  {
                        x [ NT_atom ] /= i;

                        y [ NT_atom ] /= i;

                        z [ NT_atom ] /= i;

                        /////////////////////////////////////////////////////////////////
                        
                        NT_atom += 1;
                        
                        part_key = (int *) realloc ( part_key, (NT_atom + 1) * sizeof(int) );

                        maxi_key = (int *) realloc ( maxi_key, (NT_atom + 1) * sizeof(int) );
                        
                        x = (double *) realloc ( x, (NT_atom + 1) * sizeof(double) );
                        
                        y = (double *) realloc ( y, (NT_atom + 1) * sizeof(double) );

                        z = (double *) realloc ( z, (NT_atom + 1) * sizeof(double) );
                        
                        /////////////////////////////////////////////////////////////////

                        atom_key [k][1] = NT_atom;
                              
                        part_key [ NT_atom ] = 2;
                              
                        if ( amino_key [k] == 'A' ) maxi_key [ NT_atom ] = 3;

                        else if ( amino_key [k] == 'G' ) maxi_key [ NT_atom ] = 4;

                        else if ( amino_key [k] == 'C' ) maxi_key [ NT_atom ] = 5;

                        else if ( amino_key [k] == 'U' ) maxi_key [ NT_atom ] = 6;

                        /////////////////////////////////////////////////////////////////

                        i = 1;

                        x [ NT_atom ] = 0;

                        y [ NT_atom ] = 0;

                        z [ NT_atom ] = 0;
                  }

                  /////////////////////////////////////////////////////////////////
                        
                  else if ( str [0] != 'H' ) i++;
                  
                  /////////////////////////////////////////////////////////////////

                  else j = 1;

                  /////////////////////////////////////////////////////////////////
                  
                  break;

            case 7: if (!j) x [ NT_atom ] += atof ( str ); break;

            case 8: if (!j) y [ NT_atom ] += atof ( str ); break;
            
            case 9: if (!j) z [ NT_atom ] += atof ( str ); fgets ( junk, 200, f1 ); entry = 0; break;
                        
            default: break;

            }
      }
      while ( !feof(f1) );
      
      fclose (f1);

      /////////////////////////////////////////////////////////////////
      
      x [ NT_atom ] /= i;

      y [ NT_atom ] /= i;

      z [ NT_atom ] /= i;
            
      /////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////
      
      IR1 = (double *) calloc ( NB_atom + 1, sizeof(double) );

      IR2 = (double *) calloc ( NB_atom + 1, sizeof(double) );

      IR3 = (double *) calloc ( NB_atom + 1, sizeof(double) );

      /////////////////////////////////////////////////////////////////

      CMS = (double *) calloc ( 3 * NB_atom + 1, sizeof(double) );
      
      VCMS = (double *) calloc ( 3 * NB_atom + 1, sizeof(double) );

      /////////////////////////////////////////////////////////////////
      
      AXES = (double *) calloc ( 9 * NB_atom + 1, sizeof(double) );
      
      W = (double *) calloc ( 3 * NB_atom + 1, sizeof(double) );

      /////////////////////////////////////////////////////////////////
      
      x_old = (double *) calloc ( NT_atom + 1, sizeof(double) );

      y_old = (double *) calloc ( NT_atom + 1, sizeof(double) );

      z_old = (double *) calloc ( NT_atom + 1, sizeof(double) );

      vx = (double *) calloc ( NT_atom + 1, sizeof(double) );

      vy = (double *) calloc ( NT_atom + 1, sizeof(double) );

      vz = (double *) calloc ( NT_atom + 1, sizeof(double) );

      fx = (double *) calloc ( NT_atom + 1, sizeof(double) );

      fy = (double *) calloc ( NT_atom + 1, sizeof(double) );

      fz = (double *) calloc ( NT_atom + 1, sizeof(double) );

      flx = (double *) calloc ( NT_atom + 1, sizeof(double) );

      fly = (double *) calloc ( NT_atom + 1, sizeof(double) );

      flz = (double *) calloc ( NT_atom + 1, sizeof(double) );

      /////////////////////////////////////////////////////////////////
      
      INDX = (int *) calloc ( NT_atom + 1, sizeof(int) );

      JNDX = (int *) calloc ( NT_atom + 1, sizeof(int) );

      for (i = 1; i < NB_atom + 1; i++)
      {
            for (j = 0; j < NS_atom [i] + 1; j++)
            {
                  k = atom_key [i][j];

                  INDX [k] = i;

                  JNDX [k] = j;
            }
      }
      
      /////////////////////////////////////////////////////////////////

      RMASS = (double *) calloc ( NB_atom + 1, sizeof(double) );

      RVISC = (double *) calloc ( NB_atom + 1, sizeof(double) );

      RXYZ = (double *) calloc ( 3 * NB_atom + 1, sizeof(double) );
      
      RX = (double **) calloc ( NB_atom + 1, sizeof(double *) );

      RY = (double **) calloc ( NB_atom + 1, sizeof(double *) ); 

      RZ = (double **) calloc ( NB_atom + 1, sizeof(double *) );

      for (i = 1; i < NB_atom + 1; i++)
      {
            RX [i] = (double *) calloc ( NS_atom [i] + 1, sizeof(double) );

            RY [i] = (double *) calloc ( NS_atom [i] + 1, sizeof(double) );

            RZ [i] = (double *) calloc ( NS_atom [i] + 1, sizeof(double) );
      }
      
      AV = (double *) calloc ( NB_atom + 1, sizeof(double) );

      BV = (double *) calloc ( NB_atom + 1, sizeof(double) );

      CV = (double *) calloc ( NB_atom + 1, sizeof(double) );

      FV = (double *) calloc ( NB_atom + 1, sizeof(double) );

      GV = (double *) calloc ( NB_atom + 1, sizeof(double) );

      HV = (double *) calloc ( NB_atom + 1, sizeof(double) );

      /////////////////////////////////////////////////////////////////
      
      for (i = 1; i < NB_atom + 1; i++)
            
            for (j = 0; j < NS_atom [i] + 1; j++) 
                  
                  generate_atom_velocity ( i, j );
            
      /////////////////////////////////////////////////////////////////
            
      NBP_atom = NB_atom;

      NTP_atom = NT_atom;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int return_HB_JNDX ( char * ATOM, char NUCLEOTIDE )
{
      int j;

      if ( strcmp ( ATOM, "N1" ) == 0 ) j = 2;

      if ( strcmp ( ATOM, "N2" ) == 0 ) j = 3;

      if ( strcmp ( ATOM, "N3" ) == 0 ) { if ( NUCLEOTIDE == 'G' ) j = 4; else j = 3; }

      if ( strcmp ( ATOM, "N4" ) == 0 ) j = 4;

      if ( strcmp ( ATOM, "N6" ) == 0 ) j = 4;

      if ( strcmp ( ATOM, "N7" ) == 0 ) { if ( NUCLEOTIDE == 'A' ) j = 5; if ( NUCLEOTIDE == 'G' ) j = 6; }
            
      if ( strcmp ( ATOM, "O2" ) == 0 ) j = 2;

      if ( strcmp ( ATOM, "O4" ) == 0 ) j = 4;

      if ( strcmp ( ATOM, "O6" ) == 0 ) j = 5;

      if ( strcmp ( ATOM, "O2'" ) == 0 ) j = 0;

      if ( strcmp ( ATOM, "O4'" ) == 0 ) j = 1;
                  
      if ( strcmp ( ATOM, "OP1" ) == 0 ) j = 0;

      if ( strcmp ( ATOM, "OP2" ) == 0 ) j = 1;

      return (j);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initialize_HB_atoms ( void )
{
      int i, j, k, l, m;
      
      ///////////////////////////////////////////////////////

      HB_NS_atom = (int *) calloc ( NBP_atom + 1, sizeof(int) );
      
      HB_atom_key = (int **) calloc ( NBP_atom + 1, sizeof(int *) );

      for (i = 1; i < NBP_atom + 1; i++)
      {
            switch ( amino_key [i] )
            {
                  case 'A': HB_NS_atom [i] = 5; k = 4; l = 0; m = -1; break;

                  case 'C': HB_NS_atom [i] = 4; k = 4; l = 0; m = 2; break;
            
                  case 'G': HB_NS_atom [i] = 6; k = 3; l = 0; m = 5; break;

                  case 'U': HB_NS_atom [i] = 4; k = -1; l = 0; m = 2; break; //case 'U': HB_NS_atom [i] = 4; k = -1; l = 0; m = -1; break; //original

                  case 'P': HB_NS_atom [i] = 1; k = -1; l = -1; m = -1; break;
                  
                  default: break;
            }

            HB_atom_key [i] = (int *) calloc ( HB_NS_atom [i] + 1, sizeof(int) );

            for (j = 0; j < HB_NS_atom [i] + 1; j++)
            {
                  HB_NT_atom += 1;

                  HB_INDX = (int *) realloc ( HB_INDX, ( HB_NT_atom + 1 ) * sizeof(int) );

                  HB_JNDX = (int *) realloc ( HB_JNDX, ( HB_NT_atom + 1 ) * sizeof(int) );
                              
                  VALENCE = (int *) realloc ( VALENCE, ( HB_NT_atom + 1 ) * sizeof(int) );

                  HB_atom_key [i][j] = HB_NT_atom;

                  HB_INDX [ HB_NT_atom ] = i;
                              
                  HB_JNDX [ HB_NT_atom ] = j;

                  if (j == k || j == l || j == m) VALENCE [ HB_NT_atom ] = 2; else VALENCE [ HB_NT_atom ] = 1;
            }
      }

      HB_EXCESS = (int *) calloc ( HB_NT_atom + 1, sizeof(int) );

      HB_ATOM_N = (int *) calloc ( HB_NT_atom + 1, sizeof(int) );

      HB_ATOM = (int **) calloc ( HB_NT_atom + 1, sizeof(int *) );

      HB_A = (int *) calloc ( HB_NT_atom + 1, sizeof(int) );

      EXCESS = (int *) calloc ( HB_NT_atom + 1, sizeof(int) );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initialize_unprocessed_bonds ( char * file_unprocessed_bonds )
{
      FILE * f1;

      int i, k, k_min;
      
      int k10, k11, k12, k20, k21, k22;

      int bp_N;

      int status, RES1, RES2;

      char A0 [10], A01 [10], A02 [10], A1 [10], A2 [10];
      
      ///////////////////////////////////////////////////////
      fprintf( f_trace, "hb_N= %d\n", hb_N);
      
      f1 = fopen ( file_unprocessed_bonds, "r" );
      
      ///////////////////////////////////////////////////////
      
      bp_N = 0;
      
      do
      {
            bp_N += 1;

            do fscanf ( f1, "%s", A1 ); while ( strcmp ( A1, "TER" ) != 0 );
      }
      while ( !feof(f1) );
      
      /////////////////////////////////////////////////////////////////
      
      rewind ( f1 );

      status = 0;
      
      for (i = 1; i < bp_N + 1; i++)
      {
            k_min = hb_N;

            /////////////////////////////////////////////////////////////////
            
            fscanf ( f1, "%s", A0 );

            /////////////////////////////////////////////////////////////////

            fscanf ( f1, "%s", A01 );

            RES1 = atoi ( A01 );

            if ( RES1 < 842 ) RES1 -= 561;

            else RES1 -= 567;

            /////////////////////////////////////////////////////////////////

            fscanf ( f1, "%s", A02 );

            RES2 = atoi ( A02 );

            if ( RES2 < 842 ) RES2 -= 561;

            else RES2 -= 567;

            /////////////////////////////////////////////////////////////////

            fscanf ( f1, "%s", A1 );
            
            /////////////////////////////////////////////////////////////////

            do
            {
                  fscanf ( f1, "%s", A2 );
                  
                  /////////////////////////////////////////////////////////////////
                  
                  if ( strcmp ( A1, "OP1" ) == 0 || strcmp ( A1, "OP2" ) == 0 )
                  {
                        k10 = atom_key [ 2 * RES1 - 1 ][0];
                        
                        k11 = atom_key [ 2 * RES1 - 2 ][0];

                        k12 = atom_key [ 2 * RES1 ][0];
                  }

                  else if ( strcmp ( A1, "O2'" ) == 0 || strcmp ( A1, "O4'" ) == 0 ) 
                  {
                        k10 = atom_key [ 2 * RES1 ][0];
                        
                        k11 = atom_key [ 2 * RES1 - 1 ][0];

                        k12 = atom_key [ 2 * RES1 + 1 ][0];
                  }

                  else
                  {
                        k10 = atom_key [ 2 * RES1 - 1 ][0];

                        k11 = atom_key [ 2 * RES1 - 1 ][1];

                        k12 = atom_key [ 2 * RES1 ][0];
                  }
                  
                  /////////////////////////////////////////////////////////////////
                  
                  if ( strcmp ( A2, "OP1" ) == 0 || strcmp ( A2, "OP2" ) == 0 )
                  { 
                        k20 = atom_key [ 2 * RES2 - 1 ][0];
                        
                        k21 = atom_key [ 2 * RES2 - 2 ][0];

                        k22 = atom_key [ 2 * RES2 ][0];
                  }

                  else if ( strcmp ( A2, "O2'" ) == 0 || strcmp ( A2, "O4'" ) == 0 ) 
                  {
                        k20 = atom_key [ 2 * RES2 ][0];
                        
                        k21 = atom_key [ 2 * RES2 - 1 ][0];

                        k22 = atom_key [ 2 * RES2 + 1 ][0];
                  }

                  else
                  {
                        k20 = atom_key [ 2 * RES2 - 1 ][0];

                        k21 = atom_key [ 2 * RES2 - 1 ][1];

                        k22 = atom_key [ 2 * RES2 ][0];
                  }
                  
                  /////////////////////////////////////////////////////////////////

                  for (k = k_min + 1; k < hb_N + 1; k++)
                  {
                        if ( k10 == hb_k10 [k] && k11 == hb_k11 [k] && k12 == hb_k12 [k] && k20 == hb_k20 [k] && k21 == hb_k21 [k] && k22 == hb_k22 [k] )
                        {
                              status = 1;
                              
                              break;
                        }
                  }

                  /////////////////////////////////////////////////////////////////
                  
                  if ( status )
                  {
                        hb_E [k] -= hs_D * 0.73753;

                        status = 0;
                  }

                  /////////////////////////////////////////////////////////////////

                  else if ( strcmp ( A0, "CAN51" ) == 0 && strcmp ( A2, "O2'" ) == 0 || strcmp ( A0, "CAN52" ) == 0 && strcmp ( A1, "O2'" ) == 0 )
                  {
                        k = hb_N;
                        
                        //hb_E [k] -= hs_D * 0.73753;
                  }
                  
                  /////////////////////////////////////////////////////////////////
                  
                  else
                  {
                        hb_N++;

                        ///////////////////////////////////////////////////////

                        hb_E = (double *) realloc ( hb_E, ( hb_N + 1 ) * sizeof(double) );

                        ///////////////////////////////////////////////////////
            
                        hb_k10 = (int *) realloc ( hb_k10, ( hb_N + 1 ) * sizeof(int) );

                        hb_k11 = (int *) realloc ( hb_k11, ( hb_N + 1 ) * sizeof(int) );

                        hb_k12 = (int *) realloc ( hb_k12, ( hb_N + 1 ) * sizeof(int) );

                        hb_k20 = (int *) realloc ( hb_k20, ( hb_N + 1 ) * sizeof(int) );

                        hb_k21 = (int *) realloc ( hb_k21, ( hb_N + 1 ) * sizeof(int) );

                        hb_k22 = (int *) realloc ( hb_k22, ( hb_N + 1 ) * sizeof(int) );

                        ///////////////////////////////////////////////////////
            
                        hb_r = (double *) realloc ( hb_r, ( hb_N + 1 ) * sizeof(double) );

                        hb_theta1 = (double *) realloc ( hb_theta1, ( hb_N + 1 ) * sizeof(double) );

                        hb_theta2 = (double *) realloc ( hb_theta2, ( hb_N + 1 ) * sizeof(double) );

                        hb_psi = (double *) realloc ( hb_psi, ( hb_N + 1 ) * sizeof(double) );

                        hb_psi1 = (double *) realloc ( hb_psi1, ( hb_N + 1 ) * sizeof(double) );

                        hb_psi2 = (double *) realloc ( hb_psi2, ( hb_N + 1 ) * sizeof(double) );
            
                        ///////////////////////////////////////////////////////

                        hb_k10 [ hb_N ] = k10;

                        hb_k11 [ hb_N ] = k11;

                        hb_k12 [ hb_N ] = k12;

                        hb_k20 [ hb_N ] = k20;

                        hb_k21 [ hb_N ] = k21;

                        hb_k22 [ hb_N ] = k22;

                        ///////////////////////////////////////////////////////

                        hb_r [ hb_N ] = CC_distance ( k11, k21 ); 
                                                
                        hb_theta1 [ hb_N ] = valence_angle ( k10, k11, k21 ); 
                                                
                        hb_theta2 [ hb_N ] = valence_angle ( k20, k21, k11 ); 
                                                
                        hb_psi [ hb_N ] = dihedral_angle ( k10, k11, k21, k20 ); 
                                                
                        hb_psi1 [ hb_N ] = dihedral_angle ( k21, k11, k10, k12 ); 
                                                
                        hb_psi2 [ hb_N ] = dihedral_angle ( k11, k21, k20, k22 ); 

                        ///////////////////////////////////////////////////////

                        hb_E [ hb_N ] = - hs_D * 0.73753;
                        
                        /////////////////////////////////////////////////////////////////

                        if ( strcmp ( A0, "CAN11" ) == 0 ) { hb_r [ hb_N ] = 5.8815; hb_theta1 [ hb_N ] = 2.7283; hb_theta2 [ hb_N ] = 2.5117; hb_psi [ hb_N ] = 1.2559; hb_psi1 [ hb_N ] = 0.9545; hb_psi2 [ hb_N ] = 1.1747; }
                              
                        if ( strcmp ( A0, "CAN12" ) == 0 ) { hb_r [ hb_N ] = 5.8815; hb_theta1 [ hb_N ] = 2.5117; hb_theta2 [ hb_N ] = 2.7283; hb_psi [ hb_N ] = 1.2559; hb_psi1 [ hb_N ] = 1.1747; hb_psi2 [ hb_N ] = 0.9545; }
                                          
                        if ( strcmp ( A0, "CAN21" ) == 0 ) { hb_r [ hb_N ] = 5.6550; hb_theta1 [ hb_N ] = 2.4837; hb_theta2 [ hb_N ] = 2.8230; hb_psi [ hb_N ] = 1.3902; hb_psi1 [ hb_N ] = 1.2174; hb_psi2 [ hb_N ] = 0.7619; }
                                                
                        if ( strcmp ( A0, "CAN22" ) == 0 ) { hb_r [ hb_N ] = 5.6550; hb_theta1 [ hb_N ] = 2.8230; hb_theta2 [ hb_N ] = 2.4837; hb_psi [ hb_N ] = 1.3902; hb_psi1 [ hb_N ] = 0.7619; hb_psi2 [ hb_N ] = 1.2174; }

                        if ( strcmp ( A0, "CAN31" ) == 0 ) { hb_r [ hb_N ] = 5.9234; hb_theta1 [ hb_N ] = 2.9051; hb_theta2 [ hb_N ] = 2.0544; hb_psi [ hb_N ] = 2.1469; hb_psi1 [ hb_N ] = -0.5366; hb_psi2 [ hb_N ] = 1.4288; }
                                                
                        if ( strcmp ( A0, "CAN32" ) == 0 ) { hb_r [ hb_N ] = 5.9234; hb_theta1 [ hb_N ] = 2.0544; hb_theta2 [ hb_N ] = 2.9051; hb_psi [ hb_N ] = 2.1469; hb_psi1 [ hb_N ] = 1.4288; hb_psi2 [ hb_N ] = -0.5366; }

                        if ( strcmp ( A0, "CAN41" ) == 0 ) { hb_r [ hb_N ] = 5.3737; hb_theta1 [ hb_N ] = 2.8295; hb_theta2 [ hb_N ] = 2.0314; hb_psi [ hb_N ] = 0.6063; hb_psi1 [ hb_N ] = 1.0893; hb_psi2 [ hb_N ] = 1.4835; }
                                                
                        if ( strcmp ( A0, "CAN42" ) == 0 ) { hb_r [ hb_N ] = 5.3737; hb_theta1 [ hb_N ] = 2.0314; hb_theta2 [ hb_N ] = 2.8295; hb_psi [ hb_N ] = 0.6063; hb_psi1 [ hb_N ] = 1.4835; hb_psi2 [ hb_N ] = 1.0893; }

                        if ( strcmp ( A0, "CAN51" ) == 0 ) { hb_r [ hb_N ] = 6.4922; hb_theta1 [ hb_N ] = 2.0579; hb_theta2 [ hb_N ] = 1.3569; hb_psi [ hb_N ] = 3.1147; hb_psi1 [ hb_N ] = -1.4796; hb_psi2 [ hb_N ] = 1.5278; }
                        
                        if ( strcmp ( A0, "CAN52" ) == 0 ) { hb_r [ hb_N ] = 6.4922; hb_theta1 [ hb_N ] = 1.3569; hb_theta2 [ hb_N ] = 2.0579; hb_psi [ hb_N ] = 3.1147; hb_psi1 [ hb_N ] = 1.5278; hb_psi2 [ hb_N ] = -1.4796; }
                        
                        /////////////////////////////////////////////////////////////////

                        ATOM_HB_N = (int *) realloc ( ATOM_HB_N, (hb_N + 1) * sizeof(int) ); ATOM_HB_N [ hb_N ] = 0;

                        ATOM_HB = (int **) realloc ( ATOM_HB, (hb_N + 1) * sizeof(int *) ); ATOM_HB [ hb_N ] = NULL;

                        /////////////////////////////////////////////////////////////////

                        hb_dode = (int *) realloc ( hb_dode, ( hb_N + 1 ) * sizeof(int) );

                        if ( strcmp ( A0, "NON00" ) == 0 || strcmp ( A0, "APLAT" ) == 0 ) hb_dode [ hb_N ] = 1;

                        else hb_dode [ hb_N ] = 0;

                        /////////////////////////////////////////////////////////////////

                        hb_code = (char **) realloc ( hb_code, (hb_N + 1) * sizeof(char *) );

                        hb_code [ hb_N ] = (char *) calloc ( 10, sizeof(char) );

                        strcpy ( hb_code [ hb_N ], A0 );

                        hb_RES1 = (int *) realloc ( hb_RES1, (hb_N + 1) * sizeof(int) ); 
                        
                        hb_RES1 [ hb_N ] = atoi ( A01 );
                        
                        hb_RES2 = (int *) realloc ( hb_RES2, (hb_N + 1) * sizeof(int) );

                        hb_RES2 [ hb_N ] = atoi ( A02 );

                        /////////////////////////////////////////////////////////////////

                        k = hb_N;
                  }
                  
                  /////////////////////////////////////////////////////////////////

                  if ( strcmp ( A1, "OP1" ) == 0 || strcmp ( A1, "OP2" ) == 0 ) k10 = 2 * RES1 - 2; 
                  
                  else k10 = 2 * RES1 - 1;

                  k11 = return_HB_JNDX ( A1, amino_key [ k10 ] );

                  k12 = HB_atom_key [ k10 ][ k11 ];

                  /////////////////////////////////////////////////////////////////

                  if ( strcmp ( A2, "OP1" ) == 0 || strcmp ( A2, "OP2" ) == 0 ) k20 = 2 * RES2 - 2; 
                  
                  else k20 = 2 * RES2 - 1;

                  k21 = return_HB_JNDX ( A2, amino_key [ k20 ] );
                  
                  k22 = HB_atom_key [ k20 ][ k21 ];
                  
                  /////////////////////////////////////////////////////////////////
                  
                  ATOM_HB_N [k] += 2;

                  ATOM_HB [k] = (int *) realloc ( ATOM_HB [k], ( ATOM_HB_N [k] + 1 ) * sizeof(int) );

                  ATOM_HB [k][ ATOM_HB_N [k] - 1 ] = k12;
                  
                  ATOM_HB [k][ ATOM_HB_N [k] ] = k22;

                  /////////////////////////////////////////////////////////////////

                  HB_ATOM_N [ k12 ] += 1;

                  HB_ATOM [ k12 ] = (int *) realloc ( HB_ATOM [ k12 ], ( HB_ATOM_N [ k12 ] + 1 ) * sizeof(int) );

                  HB_ATOM [ k12 ][ HB_ATOM_N [ k12 ] ] = k;

                  /////////////////////////////////////////////////////////////////
                  
                  HB_ATOM_N [ k22 ] += 1;

                  HB_ATOM [ k22 ] = (int *) realloc ( HB_ATOM [ k22 ], ( HB_ATOM_N [ k22 ] + 1 ) * sizeof(int) );
                  
                  HB_ATOM [ k22 ][ HB_ATOM_N [ k22 ] ] = k;

                  /////////////////////////////////////////////////////////////////

                  fscanf ( f1, "%s", A1 );
            }
            while ( strcmp ( A1, "TER" ) != 0 );
      }

      fclose (f1);

      /////////////////////////////////////////////////////////////////
      
      hb_psi0 = (double *) calloc ( hb_N + 1, sizeof(double) );

      hb_psi10 = (double *) calloc ( hb_N + 1, sizeof(double) );

      hb_psi20 = (double *) calloc ( hb_N + 1, sizeof(double) );
      
      hb_energy = (double *) calloc ( hb_N + 1, sizeof(double) );

      hb_status = (int *) calloc ( hb_N + 1, sizeof(int) );
      
      hb_K = (int *) calloc ( hb_N + 1, sizeof(int) );

      /////////////////////////////////////////////////////////////////

      NpTP_atom = NTP_atom * ( NTP_atom + 1 ) / 2;

      HB_PAIR_N = (int *) calloc ( NpTP_atom + 1, sizeof(int) );

      HB_PAIR = (int **) calloc ( NpTP_atom + 1, sizeof(int *) );

      for (i = 1; i < hb_N + 1; i++)
      {
            k11 = hb_k11 [i];

            k21 = hb_k21 [i];

            /////////////////////////////////////////////////////////////////

            if (k11 < k21) k = (k11 - 1) * NTP_atom - (k11 - 1) * k11 / 2 + k21;

            else k = (k21 - 1) * NTP_atom - (k21 - 1) * k21 / 2 + k11;

            /////////////////////////////////////////////////////////////////

            HB_PAIR_N [k] += 1;

            HB_PAIR [k] = (int *) realloc ( HB_PAIR [k], ( HB_PAIR_N [k] + 1 ) * sizeof(int) );

            HB_PAIR [k] [ HB_PAIR_N [k] ] = i;
      }

      fprintf( f_trace, "hb_N= %d\n", hb_N);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void free_hydrogen_bonds ( void )
{
      int i;

      ///////////////////////////////////////////////////////

      for (i = 1; i < hb_N + 1; i++) free ( hb_code [i] );

      free ( hb_code );

      free ( hb_RES1 ); 
                        
      free ( hb_RES2 );

      ///////////////////////////////////////////////////////
      
      free ( ATOM_HB_N );
      
      for (i = 1; i < hb_N + 1; i++) free ( ATOM_HB [i] );

      free ( ATOM_HB );

      ///////////////////////////////////////////////////////

      free ( HB_EXCESS );

      free ( EXCESS );

      free ( HB_A );

      ///////////////////////////////////////////////////////
      
      free ( HB_ATOM_N );

      for (i = 1; i < HB_NT_atom + 1; i++) free ( HB_ATOM [i] );
      
      free ( HB_ATOM );

      ///////////////////////////////////////////////////////

      free ( HB_PAIR_N );

      for (i = 1; i < NpTP_atom + 1; i++) free ( HB_PAIR [i] );
      
      free ( HB_PAIR );

      ///////////////////////////////////////////////////////

      HB_NT_atom = 0;

      free ( HB_NS_atom ); HB_NS_atom = NULL;

      for (i = 1; i < NBP_atom + 1; i++) free ( HB_atom_key [i] );
      
      free ( HB_atom_key ); HB_atom_key = NULL;

      free ( HB_INDX ); HB_INDX = NULL; 

      free ( HB_JNDX ); HB_JNDX = NULL;

      free ( VALENCE ); VALENCE = NULL;

      ///////////////////////////////////////////////////////
      
      hb_N = 0;

      free ( hb_dode ); hb_dode = NULL;
      
      free ( hb_k10 ); hb_k10 = NULL;

      free ( hb_k11 ); hb_k11 = NULL;

      free ( hb_k12 ); hb_k12 = NULL;

      free ( hb_k20 ); hb_k20 = NULL;

      free ( hb_k21 ); hb_k21 = NULL;

      free ( hb_k22 ); hb_k22 = NULL;
      
      free ( hb_r ); hb_r = NULL;
      
      free ( hb_theta1 ); hb_theta1 = NULL;

      free ( hb_theta2 ); hb_theta2 = NULL;

      free ( hb_psi ); hb_psi = NULL;

      free ( hb_psi1 ); hb_psi1 = NULL;

      free ( hb_psi2 ); hb_psi2 = NULL;

      free ( hb_E ); hb_E = NULL;

      free ( hb_psi0 ); hb_psi0 = NULL;

      free ( hb_psi10 ); hb_psi10 = NULL;

      free ( hb_psi20 ); hb_psi20 = NULL;

      free ( hb_energy ); hb_energy = NULL;
      
      free ( hb_status ); hb_status = NULL;

      free ( hb_K ); hb_K = NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initialize_unprocessed_stacks ( char * file_unprocessed_stacks )
{
      FILE * f1;

      int i, j;
      
      int k10, k11, k12, k20, k21, k22;
      
      int RES1, RES2;

      char A1 [10], A2 [10];

      long position, finish;

      ///////////////////////////////////////////////////////

      ST_EXCESS = (int *) calloc ( NBP_atom + 1, sizeof(int) );
      
      ///////////////////////////////////////////////////////

      f1 = fopen ( file_unprocessed_stacks, "r" );
      
      fseek ( f1, 0L, SEEK_END ); 
            
      finish = ftell ( f1 );
      
      fseek ( f1, 0L, SEEK_SET );
      
      do
      {
            st_N += 1;

            st_i = (int *) realloc ( st_i, ( st_N + 1 ) * sizeof(int) );

            st_j = (int *) realloc ( st_j, ( st_N + 1 ) * sizeof(int) );

            st_r = (double *) realloc ( st_r, ( st_N + 1 ) * sizeof(double) );

            st_theta1 = (double *) realloc ( st_theta1, ( st_N + 1 ) * sizeof(double) );

            st_theta2 = (double *) realloc ( st_theta2, ( st_N + 1 ) * sizeof(double) );

            st_psi = (double *) realloc ( st_psi, ( st_N + 1 ) * sizeof(double) );

            st_psi1 = (double *) realloc ( st_psi1, ( st_N + 1 ) * sizeof(double) );

            st_psi2 = (double *) realloc ( st_psi2, ( st_N + 1 ) * sizeof(double) );

            st_E = (double *) realloc ( st_E, ( st_N + 1 ) * sizeof(double) );
            
            /////////////////////////////////////////////////////////////////

            fscanf ( f1, "%s %s\n", A1, A2 );

            /////////////////////////////////////////////////////////////////

            RES1 = atoi ( A1 + 1 );

            if ( RES1 < 842 ) RES1 -= 561;

            else RES1 -= 567;

            /////////////////////////////////////////////////////////////////

            RES2 = atoi ( A2 + 1 );

            if ( RES2 < 842 ) RES2 -= 561;

            else RES2 -= 567;

            /////////////////////////////////////////////////////////////////
                  
            i = 2 * RES1 - 1;
            
            j = 2 * RES2 - 1;
            
            /////////////////////////////////////////////////////////////////
                        
            k10 = atom_key [i][0];

            k11 = atom_key [i][1];

            k12 = atom_key [i + 1][0];
            
            /////////////////////////////////////////////////////////////////
                  
            k20 = atom_key [j][0];

            k21 = atom_key [j][1];

            k22 = atom_key [j + 1][0];

            /////////////////////////////////////////////////////////////////

            st_i [ st_N ] = i;

            st_j [ st_N ] = j;
            
            st_r [ st_N ] = CC_distance ( k11, k21 );

            st_theta1 [ st_N ] = valence_angle ( k10, k11, k21 );
            
            st_theta2 [ st_N ] = valence_angle ( k20, k21, k11 );

            st_psi [ st_N ] = dihedral_angle ( k10, k11, k21, k20 );

            st_psi1 [ st_N ] = dihedral_angle ( k21, k11, k10, k12 );
            
            st_psi2 [ st_N ] = dihedral_angle ( k11, k21, k20, k22 );
            
            /////////////////////////////////////////////////////////////////

            st_E [ st_N ] = - st_D;

            /////////////////////////////////////////////////////////////////
            
            ATOM_ST = (int **) realloc ( ATOM_ST, (st_N + 1) * sizeof(int *) ); 
            
            ATOM_ST [ st_N ] = (int *) calloc ( 2, sizeof(int) );

            /////////////////////////////////////////////////////////////////
            //just to be able to record "-" and "+" sites separately///////// 

            if ( A1 [0] == '-' ) i--; 

            if ( A2 [0] == '-' ) j--;

            /////////////////////////////////////////////////////////////////

            ATOM_ST [ st_N ][ 0 ] = i;

            ATOM_ST [ st_N ][ 1 ] = j;

            /////////////////////////////////////////////////////////////////

            ST_EXCESS [i] += 1;

            ST_EXCESS [j] += 1;

            /////////////////////////////////////////////////////////////////

            position = ftell ( f1 );
      }
      while ( position < finish );

      fclose (f1);
      
      /////////////////////////////////////////////////////////////////

      s3_N = st_N;

      /////////////////////////////////////////////////////////////////

      for (i = 3; i < NBP_atom - 2; i += 2)
      {
            st_N += 1;

            st_i = (int *) realloc ( st_i, ( st_N + 1 ) * sizeof(int) );

            st_r = (double *) realloc ( st_r, ( st_N + 1 ) * sizeof(double) );

            st_E = (double *) realloc ( st_E, ( st_N + 1 ) * sizeof(double) );
            
            /////////////////////////////////////////////////////////////////

            st_i [ st_N ] = i;

            if ( amino_key [i] == 'A' && amino_key [i + 2] == 'A' ) { st_r [ st_N ] = 4.1806530; st_E [ st_N ] = - 5.194 + ss_D / 0.7093838769 - 0.319 * (T - 0.594); }

            if ( amino_key [i] == 'A' && amino_key [i + 2] == 'C' ) { st_r [ st_N ] = 3.8260185; st_E [ st_N ] = - 5.146 + ss_D / 0.7185955600 - 0.319 * (T - 0.594); }

            if ( amino_key [i] == 'A' && amino_key [i + 2] == 'G' ) { st_r [ st_N ] = 4.4255305; st_E [ st_N ] = - 5.977 + ss_D / 0.6968019829 + 5.301 * (T - 0.678); }

            if ( amino_key [i] == 'A' && amino_key [i + 2] == 'U' ) { st_r [ st_N ] = 3.8260185; st_E [ st_N ] = - 5.146 + ss_D / 0.7185955600 - 0.319 * (T - 0.594); }

            if ( amino_key [i] == 'C' && amino_key [i + 2] == 'A' ) { st_r [ st_N ] = 4.7010580; st_E [ st_N ] = - 5.163 + ss_D / 0.6847830171 - 0.319 * (T - 0.594); }

            if ( amino_key [i] == 'C' && amino_key [i + 2] == 'C' ) { st_r [ st_N ] = 4.2500910; st_E [ st_N ] = - 4.873 + ss_D / 0.6991615586 - 1.567 * (T - 0.568); }
            
            if ( amino_key [i] == 'C' && amino_key [i + 2] == 'G' ) { st_r [ st_N ] = 4.9790760; st_E [ st_N ] = - 5.482 + ss_D / 0.6816268897 + 0.774 * (T - 0.627); }

            if ( amino_key [i] == 'C' && amino_key [i + 2] == 'U' ) { st_r [ st_N ] = 4.2273615; st_E [ st_N ] = - 4.873 + ss_D / 0.6832570771 - 1.567 * (T - 0.568); }

            if ( amino_key [i] == 'G' && amino_key [i + 2] == 'A' ) { st_r [ st_N ] = 4.0128560; st_E [ st_N ] = - 5.948 + ss_D / 0.6903176657 + 5.301 * (T - 0.678); }

            if ( amino_key [i] == 'G' && amino_key [i + 2] == 'C' ) { st_r [ st_N ] = 3.6784360; st_E [ st_N ] = - 5.927 + ss_D / 0.7042060343 + 4.370 * (T - 0.682); }
            
            if ( amino_key [i] == 'G' && amino_key [i + 2] == 'G' ) { st_r [ st_N ] = 4.2427250; st_E [ st_N ] = - 6.416 + ss_D / 0.6971421514 + 7.346 * (T - 0.728); }

            if ( amino_key [i] == 'G' && amino_key [i + 2] == 'U' ) { st_r [ st_N ] = 3.6616930; st_E [ st_N ] = - 5.836 + ss_D / 0.6984426543 + 2.924 * (T - 0.672); }

            if ( amino_key [i] == 'U' && amino_key [i + 2] == 'A' ) { st_r [ st_N ] = 4.7010580; st_E [ st_N ] = - 5.163 + ss_D / 0.6847830171 - 0.319 * (T - 0.594); }

            if ( amino_key [i] == 'U' && amino_key [i + 2] == 'C' ) { st_r [ st_N ] = 4.2679180; st_E [ st_N ] = - 4.880 + ss_D / 0.6758595771 - 1.567 * (T - 0.568); }

            if ( amino_key [i] == 'U' && amino_key [i + 2] == 'G' ) { st_r [ st_N ] = 4.9977560; st_E [ st_N ] = - 5.886 + ss_D / 0.7025528229 + 2.924 * (T - 0.672); }

            if ( amino_key [i] == 'U' && amino_key [i + 2] == 'U' ) { st_r [ st_N ] = 4.2453650; st_E [ st_N ] = - 4.267 + ss_D / 0.6686014771 - 3.563 * (T - 0.500); }
            
            /////////////////////////////////////////////////////////////////

            ATOM_ST = (int **) realloc ( ATOM_ST, (st_N + 1) * sizeof(int *) ); 
            ATOM_ST [ st_N ] = (int *) calloc ( 2, sizeof(int) );
            
            ATOM_ST [ st_N ][ 0 ] = i;
            ATOM_ST [ st_N ][ 1 ] = i + 1;
                  
            /////////////////////////////////////////////////////////////////

            ST_EXCESS [i] += 1;

            ST_EXCESS [i + 1] += 1;
      }

      /////////////////////////////////////////////////////////////////

      st_energy = (double *) calloc ( st_N + 1, sizeof(double) );

      st_status = (int *) calloc ( st_N + 1, sizeof(int) );

      /////////////////////////////////////////////////////////////////

      for (i = 1; i < NBP_atom + 1; i++) 
      {
            if ( i == 2 * (583 - 561) - 2 ) ST_EXCESS [i] -= 2;
            
            else if ( i == 2 * (587 - 561) - 2 ) ST_EXCESS [i] -= 2;
            
            else if ( i == 2 * (658 - 561) - 2 ) ST_EXCESS [i] -= 2;

            else if ( i == 2 * (687 - 561) - 2 ) ST_EXCESS [i] -= 2;

            else if ( i == 2 * (700 - 561) - 2 ) ST_EXCESS [i] -= 2;

            else if ( i == 2 * (701 - 561) - 2 ) ST_EXCESS [i] -= 2;
            
            else if ( i == 2 * (721 - 561) - 2 ) ST_EXCESS [i] -= 2;

            else if ( i == 2 * (727 - 561) - 2 ) ST_EXCESS [i] -= 2;
            
            else if ( i == 2 * (760 - 561) - 2 ) ST_EXCESS [i] -= 2;
            
            else if ( i == 2 * (782 - 561) - 2 ) ST_EXCESS [i] -= 2;
            
            else if ( i == 2 * (802 - 561) - 2 ) ST_EXCESS [i] -= 2;
            
            else if ( i == 2 * (814 - 561) - 2 ) ST_EXCESS [i] -= 2;
            
            else if ( i == 2 * (881 - 567) - 2 ) ST_EXCESS [i] -= 2;

            else if ( i == 2 * (892 - 567) - 2 ) ST_EXCESS [i] -= 2;

            else if ( i == 2 * (898 - 567) - 2 ) ST_EXCESS [i] -= 2;

            else ST_EXCESS [i] -= 1;
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void free_stacks ( void )
{
      int k;
      
      /////////////////////////////////////////////////

      for (k = 1; k < st_N + 1; k++) free ( ATOM_ST [k] );

      free ( ATOM_ST );

      ///////////////////////////////////////////////////////

      free ( ST_EXCESS );

      ///////////////////////////////////////////////////////

      s3_N = 0; st_N = 0;

      free ( st_i ); st_i = NULL;

      free ( st_j ); st_j = NULL;

      free ( st_r ); st_r = NULL;
      
      free ( st_theta1 ); st_theta1 = NULL;

      free ( st_theta2 ); st_theta2 = NULL;

      free ( st_psi ); st_psi = NULL;

      free ( st_psi1 ); st_psi1 = NULL;

      free ( st_psi2 ); st_psi2 = NULL;

      free ( st_E ); st_E = NULL;
      
      free ( st_energy ); st_energy = NULL;
      
      free ( st_status ); st_status = NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////ADJUST BOX ROUTINES//////////////////////////////////////////////////

void set_box ( int kh )
{
      int i;
      
      double x_min, y_min, z_min;
      
      double x_max, y_max, z_max;
      
      double A, B, C, F, G, H;

      double a0, a1, a2;
      
      double * ep1 = NULL, * ep2 = NULL, * ep3 = NULL;

      double n1 [3], n2 [3], n3 [3];
      
      /////////////RNA dimensions/////////////
      ////////////////////////////////////////
      
      a0 = 0; a1 = 0; a2 = 0;
            
      for (i = 1; i < NTP_atom + 1; i++) 
      {
            a0 += x [i]; 
                  
            a1 += y [i]; 
                  
            a2 += z [i];
      }

      a0 /= NTP_atom; 
      
      a1 /= NTP_atom; 
      
      a2 /= NTP_atom;
      
      ///////////////////////////////////////////////

      A = 0; B = 0; C = 0; F = 0; G = 0; H = 0; 
      
      for (i = 1; i < NTP_atom + 1; i++) 
      {
            x_min = x [i] - a0;
                  
            y_min = y [i] - a1;
                  
            z_min = z [i] - a2;

            ////////////////////////////////////////
                  
            A += y_min * y_min + z_min * z_min;

            B += x_min * x_min + z_min * z_min;

            C += x_min * x_min + y_min * y_min;

            F += y_min * z_min;

            G += x_min * z_min;

            H += x_min * y_min;
      }

      ///////////////////////////////////////////////

      if ( fabs (F) < 1.0e-6 && fabs (G) < 1.0e-6  )
      {
            a0 = sqrt (B * B - 2.0 * A * B + A * A + 4.0 * H * H);

            /////////////////////////////////
            
            x_max = 0.5 * (B + A + a0);

            n1 [0] = 0.5 * (B - A - a0) / H;

            n1 [1] = 1;

            n1 [2] = 0;

            a1 = sqrt ( n1 [0] * n1 [0] + n1 [1] * n1 [1] + n1 [2] * n1 [2] );

            n1 [0] /= a1;

            n1 [1] /= a1;

            n1 [2] /= a1;

            /////////////////////////////////
            
            y_max = 0.5 * (B + A - a0);

            n2 [0] = 0.5 * (B - A + a0) / H;

            n2 [1] = 1;

            n2 [2] = 0;

            a1 = sqrt ( n2 [0] * n2 [0] + n2 [1] * n2 [1] + n2 [2] * n2 [2] );

            n2 [0] /= a1;

            n2 [1] /= a1;

            n2 [2] /= a1;

            /////////////////////////////////
            
            z_max = C;

            n3 [0] = 0;

            n3 [1] = 0;

            n3 [2] = 1;
            
            /////////////////////////////////

            if ( x_max < y_max ) { x_min = x_max; y_min = y_max; ep1 = n1; ep2 = n2; }
      
            else { x_min = y_max; y_min = x_max; ep1 = n2; ep2 = n1; }

            if ( z_max < x_min ) { z_min = y_min; y_min = x_min; x_min = z_max; ep3 = ep2; ep2 = ep1; ep1 = n3; }

            else if ( z_max < y_min ) { z_min = y_min; y_min = z_max; ep3 = ep2; ep2 = n3; }

            else { z_min = z_max; ep3 = n3; }

            /////////////////////////////////

            if ( ep1 [0] < 0 ) { ep1 [0] = - ep1 [0]; ep1 [1] = - ep1 [1]; ep1 [2] = - ep1 [2]; } 

            if ( ep2 [1] < 0 ) { ep2 [0] = - ep2 [0]; ep2 [1] = - ep2 [1]; ep2 [2] = - ep2 [2]; }
      }

      ///////////////////////////////////////////////
      
      else if ( fabs (F) < 1.0e-6 && fabs (H) < 1.0e-6  )
      {
            a0 = sqrt (C * C - 2.0 * A * C + A * A + 4.0 * G * G);

            /////////////////////////////////
            
            x_max = 0.5 * (C + A + a0);

            n1 [0] = 0.5 * (C - A - a0) / G;

            n1 [1] = 0;

            n1 [2] = 1;

            a1 = sqrt ( n1 [0] * n1 [0] + n1 [1] * n1 [1] + n1 [2] * n1 [2] );

            n1 [0] /= a1;

            n1 [1] /= a1;

            n1 [2] /= a1;

            /////////////////////////////////
            
            y_max = 0.5 * (C + A - a0);

            n2 [0] = 0.5 * (C - A + a0) / G;

            n2 [1] = 0;

            n2 [2] = 1;

            a1 = sqrt ( n2 [0] * n2 [0] + n2 [1] * n2 [1] + n2 [2] * n2 [2] );

            n2 [0] /= a1;

            n2 [1] /= a1;

            n2 [2] /= a1;

            /////////////////////////////////
            
            z_max = B;

            n3 [0] = 0;

            n3 [1] = 1;

            n3 [2] = 0;
            
            /////////////////////////////////

            if ( x_max < y_max ) { x_min = x_max; y_min = y_max; ep1 = n1; ep2 = n2; }
      
            else { x_min = y_max; y_min = x_max; ep1 = n2; ep2 = n1; }

            if ( z_max < x_min ) { z_min = y_min; y_min = x_min; x_min = z_max; ep3 = ep2; ep2 = ep1; ep1 = n3; }

            else if ( z_max < y_min ) { z_min = y_min; y_min = z_max; ep3 = ep2; ep2 = n3; }

            else { z_min = z_max; ep3 = n3; }

            /////////////////////////////////
            
            if ( ep1 [0] < 0 ) { ep1 [0] = - ep1 [0]; ep1 [1] = - ep1 [1]; ep1 [2] = - ep1 [2]; } 

            if ( ep2 [1] < 0 ) { ep2 [0] = - ep2 [0]; ep2 [1] = - ep2 [1]; ep2 [2] = - ep2 [2]; }
      }

      ///////////////////////////////////////////////
      
      else if ( fabs (G) < 1.0e-6 && fabs (H) < 1.0e-6  )
      {
            a0 = sqrt (C * C - 2.0 * B * C + B * B + 4.0 * F * F);

            /////////////////////////////////
            
            x_max = 0.5 * (C + B + a0);

            n1 [0] = 0;

            n1 [1] = 0.5 * (C - B - a0) / F;

            n1 [2] = 1;

            a1 = sqrt ( n1 [0] * n1 [0] + n1 [1] * n1 [1] + n1 [2] * n1 [2] );

            n1 [0] /= a1;

            n1 [1] /= a1;

            n1 [2] /= a1;

            /////////////////////////////////
            
            y_max = 0.5 * (C + B - a0);

            n2 [0] = 0;

            n2 [1] = 0.5 * (C - B + a0) / F;

            n2 [2] = 1;

            a1 = sqrt ( n2 [0] * n2 [0] + n2 [1] * n2 [1] + n2 [2] * n2 [2] );

            n2 [0] /= a1;

            n2 [1] /= a1;

            n2 [2] /= a1;

            /////////////////////////////////

            z_max = A;

            n3 [0] = 1;

            n3 [1] = 0;

            n3 [2] = 0;
            
            /////////////////////////////////

            if ( x_max < y_max ) { x_min = x_max; y_min = y_max; ep1 = n1; ep2 = n2; }
      
            else { x_min = y_max; y_min = x_max; ep1 = n2; ep2 = n1; }

            if ( z_max < x_min ) { z_min = y_min; y_min = x_min; x_min = z_max; ep3 = ep2; ep2 = ep1; ep1 = n3; }

            else if ( z_max < y_min ) { z_min = y_min; y_min = z_max; ep3 = ep2; ep2 = n3; }

            else { z_min = z_max; ep3 = n3; }

            /////////////////////////////////

            if ( ep1 [0] < 0 ) { ep1 [0] = - ep1 [0]; ep1 [1] = - ep1 [1]; ep1 [2] = - ep1 [2]; } 

            if ( ep2 [1] < 0 ) { ep2 [0] = - ep2 [0]; ep2 [1] = - ep2 [1]; ep2 [2] = - ep2 [2]; }
      }
      
      ///////////////////////////////////////////////
      
      else if ( fabs (F) < 1.0e-6 && fabs (G) < 1.0e-6 && fabs (H) < 1.0e-6  )
      {
            x_max = A;

            n1 [0] = 1;

            n1 [1] = 0;

            n1 [2] = 0;

            /////////////////////////////////
            
            y_max = B;

            n2 [0] = 0;

            n2 [1] = 1;

            n2 [2] = 0;

            /////////////////////////////////

            z_max = C;

            n3 [0] = 0;

            n3 [1] = 0;

            n3 [2] = 1;
            
            /////////////////////////////////

            if ( x_max < y_max ) { x_min = x_max; y_min = y_max; ep1 = n1; ep2 = n2; }
      
            else { x_min = y_max; y_min = x_max; ep1 = n2; ep2 = n1; }

            if ( z_max < x_min ) { z_min = y_min; y_min = x_min; x_min = z_max; ep3 = ep2; ep2 = ep1; ep1 = n3; }

            else if ( z_max < y_min ) { z_min = y_min; y_min = z_max; ep3 = ep2; ep2 = n3; }

            else { z_min = z_max; ep3 = n3; }

            /////////////////////////////////

            if ( ep1 [0] < 0 ) { ep1 [0] = - ep1 [0]; ep1 [1] = - ep1 [1]; ep1 [2] = - ep1 [2]; } 

            if ( ep2 [1] < 0 ) { ep2 [0] = - ep2 [0]; ep2 [1] = - ep2 [1]; ep2 [2] = - ep2 [2]; }
      }

      ///////////////////////////////////////////////

      else
      {
            a2 = - A - B - C;

            a1 = A * B + A * C + B * C - H * H - G * G - F * F;

            a0 = A * F * F + B * G * G + C * H * H + 2.0 * H * G * F - A * B * C;

            /////////////////////////////////

            x_min = - a2 / 3.0;

            x_max = ( 9.0 * a1 * a2 - 27 * a0 - 2.0 * a2 * a2 * a2 ) / 54.0;
      
            y_max = x_min * x_min - a1 / 3.0;

            y_min = 2.0 * sqrt ( y_max );
      
            z_min = acos ( x_max * pow ( y_max, - 1.5 ) ) / 3.0;

            a0 = pi / 3.0;
      
            /////////////////////////////////

            x_max = x_min + y_min * cos ( z_min );
      
            y_max = x_min - y_min * cos ( z_min - a0 );

            z_max = x_min - y_min * cos ( z_min + a0 );
      
            /////////////////////////////////

            if ( x_max < y_max ) { x_min = x_max; y_min = y_max; }
      
            else { x_min = y_max; y_min = x_max; }

            if ( z_max < x_min ) { z_min = y_min; y_min = x_min; x_min = z_max; }

            else if ( z_max < y_min ) { z_min = y_min; y_min = z_max; }

            else z_min = z_max;

            /////////////////////////////////
      
            x_max = (A - x_min) * F + G * H;

            y_max = (B - x_min) * G + F * H;
      
            z_max = (C - x_min) * H + G * F;

            a0 = x_max / y_max;

            a1 = x_max / z_max;

            a2 = sqrt ( 1.0 + a0 * a0 + a1 * a1 );

            n1 [0] = 1.0 / a2; 
      
            n1 [1] = a0 / a2; 
      
            n1 [2] = a1 / a2;

            ep1 = n1;
            
            /////////////////////////////////
      
            x_max = (A - y_min) * F + G * H;

            y_max = (B - y_min) * G + F * H;
      
            z_max = (C - y_min) * H + G * F;

            a0 = y_max / x_max;

            a1 = y_max / z_max;

            a2 = sqrt ( 1.0 + a0 * a0 + a1 * a1 );

            n2 [0] = a0 / a2; 
      
            n2 [1] = 1.0 / a2; 
      
            n2 [2] = a1 / a2;

            ep2 = n2;

            /////////////////////////////////

            ep3 = n3;
      }
      
      ///////////////////////////////////////////////
      
      cross_product ( ep1, ep2, ep3 );
      
      /////////////////////////////////   
      /////////////////Set box////////////////
      ////////////////////////////////////////
      
      for (i = 1; i < NTP_atom + 1; i++)
      {
            x_max = x [i] * ep1 [0] + y [i] * ep1 [1] + z [i] * ep1 [2]; 
            
            y_max = x [i] * ep2 [0] + y [i] * ep2 [1] + z [i] * ep2 [2]; 
            
            z_max = x [i] * ep3 [0] + y [i] * ep3 [1] + z [i] * ep3 [2];
            
            x [i] = x_max; 
            
            y [i] = y_max; 
            
            z [i] = z_max;
      }

      /////////////////////////////////////////////////

      //side_X = 250.0;

      //side_Y = 250.0;

      //side_Z = 250.0;

      /////////////////////////////////////////////////

      side_XYZ = side_X * side_Y * side_Z;
      
      /////////////////////////////////////////////////
      
      side_hX = 0.5 * side_X;
      
      side_hY = 0.5 * side_Y;
      
      side_hZ = 0.5 * side_Z;

      /////////////////////////////////////////////////
      
      x_min = side_hX - x [kh];

      y_min = side_hY - y [kh];

      z_min = side_hZ - z [kh];

      /////////////////////////////////////////////////

      for (i = 1; i < NTP_atom + 1; i++)
      {
            x [i] += x_min;

            y [i] += y_min;

            z [i] += z_min;

            /////////////////////////////////////////////////

            x_max = x_old [i] * ep1 [0] + y_old [i] * ep1 [1] + z_old [i] * ep1 [2] + x_min;

            y_max = x_old [i] * ep2 [0] + y_old [i] * ep2 [1] + z_old [i] * ep2 [2] + y_min;
            
            z_max = x_old [i] * ep3 [0] + y_old [i] * ep3 [1] + z_old [i] * ep3 [2] + z_min;
            
            x_old [i] = x_max; 
            
            y_old [i] = y_max; 
            
            z_old [i] = z_max;

            /////////////////////////////////////////////////

            x_max = vx [i] * ep1 [0] + vy [i] * ep1 [1] + vz [i] * ep1 [2];

            y_max = vx [i] * ep2 [0] + vy [i] * ep2 [1] + vz [i] * ep2 [2];
            
            z_max = vx [i] * ep3 [0] + vy [i] * ep3 [1] + vz [i] * ep3 [2];
            
            vx [i] = x_max; 
            
            vy [i] = y_max; 
            
            vz [i] = z_max;
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void add_crowder ( int tp )
{
      int i, j, k;
      
      int i2, j2, k2;
      
      int J1, J2;

      int In1, In2, In12;

      int round;

      double tmp1, tmp2, tmp3, xx [3]; 

      double A, B, C, F, G, H;

      double a0, a1, a2;
      
      double * ep1 = NULL, * ep2 = NULL, * ep3 = NULL;

      double n1 [3], n2 [3], n3 [3];
      
      /////////////////////////////////////////////////////////////////
      
      NB_atom += 1;
      
      i = NB_atom;

      /////////////////////////////////////////////////////////////////

      crowder_mass [tp] += 1;

      crowder_content [tp] = (int *) realloc ( crowder_content [tp], (crowder_mass [tp] + 1) * sizeof(int) ); 
                              
      crowder_content [tp][ crowder_mass [tp] ] = i;

      /////////////////////////////////////////////////////////////////

      J1 = 3 * (i - 1);

      J2 = 9 * (i - 1);

      /////////////////////////////////////////////////////////////////

      crowder_key = (int *) realloc ( crowder_key, (i + 1) * sizeof(int) );

      crowder_key [i] = tp;

      /////////////////////////////////////////////////////////////////

      NS_atom = (int *) realloc ( NS_atom, (i + 1) * sizeof(int) );

      NS_atom [i] = m_crwd [tp] - 1;

      /////////////////////////////////////////////////////////////////
      
      atom_key = (int **) realloc ( atom_key, (i + 1) * sizeof(int *) );
      
      atom_key [i] = NULL;

      atom_key [i] = (int *) calloc ( NS_atom [i] + 1, sizeof(int) );

      /////////////////////////////////////////////////////////////////

      RIGID_SET = (int *) realloc ( RIGID_SET, (i + 1) * sizeof(int) );

      RIGID_END = (int *) realloc ( RIGID_END, (i + 1) * sizeof(int) );

      RIGID_SET [i] = RIGID_SET_crwd [tp];
      
      RIGID_END [i] = RIGID_END_crwd [tp];

      /////////////////////////////////////////////////////////////////
      
      for (j = 0; j < NS_atom [i] + 1; j++)
      {
            NT_atom += 1;
                  
            atom_key [i][j] = NT_atom;

            /////////////////////////////////////////////////////////////////

            part_key = (int *) realloc ( part_key, (NT_atom + 1) * sizeof(int) );

            part_key [ NT_atom ] = part_key_crwd [tp][j];

            /////////////////////////////////////////////////////////////////

            maxi_key = (int *) realloc ( maxi_key, (NT_atom + 1) * sizeof(int) );

            maxi_key [ NT_atom ] = maxi_key_crwd [tp][j];

            /////////////////////////////////////////////////////////////////

            INDX = (int *) realloc ( INDX, (NT_atom + 1) * sizeof(int) );
            
            JNDX = (int *) realloc ( JNDX, (NT_atom + 1) * sizeof(int) );

            INDX [ NT_atom ] = i;

            JNDX [ NT_atom ] = j;
      }
      
      /////////////////////////////////////////////////////////////////
      
      IR1 = (double *) realloc ( IR1, (i + 1) * sizeof(double) );

      IR1 [i] = 0;

      IR2 = (double *) realloc ( IR2, (i + 1) * sizeof(double) );
      
      IR2 [i] = 0;

      IR3 = (double *) realloc ( IR3, (i + 1) * sizeof(double) ); 

      IR3 [i] = 0;
      
      /////////////////////////////////////////////////////////////////
      
      CMS = (double *) realloc ( CMS, (3 * i + 1) * sizeof(double) );

      CMS [J1 + 1] = 0; 
      
      CMS [J1 + 2] = 0; 
      
      CMS [J1 + 3] = 0;

      /////////////////////////////////////////////////////////////////

      VCMS = (double *) realloc ( VCMS, (3 * i + 1) * sizeof(double) );

      VCMS [J1 + 1] = 0; 
      
      VCMS [J1 + 2] = 0; 
      
      VCMS [J1 + 3] = 0;
      
      /////////////////////////////////////////////////////////////////
      
      AXES = (double *) realloc ( AXES, (9 * i + 1) * sizeof(double) );

      AXES [J2 + 1] = 0; 
      
      AXES [J2 + 2] = 0; 
      
      AXES [J2 + 3] = 0;

      AXES [J2 + 4] = 0; 
      
      AXES [J2 + 5] = 0; 
      
      AXES [J2 + 6] = 0;

      AXES [J2 + 7] = 0; 
      
      AXES [J2 + 8] = 0; 
      
      AXES [J2 + 9] = 0;

      /////////////////////////////////////////////////////////////////
      
      W = (double *) realloc ( W, (3 * i + 1) * sizeof(double) );

      W [J1 + 1] = 0; 
      
      W [J1 + 2] = 0; 
      
      W [J1 + 3] = 0;
      
      /////////////////////////////////////////////////////////////////

      x = (double *) realloc ( x, (NT_atom + 1) * sizeof(double) );
                        
      y = (double *) realloc ( y, (NT_atom + 1) * sizeof(double) );

      z = (double *) realloc ( z, (NT_atom + 1) * sizeof(double) );

      x [ NT_atom ] = 0;

      y [ NT_atom ] = 0;

      z [ NT_atom ] = 0;

      /////////////////////////////////////////////////////////////////

      x_old = (double *) realloc ( x_old, (NT_atom + 1) * sizeof(double) );

      y_old = (double *) realloc ( y_old, (NT_atom + 1) * sizeof(double) );

      z_old = (double *) realloc ( z_old, (NT_atom + 1) * sizeof(double) );

      x_old [ NT_atom ] = 0;

      y_old [ NT_atom ] = 0;

      z_old [ NT_atom ] = 0;

      /////////////////////////////////////////////////////////////////

      vx = (double *) realloc ( vx, (NT_atom + 1) * sizeof(double) );
                        
      vy = (double *) realloc ( vy, (NT_atom + 1) * sizeof(double) );

      vz = (double *) realloc ( vz, (NT_atom + 1) * sizeof(double) );

      vx [ NT_atom ] = 0;

      vy [ NT_atom ] = 0;

      vz [ NT_atom ] = 0;

      /////////////////////////////////////////////////////////////////

      fx = (double *) realloc ( fx, (NT_atom + 1) * sizeof(double) );

      fy = (double *) realloc ( fy, (NT_atom + 1) * sizeof(double) );

      fz = (double *) realloc ( fz, (NT_atom + 1) * sizeof(double) );

      fx [ NT_atom ] = 0;

      fy [ NT_atom ] = 0;

      fz [ NT_atom ] = 0;
      
      /////////////////////////////////////////////////////////////////

      flx = (double *) realloc ( flx, (NT_atom + 1) * sizeof(double) );

      fly = (double *) realloc ( fly, (NT_atom + 1) * sizeof(double) );

      flz = (double *) realloc ( flz, (NT_atom + 1) * sizeof(double) );

      flx [ NT_atom ] = 0;

      fly [ NT_atom ] = 0;

      flz [ NT_atom ] = 0;
      
      /////////////////////////////////////////////////////////////////

      RMASS = (double *) realloc ( RMASS, (i + 1) * sizeof(double) ); 

      RMASS [i] = 0;

      /////////////////////////////////////////////////////////////////
      
      RVISC = (double *) realloc ( RVISC, (i + 1) * sizeof(double) ); 
      
      RVISC [i] = 0;

      /////////////////////////////////////////////////////////////////
      
      RXYZ = (double *) realloc ( RXYZ, (3 * i + 1) * sizeof(double) ); 
      
      RXYZ [J1 + 1] = 0; 
      
      RXYZ [J1 + 2] = 0; 
      
      RXYZ [J1 + 3] = 0;

      /////////////////////////////////////////////////////////////////
      
      RX = (double **) realloc ( RX, (i + 1) * sizeof(double *) ); 
      
      RX [i] = NULL;

      RX [i] = (double *) calloc ( NS_atom [i] + 1, sizeof(double) );

      /////////////////////////////////////////////////////////////////
      
      RY = (double **) realloc ( RY, (i + 1) * sizeof(double *) ); 
      
      RY [i] = NULL;

      RY [i] = (double *) calloc ( NS_atom [i] + 1, sizeof(double) );

      /////////////////////////////////////////////////////////////////
      
      RZ = (double **) realloc ( RZ, (i + 1) * sizeof(double *) ); 
      
      RZ [i] = NULL;
      
      RZ [i] = (double *) calloc ( NS_atom [i] + 1, sizeof(double) );

      /////////////////////////////////////////////////////////////////

      AV = (double *) realloc ( AV, (i + 1) * sizeof(double) ); 
      
      AV [i] = 0;

      BV = (double *) realloc ( BV, (i + 1) * sizeof(double) ); 
      
      BV [i] = 0;

      CV = (double *) realloc ( CV, (i + 1) * sizeof(double) ); 
      
      CV [i] = 0;

      FV = (double *) realloc ( FV, (i + 1) * sizeof(double) ); 
      
      FV [i] = 0;

      GV = (double *) realloc ( GV, (i + 1) * sizeof(double) ); 
      
      GV [i] = 0;

      HV = (double *) realloc ( HV, (i + 1) * sizeof(double) ); 
      
      HV [i] = 0;
      
      /////////////////////////////////////////////////////////////////

      round = 0;

gencoords:

      round ++;

      //printf ( "%d\n", round );

      /////////////////////////////////////////////////////////////////
      
      a1 = (double) rand () / RAND_MAX;

      tmp1 = a1 * side_X;
      
      a1 = (double) rand () / RAND_MAX; 
      
      tmp2 = a1 * side_Y;

      a1 = (double) rand () / RAND_MAX; 
      
      tmp3 = a1 * side_Z;
      
      /////////////////////////////////////////////////////////////////

      n1 [0] = rand (); 
            
      a1 = (double) rand () / RAND_MAX; 
            
      if (a1 > 0.5) n1 [0] = - n1 [0];

      /////////////////////////////////////////////////
            
      n1 [1] = rand (); 
            
      a1 = (double) rand () / RAND_MAX; 
            
      if (a1 > 0.5) n1 [1] = - n1 [1];

      /////////////////////////////////////////////////
            
      n1 [2] = rand (); 
            
      a1 = (double) rand () / RAND_MAX; 
            
      if (a1 > 0.5) n1 [2] = - n1 [2];

      /////////////////////////////////////////////////
            
      a1 = sqrt ( n1 [0] * n1 [0] + n1 [1] * n1 [1] + n1 [2] * n1 [2] );

      n1 [0] /= a1;

      n1 [1] /= a1;

      n1 [2] /= a1;

      /////////////////////////////////////////////////
            
      generate_perpendicular_vector ( n1, n2 );

      /////////////////////////////////////////////////
                  
      cross_product ( n1, n2, n3 );
      
      /////////////////////////////////////////////////

      for (j = 0; j < NS_atom [i] + 1; j++)
      {
            k = atom_key [i][j];

            In1 = maxi_key [k];

            /////////////////////////////////////////////////////////////////

            x [k] = tmp1 + n1 [0] * RX_crwd [tp][j] + n2 [0] * RY_crwd [tp][j] + n3 [0] * RZ_crwd [tp][j];

            y [k] = tmp2 + n1 [1] * RX_crwd [tp][j] + n2 [1] * RY_crwd [tp][j] + n3 [1] * RZ_crwd [tp][j];

            z [k] = tmp3 + n1 [2] * RX_crwd [tp][j] + n2 [2] * RY_crwd [tp][j] + n3 [2] * RZ_crwd [tp][j];
            
            /////////////////////////////////////////////////////////////////
            
            for (i2 = 1; i2 < i; i2++)
            {
                  for (j2 = 0; j2 < NS_atom [i2] + 1; j2++)
                  {
                        k2 = atom_key [i2][j2];

                        In2 = maxi_key [k2];

                        if (In1 < In2) In12 = (In1 - 1) * na - (In1 - 1) * In1 / 2 + In2;

                        else In12 = (In2 - 1) * na - (In2 - 1) * In2 / 2 + In1;
                        
                        /////////////////////////////////////////////////////////////////
                        
                        xx [0] = x [k] - x [k2]; 
                        
                        xx [1] = y [k] - y [k2]; 
                        
                        xx [2] = z [k] - z [k2];

                        half_shift ( xx );

                        a2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

                        if ( a2 < D2_LJ_OVERLAP [ In12 ] ) goto gencoords;
                  }
            }
      }
      
      /////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////

      for (j = 0; j < RIGID_SET [i]; j++) generate_atom_velocity ( i, j );

      for (j = RIGID_END [i] + 1; j < NS_atom [i] + 1; j++) generate_atom_velocity ( i, j );

      /////////////////////////////////////////////////////////////////
                  
      if ( RIGID_SET [i] <= RIGID_END [i] )
      {
            for (j = RIGID_SET [i]; j < RIGID_END [i] + 1; j++)
            {
                  k = atom_key [i][j];

                  In1 = maxi_key [k];

                  ///////////////////////////////////////////////////////////////// 
                  
                  RMASS [i] += MASS [ In1 ];

                  RVISC [i] += VISC [ In1 ];

                  ///////////////////////////////////////////////////////////////// 

                  CMS [J1 + 1] += MASS [ In1 ] * x [k];
                  
                  CMS [J1 + 2] += MASS [ In1 ] * y [k];
                  
                  CMS [J1 + 3] += MASS [ In1 ] * z [k];
            }
            
            CMS [J1 + 1] /= RMASS [i]; 
            
            CMS [J1 + 2] /= RMASS [i]; 
            
            CMS [J1 + 3] /= RMASS [i];

            /////////////////////////////////////////////////////////////////
                  
            A = 0; B = 0; C = 0; F = 0; G = 0; H = 0;
                  
            for (j = RIGID_SET [i]; j < RIGID_END [i] + 1; j++)
            {
                  k = atom_key [i][j];

                  In1 = maxi_key [k];

                  /////////////////////////////////////////////////////////////////

                  xx [0] = x [k] - CMS [J1 + 1];
                  
                  xx [1] = y [k] - CMS [J1 + 2];
                  
                  xx [2] = z [k] - CMS [J1 + 3];

                  /////////////////////////////////////////////////////////////////

                  A += MASS [ In1 ] * ( xx [1] * xx [1] + xx [2] * xx [2] );
                  
                  B += MASS [ In1 ] * ( xx [0] * xx [0] + xx [2] * xx [2] );
                  
                  C += MASS [ In1 ] * ( xx [0] * xx [0] + xx [1] * xx [1] );
                  
                  F += MASS [ In1 ] * xx [1] * xx [2];
                  
                  G += MASS [ In1 ] * xx [0] * xx [2];
                  
                  H += MASS [ In1 ] * xx [0] * xx [1];
            }

            /////////////////////////////////////////////////////////////////
            
            if ( fabs (F) < 1.0e-6 && fabs (G) < 1.0e-6  )
            {
                  a0 = sqrt (B * B - 2.0 * A * B + A * A + 4.0 * H * H);

                  /////////////////////////////////
            
                  xx [0] = 0.5 * (B + A + a0);

                  n1 [0] = 0.5 * (B - A - a0) / H;

                  n1 [1] = 1;

                  n1 [2] = 0;

                  a1 = sqrt ( n1 [0] * n1 [0] + n1 [1] * n1 [1] + n1 [2] * n1 [2] );

                  n1 [0] /= a1;

                  n1 [1] /= a1;

                  n1 [2] /= a1;

                  /////////////////////////////////
            
                  xx [1] = 0.5 * (B + A - a0);

                  n2 [0] = 0.5 * (B - A + a0) / H;

                  n2 [1] = 1;

                  n2 [2] = 0;

                  a1 = sqrt ( n2 [0] * n2 [0] + n2 [1] * n2 [1] + n2 [2] * n2 [2] );

                  n2 [0] /= a1;

                  n2 [1] /= a1;

                  n2 [2] /= a1;

                  /////////////////////////////////
            
                  xx [2] = C;

                  n3 [0] = 0;

                  n3 [1] = 0;

                  n3 [2] = 1;
            
                  /////////////////////////////////

                  if ( xx [0] < xx [1] ) { tmp1 = xx [0]; tmp2 = xx [1]; ep1 = n1; ep2 = n2; }
      
                  else { tmp1 = xx [1]; tmp2 = xx [0]; ep1 = n2; ep2 = n1; }

                  if ( xx [2] < tmp1 ) { tmp3 = tmp2; tmp2 = tmp1; tmp1 = xx [2]; ep3 = ep2; ep2 = ep1; ep1 = n3; }

                  else if ( xx [2] < tmp2 ) { tmp3 = tmp2; tmp2 = xx [2]; ep3 = ep2; ep2 = n3; }

                  else { tmp3 = xx [2]; ep3 = n3; }
            }

            ///////////////////////////////////////////////
      
            else if ( fabs (F) < 1.0e-6 && fabs (H) < 1.0e-6  )
            {
                  a0 = sqrt (C * C - 2.0 * A * C + A * A + 4.0 * G * G);

                  /////////////////////////////////
            
                  xx [0] = 0.5 * (C + A + a0);

                  n1 [0] = 0.5 * (C - A - a0) / G;

                  n1 [1] = 0;

                  n1 [2] = 1;

                  a1 = sqrt ( n1 [0] * n1 [0] + n1 [1] * n1 [1] + n1 [2] * n1 [2] );

                  n1 [0] /= a1;

                  n1 [1] /= a1;

                  n1 [2] /= a1;

                  /////////////////////////////////
            
                  xx [1] = 0.5 * (C + A - a0);

                  n2 [0] = 0.5 * (C - A + a0) / G;

                  n2 [1] = 0;

                  n2 [2] = 1;

                  a1 = sqrt ( n2 [0] * n2 [0] + n2 [1] * n2 [1] + n2 [2] * n2 [2] );

                  n2 [0] /= a1;

                  n2 [1] /= a1;

                  n2 [2] /= a1;

                  /////////////////////////////////
            
                  xx [2] = B;

                  n3 [0] = 0;

                  n3 [1] = 1;

                  n3 [2] = 0;
            
                  /////////////////////////////////

                  if ( xx [0] < xx [1] ) { tmp1 = xx [0]; tmp2 = xx [1]; ep1 = n1; ep2 = n2; }
      
                  else { tmp1 = xx [1]; tmp2 = xx [0]; ep1 = n2; ep2 = n1; }

                  if ( xx [2] < tmp1 ) { tmp3 = tmp2; tmp2 = tmp1; tmp1 = xx [2]; ep3 = ep2; ep2 = ep1; ep1 = n3; }

                  else if ( xx [2] < tmp2 ) { tmp3 = tmp2; tmp2 = xx [2]; ep3 = ep2; ep2 = n3; }

                  else { tmp3 = xx [2]; ep3 = n3; }
            }

            ///////////////////////////////////////////////
      
            else if ( fabs (G) < 1.0e-6 && fabs (H) < 1.0e-6  )
            {
                  a0 = sqrt (C * C - 2.0 * B * C + B * B + 4.0 * F * F);

                  /////////////////////////////////
            
                  xx [0] = 0.5 * (C + B + a0);

                  n1 [0] = 0;

                  n1 [1] = 0.5 * (C - B - a0) / F;

                  n1 [2] = 1;

                  a1 = sqrt ( n1 [0] * n1 [0] + n1 [1] * n1 [1] + n1 [2] * n1 [2] );

                  n1 [0] /= a1;

                  n1 [1] /= a1;

                  n1 [2] /= a1;

                  /////////////////////////////////
            
                  xx [1] = 0.5 * (C + B - a0);

                  n2 [0] = 0;

                  n2 [1] = 0.5 * (C - B + a0) / F;

                  n2 [2] = 1;

                  a1 = sqrt ( n2 [0] * n2 [0] + n2 [1] * n2 [1] + n2 [2] * n2 [2] );

                  n2 [0] /= a1;

                  n2 [1] /= a1;

                  n2 [2] /= a1;

                  /////////////////////////////////

                  xx [2] = A;

                  n3 [0] = 1;

                  n3 [1] = 0;

                  n3 [2] = 0;
            
                  /////////////////////////////////

                  if ( xx [0] < xx [1] ) { tmp1 = xx [0]; tmp2 = xx [1]; ep1 = n1; ep2 = n2; }
      
                  else { tmp1 = xx [1]; tmp2 = xx [0]; ep1 = n2; ep2 = n1; }

                  if ( xx [2] < tmp1 ) { tmp3 = tmp2; tmp2 = tmp1; tmp1 = xx [2]; ep3 = ep2; ep2 = ep1; ep1 = n3; }

                  else if ( xx [2] < tmp2 ) { tmp3 = tmp2; tmp2 = xx [2]; ep3 = ep2; ep2 = n3; }

                  else { tmp3 = xx [2]; ep3 = n3; }
            }
      
            ///////////////////////////////////////////////
      
            else if ( fabs (F) < 1.0e-6 && fabs (G) < 1.0e-6 && fabs (H) < 1.0e-6  )
            {
                  xx [0] = A;

                  n1 [0] = 1;

                  n1 [1] = 0;

                  n1 [2] = 0;

                  /////////////////////////////////
            
                  xx [1] = B;

                  n2 [0] = 0;

                  n2 [1] = 1;

                  n2 [2] = 0;

                  /////////////////////////////////

                  xx [2] = C;

                  n3 [0] = 0;

                  n3 [1] = 0;

                  n3 [2] = 1;
            
                  /////////////////////////////////

                  if ( xx [0] < xx [1] ) { tmp1 = xx [0]; tmp2 = xx [1]; ep1 = n1; ep2 = n2; }
      
                  else { tmp1 = xx [1]; tmp2 = xx [0]; ep1 = n2; ep2 = n1; }

                  if ( xx [2] < tmp1 ) { tmp3 = tmp2; tmp2 = tmp1; tmp1 = xx [2]; ep3 = ep2; ep2 = ep1; ep1 = n3; }

                  else if ( xx [2] < tmp2 ) { tmp3 = tmp2; tmp2 = xx [2]; ep3 = ep2; ep2 = n3; }

                  else { tmp3 = xx [2]; ep3 = n3; }
            }

            ///////////////////////////////////////////////

            else
            {
                  a2 = - A - B - C;

                  a1 = A * B + A * C + B * C - H * H - G * G - F * F;

                  a0 = A * F * F + B * G * G + C * H * H + 2.0 * H * G * F - A * B * C;

                  /////////////////////////////////

                  tmp1 = - a2 / 3.0;

                  xx [0] = ( 9.0 * a1 * a2 - 27 * a0 - 2.0 * a2 * a2 * a2 ) / 54.0;
      
                  xx [1] = tmp1 * tmp1 - a1 / 3.0;

                  tmp2 = 2.0 * sqrt ( xx [1] );
      
                  tmp3 = acos ( xx [0] * pow ( xx [1], - 1.5 ) ) / 3.0;

                  a0 = pi / 3.0;
      
                  /////////////////////////////////

                  xx [0] = tmp1 + tmp2 * cos ( tmp3 );
      
                  xx [1] = tmp1 - tmp2 * cos ( tmp3 - a0 );

                  xx [2] = tmp1 - tmp2 * cos ( tmp3 + a0 );
                  
                  /////////////////////////////////

                  if ( xx [0] < xx [1] ) { tmp1 = xx [0]; tmp2 = xx [1]; }
      
                  else { tmp1 = xx [1]; tmp2 = xx [0]; }

                  if ( xx [2] < tmp1 ) { tmp3 = tmp2; tmp2 = tmp1; tmp1 = xx [2]; }

                  else if ( xx [2] < tmp2 ) { tmp3 = tmp2; tmp2 = xx [2]; }

                  else tmp3 = xx [2];

                  /////////////////////////////////
                  
                  xx [0] = (A - tmp1) * F + G * H;

                  xx [1] = (B - tmp1) * G + F * H;
      
                  xx [2] = (C - tmp1) * H + G * F;

                  a0 = xx [0] / xx [1];

                  a1 = xx [0] / xx [2];

                  a2 = sqrt ( 1.0 + a0 * a0 + a1 * a1 );

                  n1 [0] = 1.0 / a2; 
      
                  n1 [1] = a0 / a2; 
      
                  n1 [2] = a1 / a2;

                  ep1 = n1;
            
                  /////////////////////////////////
      
                  xx [0] = (A - tmp2) * F + G * H;

                  xx [1] = (B - tmp2) * G + F * H;
      
                  xx [2] = (C - tmp2) * H + G * F;

                  a0 = xx [1] / xx [0];

                  a1 = xx [1] / xx [2];

                  a2 = sqrt ( 1.0 + a0 * a0 + a1 * a1 );

                  n2 [0] = a0 / a2; 
      
                  n2 [1] = 1.0 / a2; 
      
                  n2 [2] = a1 / a2;

                  ep2 = n2;

                  /////////////////////////////////

                  ep3 = n3;
            }
            
            /////////////////////////////////////////////////////////////////

            IR1 [i] = tmp1;
            
            IR2 [i] = tmp2;
            
            IR3 [i] = tmp3;

            /////////////////////////////////////////////////////////////////
            
            if ( ep1 [0] < 0 ) { ep1 [0] = - ep1 [0]; ep1 [1] = - ep1 [1]; ep1 [2] = - ep1 [2]; } 

            if ( ep2 [1] < 0 ) { ep2 [0] = - ep2 [0]; ep2 [1] = - ep2 [1]; ep2 [2] = - ep2 [2]; }
            
            cross_product ( ep1, ep2, ep3 );

            /////////////////////////////////////////////////////////////////

            AXES [J2 + 1] = ep1 [0]; AXES [J2 + 2] = ep1 [1]; AXES [J2 + 3] = ep1 [2];

            AXES [J2 + 4] = ep2 [0]; AXES [J2 + 5] = ep2 [1]; AXES [J2 + 6] = ep2 [2];

            AXES [J2 + 7] = ep3 [0]; AXES [J2 + 8] = ep3 [1]; AXES [J2 + 9] = ep3 [2];
                  
            ///////////////////////////////////////////////////////////////// 
                  
            for (j = RIGID_SET [i]; j < RIGID_END [i] + 1; j++)
            {
                  k = atom_key [i][j];

                  In1 = maxi_key [k];

                  ///////////////////////////////////////////////////////////////// 
                        
                  xx [0] = x [k] - CMS [J1 + 1];
                  
                  xx [1] = y [k] - CMS [J1 + 2];
                  
                  xx [2] = z [k] - CMS [J1 + 3];

                  ///////////////////////////////////////////////////////////////// 
                  
                  RX [i][j] = xx [0] * ep1 [0] + xx [1] * ep1 [1] + xx [2] * ep1 [2];
            
                  RY [i][j] = xx [0] * ep2 [0] + xx [1] * ep2 [1] + xx [2] * ep2 [2];

                  RZ [i][j] = xx [0] * ep3 [0] + xx [1] * ep3 [1] + xx [2] * ep3 [2];

                  ///////////////////////////////////////////////////////////////// 

                  RXYZ [J1 + 1] += VISC [ In1 ] * RX [i][j];

                  RXYZ [J1 + 2] += VISC [ In1 ] * RY [i][j];

                  RXYZ [J1 + 3] += VISC [ In1 ] * RZ [i][j];
                        
                  ///////////////////////////////////////////////////////////////// 

                  AV [i] += VISC [ In1 ]  * ( RY [i][j] * RY [i][j] + RZ [i][j] * RZ [i][j] );

                  BV [i] += VISC [ In1 ] * ( RX [i][j] * RX [i][j] + RZ [i][j] * RZ [i][j] );

                  CV [i] += VISC [ In1 ] * ( RX [i][j] * RX [i][j] + RY [i][j] * RY [i][j] );

                  FV [i] += VISC [ In1 ] * RY [i][j] * RZ [i][j];

                  GV [i] += VISC [ In1 ] * RX [i][j] * RZ [i][j];

                  HV [i] += VISC [ In1 ] * RX [i][j] * RY [i][j];
            }
                  
            /////////////////////////////////////////////////////////////////

genw1:
            if ( IR1 [i] == 0 ) goto genw2;
            
            a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genw1;
      
            a2 = (double) rand () / RAND_MAX;
      
            tmp1 = sqrt ( T / IR1 [i] ) * sqrt ( -2.0 * log (a1) ); 
            
            tmp2 = 2.0 * pi * a2;

            W [J1 + 1] = tmp1 * cos ( tmp2 );

genw2:
            if ( IR2 [i] == 0 ) goto genw3;
            
            a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genw2;
      
            a2 = (double) rand () / RAND_MAX;
      
            tmp1 = sqrt ( T / IR2 [i] ) * sqrt ( -2.0 * log (a1) ); 
            
            tmp2 = 2.0 * pi * a2;

            W [J1 + 2] = tmp1 * cos ( tmp2 );

genw3:
            if ( IR3 [i] == 0 ) goto genvcms1;
            
            a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genw3;
      
            a2 = (double) rand () / RAND_MAX;
      
            tmp1 = sqrt ( T / IR3 [i] ) * sqrt ( -2.0 * log (a1) ); 
            
            tmp2 = 2.0 * pi * a2;

            W [J1 + 3] = tmp1 * cos ( tmp2 );

genvcms1:
            
            a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvcms1;
      
            a2 = (double) rand () / RAND_MAX;
      
            tmp1 = sqrt ( T / RMASS [i] ) * sqrt ( -2.0 * log (a1) ); 
            
            tmp2 = 2.0 * pi * a2;

            VCMS [J1 + 1] = tmp1 * cos ( tmp2 );

genvcms2:
            
            a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvcms2;
      
            a2 = (double) rand () / RAND_MAX;
      
            tmp1 = sqrt ( T / RMASS [i] ) * sqrt ( -2.0 * log (a1) ); 
            
            tmp2 = 2.0 * pi * a2;

            VCMS [J1 + 2] = tmp1 * cos ( tmp2 );

genvcms3:
            
            a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvcms3;
      
            a2 = (double) rand () / RAND_MAX;
      
            tmp1 = sqrt ( T / RMASS [i] ) * sqrt ( -2.0 * log (a1) ); 
            
            tmp2 = 2.0 * pi * a2;

            VCMS [J1 + 3] = tmp1 * cos ( tmp2 );
      }
}

/////////////////////////////////////////////////////////////////////////
////////////////////////////// Nativeinfo ///////////////////////////////
/////////////////////////////////////////////////////////////////////////
void write_ninfo_basestack( FILE* f )
{
      int i, j, k, k_10, k00, k01, k10, k11, k12, k20, k21, k22, k30;
      int ii, jj;
      double r, theta1, theta2, psi, psi1, psi2;
      double psi0, psi10, psi20;
      int istack = 0;
      int ibsdih = 0;
      int itstack = 0;
      int itbsdih = 0;
      int itbsangl = 0;
      
      for (k = 1; k < s3_N + 1; k++)
      {
            i = st_i [k];
            j = st_j [k];

            //ii = ATOM_ST [ k ][ 0 ];
            //jj = ATOM_ST [ k ][ 1 ];
            if (ST_EXCESS[ ATOM_ST[k][0] ] == 1) ii = 2;
            else ii = 1;
            if (ST_EXCESS[ ATOM_ST[k][1] ] == 1) jj = 2;
            else jj = 1;

            k10 = atom_key [i][0];
            k11 = atom_key [i][1];
            k12 = atom_key [i + 1][0];
                  
            k20 = atom_key [j][0];
            k21 = atom_key [j][1];
            k22 = atom_key [j + 1][0];
                        
            if ( fabs ( psi - st_psi [k] ) <= pi ) psi0 = st_psi [k];
            else if ( psi - st_psi [k] > pi ) psi0 = st_psi [k] + 2.0 * pi;
            else psi0 = st_psi [k] - 2.0 * pi;

            if ( fabs ( psi1 - st_psi1 [k] ) <= pi ) psi10 = st_psi1 [k];
            else if ( psi1 - st_psi1 [k] > pi ) psi10 = st_psi1 [k] + 2.0 * pi;
            else psi10 = st_psi1 [k] - 2.0 * pi;

            if ( fabs ( psi2 - st_psi2 [k] ) <= pi ) psi20 = st_psi2 [k];
            else if ( psi2 - st_psi2 [k] > pi ) psi20 = st_psi2 [k] + 2.0 * pi;
            else psi20 = st_psi2 [k] - 2.0 * pi;

            //stacking_bond_force ( k11, k21, 5.00 * st_r2 [k], st_r [k] );
            ++ itstack;
            fprintf(f, "tbs-dist %6i %2i %2i %4i %4i %4i %4i %9.3f %7.2f 5.0 %3i %3i\n",
                                 itstack, 1, 1, k11, k21, k11, k21,
                                 st_E[k], st_r[k], ii, jj);

            //stacking_valence_force ( k10, k11, k21, 1.50 * st_r2 [k], st_theta1 [k] );
            ++ itbsangl;
            fprintf(f, "tbs-angl %6i %6i %2i %2i %4i %4i %4i %4i %4i %4i",
                             itstack, itbsangl, 1,1,
                             k10, k11, k21, k10, k11, k21);
            fprintf(f, " %7.2f %3.1f\n",
                        st_theta1[k]/pi*180.0, 1.5);

            //stacking_valence_force ( k20, k21, k11, 1.50 * st_r2 [k], st_theta2 [k] );
            ++ itbsangl;
            fprintf(f, "tbs-angl %6i %6i %2i %2i %4i %4i %4i %4i %4i %4i",
                             itstack, itbsangl, 1,1,
                             k20, k21, k11, k20, k21, k11);
            fprintf(f, " %7.2f %3.1f\n",
                        st_theta2[k]/pi*180.0, 1.5);

            //stacking_dihedral_force ( k10, k11, k21, k20, 0.15 * st_r2 [k], st_psi0 [k] );
            ++ itbsdih;
            fprintf(f, "tbs-dihd %6i %6i %2i %2i %4i %4i %4i %4i %4i %4i %4i %4i",
                             itstack, itbsdih,1, 1, k10, k11, k21, k20, k10, k11, k21, k20);
            fprintf(f, " %7.2f %4.2f\n",
                        psi0/pi*180.0, 0.15);

            //stacking_dihedral_force ( k21, k11, k10, k12, 0.15 * st_r2 [k], st_psi10 [k] );
            ++ itbsdih;
            fprintf(f, "tbs-dihd %6i %6i %2i %2i %4i %4i %4i %4i %4i %4i %4i %4i",
                             itstack, itbsdih,1, 1, k21, k11, k10, k12, k21, k11, k10, k12);
            fprintf(f, " %7.2f %4.2f\n",
                        psi10/pi*180.0, 0.15);

            //stacking_dihedral_force ( k11, k21, k20, k22, 0.15 * st_r2 [k], st_psi20 [k] );
            ++ itbsdih;
            fprintf(f, "tbs-dihd %6i %6i %2i %2i %4i %4i %4i %4i %4i %4i %4i %4i",
                             itstack, itbsdih,1, 1, k11, k21, k20, k22, k11, k21, k20, k22);
            fprintf(f, " %7.2f %4.2f\n",
                        psi20/pi*180.0, 0.15);
      }
            
      /*
      for (k = 1; k < s3_N + 1; k++)
      {
            if ( e3_stack [k] < 10.0 )
            {
                  i = ATOM_ST [k][0];

                  if ( ST_EXCESS [i] == 1 )
                  {
                        k00 = s3_N + (int) ceil (0.5 * i) - 1;
                        st_status [ k00 ] = 1;
                  }

                  i = ATOM_ST [k][1];

                  if ( ST_EXCESS [i] == 1 )
                  {
                        k00 = s3_N + (int) ceil (0.5 * i) - 1;

                        st_status [ k00 ] = 1;
                  }
            }
      } */

      
      for (k = s3_N + 1; k < st_N + 1; k++)
      {
            //if ( st_status [k] ) { st_status [k] = 0; continue; }

            i = st_i [k];
                  
            k_10 = atom_key [i - 1][0]; // P(n)
            
            k00 = atom_key [i][0];  // S(n)
            k01 = atom_key [i][1];  // B(n)
            k10 = atom_key [i + 1][0];  // P(n+1)
            k20 = atom_key [i + 2][0];  // S(n+1)
            k21 = atom_key [i + 2][1];  // B(n+1)
            k30 = atom_key [i + 3][0];  // P(n+2)

            psi1 = dihedral_angle ( k_10, k00, k10, k20 );
            psi2 = dihedral_angle ( k30, k20, k10, k00 );
            
            if ( fabs ( psi1 + 2.58684 ) <= pi ) psi10 = - 2.58684;
            else if ( psi1 + 2.58684 > pi ) psi10 = - 2.58684 + 2.0 * pi;
            else psi10 = - 2.58684 - 2.0 * pi;

            if ( fabs ( psi2 - 3.07135 ) <= pi ) psi20 = 3.07135;
            else if ( psi2 - 3.07135 > pi ) psi20 = 3.07135 + 2.0 * pi;
            else psi20 = 3.07135 - 2.0 * pi;
            
            ++ istack;
            fprintf(f, "bs-dist %6i %2i %2i %4i %4i %4i %4i %9.3f %7.2f 1.40 %c-%c\n",
                                 istack, 1, 1, k01, k21, k01, k21, st_E[k], st_r[k], 
                                 amino_key [ i ], amino_key[i+2]);
            fflush(f);
            ++ ibsdih;
            fprintf(f, "bs-dihd %6i %6i %2i %2i %4i %4i %4i %4i %4i %4i %4i %4i",
                             istack, ibsdih, 1, 1,  k_10,k00,k10,k20, k_10,k00,k10,k20);
            fprintf(f, " %7.2f %3.1f PSPS\n",
                        psi10/pi*180.0, 4.0);
            ++ ibsdih;
            fprintf(f, "bs-dihd %6i %6i %2i %2i %4i %4i %4i %4i %4i %4i %4i %4i",
                             istack, ibsdih, 1, 1,  k00,k10,k20,k30, k00,k10,k20,k30);
            fprintf(f, " %7.2f %3.1f SPSP\n",
                        psi20/pi*180.0, 4.0);
            fflush(f);
      }
}

///////////////////////////////////////////////////////////////////////////////////

void write_ninfo_basepair ( FILE* f )
{
      int k;
      int ihb = 0;
      int ihbdih = 0;
      int ihbangl = 0;
      
      //fprintf(f, "hb_N= %d\n", hb_N);
      
      ///////////////////////////////////////////
      for (k = 1; k < hb_N + 1; k++)
      {
            ++ ihb;
            if ( hb_dode [k] )
            {
               fprintf(f, "hb-dist %6i %2i %2i %4i %4i %4i %4i %9.3f %7.2f 5.0 T\n",
                                    ihb,1, 1, hb_k11[k], hb_k21[k], hb_k11[k], hb_k21[k],
                                    hb_E[k], hb_r[k]);
            }
            else
            {
               fprintf(f, "hb-dist %6i %2i %2i %4i %4i %4i %4i %9.3f %7.2f 5.0 S\n",
                                    ihb,1, 1, hb_k11[k], hb_k21[k], hb_k11[k], hb_k21[k],
                                    hb_E[k], hb_r[k]);
            }


            //theta1 = valence_angle ( hb_k10 [k], hb_k11 [k], hb_k21 [k] );
            //r2 += 1.50 * ( theta1 - hb_theta1 [k] ) * ( theta1 - hb_theta1 [k] );
            ++ ihbangl;
            fprintf(f, "hb-angl %6i %6i %2i %2i %4i %4i %4i %4i %4i %4i",
                             ihb, ihbangl, 1,1,
                             hb_k10[k], hb_k11[k], hb_k21[k], hb_k10[k], hb_k11[k], hb_k21[k]);
            fprintf(f, " %7.2f %3.1f\n",
                        hb_theta1[k]/pi*180.0,  1.5);

            //theta2 = valence_angle ( hb_k20 [k], hb_k21 [k], hb_k11 [k] );
            //r2 += 1.50 * ( theta2 - hb_theta2 [k] ) * ( theta2 - hb_theta2 [k] );
            ++ ihbangl;
            fprintf(f, "hb-angl %6i %6i %2i %2i %4i %4i %4i %4i %4i %4i",
                             ihb, ihbangl, 1,1,
                             hb_k20[k], hb_k21[k], hb_k11[k], hb_k20[k], hb_k21[k], hb_k11[k]);
            fprintf(f, " %7.2f %3.1f\n",
                        hb_theta2[k]/pi*180.0,  1.5);

            //psi  = dihedral_angle ( hb_k10 [k], hb_k11 [k], hb_k21 [k], hb_k20 [k] );
            //if ( fabs ( psi - hb_psi [k] ) <= pi ) psi0 = hb_psi [k];
            //else if ( psi - hb_psi [k] > pi ) psi0 = hb_psi [k] + 2.0 * pi;
            //else psi0 = hb_psi [k] - 2.0 * pi;
            //r2 += 0.15 * ( psi - psi0 ) * ( psi - psi0 );
            ++ ihbdih;
            fprintf(f, "hb-dihd %6i %6i %2i %2i %4i %4i %4i %4i %4i %4i %4i %4i",
                             ihb, ihbdih,1, 1, hb_k10[k], hb_k11[k], hb_k21[k], hb_k20[k],
                                        hb_k10[k], hb_k11[k], hb_k21[k], hb_k20[k]);
            fprintf(f, " %7.2f %4.2f\n",
                        hb_psi[k]/pi*180.0, 0.15);

            //psi1 = dihedral_angle ( hb_k21 [k], hb_k11 [k], hb_k10 [k], hb_k12 [k] );
            //if ( fabs ( psi1 - hb_psi1 [k] ) <= pi ) psi10 = hb_psi1 [k];
            //else if ( psi1 - hb_psi1 [k] > pi ) psi10 = hb_psi1 [k] + 2.0 * pi;
            //else psi10 = hb_psi1 [k] - 2.0 * pi;
            //r2 += 0.15 * ( psi1 - psi10 ) * ( psi1 - psi10 );
            ++ ihbdih;
            fprintf(f, "hb-dihd %6i %6i %2i %2i %4i %4i %4i %4i %4i %4i %4i %4i",
                             ihb, ihbdih,1, 1, hb_k21[k], hb_k11[k], hb_k10[k], hb_k12[k],
                                        hb_k21[k], hb_k11[k], hb_k10[k], hb_k12[k]);
            fprintf(f, " %7.2f %4.2f\n",
                        hb_psi1[k]/pi*180.0, 0.15);

            //psi2 = dihedral_angle ( hb_k11 [k], hb_k21 [k], hb_k20 [k], hb_k22 [k] );
            //if ( fabs ( psi2 - hb_psi2 [k] ) <= pi ) psi20 = hb_psi2 [k];
            //else if ( psi2 - hb_psi2 [k] > pi ) psi20 = hb_psi2 [k] + 2.0 * pi;
            //else psi20 = hb_psi2 [k] - 2.0 * pi;
            //r2 += 0.15 * ( psi2 - psi20 ) * ( psi2 - psi20 );
            ++ ihbdih;
            fprintf(f, "hb-dihd %6i %6i %2i %2i %4i %4i %4i %4i %4i %4i %4i %4i",
                             ihb, ihbdih,1, 1, hb_k11[k], hb_k21[k], hb_k20[k], hb_k22[k],
                                        hb_k11[k], hb_k21[k], hb_k20[k], hb_k22[k]);
            fprintf(f, " %7.2f %4.2f\n",
                        hb_psi2[k]/pi*180.0, 0.15);
      }
}

void write_ninfo_bond_angle ( FILE* f)
{
    int i, j, k, entry;
    int ibd = 0;
    int iangl = 0;

    for (entry = 1; entry < NBP_atom + 1; entry += 2)
    {
       i = atom_key [ entry ][0];
       j = atom_key [ entry ][1];

       switch ( amino_key [ entry ] )
       {
       case 'A': 
          ++ ibd;
          fprintf(f, "bond %6i %2i %2i %4i %4i %4i %4i %12.4f %3.1f %3.1f %12.4f SA\n",
                      ibd,1,  1,  i,  j,  i,  j,  4.8515,1.0,  1.0,  10.0);

          if ( entry < NBP_atom ) 
          {
             k = atom_key [ entry + 1 ][0];
             ++ iangl;
             fprintf(f, "angl %6i %2i %2i %4i %4i %4i %4i %4i %4i",
                          iangl,1,1,  j,  i,  k,  j,  i,  k);
             fprintf(f, "%12.4f %3.1f %3.1f %12.4f ASP\n",
                         1.9259/pi*180.0,1.0,  1.0,  5.0);
          }
                        
          if ( entry > 1 ) 
          {
             k = atom_key [ entry - 1 ][0];
             ++ iangl;
             fprintf(f, "angl %6i %2i %2i %4i %4i %4i %4i %4i %4i",
                          iangl,1,1,  k,  i,  j,  k,  i,  j);
             fprintf(f, "%12.4f %3.1f %3.1f %12.4f PSA\n",
                         1.7029/pi*180.0,1.0,  1.0,  5.0);
          }

          break;

       case 'C': 
          ++ ibd;
          fprintf(f, "bond %6i %2i %2i %4i %4i %4i %4i %12.4f %3.1f %3.1f %12.4f SC\n",
                           ibd,1,  1,  i,  j,  i,  j,  4.2738,1.0,  1.0,  10.0);

          if ( entry < NBP_atom ) 
          {
             k = atom_key [ entry + 1 ][0];
             ++ iangl;
             fprintf(f, "angl %6i %2i %2i %4i %4i %4i %4i %4i %4i",
                             iangl,1,1,  j,  i,  k,  j,  i,  k);
             fprintf(f, "%12.4f %3.1f %3.1f %12.4f CSP\n",
                          1.9655/pi*180.0,1.0,  1.0,  5.0);
          }
                        
          if ( entry > 1 ) 
          {
              k = atom_key [ entry - 1 ][0];
              ++ iangl;
              fprintf(f, "angl %6i %2i %2i %4i %4i %4i %4i %4i %4i",
                           iangl,1,1,  k,  i,  j,  k,  i,  j);
              fprintf(f, "%12.4f %3.1f %3.1f %12.4f PSC\n",
                         1.5803/pi*180.0,1.0,  1.0,   5.0);
          }

          break;

       case 'G': 
          ++ ibd;
          fprintf(f, "bond %6i %2i %2i %4i %4i %4i %4i %12.4f %3.1f %3.1f %12.4f SG\n",
                         ibd,1,  1,  i,  j,  i,  j,  4.9659,1.0,  1.0,  10.0);

          if ( entry < NBP_atom ) 
          {
             k = atom_key [ entry + 1 ][0];
             ++ iangl;
             fprintf(f, "angl %6i %2i %2i %4i %4i %4i %4i %4i %4i",
                           iangl,1,1,  j,  i,  k,  j,  i,  k);
             fprintf(f, "%12.4f %3.1f %3.1f %12.4f GSP\n",
                           1.9150/pi*180.0,1.0,  1.0,  5.0);
          }
                        
          if ( entry > 1 ) 
          {
             k = atom_key [ entry - 1 ][0];
             ++ iangl;
             fprintf(f, "angl %6i %2i %2i %4i %4i %4i %4i %4i %4i",
                             iangl,1,1,  k,  i,  j,  k,  i,  j);
             fprintf(f, "%12.4f %3.1f %3.1f %12.4f PSG\n",
                          1.7690/pi*180.0,1.0,  1.0,  5.0);
          }

          break;

       case 'U': 
          ++ ibd;
          fprintf(f, "bond %6i %2i %2i %4i %4i %4i %4i %12.4f %3.1f %3.1f %12.4f SU\n",
                         ibd,1,  1,  i,  j,  i,  j,  4.2733,1.0,  1.0,  10.0);

          if ( entry < NBP_atom ) 
          {
              k = atom_key [ entry + 1 ][0];
              ++ iangl;
              fprintf(f, "angl %6i %2i %2i %4i %4i %4i %4i %4i %4i",
                             iangl,1,1,  j,  i,  k,  j,  i,  k);
              fprintf(f, "%12.4f %3.1f %3.1f %12.4f USP\n",
                          1.9663/pi*180.0,1.0,  1.0,  5.0);
           }
                        
           if ( entry > 1 ) 
           {
               k = atom_key [ entry - 1 ][0];
               ++ iangl;
               fprintf(f, "angl %6i %2i %2i %4i %4i %4i %4i %4i %4i",
                             iangl,1,1,  k,  i,  j,  k,  i,  j);
               fprintf(f, "%12.4f %3.1f %3.1f %12.4f PSU\n",
                          1.5735/pi*180.0,1.0,  1.0,  5.0);
            }

            break;
             
         default: 
            break;
         }
      }
            
      /////////////////////////////////////////

      for (entry = 2; entry < NBP_atom + 1; entry += 2)
      {
            i = atom_key [ entry ][0];
            j = atom_key [ entry - 1 ][0];

        ++ ibd;
        fprintf(f, "bond %6i %2i %2i %4i %4i %4i %4i %12.4f %3.1f %3.1f %12.4f SP\n",
                         ibd,1,  1,  j,  i,  j,  i,  3.8157,1.0,  1.0,  64.0);

            if ( entry > 2 )
            {
                  k = atom_key [ entry - 2 ][0];
            ++ iangl;
            fprintf(f, "angl %6i %2i %2i %4i %4i %4i %4i %4i %4i",
                             iangl,1,1,  k,  j,  i,  k,  j,  i);
            fprintf(f, "%12.4f %3.1f %3.1f %12.4f PSP\n",
                        1.4440/pi*180.0, 1.0,   1.0,  20.0);
            }
      }

      /////////////////////////////////////////

      for (entry = 3; entry < NBP_atom + 1; entry += 2)
      {
            i = atom_key [ entry ][0];
            j = atom_key [ entry - 1 ][0];
            k = atom_key [ entry - 2 ][0];

        ++ ibd;
        fprintf(f, "bond %6i %2i %2i %4i %4i %4i %4i %12.4f %3.1f %3.1f %12.4f PS\n",
                         ibd,1,  1,  j,  i,  j,  i,  4.6010,1.0,  1.0,   23.0);
        ++ iangl;
        fprintf(f, "angl %6i %2i %2i %4i %4i %4i %4i %4i %4i",
                         iangl,1,1,  k,  j,  i,  k,  j,  i);
        fprintf(f, "%12.4f %3.1f %3.1f %12.4f SPS\n",
                    1.5256/pi*180.0,1.0,  1.0,  20.0);
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////MAIN ROUTINE//////////////////////////////////////////////////////

int main ( int argc, char ** argv )
{
      FILE * f_con, * f_ene, *f_dcd, *f_ninfo, *f_cgpdb, *f_psf;
      FILE * f1, * f2, * f3, * f4;
      
      int In1, In2, In12, InInt, seed;
      
      long klok_old = 0;

      char syscall [500];

      double * erna = NULL;

      /////////////////////////////////////////////////////////////
      
      char * file_PDB = "16SCD.pdb";
      char * file_ninfo = "16SCD_ninfo.ninfo";
      char * file_cgpdb = "16SCD_ninfo.cg.pdb";
      char * file_psf   = "16SCD_ninfo.psf";

      char * file_unprocessed_bonds_16SCD_16SCD = "16SCD_16SCD_unprocessed_bonds.dat";

      char * file_unprocessed_stacks = "16SCD_16SCD_unprocessed_stacks.dat";
      
      /////////////////////////////////////////////////////////////
      
      char file_contacts [200] = "Contacts_16SCD_EPS_Mg_Cl_K_";

      char file_coordinates [200] = "Coordinates_16SCD_EPS_Mg_Cl_K_";
      
      char file_energies [200] = "Energies_16SCD_EPS_Mg_Cl_K_";
      
      /////////////////////////////////////////////////////////////

      char tmp_contacts [200] = "/tmp/Contacts_16SCD_EPS_Mg_Cl_K_";

      char tmp_coordinates [200] = "/tmp/Coordinates_16SCD_EPS_Mg_Cl_K_";
      
      char tmp_energies [200] = "/tmp/Energies_16SCD_EPS_Mg_Cl_K_";
      
      f_trace = fopen( "trace", "w");
      f_atoms = fopen( "atoms", "w");
      /////////////////////////////////////////////////////////////

      side_X = 350.0; side_Y = 350.0; side_Z = 350.0;

      /////////////////////////////////////////////////////////////
      
      for (In1 = 1; In1 < argc; In1++) 
      {
            if ( argv [ In1 ][ 0 ] == '-' ) 
            {
                  switch ( argv [ In1 ][ 1 ] ) 
                  {
                        case 'D': ss_D = atof ( argv [ In1 ] + 2 ); break;
                        
                        case 'E': hs_D = atof ( argv [ In1 ] + 2 ); break;

                        case 'd': st_D = atof ( argv [ In1 ] + 2 ); break;

                        case 'T': T = atof ( argv [ In1 ] + 2 ); break;

                        case 'c': cM_crwd [0] = atof ( argv [ In1 ] + 2 ); break;

                        case 'f': cM_crwd [2] = atof ( argv [ In1 ] + 2 ); break;

                        case 's': seed = atoi ( argv [ In1 ] + 2 ); break;
                        
                        default: printf ( "unrecognized argument %s\n", argv [ In1 ] ); break;
                  }
            }
      }
      
      srand ( (unsigned) ( T * 1000 * seed ) );

      //////////////////////////////////////////////////////////////

      lB = 332.0637090 * ( 0.03273600947 * T - 0.00669625750 );
      fprintf (f_trace, "lB = %d\n", lB);

      //////////////////////////////////////////////////////////////

      sprintf ( file_contacts + strlen ( file_contacts ), "D%.3f_E%.3f_d%.3f_T%.3f_cM%.4f_f%.3f_s%d.dat", ss_D, hs_D, st_D, T, cM_crwd [0], cM_crwd [2], seed );
      
      sprintf ( file_coordinates + strlen ( file_coordinates ), "D%.3f_E%.3f_d%.3f_T%.3f_cM%.4f_f%.3f_s%d.dat", ss_D, hs_D, st_D, T, cM_crwd [0], cM_crwd [2], seed );

      sprintf ( file_energies + strlen ( file_energies ), "D%.3f_E%.3f_d%.3f_T%.3f_cM%.4f_f%.3f_s%d.dat", ss_D, hs_D, st_D, T, cM_crwd [0], cM_crwd [2], seed );
      
      //////////////////////////////////////////////////////////////

      sprintf ( tmp_contacts + strlen ( tmp_contacts ), "D%.3f_E%.3f_d%.3f_T%.3f_cM%.4f_f%.3f_s%d.dat", ss_D, hs_D, st_D, T, cM_crwd [0], cM_crwd [2], seed );
      
      sprintf ( tmp_coordinates + strlen ( tmp_coordinates ), "D%.3f_E%.3f_d%.3f_T%.3f_cM%.4f_f%.3f_s%d.dat", ss_D, hs_D, st_D, T, cM_crwd [0], cM_crwd [2], seed );

      sprintf ( tmp_energies + strlen ( tmp_energies ), "D%.3f_E%.3f_d%.3f_T%.3f_cM%.4f_f%.3f_s%d.dat", ss_D, hs_D, st_D, T, cM_crwd [0], cM_crwd [2], seed );

      //////////////////////////////////////////////////////////////

      for (In1 = 0; In1 < ns; In1++)
      {
            for (In2 = In1; In2 < ns; In2++)
            {
                  In12 = ns * In1 - In1 * (In1 + 1) / 2 + In2;

                  PD [ In12 ] = ( D [ In1 ] + D [ In2 ] ) / 2;

                  PD_sq [ In12 ] = PD [ In12 ] * PD [ In12 ];
            }
      }
      
      //////////////////////////////////////////////////////////////
      
      read_maxi_key ( "MAXIKEY_Mg_Cl_K.dat" );
      printf ("read_maxi_key\n");
      
      //////////////////////////////////////////////////////////////
      
      list_crowder_types ();

      /////////////////////////////////////////////////////////////////

      read_PDB ( file_PDB );
      printf ("read_PDB\n");

      /////////////////////////////////////////////////////////////////

      initialize_HB_atoms ();

      initialize_unprocessed_bonds ( file_unprocessed_bonds_16SCD_16SCD );

      initialize_unprocessed_stacks ( file_unprocessed_stacks );

      /////////////////////////////////////////////////////////////////

      set_box ( 500 );
      printf ("set_box\n");
            
      /////////////////////////////////////////////////////////////////

      N_crwd [0] = (int) ( cM_crwd [0] * 6.022 * 0.0001 * side_XYZ );
      N_crwd [2] = (int) ( cM_crwd [2] * 6.022 * 0.0001 * side_XYZ );
      N_crwd [1] = 2 * N_crwd [0] + N_crwd [2] - ( NBP_atom - 1 ) / 2;

      fprintf( f_trace, "cM_crwd[0]= %f\n", cM_crwd[0]);
      fprintf( f_trace, "cM_crwd[1]= %f\n", cM_crwd[1]);
      fprintf( f_trace, "cM_crwd[2]= %f\n", cM_crwd[2]);
      fprintf( f_trace, "N_crwd[0]= %d\n", N_crwd[0]);
      fprintf( f_trace, "N_crwd[1]= %d\n", N_crwd[1]);
      fprintf( f_trace, "N_crwd[2]= %d\n", N_crwd[2]);
                  
      for (In1 = 0; In1 < n_crwd; In1++) 
            
            for (In2 = 1; In2 < N_crwd [ In1 ] + 1; In2++) 
                        
                  add_crowder ( In1 );
            
      /////////////////////////////////////////////////////////////////
            
      f1 = fopen ( file_coordinates, "r" );

      if (f1)
      {
            do
            {
                  klok_old += 1000000;
                  
                  fread ( x, sizeof(double), NT_atom + 1, f1 ); 
      
                  fread ( y, sizeof(double), NT_atom + 1, f1 );
      
                  fread ( z, sizeof(double), NT_atom + 1, f1 );
            }
            while ( !feof(f1) );

            fclose (f1);

            klok_old -= 1000000;

            for (In1 = 1; In1 < NT_atom + 1; In1++) 
            {
                  x_old [ In1 ] = x [ In1 ] - vx [ In1 ] * h;
      
                  y_old [ In1 ] = y [ In1 ] - vy [ In1 ] * h;

                  z_old [ In1 ] = z [ In1 ] - vz [ In1 ] * h;
            }
      }
      
      /////////////////////////////////////////////////////////////////

      fsend = (double *) calloc ( 3 * NT_atom, sizeof(double) );

      freceive = (double *) calloc ( 3 * NT_atom, sizeof(double) );

      /////////////////////////////////////////////////////////////////
      
      //MPI_Init ( &argc, &argv );
      //MPI_Comm_rank ( MPI_COMM_WORLD, &my_rank );
      //MPI_Comm_size ( MPI_COMM_WORLD, &comm_size );
      my_rank = 0; 
      comm_size = 1;

      //////////////////////////////////////////////////////////////
      dr1 = 3.0;
      R1_LIST [0] = 12.0 + dr1; R2_LIST [0] = R1_LIST [0] * R1_LIST [0];
      R1_LIST [1] = D_CT_CUTOFF + dr1; R2_LIST [1] = R1_LIST [1] * R1_LIST [1];

      list_content [0] = (int *) calloc ( 100000 + 1, sizeof(int) );
      list_content [1] = (int *) calloc ( 500000 + 1, sizeof(int) );

      populate_lists ();
      PROGRESS = 0;
      
      /////////////////////////////////////////////////////////////////
      
      //klok_old = 0;
     
      // ninfo
      f_ninfo = fopen( file_ninfo, "w");
      write_ninfo_bond_angle( f_ninfo );
      write_ninfo_basestack( f_ninfo );
      write_ninfo_basepair( f_ninfo );
      fclose(f_ninfo);

      // PDB
      f_cgpdb = fopen( file_cgpdb, "w");
      fprintf( f_cgpdb, "RECORD   N_amino   %5i\n", N_amino);
      fprintf( f_cgpdb, "RECORD   NBP_atom  %5i\n", NBP_atom);
      fprintf( f_cgpdb, "RECORD   NB_atom   %5i\n", NB_atom);
      fprintf( f_cgpdb, "RECORD   NTP_atom  %5i\n", NTP_atom);
      fprintf( f_cgpdb, "RECORD   NpTP_atom %5i\n", NpTP_atom);
      fprintf( f_cgpdb, "RECORD   NT_atom   %5i\n", NT_atom);

      for (In1 = 1; In1 < NTP_atom + 1; In1++)
      {
          fprintf ( f_cgpdb, "ATOM  %5i", In1);
          fprintf ( f_cgpdb, " ");

          int i_amino = (In1 / 3) * 2 + 1;
          if (In1 % 3 == 1) 
          {
              fprintf ( f_cgpdb, " S  ");
              fprintf ( f_cgpdb, " R%c ", amino_key[i_amino]);
              fprintf ( f_cgpdb, " A");
              fprintf ( f_cgpdb, "%4i", In1/3+1);
          }
          else if (In1 % 3 == 2)
          {
              fprintf ( f_cgpdb, " %cb ", amino_key[i_amino]);
              fprintf ( f_cgpdb, " R%c ", amino_key[i_amino]);
              fprintf ( f_cgpdb, " A");
              fprintf ( f_cgpdb, "%4i", In1/3+1);
          }
          else
          {
              fprintf ( f_cgpdb, " P  ");
              fprintf ( f_cgpdb, " R%c ", amino_key[i_amino]);
              fprintf ( f_cgpdb, " A");
              fprintf ( f_cgpdb, "%4i", In1/3+1);
          }
          fprintf ( f_cgpdb , "    ");
          fprintf ( f_cgpdb , "%8.3f", x[In1]);
          fprintf ( f_cgpdb , "%8.3f", y[In1]);
          fprintf ( f_cgpdb , "%8.3f", z[In1]);
          fprintf ( f_cgpdb , "\n");
      }
      fprintf ( f_cgpdb, "TER\n");

      for (In1 = NTP_atom + 1; In1 < NT_atom + 1; In1++)
      {
          fprintf ( f_cgpdb, "ATOM  %5i", In1);
          fprintf ( f_cgpdb, " ");

          int i_tp = crowder_key[ N_amino + (In1 - NTP_atom) ];
          if (i_tp == 0)
          {
              fprintf ( f_cgpdb, "MG  ");
              fprintf ( f_cgpdb, "  MG");
              fprintf ( f_cgpdb, " B");
              fprintf ( f_cgpdb, "%4i", In1);
          }
          else if (i_tp == 1)
          {
              fprintf ( f_cgpdb, " Cl ");
              fprintf ( f_cgpdb, "  Cl");
              fprintf ( f_cgpdb, " B");
              fprintf ( f_cgpdb, "%4i", In1);
          }
          else if (i_tp == 2)
          {
              fprintf ( f_cgpdb, " K  ");
              fprintf ( f_cgpdb, "  K ");
              fprintf ( f_cgpdb, " B");
              fprintf ( f_cgpdb, "%4i", In1);
          }
          else
          {
              std::cout << "Error";
          }
          fprintf ( f_cgpdb , "    ");
          fprintf ( f_cgpdb , "%8.3f", x[In1]);
          fprintf ( f_cgpdb , "%8.3f", y[In1]);
          fprintf ( f_cgpdb , "%8.3f", z[In1]);
          fprintf ( f_cgpdb , "\n");
      }
      fprintf ( f_cgpdb, "TER\n");
      fclose(f_cgpdb);

      f_psf = fopen( file_psf, "w");
      fprintf ( f_psf , "%8i !NATOM\n", NT_atom );

      for (In1 = 1; In1 < NTP_atom + 1; In1++)
      {
          fprintf ( f_psf, "%8i", In1);

          int i_amino = (In1 / 3) * 2 + 1;
          if (In1 % 3 == 1) 
          {
              fprintf ( f_psf, "R%c ", amino_key[i_amino]);
              fprintf ( f_psf, "%4i",  In1/3+1);
              fprintf ( f_psf, "      R%c", amino_key[i_amino]);
              fprintf ( f_psf, "    S ");
              fprintf ( f_psf, "  S ");
              fprintf ( f_psf, "    ");
              fprintf ( f_psf, "%11.6f", 0.0);
              fprintf ( f_psf, "%14.4f", 10.0);
          }
          else if (In1 % 3 == 2)
          {
              fprintf ( f_psf, "R%c ", amino_key[i_amino]);
              fprintf ( f_psf, "%4i",  In1/3+1);
              fprintf ( f_psf, "      R%c", amino_key[i_amino]);
              fprintf ( f_psf, "    %cb", amino_key[i_amino]);
              fprintf ( f_psf, "  B ");
              fprintf ( f_psf, "    ");
              fprintf ( f_psf, "%11.6f", 0.0);
              fprintf ( f_psf, "%14.4f", 10.0);
          }
          else
          {
              fprintf ( f_psf, "R%c ", amino_key[i_amino]);
              fprintf ( f_psf, "%4i",  In1/3+1);
              fprintf ( f_psf, "      R%c", amino_key[i_amino]);
              fprintf ( f_psf, "    P ", amino_key[i_amino]);
              fprintf ( f_psf, "  P ");
              fprintf ( f_psf, "    ");
              fprintf ( f_psf, "%11.6f", -1.0);
              fprintf ( f_psf, "%14.4f", 10.0);
          }
          fprintf ( f_psf, "           0\n");
      }

      for (In1 = NTP_atom + 1; In1 < NT_atom + 1; In1++)
      {
          fprintf ( f_psf, "%8i", In1);

          int i_tp = crowder_key[ N_amino + (In1 - NTP_atom) ];
          if (i_tp == 0)
          {
              fprintf ( f_psf, "Mg ");
              fprintf ( f_psf, "%4i",  NTP_atom/3 + In1 - NTP_atom +1);
              fprintf ( f_psf, "      Mg");
              fprintf ( f_psf, "    Mg");
              fprintf ( f_psf, "  Mg");
              fprintf ( f_psf, "    ");
              fprintf ( f_psf, "%11.6f", 2.0);
              fprintf ( f_psf, "%14.4f", 10.0);
          }
          else if (i_tp == 1)
          {
              fprintf ( f_psf, "Cl ");
              fprintf ( f_psf, "%4i",  NTP_atom/3 + In1 - NTP_atom +1);
              fprintf ( f_psf, "      Cl");
              fprintf ( f_psf, "    Cl");
              fprintf ( f_psf, "  Cl");
              fprintf ( f_psf, "    ");
              fprintf ( f_psf, "%11.6f", -1.0);
              fprintf ( f_psf, "%14.4f", 10.0);
          }
          else if (i_tp == 2)
          {
              fprintf ( f_psf, "K  ");
              fprintf ( f_psf, "%4i",  NTP_atom/3 + In1 - NTP_atom +1);
              fprintf ( f_psf, "      K ");
              fprintf ( f_psf, "    K ");
              fprintf ( f_psf, "  K ");
              fprintf ( f_psf, "    ");
              fprintf ( f_psf, "%11.6f", 1.0);
              fprintf ( f_psf, "%14.4f", 10.0);
          }
          else
          {
              std::cout << "Error";
          }
          fprintf ( f_psf, "           0\n");
      }

      fprintf ( f_psf , "%8i !NBOND\n", NTP_atom-1 );
      int icount = 0;
      for (In1 = 2; In1<=NTP_atom ; ++In1 )
      {
          ++ icount;
          if (In1 % 3 == 1)  // S
          {
              fprintf( f_psf, "%8i%8i", In1-1, In1);
          }
          else if (In1 % 3 == 2) // B
          {
              fprintf( f_psf, "%8i%8i", In1-1, In1);
          }
          else  // P
          {
              fprintf( f_psf, "%8i%8i", In1-2, In1);
          }

          if (icount == 4) 
          {
              fprintf ( f_psf, "\n" );
              icount = 0;
          }
      }
                  
      fclose( f_psf );
      fclose( f_trace) ;
      fclose( f_atoms) ;
      return(0);

      if ( my_rank == 0 ) 
      {
            erna = (double *) calloc ( s3_N + 3, sizeof(double) );

            /////////////////////////////////////////////////////////////////
            
            sprintf ( syscall, "cp %s %s", file_contacts, tmp_contacts );
            
            system ( syscall );

            sprintf ( syscall, "cp %s %s", file_coordinates, tmp_coordinates );
            
            system ( syscall );

            sprintf ( syscall, "cp %s %s", file_energies, tmp_energies );
            
            system ( syscall );

            //////////////////////////////////////////////////////////////
            
            f1 = fopen ( tmp_contacts, "a" );

            f2 = fopen ( tmp_coordinates, "a" );

            f3 = fopen ( tmp_energies, "a" );
      }

      //////////////////////////////////////////////////////////////
      
      for (klok = klok_old + 1; klok < 3000000001; klok++)
      {
            // printf ("klok %i\n", klok);

            if ( (klok % 1000000) == 0 ) 
            {
                  st_energy [ 0 ] = E2_STACK;

                  st_energy [ s3_N + 1 ] = E2_HB;

                  st_energy [ s3_N + 2 ] = E3_HB;

                  //MPI_Reduce ( st_energy, erna, s3_N + 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

                  //////////////////////////////////////////////////////////////
                  
                  if ( my_rank == 0 )
                  {
                        for (In1 = 1; In1 < s3_N + 1; In1++) fprintf ( f1, "%+.3f ", erna [ In1 ] );
                        
                        fprintf ( f1, "\n" );
                        
                        fflush (f1);

                        //////////////////////////////////////////////////////////////
                        
                        fwrite ( x, sizeof(double), NT_atom + 1, f2 ); 
      
                        fwrite ( y, sizeof(double), NT_atom + 1, f2 );
      
                        fwrite ( z, sizeof(double), NT_atom + 1, f2 );

                        fflush (f2);

                        //////////////////////////////////////////////////////////////

                        E3_STACK = 0;
                        
                        for (In1 = 1; In1 < s3_N + 1; In1++) E3_STACK += erna [ In1 ];

                        fprintf ( f3, "%le %le %le %le %le\n", erna [ s3_N + 1 ], erna [ s3_N + 2 ], erna [ 0 ], E3_STACK, radius_of_gyration ( 1, NTP_atom ) );

                        fflush (f3);

                        //////////////////////////////////////////////////////////////

                        if ( klok == 3000000000 )
                        {
                              sprintf ( syscall, "cp %s Contacts_16SCD_EPS_Mg_Cl_K_D%.3f_E%.3f_d%.3f_T%.3f_cM%.4f_f%.3f_s%d_%d.dat", tmp_contacts, ss_D, hs_D, st_D, T, cM_crwd [0], cM_crwd [2], seed, (int) (klok / 5000000) );
            
                              system ( syscall );

                              //////////////////////////////////////////////////////////////
                              
                              sprintf ( syscall, "cp %s Coordinates_16SCD_EPS_Mg_Cl_K_D%.3f_E%.3f_d%.3f_T%.3f_cM%.4f_f%.3f_s%d_%d.dat", tmp_coordinates, ss_D, hs_D, st_D, T, cM_crwd [0], cM_crwd [2], seed, (int) (klok / 5000000) );
            
                              system ( syscall );

                              //////////////////////////////////////////////////////////////
                              
                              sprintf ( syscall, "cp %s Energies_16SCD_EPS_Mg_Cl_K_D%.3f_E%.3f_d%.3f_T%.3f_cM%.4f_f%.3f_s%d_%d.dat", tmp_energies, ss_D, hs_D, st_D, T, cM_crwd [0], cM_crwd [2], seed, (int) (klok / 5000000) );
            
                              system ( syscall );
                        }

                        //////////////////////////////////////////////////////////////

                        else if ( (klok % 5000000) == 0 )
                        {
                              sprintf ( syscall, "cp %s Contacts_16SCD_EPS_Mg_Cl_K_D%.3f_E%.3f_d%.3f_T%.3f_cM%.4f_f%.3f_s%d_%d.dat &", tmp_contacts, ss_D, hs_D, st_D, T, cM_crwd [0], cM_crwd [2], seed, (int) (klok / 5000000) );
            
                              system ( syscall );
                              
                              //////////////////////////////////////////////////////////////

                              sprintf ( syscall, "cp %s Coordinates_16SCD_EPS_Mg_Cl_K_D%.3f_E%.3f_d%.3f_T%.3f_cM%.4f_f%.3f_s%d_%d.dat &", tmp_coordinates, ss_D, hs_D, st_D, T, cM_crwd [0], cM_crwd [2], seed, (int) (klok / 5000000) );
            
                              system ( syscall );

                              //////////////////////////////////////////////////////////////
                              
                              sprintf ( syscall, "cp %s Energies_16SCD_EPS_Mg_Cl_K_D%.3f_E%.3f_d%.3f_T%.3f_cM%.4f_f%.3f_s%d_%d.dat &", tmp_energies, ss_D, hs_D, st_D, T, cM_crwd [0], cM_crwd [2], seed, (int) (klok / 5000000) );
            
                              system ( syscall );

                              //////////////////////////////////////////////////////////////
                              
                              sprintf ( syscall, "rm -f Contacts_16SCD_EPS_Mg_Cl_K_D%.3f_E%.3f_d%.3f_T%.3f_cM%.4f_f%.3f_s%d_%d.dat &", ss_D, hs_D, st_D, T, cM_crwd [0], cM_crwd [2], seed, (int) (klok / 5000000) - 3 );
            
                              system ( syscall );
                              
                              //////////////////////////////////////////////////////////////

                              sprintf ( syscall, "rm -f Coordinates_16SCD_EPS_Mg_Cl_K_D%.3f_E%.3f_d%.3f_T%.3f_cM%.4f_f%.3f_s%d_%d.dat &", ss_D, hs_D, st_D, T, cM_crwd [0], cM_crwd [2], seed, (int) (klok / 5000000) - 3 );
            
                              system ( syscall );

                              //////////////////////////////////////////////////////////////
                              
                              sprintf ( syscall, "rm -f Energies_16SCD_EPS_Mg_Cl_K_D%.3f_E%.3f_d%.3f_T%.3f_cM%.4f_f%.3f_s%d_%d.dat &", ss_D, hs_D, st_D, T, cM_crwd [0], cM_crwd [2], seed, (int) (klok / 5000000) - 3 );
            
                              system ( syscall );
                        }
                  }
            }

            //////////////////////////////////////////////////////////////

            deterministic_forces_MPI ();

            full_forces ();

            move_rigid_units ();

            check_shifts ();
      }
      
      //////////////////////////////////////////////////////////////

      if ( my_rank == 0 ) 
      {
            fclose ( f1 );

            fclose ( f2 );

            fclose ( f3 );

            free ( erna );
      }

      //////////////////////////////////////////////////////////////

      //MPI_Finalize ();

      //////////////////////////////////////////////////////////////
      
      free ( NS_atom );

      free ( RIGID_SET );

      free ( RIGID_END );
      
      for (In1 = 1; In1 < NB_atom + 1; In1++) free ( atom_key [ In1 ] );
            
      free ( atom_key );

      free ( amino_key );

      free ( crowder_key );
      
      free ( part_key );

      free ( maxi_key );

      free ( residue_mass );

      for (In1 = 0; In1 < 5; In1++) free ( residue_content [In1] );

      free ( residue_content );
      
      free_stacks ();
      
      free_hydrogen_bonds ();

      ////////////////////////////////////////
      
      for (In1 = 0; In1 < n_crwd; In1++)
      {
            free ( RX_crwd [ In1 ] );
            
            free ( RY_crwd [ In1 ] );
            
            free ( RZ_crwd [ In1 ] );
            
            free ( part_key_crwd [ In1 ] );

            free ( maxi_key_crwd [ In1 ] );
      }

      free ( crowder_mass );

      for (In1 = 0; In1 < n_crwd; In1++) free ( crowder_content [In1] );

      free ( crowder_content );

      ////////////////////////////////////////
      
      free ( CMS );
      
      free ( VCMS );
      
      ////////////////////////////////////////

      free ( AXES );
      
      free ( W );
      
      ////////////////////////////////////////
      
      free ( x ); free ( y ); free ( z );

      free ( x_old ); free ( y_old ); free ( z_old );

      free ( vx ); free ( vy ); free ( vz );
      
      free ( fx ); free ( fy ); free ( fz );

      free ( flx ); free ( fly ); free ( flz );
      
      free ( fsend ); free ( freceive );
      
      ////////////////////////////////////////
      
      free ( INDX ); free ( JNDX );
      
      free ( RMASS );
      
      free ( RVISC );
      
      free ( RXYZ );

      for (In1 = 1; In1 < NB_atom + 1; In1++) { free ( RX [ In1 ] ); free ( RY [ In1 ] ); free ( RZ [ In1 ] ); }
      
      free ( RX ); free ( RY ); free ( RZ );

      free ( IR1 ); free ( IR2 ); free ( IR3 );

      free ( AV ); free ( BV ); free ( CV );
      
      free ( FV ); free ( GV ); free ( HV );

      ////////////////////////////////////////
      
      for (InInt = 0; InInt < npi; InInt++) free ( list_content [ InInt ] );
            
      ////////////////////////////////////////
            
      for (In1 = 1; In1 < McT + 1; In1++) free ( cell_content [ In1 ] ); 
      
      free (cell_content); free (cell_mass);

      ////////////////////////////////////////

      return (0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//struct timeb t0;
      
      //struct timeb t1;

      /////////////////////////////////////////////////////////////
      
      //double sc [1];
      
      //double maxwell_tail;
      /*ftime (&t0);
      
      populate_lists (); //populate_lists_MPI (0);

      ftime (&t1);
      
      f4 = fopen ( "Clock.dat", "a" );

      fprintf ( f4, "%d: %u\n\n", my_rank, t1.millitm - t0.millitm );

      fclose (f4);
      
      PROGRESS = 0;*/
      
/////////////////////////////////////////////////////////////////

/*void populate_lists_MPI ( int root )
{
      FILE * f1;
      
      int k, InInt;
      
      int m1, m2, m3;

      int n1, n2, n3;

      int J1, J2; 
            
      int i1, i2;
      
      int j1, j2; 

      int tmp, McYZ;

      double r2, xx [3];
      
      int index [13][3] = { { -1, -1, -1 }, { -1, -1, 0 },  { -1, -1, 1 },  { -1, 0, -1 }, { -1, 0, 0 }, { -1, 0, 1 }, 
      
      { -1, 1, -1 }, { -1, 1, 0 }, { -1, 1, 1 }, { 0, -1, -1 }, { 0, -1, 0 }, { 0, -1, 1 }, { 0, 0, -1 } };
      
      int np_busy, load, offset;

      int temp_mass = 0, * temp_content = NULL;

      int * rcounts = NULL, * displs = NULL;

      //////////////////////////////////////////
      
      list_mass [0] = 0;

      list_mass [1] = 0;
      
      //////////////////////////////////////////

      separate_in_cells ();
      
      //////////////////////////////////////////

      McYZ = McY * McZ;

      np_busy = McT % comm_size;

      if ( my_rank < np_busy ) { load = McT / comm_size + 1; offset = my_rank * load; }
      
      else { load = McT / comm_size; offset = my_rank * load + np_busy; }

      //////////////////////////////////////////

      for (J1 = offset + 1; J1 < offset + load + 1; J1++)
      {
            J2 = J1 % McYZ;
      
            if ( J2 == 0 ) { m1 = J1 / McYZ; m2 = McY; m3 = McZ; }

            else 
            { 
                  m1 = J1 / McYZ + 1; tmp = J2 % McZ; 
            
                  if ( tmp == 0 ) { m2 = J2 / McZ; m3 = McZ; } 
                  
                  else { m2 = 1 + J2 / McZ; m3 = tmp; } 
            }
            
            for (k = 0; k < 13; k++)
            {
                  n1 = m1 + index [k][0]; 
            
                  if (n1 > McX) { n1 -= McX; if ( n1 == m1 || n1 == m1 - 1 ) continue; }
            
                  else if (n1 < 1) { n1 += McX; if ( n1 == m1 || n1 == m1 + 1 ) continue; }

                  n2 = m2 + index [k][1]; 
            
                  if (n2 > McY) { n2 -= McY; if ( n2 == m2 || n2 == m2 - 1 ) continue; }
            
                  else if (n2 < 1) { n2 += McY; if ( n2 == m2 || n2 == m2 + 1 ) continue; }
            
                  n3 = m3 + index [k][2]; 
            
                  if (n3 > McZ) { n3 -= McZ; if ( n3 == m3 || n3 == m3 - 1 ) continue; } 
            
                  else if (n3 < 1) { n3 += McZ; if ( n3 == m3 || n3 == m3 + 1 ) continue; }
                              
                  J2 = (n1 - 1) * McY * McZ + (n2 - 1) * McZ + n3;
                              
                  for (j1 = 1; j1 < cell_mass [J1] + 1; j1++)
                  {
                        i1 = cell_content [J1][j1];

                        n1 = part_key [i1];
                        
                        for (j2 = 1; j2 < cell_mass [J2] + 1; j2++)
                        {
                              i2 = cell_content [J2][j2];

                              n2 = part_key [i2];
                                          
                              ///////////////////////////////////////////////////////

                              CC_vector ( i1, i2, xx );
                        
                              r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

                              ///////////////////////////////////////////////////////

                              if ( r2 <= R2_LIST [0] )
                              {
                                    list_mass [0] += 2;
                                                
                                    list_content [0][ list_mass [0] - 1 ] = i1; 
                                                
                                    list_content [0][ list_mass [0] ] = i2;
                              }

                              ///////////////////////////////////////////////////////

                              else if ( r2 <= R2_LIST [1] && n1 != 1 && n1 != 2 && n2 != 1 && n2 != 2 )
                              {
                                    list_mass [1] += 2;
                                                
                                    list_content [1][ list_mass [1] - 1 ] = i1; 
                                                
                                    list_content [1][ list_mass [1] ] = i2;
                              }
                        }
                  }
            }

            /////////////////////////////////////////
                        
            for (j1 = 1; j1 < cell_mass [J1] + 1; j1++)
            {
                  i1 = cell_content [J1][j1];

                  n1 = part_key [i1];
                  
                  for (j2 = j1 + 1; j2 < cell_mass [J1] + 1; j2++)
                  {
                        i2 = cell_content [J1][j2];

                        n2 = part_key [i2];

                        ///////////////////////////////////////////////////////
                        
                        CC_vector ( i1, i2, xx );
                        
                        r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

                        ///////////////////////////////////////////////////////

                        if ( r2 <= R2_LIST [0] )
                        {
                              list_mass [0] += 2;
                                                
                              list_content [0][ list_mass [0] - 1 ] = i1; 
                                                
                              list_content [0][ list_mass [0] ] = i2;
                        }

                        ///////////////////////////////////////////////////////

                        else if ( r2 <= R2_LIST [1] && n1 != 1 && n1 != 2 && n2 != 1 && n2 != 2 )
                        {
                              list_mass [1] += 2;
                                                
                              list_content [1][ list_mass [1] - 1 ] = i1; 
                                                
                              list_content [1][ list_mass [1] ] = i2;
                        }
                  }
            }
      }

      //////////////////////////////////////////

      if ( klok == 0 )
      {
            f1 = fopen ( "List.dat", "a" );

            fprintf ( f1, "list_mass [0] = %d, list_mass [1] = %d\n\n", list_mass [0], list_mass [1] );

            fclose (f1);
      }

      //////////////////////////////////////////

      for (InInt = 0; InInt < 2; InInt++) 
      {
            if ( my_rank == root ) 
            { 
                  rcounts = (int *) calloc ( comm_size, sizeof(int) );
      
                  displs = (int *) calloc ( comm_size, sizeof(int) );
            
                  MPI_Gather ( &list_mass [ InInt ], 1, MPI_INT, rcounts, 1, MPI_INT, root, MPI_COMM_WORLD );
            
                  displs [0] = 1; for (k = 1; k < comm_size; k++) displs [k] = displs [k - 1] + rcounts [k - 1];
            
                  for (k = 0; k < comm_size; k++) temp_mass += rcounts [k];
            
                  temp_content = (int *) calloc ( temp_mass + 1, sizeof(int) );
            
                  MPI_Gatherv ( list_content [ InInt ] + 1, list_mass [ InInt ], MPI_INT, temp_content, rcounts, displs, MPI_INT, root, MPI_COMM_WORLD );

                  J1 = temp_mass >> 1;

                  np_busy = J1 % comm_size;

                  for (k = 0; k < comm_size; k++)
                  {
                        if ( k < np_busy ) rcounts [k] = J1 / comm_size + 1 << 1;

                        else rcounts [k] = J1 / comm_size << 1;
                  }

                  displs [0] = 1; for (k = 1; k < comm_size; k++) displs [k] = displs [k - 1] + rcounts [k - 1];
            
                  MPI_Scatter ( rcounts, 1, MPI_INT, &list_mass [ InInt ], 1, MPI_INT, root, MPI_COMM_WORLD);

                  MPI_Scatterv ( temp_content, rcounts, displs, MPI_INT, list_content [ InInt ] + 1, list_mass [ InInt ], MPI_INT, root, MPI_COMM_WORLD);
            
                  temp_mass = 0;
                  
                  free ( rcounts ); rcounts = NULL;
                  
                  free ( displs ); displs = NULL; 
            
                  free ( temp_content ); temp_content = NULL;
            }

            else
            {
                  MPI_Gather ( &list_mass [ InInt ], 1, MPI_INT, rcounts, 1, MPI_INT, root, MPI_COMM_WORLD );

                  MPI_Gatherv ( list_content [ InInt ] + 1, list_mass [ InInt ], MPI_INT, temp_content, rcounts, displs, MPI_INT, root, MPI_COMM_WORLD );

                  MPI_Scatter ( rcounts, 1, MPI_INT, &list_mass [ InInt ], 1, MPI_INT, root, MPI_COMM_WORLD);

                  MPI_Scatterv ( temp_content, rcounts, displs, MPI_INT, list_content [ InInt ] + 1, list_mass [ InInt ], MPI_INT, root, MPI_COMM_WORLD);
            }
      }

      ///////////////////////////////////////////////////////

      for (k = 1; k < hb_N + 1; k++) hb_status_proxy [k] = -1;
            
      ///////////////////////////////////////////////////////
            
      for (k = 1; k < list_mass [0] + 1; k += 2)
      {
            i1 = list_content [0][k];

            n1 = part_key [i1];

            if ( n1 > 2 ) continue;

            ///////////////////////////////////////////////////////

            i2 = list_content [0][k + 1];

            n2 = part_key [i2];

            if ( n2 > 2 ) continue;

            ///////////////////////////////////////////////////////
                  
            if (i1 < i2) J1 = (i1 - 1) * NTP_atom - (i1 - 1) * i1 / 2 + i2;

            else J1 = (i2 - 1) * NTP_atom - (i2 - 1) * i2 / 2 + i1;
                  
            for (j1 = 1; j1 < HB_PAIR_N [J1] + 1; j1++) hb_status_proxy [ HB_PAIR [J1][j1] ] = 0;
      }
            
      ///////////////////////////////////////////////////////

      MPI_Allreduce ( hb_status_proxy, hb_status, hb_N + 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

      ///////////////////////////////////////////////////////

      for (k = 1; k < HB_NT_atom + 1; k++) HB_EXCESS [k] = HB_ATOM_N [k] - VALENCE [k];

      for (k = 1; k < hb_N + 1; k++)
      {
            if ( hb_status [k] == 0 ) continue;
            
            for (i1 = 1; i1 < ATOM_HB_N [k] + 1; i1++) HB_EXCESS [ ATOM_HB [k][i1] ] -= 1;
      }
}*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*void SR_interactions_lists ( void )
{
      int i, j, k, i1, i2, j1;

      double xx [3], r2, tmp;
      
      ///////////////////////////////////////////////////////
      
      for (k = 1; k < list_mass [0] + 1; k += 2)
      {
            i = list_content [0][k];

            j = list_content [0][k + 1];

            ///////////////////////////////////////////////////////
                  
            i1 = maxi_key [i];

            i2 = maxi_key [j];

            if (i1 < i2) j1 = (i1 - 1) * na - (i1 - 1) * i1 / 2 + i2;

            else j1 = (i2 - 1) * na - (i2 - 1) * i2 / 2 + i1;

            /////////////////////////////////////////////////////////////////

            xx [0] = x [i] - x [j]; 
                  
            xx [1] = y [i] - y [j]; 
                  
            xx [2] = z [i] - z [j]; 
      
            half_shift (xx);

            r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

            /////////////////////////////////////////////////////////////////

            tmp = LJ_FORCE ( r2, j1 );

            /////////////////////////////////////////////////////////////////

            i1 = part_key [i];

            i2 = part_key [j];

            if ( i1 != 1 && i1 != 2 && i2 != 1 && i2 != 2 ) tmp += COULOMB_TRUNCATED_FORCE ( r2, j1 );
      
            /////////////////////////////////////////////////////////////////
      
            xx [0] = xx [0] * tmp; 
      
            xx [1] = xx [1] * tmp; 
      
            xx [2] = xx [2] * tmp;

            /////////////////////////////////////////////////////////////////
      
            fx [i] += xx [0]; fx [j] -= xx [0]; 
                  
            fy [i] += xx [1]; fy [j] -= xx [1]; 
                  
            fz [i] += xx [2]; fz [j] -= xx [2];
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void LR_interactions_lists ( void )
{
      int i, j, k, i1, i2, j1;

      double xx [3], r2, tmp;
      
      ///////////////////////////////////////////////////////
      
      for (k = 1; k < list_mass [1] + 1; k += 2)
      {
            i = list_content [1][k];

            j = list_content [1][k + 1];

            ///////////////////////////////////////////////////////
                  
            i1 = maxi_key [i];

            i2 = maxi_key [j];

            if (i1 < i2) j1 = (i1 - 1) * na - (i1 - 1) * i1 / 2 + i2;

            else j1 = (i2 - 1) * na - (i2 - 1) * i2 / 2 + i1;

            /////////////////////////////////////////////////////////////////

            xx [0] = x [i] - x [j]; 
                  
            xx [1] = y [i] - y [j]; 
                  
            xx [2] = z [i] - z [j]; 
      
            half_shift (xx);

            r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

            /////////////////////////////////////////////////////////////////

            tmp = COULOMB_TRUNCATED_FORCE ( r2, j1 );
      
            /////////////////////////////////////////////////////////////////
      
            xx [0] = xx [0] * tmp; 
      
            xx [1] = xx [1] * tmp; 
      
            xx [2] = xx [2] * tmp;

            /////////////////////////////////////////////////////////////////
      
            flx [i] += xx [0]; flx [j] -= xx [0]; 
                  
            fly [i] += xx [1]; fly [j] -= xx [1]; 
                  
            flz [i] += xx [2]; flz [j] -= xx [2];
      }
}*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*void polymer_constraints ( void )
{
      int i, j, k, entry;
      
      /////////////////////////////////////////////////////////////////

      for (entry = 1; entry < NBP_atom + 1; entry += 2)
      {
            i = atom_key [ entry ][0];

            j = atom_key [ entry ][1];

            /////////////////////////////////////////////////////////////////

            switch ( amino_key [ entry ] )
            {
                  case 'A': 

                        bead_bond_force ( i, j, 10.0, 4.8515 ); //r3
                        
                        //E_BOND_B += bead_bond_energy ( i, j, 10.0, 4.8515 );

                        if ( entry < NBP_atom ) 
                        {
                              k = atom_key [ entry + 1 ][0];

                              bead_valence_force ( j, i, k, 5.0, 1.9259 ); //alpha1
                              
                              //E_VALENCE_B += bead_valence_energy ( j, i, k, 5.0, 1.9259 );
                        }
                        
                        if ( entry > 1 ) 
                        {
                              k = atom_key [ entry - 1 ][0];

                              bead_valence_force ( j, i, k, 5.0, 1.7029 ); //alpha2
                              
                              //E_VALENCE_B += bead_valence_energy ( j, i, k, 5.0, 1.7029 );
                        }

                        break;

                  case 'C': 

                        bead_bond_force ( i, j, 10.0, 4.2738 ); //r3
                        
                        //E_BOND_B += bead_bond_energy ( i, j, 10.0, 4.2738 );

                        if ( entry < NBP_atom ) 
                        {
                              k = atom_key [ entry + 1 ][0];

                              bead_valence_force ( j, i, k, 5.0, 1.9655 ); //alpha1
                              
                              //E_VALENCE_B += bead_valence_energy ( j, i, k, 5.0, 1.9655 );
                        }
                        
                        if ( entry > 1 ) 
                        {
                              k = atom_key [ entry - 1 ][0];

                              bead_valence_force ( j, i, k, 5.0, 1.5803 ); //alpha2
                              
                              //E_VALENCE_B += bead_valence_energy ( j, i, k, 5.0, 1.5803 );
                        }

                        break;

                  case 'G': 

                        bead_bond_force ( i, j, 10.0, 4.9659 ); //r3
                        
                        //E_BOND_B += bead_bond_energy ( i, j, 10.0, 4.9659 );

                        if ( entry < NBP_atom ) 
                        {
                              k = atom_key [ entry + 1 ][0];

                              bead_valence_force ( j, i, k, 5.0, 1.9150 ); //alpha1
                              
                              //E_VALENCE_B += bead_valence_energy ( j, i, k, 5.0, 1.9150 );
                        }
                        
                        if ( entry > 1 ) 
                        {
                              k = atom_key [ entry - 1 ][0];

                              bead_valence_force ( j, i, k, 5.0, 1.7690 ); //alpha2
                              
                              //E_VALENCE_B += bead_valence_energy ( j, i, k, 5.0, 1.7690 );
                        }

                        break;

                  case 'U': 

                        bead_bond_force ( i, j, 10.0, 4.2733 ); //r3
                        
                        //E_BOND_B += bead_bond_energy ( i, j, 10.0, 4.2733 );

                        if ( entry < NBP_atom ) 
                        {
                              k = atom_key [ entry + 1 ][0];

                              bead_valence_force ( j, i, k, 5.0, 1.9663 ); //alpha1
                              
                              //E_VALENCE_B += bead_valence_energy ( j, i, k, 5.0, 1.9663 );
                        }
                        
                        if ( entry > 1 ) 
                        {
                              k = atom_key [ entry - 1 ][0];

                              bead_valence_force ( j, i, k, 5.0, 1.5735 ); //alpha2
                              
                              //E_VALENCE_B += bead_valence_energy ( j, i, k, 5.0, 1.5735 );
                        }

                        break;
                  
                  default: 
                        
                        break;
            }
      }
            
      /////////////////////////////////////////////////////////////////

      for (entry = 2; entry < NBP_atom + 1; entry += 2)
      {
            i = atom_key [ entry ][0];

            j = atom_key [ entry - 1 ][0];

            bead_bond_force ( i, j, 64.0, 3.8157 ); //r1 
            
            //E_BOND_B += bead_bond_energy ( i, j, 64.0, 3.8157 );

            if ( entry > 2 )
            {
                  k = atom_key [ entry - 2 ][0];
            
                  bead_valence_force ( i, j, k, 20.0, 1.4440 ); //beta2 
                  
                  //E_VALENCE_B += bead_valence_energy ( i, j, k, 20.0, 1.4440 );
            }
      }

      /////////////////////////////////////////////////////////////////

      for (entry = 3; entry < NBP_atom + 1; entry += 2)
      {
            i = atom_key [ entry ][0];

            j = atom_key [ entry - 1 ][0];

            k = atom_key [ entry - 2 ][0];

            bead_bond_force ( i, j, 23.0, 4.6010 ); //r2 
            
            //E_BOND_B += bead_bond_energy ( i, j, 23.0, 4.6010 );

            bead_valence_force ( i, j, k, 20.0, 1.5256 ); //beta1 
            
            //E_VALENCE_B += bead_valence_energy ( i, j, k, 20.0, 1.5256 );
      }
}*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*void stacking_interactions ( void )
{
      int i, j, k, k_10, k00, k01, k10, k11, k12, k20, k21, k22, k30;
      
      double r, theta1, theta2, psi, psi1, psi2;

      double psi0, psi10, psi20, r2;
      
      ///////////////////////////////////////////////////////

      for (k = 1; k < s3_N + 1; k++)
      {
            i = st_i [k];

            j = st_j [k];

            /////////////////////////////////////////////////////////////////
                  
            k10 = atom_key [i][0];

            k11 = atom_key [i][1];

            k12 = atom_key [i + 1][0];
            
            /////////////////////////////////////////////////////////////////
                  
            k20 = atom_key [j][0];

            k21 = atom_key [j][1];
                        
            k22 = atom_key [j + 1][0];
            
            /////////////////////////////////////////////////////////////////
            
            r = CC_distance ( k11, k21 );
                        
            theta1 = valence_angle ( k10, k11, k21 );
            
            theta2 = valence_angle ( k20, k21, k11 );

            psi = dihedral_angle ( k10, k11, k21, k20 );

            psi1 = dihedral_angle ( k21, k11, k10, k12 );
            
            psi2 = dihedral_angle ( k11, k21, k20, k22 );
                        
            ///////////////////////////////////////////////////////

            if ( fabs ( psi - st_psi [k] ) <= pi ) psi0 = st_psi [k];
                  
            else if ( psi - st_psi [k] > pi ) psi0 = st_psi [k] + 2.0 * pi;

            else psi0 = st_psi [k] - 2.0 * pi;

            ///////////////////////////////////////////////////////

            if ( fabs ( psi1 - st_psi1 [k] ) <= pi ) psi10 = st_psi1 [k];
            
            else if ( psi1 - st_psi1 [k] > pi ) psi10 = st_psi1 [k] + 2.0 * pi;

            else psi10 = st_psi1 [k] - 2.0 * pi;

            ///////////////////////////////////////////////////////

            if ( fabs ( psi2 - st_psi2 [k] ) <= pi ) psi20 = st_psi2 [k];
                  
            else if ( psi2 - st_psi2 [k] > pi ) psi20 = st_psi2 [k] + 2.0 * pi;

            else psi20 = st_psi2 [k] - 2.0 * pi;

            ///////////////////////////////////////////////////////

            r -= st_r [k];

            if ( r < 10.0 )
            {
                  i = ATOM_ST [k][0];

                  if ( ST_EXCESS [i] == 1 )
                  {
                        k00 = s3_N + (int) ceil (0.5 * i) - 1;

                        st_status [ k00 ] = 1;
                  }

                  ///////////////////////////////////////////////////////

                  i = ATOM_ST [k][1];

                  if ( ST_EXCESS [i] == 1 )
                  {
                        k00 = s3_N + (int) ceil (0.5 * i) - 1;

                        st_status [ k00 ] = 1;
                  }
            }
            
            ///////////////////////////////////////////////////////
            
            r2 = 1.0 + 5.00 * r * r;

            r2 += 1.50 * ( theta1 - st_theta1 [k] ) * ( theta1 - st_theta1 [k] );

            r2 += 1.50 * ( theta2 - st_theta2 [k] ) * ( theta2 - st_theta2 [k] );

            r2 += 0.15 * ( psi - psi0 ) * ( psi - psi0 );

            r2 += 0.15 * ( psi1 - psi10 ) * ( psi1 - psi10 );

            r2 += 0.15 * ( psi2 - psi20 ) * ( psi2 - psi20 );
            
            ///////////////////////////////////////////////////////

            st_energy [k] = st_E [k] / r2;

            /////////////////////////////////////////////////////////////////

            r2 = - st_energy [k] / r2;

            /////////////////////////////////////////////////////////////////
                        
            stacking_bond_force ( k11, k21, 5.00 * r2, st_r [k] );

            stacking_valence_force ( k10, k11, k21, 1.50 * r2, st_theta1 [k] );

            stacking_valence_force ( k20, k21, k11, 1.50 * r2, st_theta2 [k] );
            
            stacking_dihedral_force ( k10, k11, k21, k20, 0.15 * r2, psi0 );

            stacking_dihedral_force ( k21, k11, k10, k12, 0.15 * r2, psi10 );

            stacking_dihedral_force ( k11, k21, k20, k22, 0.15 * r2, psi20 );

            ///////////////////////////////////////////////////////
                  
            E3_STACK += st_energy [k];
      }
            
      ///////////////////////////////////////////////////////
      
      for (k = s3_N + 1; k < st_N + 1; k++)
      {
            if ( st_status [k] ) { st_status [k] = 0; st_energy [k] = 0; continue; }

            ///////////////////////////////////////////////////////
            
            i = st_i [k];
                  
            ///////////////////////////////////////////////////////
            
            k_10 = atom_key [i - 1][0]; 
            
            k00 = atom_key [i][0]; 

            k01 = atom_key [i][1];

            k10 = atom_key [i + 1][0]; 

            k20 = atom_key [i + 2][0]; 

            k21 = atom_key [i + 2][1]; 

            k30 = atom_key [i + 3][0];
            
            ///////////////////////////////////////////////////////
            
            r = CC_distance ( k01, k21 );

            psi1 = dihedral_angle ( k_10, k00, k10, k20 );
            
            psi2 = dihedral_angle ( k30, k20, k10, k00 );
            
            ///////////////////////////////////////////////////////

            if ( fabs ( psi1 + 2.58684 ) <= pi ) psi10 = - 2.58684;
                  
            else if ( psi1 + 2.58684 > pi ) psi10 = - 2.58684 + 2.0 * pi;

            else psi10 = - 2.58684 - 2.0 * pi;

            ///////////////////////////////////////////////////////

            if ( fabs ( psi2 - 3.07135 ) <= pi ) psi20 = 3.07135;
                  
            else if ( psi2 - 3.07135 > pi ) psi20 = 3.07135 + 2.0 * pi;

            else psi20 = 3.07135 - 2.0 * pi;

            ///////////////////////////////////////////////////////

            r2 = 1.0 + 1.40 * ( r - st_r [k] ) * ( r - st_r [k] );

            r2 += 4.00 * ( psi1 - psi10 ) * ( psi1 - psi10 );
            
            r2 += 4.00 * ( psi2 - psi20 ) * ( psi2 - psi20 );

            ///////////////////////////////////////////////////////
            
            st_energy [k] = st_E [k] / r2;

            /////////////////////////////////////////////////////////////////

            r2 = - st_energy [k] / r2;

            /////////////////////////////////////////////////////////////////
                  
            stacking_bond_force ( k01, k21, 1.40 * r2, st_r [k] );

            stacking_dihedral_force ( k_10, k00, k10, k20, 4.00 * r2, psi10 );
      
            stacking_dihedral_force ( k30, k20, k10, k00, 4.00 * r2, psi20 );

            /////////////////////////////////////////////////////////////////
            
            E2_STACK += st_energy [k];
      }
}*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

