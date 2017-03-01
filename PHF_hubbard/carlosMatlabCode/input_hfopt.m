%
%  input_hfopt :
%
%
%    Input file to be used in performing (symmetry-projected)
%    Hartree-Fock calculations on the 1D/2D Hubbard Hamiltonian
%    with periodic boundary conditions.
%

%  make sure that the path of the code is correct

addpath ('/shared.scratch/caj5/code/');
addpath ('/shared.scratch/caj5/code/minFunc/');


%  Hubbard Hamiltonian parameters

%    [ the hopping interaction t is set to 1 ]

%    U  - on-site repulsion strength
%    N  - number of electrons in the system

%    nx - number of unit cells in x-direction
%    ny - number of unit cells in y-direction

U  = 4.0;
N  = 8;

nx = 4;
ny = 2;


%  type of SCF to use

%    iscfty - type of SCF to use
%               iscfty  =  0,  Hartree-Fock
%                       =  1,  spin-projected Hartree-Fock
%                       =  2,  LM projected Hartree-Fock
%                       =  3,  LM + spin-projected Hartree-Fock

%    igues  - type of guess to use
%               igues   =  0,  prepare random initial guess
%                       =  1,  read initial guess from file

iscfty = 0;
igues  = 0;


%  LM projection parameters

%    ikx - momentum index in x-direction
%    iky - momentum index in y-direction

ikx = 0;
iky = 0;


%  spin projection parameters

%    imult - multiplicity of state to recover **
%    Sgrd  - number of grid points in spin projection

%      **  s  =   0   =>   imult  =  1
%             = 1/2               =  2
%             =   1               =  3,  etc.

imult = 1;
Sgrd  = [ 9 16 9 ];


%  SCF convergence

%    maxit  - maximum number of SCF iterations
%    maxerr - convergence criterion
%             ( in terms of norm of the gradient )

maxit  = 2000;
maxerr = 1.0e-6;


%  output files

%    outstr  - string to use in producing output files **
%    loadmat - name of MAT file to be used in reading guess ***
%    savemat - whether to save optimized state into MAT file

%      **  read from the environment variable HUBOUT
%      *** read from the environment variable HUBLOAD

outstr  = getenv ('HUBOUT');
loadmat = getenv ('HUBLOAD');
savemat = 1;


%  run calculation

  %if ( iscfty ~= 0 )
  %  matlabpool local 8
  %end


run_hfopt


  %if ( iscfty ~= 0 )
  %  matlabpool close
  %end


quit force

