#include "stdint.h"
#define __STDC_FORMAT_MACROS
#include "inttypes.h"

#ifndef PRId64
#define PRId64 "ld"
#endif
#if !defined(LAMMPS_SMALLSMALL) && !defined(LAMMPS_BIGBIG) && !defined(LAMMPS_SMALLBIG)
#define LAMMPS_SMALLBIG
#endif

#if defined(LAMMPS_SMALLBIG)
typedef int tagint;
typedef int64_t bigint;
#define BIGINT_FORMAT "%" PRId64
#elif defined(LAMMPS_SMALLSMALL)
typedef int tagint;
typedef int bigint;
#define BIGINT_FORMAT "%d"
#else /* LAMMPS_BIGBIG */
typedef int64_t tagint;
typedef int64_t bigint;
#define BIGINT_FORMAT "%" PRId64
#endif

inline bool fexists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

int n;
bigint ntimestep,natoms;
int size_one,nchunk,triclinic;
double xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz;

double *DATA = NULL;

char *lineS = NULL;
char *lineD = NULL;

double ZLO,ZHI;
int zloflag, zhiflag;

int boundary[3][2];
int bin_changed = 0;

int flag_VD,flag_GR,flag_VACF,flag_MSD;

double vmin,vmax,occupation;
int iV,start,stop,msd,vacf;
int nstep = 0;
int N,M,output,step_size,binflag,equalbin;
int run,plot;
int step = 0;
double dt;
int post;

FILE *screen = stdout;  
FILE *ftemp = NULL;

char *filedatgr,*veldat,*msddat,*vacfdat,*filedat,*vardat;
int *vel = NULL;

double *variance = NULL;
double *vector_msd = NULL;
double *vector_vacf = NULL;
double *xoriginal = NULL;
double *voriginal = NULL;

double *varT = NULL;
double *varT_orth = NULL;
double *varT_paral = NULL;

double *varT1 = NULL;
double *varT_orth1 = NULL;
double *varT_paral1 = NULL;

double *prop = NULL;
double *prop_U = NULL;

clock_t tbins2;
double tbin2=0.0;

double *bin = NULL;
double *Nocpt = NULL;
double *Nsize = NULL;
double *time_div = NULL;
double *den = NULL;	
double *bulk_vel = NULL;
double *zero_bulk = NULL;

double *bin_w = NULL;
double *bin_vol = NULL;
double *bin_mean_z = NULL;

int maxbuf = 0;
double *buf = NULL;
//double *Tbuf = NULL;
//double *tmp = NULL;
//double *buf1 = NULL;

const double Na = 6.02214086e23; // 1/mol
const double kb = 1.38064852e-23; // J/K
const double R = Na * kb;
const double pi = 3.141592653589793;


