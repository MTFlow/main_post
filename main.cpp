#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sys/stat.h>
#include "stdlib.h"
#include "math.h"
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <list>
#include <time.h>
//#include <mpi.h>
#include <omp.h>
#include <ctime>
#include <cmath>
#include <string>
//#include <tuple/tuple.hpp>
//#include "gnuplot-iostream.h"
#include "main.h"
#include "func.h"



void help() {
	printf(
		"\nList of command line options which can be specified\n\n"
		"start(>=0)\t: start at this timestep (-start) (default: %d)\n"
		"stop(>=0)\t: stop at this timestep (-stop) (default: %d)\n"
		"N(>0)\t\t: number of bins (-N) (default: %d)\n"
		"output every(>=0)\t: when to print timestep to screen (-output_every) (default: %d)\n"
		"stepsize(>=0)\t: timestep between data saved in file (-stepsize) (default: %d)\n"
		"binflag(0|1)\t: find best distribution of bins with minimum occupation per bin (-binflag) (default: %d)\n"
		"occupation(0:1)\t: occupation per bin (-occupation) (default: %f)\n"
		"iV(>=0)\t\t: number of bins for velocity distribution (-iV) (default: %d)\n"
		"vmin\t\t: minimum velocity for distribution (-vmin) (default: %f)\n"
		"vmax\t\t: maximum velocity for distribution (-vmax) (default: %f)\n"
		"post(0|1)\t: print task timing breakdown table (-post) (default: %d)\n"
		"zlo\t\t: define zlo (-zlo) (default zloflag: %d)\n"
		"zhi\t\t: define zhi (-zhi) (default zhiflag: %d)\n"
		"Bin spacing(0|1): bins are equally spaced (-equalbin) (default: %d)\n"
		"VD(0|1)\t\t: write velocity distribution (-flag_VD) (default: %d)\n"
		"GR(0|1)\t\t: write radial distribution function (-flag_GR) (default: %d)\n"
		"VACF(0|1)\t: write velocity autocorrelation function (-flag_VACF) (default: %d)\n"
		"MSD(0|1)\t: write mean square displacement (-flag_MSD) (default: %d)\n\n",
		start,stop,N,output_every,step_size,binflag,occupation,iV,vmin,vmax,post,zloflag,zhiflag,equalbin,flag_VD,flag_GR,flag_VACF,flag_MSD );
};

void input_check(FILE *fp) {
	fprintf(screen,"\n%s\n",lineD);
	fprintf(screen,"CHECK INPUT ARGUMENTS\n");

	int I = 0;	
	int NCHUNK = 0;	
	int TIME[2];
	int flag = 0;

	while (I < 2) {
		fread(&ntimestep,sizeof(bigint),1,fp);
		fread(&natoms,sizeof(bigint),1,fp);		
		fread(&triclinic,sizeof(int),1,fp);		
		fread(&boundary[0][0],6*sizeof(int),1,fp);	
		fread(&xlo,sizeof(double),1,fp);
		fread(&xhi,sizeof(double),1,fp);
		fread(&ylo,sizeof(double),1,fp);
		fread(&yhi,sizeof(double),1,fp);
		fread(&zlo,sizeof(double),1,fp);
		fread(&zhi,sizeof(double),1,fp);
		if (triclinic) {
			fread(&xy,sizeof(double),1,fp);
			fread(&xz,sizeof(double),1,fp);
			fread(&yz,sizeof(double),1,fp);
		}
		fread(&size_one,sizeof(int),1,fp);
		fread(&nchunk,sizeof(int),1,fp);

		for (int i = 0; i < nchunk; i++) {
			fread(&n,sizeof(int),1,fp);		
			if (buf) delete [] buf;
			buf = new double[n];
			fread(buf,sizeof(double),n,fp);

			NCHUNK += n;
		}

		mass = buf[0] / (Na * 1000.0); //mass [kg]

		TIME[I] = (int)ntimestep;
		I += 1;
	}

	char boundstr[9];
  int m = 0;
  for (int idim = 0; idim < 3; idim++) {
		for (int iside = 0; iside < 2; iside++) {
		  if (boundary[idim][iside] == 0) boundstr[m++] = 'p';
	  	else if (boundary[idim][iside] == 1) boundstr[m++] = 'f';
	  	else if (boundary[idim][iside] == 2) boundstr[m++] = 's';
	  	else if (boundary[idim][iside] == 3) boundstr[m++] = 'm';
		}
		boundstr[m++] = ' ';
  }
  boundstr[8] = '\0';

	if (start < TIME[0] || (start%(TIME[1]-TIME[0]))!=0 ) {
		fprintf(screen,"WARNING >> Changed start value from %d to %d!\n",start,TIME[0]);
		start = TIME[0];
		flag = 1;
	}

	if (step_size < (TIME[1]-TIME[0]) || (step_size%(TIME[1]-TIME[0]))!=0) {
		fprintf(screen,"WARNING >> Changed stepsize value from %d to %d!\n",step_size,(TIME[1]-TIME[0]));
		step_size = (TIME[1] - TIME[0]);
		flag = 1;
	};

	if (stop < start) {
		fprintf(screen,"PROGRAM TERMINATED >> Start value %d larger than stop value %d\n",start,stop);
		fprintf(screen,"CHECK INPUT ARGUMENTS\n");
		fprintf(screen,"%s\n\n",lineD);
		fclose(fp);	
		exit(1);
	}

	if (stop%step_size!=0) {
		fprintf(screen,"WARNING >> Stop value %d NOT multiple of stepsize %d\n",stop,step_size);
		int stop_new = stop - stop%step_size;
		fprintf(screen,"           Changed stop value %d to %d\n",stop,stop_new);
		stop = stop_new;
		flag = 1;
	}

	
//	int BYTES = (stop/(TIME[1]-TIME[0]))*(2*sizeof(bigint)+(9+nchunk)*sizeof(int)+(0.5*NCHUNK+6)*sizeof(double));
//	fseek(fp,BYTES,SEEK_SET);

//	bigint stop_check;
//	fread(&stop_check,sizeof(bigint),1,fp);

//	fprintf(screen,"check:" BIGINT_FORMAT " BYTES: %d\n",stop_check,BYTES);

//	if (stop_check < stop) {
//		fprintf(screen,"WARNING >> Stop value is larger than max # timesteps in data file\n");
//		flag = 1;
//	};

	if (!flag)	fprintf(screen,"OKE\n");
//	fprintf(screen,"CHECK INPUT ARGUMENTS\n");
	fprintf(screen,"%s\n",lineS);
	fprintf(screen,"ATOM ATTRIBUTES\n%d\n",size_one);
	fprintf(screen,"1:Mass | 2:Mol | 3:Id | 4:Type | \n5:x | 6:y | 7:z | 8:vx | 9:vy | 10:vz | \n11:Ke | 12:Pe | \n13:Sxx | 14:Syy | 15:Szz | 16:Sxy | 17:Sxz | 18:Syz\n");
	fprintf(screen,"%s\n",lineS);
	fprintf(screen,"BOX BOUNDARIES\n");
	fprintf(screen,"X: %c%c | Y: %c%c | Z: %c%c\n",boundstr[0],boundstr[1],boundstr[3],boundstr[4],boundstr[6],boundstr[7]);
	fprintf(screen,"%s\n",lineS);
	fprintf(screen,"BOX DIMENSIONS\n");
	fprintf(screen,"Xlo: %-10.2f | %-5s Xhi: %-10.2f | %-5s dX: %.2f\n",xlo,"",xhi,"",xhi-xlo);
	fprintf(screen,"Ylo: %-10.2f | %-5s Yhi: %-10.2f | %-5s dY: %.2f\n",ylo,"",yhi,"",yhi-ylo);
	fprintf(screen,"Zlo: %-10.2f | %-5s Zhi: %-10.2f | %-5s dZ: %.2f\n",zlo,"",zhi,"",zhi-zlo);
	fprintf(screen,"%s\n",lineS);
	fprintf(screen,"NUMBER OF ATOMS\n");
	fprintf(screen, BIGINT_FORMAT "\n",natoms);
	fprintf(screen,"%s\n",lineS);
	fprintf(screen,"NUMBER OF PROCESSORS USED\n");
	fprintf(screen,"%d\n",nchunk);
	fprintf(screen,"%s\n",lineS);
	fprintf(screen,"STEPS BETWEEN OUTPUT DATA\n");
	fprintf(screen,"%d\n",TIME[1]-TIME[0]);

	fprintf(screen,"%s\n",lineD);

}

inline double fh(double sec) {
	double scd = sec/(double)CLOCKS_PER_SEC;
	return floor(scd/3600.0);
}

inline double fm(double sec) {
	double scd = sec/(double)CLOCKS_PER_SEC;
	return floor((scd - 3600*floor(scd/3600.0))/60.0);
}

inline double fs(double sec) {
	double scd = sec/(double)CLOCKS_PER_SEC;
	return scd - 3600.0*floor(scd/3600.0) - 60.0*floor((scd - 3600*floor(scd/3600.0))/60.0);
}



int main(int narg, char **arg)
{

//	Gnuplot gp;

//	fprintf(screen,"Working directory: %s\n");

//set initial values
	start = 0;
	stop = 0;
	N = 1;
	output_every = 100;
	step_size = 500;
	binflag = 1;
	occupation = 1.0;
	dt = 1.0;
	iV = 500;
	vmin = -0.01;
	vmax = 0.01;
	post = 1;
	M = 33; //SIZE OF PROP ARRAY;
  zloflag = 0;
	zhiflag = 0;
	equalbin = 0;
	flag_VD = 0;
	flag_GR = 0;
	flag_VACF = 0;
	flag_MSD = 0;
	flag_VAR = 0;
	flag_REGION = 0;
	flag_output_every = 0;
	plot = 0;
	dregionLO = 0.0;
	dregionHI = 0.0;

	lineD = new char[91];
	strcpy(lineD,"==========================================================================================");
	lineS = new char[91];
	strcpy(lineS,"------------------------------------------------------------------------------------------");


	if (strcmp(arg[1],"-help") == 0 || strcmp(arg[1],"-h") == 0) {
		help();
		return 1;
	};

	int iarg = 2;
	while (iarg < narg) {
		if ( strcmp(arg[iarg],"-help") == 0 || strcmp(arg[iarg],"-h") == 0) {
			help();
			return 1;
		} else if ( strcmp(arg[iarg],"-start") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			double starti = atof(arg[iarg+1]);
			if (starti!=floor(starti)) {
				fprintf(screen,"WARNING: Start value changed to %d!\n",(int)floor(starti));
 				start = (int)floor(starti);
			} else {
				start = (int)starti;
			};
		} else if ( strcmp(arg[iarg],"-stop") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			double stopi = atof(arg[iarg+1]);
			if (stopi!=floor(stopi)) {
				fprintf(screen,"WARNING: Stop value changed to %d!\n",(int)floor(stopi));
 				stop = (int)floor(stopi);
			} else {
				stop = (int)stopi;
			};
		} else if ( strcmp(arg[iarg],"-N") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			double Ni = atof(arg[iarg+1]);
			if (Ni < 0) {fprintf(screen,"N should be larger than '0'\n "); return 1;}
			if (Ni!=floor(Ni)) {
				fprintf(screen,"WARNING: N value changed to %d!\n",(int)floor(Ni));
 				N = (int)floor(Ni);
			} else {
				N = (int)Ni;
			};
		} else if ( strcmp(arg[iarg],"-output_every") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			double output_everyi = atof(arg[iarg+1]);
			if (output_everyi!=floor(output_everyi)) {
				fprintf(screen,"WARNING: Output value changed to %d!\n",(int)floor(output_everyi));
 				output_every = (int)floor(output_everyi);
			} else {
				output_every = (int)output_everyi;
			};
			flag_output_every = 1;
		} else if ( strcmp(arg[iarg],"-stepsize") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			double step_sizei = atof(arg[iarg+1]);
			if (step_sizei!=floor(step_sizei)) {
				fprintf(screen,"WARNING: step_size value changed to %d!\n",(int)floor(step_sizei));
 				step_size = (int)floor(step_sizei);
			} else {
				step_size = (int)step_sizei;
			};
		} else if ( strcmp(arg[iarg],"-binflag") == 0 ) {
				if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
				binflag = atoi(arg[iarg+1]);
				if (binflag > 1) {fprintf(screen,"Binflag should be '0' or '1'\n "); return 1;}
		} else if ( strcmp(arg[iarg],"-occupation") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			occupation = atof(arg[iarg+1]);
			if (occupation > 1 | occupation < 0) {fprintf(screen,"Occupation should be between '0' and '1'\n "); return 1;}
		} else if ( strcmp(arg[iarg],"-iV") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			iV = atof(arg[iarg+1]);
			if (iV < 0) {fprintf(screen,"iV should be larger than '0'\n "); return 1;}
			double iVi = atof(arg[iarg+1]);
			if (iVi!=floor(iVi)) {
				fprintf(screen,"WARNING: iV value changed to %d!\n",(int)floor(iVi));
 				iV = (int)floor(iVi);
			} else {
				iV = (int)iV;
			};
		} else if ( strcmp(arg[iarg],"-vmin") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			vmin = atof(arg[iarg+1]);
		} else if ( strcmp(arg[iarg],"-vmax") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			vmax = atof(arg[iarg+1]);
		} else if ( strcmp(arg[iarg],"-dt") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			dt = atof(arg[iarg+1]);
			if (dt <= 0) {fprintf(screen,"Timestep should be larger than '0'\n "); return 1;}
		} else if ( strcmp(arg[iarg],"-post") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			post = atoi(arg[iarg+1]);
		} else if ( strcmp(arg[iarg],"-zlo") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			zloflag = 1;
			ZLO = atoi(arg[iarg+1]);
		} else if ( strcmp(arg[iarg],"-zhi") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			zhiflag = 1;
			ZHI = atoi(arg[iarg+1]);
		} else if ( strcmp(arg[iarg],"-equalbin") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			equalbin = atoi(arg[iarg+1]);
		} else if ( strcmp(arg[iarg],"-flag_GR") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			flag_GR = atoi(arg[iarg+1]);
		} else if ( strcmp(arg[iarg],"-flag_VD") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			flag_VD = atoi(arg[iarg+1]);
		} else if ( strcmp(arg[iarg],"-flag_VACF") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			flag_VACF = atoi(arg[iarg+1]);
		} else if ( strcmp(arg[iarg],"-flag_MSD") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			flag_MSD = atoi(arg[iarg+1]);
		} else if ( strcmp(arg[iarg],"-flag_VAR") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			flag_VAR = atoi(arg[iarg+1]);
		} else if ( strcmp(arg[iarg],"-plot") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			plot = atoi(arg[iarg+1]);
		} else if ( strcmp(arg[iarg],"-region") == 0 ) {
			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
			flag_REGION = 1;
			if (arg[iarg+1]!=NULL) dregionLO = atof(arg[iarg+1]);
			if (arg[iarg+2]!=NULL) dregionHI = atof(arg[iarg+2]);
			iarg += 1;
//		} else if ( strcmp(arg[iarg],"-no_ocp") == 0 ) {
//			if (iarg+2 > narg) {fprintf(screen,"Missing argument value of '%s'\n",arg[iarg]); return 1;}
//			char *str = arg[iarg+1];
//			fprintf(screen,"%lu \n",strlen(str));	
//		no_array = atoi(arg[iarg+1]);
		} else {fprintf(screen,"Illegal name: '%s'\n",arg[iarg]); return 1;};
		iarg += 2;
	}

	time_t now = time(0);
	tm *ltm = localtime(&now);
	int year = 1900 + ltm->tm_year;
	int month = 1 + ltm->tm_mon;
	int day = ltm->tm_mday;
	int hour = ltm->tm_hour;
	int minute = ltm->tm_min;
	int second = ltm->tm_sec;

	fprintf(screen,"\n%s\n",lineD);
	fprintf(screen,"Date: %02d/%02d/%04d (d/m/y)\n",day,month,year);
	fprintf(screen,"Time: %02d:%02d:%02d (hr:min:sec)",hour,minute,second);
	fprintf(screen,"\n%s\n",lineD);


	FILE *fp = fopen(arg[1],"rb");
	if (!fp) {
		fprintf(screen,"ERROR: could not open file\n");
		exit(1);	
	}

	input_check(fp);
	fseek(fp,0,SEEK_SET);


	DATA = new double[natoms*3];


//	msd = 1;
//	vacf = 1;

	clock_t tinits,ttotals;
	clock_t tbins,tbin_reads,tbin_assigns,tbin_occups,tbin_props,tbin_bins;
	clock_t trun1s,trun1_reads,trun1_assigns,trun1_props;
	clock_t twrite_grs,twrite_vels,twrite_vars,twrite_bins;


	double tinit=0.0,tbin=0.0,tbin_read=0.0,tbin_assign=0.0,tbin_occup=0.0,tbin_prop=0.0,tbin_bin=0.0;
	double trun1=0.0,trun1_read=0.0,trun1_assign=0.0,trun1_prop=0.0,twrite_gr=0.0,twrite_vel=0.0,twrite_var=0.0,twrite_bin=0.0,ttotal=0.0;

//	double totaltime;
//	starttime = clock();
//	double start_omp = omp_get_wtime();
//	int threads = 4;
//	omp_set_num_threads(threads);

	ttotals = clock();

//====INIT====//
	tinits = clock();

	double coord[7] = {50.0, 50.0, 10.0, 10.0, 40.0, 15.0, 60.0}; //x,y,z,dxyz,Ra,Rrho,Nshells

	double Nstep = ((stop-start)/step_size) + 1;

	if (!flag_output_every) output_every = (int)((stop-start)/10.0);

	fprintf(screen,"\n%s\n",lineD);
	fprintf(screen,"INPUT");
	fprintf(screen,"\n%s\n",lineS);
	fprintf(screen,"%-16s%-25s%s\n","NAME","COMMAND","VALUE");
	fprintf(screen,"%-16s%-25s%d\n","Start","-start",start);
	fprintf(screen,"%-16s%-25s%d\n","Stop","-stop",stop);
	fprintf(screen,"%-16s%-25s%d\n","Nbins","-N",N);
	fprintf(screen,"%-16s%-25s%d\n","Binflag","-binflag",binflag);
	fprintf(screen,"%-16s%-25s%d\n","Output every","-output_every",output_every);
	fprintf(screen,"%-16s%-25s%d\n","Step size","-stepsize",step_size);
	fprintf(screen,"%-16s%-25s%f\n","Timestep","-dt",dt);
	fprintf(screen,"%-16s%-25s%d\n","iV","-iV",iV);
	fprintf(screen,"%-16s%-24s%f\n","vmin","-vmin",vmin);
	fprintf(screen,"%-16s%-25s%f\n","vmax","-vmax",vmax);
	if (zloflag) fprintf(screen,"%-16s%-25s%f\n","ZLO","-zlo",ZLO);
	if (zhiflag) fprintf(screen,"%-16s%-25s%f\n","ZHI","-zhi",ZHI);
	if (flag_REGION) fprintf(screen,"%-16s%-25s%f\n","Region lower"," ",dregionLO);
	if (flag_REGION) fprintf(screen,"%-16s%-25s%f\n","Region higher"," ",dregionHI);
	fprintf(screen,"%s\n",lineD);


	double temp[N][3];
	double temp1[N][4];

//	double den[N];
	double press[N][6];
	double hflux[N][3];
//	double bulk_vel[N][3];
	double P[N];
//	double zero_bulk[N][3];
	double Pe[N];
	double enthlp[N];

	double eflux[N][3];

	double gr[(int)coord[6]];
	double Nra[(int)coord[6]];

	int sizeY = (((stop-start)/step_size)+1)+1;
	double *variance = new double[N*sizeY];

	varT = new double[N*sizeY];
	varT_orth = new double[N*sizeY];
	varT_paral = new double[N*sizeY];

	varT1 = new double[N*sizeY];
	varT_orth1 = new double[N*sizeY];
	varT_paral1 = new double[N*sizeY];

	zero_bulk = new double[3*N];
	bulk_vel = new double[3*N];

//	vel = new int[3*N*iV]; //Vx, Vy, Vz

	int indexrdf[10];

	for (int i = 0; i < 10; i++) indexrdf[i] = 0;

	for (int i1 = 0; i1 < coord[6]; i1 ++) {
		gr[i1]=0.0;
		Nra[i1]=0.0;
	}

//	for (int v1 = 0; v1 < 3*N*iV; v1++) vel[v1] = 0;

	for (int ii = 0; ii < N; ii++) {
		temp[ii][0]=0.0; temp[ii][1]=0.0; temp[ii][2]=0.0;
		temp1[ii][0]=0.0; temp1[ii][1]=0.0; temp1[ii][2]=0.0; temp1[ii][3]=0.0;
	//	den[ii]=0.0;
		press[ii][0]=0.0; press[ii][1]=0.0; press[ii][2]=0.0; 
		press[ii][3]=0.0; press[ii][4]=0.0; press[ii][5]=0.0;
		hflux[ii][0]=0.0; hflux[ii][1]=0.0; hflux[ii][2]=0.0;
		bulk_vel[ii*3]=0.0; bulk_vel[ii*3+1]=0.0; bulk_vel[ii*3+2]=0.0;
		P[ii]=0.0;
		zero_bulk[ii*3]=0.0;zero_bulk[ii*3+1]=0.0;zero_bulk[ii*3+2]=0.0;
		eflux[ii][0]=0.0; eflux[ii][1]=0.0; eflux[ii][2]=0.0;
		Pe[ii]=0.0;
		enthlp[ii]=0.0;

		for (int jj = 0; jj < sizeY; jj ++) {
			varT[ii*sizeY+jj] = 0.0;
			varT_orth[ii*sizeY+jj] = 0.0;
			varT_paral[ii*sizeY+jj] = 0.0;

			varT1[ii*sizeY+jj] = 0.0;
			varT_orth1[ii*sizeY+jj] = 0.0;
			varT_paral1[ii*sizeY+jj] = 0.0;
		}
	}

	tinit += clock() - tinits;
	double value = 1.0;
		
	// Gnuplot vectors (i.e. arrows) require four columns: (x,y,dx,dy)
//	std::vector<boost::tuple<double, double, double, double> > pts_A;
//	gp << "set terminal wxt size 1500,900\n";

//==========BIN==========//

		tbins = clock();
		double min_occ;
		int min_bin;
	
		double time_prv = 0.0;


		if (flag_REGION) {
			fnameregion = new char[50];
			sprintf(fnameregion,"region_%s_%s_%d.dat",dot2underscore(dregionLO,2),dot2underscore(dregionHI,2),N);
			if (fexists(fnameregion)) sprintf(fnameregion,"region_%s_%s_%d_%02d_%02d_%04d_%02d%02d%02d.dat",dot2underscore(dregionLO,2),dot2underscore(dregionHI,2),N,day,month,year,hour,minute,second);
			fregion = fopen(fnameregion,"w");
		};

		if (binflag) {

			int Nnew = 0;
			int once = 1;

			bin = new double[N+1];
			Nocpt = new double[N];
			time_div = new double[N];
			den = new double[N];
			bin_mean_z = new double[N];

			fprintf(screen,"\n%s\n",lineD);
			fprintf(screen,"BINFLAG");
			fprintf(screen,"\n%s\n",lineS);
			fprintf(screen,"Occupation: %.2f",occupation);
			fprintf(screen,"\n%s",lineS);
			run = 0; //binflag


			for(int i1 = 0; i1 < N; i1++ ) {
				time_div[i1] = 0.0;
				den[i1] = 0.0;
				bulk_vel[i1*3] = 0.0; bulk_vel[i1*3+1] = 0.0; bulk_vel[i1*3+2] = 0.0;
				bin_mean_z[i1] = 0.0;
			};

			while(1) {

				vel = new int[3*N*iV]; //Vx, Vy, Vz
				for (int v1 = 0; v1 < 3*N*iV; v1++) vel[v1] = 0;

				step = 0;
				fprintf(screen,"\nN: %d\nTime steps:",N);

				while(1) {

					fread(&ntimestep,sizeof(bigint),1,fp);	

					if (feof(fp) || (ntimestep > stop) ) {
						fseek(fp, 0, SEEK_SET );
						break;
					};
	
					if ((step%output_every)==0) fprintf(screen," " BIGINT_FORMAT, ntimestep);
					fflush(stdout);		
					step += step_size;

					tbin_reads = clock();		
					fread(&natoms,sizeof(bigint),1,fp);		
					fread(&triclinic,sizeof(int),1,fp);		
					fread(&boundary[0][0],6*sizeof(int),1,fp);	
      		fread(&xlo,sizeof(double),1,fp);
      		fread(&xhi,sizeof(double),1,fp);
  	  		fread(&ylo,sizeof(double),1,fp);
	    		fread(&yhi,sizeof(double),1,fp);
      		fread(&zlo,sizeof(double),1,fp);
      		fread(&zhi,sizeof(double),1,fp);
      		if (triclinic) {
						fread(&xy,sizeof(double),1,fp);
						fread(&xz,sizeof(double),1,fp);
						fread(&yz,sizeof(double),1,fp);
      		}
      		fread(&size_one,sizeof(int),1,fp);
      		fread(&nchunk,sizeof(int),1,fp);
					tbin_read += clock() - tbin_reads;
	

					if (ntimestep >= start && ntimestep <= stop) { 

						if (once) {
							if (!zloflag) ZLO = zlo;
							if (!zhiflag) ZHI = zhi;
							chunkdomain(bin,N);
							once = 0;
						};


// prop = dof, ke, ke_orth, ke_paral, Sxx, Syy, Szz, Sxy, Sxz, Syz, Jcx, Jcy, Jcz, Jvx, Jvy, Jvz

						if (prop) delete [] prop;
						prop = new double[N*M];

						for(int i1 = 0; i1 < N*M; i1++ ) prop[i1] = 0.0;

						int temp = 0;
						int tmp = 0;

						for (int i = 0; i < nchunk; i++) {
							tbin_reads = clock();
							fread(&n,sizeof(int),1,fp);
				
							if (n > maxbuf) {
								if (buf) delete [] buf;
								buf = new double[n];
								maxbuf = n;
							}

							fread(buf,sizeof(double),n,fp);
							tbin_read += clock() - tbin_reads;	
			
							for (int s = 0; s<n/size_one; s++) {

								DATA[tmp*3+s*3+0] = buf[4+size_one*s];
								DATA[tmp*3+s*3+1] = buf[5+size_one*s];
								DATA[tmp*3+s*3+2] = buf[6+size_one*s];

							};
		
							tmp += n/size_one;

							tbin_assigns = clock();
				  		assignparticle(indexrdf,bin,N,n/size_one,zero_bulk,size_one,variance,sizeY,coord,vel,0);
							tbin_assign += clock() - tbin_assigns;

//							if (plot) {
//								for (int t = 0; t < n/size_one; t++) {
//									pts_A.push_back(boost::make_tuple(buf[5+size_one*t],buf[6+size_one*t],0.01,0.01));
//								};
//							};
						}

//						if (plot) {
//							gp << "plot '-' with points\n";
//							gp.send1d(pts_A);	
//							gp.flush();
//							pts_A.clear();
//						};
							

//						tbin_occups = clock();
//						for (int l = 0; l < N; l++) timediv(l,time_div);
//						tbin_occup += clock() - tbin_occups;

						tbin_occups = clock();
						for (int l = 0; l < N; l++) {		
							timediv(l,time_div);
							density(l,den,variance,sizeY);
							bulk_velocity(l,bulk_vel);
							pe_energy(l,Pe);
							mean_z(l,bin_mean_z);
						}
						tbin_occups += clock() - tbin_occups;



					} else {

						tbin_reads = clock();	
						for (int i = 0; i < nchunk; i++) {
							fread(&n,sizeof(int),1,fp);				

							if (n > maxbuf) {
								if (buf) delete [] buf;
								buf = new double[n];
								maxbuf = n;
							}

							fread(buf,sizeof(double),n,fp);
						}
						tbin_read += clock() - tbin_reads;

					}; // end of start/stop
				}; // end of while (1) over time steps

				
				if (flag_REGION) fclose(fregion);
	
				fprintf(screen," | time: %02.0f:%02.0f:%02f",fh(tbin_assign-time_prv),fm(tbin_assign-time_prv),fs(tbin_assign-time_prv));

				time_prv = tbin_assign;

				tbin_bins = clock();

				Nnew = 0;
				int extraN = 0;
				min_occ = 1.0;
				min_bin = 0.0;
				for (int p = 0; p < N; p++) {
					Nocpt[p] = time_div[p]/Nstep;
//					fprintf(screen,"bin: %d - OCC: %f \n",p,Nocpt[p]);

					if (Nocpt[p] < min_occ) {
						min_occ = Nocpt[p];
						min_bin = p;
					};
					if (Nocpt[p] >= occupation) Nnew += 1; 
//else {
//						if (((time_div[p+1]/Nstep) >= occupation) && p < (N-1)) extraN += 1;
//					};
				};	
	
				if (time_div[N-1]/Nstep < occupation) extraN = 1;

				Nnew = floor((N-Nnew)/2) + Nnew - floor(extraN/2);
//				Nnew = floor((N-Nnew)/2) + Nnew - extraN;
//				Nnew -= extraN;

				fprintf(screen,"\nZLO:%.2f - ZHI: %.2f | minimum occupation: %E - bin: %d",zlo,zhi,min_occ,min_bin);

				if (Nnew == N || Nnew == 0) break;

				double binnew[Nnew+1];
				int B = 0;
				binnew[0] = bin[0];
//				fprintf(screen,"\n %f - %f\n",binnew[Nnew],bin[N]);

				for (int b = 0; b < N-1; b++) {
						if (Nocpt[b] >= occupation) {
							binnew[B+1] = bin[b+1];
						}	else { 
							binnew[B+1] = bin[b+2];
							b += 1;
						}	
					B += 1;
				};

				binnew[Nnew] = bin[N];

/*
				for (int b = 0; b < N; b++) {
						if (Nocpt[b] >= occupation) {
							binnew[B+1] = bin[b+1];
						}	else { 
							binnew[B+1] = bin[b+2];
							b += 1;
						}	
					B += 1;
				};
*/

//				fprintf(screen,"\n");
//				for (int t1=0; t1<N+1;t1++) {
//					fprintf(screen,"bin-%d:%f - %f\n",t1,bin[t1],Nocpt[t1]);
//				};
//				for (int t1=0; t1<Nnew+1;t1++) {
//					fprintf(screen,"binnew-%d:%f\n",t1,binnew[t1]);
//				};


				N = Nnew;
				Nocpt = new double[N];
				bin = new double[N+1];

				time_div = new double[N];
				den = new double[N];
				bulk_vel = new double[3*N];


				for (int k = 0; k < N; k++)	{
					bin[k] = binnew[k];
					time_div[k] = 0.0;
					den[k] = 0.0;
					bulk_vel[k*3] = 0.0;
					bulk_vel[k*3+1] = 0.0;
					bulk_vel[k*3+2] = 0.0;
					bin_mean_z[k] = 0.0;
				};
				bin[N] = binnew[N];

			
				tbin_bin += clock() - tbin_bins;		

			}; //end of while (N......)

		}; //end of binflag
	
		bin_w = new double[N];
		bin_vol = new double[N];

		tbin_props = clock();
		for (int t = 0; t < N; t++) {
			bin_w[t] = (bin[t]+bin[t+1]) * 1e-10 / 2; //m
			bin_vol[t] = (bin[t+1]-bin[t]) * (xhi - xlo) * (yhi - ylo) * 1e-30; //m3
			den[t] /= (bin_vol[t]*Nstep); 
			Pe[t] /= Nstep;
			bulk_vel[t*3] /= Nstep;
			bulk_vel[t*3+1] /= Nstep;
			bulk_vel[t*3+2] /= Nstep;
			bin_mean_z[t] /= Nstep;
		};

		tbin_prop += clock() - tbin_props;

		tbin += clock() - tbins;

//		fprintf(screen,"\n-------------------------------------------------------\n");
//		if (binflag) fprintf(screen,"Minimum occupation: %.2f - bin: %d\n",min_occ,min_bin);
 		fprintf(screen,"\n%s\n",lineD);


		int ll = 1;
		while (N/(1.0+pow(10.0,(double)ll-1.0)) >= 10.0) ll++;

//======write output_bin_gr======//
		if (flag_GR) {
		twrite_grs = clock();

		filedatgr = new char[ll+100];
		sprintf(filedatgr,"output_bin_gr_%d.dat",N);
		if (fexists(filedatgr)) sprintf(filedatgr,"output_bin_gr_%d_%02d_%02d_%04d_%02d%02d%02d.dat",N,day,month,year,hour,minute,second);
		FILE *fpdatgr = fopen(filedatgr,"w");

		for (int iw = 0; iw < coord[6]; iw++) {
			fprintf(fpdatgr,"%E %E\n",(iw+1)*coord[4]/coord[6],gr[iw]);
		};

		fclose(fpdatgr);
		twrite_gr = clock() - twrite_grs;
		};
//=========end=========//


//======write vel_dist_bin======//
		if (flag_VD) {
		twrite_vels = clock();
		veldat = new char[ll+100];
		sprintf(veldat,"vel_dist_bin_%d.dat",N);
		if (fexists(veldat)) sprintf(veldat,"vel_dist_bin_%d_%02d_%02d_%04d_%02d%02d%02d.dat",N,day,month,year,hour,minute,second);
		FILE *fpdatvel = fopen(veldat,"w");

		fprintf(fpdatvel,"%E %E %d\n",vmin,vmax,iV);
		for(int iv1 = 0; iv1 < N; iv1++) {
			for (int iv2 = 0; iv2 < iV; iv2++) {
				fprintf(fpdatvel,"%d %d %d\n",vel[iv1*3*iV+iv2],vel[iv1*3*iV+iv2+iV],vel[iv1*3*iV+iv2+2*iV]);
			};
//			fprintf(fpdatvel,"\n");
		};
		fclose(fpdatvel);
		twrite_vel = clock() - twrite_vels;
		};
//=========end=========//	


//======write MSD======//	
		if (flag_MSD) {
		msddat = new char[ll+100];
		sprintf(msddat,"msd_bin_%d.dat",N);
		if (fexists(msddat)) sprintf(msddat,"msd_bin_%d_%02d_%02d_%04d_%02d%02d%02d.dat",N,day,month,year,hour,minute,second);
		FILE *fpdatmsd = fopen(msddat,"w");

		fprintf(fpdatmsd,"# Mean Square Displacement:\n");
		fprintf(fpdatmsd,"# (x-x0)^2 (y-y0)^2 (z-z0)^2 (x-x0)^2+(y-y0)^2+(z-z0)^2\n");
		for (int y1 = 0; y1 <N*(sizeY-1); y1++) {
			fprintf(fpdatmsd,"%.8f %.8f %.8f %.8f\n",vector_msd[0+y1*4],vector_msd[1+y1*4],vector_msd[2+y1*4],vector_msd[3+y1*4]);
		};
		fclose(fpdatmsd);
		};
//=========end MSD=========//


//======write VACF======//	
		if (flag_VACF) {
		vacfdat = new char[ll+100];
		sprintf(vacfdat,"vacf_bin_%d.dat",N);
		if (fexists(vacfdat)) sprintf(vacfdat,"vacf_bin_%d_%02d_%02d_%04d_%02d%02d%02d.dat",N,day,month,year,hour,minute,second);
		FILE *fpdatvacf = fopen(vacfdat,"w");

		fprintf(fpdatvacf,"vacf: (vx*vx0) (vy*vy0) (vz*vz0) (vx*vx0)+(vy*vy0)+(vz*vz0)\n");
		for (int y1 = 0; y1 <N*(sizeY-1); y1++) {
			fprintf(fpdatvacf,"%E %E %E %E\n",vector_vacf[0+y1*4],vector_vacf[1+y1*4],vector_vacf[2+y1*4],vector_vacf[3+y1*4]);
		};
		fclose(fpdatvacf);
		};
//=========end VACF=========//


//		msd = 0;
//		vacf = 0;

		trun1s = clock();

		fprintf(screen,"\n%s\n",lineD);
		fprintf(screen,"RUN 1");
		fprintf(screen,"\n%s\n",lineS);
		fprintf(screen,"Time steps:");
		run = 1;

		step = 0;
//  	double start_omp2 = omp_get_wtime();
//		start2 = clock();

		ftemp = fopen("temp.txt","w");
		fprintf(ftemp,"%d %d %d\n",start,stop,step_size);
		fprintf(ftemp,"%d %E %E\n",N,ZLO,ZHI);

		while(1) {
	
			fread(&ntimestep,sizeof(bigint),1,fp);
	
			if (feof(fp) || (ntimestep > stop) ) {
				fclose(fp);
				break;
			}
	
			if ((step%output_every)==0) fprintf(screen," " BIGINT_FORMAT, ntimestep);
			fflush(stdout);		
			step += step_size;
	
			trun1_reads = clock();
			fread(&natoms,sizeof(bigint),1,fp);		
			fread(&triclinic,sizeof(int),1,fp);		
			fread(&boundary[0][0],6*sizeof(int),1,fp);	
      fread(&xlo,sizeof(double),1,fp);
      fread(&xhi,sizeof(double),1,fp);
      fread(&ylo,sizeof(double),1,fp);
	    fread(&yhi,sizeof(double),1,fp);
      fread(&zlo,sizeof(double),1,fp);
      fread(&zhi,sizeof(double),1,fp);
      if (triclinic) {
				fread(&xy,sizeof(double),1,fp);
				fread(&xz,sizeof(double),1,fp);
				fread(&yz,sizeof(double),1,fp);
      }
      fread(&size_one,sizeof(int),1,fp);
      fread(&nchunk,sizeof(int),1,fp);
			trun1_read += clock() - trun1_reads;

			if (ntimestep >= start && ntimestep <= stop) { 
	
// prop = dof, ke, ke_orth, ke_paral, Sxx, Syy, Szz, Sxy, Sxz, Syz, Jcx, Jcy, Jcz, Jvx, Jvy, Jvz

				if (prop) delete [] prop;
				prop = new double[N*M];

				for(int i1 = 0; i1 < N*M; i1++ ) prop[i1] = 0.0;

				for (int i = 0; i < nchunk; i++) {
					trun1_reads = clock();
					fread(&n,sizeof(int),1,fp);				

					if (n > maxbuf) {
						if (buf) delete [] buf;
						buf = new double[n];
						maxbuf = n;
					}

					fread(buf,sizeof(double),n,fp);
					trun1_read += clock() - trun1_reads;
			
					trun1_assigns = clock();
  				assignparticle(indexrdf,bin,N,n/size_one,bulk_vel,size_one,variance,sizeY,coord,vel,0);
					trun1_assign += clock() - trun1_assigns;

//				rdf(Tsize,buf,gr,coord,Nra,size_one,indexrdf);
				}

//				#pragma omp parallel for	
				trun1_props = clock();
				for (int l = 0; l < N; l++) {
					temperature(l,temp,varT,varT_orth,varT_paral,sizeY);		
					pressure(l,press,bin_vol);
					heatflux(l,hflux,bin_vol);
					enthalpy(l,enthlp);
					energyflux(l,eflux,bin_vol,bulk_vel,hflux,press,Nstep);
				}
				trun1_prop += clock() - trun1_props;

//				fprintf(ftemp,"\n");

			} else {
				trun1_reads = clock();
				for (int i = 0; i < nchunk; i++) {
					fread(&n,sizeof(int),1,fp);

					// extend buffer to fit chunk size
					if (n > maxbuf) {
						if (buf) delete [] buf;
						buf = new double[n];
						maxbuf = n;
					}

					// read chunk and write as size_one values per line
					fread(buf,sizeof(double),n,fp);
				} // end of chunk
				trun1_read += clock() - trun1_reads;

			} // end of start
		} // end of while
		trun1 += clock() - trun1s;
		fprintf(screen,"\n%s\n",lineD);


	fclose(ftemp);

//=====WRITING====//
	twrite_bins = clock();
	int l = 1;
	while (N/(1.0+pow(10.0,(double)l-1.0)) >= 10.0) l++;

	filedat = new char[l+100];
	sprintf(filedat,"output_bin_%d.dat",N);
	if (fexists(filedat)) sprintf(filedat,"output_bin_%d_%02d_%02d_%04d_%02d%02d%02d.dat",N,day,month,year,hour,minute,second);
	FILE *fpdat = fopen(filedat,"w");

	fprintf(fpdat,"z[m] density[#/m3] Ux[m/s] Uy[m/s] Uz[m/s] T[K] T_orth[K] T_paral[K] P[N/m2] Pxx[N/m2] Pyy[N/m2] Pzz[N/m2] Pxy[N/m2] Pxz[N/m2] Pyz[N/m2] Jx[W/m2] Jy[W/m2] Jz[W/m2] Ex[W/m2] Ey[W/m2] Ez[W/m2] Pe[J] Enthalpy[J/kg] Occupation VarT varT_orth varT_paral bin_mean_z[m]\n");

	double NvarT;
	double NvarT_orth;
	double NvarT_paral;

	for (int iw = 0; iw < N; iw++) {
//		if (time_div[iw] == 0.0) time_div[iw] = 1.0;
		fprintf(fpdat,"%E ",bin_w[iw]);
		fprintf(fpdat,"%E ",den[iw]); //already divided by time_div
		fprintf(fpdat,"%E %E %E ",bulk_vel[iw*3],bulk_vel[iw*3+1],bulk_vel[iw*3+2]);
		fprintf(fpdat,"%E %E %E ",
						temp[iw][0]/Nstep,temp[iw][1]/Nstep,temp[iw][2]/Nstep);
		fprintf(fpdat,"%E %E %E %E %E %E %E ",
						(press[iw][0] + press[iw][1] + press[iw][2]) / (3 * Nstep),
						press[iw][0]/Nstep,press[iw][1]/Nstep,press[iw][2]/Nstep,
						press[iw][3]/Nstep,press[iw][4]/Nstep,press[iw][5]/Nstep);
		fprintf(fpdat,"%E %E %E ",
						hflux[iw][0]/Nstep,hflux[iw][1]/Nstep,hflux[iw][2]/Nstep);
		fprintf(fpdat,"%E %E %E ",
						eflux[iw][0]/Nstep,eflux[iw][1]/Nstep,eflux[iw][2]/Nstep);
		fprintf(fpdat,"%E ",Pe[iw]); //already divided by time_div
		fprintf(fpdat,"%E ",enthlp[iw]/Nstep);
		fprintf(fpdat,"%E ",time_div[iw]/Nstep); //fraction of total steps in which there is an atom in the bin


		NvarT = 0.0;
		NvarT_orth = 0.0;
		NvarT_paral = 0.0;

		for (int jw = 1; jw < sizeY; jw ++) {
			NvarT += (varT[sizeY*iw+jw] - (temp[iw][0]/Nstep))*(varT[sizeY*iw+jw] - (temp[iw][0]/Nstep));
		  NvarT_orth += (varT_orth[sizeY*iw+jw] - (temp[iw][1]/Nstep))*(varT_orth[sizeY*iw+jw] - (temp[iw][1]/Nstep));
			NvarT_paral += (varT_paral[sizeY*iw+jw] - (temp[iw][2]/Nstep))*(varT_paral[sizeY*iw+jw] - (temp[iw][2]/Nstep));
		};

		fprintf(fpdat,"%E %E %E ",NvarT/varT[sizeY*iw],NvarT_orth/varT_orth[sizeY*iw],NvarT_paral/varT_paral[sizeY*iw]);


		fprintf(fpdat,"%E ",bin_mean_z[iw]); //mean z-value per bin

//		fprintf(fpdat,"%E %E",Pe[iw],Ke[iw]); //already divided by time_div


		fprintf(fpdat,"\n");

	};

	fclose(fpdat);
	twrite_bin = clock() - twrite_bins;


	fprintf(screen,"\n");

	if (flag_VAR) {
		twrite_vars = clock();
		vardat = new char[l+100];
		sprintf(vardat,"var_bin_%d.dat",N);
		if (fexists(vardat)) sprintf(vardat,"var_bin_%d_%02d_%02d_%04d_%02d%02d%02d.dat",N,day,month,year,hour,minute,second);
		FILE *fpdatvar = fopen(vardat,"w");

		for(int idx = 0; idx < sizeY; idx++) {
			for (int iw = 0; iw < N; iw++) {
				fprintf(fpdatvar,"%E ",variance[idx+iw*sizeY]);
			};
			fprintf(fpdatvar,"\n");
		};

		fclose(fpdatvar);
		twrite_var = clock() - twrite_vars;
	};


//	totaltime = (clock()-starttime)/((double)CLOCKS_PER_SEC);
//	totaltime = omp_get_wtime()-start_omp;
	ttotal = (clock() - ttotals);
	double total_assign = (tbin_assign+trun1_assign);
	double total_prop = (tbin_prop+trun1_prop);
	double total_read = (tbin_read+trun1_read);
	double total_write = (twrite_gr+twrite_vel+twrite_var+twrite_bin);

	double total_time = ttotal/(double)CLOCKS_PER_SEC;
	int hours = floor(total_time/3600.0);
	int minutes = floor((total_time - (hours*3600))/60.0);
	double seconds = total_time - (hours*3600.0) - (minutes*60.0);

	if (post) {	
		fprintf(screen,"%s\n",lineD);
		fprintf(screen,"TASK TIMING BREAKDOWN:\n");
		fprintf(screen,"Section | total [hr:min:sec]\n");
		fprintf(screen,"%s\n",lineS);
		fprintf(screen,"INIT    | %02.0f:%02.0f:%02.5f\n",fh(tinit),fm(tinit),fs(tinit));
		fprintf(screen,"%s\n",lineS);
		fprintf(screen,"BIN     | %02.0f:%02.0f:%02.5f\n",fh(tbin),fm(tbin),fs(tbin));
		fprintf(screen,"-Read   | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(tbin_read),fm(tbin_read),fs(tbin_read),100*tbin_read/tbin);
		fprintf(screen,"-Assign | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(tbin_assign),fm(tbin_assign),fs(tbin_assign),100*tbin_assign/tbin);
		fprintf(screen,"-Occup. | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(tbin_occup),fm(tbin_occup),fs(tbin_occup),100*tbin_occup/tbin);
		fprintf(screen,"-Prop.  | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(tbin_prop),fm(tbin_prop),fs(tbin_prop),100*tbin_prop/tbin);
		fprintf(screen,"-Bin    | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(tbin_bin),fm(tbin_bin),fs(tbin_bin),100*tbin_bin/tbin);
		fprintf(screen,"%s\n",lineS);
		fprintf(screen,"RUN 1   | %02.0f:%02.0f:%02.5f\n",fh(trun1),fm(trun1),fs(trun1));
		fprintf(screen,"-Read   | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(trun1_read),fm(trun1_read),fs(trun1_read),100*trun1_read/trun1);
		fprintf(screen,"-Assign | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(trun1_assign),fm(trun1_assign),fs(trun1_assign),100*trun1_assign/trun1);
		fprintf(screen,"-Prop.  | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(trun1_prop),fm(trun1_prop),fs(trun1_prop),100*trun1_prop/trun1);
		fprintf(screen,"%s\n",lineS);
		fprintf(screen,"WRITE   | %02.0f:%02.0f:%02.5f\n",fh(total_write),fm(total_write),fs(total_write));
		fprintf(screen,"-g(r)   | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(twrite_gr),fm(twrite_gr),fs(twrite_gr),100*twrite_gr/total_write);
		fprintf(screen,"-Vel.   | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(twrite_vel),fm(twrite_vel),fs(twrite_vel),100*twrite_vel/total_write);
		fprintf(screen,"-Var.   | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(twrite_var),fm(twrite_var),fs(twrite_var),100*twrite_var/total_write);
		fprintf(screen,"-Bin    | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(twrite_bin),fm(twrite_bin),fs(twrite_bin),100*twrite_bin/total_write);
		fprintf(screen,"%s\n",lineS);
		fprintf(screen,"TOTAL   | %02.0f:%02.0f:%02.5f\n",fh(ttotal),fm(ttotal),fs(ttotal));
		fprintf(screen,"-Assign | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(total_assign),fm(total_assign),fs(total_assign),100*total_assign/ttotal);
		fprintf(screen,"--Bin   | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(tbin2),fm(tbin2),fs(tbin2),100*tbin2/total_assign);
		fprintf(screen,"-Prop.  | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(total_prop),fm(total_prop),fs(total_prop),100*total_prop/ttotal);
		fprintf(screen,"-Read   | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(total_read),fm(total_read),fs(total_read),100*total_read/ttotal);
		fprintf(screen,"-Write  | %02.0f:%02.0f:%02.5f (%05.2f%%)\n",fh(total_write),fm(total_write),fs(total_write),100*total_write/ttotal);
		fprintf(screen,"%s\n\n",lineD);
	}


	fprintf(screen,"Created output file(s):\n");
	if (flag_GR) fprintf(screen,"> %s\n",filedatgr);
	if (flag_VD) fprintf(screen,"> %s\n",veldat);
	if (flag_MSD) fprintf(screen,"> %s\n",msddat);
	if (flag_VACF) fprintf(screen,"> %s\n",vacfdat);
	if (flag_VAR) fprintf(screen,"> %s\n",vardat);
	if (flag_REGION) fprintf(screen,"> %s\n",fnameregion);
	fprintf(screen,"> %s\n\n",filedat);
	
	fprintf(screen,"Total time: %02d:%02d:%f (hr:min:sec)]\n",hours,minutes,seconds);
	fprintf(screen,"\n");


	if (variance) delete [] variance;
	if (vector_msd) delete [] vector_msd;
	if (vector_vacf) delete [] vector_vacf;
	if (xoriginal) delete [] xoriginal;
	if (voriginal) delete [] voriginal;

	if (varT) delete [] varT;
	if (varT_orth) delete [] varT_orth;
	if (varT_paral) delete [] varT_paral;

	if (varT1) delete [] varT1;
	if (varT_orth1) delete [] varT_orth1;
	if (varT_paral1) delete [] varT_paral1;

	if (prop) delete [] prop;
	if (prop_U) delete [] prop_U;

	if (bin) delete [] bin;  
	if (Nocpt) delete [] Nocpt; 
	if (Nsize) delete [] Nsize; 
	if (time_div) delete [] time_div; 
	if (den) delete [] den; 
	if (bulk_vel) delete [] bulk_vel; 
	if (zero_bulk) delete [] zero_bulk; 

	if (bin_w) delete [] bin_w; 
	if (bin_vol) delete [] bin_vol; 

	if (buf) delete [] buf;

	if (filedatgr) delete [] filedatgr;
	if (veldat)	delete [] veldat;
	if (msddat)	delete [] msddat;
	if (vacfdat)	delete [] vacfdat;
	if (filedat) delete [] filedat;
	if (vardat) delete [] vardat;
	if (fnameregion) delete [] fnameregion;

	if (DATA) delete [] DATA;

	return 0;
} // end of main
