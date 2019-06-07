//============================================================================
//TODO: enthalpy, surface tension, green-kubo relations,velocity distr
//============================================================================

//=============Chunk Domain=================
int chunkdomain(double *bin, int nz)
{
	double dz = 0.0;

	dz = (ZHI-ZLO) / (double)nz;

	for (int i = 0; i <= nz; i++) {
		bin[i] = ZLO + i * dz;
	};
//	for (int i = 0; i < nz; i++) {
//		bin_vol[i] = (xhi-xlo) * (yhi-ylo) * (bin[i+1] - bin[i]) * 1e-30; //m^3
//		printf("Dx:%f - Dy:%f - Dz:%f - Vol:%E\n",xhi-xlo,yhi-ylo,bin[i+1]-bin[i],bin_vol[i]);
//	};
}

//=============AssignParticle=================
void assignparticle(int *indexrdf, double *bin, int N, int nat,double *ibulk_vel,int size_one, double *ivar, int sizeY, double *coord, int *vel, int J)
{
	double m 			= 0.0;
	double vx 		= 0.0;
	double vy 		= 0.0;
	double vz 		= 0.0;
	double Keap 	= 0.0;
	double Kea 		= 0.0;
	double kea		= 0.0;
	double indexv	= 0.0;
	int jj 				= 0;
//	double iV = 100.0;
//	double vmin = -0.02;
//	double vmax = 0.02;
	double dv = (vmax-vmin)/iV;

//	FILE *fpdatvar = fopen("interface.dat","w");

//			0			1			2		3			4		5		6		7		8		9		10		11		12							17
// buf: mass 	mol		id	type 	x 	y 	z 	vx 	vy 	vz 	kea 	pea 	stressa[1] .... stressa[6] 

//	#pragma omp parallel for
	for (int j = 0; j < nat; j++) {

//		int ij=0;
		int i = 0;

		double z_tmp = buf[6+size_one*j];
		
		if (equalbin) {
			 double p = floor((z_tmp-ZLO)*N/(ZHI-ZLO));
				i = (int)sqrt(p*p);
//			 printf("i: %d\n",i);
		} else {
			for (int ii = 0; ii < N; ii++){
				if (z_tmp < bin[ii+1]) {
					i = ii;
					break;
				};
				i = ii;
			};
		};

		if (i>=N) i = N-1;


		if (z_tmp > 25 && z_tmp < 32) {
			fprintf(fregion,"%E %E %E %E\n",
				buf[2+size_one*j],buf[4+size_one*j],buf[5+size_one*j],buf[6+size_one*j]);	
		};




//		int i = 0;	
//		while (buf[6+size_one*j] >= bin[i] && buf[6+size_one*j] >= bin[i+1] && i < N) i++;
//		if (i == N) i = N - 1;
//		if (i == 0) {
//			if (buf[6+size_one*j] < bin[0]) i = 0;
//		};

		int iVX,iVY,iVZ;

		if (buf[7+size_one*j] <= vmin) {
			iVX = 0;
		} else if (buf[7+size_one*j] >= vmax) {
			iVX = iV-1;
		} else {
			iVX = (int)floor((buf[7+size_one*j]-vmin)*iV/(vmax-vmin));
		};

		if (buf[8+size_one*j] <= vmin) {
			iVY = 0;
		} else if (buf[8+size_one*j] >= vmax) {
			iVY = iV-1;
		} else {
			iVY = (int)floor((buf[8+size_one*j]-vmin)*iV/(vmax-vmin));
		};

		if (buf[9+size_one*j] <= vmin) {
			iVZ = 0;
		} else if (buf[9+size_one*j] >= vmax) {
			iVZ = iV-1;
		} else {
			iVZ = (int)floor((buf[9+size_one*j]-vmin)*iV/(vmax-vmin));
		};

//		int iVX = (int)sqrt(ivx*ivx);
//		int iVY = (int)sqrt(ivy*ivy);
//		int iVZ = (int)sqrt(ivz*ivz);

		vel[iVX+i*3*iV] += 1.0;
		vel[iVY+iV+i*3*iV] += 1.0;
		vel[iVZ+2*iV+i*3*iV] += 1.0;

		//index list for rdf
//		if (buf[4+size_one*j] >= coord[0] && buf[4+size_one*j] <= coord[0]+coord[3] &&
//			 	buf[5+size_one*j] >= coord[1] && buf[5+size_one*j] <= coord[1]+coord[3] &&
//				buf[6+size_one*j] >= coord[2] && buf[6+size_one*j] <= coord[2]+coord[3] && jj < 10) {
//				printf("\ntest:%d\n",jj);
//				indexrdf[jj+1] = j;
//				jj += 1;
//		};
//	
//		indexrdf[0] = jj;

	prop[i*M+29] += z_tmp * 1e-10;


	// Molar mass [kg]
	m  = buf[size_one*j] / (Na * 1000.0);
//	prop[i][29] = buf[size_one*j]/(1000.0 * Na);

	// Peculiar velocity [m/s]
	vx = ( buf[7+size_one*j] - ibulk_vel[i*3] * 1.0e-5 ) * 1e5;
	vy = ( buf[8+size_one*j] - ibulk_vel[i*3+1] * 1.0e-5 ) * 1e5;
	vz = ( buf[9+size_one*j] - ibulk_vel[i*3+2] * 1.0e-5 ) * 1e5;

	// Degrees of freedom
//	prop[i][0] += 3.0;
	prop[i*M] += 3.0;

	// Kinetic energy [J]
	prop[i*M+1] += m * ( vx * vx + vy * vy + vz * vz	);
	prop[i*M+2] += m * ( vx * vx + vy * vy	); 
	prop[i*M+3] += m * ( vz * vz	); 

//	double indexv = ivar[i*sizeY];
//	double kea = m * ( vx * vx + vy * vy + vz * vz	) * 1.0e10 / (Na * 1000.0);
//	printf("index: %d \n",(int)indexv);
//	ivar[i*sizeY+(int)indexv+1] = kea / (kb * 3.0);
//	ivar[i*sizeY] += 1.0;	

	// stress components [J] NOT divided by volume yet
	prop[i*M+4] += buf[12+size_one*j] * 101325.0 * 1e-30 - 
								m * vx * vx;
	prop[i*M+5] += buf[13+size_one*j] * 101325.0 * 1e-30 - 
								m * vy * vy;
	prop[i*M+6] += buf[14+size_one*j] * 101325.0 * 1e-30 - 
								m * vz * vz;
	prop[i*M+7] += buf[15+size_one*j] * 101325.0 * 1e-30 - 
								m * vx * vy;
	prop[i*M+8] += buf[16+size_one*j] * 101325.0 * 1e-30 - 
								m * vx * vz;
	prop[i*M+9] += buf[17+size_one*j] * 101325.0 * 1e-30 - 
								m * vy * vz;


	Keap =  0.5 * m * ( vx * vx + vy * vy + vz * vz ); // [J]

	// Convective part heat flux, NOT divided by volume yet [Wm]
	prop[i*M+10] += ( Keap + (buf[11+size_one*j]*4184/Na) ) * vx;	
	prop[i*M+11] += ( Keap + (buf[11+size_one*j]*4184/Na) ) * vy;	
	prop[i*M+12] += ( Keap + (buf[11+size_one*j]*4184/Na) ) * vz;	

	// Virial part heat flux, NOT divided by volume yet [Wm]
	prop[i*M+13] += ( buf[12+size_one*j] * vx + buf[15+size_one*j] * vy + 
									 buf[16+size_one*j] * vz ) * 101325.0 * 1.0e-30; 
	prop[i*M+14] += ( buf[15+size_one*j] * vx + buf[13+size_one*j] * vy + 
									 buf[17+size_one*j] * vz ) * 101325.0 * 1.0e-30; 
	prop[i*M+15] += ( buf[16+size_one*j] * vx + buf[17+size_one*j] * vy + 
									 buf[14+size_one*j] * vz ) * 101325.0 * 1.0e-30; 

	// Bulk velocities [m/s]
	prop[i*M+16] += buf[7+size_one*j] * 1.0e5;
	prop[i*M+17] += buf[8+size_one*j] * 1.0e5;
	prop[i*M+18] += buf[9+size_one*j] * 1.0e5;
		
	Kea =  0.5 * m * ( 	buf[7+size_one*j] * buf[7+size_one*j] + 
											buf[8+size_one*j] * buf[8+size_one*j] + 
											buf[9+size_one*j] * buf[9+size_one*j] ) * 
				1.0e10; // [J]

	// Convective part heat flux, NOT divided by volume yet [J]
	prop[i*M+19] += ( Kea + (buf[11+size_one*j]*4184/Na) );	

	//Potential energy [J]
	prop[i*M+20] += ( buf[11+size_one*j] ) * 4184.0 / Na;	

	// For enthalpy [m/s]
	prop[i*M+21] += vx * vx + vy * vy + vz * vz;

	// For enthalpy [J]
	prop[i*M+22] += ( buf[12+size_one*j] + buf[13+size_one*j] + buf[14+size_one*j]) * 101325.0 * 1e-30;	

	// Convective part heat flux, NOT divided by volume yet [Wm]
	prop[i*M+23] += ( Kea + (buf[11+size_one*j]*4184/Na) ) * buf[7+size_one*j] * 1.0e5;	
	prop[i*M+24] += ( Kea + (buf[11+size_one*j]*4184/Na) ) * buf[8+size_one*j] * 1.0e5;	
	prop[i*M+25] += ( Kea + (buf[11+size_one*j]*4184/Na) ) * buf[9+size_one*j] * 1.0e5;	

	// Virial part heat flux, NOT divided by volume yet [Wm]
	prop[i*M+26] += ( buf[12+size_one*j] * buf[7+size_one*j] + buf[15+size_one*j] * buf[8+size_one*j] + 
									 buf[16+size_one*j] * buf[9+size_one*j] ) * 1.0e5 * 101325.0 * 1.0e-30; 
	prop[i*M+27] += ( buf[15+size_one*j] * buf[7+size_one*j] + buf[13+size_one*j] * buf[8+size_one*j] + 
									 buf[17+size_one*j] * buf[9+size_one*j] ) * 1.0e5 * 101325.0 * 1.0e-30; 
	prop[i*M+28] += ( buf[16+size_one*j] * buf[7+size_one*j] + buf[17+size_one*j] * buf[8+size_one*j] + 
									 buf[14+size_one*j] * buf[9+size_one*j] ) * 1.0e5 * 101325.0 * 1.0e-30; 


	//Potential energy [J]
	prop[i*M+30] += ( buf[11+size_one*j] ) * 4184.0 / Na;	

	//Kinetic energy [J]
	prop[i*M+31] += ( buf[10+size_one*j] ) * 4184.0 / Na;	


/*	
	if (ntimestep==start && run == 1) {
		int JJ = J + j;
		xoriginal[0+JJ*4] = buf[2+size_one*j]; //id
		xoriginal[1+JJ*4] = buf[4+size_one*j]; //x
		xoriginal[2+JJ*4] = buf[5+size_one*j]; //y
		xoriginal[3+JJ*4] = buf[6+size_one*j]; //z

//		voriginal[0+JJ*4] = buf[2+size_one*j]; //id
		voriginal[1+JJ*4] = buf[7+size_one*j]; //vx
		voriginal[2+JJ*4] = buf[8+size_one*j]; //vy
		voriginal[3+JJ*4] = buf[9+size_one*j]; //vz
		JJ += 1;
	};

	int ka = 0;
	if (run == 1 && (msd | vacf)) {
		while((buf[2+size_one*j] != xoriginal[0+ka*4]) && (ka < natoms)) ka++; //find id in xoriginal
	};

	if (msd && run == 1) {
		double dx = buf[4+size_one*j] - xoriginal[1+ka*4];
		double dy = buf[5+size_one*j] - xoriginal[2+ka*4];
		double dz = buf[6+size_one*j] - xoriginal[3+ka*4];

//		double dxtmp = buf[4+size_one*j] - xoriginal[1+ka*4];
//		double dytmp = buf[5+size_one*j] - xoriginal[2+ka*4];
//		double dztmp = buf[6+size_one*j] - xoriginal[3+ka*4];

//		int MAXDX = (dxtmp > (xhi-xlo)*0.5) - (dxtmp < -(xhi-xlo)/2.0); //20*abs(buf[7+size_one*j]*step_size*dt));
//		int MAXDY = (dytmp > (yhi-ylo)*0.5) - (dytmp < -(yhi-ylo)/2.0); //20*abs(buf[8+size_one*j]*step_size*dt));
//		int MAXDZ = (dztmp > (zhi-zlo)*0.5) - (dztmp < -(zhi-zlo)/2.0); //20*abs(buf[9+size_one*j]*step_size*dt));

//		double dx = dxtmp - (xhi-xlo) * (MAXDX); //unwrap x-position
//		double dy = dytmp - (yhi-ylo) * (MAXDY); //unwrap y-position
//		double dz = dztmp - (zhi-zlo) * (MAXDZ); //unwrap z-position


//		if (MAXDY==1) {
//			printf("\nstep: %d - id: %f - k:%d - dytmp:%f - vel: %f - MAXDY: %d - dy:%f - yor:%f - y:%f\n",nstep*step_size,buf[2+size_one*j],ka,dytmp,buf[8+size_one*j]*step_size*dt,MAXDY,dy,xoriginal[2+ka*4],buf[5+size_one*j]);
//		};

		vector_msd[0+i*4+nstep*4*N] += dx*dx; 
		vector_msd[1+i*4+nstep*4*N] += dy*dy; 
		vector_msd[2+i*4+nstep*4*N] += dz*dz; 
		vector_msd[3+i*4+nstep*4*N] += dx*dx + dy*dy + dz*dz; 
	};

	//VACF
	if (vacf && run == 1) {
		double vxsq = buf[7+size_one*j] * voriginal[1+ka*4];
		double vysq = buf[8+size_one*j] * voriginal[2+ka*4];
		double vzsq = buf[9+size_one*j] * voriginal[3+ka*4];

		vector_vacf[0+i*4+nstep*4*N] += vxsq;
		vector_vacf[1+i*4+nstep*4*N] += vysq;
		vector_vacf[2+i*4+nstep*4*N] += vzsq;
		vector_vacf[3+i*4+nstep*4*N] += vxsq + vysq + vzsq;
	};

*/

	} // end of fore loop over natoms
//	fclose(fpdatvar); 
} // end of function assignparticle()


//===============MSD=======================
void MSD(int l, int N) 
{
		double nmsd;
		if (prop[l*M] < 3.0) {
			nmsd = 1.0;
		} else {
			nmsd = prop[l*M] / 3.0;
		}
		vector_msd[0+l*4+nstep*4*N] /= nmsd;
		vector_msd[1+l*4+nstep*4*N] /= nmsd;
		vector_msd[2+l*4+nstep*4*N] /= nmsd;
		vector_msd[3+l*4+nstep*4*N] /= nmsd;
}

//===============VACF=======================
void VACF(int l, int N) 
{
		double nvacf;
		if (prop[l*M] < 3.0) {
			nvacf = 1.0;
		} else {
			nvacf = prop[l*M] / 3.0;
		}
		vector_vacf[0+l*4+nstep*4*N] /= nvacf;
		vector_vacf[1+l*4+nstep*4*N] /= nvacf;
		vector_vacf[2+l*4+nstep*4*N] /= nvacf;
		vector_vacf[3+l*4+nstep*4*N] /= nvacf;

}

//===============timediv=======================
void timediv(int N, double *itime_div)
{
	if (prop[N*M] > 1.0) itime_div[N] += 1.0;
}

//===============Temperature[K]====================
void temperature(int N, double itemp[][3], double *ivarT, double *ivarT_orth, double *ivarT_paral, int sizeY)
{
	double div,div_orth,div_paral;
	if (prop[N*M] < 3.0) {
		div = 1.0; // dividing by zero
		div_orth = 1.0;
		div_paral = 1.0;	
	} else {
		div = prop[N*M];
		div_orth = 2.0 * prop[N*M] / 3.0;
		div_paral = prop[N*M] / 3.0;
	};
	itemp[N][0] += prop[N*M+1] / ( kb * div );
 	itemp[N][1] += prop[N*M+2] / ( kb * div_orth ); 
	itemp[N][2] += prop[N*M+3] / ( kb * div_paral );	


//	fprintf(ftemp,"%E %E %E %E\n",bin_w[N],prop[N*M+1] / ( kb * div ),prop[N*M+2] / ( kb * div_orth ), prop[N*M+3] / ( kb * div_paral )); 
	fprintf(ftemp,"%E %E %E\n",prop[N*M+1] / ( kb * div ),prop[N*M+2] / ( kb * div_orth ), prop[N*M+3] / ( kb * div_paral )); 



	double indexv = ivarT[N*sizeY];
	ivarT[N*sizeY+(int)indexv+1] = prop[N*M+1] / ( kb * div );
	ivarT[N*sizeY] += 1.0;	

	double indexv1 = ivarT_orth[N*sizeY];
	ivarT_orth[N*sizeY+(int)indexv1+1] = prop[N*M+2] / ( kb * div_orth );
	ivarT_orth[N*sizeY] += 1.0;	

	double indexv2 = ivarT_paral[N*sizeY];
	ivarT_paral[N*sizeY+(int)indexv2+1] = prop[N*M+3] / ( kb * div_paral );
	ivarT_paral[N*sizeY] += 1.0;	


}

//================Number density[#/m3]================
void density(int N, double *den, double *ivar, int sizeY)
{
//	den[N] += prop[N*M] / (3.0 * ibin_vol[N]);
	den[N] += prop[N*M] / (3.0);


//	double indexv = ivar[N*sizeY];
//	ivar[N*sizeY+(int)indexv+1] = iprop[N][0] / ( 3.0 * ibin_vol[N] );
//	ivar[N*sizeY] += 1.0;	

}

//================Potential energy[J]================
void pe_energy(int N, double *Pe)
{
	Pe[N] += prop[N*M+20];	
}

//===============Pressure[N/m2]=======================
void pressure(int N, double ipress[][6], double *ibin_vol)
{
	for (int j = 0; j < 6; j++) {
		// pxx pyy pzz pxy pxz pyz
		ipress[N][j] += -1 * prop[N*M+j+4] / ibin_vol[N];
	}
}

//===============Heat flux[W/m2]=========================
void heatflux(int N, double ihflux[][3], double *ibin_vol)
{
	ihflux[N][0] += ( prop[N*M+10] - prop[N*M+13] ) / ibin_vol[N]; 
	ihflux[N][1] += ( prop[N*M+11] - prop[N*M+14] ) / ibin_vol[N];
	ihflux[N][2] += ( prop[N*M+12] - prop[N*M+15] ) / ibin_vol[N];
}

//===============Energy flux[W/m2]=========================
void energyflux1(int N, double ieflux[][3], double *ibin_vol, double *bulk_vel, double ihflux[][3], double ipress[][6], double itime_div) 
{
	ieflux[N][0] += (	prop[N*M+19] / ibin_vol[N]) * bulk_vel[N*3] + 
										ihflux[N][0] + 
										( ipress[N][0] ) * bulk_vel[N*3] +
										( ipress[N][3] ) * bulk_vel[N*3+1] + 
										( ipress[N][4] ) * bulk_vel[N*3+2];
	ieflux[N][1] += (	prop[N*M+19] / ibin_vol[N]) * bulk_vel[N*3+1] + 
										ihflux[N][1] + 
										( ipress[N][3] ) * bulk_vel[N*3] +
										( ipress[N][1] ) * bulk_vel[N*3+1] + 
										( ipress[N][5] ) * bulk_vel[N*3+2];
	ieflux[N][2] += (	prop[N*M+19] / ibin_vol[N]) * bulk_vel[N*3+2] + 
										ihflux[N][2] + 
										( ipress[N][4] ) * bulk_vel[N*3] +
										( ipress[N][5] ) * bulk_vel[N*3+1] + 
										( ipress[N][2] ) * bulk_vel[N*3+2]; 
}

//===============Energy flux[W/m2]=========================
void energyflux(int N, double ieflux[][3], double *ibin_vol) 
{
	ieflux[N][0] += ( prop[N*M+23] - prop[N*M+26] ) / ibin_vol[N]; 
	ieflux[N][1] += ( prop[N*M+24] - prop[N*M+27] ) / ibin_vol[N];
	ieflux[N][2] += ( prop[N*M+25] - prop[N*M+28] ) / ibin_vol[N];
}

//==============Enthalpy[J/kg]========================
void enthalpy(int N, double *ienthlp, double *ibin_vol)
{
	double m = prop[N*M+29];
	double Nj;
	if (prop[N*M] < 3.0) {Nj = 1.0;}
	else {Nj = prop[N*M];}

	ienthlp[N] = (3.0/Nj) * ( (5.0/6.0) * (prop[N*M+21]) + 
																			(1/m) * ( prop[N*M+20] + (1.0/3.0) * prop[N*M+22] ) );
}


/*
//==============Surface tension========================
void surface_tension()
{

}
*/

/*
//==============Green-kubo relations========================
void green_kubo()
{
thermal cond. -> lambda [W/mK]
mass transport -> Kc [m/s]
heat trnas. coeff. -> h [W/m2K]
Shear viscosity
}
*/

/*
//==============Velocity distribution========================
void velocity_dist()
{

}
*/

//==============Bulk velocity[m/s]========================
void mean_z(int N, double *bin_mean_z)
{
	double Vdiv;
	if (prop[N*M] < 3.0) {
		Vdiv = 1.0;
	} else {
		Vdiv = prop[N*M] / 3.0;
	};
	bin_mean_z[N] += prop[N*M+29] / Vdiv;
}

//==============Bulk velocity[m/s]========================
void bulk_velocity(int N, double *bulk_vel)
{
	double Vdiv;
	if (prop[N*M] < 3.0) {
		Vdiv = 1.0;
	} else {
		Vdiv = prop[N*M] / 3.0;
	};
	bulk_vel[N*3] += prop[N*M+16] / Vdiv;
	bulk_vel[N*3+1] += prop[N*M+17] / Vdiv;
	bulk_vel[N*3+2] += prop[N*M+18] / Vdiv;
}

/*
//==============Bulk velocity========================
void bulk_velocity_test(int N, double iprop[][22], double ibulk_vel[][3], double *ie_bv, double *itime_div)
{
	double Vdiv, vnew, vold, idiv, idiv_new, vold_div;
	double iold_bulk_vel[N][3];
	double old_time_div[N];
	double ierror = 0.0;
	double error = 0.0;
	
	for (int i = 0; i < N; i++) {
		if (iprop[i][0] < 3.0) {
			Vdiv = 1.0;
		} else {
			Vdiv = iprop[i][0] / 3.0;
		};
		if (itime_div[i] == 0.0) {
			idiv = 1.0;
		} else {
			idiv = itime_div[i];
		};
		iold_bulk_vel[i][0] = ibulk_vel[i][0] / idiv;
 		iold_bulk_vel[i][1] = ibulk_vel[i][1] / idiv;
		iold_bulk_vel[i][2] = ibulk_vel[i][2] / idiv;

//		old_time_div[i] = itime_div[i];

		ibulk_vel[i][0] += iprop[i][16] / Vdiv;
		ibulk_vel[i][1] += iprop[i][17] / Vdiv;
		ibulk_vel[i][2] += iprop[i][18] / Vdiv;
	};

	timediv(N,iprop,itime_div);

	for (int i = 0; i < N; i++) {
		if (itime_div[i] == 0.0) {
			idiv_new = 1.0;
		} else {
			idiv_new = itime_div[i];
		}; 
		vnew = sqrt(ibulk_vel[i][0] * ibulk_vel[i][0] + 
								ibulk_vel[i][1] * ibulk_vel[i][1] + 
								ibulk_vel[i][2] * ibulk_vel[i][2]) / idiv_new;
		vold = sqrt(iold_bulk_vel[i][0] * iold_bulk_vel[i][0] + 
								iold_bulk_vel[i][1] * iold_bulk_vel[i][1] +
								iold_bulk_vel[i][2] * iold_bulk_vel[i][2]);

		if (vold == 0.0) {
			vold_div = 1.0;
		} else {
			vold_div = vold;
		};
	
		error =  sqrt(  ((vnew - vold) / vold_div) * ((vnew - vold) / vold_div) );
	
		if (error > ierror) {
			ierror = error;
		};
	};

	*ie_bv = ierror;
//	printf("e_bv: %E\n\n",*ie_bv);
}
*/

//==============Radial Distribution Function==============
void rdf(int n, double *buf, double *gr, double *coord, double *Nra, int size_one, int *index)
{
	double x,y,z,dxyz,Ra;
	double Nrho = 0.0;
	int i = 0;
//	int index[10];
//	index[0]=-1;
//	int jj = 0;
//			0			1			2		3			4		5		6		7		8		9		10		11		12							17
// buf: mass 	mol		id	type 	x 	y 	z 	vx 	vy 	vz 	kea 	pea 	stressa[1] .... stressa[6] 

//	printf("\n\nRadial distribution function");
//	printf("\n===================================================");			

/*
	for (int j = 0; j < n; j++) {
//		if (jj > 9) break;
		if (buf[4+size_one*j] >= coord[0] && buf[4+size_one*j] <= coord[0]+coord[3] &&
			 	buf[5+size_one*j] >= coord[1] && buf[5+size_one*j] <= coord[1]+coord[3] &&
				buf[6+size_one*j] >= coord[2] && buf[6+size_one*j] <= coord[2]+coord[3] && jj < 10) {
//				printf("\natom:%d - x:%f - y:%f - z:%f",j,buf[4+size_one*j],buf[5+size_one*j],buf[6+size_one*j] );
				index[jj] = j;
				jj += 1;
		};
	}
*/

	if (index[0]==-1) {
		printf("\nNo atom selected with dxyz:%f",coord[3]);	
		return;
	}
	
	srand(time(NULL));
	int id = index[rand() % index[0]];
//	printf("\n---------------------------------------------------");		
//		printf("Radial distribution function: %d\n Stop: %d\n Nbins: %d\n",coord[0],coord[1],coord[2],coord[3],coord[4],coord[5],coord[6]);
//	printf("\nRandomly selected atom number:%d",id);
//	printf("\n===================================================\n");			

	for (int j = 0; j < n; j++) {
		if ( sqrt((buf[4+size_one*j]-buf[4+size_one*id]) * (buf[4+size_one*j]-buf[4+size_one*id]) + 
							(buf[5+size_one*j]-buf[5+size_one*id]) * (buf[5+size_one*j]-buf[5+size_one*id]) + 
							(buf[6+size_one*j]-buf[6+size_one*id]) * (buf[6+size_one*j]-buf[6+size_one*id])) <= coord[5] 
						&& (buf[6+size_one*j]-buf[6+size_one*id]) >= 0.0 ) {
			Nrho += 1.0;
		};

		i = 0;
		for (int p = 0; p < coord[6]; p++ ) {
			double r = (p+1)*coord[4]/coord[6];
			if ( sqrt(	(buf[4+size_one*j]-buf[4+size_one*id]) * 
									(buf[4+size_one*j]-buf[4+size_one*id]) + 
									(buf[5+size_one*j]-buf[5+size_one*id]) * 
									(buf[5+size_one*j]-buf[5+size_one*id]) + 
									(buf[6+size_one*j]-buf[6+size_one*id]) * 
									(buf[6+size_one*j]-buf[6+size_one*id])) <= r 
								&& (buf[6+size_one*j]-buf[6+size_one*id]) >= 0.0 ) {
				if (j != id) {Nra[i] += 1.0;};
				break;
			};
			i += 1;
		};

	}
	
	double rho = Nrho / ( (4.0/6.0) * pi * coord[5] * coord[5] * coord[5]); //volume of half of sphere
//	printf("rho:%f\n",rho);

	for (int i = 0; i < coord[6]; i++) {
		gr[i] = Nra[i] / (0.5 * 4 * pi * (coord[4]/coord[6]) * (coord[4]/coord[6])*(i+1) * (coord[4]/coord[6])*(i+1) * rho );
}
	 
	
}


