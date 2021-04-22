#include "udf.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
/*User Defined Function for CFB Simulation using Fluent®. */

float Us[1001],e[201],Hd[1001][201];

DEFINE_EXECUTE_ON_LOADING(Hd_init, libname)
/*Load the data matrix of drag correction factor. */
{
#if !RP_HOST
	const int N=1001;
	const int M=201;
	int i,j;

	FILE *fp;

	if ((fp=fopen("Hd.dat","r"))==NULL)
	{
		Message("cannot open file!/n");
	}

	for(i=0;i<N;i++)
		for(j=0;j<M;j++)
		{
			fscanf(fp,"%f%f%f",&Us[i],&e[j],&Hd[i][j]); 
		}

	Message("Hd load successfully!\n");
	Message("%f,%f,%f\n",Us[0],e[0],Hd[0][0]);
	Message("%f,%f,%f\n",Us[1000],e[200],Hd[1000][200]);

	fclose(fp);
	#endif		
}

DEFINE_EXCHANGE_PROPERTY(qc_emms_drag, cell, mix_thread, s_col, f_col)
/*Calculate drag correction factor using linear interpolation. */
{
	#if !RP_HOST
	Thread *thread_g, *thread_s;

	real vel_g, vel_s, vel_slip, void_g;

	real logq, Uslip;
	real de;

	real Hd0,Hd1,Hd2,Hd3,Hd4,Hd5,Hd6;
	int N,M;
	real n,m;


	thread_g = THREAD_SUB_THREAD(mix_thread, s_col);
	thread_s = THREAD_SUB_THREAD(mix_thread, f_col);

	/*Velocities in CFB axial direction. May need revised to C_U, C_V according to CFB geometry. */
	vel_g = C_W(cell, thread_g);
    vel_s = C_W(cell, thread_s);
    vel_slip = vel_g - vel_s;

	void_g = C_VOF(cell, thread_g);

	
	Uslip=vel_slip*void_g;
	
	de=e[1]-e[0];

	logq=log(1.01);
	

	if ((void_g<e[0])||(void_g>e[200])||(Uslip<Us[0])||(Uslip>Us[1000]))
		Hd0=1;
	else
	{
		n=log(Uslip/Us[0])/logq;
		m=(void_g-e[0])/de;

		N=floor(n);
		M=floor(m);

		if ((N==1000)||(M==200))
			Hd0=1;
		else
		{
			Hd1=Hd[N][M];Hd2=Hd[N][M+1];Hd3=Hd[N+1][M];Hd4=Hd[N+1][M+1];

			Hd5=(Hd1*(e[M+1]-void_g)+Hd2*(void_g-e[M]))/de;
			Hd6=(Hd3*(e[M+1]-void_g)+Hd4*(void_g-e[M]))/de;

			Hd0=(Hd5*(Us[N+1]-Uslip)+Hd6*(Uslip-Us[N]))/(Us[N+1]-Us[N]);
		}
	}
	

	if (Hd0>1||Hd0<=0) 
	{
		Hd0=1;
	}

	C_UDMI(cell, mix_thread, 0) = Hd0;
	/*Storage Hd in Fluent® data. A User Defined Memory needs set in Fluent®. */


	return Hd0;

	#endif
}
