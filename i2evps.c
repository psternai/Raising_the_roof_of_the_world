#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mkl.h>

/* --------------------- INCLUDE PARTS --------------------- */
#include"headevps.c"
#include"loadevps.c"
#include"moveevps.c"
#include"markevps.c"
#include"gausevps.c"
#include"heatevps.c"
/* --------------------- INCLUDE PARTS --------------------- */

/* Solve differencial equations by MARKER-IN-CELL + FINITE-DIFFERENCES method */
int main()
{
/* Counters */
int n0,n1,f0,n3,fln3;
long int pos0cur0,m1,m2,m3;
FILE * tfl;
/**/
/**/
/**/
/* Load configuration from mode.t3c */
fln3=loadconf()+1;
/**/
/**/
/**/
/* Load data from input file */
loader();
/**/
/**/
/**/
/* Output File Cycle */
for (f0=fln3;f0<fl0num;f0++)
{
/* Reload Cur Output File Name, Type */
for (n1=0;n1<50;n1++) fl1out[n1]=fl0out[f0][n1];
fl1otp=fl0otp[f0];
/**/
/* Reload cyc0max_maxxystep_maxtkstep_maxtmstep_maxdsstep */
cyc0max=fl0cyc[f0];
maxxystep=fl0stp[f0][0];
maxtkstep=fl0stp[f0][1];
maxtmstep=fl0stp[f0][2];
xelvismin=fl0stp[f0][3];
/**/
/* General Cycle */
for (n0=0;n0<cyc0max;n0++)
	{
	if (printmod) printf("\n! FILE %s  KRUG ! %d\n",fl1out,n0+1);
	/**/
	/* Regrid due to the displacement */
/*
for (n3=1;n3<xnumx-1;n3++)
	{
	if(n3<21 || n3>170)
		{
		gx[n3]=gx[n3-1]+(0.40+timesum*1.5)/20.0;
		}
	else
		{
		gx[n3]=gx[n3-1]+(3.00-timesum*3.0)/150.0;
		}
	}
*/
	/**/
	/* Set initial time step */
	timestep=maxtmstep;
	timestepe=xelvismin;
	if (printmod) printf("\n !!! MAX VALID TIME STEP FOR CYCLE %e YEARS !!! \n",timestep/3.15576e+7);
	/**/
	/**/
	/**/
	/* vX,vY recalc after Stokes+Contin equation */
	if(movemod)
		{
		if (printmod) printf("\n EPS, SIG, P, VX, VY CALCULATION...\n");
		viterate(n0);
		if (printmod) printf("EPS, SIG, P, VX, VY  OK!\n");
		}
	/**/
	/**/
	/**/
	/* Time step for markers definition */
	if(markmod && !movemod)
		{
		maxvelstep();
		}
	/**/
	/**/
	/**/
	/* Tk recalc after Heat transport equation */
	if(timedir>0 && tempmod && timestep)
		{
		if (printmod) printf("\n TEMPERATURE CALCULATION...\n");
		titerate(n0);
		if (printmod) printf("TEMPERATURE OK!\n");
		}
	/**/
	/**/
	/**/
	/* Move marker */
	if(markmod && timestep)
		{
		/* Time step for markers definition */
		if(timedir<0 && timestep>0) timestep=-timestep;
		if (printmod) printf("\n CURRENT VALID TIME STEP %e YEARS IN CYCLE\n",timestep/3.15576e+7);
		if (printmod) printf("MARKERS MOVE...");
		movemark();
		if (printmod) printf("MARKERS OK!\n");
		}
	/**/
	/**/
	/**/
	/* Reset Vx, Vy, P values */
	if(ratemod)
		{
		if (printmod) printf("VX, VY, P RESET ...");
		/* Vx, Vy Reset Cycle */
		for (m1=0;m1<nodenum;m1++)
			{
			vx[m1]=vy[m1]=0;
			}
		/**/
		/* Pressure in cells Reset Cycle */
		for (m1=0;m1<xnumx1;m1++)
		for (m2=0;m2<ynumy1;m2++)
			{
			/* Pos in pr[] */
			m3=m1*ynumy1+m2;
			/* Recalc P */
			pr[m3]=pinit+((double)(m2)+0.5)*ystpy*GYKOEF*pkf[0];
			}
		if (printmod) printf("VX, VY, P OK!");
		}
	/**/
	/**/
	/**/
	/* ro[],nu[] Recalc */
	if(gridmod)
	if(!stoksmod || !movemod)
		{
		if (printmod) printf("\n RO, NU, CP etc  RECALC AFTER NEW MARKERS POSITIONS...");
		ronurecalc();
		if (printmod) printf("RO, NU, CP etc OK!\n");
		}
	/**/
	/**/
	/**/
	/* Increse Timesum */
	timesum+=timestep;
	/* 1year=3.15576*10^7sek */
	if (printmod) printf("\n %e YEARS IN CYCLE     %e YEARS FROM START\n\n",timestep/3.15576e+7,timesum/3.15576e+7);
	/**/
	/**/
	/**/
	/* Print Results */
/*
	if(n0<cyc0max-1) saver(f0+1,n0);
*/
	/* Print Results */
	}
/**/
/**/
/**/
/* Print Results */
saver(f0+1,n0-1);


/*** ARMEL'S UPDATE ***/

/*Print topo only */
	/* Load and open new TOPO file */	
	
	for (n1=0;n1<50;n1++) topo1out[n1]=topo0out[f0][n1];	
	fl = fopen(topo1out,"wt");

 	if (printmod) printf("\n NEW TOPO FILE %s ...\n",topo1out);

	/*Write TOPO information in this new file*/

	for (m1=0; m1<xnumx; m1++) 
	fprintf (fl, "%e %e %e \n", timesum/3.15576e+7, gx[m1], ep[m1]);
        if (printmod) printf("\n NEW TOPO FILE %s OK \n",topo1out);

	fclose(fl);


/* Print Density, Velocity, EPS, SIGMA, water content wa1 */
	/* Load and open new 2D data file */

	for (n1=0;n1<50;n1++) data1out[n1]=data0out[f0][n1];
        fl = fopen(data1out,"wt");

        if (printmod) printf("\n NEW 2D DATA FILE %s ...\n",data1out);
	
	/* Write in BINARY format */
		/* Define counters and buffers */
		char szlong,szfloat,szdouble;
		float ival0;
		double ival1;
			
		char nn1;
		long int mm1;

		szlong = sizeof(m1);
		szfloat = sizeof(ival0);
 		szdouble = sizeof(ival1);
			
		/* Write informations */ 
		fwrite(&xnumx,szlong,1,fl);
		fwrite(&ynumy,szlong,1,fl);
			
		ival1 = timesum/3.15576e+7; fwrite(&ival1,szdouble,1,fl);
		
		for (m1=0;m1<xnumx;m1++)
			{	
			ival0=(float)(gx[m1]); fwrite(&ival0,szfloat,1,fl);
			}
		for (m2=0;m2<ynumy;m2++)
			{
			ival0=(float)(gy[m2]); fwrite(&ival0,szfloat,1,fl);
			}			
			
		for (m1=0; m1<nodenum; m1++)
			{	
			ival0=(float)(nu[m1]); fwrite(&ival0,szfloat,1,fl);
			ival0=(float)(ro[m1]); fwrite(&ival0,szfloat,1,fl);
			ival0=(float)(vx[m1]); fwrite(&ival0,szfloat,1,fl);
			ival0=(float)(vy[m1]); fwrite(&ival0,szfloat,1,fl);
			ival0=(float)(exx[m1]); fwrite(&ival0,szfloat,1,fl);
			ival0=(float)(exy[m1]); fwrite(&ival0,szfloat,1,fl);			
			ival0=(float)(sxx[m1]); fwrite(&ival0,szfloat,1,fl);
			ival0=(float)(sxy[m1]); fwrite(&ival0,szfloat,1,fl);
			}
			
		if (printmod) printf("\n VELOCITY, EPS, ETC WRITTEN IN %s \n",data1out);
	
	fclose(fl);

/*** END OF ARMEL'S UPDATE ***/


/* Print Results */
/**/
/* End Output file Names Cycle */
/**/
// time tracking
    tfl = fopen("times.txt","a+");
    fprintf(tfl,"%8.4f \n", timesum/3.15576e+7);
    fclose(tfl);
// end of time tracking
}
/* End Program */
return 0;
}
/* Solve differencial equations by MARKER-IN-CELL + FINITE-DIFFERENCES method */
