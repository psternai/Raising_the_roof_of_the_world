#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mkl.h>

/* --------------------- INCLUDE PARTS --------------------- */
#include"headjevps.c"
#include"loadjevps.c"
#include"moveevps.c"
#include"markevps.c"
#include"gausevps.c"
#include"heatevps.c"
/* --------------------- INCLUDE PARTS --------------------- */

/* Solve differencial equations by MARKER-IN-CELL + FINITE-DIFFERENCES method */
int main()
{
/* Counters */
int n0,n1,f0;
long int pos0cur0,m1,m2,m3;
double dx1,dx2;
/**/
/**/
/**/
/* Load configuration from mode.t3c */
loadconf();
/**/
/**/
/* Topography data print */
if(fl0otp[0]>=13)
	{
	/* Reload Cur Output File Name, Type */
	for (n1=0;n1<50;n1++) fl1out[n1]=fl0out[0][n1];
	fl1otp=fl0otp[0];
	if (printmod && fl0otp[0]==13) printf("Print topography data to %s... \n",fl1out);
	if (printmod && fl0otp[0]==14) printf("Print overriding plate length to %s... \n",fl1out);
	/* Open file */
	fl1 = fopen(fl1out,"wt");
	}
/**/
/**/
/**/
/* Output File Cycle */
for (f0=0;f0<fl0num;f0++)
{
/* Reload Cur Input File Name, Type */
for (n1=0;n1<50;n1++) fl1in[n1]=fl0in[f0][n1];
fl1itp=fl0itp[f0];
/**/
/* Load data from input file */
loader();
/**/
/* Normal postprocessing */
if(fl0otp[0]<13)
	{
	/* Reload Cur Output File Name, Type */
	for (n1=0;n1<50;n1++) fl1out[n1]=fl0out[f0][n1];
	fl1otp=fl0otp[f0];
	/**/
	/* Save data to output file */
	saver(f0,0);
	}
else
	{
	/* Topography postprocessing */
	if(fl0otp[0]==13)
		{
		if (f0==0) fprintf(fl1,"%ld \n",xnumx);
		for (m1=0;m1<xnumx;m1++) fprintf(fl1,"%e %e %e \n",timesum/3.15576e+7,gx[m1],ep[m1]);
		}
	/* Plate length postprocessing */
	if(fl0otp[0]==14)
		{
		/* Define nearest crust  marker to 2000 km */
		if (f0==0)
			{
			m2=-1;
			dx2=1e+30;
			for (m1=0;m1<marknum;m1++) 
			if((markt[m1]==5 && marky[m1]<1.2e+4 && marky[m1]>1.1e+4) && (m2<0 || ABSV(markx[m1]-2e+6)<dx2)) 
				{
				dx2=ABSV(markx[m1]-2e+6);
				m2=m1;
				}
			}
		/* Define leftmost crust  marker */
		dx1=-1e+30;
		for (m1=0;m1<marknum;m1++)
		if(markt[m1]==5 && (markx[m2]-markx[m1])>dx1) 
			{
			dx1=(markx[m2]-markx[m1]);
			m3=m1;
			}
		fprintf(fl1,"%e %e   %ld %d %e %e   %ld %d %e %e\n",timesum/3.15576e+7,dx1,m3,markt[m3],markx[m3],marky[m3],m2,markt[m2],markx[m2],marky[m2]);
		}
	}
}
/* Topography data print */
if(fl0otp[0]==13)
	{
	fclose(fl1);
	}
/* End Program */
return 0;
}
/* Solve differencial equations by MARKER-IN-CELL + FINITE-DIFFERENCES method */
