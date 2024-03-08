/* Load information from configuration file mode.t3c ============== */
int loadconf()
{
/* Counter */
int n1,n2,fln3=0;
/**/
/**/
/**/
/* Open File file.t3c */
fl = fopen("file.t3c","rt");
ffscanf(); fln3=atoi(sa)-1;
fclose(fl);
/**/
/**/
/**/
/* Open File mode.t3c */
fl = fopen("mode.t3c","rt");
/**/
/* Data File name */
ffscanf();
for (n1=0;n1<50;n1++) fl1in[n1]=sa[n1];
ffscanf(); if(sa[0] == 'b') fl1itp=1;
/**/
/* Load first Results File names */
ffscanf();
fl0num=0;
while(sa[0]!='~')
	{
	/* Check file Counter */
	if(fl0num>=MAXFLN) {printf("Space out in fl0out[]"); exit(0);}
	/**/
	/* Save results file name */
	for (n1=0;n1<50;n1++) fl0out[fl0num][n1]=sa[n1];
	/**/
	/* Load TYPE_cyc0max_maxxystep_maxtkstep_maxtmstep_maxdsstep */
	ffscanf(); if(sa[0] == 'b') fl0otp[fl0num]=1;
	ffscanf();fl0cyc[fl0num]=atoi(sa);
	ffscanf();fl0stp[fl0num][0]=atof(sa);
	ffscanf();fl0stp[fl0num][1]=atof(sa);
	ffscanf();fl0stp[fl0num][2]=atof(sa)*3.15576e+7;
	ffscanf();fl0stp[fl0num][3]=atof(sa)*3.15576e+7;


	/*** ARMEL'S UPDATE ***/
	
	/* New topography output file topo0out */
	/* Load new outpout file for topography only */
	ffscanf();
	for (n1=0;n1<50;n1++)  topo0out[fl0num][n1]=sa[n1]; 
	
	/* New 2D data file 2D0out */
	ffscanf();
	for (n1=0;n1<50;n1++) data0out[fl0num][n1]=sa[n1];
	
	/*** END OF ARMEL'S UPDATE ***/


	/**/
	/* Incr File Counters */
	fl0num++;
	/**/
	/* Load Next Results File names */
	ffscanf();
	}
/**/
/* Data File name change after number */
if(fln3>=0 && fln3<fl0num)
	{
	for (n1=0;n1<50;n1++) fl1in[n1]=fl0out[fln3][n1];
	fl1itp=fl0otp[fln3];
	}
else
	{
	fln3=-1;
	}
/**/
/* Service */
ffscanf();printmod=atoi(sa);
ffscanf();timedir=atof(sa);
ffscanf();movemod=atoi(sa);
ffscanf();tempmod=atoi(sa);
ffscanf();markmod=atoi(sa);
ffscanf();ratemod=atoi(sa);
ffscanf();gridmod=atoi(sa);
ffscanf();intermod=atoi(sa);
ffscanf();intermod1=atoi(sa);
ffscanf();outgrid=atoi(sa);
ffscanf();densimod=atoi(sa);
ffscanf();timebond=atof(sa)*3.15576e+7;
/**/
/* Errosion/Sedimentation */
ffscanf();erosmod=atoi(sa);
ffscanf();eroslev=atof(sa);
ffscanf();eroscon=atof(sa);
ffscanf();eroskoe=atof(sa);
ffscanf();sedilev=atof(sa);
ffscanf();sedicon=atof(sa);
ffscanf();sedikoe=atof(sa);
ffscanf();sedimcyc=atoi(sa);
ffscanf();waterlev=atof(sa);
ffscanf();slopemax=atof(sa);
/**/
/* V */
ffscanf();DIVVMIN=atof(sa);
ffscanf();STOKSMIN=atof(sa);
ffscanf();stoksmod=atoi(sa);
ffscanf();viscmod=atof(sa);
ffscanf();stoksfd=atoi(sa);
ffscanf();nubeg=atof(sa);
ffscanf();nuend=atof(sa);
ffscanf();nucontr=atof(sa);
ffscanf();hidry=atof(sa);
ffscanf();hidrl=atof(sa);
ffscanf();strmin=atof(sa);
ffscanf();strmax=atof(sa);
ffscanf();stredif=atof(sa);
ffscanf();plastmax=atof(sa);
/**/
/**/
/**/
/* T */
ffscanf();HEATMIN=atof(sa);
ffscanf();heatmod=atoi(sa);
ffscanf();heatfd=atoi(sa);
ffscanf();heatdif=atof(sa);
ffscanf();frictyn=atoi(sa);
ffscanf();adiabyn=atoi(sa);
/**/
/**/
/* Water */
ffscanf();tkpor=atof(sa);
ffscanf();zmpor=atof(sa);
ffscanf();vyfluid=atof(sa);
ffscanf();vymelt=atof(sa);
ffscanf();dmwamin=atof(sa);
ffscanf();tdeep=atof(sa);
ffscanf();zdeep=atof(sa);
ffscanf();dxwater=atof(sa);
ffscanf();dywater=atof(sa);
ffscanf();deserp=atof(sa);
ffscanf();dyserp=atof(sa);
ffscanf();lambfld=atof(sa);


/* ARMEL'S UPDATES */

	/* Cooling zone */
	ffscanf();coolzone=atof(sa);
	ffscanf();coolkoef=atof(sa);

/* END OF ARMEL'S UPDATES */

/*
printf("%e %e %e %e %e",tkpor,zmpor,vyfluid,vymelt,dmwamin);getchar();
*/
/**/
fclose(fl);
/* End Load information from configuration file mode.t3c */
/**/
/* stop.yn file creation */
fl = fopen("stop.yn","wt");
fprintf(fl,"n \n");
fclose(fl);
/**/
/*
return fln3;
*/
/**/
/* Load thermodynamic database */
/* Dry peridotite */
/* RO - density */
fl = fopen("pdry_rho","rt");
ffscanf();
ffscanf(); tknum=atoi(sa);
ffscanf(); pbnum=atoi(sa);
ffscanf(); tkmin=atof(sa);
ffscanf(); pbmin=atof(sa);
ffscanf(); tkstp=atof(sa);
ffscanf(); pbstp=atof(sa);
ffscanf();
ffscanf();
/*
printf("%d %d %e %e %e %e",tknum,pbnum,tkmin,pbmin,tkstp,pbstp);getchar();
*/
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][0][0]=atof(sa)/1000.0;
/*
printf("%d %d %e",n1,n2,td[n1][n2][0][0]);getchar();
*/
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("pdry_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][0][1]=atof(sa)/1000.0/4.1837;
/*
	ffscanf(); td[n2][n1][0][1]=atof(sa)/1000.0/4.1837/td[n2][n1][0][0]/1000.0;
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][0][1]-td[n2-1][n1][0][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("PDRY %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Vp - seismic velosity, km/s */
fl = fopen("pdry_vp","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][0][2]=MAXV(0,atof(sa));
/*
if (n1>160) {printf("%d %d %e",n1,n2,td[n2][n1][0][2]);getchar();}
*/
	}
fclose(fl);
/**/
/* Vs - seismic velosity, km/s */
fl = fopen("pdry_vs","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][0][3]=MAXV(0,atof(sa));
	}
fclose(fl);
printf("pdry_vs OK \n");
/**/
/* Wet peridotite */
/* RO - density */
fl = fopen("pwet_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][1][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("pwet_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][1][1]=atof(sa)/1000.0/4.1837;
/*
{printf("PWET %d %d %e %e",n1,n2,td[n2][n1][0][1],td[n2][n1][1][1]);getchar();}
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][1][1]-td[n2-1][n1][1][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("PWET %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Vp - seismic velosity, km/s */
fl = fopen("pwet_vp","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][1][2]=MAXV(0,atof(sa));
/*
if (n1>25 && n1<75 && n2>245 && n2<275) {printf("<%s> %d %d %e",sa,n1,n2,td[n2][n1][1][2]);getchar();}
*/
	}
fclose(fl);
/**/
/* Vs - seismic velosity, km/s */
fl = fopen("pwet_vs","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][1][3]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("pwet_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][1][4]=MAXV(0,atof(sa));
	}
fclose(fl);
printf("pwet_h2o OK \n");
/**/
/**/
/* Molten peridotite */
/* RO - density */
fl = fopen("pwetmelt_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][2][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("pwetmelt_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][2][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][2][1]-td[n2-1][n1][2][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("PMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Vp - seismic velosity, km/s */
fl = fopen("pwetmelt_vp","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][2][2]=MAXV(0,atof(sa));
/*
td[n2][n1][2][2]=atof(sa);
if (n1>25 && n1<75 && n2>245 && n2<275) {printf("<%s> %d %d %e",sa,n1,n2,td[n2][n1][2][2]);getchar();}
*/
	}
fclose(fl);
/**/
/* Vs - seismic velosity, km/s */
fl = fopen("pwetmelt_vs","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][2][3]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("pwetmelt_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][2][4]=MAXV(0,atof(sa));
	}
fclose(fl);
printf("pwetmelt_h2o OK \n");
/**/
/**/
/* Wet Gabbro */
/* RO - density */
fl = fopen("gabwet_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][3][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("gabwet_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][3][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][3][1]-td[n2-1][n1][3][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("GWET %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Vp - seismic velosity, km/s */
fl = fopen("gabwet_vp","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][3][2]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Vs - seismic velosity, km/s */
fl = fopen("gabwet_vs","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][3][3]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("gabwet_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][3][4]=MAXV(0,atof(sa));
	}
printf("gabwet_h2o OK \n");
fclose(fl);
/**/
/**/
/* Molten Gabbro */
/* RO - density */
fl = fopen("gabwetmelt_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][4][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("gabwetmelt_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][4][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][4][1]-td[n2-1][n1][4][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("GMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Vp - seismic velosity, km/s */
fl = fopen("gabwetmelt_vp","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][4][2]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Vs - seismic velosity, km/s */
fl = fopen("gabwetmelt_vs","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][4][3]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("gabwetmelt_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][4][4]=MAXV(0,atof(sa));
	}
fclose(fl);
printf("gabwetmelt_h2o OK \n");
/**/
/**/
/* Wet sediments */
/* RO - density */
fl = fopen("swet_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][5][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("swet_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][5][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][5][1]-td[n2-1][n1][5][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("SWET %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Vp - seismic velosity, km/s */
fl = fopen("swet_vp","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][5][2]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Vs - seismic velosity, km/s */
fl = fopen("swet_vs","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][5][3]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("swet_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][5][4]=MAXV(0,atof(sa));
	}
fclose(fl);
printf("swet_h2o OK \n");
/**/
/**/
/* Molten sediments */
/* RO - density */
fl = fopen("swetmelt_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][6][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("swetmelt_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][6][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][6][1]-td[n2-1][n1][6][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("SMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Vp - seismic velosity, km/s */
fl = fopen("swetmelt_vp","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][6][2]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Vs - seismic velosity, km/s */
fl = fopen("swetmelt_vs","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][6][3]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("swetmelt_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][6][4]=MAXV(0,atof(sa));
	}
fclose(fl);
printf("swetmelt_h2o OK \n");
/**/
/**/
/* Wet Basalt */
/* RO - density */
fl = fopen("bwet_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][7][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("bwet_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][7][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][7][1]-td[n2-1][n1][7][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("BWET %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Vp - seismic velosity, km/s */
fl = fopen("bwet_vp","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][7][2]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Vs - seismic velosity, km/s */
fl = fopen("bwet_vs","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][7][3]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("bwet_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][7][4]=MAXV(0,atof(sa));
	}
fclose(fl);
printf("bwet_h2o OK \n");
/**/
/**/
/* Molten Basalt */
/* RO - density */
fl = fopen("bwetmelt_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][8][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("bwetmelt_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][8][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][8][1]-td[n2-1][n1][8][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("BMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Vp - seismic velosity, km/s */
fl = fopen("bwetmelt_vp","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][8][2]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Vs - seismic velosity, km/s */
fl = fopen("bwetmelt_vs","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][8][3]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("bwetmelt_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][8][4]=MAXV(0,atof(sa));
	}
fclose(fl);
printf("bwetmelt_h2o OK \n");
/**/
/**/
/* Dry Upper crust */
/* RO - density */
fl = fopen("ucdry_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][11][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("ucdry_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][11][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][11][1]-td[n2-1][n1][11][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("BMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Vp - seismic velosity, km/s */
fl = fopen("ucdry_vp","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][8][2]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Vs - seismic velosity, km/s */
fl = fopen("ucdry_vs","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][11][3]=MAXV(0,atof(sa));
	}
fclose(fl);
printf("ucdry_vs OK \n");
/**/
/* Wet Upper crust */
/* RO - density */
fl = fopen("ucwet_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][12][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("ucwet_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][12][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][12][1]-td[n2-1][n1][12][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("BMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Vp - seismic velosity, km/s */
fl = fopen("ucwet_vp","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][12][2]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Vs - seismic velosity, km/s */
fl = fopen("ucwet_vs","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][12][3]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("ucwet_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][12][4]=MAXV(0,atof(sa));
	}
fclose(fl);
printf("ucwet_h2o OK \n");
/**/
/**/
/* Dry Lower crust */
/* RO - density */
fl = fopen("lcdry_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][13][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("lcdry_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][13][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][13][1]-td[n2-1][n1][13][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("BMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Vp - seismic velosity, km/s */
fl = fopen("lcdry_vp","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][13][2]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Vs - seismic velosity, km/s */
fl = fopen("lcdry_vs","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][13][3]=MAXV(0,atof(sa));
	}
fclose(fl);
printf("lcdry_vs OK \n");
/**/
/* Wet Lower crust */
/* RO - density */
fl = fopen("lcwet_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][14][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("lcwet_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][14][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][14][1]-td[n2-1][n1][14][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("BMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Vp - seismic velosity, km/s */
fl = fopen("lcwet_vp","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][14][2]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Vs - seismic velosity, km/s */
fl = fopen("lcwet_vs","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][14][3]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("lcwet_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][14][4]=MAXV(0,atof(sa));
	}
fclose(fl);
printf("lcwet_h2o OK \n");
/**/
/**/
/*
printf("%d %d %e %e %e %e ",n1,n2,td[n1-1][n2-1][0][0],td[n1-1][n2-1][0][1],td[n1-1][n2-1][0][2],td[n1-1][n2-1][0][3]);getchar();
*/
/**/
/* Load global thermodynamic database */
/* Dry peridotite Xmg=0.895 */
/* RO - density */
fl = fopen("m895_ro","rt");
ffscanf();
ffscanf(); tknum1=atoi(sa);
ffscanf(); pbnum1=atoi(sa);
ffscanf(); tkmin1=atof(sa);
ffscanf(); pbmin1=atof(sa);
ffscanf(); tkstp1=atof(sa);
ffscanf(); pbstp1=atof(sa);
ffscanf();
ffscanf();
for (n1=0;n1<pbnum1;n1++)
for (n2=0;n2<tknum1;n2++)
	{
	ffscanf(); td[n2][n1][9][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("m895_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum1;n1++)
for (n2=0;n2<tknum1;n2++)
	{
	ffscanf(); td[n2][n1][9][1]=atof(sa)/1000.0/4.1837;
	}
fclose(fl);
/**/
/* Vp - seismic velosity, km/s */
fl = fopen("m895_vp","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum1;n1++)
for (n2=0;n2<tknum1;n2++)
	{
	ffscanf(); td[n2][n1][9][2]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Vs - seismic velosity, km/s */
fl = fopen("m895_vs","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum1;n1++)
for (n2=0;n2<tknum1;n2++)
	{
	ffscanf(); td[n2][n1][9][3]=MAXV(0,atof(sa));
	}
fclose(fl);
printf("m895_vs OK \n");
/**/
/* Wet peridotite */
/* RO - density */
fl = fopen("morn_ro","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum1;n1++)
for (n2=0;n2<tknum1;n2++)
	{
	ffscanf(); td[n2][n1][10][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("morn_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum1;n1++)
for (n2=0;n2<tknum1;n2++)
	{
	ffscanf(); td[n2][n1][10][1]=atof(sa)/1000.0/4.1837;
	}
fclose(fl);
/**/
/* Vp - seismic velosity, km/s */
fl = fopen("morn_vp","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum1;n1++)
for (n2=0;n2<tknum1;n2++)
	{
	ffscanf(); td[n2][n1][10][2]=MAXV(0,atof(sa));
	}
fclose(fl);
/**/
/* Vs - seismic velosity, km/s */
fl = fopen("morn_vs","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum1;n1++)
for (n2=0;n2<tknum1;n2++)
	{
	ffscanf(); td[n2][n1][10][3]=MAXV(0,atof(sa));
	}
fclose(fl);
printf("morn_vs OK \n");
/**/
/**/
return fln3;
}
/* Load information from configuration file mode.t3c ============== */



/* Load Information from data file ------------------------------- */
void loader()
/* bondv[] - bondary value */
/* bondm[] - bondary mode 0=Not, -1=Value, 1,2...=LinNum+1 */
/* m1,m2 - node X,Y number */
{
/* Counter */
int n1;
char nn1;
long int m1,m2,m3;
long int mm1;
char szint,szlong,szfloat,szdouble,szcur;
float ival0;
double ival1;
/**/
/**/
/**/
/* Load Past Results from data file-------------------------------- */
if (printmod) printf("Load Past results from %s ...",fl1in);
/**/
/**/
/**/
/* Load in Text Format ---------------------------- */
if(fl1itp==0)
	{
	fl = fopen(fl1in,"rt");
	/**/
	/* Grid Parameters */
	ffscanf();xnumx=atoi(sa);
	ffscanf();ynumy=atoi(sa);
	ffscanf();mnumx=atoi(sa);
	ffscanf();mnumy=atoi(sa);
	ffscanf();marknum=atoi(sa);
	ffscanf();xsize=atof(sa);
	ffscanf();ysize=atof(sa);
	ffscanf();pinit=atof(sa);
	ffscanf();pkf[0]=atof(sa);
	ffscanf();pkf[1]=atof(sa);
	ffscanf();pkf[2]=atof(sa);
	ffscanf();pkf[3]=atof(sa);
	ffscanf();GXKOEF=atof(sa);
	ffscanf();GYKOEF=atof(sa);
	ffscanf();rocknum=atoi(sa);
	ffscanf();bondnum=atoi(sa);
	ffscanf();
	ffscanf();timesum=atof(sa)*3.15576e+7;
	/**/
	/* Calc,Check Grid parameters */
	gridcheck();
	/**/
	/* Rock Types information */
	for (n1=0;n1<rocknum;n1++)
		{
		ffscanf();
		ffscanf();markim[n1]=atoi(sa);;
		ffscanf();markn0[n1]=atof(sa);;
		ffscanf();markn1[n1]=atof(sa);;
		ffscanf();marks0[n1]=atof(sa);;
		ffscanf();marks1[n1]=atof(sa);;
		ffscanf();marknu[n1]=atof(sa);;
		ffscanf();markdh[n1]=atof(sa);;
		ffscanf();markdv[n1]=atof(sa);;
		ffscanf();markss[n1]=atof(sa);;
		ffscanf();markmm[n1]=atof(sa);;
		ffscanf();markgg[n1]=atof(sa);;
		ffscanf();markll[n1]=atof(sa);;
		ffscanf();marka0[n1]=atof(sa);;
		ffscanf();marka1[n1]=atof(sa);;
		ffscanf();markb0[n1]=atof(sa);;
		ffscanf();markb1[n1]=atof(sa);;
		ffscanf();marke0[n1]=atof(sa);;
		ffscanf();marke1[n1]=atof(sa);;
		ffscanf();markf0[n1]=atof(sa);;
		ffscanf();markf1[n1]=atof(sa);;
		ffscanf();markro[n1]=atof(sa);;
		ffscanf();markbb[n1]=atof(sa);;
		ffscanf();markaa[n1]=atof(sa);;
		ffscanf();markcp[n1]=atof(sa);;
		ffscanf();markkt[n1]=atof(sa);;
		ffscanf();markkf[n1]=atof(sa);;
		ffscanf();markkp[n1]=atof(sa);;
		ffscanf();markht[n1]=atof(sa);;
		}
	/**/
	/* Nodes information */
	/* Vx,Vy,bondm[],ro[],nu[],ep[],et[],tk[],cp[],kt[],ht[] */
	for (m1=0;m1<xnumx;m1++)
	for (m2=0;m2<ynumy;m2++)
		{
		m3=m1*ynumy+m2;
		ffscanf();
		ffscanf();
		ffscanf();pr[m3]=atof(sa);
		ffscanf();vx[m3]=atof(sa);
		ffscanf();vy[m3]=atof(sa);
		ffscanf();bondm[m3*3+0]=atoi(sa);
		ffscanf();bondm[m3*3+1]=atoi(sa);
		ffscanf();bondm[m3*3+2]=atoi(sa);
		ffscanf();exx[m3]=atof(sa);
		ffscanf();dro[m3]=atof(sa);
		ffscanf();exy[m3]=atof(sa);
		ffscanf();sxx[m3]=atof(sa);
		ffscanf();drp[m3]=atof(sa);
		ffscanf();sxy[m3]=atof(sa);
		ffscanf();ro[m3]=atof(sa);
		ffscanf();nu[m3]=atof(sa);
		ffscanf();nd[m3]=atof(sa);
		ffscanf();sxxe[m3]=atof(sa);
		ffscanf();sppe[m3]=atof(sa);
		ffscanf();sxye[m3]=atof(sa);
		ffscanf();exxe[m3]=atof(sa);
		ffscanf();exye[m3]=atof(sa);
		ffscanf();esp[m3]=atof(sa);
		ffscanf();mu[m3]=atof(sa);
		ffscanf();gg[m3]=atof(sa);
		ffscanf();gd[m3]=atof(sa);
		ffscanf();ep[m3]=atof(sa);
		ffscanf();et[m3]=atof(sa);
		ffscanf();tk[m3]=atof(sa);
		ffscanf();bondm[nodenum3+m3]=atoi(sa);
		ffscanf();cp[m3]=atof(sa);
		ffscanf();kt[m3]=atof(sa);
		ffscanf();ht[m3]=atof(sa);
		}
	/**/
	/* Gridlines positions */
	for (m1=0;m1<xnumx;m1++)
		{
		ffscanf();
		ffscanf();gx[m1]=atof(sa);
		}
	for (m2=0;m2<ynumy;m2++)
		{
		ffscanf();
		ffscanf();gy[m2]=atof(sa);
		}
	/**/
	/* Bondary Conditions Equations */
	for (m1=1;m1<bondnum;m1++)
		{
		ffscanf();
		ffscanf();bondv[m1][0]=atof(sa);
		ffscanf();bondv[m1][1]=atof(sa);
		ffscanf();bondv[m1][2]=atof(sa);
		ffscanf();bondv[m1][3]=atof(sa);
		ffscanf();bondn[m1][0]=atoi(sa);
		ffscanf();bondn[m1][1]=atoi(sa);
		ffscanf();bondn[m1][2]=atoi(sa);
		}
	/**/
	/* Markers X,Y,types */
	for (mm1=0;mm1<marknum;mm1++)
		{
		ffscanf();markx[mm1]=atof(sa);
		ffscanf();marky[mm1]=atof(sa);
		ffscanf();markk[mm1]=atof(sa);
		ffscanf();marke[mm1]=atof(sa);
		ffscanf();markxx[mm1]=atof(sa); 
		ffscanf();markv[mm1]=atof(sa); 
		ffscanf();markxy[mm1]=atof(sa); 
		ffscanf();markp[mm1]=atof(sa);
		ffscanf();markexx[mm1]=atof(sa);
		ffscanf();markexy[mm1]=atof(sa);
		ffscanf();markd[mm1]=atof(sa);
		ffscanf();markw[mm1]=atof(sa);
		ffscanf();markt[mm1]=atoi(sa); 
		}
	}
/* Load in Text Format ---------------------------- */
/**/
/**/
/**/
/* Load in Binary Format ---------------------------- */
else
	{
	fl = fopen(fl1in,"rb");
	/**/
	/* Sizes of var definition */
	szint=sizeof(n1);
	szlong=sizeof(m1);
	szfloat=sizeof(ival0);
	szdouble=sizeof(ival1);
	/* Check sizes of variables */
	fread(&szcur,1,1,fl);
	if (szcur!=szint) {printf("Current INT size <%d> is different from given in file <%d> \n",szint,szcur); exit(0);}
	fread(&szcur,1,1,fl);
	if (szcur!=szlong) {printf("Current LONG INT size <%d> is different from given in file <%d> \n",szlong,szcur); exit(0);}
	fread(&szcur,1,1,fl);
	if (szcur!=szfloat) {printf("Current FLOAT size <%d> is different from given in file <%d> \n",szfloat,szcur); exit(0);}
	fread(&szcur,1,1,fl);
	if (szcur!=szdouble) {printf("Current DOUBLE size <%d> is different from given in file <%d> \n",szdouble,szcur); exit(0);}
	/**/
	/* Grid Parameters */
	fread(&xnumx,szlong,1,fl);
	fread(&ynumy,szlong,1,fl);
	fread(&mnumx,szlong,1,fl);
	fread(&mnumy,szlong,1,fl);
	fread(&marknum,szlong,1,fl);
	fread(&xsize,szdouble,1,fl);
	fread(&ysize,szdouble,1,fl);
	fread(&pinit,szdouble,1,fl);
	fread(pkf,szdouble,4,fl);
	fread(&GXKOEF,szdouble,1,fl);
	fread(&GYKOEF,szdouble,1,fl);
	fread(&rocknum,szint,1,fl);
	fread(&bondnum,szlong,1,fl);
	fread(&n1,szint,1,fl);
	fread(&timesum,szdouble,1,fl);timesum*=3.15576e+7;
	/**/
	/* Calc,Check Grid parameters */
	gridcheck();
	/**/
	/* Rock Types information */
	fread(markim,szint,rocknum,fl);
	fread(markn0,szdouble,rocknum,fl);
	fread(markn1,szdouble,rocknum,fl);
	fread(marks0,szdouble,rocknum,fl);
	fread(marks1,szdouble,rocknum,fl);
	fread(marknu,szdouble,rocknum,fl);
	fread(markdh,szdouble,rocknum,fl);
	fread(markdv,szdouble,rocknum,fl);
	fread(markss,szdouble,rocknum,fl);
	fread(markmm,szdouble,rocknum,fl);
	fread(markgg,szdouble,rocknum,fl);
	fread(markll,szdouble,rocknum,fl);
	fread(marka0,szdouble,rocknum,fl);
	fread(marka1,szdouble,rocknum,fl);
	fread(markb0,szdouble,rocknum,fl);
	fread(markb1,szdouble,rocknum,fl);
	fread(marke0,szdouble,rocknum,fl);
	fread(marke1,szdouble,rocknum,fl);
	fread(markf0,szdouble,rocknum,fl);
	fread(markf1,szdouble,rocknum,fl);
	fread(markro,szdouble,rocknum,fl);
	fread(markbb,szdouble,rocknum,fl);
	fread(markaa,szdouble,rocknum,fl);
	fread(markcp,szdouble,rocknum,fl);
	fread(markkt,szdouble,rocknum,fl);
	fread(markkf,szdouble,rocknum,fl);
	fread(markkp,szdouble,rocknum,fl);
	fread(markht,szdouble,rocknum,fl);
	/**/
	/* Nodes information */
	/* Vx,Vy,bondm[],ro[],nu[],ep[],et[],tk[],cp[],kt[],ht[] */
	for (m1=0;m1<nodenum;m1++)
		{
		fread(&ival0,szfloat,1,fl);pr[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);vx[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);vy[m1]=(double)(ival0);
		fread(&m2,szlong,1,fl);bondm[m1*3+0]=m2;
		fread(&m2,szlong,1,fl);bondm[m1*3+1]=m2;
		fread(&m2,szlong,1,fl);bondm[m1*3+2]=m2;
		fread(&ival0,szfloat,1,fl);exx[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);dro[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);exy[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);sxx[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);drp[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);sxy[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);ro[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);nu[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);nd[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);sxxe[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);sppe[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);sxye[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);exxe[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);exye[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);esp[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);mu[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);gg[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);gd[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);ep[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);et[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);tk[m1]=(double)(ival0);
		fread(&m2,szlong,1,fl);bondm[nodenum3+m1]=m2;
		fread(&ival0,szfloat,1,fl);cp[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);kt[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);ht[m1]=(double)(ival0);
		}
	/**/
	/* Gridlines positions */
	for (m1=0;m1<xnumx;m1++)
		{
		fread(&ival0,szfloat,1,fl);gx[m1]=(double)(ival0);
/*
printf("%ld %e",m1,gx[m1]);getchar();
*/
		}
	for (m2=0;m2<ynumy;m2++)
		{
		fread(&ival0,szfloat,1,fl);gy[m2]=(double)(ival0);
/*
printf("%ld %e",m2,gy[m2]);getchar();
*/
		}
	/**/
	/* Bondary Conditions Equations */
	for (m1=1;m1<bondnum;m1++)
		{
		fread(&ival0,szfloat,1,fl);bondv[m1][0]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);bondv[m1][1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);bondv[m1][2]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);bondv[m1][3]=(double)(ival0);
		fread(&m2,szlong,1,fl);bondn[m1][0]=m2;
		fread(&m2,szlong,1,fl);bondn[m1][1]=m2;
		fread(&m2,szlong,1,fl);bondn[m1][2]=m2;
		}
	/**/
	/* Markers X,Y,types */
	for (mm1=0;mm1<marknum;mm1++)
		{
		fread(&ival0,szfloat,1,fl);markx[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);marky[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markk[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);marke[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markxx[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markv[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markxy[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markp[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markexx[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markexy[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markd[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markw[mm1]=ival0;
		fread(&nn1,1,1,fl);markt[mm1]=nn1;
		}
	}
/* Load in Binary Format ---------------------------- */
/**/
/**/
/**/
fclose(fl);
if (printmod) printf("OK!\n");
}
/* Load Information from data file ------------------------------- */



/* Print Results to data file ----------------------------------- */
void saver(int f0, int n0)
/* n0 - circle number */
{
/* Counters */
int n1,mm2;
char nn1;
long int m1,m2,m3;
long int mm1;
/* Buffers */
char szint,szlong,szfloat,szdouble;
float ival0;
double ival1;
/**/
/**/
/**/
if (printmod) printf("Print %d circle results to %s...",n0+1,fl1out);
/**/
/**/
/**/
/* Save data in text format ---------------------------- */
if (fl1otp==0)
	{
	fl = fopen(fl1out,"wt");
	/**/
	/* Grid Parameters */
	fprintf(fl,"%ld-xnumx\n",xnumx);
	fprintf(fl,"%ld-ynumy\n",ynumy);
	fprintf(fl,"%ld-mnumx\n",mnumx);
	fprintf(fl,"%ld-mnumy\n",mnumy);
	fprintf(fl,"%ld-marknum\n",marknum);
	fprintf(fl,"% 9.8e-xsize\n",xsize);
	fprintf(fl,"% 9.8e-ysize\n",ysize);
	fprintf(fl,"% 9.8e-pinit\n",pinit);
	fprintf(fl,"% 9.8e-p0\n",pkf[0]);
	fprintf(fl,"% 9.8e-px\n",pkf[1]);
	fprintf(fl,"% 9.8e-py\n",pkf[2]);
	fprintf(fl,"% 9.8e-pxy\n",pkf[3]);
	fprintf(fl,"% 9.8e-GXKOEF\n",GXKOEF);
	fprintf(fl,"% 9.8e-GYKOEF\n",GYKOEF);
	fprintf(fl,"%d-rocknum\n",rocknum);
	fprintf(fl,"%ld-bondnum\n",bondnum);
	fprintf(fl,"%d-n0cycle\n",n0+1);
	fprintf(fl,"% 9.8e-timesum",timesum/3.15576e+7);
	fprintf(fl,"\n\n\n");
	/**/
	/* Rock Types information */
	for (n1=0;n1<rocknum;n1++)
		{
		fprintf(fl,"% 3d % 1d % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e\n",n1,markim[n1],markn0[n1],markn1[n1],marks0[n1],marks1[n1],marknu[n1],markdh[n1],markdv[n1],markss[n1],markmm[n1],markgg[n1],markll[n1],marka0[n1],marka1[n1],markb0[n1],markb1[n1],marke0[n1],marke1[n1],markf0[n1],markf1[m1],markro[n1],markbb[n1],markaa[n1],markcp[n1],markkt[n1],markkf[n1],markkp[n1],markht[n1]);
		}
	fprintf(fl,"\n\n\n");
	/**/
	/* Nodes information */
	/* Vx,Vy,bondm[],ro[],nu[],ep[],et[],tk[],cp[],kt[],ht[] */
	for (m1=0;m1<xnumx;m1++)
		{
		for (m2=0;m2<ynumy;m2++)
			{
			m3=m1*ynumy+m2;
			/* X, Y */
			fprintf(fl,"% 5.4e % 5.4e ",gx[m1],gy[m2]);
			/* P,Vx,Vy, Bond P,Vx,Vy */
			fprintf(fl,"% 9.8e % 9.8e % 9.8e %ld %ld %ld ",pr[m3],vx[m3],vy[m3],bondm[m3*3+0],bondm[m3*3+1],bondm[m3*3+2]);
			/* EPS, SIG */
			fprintf(fl,"% 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e ",exx[m3],dro[m3],exy[m3],sxx[m3],drp[m3],sxy[m3]);
			/* ro, Nu koef */
			fprintf(fl,"% 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e ",ro[m3],nu[m3],nd[m3],sxxe[m3],sppe[m3],sxye[m3],exxe[m3],exye[m3],esp[m3],mu[m3],gg[m3],gd[m3],ep[m3],et[m3]);
			/* T, Bond T, Cp, Kt, Ht */
			fprintf(fl,"% 9.8e %ld % 9.8e % 9.8e % 9.8e\n",tk[m3],bondm[nodenum3+m3],cp[m3],kt[m3],ht[m3]);
			}
		fprintf(fl,"\n");
		}
	fprintf(fl,"\n\n\n");
	/**/
	/* Gridlines positions */
	for (m1=0;m1<xnumx;m1++)
		{
		fprintf(fl,"%ld % 9.8e\n",m1,gx[m1]);
		}
	fprintf(fl,"\n");
	for (m2=0;m2<ynumy;m2++)
		{
		fprintf(fl,"%ld % 9.8e\n",m2,gy[m2]);
		}
	fprintf(fl,"\n\n\n");
	/**/
	/* Bondary Conditions Equations */
	for (m1=1;m1<bondnum;m1++)
		{
		fprintf(fl,"%6ld % 9.8e % 9.8e % 9.8e % 9.8e %ld %ld %ld \n",m1,bondv[m1][0],bondv[m1][1],bondv[m1][2],bondv[m1][3],bondn[m1][0],bondn[m1][1],bondn[m1][2]);
		}
	fprintf(fl,"\n\n\n");
	/**/
	/* Markers X,Y,Type */
	for (m1=0;m1<marknum;m1++)
		{
		mm2=markt[m1];
		fprintf(fl,"% 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e % 9.8e  % 9.8e % 9.8e % 9.8e %d \n",markx[m1],marky[m1],markk[m1],marke[m1],markxx[m1],markv[m1],markxy[m1],markp[m1],markexx[m1],markexy[m1],markd[m1],markw[m1],mm2);
		}
	}
/* Save data in text format ---------------------------- */
/**/
/**/
/**/
/* Save data in binary format ---------------------------- */
else
	{
	fl = fopen(fl1out,"wb");
	/**/
	/* Sizes of var definition */
	szint=sizeof(n1);
	szlong=sizeof(m1);
	szfloat=sizeof(ival0);
	szdouble=sizeof(ival1);
	fwrite(&szint,1,1,fl);
	fwrite(&szlong,1,1,fl);
	fwrite(&szfloat,1,1,fl);
	fwrite(&szdouble,1,1,fl);
	/**/
	/* Grid Parameters */
	fwrite(&xnumx,szlong,1,fl);
	fwrite(&ynumy,szlong,1,fl);
	fwrite(&mnumx,szlong,1,fl);
	fwrite(&mnumy,szlong,1,fl);
	fwrite(&marknum,szlong,1,fl);
	fwrite(&xsize,szdouble,1,fl);
	fwrite(&ysize,szdouble,1,fl);
	fwrite(&pinit,szdouble,1,fl);
	fwrite(pkf,szdouble,4,fl);
	fwrite(&GXKOEF,szdouble,1,fl);
	fwrite(&GYKOEF,szdouble,1,fl);
	fwrite(&rocknum,szint,1,fl);
	fwrite(&bondnum,szlong,1,fl);
	fwrite(&n0,szint,1,fl);
	ival1=timesum/3.15576e+7;fwrite(&ival1,szdouble,1,fl);
	/**/
	/* Rock Types information */
	fwrite(markim,szint,rocknum,fl);
	fwrite(markn0,szdouble,rocknum,fl);
	fwrite(markn1,szdouble,rocknum,fl);
	fwrite(marks0,szdouble,rocknum,fl);
	fwrite(marks1,szdouble,rocknum,fl);
	fwrite(marknu,szdouble,rocknum,fl);
	fwrite(markdh,szdouble,rocknum,fl);
	fwrite(markdv,szdouble,rocknum,fl);
	fwrite(markss,szdouble,rocknum,fl);
	fwrite(markmm,szdouble,rocknum,fl);
	fwrite(markgg,szdouble,rocknum,fl);
	fwrite(markll,szdouble,rocknum,fl);
	fwrite(marka0,szdouble,rocknum,fl);
	fwrite(marka1,szdouble,rocknum,fl);
	fwrite(markb0,szdouble,rocknum,fl);
	fwrite(markb1,szdouble,rocknum,fl);
	fwrite(marke0,szdouble,rocknum,fl);
	fwrite(marke1,szdouble,rocknum,fl);
	fwrite(markf0,szdouble,rocknum,fl);
	fwrite(markf1,szdouble,rocknum,fl);
	fwrite(markro,szdouble,rocknum,fl);
	fwrite(markbb,szdouble,rocknum,fl);
	fwrite(markaa,szdouble,rocknum,fl);
	fwrite(markcp,szdouble,rocknum,fl);
	fwrite(markkt,szdouble,rocknum,fl);
	fwrite(markkf,szdouble,rocknum,fl);
	fwrite(markkp,szdouble,rocknum,fl);
	fwrite(markht,szdouble,rocknum,fl);
	/**/
	/* Nodes information */
	/* Vx,Vy,bondm[],ro[],nu[],ep[],et[],tk[],cp[],kt[],ht[] */
	for (m1=0;m1<nodenum;m1++)
		{
		ival0=(float)(pr[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(vx[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(vy[m1]);fwrite(&ival0,szfloat,1,fl);
		m2=bondm[m1*3+0];fwrite(&m2,szlong,1,fl);
		m2=bondm[m1*3+1];fwrite(&m2,szlong,1,fl);
		m2=bondm[m1*3+2];fwrite(&m2,szlong,1,fl);
		ival0=(float)(exx[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(dro[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(exy[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(sxx[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(drp[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(sxy[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(ro[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(nu[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(nd[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(sxxe[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(sppe[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(sxye[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(exxe[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(exye[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(esp[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(mu[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(gg[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(gd[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(ep[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(et[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(tk[m1]);fwrite(&ival0,szfloat,1,fl);
		m2=bondm[nodenum3+m1];fwrite(&m2,szlong,1,fl);
		ival0=(float)(cp[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(kt[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(ht[m1]);fwrite(&ival0,szfloat,1,fl);
		}
	/**/
	/* Gridlines positions */
	for (m1=0;m1<xnumx;m1++)
		{
		ival0=(float)(gx[m1]);fwrite(&ival0,szfloat,1,fl);
		}
	for (m2=0;m2<ynumy;m2++)
		{
		ival0=(float)(gy[m2]);fwrite(&ival0,szfloat,1,fl);
		}
	/**/
	/* Bondary Conditions Equations */
	for (m1=1;m1<bondnum;m1++)
		{
		ival0=(float)(bondv[m1][0]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(bondv[m1][1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(bondv[m1][2]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(bondv[m1][3]);fwrite(&ival0,szfloat,1,fl);
		m2=bondn[m1][0];fwrite(&m2,szlong,1,fl);
		m2=bondn[m1][1];fwrite(&m2,szlong,1,fl);
		m2=bondn[m1][2];fwrite(&m2,szlong,1,fl);
		}
	/**/
	/* Markers X,Y,types */
	for (mm1=0;mm1<marknum;mm1++)
		{
		ival0=markx[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=marky[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markk[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=marke[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markxx[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markv[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markxy[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markp[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markexx[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markexy[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markd[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markw[mm1];fwrite(&ival0,szfloat,1,fl);
		nn1=markt[mm1];fwrite(&nn1,1,1,fl);
		}
	}
/* Save data in binary format ---------------------------- */
/**/
/**/
/**/
fclose(fl);
if (printmod) printf("OK!\n");
/**/
/**/
/* file.t3c file creation */
fl = fopen("file.t3c","wt");
fprintf(fl,"%d \n",f0);
fclose(fl);
/**/
/* stop.yn file information read */
fl = fopen("stop.yn","rt");
/* Read String */
ffscanf();
/**/
/* Stop Y/N */
if (sa[0]=='y' || sa[0]=='Y')
	{
	fclose(fl);
	printf("PROGRAM TERMINATED FROM stop.yn \n");
	exit(0);
	}
/**/
/* Change printmod */
if (sa[0]>='0' && sa[0]<='9')
	{
	printmod=atoi(sa);
	}
fclose(fl);
}
/* Print Results to data file ----------------------------------- */



/* LOAD WITHOUT EMPTY LINES from fl =================================== */
void ffscanf()
{
/* Counter */
int n1;
/**/
/* Read cycle */
do
	{
	/* Input string */
	n1=fscanf(fl,"%s",sa);
	/* Check end of file */
	if (n1<1)
		{
		printf("\n Unexpected end of file\n");
		fclose(fl);
		exit(0);
		}
	/* Delete last symbol <32 */
	for(n1=strlen(sa)-1;n1>=0;n1--)
	if (*(sa+n1)<=32)
	*(sa+n1)=0;
	else
	break;
	}
while (strlen(sa)<1 || sa[0]=='/');
}
/* End LOAD WITHOUT EMPTY LINES from fl =================================== */



/* LOAD WITHOUT EMPTY LINES from fl1 =================================== */
void ffscanf1()
{
/* Counter */
int n1;
/**/
/* Read cycle */
do
	{
	/* Input string */
	n1=fscanf(fl1,"%s",sa);
	/* Check end of file */
	if (n1<1)
		{
		printf("\n Unexpected end of file\n");
		fclose(fl1);
		exit(0);
		}
	/* Delete last symbol <32 */
	for(n1=strlen(sa)-1;n1>=0;n1--)
	if (*(sa+n1)<=32)
	*(sa+n1)=0;
	else
	break;
	}
while (strlen(sa)<1 || sa[0]=='/');
}
/* End LOAD WITHOUT EMPTY LINES from fl1 =================================== */


/* Calc,Check Parameters of Grid */
void gridcheck()
{
/* Gridlines NUM */
if(xnumx>MAXXLN) {printf("Space out in gx[] %ld",xnumx); exit(0);}
if(ynumy>MAXYLN) {printf("Space out in gy[] %ld",ynumy); exit(0);}
/**/
/* Nodes Num */
nodenum=xnumx*ynumy;
if(nodenum>MAXNOD) {printf("Space out in vx[],vy[] %ld",nodenum); exit(0);}
/**/
/* Cells Num */
cellnum=(xnumx-1)*(ynumy-1);
if(cellnum>MAXCEL) {printf("Space out in pr[]"); exit(0);}
/**/
/* Mark num */
if(marknum>MAXMRK+1) {printf("Space out in markx[]"); exit(0);}
/**/
/* Rock types Num */
if (rocknum>MAXTMR){printf("Space out in marknu[]"); exit(0);}
/**/
/* Bondary condit Equations Num */
if (bondnum>MAXBON){printf("Space out in bondv[]"); exit(0);}
/**/
/* Koef for processing */
xstpx=xsize/(double)(xnumx-1);
ystpy=ysize/(double)(ynumy-1);
kfx=1.0/xstpx;
kfy=1.0/ystpy;
kfxx=kfx*kfx;
kfyy=kfy*kfy;
kfxy=kfx*kfy;
/* Marker size */
mardx=xstpx/(double)(mnumx);
mardy=ystpy/(double)(mnumy);
/* Step for numerical differentiation */
numdx=5e-1*mardx;
numdy=5e-1*mardy;
/**/
/* Spec counters */
nodenum2=nodenum*2;
nodenum3=nodenum*3;
xnumx1=xnumx-1;
ynumy1=ynumy-1;
/**/
}
/* Calc,Check Parameters of Grid */
