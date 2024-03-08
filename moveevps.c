/* LEFT+Right Side or Err for X-Stokes Equation */
/* Stoks equation initial form */
/* dSIGxx/dX + dSIGxy/dY - dP/dX = - RO*Gx */
double xstokserr(long int m1, long int m2, int ynerr)
/* m1,m2 - node X,Y number */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Counters */
int n1,n2;
long int v[4];
/* Err Buf */
double leftx,rightx,nueff;
/* Distances */
double xkf=(gx[m1+1]-gx[m1-1])/2.0,ykf=gy[m2+1]-gy[m2];
/**/
/**/
/**/
/* Staggered Nodes num */
/*                      [0]                [2]   */
/*                  Ro0,Sxy0,Nu0                 */
/*                                               */
/*         Pr1,Sxx1    <Vx0>  Pr3,Sxx3           */
/*                                               */
/*                      [1]                [3]   */
/*                  RO1,Sxy1,Nu1                 */
/*                                               */
v[0]=m1*ynumy+m2;v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/**/
/* RIGHT parts of X-Stokes */
/* dSIGik/dXk - dP/dXi = - Gi*RO */
rightx  = -GXKOEF*(ro[v[0]]+ro[v[1]])/2.0;
/**/
/**/
/**/
/* Return val for LSQ err ----------------------------*/
if(ynerr==1)
	{
	/* LEFT part of X-Stokes */
	/* dSIGxx/dX + dSIGxy/dY - dP/dX = - RO*Gx */
	/**/
	/* dSIGxx/dX */
	leftx =(sxx[v[3]]-sxx[v[1]])/xkf;
	/* dSIGxy/dY */
	leftx+=(sxy[v[1]]-sxy[v[0]])/ykf;
	/* -dP/dX */
	leftx-=(pr[v[3]]-pr[v[1]])/xkf;
	/**/
	/* Effective NU calc */
	nueff=MAXV(ABSV(nu[v[0]]),ABSV(nu[v[1]]));
	/**/
	/* X-STOKS Error */
	leftx=(leftx-rightx)/nueff;
	/**/
	/* Min,Max Value of P Save */
	errbuf[10]=MAXV(pr[v[1]],pr[v[3]]);
	errbuf[11]=MINV(pr[v[1]],pr[v[3]]);
	/**/
	return leftx;
	}
/**/
/**/
/**/
/* Set Initial Num of lines -------------------------------------- */
wn[0]=2;
/* Save Right part Save for X-Stokes ---------------------*/
wi[0]=rightx;
/**/
/* Add Coefficients for left parts of X-Stokes ----------------*/
/* Staggered Nodes num */
/*                      [0]                [2]   */
/*                  Ro0,Sxy0,Nu0                 */
/*                                               */
/*         Pr1,Sxx1    <Vx0>  Pr3,Sxx3           */
/*                                               */
/*                      [1]                [3]   */
/*                  RO1,Sxy1,Nu1                 */
/*                                               */
/*  0(P) 1(Vx)  2(Vy)  */
/* -dP/dX */
wn[1]=v[1]*3+0;
wi[1]=+1.0/xkf;
wn[2]=v[3]*3+0;
wi[2]=-1.0/xkf;
/* dSIGxx/dX */
sxxcalc(m1,m2+1,-1.0/xkf);
sxxcalc(m1+1,m2+1,1.0/xkf);
/* dSIGxy/dY */
sxycalc(m1,m2,-1.0/ykf);
sxycalc(m1,m2+1,1.0/ykf);
/**/
/**/
/**/
/*
for(n1=0;n1<=vn[0];n1++)printf("Vx %e %ld \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left+Right Side or Err for X-Stokes Equation */




/* LEFT+Right Side or Err for Y-Stokes Equation */
/* Stoks equation initial form */
/* dSIGyy/dY + dSIGxy/dX - dP/dY = - RO*Gy */
double ystokserr(long int m1, long int m2, int ynerr)
/* m1,m2 - node X,Y number */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Counters */
int n1;
long int v[4];
/* Err Buf */
double lefty,righty,nueff;
/* Distances */
double xkf=gx[m1+1]-gx[m1],ykf=(gy[m2+1]-gy[m2-1])/2.0;
/**/
/**/
/**/
/* Staggered Nodes num */
/*                                               */
/*                Pr2,Syy2                       */
/*                                               */
/*         [0]                  [2]              */
/*    Ro0,Sxy0,Nu0   <Vy0>   Ro2,Sxy2,Nu2        */
/*                                               */
/*                Pr3,Syy3                       */
/*                                               */
/*         [1]                  [3]              */
/*                                               */
v[0]=m1*ynumy+m2;v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/**/
/* RIGHT part of Y-Stokes */
/* dSIGik/dXk - dP/dXi = - Gi*RO */
righty  = -GYKOEF*(ro[v[0]]+ro[v[2]])/2.0;
/**/
/**/
/**/
/**/
/* Return val for LSQ err ----------------------------*/
if(ynerr==1)
	{
	/* LEFT part of Y-Stokes */
	/* dSIGyy/dY + dSIGxy/dX - dP/dY = - RO*Gy */
	/**/
	/* dSIGyy/dY */
	lefty =(-sxx[v[3]]+sxx[v[2]])/ykf;
	/* dSIGxy/dX */
	lefty+=(sxy[v[2]]-sxy[v[0]])/xkf;
	/* -dP/dY */
	lefty-=(pr[v[3]]-pr[v[2]])/ykf;
	/**/
	/* Effective NU calc */
	nueff=MAXV(ABSV(nu[v[0]]),ABSV(nu[v[2]]));
	/**/
	/* Y-STOKS Error */
	lefty=(lefty-righty)/nueff;
	/**/
	/* Min,Max Value of P Save */
	errbuf[10]=MAXV(pr[v[2]],pr[v[3]]);
	errbuf[11]=MINV(pr[v[2]],pr[v[3]]);
	/**/
	return lefty;
	}
/**/
/**/
/**/
/* Set Initial Num of lines -------------------------------------- */
wn[0]=2;
/* Save Right parts Save for Y-Stokes ---------------------*/
wi[0]=righty;
/**/
/* Add Coefficients for left parts of Y-Stokes ----------------*/
/* Staggered Nodes num */
/*                                               */
/*                Pr2,Syy2                       */
/*                                               */
/*         [0]                  [2]              */
/*    Ro0,Sxy0,Nu0   <Vy0>   Ro2,Sxy2,Nu2        */
/*                                               */
/*                Pr3,Syy3                       */
/*                                               */
/*         [1]                  [3]              */
/*                                               */
/*  0(P) 1(Vx)  2(Vy)  */
/* -dP/dY */
wn[1]=v[2]*3+0;
wi[1]=+1.0/ykf;
wn[2]=v[3]*3+0;
wi[2]=-1.0/ykf;
/* dSIGyy/dY */
sxxcalc(m1+1,m2,1.0/ykf);
sxxcalc(m1+1,m2+1,-1.0/ykf);
/* dSIGxy/dX */
sxycalc(m1,m2,-1.0/xkf);
sxycalc(m1+1,m2,1.0/xkf);
/**/
/**/
/**/
/*
for(n1=0;n1<=vn[0];n1++)printf("Vy %e %ld \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left+Right Side or Err for Y-Stokes Equation */




/* Left side or Err for vX Boundary Condition Equation */ 
/* V=CONST+KOEF*Vn */
double xbonderr(long int mcmax, int ynerr)
/* mcmax - numer of cur Vx in sol[] */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Val Buffer */
double leftx=0;
int n1;
/**/
/**/
/**/
/* Error Calc */
if (ynerr)
	{
	/* Add Const */
	leftx=x[mcmax]-bondv[bondm[mcmax]][0];
	/* Add Koef */
	for (n1=0;n1<3;n1++)
		{
		if(bondn[bondm[mcmax]][n1]) 
			{
			leftx-=bondv[bondm[mcmax]][n1+1]*x[bondn[bondm[mcmax]][n1]-1];
			}
		}
	/**/
	return leftx;
	}
/**/
/**/
/**/
/* Add X CONST */
wn[0]=1;
wi[0]=bondv[bondm[mcmax]][0];
wn[1]=mcmax;
wi[1]=1.0;
/* Add X PAR1,PAR2,PAR3 */
for (n1=0;n1<3;n1++)
	{
	if(bondn[bondm[mcmax]][n1]) 
		{
		wn[0]+=1;
		wn[wn[0]]=bondn[bondm[mcmax]][n1]-1;
		wi[wn[0]]=-bondv[bondm[mcmax]][n1+1];
		}
	}
/**/
/*
for(n1=0;n1<3;n1++)printf("%e %ld \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Err for vX Boundary Condition Equation */ 





/* Left side or Err for vY Boundary Condition Equation */ 
/* V=CONST+KOEF*Vn */
double ybonderr(long int mcmax, int ynerr)
/* mcmax - numer of cur Vx in sol[] */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Val Buffer */
double lefty=0;
int n1;
/**/
/**/
/**/
/* Error Calc */
if (ynerr)
	{
	/* Add Const */
	lefty=x[mcmax]-bondv[bondm[mcmax]][0];
	/* Add Koef */
	for (n1=0;n1<3;n1++)
		{
		if(bondn[bondm[mcmax]][n1]) 
			{
			lefty-=bondv[bondm[mcmax]][n1+1]*x[bondn[bondm[mcmax]][n1]-1];
			}
		}
	/**/
	return lefty;
	}
/**/
/**/
/**/
/* Add Y CONST */
wn[0]=1;
wi[0]=bondv[bondm[mcmax]][0];
wn[1]=mcmax;
wi[1]=1.0;
/* Add Y PAR1,PAR2,PAR3 */
for (n1=0;n1<3;n1++)
	{
	if(bondn[bondm[mcmax]][n1]) 
		{
		wn[0]+=1;
		wn[wn[0]]=bondn[bondm[mcmax]][n1]-1;
		wi[wn[0]]=-bondv[bondm[mcmax]][n1+1];
		}
	}
/*
for(n1=0;n1<3;n1++)printf("%e %d \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Err for vY Boundary Condition Equation */ 




/* Left side or Err for P Boundary Condition Equation */ 
/* P=CONST+KOEF*Pn */
double pbonderr(long int mcmax, int ynerr)
/* mcmax - numer of cur P  in sol[] */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Val Buffer */
double leftp;
int n1;
/**/
/**/
/**/
/* Error Calc */
if (ynerr)
	{
	/* Add Const */
	leftp=x[mcmax]-bondv[bondm[mcmax]][0];
	/* Add Koef */
	for (n1=0;n1<3;n1++)
		{
		if(bondn[bondm[mcmax]][n1]) 
			{
			leftp-=bondv[bondm[mcmax]][n1+1]*x[bondn[bondm[mcmax]][n1]-1];
			}
		}
	/**/
	return leftp*leftp;
	}
/**/
/**/
/**/
/* Add P CONST */
wn[0]=1;
wi[0]=bondv[bondm[mcmax]][0];
wn[1]=mcmax;
wi[1]=1.0;
/* Add P PAR1,PAR2,PAR3 */
for (n1=0;n1<3;n1++)
	{
	if(bondn[bondm[mcmax]][n1]) 
		{
		wn[0]+=1;
		wn[wn[0]]=bondn[bondm[mcmax]][n1]-1;
		wi[wn[0]]=-bondv[bondm[mcmax]][n1+1];
		}
	}
/*
for(n1=0;n1<3;n1++)printf("%e %d \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Err for P Boundary Condition Equation */ 





/* Solve XY-Stokes+Continuity equations by vX,vY,P mixed arbitrary order Finite Diff method */
void viterate(int m0)
{
/* Counters */
long int m1,m2,m3,m4,m5,m1min,m1max,m2min,m2max,mcmax,mcmax0,mcmax1,dm1,dm2,ccc=0;
double vxmin,vxmax,vymin,vymax,minvx,maxvx,minvy,maxvy,mindx,pmpa,minnu,maxnu,minpr,maxpr,mukoef;
double e,maxds,dss1,dsxx,dsxx1,dspp,dspp1,dsxy,dsxy1,dexx,dexx1,dexy,dexy1,celdx,celdy,swt1,mgg,mnu,xelvis0=1.0,xelvis1=0;
int printyn=printmod,n1;
/**/
/* Val buffer */
double ival;
/* Err koef */
double bondsum,bondnum;
double stoksum,stoknum;
double contsum,contnum;
/**/
/**/
/**/
/* Change Grid */
if(1==0)
{
/* Save old erosion/sedimentation, hydration surfaces */
for (m1=0;m1<xnumx;m1++)
	{
	ep0[m1]=ep[m1];
	ep0[xnumx+m1]=ep[xnumx+m1];
	ep0[xnumx*2+m1]=ep[xnumx*2+m1];
	ep0[xnumx*3+m1]=ep[xnumx*3+m1];
	}
/* Recomputing horizontal grid */
ival=gx[300];
for (m3=0;m3<=marknum;m3++) 
/* Check markers out of grid */
if (markx[m3]>0 && marky[m3]>0 && (double)(markx[m3])<xsize && (double)(marky[m3])<ysize && markt[m3]<50)
	{
	/* Check Asthenosphere position */
	if ((markt[m3]==10 || markt[m3]==14 || markt[m3]==34 || markt[m3]==12) && markx[m3]>gx[300] && markx[m3]<gx[301] && marky[m3]>3e+4 && marky[m3]<6e+4 && markx[m3]>ival) ival=markx[m3];
	}
if (printmod) printf("Trench position = %e km      step1 = %e km    step2 = %e km \n",ival/1000.0,(ival-3e+5)/150.0/1000.0,(3e+6-ival-7e+5)/150.0/1000.0);
ep0[xnumx*4]=0;
for (m1=1;m1<=150;m1++) 
	{
	if (ival<6e+5) 
		{
		ep0[xnumx*4+m1]=ep0[xnumx*4+m1-1]+(ival-3e+5)/150.0;
		}
	else
		{
		ep0[xnumx*4+m1]=ep0[xnumx*4+m1-1]+2e+3+(ival-3e+5-2e+3*150.0)/(150.0*151.0/2.0)*(151.0-(double)(m1));
		}
	}
for (m1=151;m1<=650;m1++) ep0[xnumx*4+m1]=ep0[xnumx*4+m1-1]+2e+3;
for (m1=651;m1<=800;m1++) ep0[xnumx*4+m1]=ep0[xnumx*4+m1-1]+2e+3+(3e+6-ival-7e+5-2e+3*150.0)/(150.0*151.0/2.0)*((double)(m1)-650.0);
/* 
printf("%e %e %e %e",ep0[xnumx*4+150],ep0[xnumx*4+650],ep0[xnumx*4+800],ep0[xnumx*4+800]-ep0[xnumx*4+799]);getchar();
*/
/* Reinterpolate erosion/sedimentation, hydration surfaces */
for (m3=0;m3<xnumx;m3++)
	{
	m1=m1serch(ep0[xnumx*4+m3]);
	/* Relativ Normalized coord Calc */
	e=(ep0[xnumx*4+m3]-gx[m1])/(gx[m1+1]-gx[m1]);
	/* Surface level elevation for marker definition */
	ep[m3]=(e*ep0[m1+1]+(1.0-e)*ep0[m1]);
	ep[xnumx+m3]=(e*ep0[xnumx+m1+1]+(1.0-e)*ep0[xnumx+m1]);
	ep[xnumx*2+m3]=(e*ep0[xnumx*2+m1+1]+(1.0-e)*ep0[xnumx*2+m1]);
	ep[xnumx*3+m3]=(e*ep0[xnumx*3+m1+1]+(1.0-e)*ep0[xnumx*3+m1]);
	}
/* Reload new gridline positions */
for (m1=0;m1<xnumx;m1++)
	{
	gx[m1]=ep0[xnumx*4+m1];
	}
}
/* End Change Grid */
/**/
/**/
/**/
/* Recomputing properties with new timestep */
if (stoksmod)
	{
	if(plastmax>0)
		{
		plastmark=plastold=plastnew=plastnew1=0;
		/* Check plastic yielding */
		plastic();
		if (printmod) printf("Elastic timestep %e  Plastic markers = %e  Max starting %e   Actual starting %e       stable  %e continuing %e \n",timestepe/3.15576e+7,plastmark,plastmark*plastmax,plastnew,plastnew1,plastold);
		/* Decrease amount of markers starting yielding */
		if(plastmark>0 && plastnew/plastmark>plastmax)
			{
			timestepe*=plastmax/(plastnew/plastmark);
			if(timestep>timestepe) timestep=timestepe;
			}
		}
	/* Recalc properties for markers Y/N */
	plastmark=plastold=plastnew=plastnew1=0;
	if (printmod) printf("\n RO, NU, CP etc  RECALC...");
	ronurecalc();
	if (printmod) printf(" RO, NU, CP etc OK!\n");
	if (printmod) printf("Elastic timestep %e  Plastic markers = %e  Max starting %e   Actual starting %e       stable  %e continuing %e \n",timestepe/3.15576e+7,plastmark,plastmark*plastmax,plastnew,plastnew1,plastold);
	}
/**/
/**/
/**/
if(stoksmod)
{
/* Movement Timestep check   */
if (timestep<=0)
	{
	printf("EXIT PROGRAM:  Movement timestep<=0 <%e>",timestep);
	exit(0);
	}
/* Elastic Timestep check   */
if (timestepe<=0)
	{
	printf("EXIT PROGRAM:  Elastic timestep<=0 <%e>",timestepe);
	exit(0);
	}
/**/
/**/
/**/
/* Check  Xelvis  */	
/* Change stoksmod for Xelvis calc mode */	
stoksmod=-stoksmod;
for (m1=1;m1<xnumx;m1++)
for (m2=1;m2<ynumy;m2++)
	{
	/* Pos in vx[], vy[], pr[], etc. */
	mcmax1=m1*ynumy+m2;
	/**/	
	/* Sxx */	
	ival=sxxcalc(m1,m2,0);
	if(ival<xelvis0) xelvis0=ival; if(ival>xelvis1) xelvis1=ival;
	/**/	
	/* Sxy,Exy */	
	if(m1<xnumx-1 && m2<ynumy-1)
		{
		ival=sxycalc(m1,m2,0);
		if(ival<xelvis0) xelvis0=ival; if(ival>xelvis1) xelvis1=ival;
		}
	}
/* Restore stoksmod */	
stoksmod=-stoksmod;
if(printmod) printf("\n !!!   VISCOELASTIC TIME STEP FOR CYCLE %e YEAR !!!\n",timestepe/3.15576e+7);
if(printmod) printf("\n !!!   XLEVIS:    min=%e    max=%e !!!\n",xelvis0,xelvis1);
}
/* End Defining effective numerical timestep for Stockes Equation */
/**/	
/**/	
/**/	
/* Clear New Solution */
/* P, Vx,Vy */
pos0cur=0;
/* Err koef */
bondsum=0;bondnum=0;
stoksum=0;stoknum=0;
contsum=0;contnum=0;
/**/
/**/
/**/
/* Add Matrix by vX-vY-Stokes, Continuity, Boundary, EPS, SIG, Equations */
mukoef=nubeg*strmin;
/* Node  Cycle */
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	/* Pos P,Vx,Vy in sol0[] */
	mcmax=(m1*ynumy+m2)*3;
/*
printf("P-Vx-Vy Cycle %ld %ld    %ld",m1,m2,mcmax); getchar();
*/
	/**/
	/**/
	/**/
	/* Add Continuity equation for Cells ------------------------------------------------ */
	if(m1 && m2)
		{
/*
printf("Pr %ld %ld    %ld",m1,m2,bondm[mcmax+0]); getchar();
*/
		if(!bondm[mcmax+0]) 
			{
			/* Continuity eq. */
	                conterr(m1,m2,0);
			}
		else
			{
			/* Add P-Boundary */
			pbonderr(mcmax+0,0);
			}
		/* Rescale coefficients */
		for (n1=0;n1<=wn[0];n1++) wi[n1]*=mukoef;
	        gausmat4(1,mcmax+0,0);
		}
	else
		{
		wn[0]=1;
		wi[0]=0;
		wn[1]=mcmax+0;
		wi[1]=mukoef;
	        gausmat4(1,mcmax+0,0);
		}

	/**/
	/**/
	/**/
	/* Add vX-Equations --------------------------------------------------- */
	if(m2<ynumy-1)
		{
/*
printf("Vx %ld %ld    %ld",m1,m2,bondm[mcmax+1]); getchar();
*/
		if(!bondm[mcmax+1] || (timesum>timebond && m1>2 && m2>2 && m1<xnumx-4 && m2<ynumy-3)) 
/*
		if(!bondm[mcmax+1]) 
*/
			{
			/* Add vX-Stokes */
	                xstokserr(m1,m2,0);
			/**/
        		/* Add matrix */
			/* vX */
			gausmat4(1,mcmax+1,0);
/*
printf("Vx %ld %ld    %ld",m1,m2,bondm[mcmax+1]); getchar();
*/
			}
		else
			{
			/* Continuity Equation Vx boundary condition */
			if(bondv[bondm[mcmax+1]][1]>1e+30)
				{
				/* m1 m2 increment definition */
				m3=(bondn[bondm[mcmax+1]][0]-1-(mcmax+1))/3;
				dm1=dm2=0;
				if(m3>=ynumy) dm1=1;
				if(m3==1 || m3>ynumy) dm2=1;
/*
printf("Vx(%ld) P(%ld) boundary conditions for X=%ld Y=%ld Node",mcmax+1,bondn[bondm[mcmax+1]][0]-2,m1,m2); getchar();
*/
				if(m3<0 || m3>ynumy+1 || !(m1+dm1) || !(m2+dm2) || !bondm[bondn[bondm[mcmax+1]][0]-2])
					{
					printf("EXIT PROGRAM Inconsistent Vx(%ld) P(%ld) boundary conditions for X=%ld Y=%ld Node",mcmax+1,bondn[bondm[mcmax+1]][0]-2,m1,m2);
					exit(0);
					}
		                conterr(m1+dm1,m2+dm2,0);
				}
			else
				{
				/* Add vX Simple Boundary */
				xbonderr(mcmax+1,0);
				}
			/**/
			/* Vx boundary condition Add */
			/* Rescale coefficients */
			for (n1=0;n1<=wn[0];n1++) wi[n1]*=mukoef;
			gausmat4(1,mcmax+1,0);
			}
		}
	/* Add Ghost parameters */
	else
		{
		wn[0]=1;
		wi[0]=0;
		wn[1]=mcmax+1;
		wi[1]=mukoef;
	        gausmat4(1,mcmax+1,0);
		}
	/**/
	/**/
	/**/
	/* Add vY-Equations --------------------------------------------------- */
	if(m1<xnumx-1)
		{
/*
printf("Vy %ld %ld    %ld",m1,m2,bondm[mcmax+2]); getchar();
*/
		if(!bondm[mcmax+2]) 
			{
			/* Add vX-Stokes */
        	        ystokserr(m1,m2,0);
			/**/
        		/* Add matrix */
			/* vY */
			gausmat4(1,mcmax+2,0);
/*
printf("Vy2 %ld %ld    %ld",m1,m2,bondm[mcmax+8]); getchar();
*/
			}
		else
			{
			/* Continuity Equation Vy boundary condition */
			if(bondv[bondm[mcmax+2]][1]>1e+30)
				{
				/* m1 m2 increment definition */
				m3=(bondn[bondm[mcmax+2]][0]-1-(mcmax+2))/3;
				dm1=dm2=0;
				if(m3>=ynumy) dm1=1;
				if(m3==1 || m3>ynumy) dm2=1;
				if(m3<0 || m3>ynumy+1 || !(m1+dm1) || !(m2+dm2) || !bondm[bondn[bondm[mcmax+2]][0]-3])
/*
printf("Vy(%ld) P(%ld) boundary conditions for X=%ld Y=%ld Node",mcmax+2,bondn[bondm[mcmax+2]][0]-3,m1,m2); getchar();
*/
					{
					printf("EXIT PROGRAM Inconsistent Vy(%ld) P(%ld) boundary conditions for X=%ld Y=%ld Node",mcmax+2,bondn[bondm[mcmax+2]][0]-3,m1,m2);
					exit(0);
					}
	        	        conterr(m1+dm1,m2+dm2,0);
				}
			/* Simple Vy boundary condition */
			else
				{
				/* Add vY Simple Boundary */
				ybonderr(mcmax+2,0);
				}
			/**/
			/* Vy boundary condition Add */
			/* Rescale coefficients */
			for (n1=0;n1<=wn[0];n1++) wi[n1]*=mukoef;
			gausmat4(1,mcmax+2,0);
			}
		}
	/* Add Ghost parameters */
	else
		{
		wn[0]=1;
		wi[0]=0;
		wn[1]=mcmax+2;
		wi[1]=mukoef;
	        gausmat4(1,mcmax+2,0);
		}
	/**/
	/**/
	/**/
	/* Print cycles counter */
	if(1==0 && printmod) 
		{
		ccc++;
		if (ccc>=printmod) 
			{
			printf("X=%ld Y=%ld\n",m1,m2);
			ccc=0;
			}
		}
	}
/* End  Add Matrix By vX-vY-Stokes, Continuity Equations */
/**/
/**/
/**/
/* Solve Matrix ================================================ */
if (printmod) printf("Number of positions in global matrix = %ld  Number of unknown = %ld \n",pos0cur,nodenum3);
gausmat4(0,nodenum3,0);
/* Solve Matrix ================================================ */
/**/
/**/
/**/
/* Reload P, Vx, Vy Results */
/* Node  Cycle */
for (m1=0;m1<xnumx;m1++)
	{
	/* Set Initial p value at upper boundary */
	pmpa=pinit;
	/**/	
	for (m2=0;m2<ynumy;m2++)
		{
		/* Pos P,Vx,Vy in sol0[] */
		mcmax0=(m1*ynumy+m2)*3;
		/* Pos in vx[], vy[], pr[], etc. */
		mcmax1=m1*ynumy+m2;
		/**/	
		/**/	
		/**/	
		/* Reload/Recalc P */	
		if(m1 && m2) 
			{
			/* Reload P */	
			pr[mcmax1]=x[mcmax0+0];
/*	
printf("\n %ld %e",mcmax1,pr[mcmax1]); getchar();
*/	
			}
		/**/	
		/**/	
		/**/	
		/* Reload Vx */	
		if(m2<ynumy-1)
			{ 
			vx[mcmax1]=x[mcmax0+1];
			}
		/**/	
		/**/	
		/**/	
		/* Reload Vy */	
		if(m1<xnumx-1)
			{
			vy[mcmax1]=x[mcmax0+2];
			}
/*	
printf("\n %ld %ld %ld %e %e",m1,m2,mcmax1,vx[mcmax1],vy[mcmax1]); getchar();
*/	
		}
	}
/* End Reload P, Vx, Vy Results */
/**/	
/**/	
/**/	
/* Recalc EPS, SIG Results */
/* Node  Cycle */
maxds=0;
minnu=1e+50;maxnu=-1e+50;
for (m1=1;m1<xnumx;m1++)
for (m2=1;m2<ynumy;m2++)
	{
	/* Pos in vx[], vy[], pr[], etc. */
	mcmax1=m1*ynumy+m2;
	/**/	
	/* Sxx,Exx */	
	sxx[mcmax1]=sxxcalc(m1,m2,0); exx[mcmax1]=eps[0];
	maxds=MAXV(maxds,ABSV(sxx[mcmax1]-sxx0[mcmax1]));
	/* Min,Max Nu value Calc */
	minnu=MINV(minnu,eps[2]);
	maxnu=MAXV(maxnu,eps[2]);
	/**/	
	/* Sxy,Exy */	
	if(m1<xnumx-1 && m2<ynumy-1)
		{
		sxy[mcmax1]=sxycalc(m1,m2,0); exy[mcmax1]=eps[0]; esp[mcmax1]=eps[1];
		maxds=MAXV(maxds,ABSV(sxy[mcmax1]-sxy0[mcmax1]));
		}
	}
/* End Recalc EPS, SIG Results */
/**/	
/**/	
/**/	
/* Vx,Vy max-min definition */
minvx=1e+30;maxvx=-1e+30;minvy=1e+30;maxvy=-1e+30;
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	/* Pos of Vx in sol0[] */
	m3=m1*ynumy+m2;
	/**/
	/* Min,Max Vx definition */
	if(m2<ynumy-1)
		{
		minvx=MINV(minvx,vx[m3]);
		maxvx=MAXV(maxvx,vx[m3]);
		}
	/* Min,Max Vy definition */
	if(m1<xnumx-1)
		{
		minvy=MINV(minvy,vy[m3]);
		maxvy=MAXV(maxvy,vy[m3]);
		}
	}
/* Max Vx,Vy Diff in Grid Calc */
vxmin=minvx;
vymin=minvy;
vxmax=maxvx;
vymax=maxvy;
maxvx-=minvx;
maxvy-=minvy;
maxvx=MAXV(maxvx,maxvy);
mindx=MINV(xstpx,ystpy);
/**/
/**/
/**/
/* Check Error */
/* Node  Cycle */
minpr=1e+50;maxpr=-1e+50;
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	/* Pos of P,Vx,Vy in sol0[] */
	mcmax0=(m1*ynumy+m2)*3;
	/**/
	/**/
	/**/
	/* Check Continuity equation for Cells =========================== */
	if(m1 && m2) 
		{
		ival=conterr(m1,m2,1);
/*	
printf("\n %ld %ld   %e %e %e     %e ",m1,m2,ival,maxvx,mindx,ival/(maxvx/mindx)); getchar();
*/	
		contsum+=ival*ival;
       	        contnum+=1.0;
		/* Print Results */
		ival/=maxvx/mindx;
		if (printmod && ABSV(ival)>DIVVMIN)
			{
			printf("\n Large Continuity err at X=%ld Y=%ld:   Err=%e",m1,m2,ival);
			}
		}
	/**/
	/**/
	/**/
	/* Check vX-Equations for nodes =========================== */
	if(m2<ynumy-1)
		{
		if(!bondm[mcmax0+1]) 
			{
			/* Add vX-Stokes */
	                ival=xstokserr(m1,m2,1);
        	        stoksum+=ival*ival;
                	stoknum+=1.0;
        	        /* Min,Max Pr value Calc */
                	maxpr=MAXV(maxpr,errbuf[10]);
	                minpr=MINV(minpr,errbuf[11]);
			/* Print Results */
			ival/=maxvx/mindx/mindx;
			if (printmod && ABSV(ival)>STOKSMIN)
				{
				printf("\n Large X stokes err at X=%ld Y=%ld:   Err=%e",m1,m2,ival);
				}
	                }
		else
			{
			if(bondv[bondm[mcmax0+1]][1]<1e+30)
				{
				/* Add vX-Boundary */
				ival=xbonderr(mcmax0+1,1);
				bondsum+=ival*ival;
	                	bondnum+=1.0;
				}
			}
		}
	/**/
	/**/
	/**/
	/* Check vY-Equations for nodes =========================== */
	if(m1<xnumx-1)
		{
		if(!bondm[mcmax0+2]) 
			{
			/* Add vX-vY-Stokes */
	                ival=ystokserr(m1,m2,1);
        	        stoksum+=ival*ival;
                	stoknum+=1.0;
        	        /* Min,Max Pr value Calc */
	                maxpr=MAXV(maxpr,errbuf[10]);
        	        minpr=MINV(minpr,errbuf[11]);
			/* Print Results */
			ival/=maxvx/mindx/mindx;
			if (printmod && ABSV(ival)>STOKSMIN)
				{
				printf("\n Large Y stokes err at X=%ld Y=%ld:   Err=%e",m1,m2,ival);
				}
                	}
		else
			{
			if(bondv[bondm[mcmax0+2]][1]<1e+30)
				{
				/* Add vX-vY-Boundary */
				ival=ybonderr(mcmax0+2,1);
				bondsum+=ival*ival;
                		bondnum+=1.0;
				}
			}
		}
/*	
printf("\n %ld %ld   %e %e %e ",m1,m2,minnu,maxnu,errbuf[9]); getchar();
*/	
	}
stoksum=pow(stoksum/stoknum,0.5)/(maxvx/mindx/mindx);
contsum=pow(contsum/contnum,0.5)/(maxvx/mindx);
bondsum=pow(bondsum/bondnum,0.5)/maxvx;
/* End Check Error */
/**/
/**/
/**/
/* Print Results */
if (printmod)
	{
	 printf("\n KRUG %2d \n MIN/MAX VELOCITY Vx = % e / % e     Vy = %e / % e\n",m0+1,vxmin,vxmax,vymin,vymax);
	 printf("VISKOS: min = %e max = %e \n",minnu,maxnu);
	 printf("PRESSURE: min = %e max = %e \n",minpr,maxpr);
	 printf("STRESSS CHANGE:    max = %e \n",maxds);
	 printf("STOKS: num = %e err = %e \n",stoknum,stoksum);
	 printf("CONT : num = %e err = %e \n",contnum,contsum);
	 printf("BOUND V: num = %e err = %e \n",bondnum,bondsum);
	 }
/**/
/**/
/**/
/* Time step for markers definition */
if(markmod)
	{
	maxvx=MAXV(ABSV(vxmin),ABSV(vxmax));
	maxvy=MAXV(ABSV(vymin),ABSV(vymax));
	if (maxvx)
		{
		maxvx=(maxxystep*xstpx)/maxvx;
		if(printmod) printf("\n !!! MAX VALID TIME STEP FOR Vx-MARKER %e YEAR !!!\n",maxvx/3.15576e+7);
		if(timestep>maxvx) timestep=maxvx;
		}
	if (maxvy)
		{
		maxvy=(maxxystep*ystpy)/maxvy;
		if(printmod) printf("\n !!! MAX VALID TIME STEP FOR Vy-MARKER %e YEAR !!!\n",maxvy/3.15576e+7);
		if(timestep>maxvy) timestep=maxvy;
		}
	if(printmod) printf("\n !!!       CURRENT TIME STEP FOR CYCLE %e YEAR !!!\n",timestep/3.15576e+7);
	}
/**/
/**/
/**/
/* Recalc EPS, SIG Results for new timestep  */
/* Node  Cycle */
if(stoksmod && timestep<timestepe)
{
timestepe=timestep;
maxds=0;
if(printmod) printf("\n !!! RECALC SIGMA FOR NEW TIMESTEP=%e !!!\n",timestep/3.15576e+7);
for (m1=1;m1<xnumx;m1++)
for (m2=1;m2<ynumy;m2++)
	{
	/* Pos in vx[], vy[], pr[], etc. */
	mcmax1=m1*ynumy+m2;
	/**/	
	/* Sxx,Exx */	
	sxx[mcmax1]=sxxcalc(m1,m2,0); exx[mcmax1]=eps[0];
	maxds=MAXV(maxds,ABSV(sxx[mcmax1]-sxx0[mcmax1]));
	/**/	
	/* Sxy,Exy */	
	if(m1<xnumx-1 && m2<ynumy-1)
		{
		sxy[mcmax1]=sxycalc(m1,m2,0); exy[mcmax1]=eps[0]; esp[mcmax1]=eps[1];
		maxds=MAXV(maxds,ABSV(sxy[mcmax1]-sxy0[mcmax1]));
		}
	}
if(printmod) printf("NEW STRESSS CHANGE:    max = %e \n",maxds);
}
/* End Recalc EPS, SIG Results for new timestep  */
/**/
/**/
/**/
/* Recalc Viscoelastic stresses for markers */
m5=nodenum*2;
/* Clear nodes wt, save increment */
for (m1=0;m1<nodenum;m1++)
	{
	sol0[m1]=sol0[nodenum+m1]=sol0[m5+m1]=sol1[m1]=sol1[nodenum+m1]=sol1[m5+m1]=0;
	}
/**/
/* Set values for newly coming marker */
for (m3=0;m3<=marknum;m3++) 
/* Check markers out of grid */
if (markx[m3]>0 && marky[m3]>0 && (double)(markx[m3])<xsize && (double)(marky[m3])<ysize && markt[m3]<50 && markk[m3]<=0)
	{
	/* Interpolate Stresses */
	allinterd((double)(markx[m3]),(double)(marky[m3]));
	/**/
	/* Reset marker stresses for newly coming markers */
	markxx[m3]=(float)(eps[32]);
	markp[m3]=(float)(eps[33]);
	markxy[m3]=(float)(eps[31]);
	markexx[m3]=(float)(eps[35]);
	markexy[m3]=(float)(eps[34]);
	/* Set newly coming marker if no temperature calculation */
	if(!tempmod) markk[m3]=(float)(eps[2]);
/*
	markexx[m3]=(float)(eps[6]);
	markexy[m3]=(float)(eps[4]);
printf("DIFFUSION1 %ld %e",m3,markk[m3]); getchar();
*/
	}
/**/
/* Recalc marker SIGij + Diffusion, Add  nodes wt */
for (m3=0;m3<=marknum;m3++) 
/* Check markers out of grid */
if (markx[m3]>0 && marky[m3]>0 && (double)(markx[m3])<xsize && (double)(marky[m3])<ysize && markt[m3]<50 && markk[m3]>0)
	{
	/* Interpolate Stresses */
	allinterd((double)(markx[m3]),(double)(marky[m3]));
	/**/
	/* Numerical diffusion add to markers -------------*/
	if(stredif)
		{
		/**/
		/* Calc difference in marker stresses */
		dsxx=eps[32]-(double)(markxx[m3]);
		dspp=eps[33]-(double)(markp[m3]);
		dsxy=eps[31]-(double)(markxy[m3]);
		dexx=eps[35]-(double)(markexx[m3]);
		dexy=eps[34]-(double)(markexy[m3]);
		/**/
		/* Interpolation of G and Nu from nodes to marker */
		/* Marker weight calculation using dimension of current Cell */
		celdx=gx[wn[0]+1]-gx[wn[0]];
		celdy=gy[wn[1]+1]-gy[wn[1]];
		swt1=1.0/celdx/celdy;
		/**/
		/* Interpolation weights calc  after Fornberg (1996) */
		nodewt(wn[0],wn[0]+1,wn[1],wn[1]+1,(double)(markx[m3]),(double)(marky[m3]),0,0);
		/**/
		/* Calc, check Viscous relaxation to marker stresses */
		mgg=markgg[markt[m3]];
		mnu=markv[m3];
		/* dS=dS0*exp(-G/Nu*dt) */
		dss1=-stredif*mgg/mnu*timestep;
		if(dss1<-150.0) dss1=-150.0;
		dss1=1.0-exp(dss1);
		dsxx1=dsxx*dss1;
		dspp1=dspp*dss1;
		dsxy1=dsxy*dss1;
		dexx1=dexx*dss1;
		dexy1=dexy*dss1;
		/**/
 	      	/* Diffuse Marker Temperature */
		markxx[m3]+=(float)(dsxx1);
		markp[m3]+=(float)(dspp1);
		markxy[m3]+=(float)(dsxy1);
		markexx[m3]+=(float)(dexx1);
		markexy[m3]+=(float)(dexy1);
/*
		if(dss1>stredif) dss1=stredif;
printf("DIFFUSION1 %ld %e  %ld %e %e %e %e %e",m3,m4,mgg,mnu); getchar();
printf("DIFFUSION2 %e %e %e ",markxx[m3],dsxx,dsxx1); getchar();
printf("DIFFUSION4 %e %e %e ",markxy[m3],dsxy,dsxy1); getchar();
*/
		/**/
		/* Wt for nodes calc, add */
		m1min=wn[2]; m1max=wn[3];
		m2min=wn[4]; m2max=wn[5];
		for (m1=m1min;m1<=m1max;m1++)
		for (m2=m2min;m2<=m2max;m2++)
			{
       	 		/* Cur Node Num, wt */
		        m4=m1*ynumy+m2;
			ival=cn[m1-m1min][5]*cn[m2-m2min][4]*swt1;
			/**/
			/* Add Node wt, T */
			sol0[m4]+=ival;
			sol1[m4]+=dsxx1*ival;
			sol1[nodenum+m4]+=dspp1*ival;
			sol1[m5+m4]+=dsxy1*ival;
			sol0[nodenum+m4]+=dexx1*ival;
			sol0[m5+m4]+=dexy1*ival;
			}
		}
	/* End Numerical diffusion add to markers -------------*/
	/**/
	/**/
	/**/
	/* Change marker Stresse after solution */
	markxx[m3]+=(float)(eps[7]-eps[32]);
	markp[m3]+=(float)(eps[10]-eps[33]);
	markxy[m3]+=(float)(eps[5]-eps[31]);
	markexx[m3]+=(float)(eps[6]-eps[35]);
	markexy[m3]+=(float)(eps[4]-eps[34]);
	/**/
	/**/
	/**/
/*
printf("%ld %e %e %e   %e %e %e   %e %e",m3,eps[32],eps[33],eps[31],eps[7],eps[9],eps[5],markxx[m3],markxy[m3]); getchar();
*/
	}
/**/
/* Numerical antidiffusion add to markers -------------*/
if(stredif)
	{
	/* Recalc changes in nodes T */
	for (m1=0;m1<nodenum;m1++)
		{
/*
printf("A %ld %e %e %e",m1,tk0[m1],sol0[m1],sol1[m1]); getchar();
*/
		if(sol0[m1]) 
			{
			/* Averaged Changes in Stresses due to smoothing */
			sol1[m1]/=sol0[m1];
			sol1[nodenum+m1]/=sol0[m1];
			sol1[m5+m1]/=sol0[m1];
			sol0[nodenum+m1]/=sol0[m1];
			sol0[m5+m1]/=sol0[m1];
			}
/*
printf("B %ld %e %e %e %e %e %e ",m1,sxx0[m1],sxy0[m1],sol1[m1],sol1[m5+m1],sol0[m1]); getchar();
*/
		}
	/**/
	/* Recalc marker SIGij + Antidiffusion) */
	for (m3=0;m3<=marknum;m3++) 
	/* Check markers out of grid */
	if (markx[m3]>0 && marky[m3]>0 && (double)(markx[m3])<xsize && (double)(marky[m3])<ysize && markt[m3]<50 && markk[m3]>0)
		{
		/* Interpolate Stresses */
		allinterd((double)(markx[m3]),(double)(marky[m3]));
		/**/
/*
printf("C %ld %e",m3,markk[m3]); getchar();
*/
		/* Wt for nodes calc, add */
		m1min=wn[2]; m1max=wn[3];
		m2min=wn[4]; m2max=wn[5];
		for (m1=m1min;m1<=m1max;m1++)
		for (m2=m2min;m2<=m2max;m2++)
			{
       		 	/* Cur Node Num, wt */
		        m4=m1*ynumy+m2;
			ival=cn[m1-m1min][5]*cn[m2-m2min][4];
			/**/
       		 	/* Antidiffuse Marker Stresses */
			if(sol1[m4]) markxx[m3]-=(float)(sol1[m4]*ival);
			if(sol1[nodenum+m4]) markp[m3]-=(float)(sol1[nodenum+m4]*ival);
			if(sol1[m5+m4]) markxy[m3]-=(float)(sol1[m5+m4]*ival);
			if(sol0[nodenum+m4]) markexx[m3]-=(float)(sol0[nodenum+m4]*ival);
			if(sol0[m5+m4]) markexy[m3]-=(float)(sol0[m5+m4]*ival);
/*
printf("D %ld %ld %e %e ",m3,m4,markxx[m3],markxy[m3],sol1[m4],sol1[m5+m4]); getchar();
*/
			/**/
			}
/*
printf("E %ld %e",m3,markk[m3]); getchar();
*/
		}
	}
/* End Numerical antidiffusion add to markers -------------*/
/**/
/**/
/**/
}
/* End Solve XY-Stokes+Continuity equations by vX,vY,P mixed arbitrary order Finite Diff method */





/* Max Vx,Vy in nodes serch time step recalc */
void maxvelstep()
{
double maxvx=0,maxvy=0,ival;
long int m1,m2,m3;
/**/
/**/
/**/
/* Node  Cycle */
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	/* Pos of in Vx,Vy */
	m3=m1*ynumy+m2;
	/**/
	if(m2<ynumy-1)
		{
		ival=ABSV(vx[m3]);
		maxvx=MAXV(maxvx,ival);
		}
	if(m1<xnumx-1)
		{
		ival=ABSV(vy[m3]);
		maxvy=MAXV(maxvy,ival);
		}
	}
if (maxvx)
	{
	maxvx=(maxxystep*xstpx)/maxvx;
	if(printmod) printf("\n !!! MAX VALID TIME STEP FOR Vx-MARKER %e YEAR !!!\n",maxvx/3.15576e+7);
	timestep=MINV(maxvx,timestep);
	}
if (maxvy)
	{
	maxvy=(maxxystep*ystpy)/maxvy;
	if(printmod) printf("\n !!! MAX VALID TIME STEP FOR Vy-MARKER %e YEAR !!!\n",maxvy/3.15576e+7);
	timestep=MINV(maxvy,timestep);
	}
}
/* Max Vx,Vy in nodes serch time step recalc */




/* Weight of FD calculation for after Fornberg (1996) */
void fdweight(int n, int m, double xi)
/* n - maximal index 0-n */
/* m - required derivative order 0-m */
/* xi - derivation point coordinate */
{
/* Counters */
int i,j,k,mn;
double c1,c2,c3,c4,c5,kk;
/**/
/**/
/**/
c1=1.0;
c4=xn[0]-xi;
for(k=0;k<=m;k++)
	{
	for(j=0;j<=n;j++)
		{
		cn[j][k]=0;
		}
	}
/**/
cn[0][0]=1.0;
for(i=1;i<=n;i++)
	{
	mn=i;if(mn>m) mn=m;
	c2=1.0;
	c5=c4;
	c4=xn[i]-xi;
	for(j=0;j<i;j++)
		{
		c3=xn[i]-xn[j];
		c2*=c3;
		for(k=mn;k>0;k--)
			{
			kk=(double)(k);
			cn[i][k]=c1*(kk*cn[i-1][k-1]-c5*cn[i-1][k])/c2;
			}
		cn[i][0]=-c1*c5*cn[i-1][0]/c2;
		for(k=mn;k>0;k--)
			{
			kk=(double)(k);
			cn[j][k]=(c4*cn[j][k]-kk*cn[j][k-1])/c3;
			}
		cn[j][0]=c4*cn[j][0]/c3;
		}
	c1=c2;
	}
/*
for(i=0;i<=n;i++)printf("FD %d %d %e %d %e %e\n",n,m,xi,i,xn[i],cn[i][m]);getchar();
*/
}
/* Weight of FD calculation after Fornberg (1996) */



/* Left side or Value for Sxx  Equation */ 
/* Sxx=2Nu*Exx*X+Sxx0*(1-X), Exx=1/2(dVx/dX-dVy/dY) */
double sxxcalc(long int m1, long int m2, double ynval)
/* m1,m2 - node X,Y number */
/* ynval - Val Sxx Calc Y(0)/N(koefficient) */
{
/* Exx horisontal position */
double xi=(gx[m1-1]+gx[m1])/2.0,leftsxx=0,nueff,sxxeeff,xelvis=1.0,ggeff=0;
long int v[4];
int n1,n;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   [0]      Vy0       [2] */
/*   Nu0                Nu2 */
/*                          */
/*   Vx0    Sxx3,Exx3   Vx2 */
/*                          */
/*   [1]                [3] */
/*   Nu1      Vy1       Nu3 */
/*                          */
v[0]=(m1-1)*ynumy+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/**/
/* Effective viscosity calc */
nueff=nd[v[3]];
/**/
/**/
/**/
/* Effective elastic stress calc */
if(stoksmod)
	{
	/* Effective shear modulus calc */
	ggeff=gd[v[3]];
	/* Effective viscoelastic factor calc */
	xelvis=ggeff*timestepe/(ggeff*timestepe+nueff);
	/* Effective old elastic stress calc, save */
	sxxeeff=sxxe[v[3]];
	/* Return viscoelasticity factor */
	if(stoksmod<0)
		{
		eps[0]=nueff;
		eps[1]=ggeff;
		return xelvis;
		}
/*
		xelvis=ggeff*timestepe/nueff;
		if(xelvis>150.0) xelvis=150.0;
		xelvis=1.0-exp(-xelvis);
printf("XX %ld %ld %e %e %e %e %e %e",m1,m2,nueff,ggeff,xelvis,1.0-xelvis,sxxeeff,timestepe);getchar();
*/
	/* Nu recalc, Check */
	nueff*=xelvis; if(nueff<nubeg) nueff=nubeg;
	}
/* Nu Save */
eps[2]=nueff;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   [0]      Vy0       [2] */
/*   Nu0                Nu2 */
/*                          */
/*   Vx0    Sxx3,Exx3   Vx2 */
/*                          */
/*   [1]                [3] */
/*   Nu1      Vy1       Nu3 */
/*                          */
/**/
/* Return Sxx,Exx val ----------------------------*/
if(ynval==0)
	{
	/* Exx=1/2(dVx/dX-dVy/dY) */
	leftsxx=0.5*((vx[v[2]]-vx[v[0]])/(gx[m1]-gx[m1-1])-(vy[v[1]]-vy[v[0]])/(gy[m2]-gy[m2-1]));
	/**/
	/* Save Exx */
	eps[0]=leftsxx;
	/**/
	/* Calc Sxx=2Nu*Exx*X+Sxx0*(1-X) */
	leftsxx=2.0*nueff*leftsxx;
	/* Effective elastic stress calc */
	if(stoksmod)
		{
		leftsxx+=sxxeeff*(1.0-xelvis);
		}
	/**/
	return leftsxx;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Sxx ----------------*/
/*  0(P) 1(Vx)  2(Vy)  */
/* Sxx=2Nu*Exx*X+Sxx0*(1-X), Exx=1/2(dVx/dX-dVy/dY) */
/* Add Vx with koefficients */
wn[wn[0]+1]=v[0]*3+1;
wi[wn[0]+1]=-ynval*nueff/(gx[m1]-gx[m1-1]);
wn[wn[0]+2]=v[2]*3+1;
wi[wn[0]+2]=+ynval*nueff/(gx[m1]-gx[m1-1]);
/* Add Vy with koefficients */
wn[wn[0]+3]=v[0]*3+2;
wi[wn[0]+3]=+ynval*nueff/(gy[m2]-gy[m2-1]);
wn[wn[0]+4]=v[1]*3+2;
wi[wn[0]+4]=-ynval*nueff/(gy[m2]-gy[m2-1]);
/**/
/* Add elastic stress to the right part */
if(stoksmod)
	{
	wi[0]-=(1.0-xelvis)*ynval*sxxeeff;
	}
/**/
/* Add total Num of lines */
wn[0]+=4;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Exx %e %ld \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Sxx  Equation */ 




/* Left side or Value for Sxy  Equation */ 
/* Sxy=2Nv*Exy, Exy=1/2(dVx/dY+dVy/dX) */
double sxycalc(long int m1, long int m2, double ynval)
/* m1,m2 - node X,Y number */
/* ynval - Val Syy Calc Y(0)/N(koefficient) */
{
/* Exy position */
double xi,leftsxy=0,leftsxy1=0,leftsxy2=0,nueff,sxyeeff,xelvis=1.0,ggeff=0;
long int v[4];
int n1;
/**/
/**/
/**/
/* Staggered Nodes num */
/*  [0]                [2]            */
/*                                    */
/*                     Vx2            */
/*                                    */
/*  [1]     Vy1        [3]       Vy3  */
/*                   Exy3,Nu3         */
/*                                    */
/*                     Vx3            */
/*                                    */
v[0]=(m1-1)*ynumy+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/**/
/* Effective viscosity calc */
nueff=nu[v[3]];
/**/
/**/
/**/
/* Effective elastic stress calc */
if(stoksmod)
	{
	/* Effective shear modulus calc */
	ggeff=gg[v[3]];
	/* Effective viscoelastic factor calc */
	xelvis=ggeff*timestepe/(ggeff*timestepe+nueff);
	/* Effective old elastic stress calc, save */
	sxyeeff=sxye[v[3]];
	/* Return viscoelasticity factor */
	if(stoksmod<0)
		{
		eps[0]=nueff;
		eps[1]=ggeff;
		return xelvis;
		}
/*
		xelvis=ggeff*timestepe/nueff;
		if(xelvis>150.0) xelvis=150.0;
		xelvis=1.0-exp(-xelvis);
printf("XY %ld %ld %e %e %e %e %e %e",m1,m2,nueff,ggeff,xelvis,1.0-xelvis,sxyeeff,timestepe);getchar();
*/
	/* Nu recalc, Check */
	nueff*=xelvis; if(nueff<nubeg) nueff=nubeg;
	}
/* Nu Save */
eps[2]=nueff;
/**/
/**/
/**/
/* Staggered Nodes num */
/*  [0]                [2]            */
/*                                    */
/*                     Vx2            */
/*                                    */
/*  [1]     Vy1        [3]       Vy3  */
/*                   Exy3,Nu3         */
/*                                    */
/*                     Vx3            */
/* Return Sxy,Exy val ----------------------------*/
if(ynval==0)
	{
	/* Exy=1/2(dVx/dY+dVy/dX)=0 */
	leftsxy1=(vx[v[3]]-vx[v[2]])/(gy[m2+1]-gy[m2-1]);
	leftsxy2=(vy[v[3]]-vy[v[1]])/(gx[m1+1]-gx[m1-1]);
	/**/
	/* Save Exy */
	eps[0]=leftsxy=leftsxy1+leftsxy2;
	/* Save Esp (rotation rate) */
	eps[1]=leftsxy1-leftsxy2;
	/**/
	/* Calc Sxy=2Nu*Exy-Sxyp */
	leftsxy=2.0*nueff*leftsxy;
	/* Effective elastic stress calc */
	if(stoksmod)
		{
		leftsxy+=sxyeeff*(1.0-xelvis);
		}
	/**/
	return leftsxy;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Sxy ----------------*/
/*  0(P) 1(Vx)  2(Vy)  */
/* Sxy=2Nu*Exy, Exy=1/2(dVx/dY+dVy/dX) */
/* Add Vx with koefficients */
wn[wn[0]+1]=v[2]*3+1;
wi[wn[0]+1]=-ynval*2.0*nueff/(gy[m2+1]-gy[m2-1]);
wn[wn[0]+2]=v[3]*3+1;
wi[wn[0]+2]=+ynval*2.0*nueff/(gy[m2+1]-gy[m2-1]);
/* Add Vy with koefficients */
wn[wn[0]+3]=v[1]*3+2;
wi[wn[0]+3]=-ynval*2.0*nueff/(gx[m1+1]-gx[m1-1]);
wn[wn[0]+4]=v[3]*3+2;
wi[wn[0]+4]=+ynval*2.0*nueff/(gx[m1+1]-gx[m1-1]);
/**/
/* Add elastic stress to the right part */
if(stoksmod)
	{
	wi[0]-=(1.0-xelvis)*ynval*sxyeeff;
	}
/**/
/* Add total Num of lines */
wn[0]+=4;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Exy %e %ld \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Sxy  Equation */ 




/* Left side or Err for Compressible Continuity Equation  */
/* div(V) = -D(ln(RO))/dt, div(V)=dVx/dX+dVy/dY */
double conterr(long int m1, long int m2, int ynerr)
/* m1,m2 - node X,Y number */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Counter */
long int v[4];
/* Val Buffer */
double leftc=0,rightc=0,dlnrodp;
int n1;
/**/
/* Staggered Nodes num */
/*   [0]       Vy0      [2] */
/*                          */
/*   Vx0        <P3>    Vx2 */
/*            Exx3,Eyy3     */
/*                          */
/*   [1]       Vy1      [3] */
v[0]=(m1-1)*ynumy+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/**/
/* -D(ln(RO))/dt add to the right part */
rightc=0;
/* D(ln(RO))/dp*dP/dt add to the right part */
dlnrodp=drp[v[3]];
if(timestepe>0)
	{
	rightc=dro[v[3]];
	rightc=(dlnrodp*sppe[v[3]]-rightc)/timestepe;
/*
if(timesum>0){printf("%ld %e %e %e",v[3],dro[v[3]],drp[v[3]],sppe[v[3]]);getchar();}
*/
	}
/**/
/**/
/**/
/* Staggered Nodes num */
/*   [0]       Vy0      [2] */
/*                          */
/*   Vx0        <P3>    Vx2 */
/*            Exx3,Eyy3     */
/*                          */
/*   [1]       Vy1      [3] */
/* Return dVx/dX+dVy/dY err ----------------------------*/
if(ynerr==1)
	{
	/* div(V)=dVx/dX+dVy/dY */
	leftc=(vx[v[2]]-vx[v[0]])/(gx[m1]-gx[m1-1])+(vy[v[1]]-vy[v[0]])/(gy[m2]-gy[m2-1]);
	/* D(ln(RO))/dp*dP/dt add to the left part */
	if(timestepe>0 && dlnrodp)
		{
		leftc+=dlnrodp*pr[v[3]]/timestepe;
		}
/*
	if (bondm[v[0]*3+1] && bondm[v[0]*3+2] && bondm[v[1]*3+2] && bondm[v[2]*3+1]) rightc=0;
*/
	/**/
	return leftc-rightc;
	}
/**/
/**/
/**/
/* Add continuity equation */
/* Set Initial Num of lines -------------------------------------- */
wn[0]=0;
/**/
/* Save Right part for Contin ---------------------*/
wi[0]=rightc;
/**/
/**/
/**/
/* Add Coefficients for left parts of div(V) ----------------*/
/*  0(P) 1(Vx)  2(Vy)  */
/* div(V)=dVx/dX+dVy/dY */
/* Add Vx with koefficients */
wn[1]=v[0]*3+1;
wi[1]=-1.0/(gx[m1]-gx[m1-1]);
wn[2]=v[2]*3+1;
wi[2]=+1.0/(gx[m1]-gx[m1-1]);
/* Add Vy with koefficients */
wn[3]=v[0]*3+2;
wi[3]=-1.0/(gy[m2]-gy[m2-1]);
wn[4]=v[1]*3+2;
wi[4]=+1.0/(gy[m2]-gy[m2-1]);
/* Add total Num of lines */
wn[0]=4;
/**/
/* D(ln(RO))/dp*dP/dt add to the left part */
if(timestepe>0 && dlnrodp)
	{
	wn[5]=v[3]*3;
	wi[5]=dlnrodp/timestepe;
	wn[0]=5;
	}
/**/
/**/
/**/
/* Check Boundary conditions around Cell */
leftc=1.0;
if (!bondm[v[0]*3+1]) leftc=0;
if (!bondm[v[0]*3+2]) leftc=0;
if (!bondm[v[1]*3+2]) leftc=0;
if (!bondm[v[2]*3+1]) leftc=0;
/**/ 
/*
for(n1=0;n1<3;n1++)printf("Cont %e %d \n",wi[n1],wn[n1]);getchar();
*/
return leftc;
}
/* Left side or Err for Continuity Equation */

