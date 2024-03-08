/* Move markers by using Simple/Runge-Kutta method */
void movemark()
{
	/* Vx, Vy buffer */
	double dvxdx,dvxdy,dvydx,dvydy,celdx,celdy,vx0,vx1,vx2,vx3,vx4,vy0,vy1,vy2,vy3,vy4,ee0,ee1,ee2,ee3,ee4,sp0,sp1,sp2,sp3,sp4;
	/* Water */
	double vxwater,vywater;
	long int mm1,marknum1,m10,m20,m30,m1,m2,m3;
	/* Erosion-Sedimentation Y/N */
	int n1;
	int mm2;
	/* Nonstabilyty for immobile markers */
	double xnonstab=0.50,ynonstab=0.60;
	double dpdx,dpdy,e,n,vxkoef,vykoef,dx,dy;
	/**/
	/**/
	/**/
	/* Calc Basalt melting depth */
	/*
 	   basalt();
	   printf("Basalt melting depth=%e \n",basalty);
	   */
	/* Dehydration depths  */
	/*
	   dehydration(0);
	   dehydration(1);
	   printf("Serpentine dehydration depths: Min = %e   Max = %e \n",dehydrmin,dehydrmax);
	   */
	/* Hydration front progress  */
	if(vyfluid!=0 && timesum>1e+11) hydration2();
	/*
	*/
	/* Save number of markers */
	marknum1=marknum;
	/**/
	/**/
	/* Surface changes */
	if (timestep) erosion();
	/*
	*/
	/**/
	/**/
	/**/
	/* Move markers */
	for (mm1=0;mm1<marknum;mm1++)
	{
		/* Marker type */
		mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
		if( ((markx[mm1]>=0 && marky[mm1]>=0 && (double)(markx[mm1])<=xsize && (double)(marky[mm1])<=ysize) || outgrid!=1) && !markim[mm2] )
		{
			/**/
			/* Erosion-Sedimentation */
			if((double)(marky[mm1])<=eroslev) n1=1; else n1=0;
			/**/
			/* Interpolate velocity */
			allinteri((double)(markx[mm1]),(double)(marky[mm1]));
			/**/
			/* Water marker move */
			vxwater=vywater=0;
			if(markt[mm1]>=50 && markt[mm1]<100) 
			{
				/* Water velocity */
				vywater=vyfluid; if(markd[mm1]>1100.0) vywater=vymelt;
				/* Fluid in rock */
				if(vyfluid>0 && (markk[mm1]==0 || markk[mm1]>298.0)) 
				{
					/* Horizontal,Vertical P-cell index */
					m1=wn[0]; if(markx[mm1]>(gx[m1]+gx[m1+1])/2.0) m1+=1;
					if(m1<1) m1=1; if(m1>xnumx-2) m1=xnumx-2;
					m2=wn[1]; if(marky[mm1]>(gy[m2]+gy[m2+1])/2.0) m2+=1;
					if(m2<1) m2=1; if(m2>ynumy-2) m2=ynumy-2;
					/* Pressure gradients */
					e=(markx[mm1]-(gx[m1-1]+gx[m1])/2.0)/((gx[m1+1]-gx[m1-1])/2.0);
					n=(marky[mm1]-(gy[m2-1]+gy[m2])/2.0)/((gy[m2+1]-gy[m2-1])/2.0);
					m3=m1*ynumy+m2;
					dpdx=2.0*((1.0-n)*(pr[m3+ynumy]-pr[m3])+n*(pr[m3+ynumy+1]-pr[m3+1]))/(gx[m1+1]-gx[m1-1]);
					dpdy=2.0*((1.0-e)*(pr[m3+1]-pr[m3])+e*(pr[m3+ynumy+1]-pr[m3+ynumy]))/(gy[m2+1]-gy[m2-1]);
					/* Recalc velocity koefficients */
					vxkoef=(1000.0*GXKOEF-dpdx)/(2300.0*9.81);
					vykoef=(1000.0*GYKOEF-dpdy)/(2300.0*9.81);
					/*
					   printf("%ld %ld %e %e   %ld %d %e %e %e %e  %e %e %e %e",m1,m2,gx[m1],gy[m2],mm1,markt[mm1],markx[mm1],marky[mm1],e,n,dpdx,dpdy,vxkoef,vykoef);getchar();
					   */
					if(vxkoef>2.0) vxkoef=2.0; if(vxkoef<-2.0) vxkoef=-2.0;
					if(vykoef>2.0) vykoef=2.0; if(vykoef<-2.0) vykoef=-2.0;
					/* Recalc velocity */
					vxwater=vywater*vxkoef;
					vywater*=vykoef;
				}
				else
					/* Fluid in water */
				{
					vxwater=0;
					vywater=-ABSV(vywater);
				}
				/*
				   printf("%ld %ld %e %e   %ld %d %e %e %e %e  %e %e  %e %e",m1,m2,gx[m1],gy[m2],mm1,markt[mm1],markx[mm1],marky[mm1],dpdx,dpdy,vxkoef,vykoef,vxwater,vywater);getchar();
				   */
				/**/
			}
			/**/
			/* Motion Calc ///////////////////////////////// */
			/**/
			/* Vx, Vy, EpsII Simple calc */
			if(markmod==1)
			{
				vx0=eps[11]+vxwater; vy0=eps[12]+vywater; sp0=eps[30]; ee0=pow(eps[6]*eps[6]+eps[4]*eps[4],0.5);
				/**/
				/*
				   printf("SIMPLE %ld %d %e %e   %e %e %e",mm1,markt[mm1],markx[mm1],marky[mm1],vx0,vy0,sp0); getchar();
				   */
			}
			/* Vx, Vy, EpsII 4 Runge-Kutta koef calc */
			else
			{
				vx1=eps[11]+vxwater; vy1=eps[12]+vywater; sp1=eps[30]; ee1=pow(eps[6]*eps[6]+eps[4]*eps[4],0.5);
				/**/
				/*
				   printf("RK4   %ld %d %e %e   %e %e %e",mm1,markt[mm1],markx[mm1],marky[mm1],vx0,vy0,sp0); getchar();
				   */
				/**/
				allinteri((double)(markx[mm1])+vx1*timestep/2.0,(double)(marky[mm1])+vy1*timestep/2.0);
				vx2=eps[11]+vxwater; vy2=eps[12]+vywater; sp2=eps[30]; ee2=pow(eps[6]*eps[6]+eps[4]*eps[4],0.5);
				/**/
				allinteri((double)(markx[mm1])+vx2*timestep/2.0,(double)(marky[mm1])+vy2*timestep/2.0);
				vx3=eps[11]+vxwater; vy3=eps[12]+vywater; sp3=eps[30]; ee3=pow(eps[6]*eps[6]+eps[4]*eps[4],0.5);
				/**/
				allinteri((double)(markx[mm1])+vx3*timestep,(double)(marky[mm1])+vy3*timestep);
				vx4=eps[11]+vxwater; vy4=eps[12]+vywater; sp4=eps[30]; ee4=pow(eps[6]*eps[6]+eps[4]*eps[4],0.5);
				/**/
				/* Vx,Vy, EpsXX, EpsYY, EpsXY calc after Runge-Kutta */
				vx0=(vx1+2.0*vx2+2.0*vx3+vx4)/6.0;
				vy0=(vy1+2.0*vy2+2.0*vy3+vy4)/6.0;
				if(markmod==2)
				{
					sp0=(sp1+2.0*sp2+2.0*sp3+sp4)/6.0;
					ee0=(ee1+2.0*ee2+2.0*ee3+ee4)/6.0;
				}
				else
				{
					sp0=sp1;
					ee0=ee1;
				}
			}
			/**/
			/* Orthogonal motion only */
			if (outgrid==2)
			{
				if(markx[mm1]<0 || (double)(markx[mm1])>xsize) vy0=0;		
				if(marky[mm1]<0 || (double)(marky[mm1])>ysize) vx0=0;		
			}
			/**/
			/**/
			/**/
			/* Normal/Immobile markers */
			if(markt[mm1]<100)
			{
				/* Markers coming from the depth */
				if(marky[mm1]>zdeep && vy0<0 && markk[mm1]<tdeep) markk[mm1]=tdeep;
				/*
 * 				   if(marky[mm1]>zdeep && vy0<0 && markk[mm1]<tdeep) markk[mm1]=tdeep;
 * 				   				   */
				/* Normal markers */
				/* X,Y calc after Runge-Kutta */
				markx[mm1]+=(float)(timestep*vx0);
				marky[mm1]+=(float)(timestep*vy0);
				if(marke[mm1]>0)
				{
					marke[mm1]+=(float)(timestep*ee0);
				}
				sp0*=timestep;
				/* Turcotte & Schubert, 1995 rotation formula */
				if(stoksmod==1)
				{
					sp1=markxx[mm1]*cos(sp0)*cos(sp0)-markxx[mm1]*sin(sp0)*sin(sp0)+markxy[mm1]*sin(2.0*sp0);
					sp3=0.5*(-markxx[mm1]-markxx[mm1])*sin(2.0*sp0)+markxy[mm1]*cos(2.0*sp0);
					markxx[mm1]=sp1;
					markxy[mm1]=sp3;
				}
				/* Jaumann corrotation formula */
				if(stoksmod==2)
				{
					sp1=markxx[mm1]+markxy[mm1]*2.0*sp0;
					sp3=markxy[mm1]+0.5*(-markxx[mm1]-markxx[mm1])*2.0*sp0;
					markxx[mm1]=sp1;
					markxy[mm1]=sp3;
				}
				/**/
				/* Out of grid marker reset */
				if(markx[mm1]<0 || marky[mm1]<0 || (double)(markx[mm1])>xsize || (double)(marky[mm1])>ysize) 
				{
					markk[mm1]=0;
					markd[mm1]=-1.0;
			markw[mm1]=-1.0;
			marke[mm1]=0;
			}
		}
	else
		{
		/* Immobile markers */
		/* X,Y calc after Runge-Kutta */
		markx[mm1]+=(float)(timestep*vx0);
		marky[mm1]+=(float)(timestep*vy0);
		/**/
/*
 printf("%ld %d %e %e",mm1,markt[mm1],xnew,ynew); getchar();
 * */
		/* Check new position, add marker */
		if(markx[mm1]>=0 && marky[mm1]>=0 && markx[mm1]<=xsize && marky[mm1]<=ysize)
			{
			/* Type save */
			markt[marknum1]=markt[mm1]-100;
			/* X,Y calc after Runge-Kutta */
			markx[marknum1]=markx[mm1];
			marky[marknum1]=marky[mm1];
			/* Temperature Reset */
			markk[marknum1]=0;
			markd[marknum1]=-1.0;
			markv[marknum1]=0;
			/* Strain Reset */
			marke[marknum1]=0;
			/* Stress Reset */
			markxx[marknum1]=0;
			markxy[marknum1]=0;
			/* Pressure Reset */
			markp[marknum1]=0;
			/* Strain rate Reset */
			markexx[marknum1]=0;
			markexy[marknum1]=0;
			/* Add aditional markers counter */
			marknum1++;
			/* X,Y reset for immobile marker */
			markx[mm1]=markk[mm1];
			marky[mm1]=markv[mm1];
			}
		/* Check,Reset old position */
		dx=markx[mm1]-markk[mm1];
		dy=marky[mm1]-markv[mm1];
		dy=pow(dx*dx+dy*dy,0.5);
/*
 * 		if(dy>ystpy || (marky[mm1]<0 && vy0<0) || (marky[mm1]>ysize && vy0>0) || (markx[mm1]<0 && vx0<0) || (markx[mm1]>xsize && vx0>0))
 * 		*/
		if(dy>ystpy)
			{
			/* X,Y reset for immobile marker */
			markx[mm1]=markk[mm1];
			marky[mm1]=markv[mm1];
			}
		/**/
		}
	/**/
	/**/
	/**/
	/* Motion Calc ///////////////////////////////// */
	}
}
/*
printf("\n A Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1-1);
*/
/**/
/**/
/**/
/* Mark num */
if(marknum1>MAXMRK) {printf("Space out in markx[]"); exit(0);}
/**/
/**/
/**/
/* Reset aditional markers */
mm1=0;
while(marknum1>marknum && mm1<marknum)
	{
	/* Reload marker */
	if((markx[mm1]<0 || marky[mm1]<0 || (double)(markx[mm1])>xsize || (double)(marky[mm1])>ysize) && markt[mm1]<100) 
		{
		/* Decrease aditional markers counter */
		marknum1--;
		/* Type save */
		markt[mm1]=markt[marknum1];
		/* Temperature Reset */
		markk[mm1]=0;
		markd[mm1]=-1.0;
		/* Strain Reset */
		marke[mm1]=0;
		/* Stress Reset */
		markxx[mm1]=0;
		markxy[mm1]=0;
		/* Pressure Reset */
		markp[mm1]=0;
		/* Strain rate Reset */
		markexx[mm1]=0;
		markexy[mm1]=0;
		/* X,Y reload  */
		markx[mm1]=markx[marknum1];
		marky[mm1]=marky[marknum1];
		}
	/* Increase markers counter */
	mm1++;
	}
printf("\n Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1);
/* Set new marker number */
marknum=marknum1;
/**/
/**/
/**/
/* Incr cycle of sedimentation */
sedimnum++;
}
/* End Move markers by using Simple/Runge-Kutta method */





/* ro[],nu[] recalc after marker positions */
void ronurecalc()
{
/* Counters */
long int m1,m2,m3,m1min,m1max,m2min,m2max;
int mm2,yn,mm3,n1,n2,ncount=0;
long int mm1;
double dx,dy,swt,swt1,celdx,celdy;
double mnu,mgg,maa,mdro,msxxe,msxye,mexxe,mexye,mro,mcp,mkt,mht,mbb,mdi0,mdi1,mwa,dmwa;
double anu,agg,aro,aaa,acp,akt,aht,abb,wmark;
double wnu,wgg,wro,waa,wcp,wkt,wht,wbb;
double sigin,epsin;
/* TD Database variables,  dTK,dPB - TK, PB step for tabulation in TD database */
double H0,H1,H2,H3,R0,R1,R2,R3,G0,G1,G2,G3,W0,W1,W2,W3,dTK=20.0,dPB=1000.0,n,e;
/* RO, NU equations var */
double mpb=1.0,mtk=300.0,numax=0,numin=0;
/**/
/**/
/**/
if (printmod) printf("\n Number of nodes = %ld  Number of markers = %ld \n",nodenum,marknum);
/**/
/**/
/**/
/* Layering on sediments */
m1=(long int)(sedimnum/sedimcyc);
m2=((long int)(m1/2))*2;
if(m2==m1) yn=3; else yn=4;
/**/
/**/
/**/
/* ADD MARKERS TO THE v-CELLS ========================== */
/* Clear ro[],nu[] wt */
for (m1=0;m1<nodenum;m1++)
	{
	ro0[m1]=0;
	et0[m1]=0;
	nu0[m1]=0;
	nd0[m1]=0;
	gg0[m1]=0;
	gd0[m1]=0;
	sxxe0[m1]=0;
	sppe0[m1]=0;
	sxye0[m1]=0;
	exxe0[m1]=0;
	exye0[m1]=0;
	dro0[m1]=0;
	drp0[m1]=0;
	cp0[m1]=0;
	kt0[m1]=0;
	ht0[m1]=0;
	tk0[m1]=0;
	sol0[m1]=0;
	sol0[nodenum+m1]=0;
	sol1[m1]=0;
	sol1[nodenum+m1]=0;
	}
/**/
/**/
/**/
/* Erosion-sedimentation, Hydration, Melting  account for all markers */
for (mm1=0;mm1<marknum;mm1++)
{
/* Check markers out of grid */
/**/
if(markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && markk[mm1]>0 && markt[mm1]<50)
	{
	/* Marker type */
	mm2=markt[mm1];
	/* Up Left Node X,Y Num */
	wn[0]=m1serch((double)markx[mm1]);
	/**/
	/* Erosion/sedimentation account */
	if(0==0 && erosmod) erosmark(mm1,yn);
	/**/
	/* Water/Air account */
	if(markt[mm1]<2)
		{
		/* Change marker type */
		if((double)(marky[mm1])>waterlev) markt[mm1]=1; else markt[mm1]=0;
/*
	if(erosmod && (double)(marky[mm1])<eroslev) {markt[mm1]=0; marke[mm1]=0; markd[mm1]=-1.0;}
		if((double)(marky[mm1])>sedilev) {markt[mm1]=yn; marke[mm1]=0; markd[mm1]=-1.0;}
		if((double)(marky[mm1])>sedilev) {markt[mm1]=yn; marke[mm1]=0; markd[mm1]=-1.0; printf("%ld %d %e %e %e",mm1,markt[mm1],markx[mm1],marky[mm1],sedilev); getchar();}
*/
		}
	/**/
	/* P, T parameters calc */
	mtk=(double)(markk[mm1]);
	allinterp((double)(markx[mm1]),(double)(marky[mm1]));
	mpb=eps[10]*1e-5;
	/**/
	/* Create oceanic plate*/
	if(0==0)
	{
	/* Create Lithosphere */
	if(mm2==10  && marky[mm1]<2e+5 && markk[mm1]<1473.0)
		{
		markt[mm1]=mm2=9;
		}
	if(0==0)
	{
	/* Create Gabbroic crust */
	if(mm2==9 || mm2==10 || mm2==29 || mm2==30 && markx[mm1]<5e+5)
		{
		/* Erosion level calc */
		e=(markx[mm1]-gx[wn[0]])/(gx[wn[0]+1]-gx[wn[0]]);
		n=(e*ep[wn[0]+1]+(1.0-e)*ep[wn[0]]);
		if((marky[mm1]-n)<=7500.0 && markk[mm1]>1000.0)
			{
			markt[mm1]=mm2=8;
			markw[mm1]=0;
			markd[mm1]=-1.0;
			}
		}
	/* Create Basaltic crust */
	if(mm2==8 && markx[mm1]<5e+5)
		{
		/* Erosion level calc */
		e=(markx[mm1]-gx[wn[0]])/(gx[wn[0]+1]-gx[wn[0]]);
		n=(e*ep[wn[0]+1]+(1.0-e)*ep[wn[0]]);
		if((marky[mm1]-n)<=2500.0)
			{
			markt[mm1]=mm2=7;
			markw[mm1]=0;
			markd[mm1]=-1.0;
			}
		}
	 /* Create pelagic sediments */
         if(mm2==7 && markx[mm1]<5e+5)
                 {
                 /* Erosion level calc */
                 e=(markx[mm1]-gx[wn[0]])/(gx[wn[0]+1]-gx[wn[0]]);
                 n=(e*ep[wn[0]+1]+(1.0-e)*ep[wn[0]]);
                 if((marky[mm1]-n)<=500.0)
                         {
                         markt[mm1]=mm2=4;
                         markw[mm1]=0;
                         markd[mm1]=-1.0;
                         }
                 }
	}
	}
/*
	mpb=(double)(markp[mm1]*1e-5);
	allinterp((double)(markx[mm1]),(double)(marky[mm1]));
	mpb=eps[10]*1e-5;
	depthp((double)(markx[mm1]),(double)(marky[mm1]));eps[10]=eps[50];
*/
	/**/
	/* Serpentinization of brittle mantle faults at sub-surface */
	if((markt[mm1]==9 || markt[mm1]==9 || markt[mm1]==9) && marke[mm1]>deserp && marky[mm1]<dyserp) 
		{
		/* Mantle to Antigorite transformation */
		markt[mm1]=13; 
		markd[mm1]=-1.0; 
/*
		antigor(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1);
		if(markt[mm1]==11) 
			{
			markt[mm1]=9; 
			}
		else
			{
			markd[mm1]=-1.0; 
			}
*/
		}
	/**/
	/* Mantle to Antigorite transformation */
	antigor(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1);
/*
*/
	/**/
	/* Rocks to rock+melt transformation */
	melting(mtk,mpb,mm1);
/*
*/
	}
}
/**/
/**/
/**/
/* Add ro[] nu[] etc. using selected markers */
for (mm1=0;mm1<marknum;mm1+=gridmod)
{
/* Marker type */
mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
/* Check markers out of grid */
/**/
if(markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && markk[mm1]>0 && markt[mm1]<50)
	{
	/**/
	/* Remove water, rocks */
	if(erosmod==0)
		{
		if(marky[mm1]>sedilev && mm2<2) 
			{
			mm2=yn; markt[mm1]=yn;
			markw[mm1]=0;
			markd[mm1]=-1.0;
			}
		if(marky[mm1]<eroslev && mm2>1) 
			{
			if((double)(marky[mm1])>waterlev) markt[mm1]=1; else markt[mm1]=0;
			mm2=markt[mm1];
			markw[mm1]=0;
			markd[mm1]=-1.0;
			}
		}
/*
*/
	/* Remove Plumes */
	if(marky[mm1]>zdeep && mm2!=10) 
		{
		mm2=10; markt[mm1]=10;
		markw[mm1]=0;
		markd[mm1]=-1.0;
		}
/*
*/
	/* P, T parameters calc */
	allintere((double)(markx[mm1]),(double)(marky[mm1]));
	mpb=eps[10]*1e-5;
/*
	mpb=(double)(markp[mm1]*1e-5);
	mpb=eps[10]*1e-5;
*/
	mtk=(double)(markk[mm1]);
	/* Reset water/air temperature */
/*
	if (mm2<2) mtk=markk[mm1]=273.0;
*/
	/**/
	/* Density parameters calc */
	mro=dencalc(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm2);
	mbb=eps[20];
	maa=eps[19];
	mcp=markcp[mm2];
	mkt=(markkt[mm2]+markkf[mm2]/(mtk+77.0))*exp(markkp[mm2]*mpb);
	/* Test Heat conductivity k=ko/(1+b*(T-To)/To) */
	if (markkt[mm2]<0) mkt=-markkt[mm2]/(1.0+markkf[mm2]*(mtk-markkp[mm2])/markkp[mm2]);
/*
ncount++;if(ncount>100000){ncount=0;printf("k %ld %d %e %e %e %e %e %e \n",mm1,mm2,mtk,markk[mm1],mkt,markkt[mm2],markkf[mm2],markkp[mm2]);}
printf("k %ld %d %e %e %e %e %e %e",mm1,mm2,mtk,markk[mm1],mkt,markkt[mm2],markkf[mm2],markkp[mm2]);getchar();
*/
	mht=markht[mm2];
	/**/
	/* Molten rocks */
	if (mm2>20) 
		{
		if (meltpart(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1,mm2));
			{
			mro=eps[23];
			mbb=eps[20];
			maa=eps[19];
			mnu=eps[24];
			mcp=eps[25];
			mkt=eps[26];
			mgg=eps[27];
			}
		}
/*
*/
	/**/
	/* Thermodynamic database use for water  */
	/* Density && Water wt% save */
	if (densimod==3)
	if(mm2>1)
	if(1==0 || timesum<=1e+11 || markd[mm1]<=0) 
		{
		tdbasecalc(mtk,mpb,mm2,mm1);
		markw[mm1]=eps[42];
		markd[mm1]=mro;
		}
	/**/
	/* Thermodynamic database use for ro, Cp */
	if (densimod==2)
/*
	if(mm2>1 && mm2!=5 && mm2!=6 && mm2!=25 && mm2!=26)
*/
	if(mm2>1)
		{
		/* Compute TD variables */
		tdbasecalc(mtk,mpb,mm2,mm1);
/*
		mgg=markgg[mm2];
*/
		mgg=markgg[mm2]=eps[40];
		mro=eps[41];
		mwa=eps[42];
		mcp=eps[43];
		mbb=eps[44];
		maa=eps[45];
		/* Density && Water wt% save */
		if(1==0 || timesum<=1e+11 || mwa<=0 || markd[mm1]<=0) 
			{
			markw[mm1]=mwa;
			markd[mm1]=mro;
			}
		else
			{
			/* Recompute rock density on the basis of water density */
			wro=1050.0;
			dmwa=markw[mm1]-mwa;
			mro=mro/(1.0+dmwa*1e-2*(mro/wro-1.0));
			markd[mm1]=mro;
			}
/*
{printf("TD1 %d %d %e %e   %e %e   %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1]);getchar();}
{printf("TD2 %d %d %e %e   %e %e   %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1]);getchar();}
*/
		}
	/**/
	/* Marker Rheology */
	mdi0=0;
	mdi1=1.0;
	if(mm2<=20) 
		{
		mnu=viscalc(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1,mm2,0);
		/**/
		/* Density correction for the dilation angle */
		mdi0=eps[18];
		mdi1=1.0;
		if(markf0[mm2]>0 && markf1[mm2]>0 && marke[mm1]>0)
			{
			/* Second invariant of viscoplastic strain calc, check */
			sigin=pow(markxx[mm1]*markxx[mm1]+markxy[mm1]*markxy[mm1],0.5);
			epsin=marke[mm1]-sigin/2.0/markgg[mm2];
			if(epsin>markf1[mm2]) epsin=markf1[mm2];
			if(epsin>0) mdi1=exp(-2.0*epsin*markf0[mm2]);
/*
if(mnu<2e+13) {printf("k %ld %d  %e %e   %e %e   %e %e   %e",mm1,mm2,markx[mm1],marky[mm1],marke[mm1]*2.0*markgg[mm2],sigin,marke[mm1],epsin,mnu);getchar();}
*/
			}
		}
	/**/
	msxxe=markxx[mm1];
	msxye=markxy[mm1];
	mexxe=markexx[mm1];
	mexye=markexy[mm1];
	/* Shear modulus */
	mgg=markgg[mm2];
	/* Min,Max NU limitation */
	if(mnu<nubeg) mnu=nubeg; if(mnu>nuend) mnu=nuend;
	/**/
	/**/
	/* End Errosion/Sedimentation account --------------------- */
	/**/
	/**/
	/**/
	/* Water/Air account ====================================*/
	if(mm2<2)
		{
		markd[mm1]=mro;
		mdi0=0;
		mdi1=1.0;
		}
	/* End Water/Air account ====================================*/
	/**/
	/**/
	/**/
	/* Calc log density derivative, save new density */
	if(markd[mm1]<=0 || densimod==0) {markd[mm1]=mro; mdi1=1.0; mdi0=0; maa=0;}
	mdro=0; 
	maa=0; 
	mdi1=1.0; 
	if(timestepe) 
		{
		mdro=mro/markd[mm1];
		mdro=log(mdro)-mdi0;
/*
if(epsin>0) {printf("d %ld %d  %e %e   %e %e   %e %e   %e %e %e %e",mm1,mm2,markx[mm1],marky[mm1],marke[mm1]*2.0*markgg[mm2],sigin,marke[mm1],epsin,-2.0*epsin*markf0[mm2],markd[mm1],mro,mdro);getchar();}
*/
		}
	/* Save new density */
	mdro=-mdi0; 
	markd[mm1]=mro;
	/* Correct new density for dilation */
	mro*=mdi1;
	/**/
	/**/
	/**/
	/* Saving marker viscosity */
	markv[mm1]=mnu;
	/**/
	/**/
	/**/
/*
printf("num=%ld type=%d  x=%e y=%e mpb=%e mtk=%e nu=%e ro=%e cp=%e kt=%e ht=%e",mm1,mm2,markx[mm1],marky[mm1],mpb,mtk,mnu,mro,mcp,mkt,mht);getchar();
*/
	/* Interpolation from markers to 4 corners of the cell ====================================*/
	/* Marker weight calculation using dimension of current Cell */
	celdx=gx[wn[0]+1]-gx[wn[0]];
	celdy=gy[wn[1]+1]-gy[wn[1]];
	swt1=1.0/celdx/celdy;
	/* Marker weights calculation using dimension of current Cell */
	celdx=(markx[mm1]-gx[wn[0]])/(gx[wn[0]+1]-gx[wn[0]]);
	celdy=(marky[mm1]-gy[wn[1]])/(gy[wn[1]+1]-gy[wn[1]]);
if (celdx<0 || celdy<0 || celdx>1.0 ||celdy>1.0) {printf("num=%ld type=%d  x=%e y=%e celdx=%e celdy=%e",mm1,mm2,markx[mm1],marky[mm1],celdx,celdy);getchar();}
	/* Interpolate ro,Nu etc to nodes using interpolation koefficients */
	for (m1=0;m1<4;m1++)
		{
		/* Marker weight calculation using dimension of current Cell */
		/* Different corners */
		/* 0  2 */
		/* 1  3 */
		switch(m1)
			{
			case  0: 
			/* Calc node number */
			m3=wn[0]*ynumy+wn[1];
			/* Add shear viscosity Nu */
			if (celdx<0.5 && celdy<0.5) 
				{
				dx=1.0-2.0*celdx;
				dy=1.0-2.0*celdy;
				swt=swt1*dx*dy;
				nu0[m3]+=mnu*swt;
				gg0[m3]+=1.0/mgg*swt;
				sxye0[m3]+=msxye*swt;
				exye0[m3]+=mexye*swt;
				sol0[nodenum+m3]+=swt;
				}
			/* Calc standard wt */
			swt=swt1*(1.0-celdx)*(1.0-celdy); 
			break;
			/**/
			case  1: 
			/* Calc node number */
			m3=wn[0]*ynumy+wn[1]+1;
			/* Add shear viscosity Nu */
			if (celdx<0.5 && celdy>0.5) 
				{
				dx=1.0-2.0*celdx;
				dy=2.0*celdy-1.0;
				swt=swt1*dx*dy;
				nu0[m3]+=mnu*swt;
				gg0[m3]+=1.0/mgg*swt;
				sxye0[m3]+=msxye*swt;
				exye0[m3]+=mexye*swt;
				sol0[nodenum+m3]+=swt;
				}
			/* Calc standard wt */
			swt=swt1*(1.0-celdx)*celdy; 
			break;
			/**/
			case  2: 
			/* Calc node number */
			m3=(wn[0]+1)*ynumy+wn[1];
			/* Add shear viscosity Nu, Sxy */
			if (celdx>0.5 && celdy<0.5) 
				{
				dx=2.0*celdx-1.0;
				dy=1.0-2.0*celdy;
				swt=swt1*dx*dy;
				nu0[m3]+=mnu*swt;
				gg0[m3]+=1.0/mgg*swt;
				sxye0[m3]+=msxye*swt;
				exye0[m3]+=mexye*swt;
				sol0[nodenum+m3]+=swt;
				}
			/* Calc standard wt */
			swt=swt1*celdx*(1.0-celdy); 
			break;
			/**/
			case  3: 
			/* Calc node number */
			m3=(wn[0]+1)*ynumy+wn[1]+1;
			/* Add shear viscosity Nu */
			if (celdx>0.5 && celdy>0.5) 
				{
				dx=2.0*celdx-1.0;
				dy=2.0*celdy-1.0;
				swt=swt1*dx*dy;
				nu0[m3]+=mnu*swt;
				gg0[m3]+=1.0/mgg*swt;
				sxye0[m3]+=msxye*swt;
				exye0[m3]+=mexye*swt;
				sol0[nodenum+m3]+=swt;
				}
			/* Add Normal viscosity Nd, Sxx, Syy */
				{
				dx=1.0-2.0*ABSV(celdx-0.5);
				dy=1.0-2.0*ABSV(celdy-0.5);
				swt=swt1*dx*dy;
				nd0[m3]+=mnu*swt;
				gd0[m3]+=1.0/mgg*swt;
				sxxe0[m3]+=msxxe*swt;
				sppe0[m3]+=markp[mm1]*swt;
				exxe0[m3]+=mexxe*swt;
				dro0[m3]+=mdro*swt;
				drp0[m3]+=maa*swt;
				sol1[nodenum+m3]+=swt;
				}
			/**/
			/* Calc standard wt */
			swt=swt1*celdx*celdy; 
			break;
			}
		/**/
		/* Add Physical Properties: ro,nu, etc. */
/*
			if (0==0) 
printf("num=%ld type=%d  x=%e y=%e cell=%ld swt=%e",mm1,mm2,markx[mm1],marky[mm1],m3,swt);getchar();
nu0[m3]+=mnu*swt;
*/
		ro0[m3]+=mro*swt;
		et0[m3]+=mbb*swt;
		cp0[m3]+=mcp*mro*swt;
		kt0[m3]+=mkt*swt;
		ht0[m3]+=mht*swt;
		sol0[m3]+=swt;
		/**/
		/* Add T */
		if(!markim[mm2]) 
			{
			tk0[m3]+=mtk*swt;
			sol1[m3]+=swt;
			}
		}
	/* End Interpolation from markers to nodes ====================================*/
	}
}
/**/
/**/
/* Recalc ro[] nu[] */
for (m3=0;m3<nodenum;m3++)
/* Recalc ro[] nu[] */
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
{
/* Current node num, wt */
m3=m1*ynumy+m2;
/* Shear viscosity recalc check */
if(sol0[nodenum+m3])
	{
	if(mu[m3] && (timesum<timebond || m1<=2 || m2<=2 || m1>=xnumx-4 || m2>=ynumy-3)) 
/*
	if(mu[m3]) 
*/
		{
		if(mu[m3]>0) 
			{
			nu0[m3]=mu[m3];
			}
		else
			{
			nu0[m3]/=sol0[nodenum+m3];
			if(nu0[m3]>-mu[m3]) nu0[m3]=-mu[m3];
			}
		} 
	else
		{
		nu0[m3]/=sol0[nodenum+m3];
		} 
	/* Min,Max NU limitation */
	if(nu0[m3]<nubeg) nu0[m3]=nubeg; if(nu0[m3]>nuend) nu0[m3]=nuend;
	/* Min,Max NU definition for nu contrast limit */
	if(numin==0 || nu0[m3]<numin) numin=nu0[m3]; if(numax==0 || nu0[m3]>numax) numax=nu0[m3];
	nu[m3]=nu0[m3];
	/**/
	/* Elastic shear stress Sxy recalc */
	sxye[m3]=sxye0[m3]/sol0[nodenum+m3];
	exye[m3]=exye0[m3]/sol0[nodenum+m3];
	/**/
	/* Shear shear modulus recalc */
	gg[m3]=1.0/(gg0[m3]/sol0[nodenum+m3]);
	/**/
	/* Reset weight */
	sol0[nodenum+m3]=0;
	} 
/**/
/* Normal viscosity recalc check */
if(sol1[nodenum+m3])
	{
	if(mu[m3] && (timesum<timebond || m1<=2 || m2<=2 || m1>=xnumx-4 || m2>=ynumy-3)) 
/*
	if(mu[m3]) 
*/
		{
		if(mu[m3]>0) 
			{
			nd0[m3]=mu[m3];
			}
		else
			{
			nd0[m3]/=sol1[nodenum+m3];
			if(nd0[m3]>-mu[m3]) nd0[m3]=-mu[m3];
			}
		} 
	else
		{
		nd0[m3]/=sol1[nodenum+m3];
		} 
	/* Min,Max NU limitation */
	if(nd0[m3]<nubeg) nd0[m3]=nubeg; if(nd0[m3]>nuend) nd0[m3]=nuend;
	/* Min,Max NU definition for nu contrast limit */
	if(numin==0 || nd0[m3]<numin) numin=nd0[m3]; if(numax==0 || nd0[m3]>numax) numax=nd0[m3];
	nd[m3]=nd0[m3];
	/**/
	/* Elastic Normal stress recalc */
	sxxe[m3]=sxxe0[m3]/sol1[nodenum+m3];
	sppe[m3]=sppe0[m3]/sol1[nodenum+m3];
	exxe[m3]=exxe0[m3]/sol1[nodenum+m3];
	/* Density changes recalc */
	dro[m3]=dro0[m3]/sol1[nodenum+m3];
	drp[m3]=drp0[m3]/sol1[nodenum+m3];
	/**/
	/* Normal shear modulus recalc */
	gd[m3]=1.0/(gd0[m3]/sol1[nodenum+m3]);
	/**/
	/* Reset weight */
	sol1[nodenum+m3]=0;
	} 
/**/
/* Other variables recalc check */
if(sol0[m3])
	{
	/* Material constants recalc */
	ro[m3]=ro0[m3]/sol0[m3];
	if(gy[m2]<waterlev && ro[m3]<1000.1) ro[m3]=1.0;
	if(gy[m2]>=waterlev && ro[m3]<1000.1) ro[m3]=1000.0;
	et[m3]=et0[m3]/sol0[m3];
	cp[m3]=(cp0[m3]/sol0[m3])/ro[m3];
	kt[m3]=kt0[m3]/sol0[m3];
	ht[m3]=ht0[m3]/sol0[m3];
	/**/
	/* Advective addition for T K in nodes recalc */
	if (sol1[m3]) 
		{
		tk[m3]=tk0[m3]/sol1[m3]; 
/*
	if (sol1[m3] && tk3[m3]) 
		tk[m3]+=tk0[m3]/sol1[m3]-tk3[m3]; 
printf("%ld %e",m3,tk[m3]);getchar();
*/
		sol1[m3]=0;
		}
	/**/
	/* Reset weight */
	sol0[m3]=0;
	}
}
if (printmod) printf("Min, Max viscosity %e %e \n",numin,numax);
/*
printf("%e %e",numin,numax);getchar();
printf("%e %e %e   %e %e ",nubeg,nuend,nucontr,numin,numax);getchar();
*/
/**/
/* Reset advective temperature */
for (m3=0;m3<nodenum;m3++) tk3[m3]=0; 
/**/
/* Set Upper/Lower limits for nu[] after given contrast */
if(nucontr>1.0 && numin>0) numax=numin*nucontr;
if(nucontr<1.0 && numax>0) numin=numax*nucontr;
for (m3=0;m3<nodenum;m3++)
	{
	if(nu[m3]<numin) nu[m3]=numin; if(nu[m3]>numax) nu[m3]=numax;
	if(nd[m3]<numin) nd[m3]=numin; if(nd[m3]>numax) nd[m3]=numax;
	}
/**/
/* Water/air density */
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	m3=m1*ynumy+m2;
	if(gy[m2]<waterlev && ro[m3]<1000.1) ro[m3]=1.0;
	if(gy[m2]>=waterlev && ro[m3]<1000.1) ro[m3]=1000.0;
	}
/**/
/* Set Boundary conditions for T */
if (printmod) printf("\n AVERAGE TEMPERATURE CORRECTION FOR BOUNDARY CONDITIONS ...\n");
tkrecalc();
if (printmod) printf("AVERAGE TEMPERATURE OK!\n");
/* Adiabate computing */
if(1==0 && timesum<3.15576e+7*1e+3) 
	{
	/* Lower boundary TK - Node Cycle */
	for (m1=0;m1<xnumx;m1++)
		{
		/* Cur Line Num in bondm[] */
		m2=(m1+1)*ynumy-1;
		m3=bondm[m2+nodenum3];
		if(m3) 
			{
			bondv[m3][0]=tk[m2-1]*2.0-tk[m2-2];
/*
			bondv[m3][0]=tk[m2]+5.0;
mtk=bondv[m3][0];
printf("%ld %e %e %e %e %e",m1,tk[m2-2],tk[m2-1],tk[m2],bondv[m3][0],mtk);getchar();
*/
			}
		}
	}
/**/
/**/
/**/
/*
printf("%e %e %e   %e %e ",nubeg,nuend,nucontr,numin,numax);getchar();
*/
/* ADD MARKERS TO THE v-CELLS ========================== */
}
/* End ro[],nu[] recalc after marker positions */



/* Calc ro for given P,T */
double dencalc(double mtk, double mpb, double x, double y, int mm2)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm2 - Rock number */
{
/* Val buffer */
double ival;
eps[19]=eps[20]=0;
/**/
/* Ro=ro0*(1-bro*(TK-298.15))*(1+aro*(Pkbar-0.001)) */
ival=markro[mm2]*(1.0-markbb[mm2]*(mtk-298.15))*(1.0+markaa[mm2]*(mpb-1.0)*1e-3);
/* Adiabatic term: al=bro/(1-bro*(Tk-298.15)) */
eps[20]=markbb[mm2]/(1.0-markbb[mm2]*(mtk-298.15));
/* Compressibility: be=aro/(1+aro*(Pkbar-0.0001) */
eps[19]=1.e-8*markaa[mm2]/(1.0+markaa[mm2]*(mpb-1.0)*1e-3);
/* Constant density */
if (densimod==0) return markro[mm2];
return  ival;
}
/* End Calc ro for given P,T */



/* Antigorite weakening of mantle */
void antigor(double mtk, double mpb, double x, double y, long int mm1)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm1 - mark number */
{
/* Val buffer */
double k1,sy1,e,hydry,yfiltr,hydryl,tsubd,vxs,vys;
int mm2=(int)markt[mm1];
long int m1;
/**/
/* Check marker type */
if(mm2!=11 &&  mm2!=13) return;
/**/
/* Up Left Node X,Y Num */
m1=wn[0];
m1=m1serch((double)markx[mm1]);
/* Relativ Normalized coord Calc */
e=(x-gx[m1])/(gx[m1+1]-gx[m1]);
/* Erosion surface */
sy1=(e*ep[m1+1]+(1.0-e)*ep[m1]);
/**/
/* Antigorite weakening of mantle above oceanic crust */
/* Atg stability field after Schmidt and Poli, 1998 */
if((y-sy1)>63000.0)
	{
	k1=1013.17699-0.060387633e-3*(y-sy1)-0.004289442e-6*(y-sy1)*(y-sy1);
	}
else
	{
	k1=751.490422+6.00773668e-3*(y-sy1)-0.034690759e-6*(y-sy1)*(y-sy1);
	}
/* Change marker Type */
/* Serpentinized (13) - to hydrated (11) */
if(k1<=mtk && markt[mm1]==13) markt[mm1]=11;
/* Hydrated(11) -  to serpentinized (13) */
if(k1>mtk && markt[mm1]==11) markt[mm1]=13;
}
/* Antigorite weakening of mantle */



/* Basalt melting in oceanic crust */
double basalt()
{
/* Val buffer */
double x0,x1,y0,y1,k0,k1,t0,t1,sy1;
long int m1;
/**/
/**/
/**/
/* Serch for melting depth along lower hydration surface */
for (m1=1;m1<xnumx;m1++)
	{
	/* TK calc along the lower hydr surface */
	x0=gx[m1-1];
	y0=ep[m1+xnumx*3-1];
	allintert(x0,y0); t0=eps[2];
	x1=gx[m1];
	y1=ep[m1+xnumx*3];
	allintert(x1,y1); t1=eps[2];
	/* Oceanic crust top */
	sy1=ep[m1];
	if((y0-sy1)>0 && (y1-sy1)>0)
		{
		/* Basalt solidus temperature calc  */
		if((y0-sy1)>48000.0)
			{
			k0=940.19-0.023633e-3*(y0-sy1)+0.007774e-6*(y0-sy1)*(y0-sy1);
			}
		else
			{
			k0=941.60+704.4157e+3/(y0-sy1)+1217.45485e+6/(y0-sy1)/(y0-sy1);
			}
		if((y1-sy1)>48000.0)
			{
			k1=940.19-0.023633e-3*(y1-sy1)+0.007774e-6*(y1-sy1)*(y1-sy1);
			}
		else
			{
			k1=941.60+704.4157e+3/(y1-sy1)+1217.45485e+6/(y1-sy1)/(y1-sy1);
			}
		/* Basalt solidus depth calc */
		if(t1>=k1)
			{
			basalty=y1-(y1-y0)/((t1-k1)+(k0-t0))*(t1-k1);
			return 0;
			}

		}
	}
basalty=ysize;
return 0;
}
/* Basalt melting in oceanic crust */




/* Nu calc after reological equation */
/* P-T-stress dependent rheology without/with brittle/ductile transition */
/* Reological equations */
/* Stress>SScr */
/* Power law dislocation creep: SSii={NU0*EEii*exp[(E+PV)/RT]}^(1/n) */
/* Effective viscosity: NU=1/2*{NU0*exp[(E+PV)/RT]}^(1/n)*EEii^[(1-n)/n] */
/* Stress<SScr */
/* Newtonian diffusion   creep: SSii=NU1*EEii*exp[(E+PV)/RT] */
/* Effective viscosity: NU=NU0/2*exp[(E+PV)/RT] */
/* NU1=NU0/SScr^(n-1) */
/* SScr - dislocation, diffusion transition stress */
/* SSii - second invariant of deviatoric stress tensor */
/* EEii - second invariant of strain rate tensor */
/* E - activation energy, J */
/* V - activation volume, J/bar */
/* R - gase constant 8.314 J/K */
/* Viscosity NU  calc after reological equations */
/* NU=SSii/(2*EEii) */
/* Brittle - Ductile transition */
/* sbrit=MINV(0.85e+5*pb,60e+6+0.6e+5*pb)*lambda;  (Schott & Scmeling, 1998) */
/* sbrit=MINV(0.667e+5*pb,51.2e+6+0.512e+5*pb)*lambda; (Brace & Kohlsstedt, 1980) */
double viscalc(double mtk, double mpb, double x, double y, long int mm1, int mm2, int yn)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm1 - Marker number */
/* mm2 - rock type */
/* yn - plastic reset yes(0)/no(1) */
{
/* Val buffer */
double xnu,nnu,e,n,rt=8.314*mtk,k1,e1,epsin,sigin,sduct,sbrit,nueff,strain,abrit,bbrit,nubrit,nunewt,nupowl,nuduct;
/* Reological Eq par */
double sy1,mubrit,lamb,xelvis,sxxnew,sxynew,siginnew,mnu0,mnu1,mnu2,siginnew0,siginnew1,siginnew2,dsiginnew0,dsiginnew1,dsiginnew2;
/* Counters */
long int m1;
int ncount=0;
/**/
/* Melted rocks */
if (mm2>20) return markn0[mm2];
/**/
/* Calc effective strain rate, stress after second strain rate Tenzor invariant EEii=(1/2SUM(EPSik^2))^(1/2) */
/*
allintere(x,y);
epsin=pow(eps[6]*eps[6]+eps[4]*eps[4],0.5);
*/
epsin=pow(markexx[mm1]*markexx[mm1]+markexy[mm1]*markexy[mm1],0.5);
sigin=pow(markxx[mm1]*markxx[mm1]+markxy[mm1]*markxy[mm1],0.5);
/**/
/* Brittle "Mohr-Coulomb" viscosity calc: sbrit=A(eps)+LAM*B(eps)*P ------------------- */
/* Lambda brittle weackening factor calc for hydrostatic pore pressure */
/* Up Left Node X,Y Num */
m1=wn[0];
/* Relativ Normalized coord Calc */
e=(x-gx[m1])/(gx[m1+1]-gx[m1]);
n=(e*ep[m1+1]+(1.0-e)*ep[m1]);
lamb=markll[mm2]; 
if ((y-n)<=0) lamb=hidrl;
if ((y-n)>0 && (y-n)<hidry) lamb=hidrl*(1.0-(y-n)/hidry)+lamb*(y-n)/hidry;
/**/
/* Lower friction in fluid/melt present areas */
if (markv[mm1]<0) lamb=lambfld;
/**/
/* Calc Second strain Tenzor invariant EEii=(1/2SUM(EPSik^2))^(1/2) */
strain=marke[mm1];
/*
if(strain>1.5) {printf("%ld %d  %e %e  %e %e  %e ",mm1,mm2,x,y,mtk,mpb,strain); getchar();}
*/
/**/
/* A,B coefficients calc depending on integral strain */
abrit=marka0[mm2]; 
bbrit=markb0[mm2]; 
if(strain>marke1[mm2]) 
	{
	abrit=marka1[mm2]; 
	bbrit=markb1[mm2]; 
	}
else
	{
	if(strain>marke0[mm2] && marke1[mm2]>marke0[mm2]) 
		{
		abrit=marka0[mm2]+(marka1[mm2]-marka0[mm2])*(strain-marke0[mm2])/(marke1[mm2]-marke0[mm2]);
		bbrit=markb0[mm2]+(markb1[mm2]-markb0[mm2])*(strain-marke0[mm2])/(marke1[mm2]-marke0[mm2]);
		}
	}
/**/
/**/
/*
if (mm2>4 && mm2<10 && epsin>0) {printf("%ld %d  %e %e  %e %e  %e %e %e %e %e %e   %e %e ",mm1,mm2,x,y,mtk,mpb,strain,lamb,abrit,bbrit,sbrit,nubrit,epsin,sigin); getchar();}
*/
/**/
/* End Brittle "Mohr-Coulomb" viscosity calc: sbrit=A(eps)+LAM*B(eps)*P ------------------- */
/**/
/**/
/**/
/* Ductile viscosity calc -------------------------------------------*/
/* Inverted value of newtonian NU set */
nunewt=0;
/**/
/* Inverted value of power-low NU set */
nupowl=0;
/**/
/* Check for the presence of ductile rheology */
if (marknu[mm2])
	{
	/* A)  Simple Newtonian rheology */
	/* Newtonian creep: SSii=NU0*2.0*EEii */
	/* Effective viscosity: NU=NU0 */
	/* Effective viscosity member in Stoks: NUs=NU */
	if(markdh[mm2]==0 && markdv[mm2]==0 && (markss[mm2]==0 || markmm[mm2]==1.0))
		{
		/* Inverted value of newtonian NU calc */
		nunewt=1.0/marknu[mm2];
		}
	/**/
	/**/
	/**/
	/* B)  P-T dependent, stress independent Newtonian rheology */
	/* Newtonian diffusion creep: SSii=NU0*EEii*exp[(E+PV)/RT] */
	/* Effective viscosity: NU=NU0*exp[(E+PV)/RT] */
	if((markdh[mm2]!=0 || markdv[mm2]!=0) && (markss[mm2]==0 || markmm[mm2]==1.0))
		{
		/* Inverted value of newtonian NU calc */
		e1=(markdh[mm2]+markdv[mm2]*mpb)/rt;
		if(e1>150.0) e1=150.0;
		/* Test creep Moresi & Solomatov (1995): SSii=NU0*exp[-a*(T-T0)] */
		if(markdh[mm2]<0 && markdv[mm2]>0) 
			{
			e1=markdh[mm2]*(mtk-markdv[mm2]);
			if(e1<-150.0) e1=-150.0;
			}
		/* Test creep Turkotte & Schubert(1982): SSii=NU0*exp[E/RTo(1-(T-T0)/T0)] */
		if(markdh[mm2]<0 && markdv[mm2]<0) 
			{
			e1=(-markdh[mm2])*(1.0-(mtk-(-markdv[mm2]))/(-markdv[mm2]))/8.314/(-markdv[mm2]);
			if(e1>150.0) e1=150.0;
			}
		nunewt=1.0/(marknu[mm2]*exp(e1));
		}
	/**/
	/**/
	/**/
	/* C)  P-T independent, stress dependent rheology without/with brittle/ductile transition */
	/* Stress>SScr */
	/* Power law creep: SSii={NU0*EEii}^(1/n) */
	/* Effective viscosity: NU=1/2*NU0^(1/n)*EEii^[(1-n)/n] */
	/* Effective viscosity member in Stoks: NUs=NU/n */
	/* Stress<SScr */
	/* Newtonian creep: SSii=NU1*EEii */
	/* Effective viscosity: NU=NU1/2 */
	/* Effective viscosity member in Stoks: NUs=NU */
	/* NU1=NU0/SScr^(n-1) */
	if((markdh[mm2]==0 && markdv[mm2]==0) && markss[mm2]!=0 && markmm[mm2]!=1.0)
		{
		/* Koef for stress independent creep NU1 calc */
		k1=marknu[mm2]/pow(markss[mm2],markmm[mm2]-1.0);
		/**/
		/* Inverted value of newtonian NU calc */
		nunewt=1.0/(0.5*k1);
		/**/
		/* Inverted value of power-low NU calc */
		if (sigin>0) nupowl=1.0/(0.5*sigin*marknu[mm2]/pow(sigin,markmm[mm2]));
		}
	/**/
	/**/
	/**/
	/* D)  P-T-stress dependent rheology without/with brittle/ductile transition */
	/* Reological equations */
	/* Stress>SScr */
	/* Power law dislocation creep: SSii={NU0*EEii*exp[(E+PV)/RT]}^(1/n) */
	/* Effective viscosity: NU=1/2*{NU0*exp[(E+PV)/RT]}^(1/n)*EEii^[(1-n)/n] */
	/* Effective viscosity member in Stoks: NUs=NU/n */
	/* Stress<SScr */
	/* Newtonian diffusion   creep: SSii=NU1*EEii*exp[(E+PV)/RT] */
	/* Effective viscosity: NU=NU0/2*exp[(E+PV)/RT] */
	/* Effective viscosity member in Stoks: NUs=NU */
	/* NU1=NU0/SScr^(n-1) */
	if(marknu[mm2]>0 && (markdh[mm2]!=0 || markdv[mm2]!=0) && markss[mm2]!=0 && markmm[mm2]!=1.0)
		{
		/* T-P exponent for effective NU calc */
		e1=(markdh[mm2]+markdv[mm2]*mpb)/rt;
		if(e1>150.0) e1=150.0;
		e1=exp(e1);
		/**/
		/* Koef for stress independent creep NU1 calc */
		k1=marknu[mm2]/pow(markss[mm2],markmm[mm2]-1.0);
		/**/
		/* Inverted value of newtonian NU calc */
		nunewt=1.0/(0.5*k1*e1);
		mnu2=nunewt;
		/**/
		/* Effective viscosity1 calc */
		siginnew1=siginnew=sigin;
		nupowl=0;
		if (siginnew>0) nupowl=1.0/(0.5*siginnew*marknu[mm2]*e1/pow(siginnew,markmm[mm2]));
		mnu1=nupowl;
		mnu0=1.0/(mnu1+mnu2);
		xelvis=markgg[mm2]*timestepe/(markgg[mm2]*timestepe+mnu0);
		siginnew2=2.0*mnu0*epsin*xelvis+sigin*(1.0-xelvis);
		dsiginnew1=siginnew2-siginnew1;
		/* Effective viscosity2 calc */
		siginnew=siginnew2;
		nupowl=0;
		if (siginnew>0) nupowl=1.0/(0.5*siginnew*marknu[mm2]*e1/pow(siginnew,markmm[mm2]));
		mnu1=nupowl;
		mnu0=1.0/(mnu1+mnu2);
		xelvis=markgg[mm2]*timestepe/(markgg[mm2]*timestepe+mnu0);
		siginnew=2.0*mnu0*epsin*xelvis+sigin*(1.0-xelvis);
		dsiginnew2=siginnew-siginnew2;
/*
if(timesum>0 && mnu0<1e+19){printf("TD %ld %d   %e %e %e   %e",mm1,mm2,mnu1,mnu2,mnu0,siginnew2);getchar();}
if(timesum>0 && mnu0<1e+19){printf("TD %ld %d   %e %e %e   %e",mm1,mm2,mnu1,mnu2,mnu0,siginnew);getchar();}
if(timesum>0 && mnu0<1e+19){printf("TD %ld %d   %e %e    %e %e   %e %e   %e %e   %e %e %e",mm1,mm2,mtk-273.15,mpb/1000.0,epsin,sigin,siginnew1,dsiginnew1,siginnew2,dsiginnew2, mnu1,mnu2,siginnew);getchar();}
*/
		/* Dislocation viscosity calc by Bisection method */
		ncount=0;
		do
			{
			dsiginnew0=ABSV(dsiginnew1)+ABSV(dsiginnew2);
			if(dsiginnew0>0)
				{
/*
				dsiginnew0=ABSV(dsiginnew1)/(ABSV(dsiginnew1)+ABSV(dsiginnew2));
*/
				dsiginnew0=0.5;
				siginnew0=siginnew=siginnew1*(1.0-dsiginnew0)+siginnew2*dsiginnew0;
/*
if(timesum>0){printf("TD %ld %d %d   %e %e    %e %e   %e %e   %e %e   %e %e  %e %e %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,epsin,sigin,siginnew1,dsiginnew1,siginnew2,dsiginnew2,siginnew0,dsiginnew0, mnu1,mnu2,mnu0);getchar();}
*/
				nupowl=0;
				if (siginnew>0) nupowl=1.0/(0.5*siginnew*marknu[mm2]*e1/pow(siginnew,markmm[mm2]));
				mnu1=nupowl;
				mnu0=1.0/(mnu1+mnu2);
				xelvis=markgg[mm2]*timestepe/(markgg[mm2]*timestepe+mnu0);
				siginnew=2.0*mnu0*epsin*xelvis+sigin*(1.0-xelvis);
				dsiginnew0=siginnew-siginnew0;
/*
if(timesum>0 && mm1>=2022){printf("NU %ld %d  %e %e    %e %e   %e %e   %e %e   %e %e  %e %e %e",mm1,mm2,mtk-273.15,mpb/1000.0,epsin,sigin,siginnew1,dsiginnew1,siginnew2,dsiginnew2,siginnew0,dsiginnew0, mnu1,mnu2,mnu0);getchar();}
*/
				if((dsiginnew0>=0 && dsiginnew1>=0) || (dsiginnew0<0 && dsiginnew1<0)) 
					{
					siginnew1=siginnew0;
					dsiginnew1=dsiginnew0;
					}
				else
					{
					siginnew2=siginnew0;
					dsiginnew2=dsiginnew0;
					}
				}
			ncount++;
			}
		while(ABSV(dsiginnew0)>10.0 && ncount<101);
/*
if(ncount==101){printf("NU %ld %d  %e %e    %e",mm1,mm2,mtk-273.15,mpb/1000.0,dsiginnew0);getchar();}
*/
		}
	/**/
	/**/
	/**/
	/* F)  P-T-stress dependent rheology of Mantle Karato & Wu (1993) */
	if(marknu[mm2]<0)
		{
		/* Diffusion viscosity calc for given grainsize, mm, Karato & Wu (1993) */
		mnu2=1.0; 
		mnu2=1.0/(0.5*pow(mnu2*1e-3/5e-10,2.5)*8e+10/8.7e+15*exp((3.06e+5+0.6*(eps[46]/28.0*markdv[mm2]))/8.314/mtk));
		mnu0=1.0/mnu2;
		/* Dislocation viscosity calc, Karato & Wu (1993)  Bisection method */
		/* Effective viscosity1 calc */
		siginnew1=siginnew=sigin;
		mnu1=0; 
/* No dislocation creep */
if(0==0)
{
		if(siginnew>0 && 0==0) mnu1=1.0/(0.5*siginnew/pow(siginnew/8e+10,3.5)/3.5e+22*exp((5.40e+5+2.0*mpb)/8.314/mtk));
		if(siginnew>0 && 1==0) mnu1=1.0/(0.5*3.98e+16*exp((5.32e+5+1.6*mpb)/8.314/mtk)/pow(siginnew,2.5));
		mnu0=1.0/(mnu1+mnu2);
		xelvis=markgg[mm2]*timestepe/(markgg[mm2]*timestepe+mnu0);
		siginnew2=2.0*mnu0*epsin*xelvis+sigin*(1.0-xelvis);
		dsiginnew1=siginnew2-siginnew1;
		/* Effective viscosity2 calc */
		siginnew=siginnew2;
		mnu1=0; 
		if(siginnew>0) mnu1=1.0/(0.5*siginnew/pow(siginnew/8e+10,3.5)/3.5e+22*exp((5.40e+5+2.0*mpb)/8.314/mtk));
		mnu0=1.0/(mnu1+mnu2);
		xelvis=markgg[mm2]*timestepe/(markgg[mm2]*timestepe+mnu0);
		siginnew=2.0*mnu0*epsin*xelvis+sigin*(1.0-xelvis);
		dsiginnew2=siginnew-siginnew2;
/*
if(timesum>0){printf("TD %ld %d %d   %e %e    %e %e   %e %e   %e %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,epsin,sigin,siginnew1,dsiginnew1,siginnew2,dsiginnew2);getchar();}
*/
		/* Dislocation viscosity calc by Bisection method */
		ncount=0;
		do
			{
			dsiginnew0=ABSV(dsiginnew1)+ABSV(dsiginnew2);
			if(dsiginnew0>0)
				{
/*
				dsiginnew0=ABSV(dsiginnew1)/(ABSV(dsiginnew1)+ABSV(dsiginnew2));
*/
				dsiginnew0=0.5;
				siginnew0=siginnew=siginnew1*(1.0-dsiginnew0)+siginnew2*dsiginnew0;
/*
if(timesum>0){printf("TD %ld %d %d   %e %e    %e %e   %e %e   %e %e   %e %e  %e %e %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,epsin,sigin,siginnew1,dsiginnew1,siginnew2,dsiginnew2,siginnew0,dsiginnew0, mnu1,mnu2,mnu0);getchar();}
*/
				mnu1=0; 
				if(siginnew>0) mnu1=1.0/(0.5*siginnew/pow(siginnew/8e+10,3.5)/3.5e+22*exp((5.40e+5+2.0*mpb)/8.314/mtk));
				mnu0=1.0/(mnu1+mnu2);
				xelvis=markgg[mm2]*timestepe/(markgg[mm2]*timestepe+mnu0);
				siginnew=2.0*mnu0*epsin*xelvis+sigin*(1.0-xelvis);
				dsiginnew0=siginnew-siginnew0;
/*
{printf("TD %ld %d   %e %e    %e %e   %e %e   %e %e   %e %e  %e %e %e",mm1,mm2,mtk-273.15,mpb/1000.0,epsin,sigin,siginnew1,dsiginnew1,siginnew2,dsiginnew2,siginnew0,dsiginnew0, mnu1,mnu2,mnu0);getchar();}
*/
				if((dsiginnew0>=0 && dsiginnew1>=0) || (dsiginnew0<0 && dsiginnew1<0)) 
					{
					siginnew1=siginnew0;
					dsiginnew1=dsiginnew0;
					}
				else
					{
					siginnew2=siginnew0;
					dsiginnew2=dsiginnew0;
					}
				}
			ncount++;
			}
		while(ABSV(dsiginnew0)>10.0 && ncount<101);
}
		mnu0=1.0/(mnu1+mnu2);
		if(0==0 && timesum<3.15576e+7*1e+5 && mnu0<1e+21) {mnu0=0; mnu1=1e-21;}
		nunewt=mnu2;
		nupowl=mnu1;
/*
if (mpb>0){printf("TD %ld %d   %e %e    %e %e   %e %e   %e %e   %e %e  %e %e %e   %e %e %e %e",mm1,mm2,mtk-273.15,mpb/1000.0,epsin,sigin,siginnew1,dsiginnew1,siginnew2,dsiginnew2,siginnew0,dsiginnew0, mnu1,mnu2,mnu0,eps[46],eps[47],eps[48],eps[46]/28.0*0.6*markdv[mm2]/mpb);getchar();}
*/
		}
	}
/* End Ductile viscosity calc -------------------------------------------*/
/**/
/**/
/**/
/* Ductile effective viscosity calc, check */
nueff=1.0/(nunewt+nupowl);
/* Mantle viscosity */
if(0==0 && (markt[mm1]==9 || markt[mm1]==10) && timesum<3.15576e+7*1e+4 && nueff<1e+20) nueff=1e+20;
/**/
if(nueff<nubeg) nueff=nubeg; if(nueff>nuend) nueff=nuend;
if(nueff<markn0[mm2]) nueff=markn0[mm2]; if(nueff>markn1[mm2]) nueff=markn1[mm2];
nuduct=nueff;
eps[18]=0;
/**/
/**/
/**/
/* Brittle falure stress calc */
mubrit=lamb*bbrit;
if((mubrit || abrit) && epsin)
	{
	/* Add number of markers with plastic rheology */
	plastmark+=1.0;
	/* Checking of yeld criterion */
	sbrit=markp[mm1]*mubrit+abrit;
/*
	sbrit=markp[mm1]*mubrit+abrit*cos(asin(mubrit));
*/
/* Depth-dependent pressure */
if(1==0)
{
depthp(x,y);
sbrit=eps[50]*mubrit+abrit*cos(asin(mubrit));
}
	if(sbrit<0) sbrit=0;
	if(sbrit>marks1[mm2]) sbrit=marks1[mm2];
	/* Viscous case */
	if(sbrit<2.0*nueff*epsin && !stoksmod)
/*
marks1[mm2]=1e+10;
	sbrit=eps[10]*mubrit+abrit*cos(asin(mubrit));
	if(sbrit<2.0*nueff*epsin && !stoksmod)
	sbrit=mpb*1e+5*mubrit+abrit*cos(asin(mubrit));
printf("%e %e %e ",sbrit,abrit,2.0*nueff*epsin);getchar();
printf("%ld %e %e %e %e %e",mm1,x,y,sbrit,epsin,nueff);getchar();
		nueff=0.5*sbrit/epsin;
*/
		{
		/* Recalculating of viscosity */
		nueff=0.5*sbrit/epsin;
		}
	/* Viscoelastic case */
	if(stoksmod && timestepe && epsin)
		{
		/* Future plastic creep */
		/* Future stresses calc */
		xelvis=markgg[mm2]*timestepe/(markgg[mm2]*timestepe+nueff);
		siginnew=2.0*nueff*epsin*xelvis+sigin*(1.0-xelvis);
/*
if(sigin>1e+10) {printf("%ld %d %e %e  %e %e    %e %e   %e %e   %e %e",mm1,mm2,markx[mm1],marky[mm1],markk[mm1]-273.15,markp[mm1]*1e-8,sbrit,sigin,siginnew,2.0*nueff*epsin,nueff,0.5*sbrit/epsin);getchar();}
*/
		if(sbrit<siginnew || sbrit<sigin)
			{
			/* Add number of markers starting/continuing yelding */
			if(sigin<sbrit && marke[mm1]<2e-20) plastnew+=1.0; else plastold+=1.0;
			if(sigin>=sbrit && marke[mm1]<2e-20) plastnew1+=1.0;
			/* Plastic reset */
			if(yn==0)
				{
				/* Density correction for the dilation angle */
				if(markf0[mm2]>0 && markf1[mm2]>0)
					{
					/* Second invariant of viscoplastic strain calc, check */
					e1=marke[mm1]-sbrit/2.0/markgg[mm2];
					/* Correction of divergence rate for plastic strain rate */
					if(e1<markf1[mm2])
						{
						e1=epsin-sbrit/2.0/nuduct;
						if(e1) eps[18]=2.0*e1*markf0[mm2]*timestepe;
						}
					}
				/* Recompute stress */
				if(sigin && sbrit<sigin)
					{
					markxx[mm1]*=sbrit/sigin;
					markxy[mm1]*=sbrit/sigin;
					sigin=sbrit;
					}
				/* New viscosity calc */
				nubrit=sbrit/(2.0*epsin+(sigin-sbrit)/timestepe/markgg[mm2]);
				if(nubrit<nueff) nueff=nubrit;
				/* Set initial plastic strain */
				if(marke[mm1]<=0) marke[mm1]=1e-20;
				}
			}
/*
if(markx[mm1]>1000000.0 && markk[mm1]>1500.0) {printf("%ld %d %e %e  %e %e    %e %e   %e %e   %e %e",mm1,mm2,markx[mm1],marky[mm1],markk[mm1]-273.15,markp[mm1]*1e-8,sbrit,sigin,siginnew,2.0*nueff*epsin,nueff,0.5*sbrit/epsin);getchar();}
			marke[mm1]=1e-50;
*/
		else
			{
			if(yn==0) marke[mm1]=0;
/*
			if (marke[mm1]>0) marke[mm1]=-ABSV(marke[mm1]);
*/
			}
		}
	}
/*
printf("%ld %e %e ",mm1,nueff,abrit);getchar();
*/
/**/
/* Set zero memory stresses for after wall */
/*
if(x>0.260) markxx[mm1]=markxy[mm1]=0;
*/
/**/
/* Check calculated viscosity */
if(nueff<nubeg) nueff=nubeg; if(nueff>nuend) nueff=nuend;
if(nueff<markn0[mm2]) nueff=markn0[mm2]; if(nueff>markn1[mm2]) nueff=markn1[mm2];
/**/
/* Set max viscosity for air/water 4 km above erosion level */
/*
if(mm2<2 && (n-y)>4000.0) nueff=markn1[mm2];
*/
/**/
/**/
/*
if(mm2<2 && x>(0.205+timesum*2.5/100.0/3600.0)) nueff=1000.0;
if(mm2<2 && (x<0.4+timesum*1.5 || x>3.4-timesum*1.5)) nueff=1000.0;
if(mm2<2 && x>(0.450-timesum*2.5/100.0/3600.0)) nueff=1000.0;
if(mm2==12 && x>0.220 && y>0.0445) nueff=10000.0;
if(mm2==10 && x>0.220) nueff=10000.0;
*/
/**/
/**/
/**/
return nueff;
}
/* Nu calc after reological equation */



/* Number of nearest left vertical line find */
long int m1serch(double x)
/* x - X coordinate */
{
/* Variables */
long int m1,m10=0,m11=xnumx-1;
/**/
/* Serch cycle */
do
	{
	m1=(m10+m11)/2;
	if (gx[m1]>x) m11=m1; else m10=m1;
	}
while((m11-m10)>1);
if(m10>xnumx-2) m10=xnumx-2;
/*
if(x<gx[m10] || x>gx[m10+1]) {printf("XXX %ld %ld %ld  %e %e  %e ",m10,m11,m1,gx[m10],gx[m11],x); getchar();}
*/
return m10;
}
/* Number of nearest left vertical line find */




/* Number of nearest upper horizontal line find */
long int m2serch(double y)
/* y - Y coordinate */
{
/* Variables */
long int m2,m20=0,m21=ynumy-1;
/**/
/* Serch cycle */
do
	{
	m2=(m20+m21)/2;
	if (gy[m2]>y) m21=m2; else m20=m2;
	}
while((m21-m20)>1);
if(m20>ynumy-2) m20=ynumy-2;
/*
if(y<gy[m20] || y>gy[m20+1]) {printf("YYY %ld %ld %ld  %e %e  %e ",m20,m21,m2,gy[m20],gy[m21],y); getchar();}
*/
return m20;
}
/* Number of nearest upper horizontal line find */



/* Erosion/Sedimentation Function for markers */
/* mardy - marker vertical size, m */
void erosmark(long int mm1, int yn)
/* mm1 - marker number */
/* yn - current sedimnts type 2,3 */
{
/* Variables */
double e,e0;
long int m1;
/**/
/**/
/**/
/* Surface level elevation definition */
/* Up Left Node X,Y Num */
m1=wn[0];
m1=m1serch((double)markx[mm1]);
/* Relativ Normalized coord Calc */
e=((double)(markx[mm1])-gx[m1])/(gx[m1+1]-gx[m1]);
/**/
/* Surface level elevation for marker definition */
e0=(e*ep[m1+1]+(1.0-e)*ep[m1]);
/**/
/*
printf("MARKER %ld %ld %e %e %e %d %d",mm1,m1,markx[mm1],marky[mm1],e0,markt[mm1],yn);getchar();
*/
/* Marker surface elevation definition */
if(markt[mm1]<2)
	{
	/* Water/Air -> Sediments conversion */
	if((double)(marky[mm1])>e0) {markt[mm1]=yn; marke[mm1]=0; markd[mm1]=-1.0;}
	}
if(markt[mm1]>1)
	{
	/* Rock->Water/Air conversion */
	if((double)(marky[mm1])<e0) {markt[mm1]=0; marke[mm1]=0; markd[mm1]=-1.0;}
	}
/*
if(mm2!=markt[mm1]) {printf("MARKER %ld %ld %e %e %e %d %d %d",mm1,m1,markx[mm1],marky[mm1],e0,mm2,markt[mm1],yn);getchar();}
*/
}
/* Erosion/Sedimentation Function for markers */



/* Rock to rock+melt transformation */
void melting(double mtk, double mpb, long int mm1)
/* mtk - T, K */
/* mpb - P, bar */
/* mm1 - mark number */
{
/* Marker type */
int mm2=(int)markt[mm1];
/**/
/* Melting related cahnge of the marker type */
/* Check marker type */
/*
if (mm2==3 || mm2==4 || mm2==5 || mm2==6 || mm2==7 || mm2==8 || mm2==23 || mm2==24 || mm2==25 || mm2==26 || mm2==27 || mm2==28)
*/
if (mm2==3 || mm2==4 || mm2==5 || mm2==6 || mm2==7 || mm2==8 || mm2==11 || mm2==16 || mm2==23 || mm2==24 || mm2==25 || mm2==26 || mm2==27 || mm2==28 || mm2==34 || mm2==36 || mm2==37 || mm2==38)
if (mpb<0) mpb=0;
switch(mm2)
	{
	/* Sediments, upper crust */
	case  3:
	case  4:
	case  5:
	case  17:
	case  23:
	case  24:
	case  25:
	case  26:
	case  37:
	/* Basalt, Gabbro */
	case  7:
	case  8:
	case  16:
	case  6:
	case  18:
	case  27:
	case  28:
	case  36:
	case  38:
	meltpart1(mtk,mpb,mm2);
	if(eps[21]>0 && mm2<20) {markt[mm1]+=20; marke[mm1]=0;}
	if(eps[21]<=0 && mm2>20) {markt[mm1]-=20; marke[mm1]=0;}
/*
	if(eps[21]>0 && mm2<20) {markt[mm1]+=20; markd[mm1]=-1.0;}
	if(eps[21]<=0 && mm2>20) {markt[mm1]-=20; markd[mm1]=-1.0;}
if (mm2!=markt[mm1]){printf("Granite/Basalt %ld %d %d  %e %e  %e",mm1,mm2,markt[mm1],mtk,mpb,eps[21]);getchar();}
	if(eps[21]>0 && mm2==11) {markt[mm1]=31; markd[mm1]=-1.0;}
	if(eps[21]<=0 && mm2==31) {markt[mm1]=12; markd[mm1]=-1.0;}
*/
 	return;
	/**/
	/* Hydrated Peridotite */
	case  11:
	case  34:
	meltpart1(mtk,mpb,mm2);
	if(eps[21]>0 && mm2==11) {markt[mm1]=34; marke[mm1]=0;}
	if(eps[21]<=0 && mm2==34) {markt[mm1]=14; marke[mm1]=0;}
/*
if (mm2!=markt[mm1]){printf("Peridotite %ld %d %d  %e %e  %e",mm1,mm2,markt[mm1],mtk,mpb,eps[21]);getchar();}
*/
 	return;
	/* Others */
	default: return;
	}
}
/* Rock to rock+melt transformation */




/* Melt fraction, density, viscosity, heat capacity calculation */
double meltpart(double mtk, double mpb, double x, double y, long int mm1, int mm2)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm1 - mark number */
/* mm2 - mark type */
{
/* Val buffer */
double xmelt=0,ival,dmpb,dmtk,epsin,sduct,nueff,smin,smax,nmin,nmax,cpadd=0;
long int m1;
/**/
/* Check marker type */
/*
if (mm2==23 || mm2==24 || mm2==25 || mm2==26 || mm2==27 || mm2==28)
*/
if (mm2==23 || mm2==24 || mm2==25 || mm2==26 || mm2==27 || mm2==28 || mm2==34 || mm2==36 || mm2==37 || mm2==38)
	{
	/* Calculate melt fraction */
	meltpart1(mtk,mpb,mm2);
	xmelt=eps[21];
	/**/
	/* Standard adiabatic term: al=bro/(1+bro*(Tk-298.15)) */
	eps[20]=(markbb[mm2]*xmelt+markbb[mm2-20]*(1.0-xmelt))/(1.0-(markbb[mm2]*xmelt+markbb[mm2-20]*(1.0-xmelt))*(mtk-298.15));
	eps[19]=(markaa[mm2]*xmelt+markaa[mm2-20]*(1.0-xmelt))/(1.0+(markaa[mm2]*xmelt+markaa[mm2-20]*(1.0-xmelt))*(mpb-1.0)*1e-3);
	/**/
	/* Density */
	/* Ro=ro0 */
	if (densimod==0) 
		{
		eps[23]=markro[mm2]*xmelt+markro[mm2-20]*(1.0-xmelt);
		}
	/* Ro=ro0*(1-bro*(TK-298.15))*(1+aro*(Pkbar-0.001)) */
	else
		{
		eps[23]=(markro[mm2]*xmelt+markro[mm2-20]*(1.0-xmelt))*(1.0-(markbb[mm2]*xmelt+markbb[mm2-20]*(1.0-xmelt))*(mtk-298.15))*(1.0+(markaa[mm2]*xmelt+markaa[mm2-20]*(1.0-xmelt))*(mpb-1.0)*1e-3);
		}
	/**/
	/* Viscosity */
	/* Effective NU calc check */
	/* Little melt */
	if(xmelt<0.1)
		{
		eps[24]=viscalc(mtk,mpb,x,y,mm1,mm2-20,0);
		/* Shear modulus */
		eps[27]=markgg[mm2-20];
		}
	else
		{
		/* Significant melt */
		/* Set viscosity and stress limits */
		nmin=MAXV(markn0[mm2],nubeg);
		nmax=MINV(markn1[mm2],nuend);
		smin=MAXV(marks0[mm2],strmin);
		smax=MINV(marks1[mm2],strmax);
		/* Calc effective strain rate after second strain rate Tenzor invariant EEii=(1/2SUM(EPSik^2))^(1/2) */
		epsin=pow(eps[6]*eps[6]+eps[4]*eps[4],0.5);
		/* Effective NU calc check */
		nueff=marknu[mm2]*exp(2.5+pow((1.0-xmelt)/xmelt,0.48)*(1.0-xmelt));
		if(nueff<nmin) nueff=nmin;
		if(nueff>nmax) nueff=nmax;
		/* Ductile stress calc check */
		sduct=nueff*2.0*epsin;
		if(sduct<smin && epsin) {nueff=0.5*smin/epsin; sduct=smin;}
		if(sduct>smax) {nueff=0.5*smax/epsin; sduct=smax;}
		eps[24]=nueff;
		/* Shear modulus */
		eps[27]=markgg[mm2];
		}
	/**/
	/* Heat capacity */
	eps[25]=markcp[mm2]*xmelt+markcp[mm2-20]*(1.0-xmelt);
	/**/
	/* heat conductivity */
	eps[26]=((markkt[mm2]+markkf[mm2]/(mtk+77.0))*exp(markkp[mm2]*mpb))*xmelt+((markkt[mm2-20]+markkf[mm2-20]/(mtk+77.0))*exp(markkp[mm2-20]*mpb))*(1.0-xmelt);
	/**/
	/* Additional melting adiabatic term, heat capacity */
	if(xmelt>0 && xmelt<1.0)
		{
		/* Melting adiabatic term: alm=-ro*(dHlat/dP)/T */
		/* Numerical differentiation */
		dmpb=mpb*0.001;
		meltpart1(mtk,mpb-dmpb,mm2);
		ival=eps[22];
		meltpart1(mtk,mpb+dmpb,mm2);
		ival-=eps[22];
		ival*=eps[23]/(mtk*2.0*dmpb*1e+5);
		eps[20]+=ival;
		/**/
		/* Melting heat capacity term: cpm=dHlat/dT */
		/* Numerical differentiation */
		dmtk=1.0;
		meltpart1(mtk+dmtk,mpb,mm2);
		ival=eps[22];
		meltpart1(mtk-dmtk,mpb,mm2);
		ival-=eps[22];
		ival/=2.0*dmtk;
		eps[25]+=ival;
/*
printf("Meltpart: %ld %d %e %e %e  %e %e %e %e %e %e",mm1,mm2,mtk,mpb,xmelt,eps[20],eps[22],eps[23],eps[24],eps[25],eps[26]);getchar();
*/
		}
	return 1.0;
	}
eps[20]=eps[21]=eps[22]=eps[23]=eps[24]=eps[25]=eps[26]=0;
return 0;
}
/* Rock to rock+melt transformation */





/* Calculation of Vx,Vy,EPSxx,EPSyy,EPSxy,SIGxx,SIGyy,SIGxy,T,T0,T1,T2,P for current location by Linear Interpolation */
/* Staggered Nodes num */
/*   [0]                [3]                [6]   */
/*  T0,xy0    Vy0     T3,xy3     Vy3             */
/*                                               */
/*   Vx0    P4,xx4,yy4  Vx3    P7,xx7,yy7        */
/*                                               */
/*   [1]                [4]                [7]   */
/*  T,xy1     Vy1     T4,xy4     Vy4             */
/*                                               */
/*   Vx1    P5,xx5,yy5  Vx4    P8,xx8,yy8        */
/*                                               */
/*   [2]                [5]                [8]   */
/*                                               */
/*                                               */
void allinter(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* Buffer clear */
for (m1=0;m1<=15;m1++) eps[m1]=0;
/**/
/**/
/**/
/* T interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* T Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[0]+=ival*tk0[m3];
	eps[1]+=ival*tk1[m3];
	eps[2]+=ival*tk[m3];
	eps[3]+=ival*tk2[m3];
	}
/**/
/* Wt for nodes save for T */
wn[2]=m1min; wn[3]=m1max;
for (m1=m1min;m1<=m1max;m1++)
	{
	cn[m1-m1min][5]=cn[m1-m1min][1];
	}
wn[4]=m2min; wn[5]=m2max;
for (m2=m2min;m2<=m2max;m2++)
	{
	cn[m2-m2min][4]=cn[m2-m2min][0];
	}
/* End T interpolation ------------------------ */
/**/
/**/
/**/
/* SIGxy,EPSxy, SIGxy*EPSxy interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(m1min<1) m1min=1; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(m2min<1) m2min=1; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* SIGxy,EPSxy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[4]+=ival*exy[m3];
	eps[5]+=ival*sxy[m3];
	eps[13]+=ival*sxy[m3]*sxy[m3]/(2.0*nu[m3]);
	}
/* End SIGxy,EPSxy interpolation ------------------------ */
/**/
/**/
/**/
/* SIGxx,EPSxx,SIGyy,EPSyy,P, SIGxx*EPSxx, SIGyy*EPSyy interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* SIGxx,EPSxx,SIGyy,EPSyy,P Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[6]+=ival*exx[m3];
	eps[7]+=ival*sxx[m3];
	eps[10]+=ival*pr[m3];
	eps[14]+=ival*sxx[m3]*sxx[m3]/(2.0*nd[m3]);
	}
/* Pressure as function of depth from erosion level calc */
/*
depthp(x,y);eps[10]=eps[50];
*/
/* End SIGxx,EPSxx,SIGyy,EPSyy,P interpolation ------------------------ */
/**/
/**/
/**/
/* Vx interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20; if(y<(gy[m20]+gy[m20+1])/2.0) m2min-=1;
if(m2min<0) m2min=0; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<0) m2min=0;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,+1);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[11]+=ival*vx[m3];
	}
/* End Vx interpolation ------------------------ */
/**/
/**/
/**/
/* Vy interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x<(gx[m10]+gx[m10+1])/2.0) m1min-=1;
if(m1min<0) m1min=0; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<0) m1min=0;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,+1,0);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[12]+=ival*vy[m3];
	}
/* End Vy interpolation ------------------------ */
/**/
/**/
/**/
/*
printf("eps %e %e ",m1,m2,e,n); getchar();
*/
}
/* Calculation of Vx,Vy,EPSxx,EPSyy,EPSxy,SIGxx,SIGyy,SIGxy,T,T0,T1,T2,P for current location by Linear Interpolation */





/* Weights for horisontal and vertical nodes calculation for marker interpolation */ 
void nodewt(long int m1min, long int m1max, long int m2min, long int m2max, double x, double y, int ynx, int yny)
/* m1min,m1max, m2min,m2max - node X,Y number limits */
/* x,y - curent pont coordinates */
/* ynx, yny - Type of shifts: No(0), Back(-1), Forw(1) */
{
/* Eyy vertical position */
long int m3;
int nx,ny;
/**/
/**/
/**/
/* Weigths in horizontal directions */
/* Load distances to xn[] */
if(ynx<0) 
	{
	for (m3=m1min;m3<=m1max;m3++)
		{
		xn[m3-m1min]=(gx[m3]+gx[m3-1])/2.0;
		}
	}
if(ynx==0) 
	{
	for (m3=m1min;m3<=m1max;m3++)
		{
		xn[m3-m1min]=gx[m3];
		}
	}
if(ynx>0) 
	{
	for (m3=m1min;m3<=m1max;m3++)
		{
		xn[m3-m1min]=(gx[m3]+gx[m3+1])/2.0;
		}
	}
/**/
/* Calc maximal position in xn[] */
nx=(int)(m1max-m1min);
/**/
/* Calc coefficients for horizontal direction */
fdweight(nx,0,x);
/**/
/* Reload horizontal coefficients to cn[] */
for (m3=0;m3<=nx;m3++)
	{
	cn[m3][1]=cn[m3][0];
	}
/**/
/**/
/**/
/* Weigths in vertical directions */
/* Load distances to xn[] */
if(yny<0) 
	{
	for (m3=m2min;m3<=m2max;m3++)
		{
		xn[m3-m2min]=(gy[m3]+gy[m3-1])/2.0;
		}
	}
if(yny==0) 
	{
	for (m3=m2min;m3<=m2max;m3++)
		{
		xn[m3-m2min]=gy[m3];
		}
	}
if(yny>0) 
	{
	for (m3=m2min;m3<=m2max;m3++)
		{
		xn[m3-m2min]=(gy[m3]+gy[m3+1])/2.0;
		}
	}
/**/
/* Calc maximal position in xn[] */
ny=(int)(m2max-m2min);
/**/
/* Calc coefficients for horizontal direction */
fdweight(ny,0,y);
/**/
/**/
/**/
}
/* Weights for horisontal and vertical nodes calculation for marker interpolation */ 



/* Calculation of Vx,Vy by Interpolation */
void allinterv1(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* Vx interpolation ------------------------ */
/* Buffer clear */
eps[11]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20; if(y<(gy[m20]+gy[m20+1])/2.0) m2min-=1;
if(m2min<0) m2min=0; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<0) m2min=0;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,+1);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[11]+=ival*vx[m3];
	}
/* End Vx interpolation ------------------------ */
/**/
/**/
/**/
/* Vy interpolation ------------------------ */
/* Buffer clear */
eps[12]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x<(gx[m10]+gx[m10+1])/2.0) m1min-=1;
if(m1min<0) m1min=0; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<0) m1min=0;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,+1,0);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[12]+=ival*vy[m3];
	}
/* End Vy interpolation ------------------------ */
/*
printf("eps %e %e ",m1,m2,e,n); getchar();
*/
}
/* Calculation of Vx,Vy by Interpolation */



/* Calculation of EPSij, SIGij, EPSij*SIGij, P by Interpolation */
void allintere(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* EPSxy interpolation ------------------------ */
/* Buffer clear */
eps[4]=eps[5]=eps[13]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(m1min<1) m1min=1; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(m2min<1) m2min=1; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* EPSxy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[4]+=ival*exy[m3];
	eps[5]+=ival*sxy[m3];
	eps[13]+=ival*sxy[m3]*sxy[m3]/(2.0*nu[m3]);
	}
/* End SIGxy,EPSxy interpolation ------------------------ */
/**/
/**/
/**/
/* EPSxx,EPSyy, P interpolation ------------------------ */
/* Buffer clear */
eps[6]=eps[7]=eps[10]=eps[14]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* EPSxx,EPSyy,P Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[6]+=ival*exx[m3];
	eps[7]+=ival*sxx[m3];
	eps[10]+=ival*pr[m3];
	eps[14]+=ival*sxx[m3]*sxx[m3]/(2.0*nd[m3]);
	}
/* End EPSxx,EPSyy,P interpolation ------------------------ */
/* Pressure as function of depth from erosion level calc */
/*
depthp(x,y);eps[10]=eps[50];
*/
/*
printf("eps %e %e ",m1,m2,e,n); getchar();
*/
}
/* Calculation of EPSij, SIGij, EPSij*SIGij, P by Interpolation */



/* Calculation of Vx,Vy, EPSxx*SIGxx,EPSyy*SIGyy,EPSxy*SIGxy by Interpolation */
void allinters(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* SIGxy*EPSxy interpolation ------------------------ */
/* Buffer clear */
eps[13]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(m1min<1) m1min=1; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(m2min<1) m2min=1; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* SIGxy,EPSxy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[13]+=ival*sxy[m3]*sxy[m3]/(2.0*nu[m3]);
	}
/* End SIGxy*EPSxy interpolation ------------------------ */
/**/
/**/
/**/
/* SIGxx*EPSxx, SIGyy*EPSyy interpolation ------------------------ */
/* Buffer clear */
eps[14]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* SIGxx,EPSxx,SIGyy,EPSyy,P Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[14]+=ival*sxx[m3]*sxx[m3]/(2.0*nd[m3]);
	}
/* End SIGxx*EPSxx,SIGyy*EPSyy interpolation ------------------------ */
/**/
/**/
/**/
/* Vx interpolation ------------------------ */
/* Buffer clear */
eps[11]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20; if(y<(gy[m20]+gy[m20+1])/2.0) m2min-=1;
if(m2min<0) m2min=0; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<0) m2min=0;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,+1);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[11]+=ival*vx[m3];
	}
/* End Vx interpolation ------------------------ */
/**/
/**/
/**/
/* Vy interpolation ------------------------ */
/* Buffer clear */
eps[12]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x<(gx[m10]+gx[m10+1])/2.0) m1min-=1;
if(m1min<0) m1min=0; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<0) m1min=0;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,+1,0);
/**/
/* Vy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[12]+=ival*vy[m3];
	}
/* End Vy interpolation ------------------------ */
/*
printf("eps %e %e ",m1,m2,e,n); getchar();
*/
}
/* Calculation of Vx,Vy, EPSxx*SIGxx,EPSyy*SIGyy,EPSxy*SIGxy by Interpolation */



/* Calculation of Vx,Vy, EPSxx,EPSyy,EPSxy, SPINxy by Interpolation */
void allinteri(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival,xrat;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Check weighting for interpolation */
xrat=2.0/3.0;
if(x<(gx[0]+gx[1])/2.0) xrat=1.0;
if(x>(gx[xnumx-2]+gx[xnumx-1])/2.0) xrat=1.0;
if(y<(gy[0]+gy[1])/2.0) xrat=1.0;
if(y>(gy[ynumy-2]+gy[ynumy-1])/2.0) xrat=1.0;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* EPSxy, SPINxy interpolation ------------------------ */
/* Buffer clear */
eps[4]=eps[30]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(m1min<1) m1min=1; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(m2min<1) m2min=1; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* SIGxy,EPSxy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[4]+=ival*exy[m3];
	eps[30]+=ival*esp[m3];
	}
/* End SIGxy*EPSxy interpolation ------------------------ */
/**/
/**/
/**/
/* EPSxx, EPSyy, P interpolation ------------------------ */
/* Buffer clear */
eps[6]=eps[10]=eps[11]=eps[12]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* SIGxx,EPSxx,SIGyy,EPSyy,P Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[6]+=ival*exx[m3];
	eps[10]+=ival*pr[m3];
	eps[11]+=ival*(vx[m3-1]+vx[m3-ynumy-1])*0.5*(1.0-xrat);
	eps[12]+=ival*(vy[m3-ynumy]+vy[m3-ynumy-1])*0.5*(1.0-xrat);
	}
/* End SIGxx*EPSxx,SIGyy*EPSyy interpolation ------------------------ */
/*
depthp(x,y);eps[10]=eps[50];
*/
/**/
/**/
/**/
/* Vx interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20; if(y<(gy[m20]+gy[m20+1])/2.0) m2min-=1;
if(m2min<0) m2min=0; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<0) m2min=0;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,+1);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[11]+=ival*vx[m3]*xrat;
	}
/* End Vx interpolation ------------------------ */
/**/
/**/
/**/
/* Vy interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x<(gx[m10]+gx[m10+1])/2.0) m1min-=1;
if(m1min<0) m1min=0; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<0) m1min=0;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,+1,0);
/**/
/* Vy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[12]+=ival*vy[m3]*xrat;
	}
/* End Vy interpolation ------------------------ */
/*
printf("eps %e %e ",m1,m2,e,n); getchar();
*/
}
/* Calculation of Vx,Vy, EPSxx,EPSyy,EPSxy, SPINxy by Interpolation */




/* Calculation of Vx,Vy by Interpolation */
void allinterv(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival,xrat;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Check weighting for interpolation */
xrat=2.0/3.0;
if(x<(gx[0]+gx[1])/2.0) xrat=1.0;
if(x>(gx[xnumx-2]+gx[xnumx-1])/2.0) xrat=1.0;
if(y<(gy[0]+gy[1])/2.0) xrat=1.0;
if(y>(gy[ynumy-2]+gy[ynumy-1])/2.0) xrat=1.0;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* Vx interpolation ------------------------ */
/* Buffer clear */
eps[11]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20; if(y<(gy[m20]+gy[m20+1])/2.0) m2min-=1;
if(m2min<0) m2min=0; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<0) m2min=0;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,+1);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[11]+=ival*vx[m3]*xrat;
	}
/* End Vx interpolation ------------------------ */
/**/
/**/
/**/
/* Vy interpolation ------------------------ */
/* Buffer clear */
eps[12]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x<(gx[m10]+gx[m10+1])/2.0) m1min-=1;
if(m1min<0) m1min=0; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<0) m1min=0;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,+1,0);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[12]+=ival*vy[m3]*xrat;
	}
/* End Vy interpolation ------------------------ */
/*
printf("eps %e %e ",m1,m2,e,n); getchar();
*/
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* Vx, Vy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	/* 0    2 */
	/* 1    3 */
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[11]+=ival*(vx[m3-1]+vx[m3-ynumy-1])*0.5*(1.0-xrat);
	eps[12]+=ival*(vy[m3-ynumy]+vy[m3-ynumy-1])*0.5*(1.0-xrat);
	}
}
/* Calculation of Vx,Vy by Interpolation */





/* Calculation of T,T0 for current location by Interpolation */
void allintert(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* T interpolation ------------------------ */
/* Buffer clear */
eps[2]=eps[3]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* T Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[2]+=ival*tk[m3];
	eps[3]+=ival*tk2[m3];
	}
/**/
/* Wt for nodes save for T */
wn[2]=m1min; wn[3]=m1max;
for (m1=m1min;m1<=m1max;m1++)
	{
	cn[m1-m1min][5]=cn[m1-m1min][1];
	}
wn[4]=m2min; wn[5]=m2max;
for (m2=m2min;m2<=m2max;m2++)
	{
	cn[m2-m2min][4]=cn[m2-m2min][0];
	}
/* End T interpolation ------------------------ */
/*
printf("eps %e %e ",m1,m2,e,n); getchar();
*/
}
/* Calculation of T,T0 for current location by Interpolation */



/* Calculation of SIGij by Interpolation */
void allinterd(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/* T interpolation ------------------------ */
/* Buffer clear */
eps[2]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* SIGij, T Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[2]+=ival*tk[m3];
	}
/**/
/* Wt for nodes save for SIGij old */
wn[2]=m1min; wn[3]=m1max;
for (m1=m1min;m1<=m1max;m1++)
	{
	cn[m1-m1min][5]=cn[m1-m1min][1];
	}
wn[4]=m2min; wn[5]=m2max;
for (m2=m2min;m2<=m2max;m2++)
	{
	cn[m2-m2min][4]=cn[m2-m2min][0];
	}
/* End SIGij old interpolation ------------------------ */
/**/
/**/
/**/
/* SIGxy interpolation ------------------------ */
/* Buffer clear */
eps[4]=eps[5]=eps[31]=eps[34]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(m1min<1) m1min=1; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(m2min<1) m2min=1; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* EPSxy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[4]+=ival*exy[m3];
	eps[5]+=ival*sxy[m3];
	eps[31]+=ival*sxye[m3];
	eps[34]+=ival*exye[m3];
	}
/* End SIGxy interpolation ------------------------ */
/**/
/**/
/**/
/* SIGxx,SIGyy interpolation ------------------------ */
/* Buffer clear */
eps[6]=eps[7]=eps[10]=eps[32]=eps[33]=eps[35]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* SIGxx,SIGyy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[6]+=ival*exx[m3];
	eps[7]+=ival*sxx[m3];
	eps[10]+=ival*pr[m3];
	eps[32]+=ival*sxxe[m3];
	eps[33]+=ival*sppe[m3];
	eps[35]+=ival*exxe[m3];
	}
/* End SIGxx,SIGyy interpolation ------------------------ */
/*
printf("eps %e %e ",m1,m2,e,n); getchar();
*/
}
/* Calculation of SIGij by Interpolation */



/* Calculation of P by Interpolation */
void allinterp(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* P interpolation ------------------------ */
/* Buffer clear */
eps[10]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* EPSxx,EPSyy,P Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[10]+=ival*pr[m3];
	}
/* End EPSxx,EPSyy,P interpolation ------------------------ */
/*
printf("eps %e %e ",m1,m2,e,n); getchar();
*/
}
/* Calculation of P by Interpolation */



/* Calculation of P from the depth below the surface */
void depthp(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m10;
/* en-NormalisedDistance */
double e,n,ival;
/**/
/* Up Left Node X,Y Num */
m10=wn[0];
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Pressure as function of depth from erosion level calc */
/* Relativ Normalized coord Calc */
e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
n=(e*ep[m10+1]+(1.0-e)*ep[m10]);
ival=3300.0*(y-n)*9.80665;
if (ival<1e+5) ival=1e+5;
eps[50]=ival;
/*
*/
/*
printf("eps %e %e ",m1,m2,e); getchar();
*/
}
/* Calculation of P from the depth below the surface */



/* Melt fraction, latent heat calculation */
void meltpart1(double mtk, double mpb, int mm2)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm1 - mark number */
/* mm2 - mark type */
/* yn  - type of calculation: 0 - Ro, 1 - Nu, 2 - Cp, 3 - kt */
{
/* Val buffer */
double xmelt=0,hlatent=0,ival;
long int m1;
double ykm=mpb*3e-3,ts=0,tl=0;
/**/
/**/
/**/
/* Calculate melt fraction using marker type */
if (ykm>0)
switch(mm2)
	{
	/* Sediments: latent heat 300 kJ/kg (Bittner & Schmeling, 1995) */
	case  3:
	case  4:
	case  5:
	case 17:
	case 23:
	case 24:
	case 25:
	case 37:
	/* Wet Solidus Temperature, Johannes, 1985, Poli & Schmidt, 2002 */
	if (ykm<36.0) 
		{
		ts=889.0+536.6/(ykm+1.609)+18.21/(ykm+1.609)/(ykm+1.609);
		}
	else
		{
		ts=831.3+2.0*ykm;
		}
	/* Dry Granite Liquidus, Johannes, 1985 */
	tl=1262.0+3.0*ykm;
	hlatent=300000.0;
	break;
	/**/
	/* Basalt, Gabbro: latent heat 380 kJ/kg (Bittner & Schmeling, 1995) */
	case 7:
	case 8:
	case 16:
	case 27:
	case 28: 
	case 36: 
	case  6:
	case 18:
	case 26:
	case 38:
	/* Wet solidus, Schmidt & Poli, 1998  */
	if (ykm<48.0) 
		{
		ts=972.6-2111.0/(ykm+10.63)+70033.0/(ykm+10.63)/(ykm+10.63);
		}
	else
		{
		ts=935.4+0.1162*ykm+0.006937*ykm*ykm;
		}
	/* Dry Toleitic Basalt Liquidus, Hess, 1989 */
	tl=1423.15+3.5*ykm;
	hlatent=380000.0;
	break;
	/**/
	/* Peridotite: latent heat 400 kJ/kg Turcotte & Schubert, 1982, p.171 */
	case 11:
	case 34:
	/* Wet solidus, Schmidt & Poli, 1998  */
	if (ykm<72.0) 
		{
		ts=1239.8+1493.0/(ykm+9.701);
		}
	else
		{
		ts=1266.3-0.3948*ykm+0.003893*ykm*ykm;
		}
	/* Dry Peridotite Liquidus, Hess, 1989 */
	tl=2073.15+3.8*ykm;
	hlatent=400000.0;
	break;
	/**/
	/* Other rocks - No melting */
	default:
	break;
	}
/**/
/* Melt fraction, latent heat calculation */
eps[21]=eps[22]=0;
if(tl)
	{
	/* Melt fraction calc, check */
	xmelt=(mtk-ts)/(tl-ts);
	if(xmelt<0) xmelt=0;
	if(xmelt>1.0) xmelt=1.0;
	eps[21]=xmelt;
	/* Latent heat calc */
	hlatent*=xmelt;
	eps[22]=hlatent;
/*
if(mm2<20 && xmelt) {printf("Meltpart1: %d %e %e  %e %e",mm2,mtk,mpb,xmelt,hlatent);getchar();}
*/
	}
/**/
}
/* Melt fraction, latent heat calculation */



/* Antigorite dehydration depth */
double dehydration(int yn)
/* yn - Upper(0)/Lower boundary */
{
/* Val buffer */
double x0,x1,y0,y1,k0,k1,t0,t1,sy1,ival,yshift;
long int m1;
/**/
/**/
/**/
/* Set Upper/Lower boundary */
if(yn) yshift=5000.0; else yshift=-5000.0;
/**/
/* Serch for melting depth along lower hydration surface */
for (m1=1;m1<xnumx;m1++)
	{
	/* TK calc along the lower hydr surface */
	x0=gx[m1-1];
	y0=ep[m1+xnumx*3-1]+yshift;
	allintert(x0,y0); t0=eps[2];
	x1=gx[m1];
	y1=ep[m1+xnumx*3]+yshift;
	allintert(x1,y1); t1=eps[2];
	/* Oceanic crust top */
	sy1=ep[m1];
	if((y0-sy1)>0 && (y1-sy1)>0)
		{
		/* Atg stability field after Schmidt and Poli, 1998 */
		if((y0-sy1)>63000.0)
			{
			k0=1013.17699-0.060387633e-3*(y0-sy1)-0.004289442e-6*(y0-sy1)*(y0-sy1);
			}
		else
			{
			k0=751.490422+6.00773668e-3*(y0-sy1)-0.034690759e-6*(y0-sy1)*(y0-sy1);
			}
		if((y1-sy1)>63000.0)
			{
			k1=1013.17699-0.060387633e-3*(y1-sy1)-0.004289442e-6*(y1-sy1)*(y1-sy1);
			}
		else
			{
			k1=751.490422+6.00773668e-3*(y1-sy1)-0.034690759e-6*(y1-sy1)*(y1-sy1);
			}
		/* Antigorite stability depth calc */
		if(t1>=k1)
			{
			ival=y1-(y1-y0)/((t1-k1)+(k0-t0))*(t1-k1);
			/* Set Upper/Lower boundary */
			if(yn) dehydrmax=ival; else dehydrmin=ival;
			return 0;
			}
		}
	}
basalty=ysize;
if(yn) dehydrmax=ysize; else dehydrmin=ysize;
return 0;
}
/* Antigorite dehydration depth */



/* Hydration front progress */
double hydration()
{
/* Val buffer */
double vsubd,tgasubd,afiltr,bfiltr,rohidr,h2ohidr,mh2o,vfiltr,yfiltr,tsubd,dydx,vxs,vys,sy1,sy2,sy3,sy4,sy5,x0,y0,x1,y1,vx1,vy1;
double hytimesum,hytimesum0;
double serp1,serp2,serp3;
long int m1,m2;
/**/
/**/
/* Vx, Vy Subduction velosity */
vxs=vx[25];
vys=vy[25];
if(vxs<=0) return 0;
/* Upper surface for hydration rate calculation: e.g., oceanic crust entrance top */
sy1=10000.0;
/* Min, Max level for material displacement account: e.g., continental crust bottom */
sy2=14000.0;
sy3=400000.0;
/* Max vertical distance between next and previose and next nodes of the upper hydration surface */
sy4=2000.0;
sy5=4000.0;
/* Oceanic crust age */
/* 1year=3.15576*10^7sek */
tsubd=0.0e+6*3.15576e+7;
/* Check subduction time */
if((timesum-tsubd)<=0) return 0;
/* Max alteration depth */
yfiltr=dehydrmin;
/*
yfiltr=125000.0;
if(yfiltr>150000.0) yfiltr=150000.0;
if(yfiltr<45000.0) yfiltr=45000.0;
*/
/* Filtration velosity constants */
vsubd=pow((vys*vys+vxs*vxs),0.5);
afiltr=vsubd*0.05;
if (yfiltr>100000.0) afiltr=vsubd*0.05*100000.0/yfiltr;
/* Depths and hydration rate due to serpentine */
serp1=dehydrmin;
serp2=dehydrmax;
serp3=vsubd*0.10;
if ((dehydrmax-dehydrmin)>50000.0) serp3=vsubd*0.10*50000.0/(dehydrmax-dehydrmin);
/**/
/*
rohidr=markro[13];
mh2o=0.50e+6;
bfiltr=0.0;
h2ohidr=(3300.0-rohidr)*0.07/420.0;
tgasubd=1.0;
vys=vsubd*pow(0.5,0.5);
afiltr=vsubd*2.0*tgasubd*mh2o/(2.0-bfiltr)/rohidr/h2ohidr/(yfiltr-sy1);
*/
printf("Hydration: Zlim=%3.2e Vs=%3.2e Vh=%3.2e Vh/Vs=%3.2e Zsermin=%3.2e Zsermax=%3.2e Vhser=%3.2e Vhser/Vs=%3.2e \n",yfiltr,vsubd,afiltr,afiltr/vsubd,serp1,serp2,serp3,serp3/vsubd);
/*
printf("Subduction %e %e",vxs,vys); getchar();
getchar();
*/
/**/
/**/
/* Hydration Solution Cycle ------------------------------------------ */
hytimesum=0;
hytimesum0=timestep;
do
{
/* Save old cycle results */
for (m1=0;m1<xnumx;m1++)
	{
	ep0[xnumx*2+m1]=ep[xnumx*2+m1];
	ep0[xnumx*3+m1]=ep[xnumx*3+m1];
	}
/**/
/**/
/**/
/* Initial timestep definition */
timestep=hytimesum0-hytimesum;
/**/
/**/
/**/
/* Hydration timestep definition using material velosity field */
for (m1=0;m1<xnumx;m1++)
	{
	/* Calc horisontal,vertical Coordinate for Lower surface and Hydration Surface */
	x0=x1=gx[m1];
	y0=ep0[m1+xnumx*3];
	y1=ep0[m1+xnumx*2];
	/**/
	/* HYDRATION SURFACE */
	/* Calc material velocity on the Lower Surface using velosity field */
	vx1=0; 
	vy1=0;
	/* Material displacement starting from min level */
	if(y1>=sy2 && y1<=sy3)
		{
		allinterv(x1,y1);
		vx1=eps[11];
		vy1=eps[12];
		}
	/* Calc x derivative of y position of the Surface using upwind differences */
	dydx=0;
	if(vx1>0 && m1>0)
		{
		timestep=MINV(timestep,(gx[m1]-gx[m1-1])/vx1);
/*
printf("111 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1-1])/vx1,hytimesum0,hytimesum,hytimesum0-hytimesum,timestep);getchar();
*/
		}
	if(vx1<0 && m1<xnumx-1)
		{
		timestep=MINV(timestep,(gx[m1]-gx[m1+1])/vx1);
/*
printf("222 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1+1])/vx1,hytimesum0,hytimesum,hytimesum0-hytimesum,timestep);getchar();
*/
		}
	/* Check vertical timestep */
	if(vy1)
		{
		/* Horizontal line num definition */
		m2=m2serch(y1);
		/* Check timestep */
		timestep=MINV(timestep,(gy[m2+1]-gy[m2])/ABSV(vy1));
/*
printf("333 %ld  %e %e %e %e %e %e %e",m2,vx1,vy1,(gy[m2+1]-gy[m2])/ABSV(vy1),hytimesum0,hytimesum,hytimesum0-hytimesum,timestep);getchar();
*/
		}
	/**/
	/**/
	/* LOWER SURFACE */
	/* Calc material velocity on hydration front using velosity field */
	vx1=0; 
	vy1=0;
	/* Material displacement starting from min level */
	if(y0>=sy2 && y0<=sy3)
		{
		allinterv(x0,y0);
		vx1=eps[11];
		vy1=eps[12];
		}
	/* Calc x derivative of y position of the Surface using upwind differences */
	dydx=0;
	if(vx1>0 && m1>0)
		{
		timestep=MINV(timestep,(gx[m1]-gx[m1-1])/vx1);
/*
printf("444 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1-1])/vx1,hytimesum0,hytimesum,hytimesum0-hytimesum,timestep);getchar();
*/
		}
	if(vx1<0 && m1<xnumx-1)
		{
		timestep=MINV(timestep,(gx[m1]-gx[m1+1])/vx1);
/*
printf("555 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1+1])/vx1,hytimesum0,hytimesum,hytimesum0-hytimesum,timestep);getchar();
*/
		}
	/* Check vertical timestep */
	if(vy1)
		{
		/* Horizontal line num definition */
		m2=m2serch(y0);
		/* Check timestep */
		timestep=MINV(timestep,(gy[m2+1]-gy[m2])/ABSV(vy1));
/*
printf("666 %ld  %e %e %e %e %e %e %e",m2,vx1,vy1,(gy[m2+1]-gy[m2])/ABSV(vy1),hytimesum0,hytimesum,hytimesum0-hytimesum,timestep);getchar();
*/
		}
	}
/*
printf("777 %e %e %e %e",hytimesum0,hytimesum,hytimesum0-hytimesum,timestep);getchar();
*/
/**/
/**/
/**/
/* Displace Lower and Upper Hydration boundaries */
for (m1=0;m1<xnumx;m1++)
	{
	/* Location of the Lower, Upper hydration boundary */
	x0=x1=gx[m1];
	y0=ep0[m1+xnumx*3];
	y1=ep0[m1+xnumx*2];
/*
printf("%e %e %e %e \n",x1,y1,timestep,ep0[xnumx*2+m1]);
printf("%e %e %e %e \n",x1,y1,timestep,ep0[xnumx*2+m1]);
getchar();
*/
	/**/
	/**/
	/**/
	/* Displacement of the Lower Surface */
	/* Calc material velocity on Lower Surface using velosity field */
	vx1=0; 
	vy1=0;
	/* Material displacement starting from min level */
	if(y0>=sy2 && y0<=sy3)
		{
		allinterv(x0,y0);
		vx1=eps[11];
		vy1=eps[12];
		}
	/* Calc x derivative of y position of hydration boundary using upwind differences */
	dydx=0;
	if(vx1>0 && m1>0)
		{
		dydx=(ep0[xnumx*3+m1]-ep0[xnumx*3+m1-1])/(gx[m1]-gx[m1-1]);
/*
printf("AAA %e ",dydx);getchar();
*/
		}
	if(vx1<0 && m1<xnumx-1)
		{
		dydx=(ep0[xnumx*3+m1+1]-ep0[xnumx*3+m1])/(gx[m1+1]-gx[m1]);
/*
printf("BBB %e ",dydx);getchar();
*/
		}
	/* Recalc, Check new hydration front position */
	ep[xnumx*3+m1]+=timestep*(vy1-dydx*vx1);
/*
if (ep[xnumx*2+m1]>y1) ep[xnumx*2+m1]=y1;
printf("%e %e %e %e %e %e %e %e ",x1,y1,timestep,vfiltr,vx1,vy1,dydx,ep[xnumx*2+m1]);
printf("%e \n",ep[xnumx*2+m1]);
printf("%e %e %e %e %e %e %e %e",x1,y1,timestep,vfiltr,vx1,vy1,dydx,ep[xnumx*2+m1]);getchar();
*/
	/**/
	/**/
	/**/
	/* Displacement of the upper Hydration Surface */
	/* Hydration velosity for current y0 location along lower hydration surface */
	vfiltr=0;
	if(y0>=sy1 && y0<=yfiltr) vfiltr=afiltr*(1.0-bfiltr*(y0-sy1)/(yfiltr-sy1));
	if(y0>=serp1 && y0<=serp2) vfiltr=serp3;
	/* Calc material velocity on hydration front using velosity field */
	vx1=0; 
	vy1=0;
	/* Material displacement starting from min level */
	if(y1>=sy2 && y1<=sy3)
		{
		allinterv(x1,y1);
		vx1=eps[11];
		vy1=eps[12];
		}
	/* Calc x derivative of y position of hydration boundary using upwind differences */
	dydx=0;
	if(vx1>0 && m1>0)
		{
		dydx=(ep0[xnumx*2+m1]-ep0[xnumx*2+m1-1])/(gx[m1]-gx[m1-1]);
/*
printf("AAA %e ",dydx);getchar();
*/
		}
	if(vx1<0 && m1<xnumx-1)
		{
		dydx=(ep0[xnumx*2+m1+1]-ep0[xnumx*2+m1])/(gx[m1+1]-gx[m1]);
/*
printf("BBB %e ",dydx);getchar();
*/
		}
	/* Recalc, Check new hydration front position */
	ep[xnumx*2+m1]+=timestep*(-vfiltr+vy1-dydx*vx1);
	if (ep[xnumx*2+m1]>ep[xnumx*3+m1]) ep[xnumx*2+m1]=ep[xnumx*3+m1];
	if (m1 && ep[xnumx*2+m1]<(ep[xnumx*2+m1-1]-sy4)) ep[xnumx*2+m1]=(ep[xnumx*2+m1-1]-sy4);
	if (m1<xnumx-1 && ep[xnumx*2+m1]<(ep[xnumx*2+m1+1]-sy5)) ep[xnumx*2+m1]=(ep[xnumx*2+m1+1]-sy5);
/*
printf("%e %e %e %e %e %e %e %e ",x1,y1,timestep,vfiltr,vx1,vy1,dydx,ep[xnumx*2+m1]);
printf("%e %e %e %e %e %e %e %e %e",vys,x1,y1,timestep,vfiltr,vx1,vy1,dydx,ep[xnumx*2+m1]);getchar();
printf("%e \n",ep[xnumx*2+m1]);
*/
	if(y0>sy3) break;
	}
/**/
/**/
/**/
/* Add Hydration step */
hytimesum+=timestep;
/**/
/**/
/**/
/* Print Results */
if (printmod)
	{
	printf("\n HYDRATION STEP = %e YEARS    HYDRATION TIME = %e YEARS \n",timestep/3.15576e+7,hytimesum/3.15576e+7);
/*
getchar();
*/
	}
}
while(hytimesum<hytimesum0);
/* Restore timestep */
timestep=hytimesum0;
return 0;
}
/* Hydration front progress */



/* Errosion Surface progress */
void erosion1()
{
/* Val buffer */
double v0,v1,dydt,x1,y1,x11,y11,x2,y2,x21,y21,vx1,vy1,vx2,vy2,dydx,dy;
double ertimesum,ertimesum0;
long int m1,m2;
/**/
/* Erosion Solution Cycle ------------------------------------------ */
ertimesum=0;
ertimesum0=timestep;
do
{
/* Save old cycle results */
for (m1=0;m1<xnumx;m1++)
	{
	ep0[m1]=ep[m1];
	ep0[xnumx+m1]=ep[xnumx+m1];
	}
/**/
/**/
/**/
/* Initial timestep definition */
timestep=ertimesum0-ertimesum;
/**/
/**/
/**/
/* Erosion timestep definition using material velosity field */
for (m1=0;m1<xnumx;m1++)
	{
	/* Calc horisontal Coordinate */
	x1=gx[m1];
	/**/
	/* EROSION SURFACE */
	/* Calc material velocity on the Surface using velosity field */
	allinterv(x1,ep0[m1]);
	vx1=eps[11];
	vy1=eps[12];
	/* Check horizontal timestep */
	/* Calc x derivative of y position of the Surface using upwind differences */
	if(vx1>0 && m1>0)
		{
		timestep=MINV(timestep,0.5*(gx[m1]-gx[m1-1])/vx1);
/*
printf("111 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1-1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
		}
	if(vx1<0 && m1<xnumx-1)
		{
		timestep=MINV(timestep,0.5*(gx[m1]-gx[m1+1])/vx1);
/*
printf("222 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1+1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
		}
	/* Check vertical timestep */
	if(vy1)
		{
		/* Horizontal line num definition */
		m2=m2serch(ep0[m1]);
		/* Check timestep */
		timestep=MINV(timestep,0.5*(gy[m2+1]-gy[m2])/ABSV(vy1));
/*
printf("333 %ld  %e %e %e %e %e %e %e",m2,vx1,vy1,(gy[m2+1]-gy[m2])/ABSV(vy1),ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
		}
	/**/
	/**/
	/* INITIAL SURFACE */
	/* Calc material velocity on the Initial Surface using velosity field */
	allinterv(x1,ep0[xnumx+m1]);
	vx1=eps[11];
	vy1=eps[12];
	/* Check horizontal timestep */
	/* Calc x derivative of y position of the Surface using upwind differences */
	if(vx1>0 && m1>0)
		{
		timestep=MINV(timestep,0.5*(gx[m1]-gx[m1-1])/vx1);
/*
printf("444 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1-1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
		}
	if(vx1<0 && m1<xnumx-1)
		{
		timestep=MINV(timestep,0.5*(gx[m1]-gx[m1+1])/vx1);
/*
printf("555 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1+1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
		}
	/* Check vertical timestep */
	if(vy1)
		{
		/* Horizontal line num definition */
		m2=m2serch(ep0[xnumx+m1]);
		/* Check timestep */
		timestep=MINV(timestep,0.5*(gy[m2+1]-gy[m2])/ABSV(vy1));
/*
printf("666 %ld  %e %e %e %e %e %e %e",m2,vx1,vy1,(gy[m2+1]-gy[m2])/ABSV(vy1),ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
		}
	}
/*
printf("777 %e %e %e %e",ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
/**/
/**/
/**/
/* Displace Surface boundary */
/*
for (m1=1;m1<xnumx-1;m1++)
*/
for (m1=0;m1<xnumx;m1++)
	{
	/* EROSION SURFACE */
	/* Calculation of errosion rate */
	v0=0;
	if(ep0[m1]<eroslev)
		{
		v0=eroscon+eroskoe*(eroslev-ep0[m1]);
		}
	/* Calculation of sedimentation rate */
	v1=0;
	if(ep0[m1]>sedilev)
		{
		v1=sedicon+sedikoe*(ep0[m1]-sedilev);
		}
	/* Calc horisontal Coordinate */
	x1=gx[m1];
	y1=ep0[m1];
	/**/
	/* Calc material velocity on the Surface using velosity field */
	allinterv(x1,ep0[m1]);
	vx1=eps[11];
	vy1=eps[12];
	/**/
	/* Erase erosion/sedimentation rate for marginal points */
	if((m1==0 && vx1>0) || (m1==xnumx-1 && vx1<0)) v0=v1=0;
	/**/
	/* Calc advective displacement of two points */
	dydt=vy1*timestep;
	if(m1==0 && vx1>0) dydt=0;
	if(m1==xnumx-1 && vx1<0) dydt=0;
	if(vx1>0 && m1>0)
		{
		x2=gx[m1-1];
		y2=ep0[m1-1];
		allinterv(x2,y2);
		vx2=eps[11];
		vy2=eps[12];
		/* Calc Displacements */
		x11=x1+timestep*vx1;
		y11=y1+timestep*vy1;
		x21=x2+timestep*vx2;
		y21=y2+timestep*vy2;
		/* Calc advected Position of erosion surface */
		dydt=(y11+(y21-y11)/(x21-x11)*(x1-x11))-y1;
/*
printf("AAA %e %e",ep0[m1],dydx);getchar();
*/
		}
	if(vx1<0 && m1<xnumx-1)
		{
		x2=gx[m1+1];
		y2=ep0[m1+1];
		allinterv(x2,y2);
		vx2=eps[11];
		vy2=eps[12];
		/* Calc Displacements */
		x11=x1+timestep*vx1;
		y11=y1+timestep*vy1;
		x21=x2+timestep*vx2;
		y21=y2+timestep*vy2;
		/* Calc advected Position of erosion surface */
		dydt=(y11+(y21-y11)/(x21-x11)*(x1-x11))-y1;
/*
printf("BBB %e %e",ep0[m1],dydx);getchar();
*/
		}
	/* Recalc new Surface position */
	ep[m1]+=dydt+timestep*(v0-v1);
/*
printf("SURFACE %ld %e %e %e %e %e %e %e %e",m1,x1,v0,v1,vx1,vy1,dydx,ep[m1]);getchar();
*/
	/**/
	/**/
	/**/
	/* INITIAL SURFACE */
	/* Initial surface displacement */
	/* Calc material velocity on the Surface using velosity field */
	y1=ep0[xnumx+m1];
	allinterv(x1,ep0[xnumx+m1]);
	vx1=eps[11];
	vy1=eps[12];
	/* Calc x derivative of y position of Initial Surface using upwind differences */
	dydt=vy1*timestep;
	if(m1==0 && vx1>0) dydt=0;
	if(m1==xnumx-1 && vx1<0) dydt=0;
	if(vx1>0 && m1>0)
		{
		x2=gx[m1-1];
		y2=ep0[xnumx+m1-1];
		allinterv(x2,y2);
		vx2=eps[11];
		vy2=eps[12];
		/* Calc Displacements */
		x11=x1+timestep*vx1;
		y11=y1+timestep*vy1;
		x21=x2+timestep*vx2;
		y21=y2+timestep*vy2;
		/* Calc advected Position of initial surface */
		dydt=(y11+(y21-y11)/(x21-x11)*(x1-x11))-y1;
/*
printf("AAA %e ",dydx);getchar();
*/
		}
	if(vx1<0 && m1<xnumx-1)
		{
		x2=gx[m1+1];
		y2=ep0[xnumx+m1+1];
		allinterv(x2,y2);
		vx2=eps[11];
		vy2=eps[12];
		/* Calc Displacements */
		x11=x1+timestep*vx1;
		y11=y1+timestep*vy1;
		x21=x2+timestep*vx2;
		y21=y2+timestep*vy2;
		/* Calc advected Position of initial surface */
		dydt=(y11+(y21-y11)/(x21-x11)*(x1-x11))-y1;
/*
printf("BBB %e ",dydx);getchar();
*/
		}
	/* Recalc new Initial Surface position */
	ep[xnumx+m1]+=dydt;
	/**/
	}
/**/
/**/
/**/
/* Relax EROSION surface */
for (m1=0;m1<xnumx-1;m1++)
	{
	/* Calc x derivative of y position */
	dydx=(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]);
	/* Relax surface for critical slope */
	if(dydx>slopemax)
		{
		dy=((ep[m1+1]-ep[m1])-slopemax*(gx[m1+1]-gx[m1]))/2.0;
		ep[m1]  +=dy;
		ep[m1+1]-=dy;
/*
dydx=(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]);
printf("AAA %ld %e %e",m1,slopemax,dydx);getchar();
*/
		}
	if(dydx<-slopemax)
		{
		dy=((ep[m1+1]-ep[m1])+slopemax*(gx[m1+1]-gx[m1]))/2.0;
		ep[m1]  +=dy;
		ep[m1+1]-=dy;
/*
dydx=(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]);
printf("BBB %ld %e %e",m1,slopemax,dydx);getchar();
*/
		}
	}
/**/
/**/
/**/
/* Add Erosion step */
ertimesum+=timestep;
/**/
/**/
/**/
/* Print Results */
if (printmod)
	{
	printf("\n EROSION STEP = %e YEARS    EROSION TIME = %e YEARS \n",timestep/3.15576e+7,ertimesum/3.15576e+7);
/*
getchar();
*/
	}
}
while(ertimesum<ertimesum0);
/* Restore timestep */
timestep=ertimesum0;
}
/* Errosion Surface progress */





/* Hydration front progress after H2O budget */
double hydration2()
{
/* Val buffer */
double ysurf,vfiltr,yfiltr,dydx,dydx1,sy1,sy2,sy3,sy4,sy5,e1,mwamin,x0,y0,x1,y1,vx1,vy1;
double hytimesum,hytimesum0;
/* TD Database variables */
double W0,W1,W2,W3,R0,R1,R2,R3,n,e,dx,dy;
double mtk,mpb,mwa,mro,dmwa,wro;
long int m1,m2,m3,mm1,marknum1=marknum;
int mm2,mm3,n1,n2;
/**/
printf("\n WATER Transport BEG \n");
/* Marker steps */
dx=dxwater;
dy=dywater;
/**/
/**/
/* Min water contents in the hydraten mantle wt% */
mwamin=0.1;
/* Min Distance from erosion surface for water release */
ysurf=8000.0;
/**/
/**/
/* Clear wa[] wt */
for (m1=0;m1<nodenum;m1++)
	{
	wa0[m1]=0;
	wa1[m1]=0;
	sol0[m1]=0;
	sol1[m1]=0;
	sol0[nodenum+m1]=1e+30;
	sol1[nodenum+m1]=-1e+30;
	sol0[nodenum2+m1]=1e+30;
	sol1[nodenum2+m1]=-1e+30;
	fre0[         m1]=1e+30;
	fre0[nodenum +m1]=-1e+30;
	fre0[nodenum2+m1]=1e+30;
	fre0[nodenum3+m1]=-1e+30;
	}
/**/
/**/
/**/
/* Fluid marker generation cycle */
for (mm1=0;mm1<marknum;mm1++)
{
/* Marker type */
mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
/**/
/* Marker cell number */
m1=m1serch((double)markx[mm1]);
m2=m2serch((double)marky[mm1]);
m3=m1*ynumy+m2;
/**/
/* Erosion surface */
e1=(markx[mm1]-gx[m1])/(gx[m1+1]-gx[m1]);
sy1=(e1*ep[m1+1]+(1.0-e1)*ep[m1]);
/**/
/*
printf("%ld %e %e %e %e %e %e",mm1,markx[mm1],marky[mm1],e1,sy1,sy2,sy3);getchar();
*/
/**/
/* Check markers out of grid and within hydration range */
if(markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && (markk[mm1]>0 || markt[mm1]>=50) && markt[mm1]<100)
if((double)(markd[mm1])>=0 && (double)(markw[mm1])>=0 && mm2>1 && mm2!=9 && mm2!=10)
	{
	if(mm2<50)
		{
		/* P, T parameters calc */
		allinterp((double)(markx[mm1]),(double)(marky[mm1]));
		mpb=eps[10]*1e-5;
		mtk=(double)(markk[mm1]);
		/**/
		/* Mantle to Antigorite transformation */
		antigor(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1);
/*
*/
		/**/
		/* Rocks to rock+melt transformation */
/*
		melting(mtk,mpb,mm1);
*/
		if (markt[mm1]>=20)
			{
			/* Check melting extent */
			if(fre0[        +m3]>markx[mm1]-dx) fre0[         m3]=markx[mm1]-dx;
			if(fre0[nodenum +m3]<markx[mm1]+dx) fre0[nodenum +m3]=markx[mm1]+dx;
			if(fre0[nodenum2+m3]>marky[mm1]-dy) fre0[nodenum2+m3]=marky[mm1]-dy;
			if(fre0[nodenum3+m3]<marky[mm1]+dy) fre0[nodenum3+m3]=marky[mm1]+dy;
			}
		/* Compute TD variables */
		tdbasecalc(mtk,mpb,mm2,mm1);
		mro=eps[41];
		mwa=eps[42];
		/**/
		/* Water changes in kg/m3 calc */
		dmwa=mro*(mwa-markw[mm1])*1e-2;
/*
{printf("H2O MARKER %ld %d %d %e %e   %e %e  %e %e   %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1],dmwa);getchar();}
{printf("H2O RELEASE %ld %d %d %e %e   %e %e  %e %e   %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1],dmwa);getchar();}
*/
		/**/
		/* Add water changes to the current cell, kg/m3 */
		/* Water release */
		if ((markw[mm1]-mwa)>dmwamin)
			{
			/* Save new water content */
			markw[mm1]=mwa;
			/* Generation of fluid marker (NO FLUID From melts */
			if (markt[mm1]<20 && marky[mm1]>sy1)
				{
				markt[marknum1]=markt[mm1]+50;
				markx[marknum1]=markx[mm1];
				marky[marknum1]=marky[mm1];
				markk[marknum1]=markk[mm1];
				markd[marknum1]=1050.0;
				markw[marknum1]=-dmwa;
				/* Add aditional markers counter */
				marknum1++;
				/* Check hydration extent */
				if(sol0[nodenum+m3]>markx[mm1]-dx) sol0[nodenum+m3]=markx[mm1]-dx;
				if(sol1[nodenum+m3]<markx[mm1]+dx) sol1[nodenum+m3]=markx[mm1]+dx;
				if(sol0[nodenum2+m3]>marky[mm1]-dy) sol0[nodenum2+m3]=marky[mm1]-dy;
				if(sol1[nodenum2+m3]<marky[mm1]+dy) sol1[nodenum2+m3]=marky[mm1]+dy;
				}
			}
		else
		/* Water consuming */
			{
			if(dmwa>0)
				{
				wa1[m3]+=dmwa;
				sol1[m3]+=1.0;
				}
			}
		}
	else
	/* Fluid marker count */
		{
		/* Check position */
		if(marky[mm1]>sy1)
			{
			/* Check hydration extent */
			if(sol0[nodenum+m3]>markx[mm1]-dx) sol0[nodenum+m3]=markx[mm1]-dx;
			if(sol1[nodenum+m3]<markx[mm1]+dx) sol1[nodenum+m3]=markx[mm1]+dx;
			if(sol0[nodenum2+m3]>marky[mm1]-dy) sol0[nodenum2+m3]=marky[mm1]-dy;
			if(sol1[nodenum2+m3]<marky[mm1]+dy) sol1[nodenum2+m3]=marky[mm1]+dy;
			}
		else
		/* Erase fluid marker */
			{
			markx[mm1]=-1.0;
			markk[mm1]=0;
			}
		}
	}
}
/**/
/**/
/*
for(m1=0;m1<xnumx-1;m1++)
for(m2=0;m2<ynumy-1;m2++)
{
m3=m1*ynumy+m2;
if(sol0[nodenum+m3]<xsize) {printf("%ld %ld  %e %e  %e %e  %e %e   %e %e %e %e",m1,m2,gx[m1],gx[m1+1],gy[m2],gy[m2+1],wa0[m3],wa1[m3],sol0[nodenum+m3],sol1[nodenum+m3],sol0[nodenum*2+m3],sol1[nodenum*2+m3]);getchar();}
}
*/
/**/
/* Rock hydration cycle */
for (mm1=0;mm1<marknum;mm1++)
if(markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && markt[mm1]<50)
{
/*
printf("%ld %e %e %e %e %e %e",mm1,markx[mm1],marky[mm1],e1,sy1,sy2,sy3);getchar();
*/
/**/
/* Marker cell number */
m1=m1serch((double)markx[mm1]);
m2=m2serch((double)marky[mm1]);
m3=m1*ynumy+m2;
/* Check markers within hydration range */
if(markx[mm1]>sol0[nodenum+m3] && marky[mm1]>sol0[nodenum2+m3] && (double)(markx[mm1])<sol1[nodenum+m3] && (double)(marky[mm1])<sol1[nodenum2+m3])
	{
	/* Fluid presence mark */
	markv[mm1]=-1.0;
if(markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==12 || markt[mm1]==14 || markt[mm1]==5 || markt[mm1]==6)
	{
	/* Mantle Hydration */
	if (markt[mm1]!=5 && markt[mm1]!=6)
		{
		mm2=markt[mm1]=11;
		}
	else
		{
		mm2=markt[mm1]=markt[mm1]+12;
		}
	/* P, T parameters calc */
	allinterp((double)(markx[mm1]),(double)(marky[mm1]));
	mpb=eps[10]*1e-5;
	mtk=(double)(markk[mm1]);
	/**/
	/* Mantle to Antigorite transformation */
	antigor(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1);
/*
*/
	/**/
	/* Rocks to rock+melt transformation */
/*
	melting(mtk,mpb,mm1);
*/
	if (markt[mm1]>=20)
		{
		/* Check melting extent */
		if(fre0[        +m3]>markx[mm1]-dx) fre0[         m3]=markx[mm1]-dx;
		if(fre0[nodenum +m3]<markx[mm1]+dx) fre0[nodenum +m3]=markx[mm1]+dx;
		if(fre0[nodenum2+m3]>marky[mm1]-dy) fre0[nodenum2+m3]=marky[mm1]-dy;
		if(fre0[nodenum3+m3]<marky[mm1]+dy) fre0[nodenum3+m3]=marky[mm1]+dy;
		}
	/**/
	/* Thermodynamic database use for Ro, Water */
	/* Compute TD variables */
	tdbasecalc(mtk,mpb,mm2,mm1);
	mro=eps[41];
	mwa=eps[42];
	/**/
	/* Water changes in kg/m3 calc */
	dmwa=mro*(mwa-markw[mm1])*1e-2;
/*
{printf("H2O HYDRATE %ld %d %d %e %e   %e %e  %e %e   %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1],dmwa);getchar();}
*/
	/**/
	/* Add water changes to the current cell, kg/m3 */
	/* Water consuming */
	if (dmwa>0)
		{
		wa1[m3]+=dmwa;
		sol1[m3]+=1.0;
		}
	}
	}
}
/**/
/**/
/*
for(m1=0;m1<xnumx-1;m1++)
for(m2=0;m2<ynumy-1;m2++)
{
m3=m1*ynumy+m2;
printf("%ld %ld  %e %e  %e %e  %e %e   %e %e %e %e",m1,m2,gx[m1],gx[m1+1],gy[m2],gy[m2+1],wa0[m3],wa1[m3],sol0[nodenum+m3],sol1[nodenum+m3],sol0[nodenum*2+m3],sol1[nodenum*2+m3]);getchar();
}
*/
/**/
/**/
/**/
/* Fluid marker computing cycle */
for (mm1=0;mm1<marknum1;mm1++)
{
/**/
/*
printf("%ld %e %e %e %e %e %e",mm1,markx[mm1],marky[mm1],e1,sy1,sy2,sy3);getchar();
*/
/**/
/* Check markers out of grid and within hydration range */
if(markt[mm1]>=50 && markt[mm1]<100 && markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize)
	{
	/* Marker cell number */
	m1=m1serch((double)markx[mm1]);
	m2=m2serch((double)marky[mm1]);
	m3=m1*ynumy+m2;
	/**/
	/* Erosion surface */
	e1=(markx[mm1]-gx[m1])/(gx[m1+1]-gx[m1]);
	sy1=(e1*ep[m1+1]+(1.0-e1)*ep[m1]);
	/* Water in melt region conversion */
	if(markd[mm1]<1100.0 && markx[mm1]>fre0[m3] && marky[mm1]>fre0[nodenum2+m3] && markx[mm1]<fre0[nodenum+m3] && marky[mm1]<fre0[nodenum3+m3]) markd[mm1]=1150.0;
/*
if(markd[mm1]>1100.0) {printf("H2O REDUCE1 %ld %d   %ld %ld %ld   %e %e    %e %e %e",mm1,mm2,m1,m2,m3,wa0[m3],wa1[m3],markx[mm1],markw[mm1],markd[mm1]);getchar();}
*/
	/* Check position, no fluid above erosion/sedimentation level, no fluid passing through the melt */
	if(marky[mm1]>sy1 && marky[mm1]<zdeep && (markd[mm1]<1100.0 || (markx[mm1]>fre0[m3] && marky[mm1]>fre0[nodenum2+m3] && markx[mm1]<fre0[nodenum+m3] && marky[mm1]<fre0[nodenum3+m3])))
		{
		wa0[m3]+=markw[mm1];
		sol0[m3]+=1.0;
		}
	else
	/* Erase fluid marker */
		{
		markx[mm1]=-1.0;
		markk[mm1]=0;
		}
/*
if(markd[mm1]>1100.0) {printf("H2O REDUCE2 %ld %d   %ld %ld %ld   %e %e    %e %e %e",mm1,mm2,m1,m2,m3,wa0[m3],wa1[m3],markx[mm1],markw[mm1],markd[mm1]);getchar();}
*/
	}
}
/**/
/**/
/**/
/*
for(m1=0;m1<xnumx-1;m1++)
for(m2=0;m2<ynumy-1;m2++)
{
m3=m1*ynumy+m2;
printf("%ld %ld  %e %e  %e %e  %e %e   %e %e %e %e",m1,m2,gx[m1],gx[m1+1],gy[m2],gy[m2+1],wa0[m3],wa1[m3],sol0[nodenum+m3],sol1[nodenum+m3],sol0[nodenum*2+m3],sol1[nodenum*2+m3]);getchar();
}
*/
/* Fluid marker consuming cycle */
for (mm1=0;mm1<marknum1;mm1++)
{
/* Marker type */
mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
/**/
/* Marker cell number */
m1=m1serch((double)markx[mm1]);
m2=m2serch((double)marky[mm1]);
m3=m1*ynumy+m2;
/**/
/*
printf("%ld %e %e %e %e %e %e",mm1,markx[mm1],marky[mm1],e1,sy1,sy2,sy3);getchar();
if(mm1>marknum){printf("%ld %d  %e %e  %e %e ",mm1,mm2,markx[mm1],marky[mm1],e1,sy1);getchar();}
*/
/**/
/* Change water consuming rocks  and fluid makers */
if(markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && (markk[mm1]>0 || markt[mm1]>=50) && markt[mm1]<100)
if((double)(markd[mm1])>=0 && (double)(markw[mm1])>=0 && mm2>1 && mm2!=9 && mm2!=10 && mm2!=12 && mm2!=14 && mm2!=5 && mm2!=6)
	{
	if(mm2<50)
		{
		/* P, T parameters calc */
		allinterp((double)(markx[mm1]),(double)(marky[mm1]));
		mpb=eps[10]*1e-5;
		mtk=(double)(markk[mm1]);
		/**/
		/* Thermodynamic database use for Ro, Water */
		/* Compute TD variables */
		tdbasecalc(mtk,mpb,mm2,mm1);
		mwa=eps[42];
		/**/
		/* Water changes in kg/m3 calc */
		dmwa=mwa-markw[mm1];
/*
{printf("TD! %ld %d %d %e %e %e %e ",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro);getchar();}
{printf("TDa %d %d %e %e %e %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,hpor,ppor,tpor,mwa);getchar();}
{printf("TDb %d %d %e %e %e %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,hpor,ppor,tpor,mwa);getchar();}
*/
		/**/
		/* Add water changes to the current cell, kg/m3 */
		/* Water consuming */
		if(dmwa>0)
			{
			if (wa1[m3]<=wa0[m3])
				{
				/* Save complete new water content */
				markw[mm1]=mwa;
				}
			else
				{
				/* COmpute, Save partial new water content */
				markw[mm1]=markw[mm1]+dmwa*wa0[m3]/wa1[m3];
				}
/*
{printf("H2O CONSUME %ld %d %d %e %e   %e %e  %e %e   %e    %e %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1],dmwa,wa0[m3],wa1[m3]);getchar();}
*/
			}
		}
	else
	/* Fluid marker change */
		{
		/* Check position, no fluid above erosion/sedimentation level, no fluid passing through the melt */
		if(wa1[m3]<wa0[m3])
			{
			/* Count water changes for fluid marker */
			markw[mm1]*=1.0-wa1[m3]/wa0[m3];
			}
		else
		/* Erase fluid marker */
			{
			markx[mm1]=-1.0;
			markk[mm1]=0;
			}
/*
if(markd[mm1]>1100.0) {printf("H2O REDUCE2 %ld %d   %ld %ld %ld   %e %e    %e %e %e",mm1,mm2,m1,m2,m3,wa0[m3],wa1[m3],markx[mm1],markw[mm1],markd[mm1]);getchar();}
*/
		}
	}
}
/*
marknum=marknum1;
return 0;
*/
/**/
/**/
/**/
/* Reset aditional markers */
printf("\n WATER BEG Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1);
mm1=0;
while(marknum1>marknum && mm1<marknum)
	{
	/* Reload marker */
	if((markx[mm1]<0 || marky[mm1]<0 || (double)(markx[mm1])>xsize || (double)(marky[mm1])>ysize) && markt[mm1]<100) 
		{
		/* Decrease aditional markers counter */
		marknum1--;
		if(markx[marknum1]>=0);
			{
			/* Type save */
			markt[mm1]=markt[marknum1];
			/* X,Y, water reload */
			markx[mm1]=markx[marknum1];
			marky[mm1]=marky[marknum1];
			markw[mm1]=markw[marknum1];
			markd[mm1]=markd[marknum1];
			markk[mm1]=markk[marknum1];
			}
		}
	/* Increase markers counter */
	mm1++;
	}
printf("\n WATER END Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1);
/* Set new marker number */
marknum=marknum1;
/**/
/**/
/**/
return 0;
}
/* Hydration front progress after H2O budget */


/* Errosion Surface progress */
void erosion()
{
/* Val buffer */
double v0,v1,dydx,x1,vx1,vy1,dy;
double ertimesum,ertimesum0;
long int m1,m2;
/**/
/* Erosion Solution Cycle ------------------------------------------ */
ertimesum=0;
ertimesum0=timestep;
do
{
/* Save old cycle results */
for (m1=0;m1<xnumx;m1++)
	{
	ep0[m1]=ep[m1];
	ep0[xnumx+m1]=ep[xnumx+m1];
	}
/**/
/**/
/**/
/* Initial timestep definition */
timestep=ertimesum0-ertimesum;
/**/
/**/
/**/
/* Erosion timestep definition using material velosity field */
for (m1=0;m1<xnumx;m1++)
	{
	/* Calc horisontal Coordinate */
	x1=gx[m1];
	/**/
	/* EROSION SURFACE */
	/* Calc material velocity on the Surface using velosity field */
	allinteri(x1,ep0[m1]);
	vx1=eps[11];
	vy1=eps[12];
	/* Check horizontal timestep */
	/* Calc x derivative of y position of the Surface using upwind differences */
	dydx=0;
	if(vx1>0 && m1>0)
		{
		timestep=MINV(timestep,(gx[m1]-gx[m1-1])/vx1);
/*
printf("111 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1-1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
		}
	if(vx1<0 && m1<xnumx-1)
		{
		timestep=MINV(timestep,(gx[m1]-gx[m1+1])/vx1);
/*
printf("222 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1+1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
		}
	/* Check vertical timestep */
	if(vy1)
		{
		/* Horizontal line num definition */
		m2=m2serch(ep0[m1]);
		/* Check timestep */
		timestep=MINV(timestep,(gy[m2+1]-gy[m2])/ABSV(vy1));
/*
printf("333 %ld  %e %e %e %e %e %e %e",m2,vx1,vy1,(gy[m2+1]-gy[m2])/ABSV(vy1),ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
		}
	/**/
	/**/
	/* INITIAL SURFACE */
	/* Calc material velocity on the Initial Surface using velosity field */
	allinteri(x1,ep0[xnumx+m1]);
	vx1=eps[11];
	vy1=eps[12];
	/* Check horizontal timestep */
	/* Calc x derivative of y position of the Surface using upwind differences */
	dydx=0;
	if(vx1>0 && m1>0)
		{
		timestep=MINV(timestep,(gx[m1]-gx[m1-1])/vx1);
/*
printf("444 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1-1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
		}
	if(vx1<0 && m1<xnumx-1)
		{
		timestep=MINV(timestep,(gx[m1]-gx[m1+1])/vx1);
/*
printf("555 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1+1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
		}
	/* Check vertical timestep */
	if(vy1)
		{
		/* Horizontal line num definition */
		m2=m2serch(ep0[xnumx+m1]);
		/* Check timestep */
		timestep=MINV(timestep,(gy[m2+1]-gy[m2])/ABSV(vy1));
/*
printf("666 %ld  %e %e %e %e %e %e %e",m2,vx1,vy1,(gy[m2+1]-gy[m2])/ABSV(vy1),ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
		}
	}
/*
printf("777 %e %e %e %e",ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
/**/
/**/
/**/
/* Displace Surface boundary */
/*
for (m1=1;m1<xnumx-1;m1++)
*/
for (m1=0;m1<xnumx;m1++)
	{
	/* EROSION SURFACE */
	/* Calculation of errosion rate */
	v0=0;
	if(ep0[m1]<eroslev)
		{
		v0=eroscon+eroskoe*(eroslev-ep0[m1]);
		}
	/* Calculation of sedimentation rate */
	v1=0;
	if(ep0[m1]>sedilev)
		{
		v1=sedicon+sedikoe*(ep0[m1]-sedilev);
		}
	/* Calc horisontal Coordinate */
	x1=gx[m1];
	/**/
	/* Calc material velocity on the Surface using velosity field */
	allinteri(x1,ep0[m1]);
	vx1=eps[11];
	vy1=eps[12];
	/**/
	/* Erase erosion/sedimentation rate for marginal points */
	if((m1==0 && vx1>0) || (m1==xnumx-1 && vx1<0)) v0=v1=0;
	/**/
	/* Calc x derivative of y position of the Surface using upwind differences */
	dydx=0;
	if(vx1>0 && m1>0)
		{
		dydx=(ep0[m1]-ep0[m1-1])/(gx[m1]-gx[m1-1]);
/*
printf("AAA %e %e",ep0[m1],dydx);getchar();
*/
		}
	if(vx1<0 && m1<xnumx-1)
		{
		dydx=(ep0[m1+1]-ep0[m1])/(gx[m1+1]-gx[m1]);
/*
printf("BBB %e %e",ep0[m1],dydx);getchar();
*/
		}
	/* Recalc new Surface position */
	ep[m1]+=timestep*(v0-v1+vy1-dydx*vx1);
/*
printf("SURFACE %ld %e %e %e %e %e %e %e %e",m1,x1,v0,v1,vx1,vy1,dydx,ep[m1]);getchar();
*/
	/**/
	/**/
	/**/
	/* INITIAL SURFACE */
	/* Initial surface displacement */
	/* Calc material velocity on the Surface using velosity field */
	allinteri(x1,ep0[xnumx+m1]);
	vx1=eps[11];
	vy1=eps[12];
	/* Calc x derivative of y position of Initial Surface using upwind differences */
	dydx=0;
	if(vx1>0 && m1>0)
		{
		dydx=(ep0[xnumx+m1]-ep0[xnumx+m1-1])/(gx[m1]-gx[m1-1]);
/*
printf("AAA %e ",dydx);getchar();
printf("AAA %e ",dydx);getchar();
*/
		}
	if(vx1<0 && m1<xnumx-1)
		{
		dydx=(ep0[xnumx+m1+1]-ep0[xnumx+m1])/(gx[m1+1]-gx[m1]);
/*
printf("BBB %e ",dydx);getchar();
*/
		}
	/* Recalc new Initial Surface position */
	ep[xnumx+m1]+=timestep*(vy1-dydx*vx1);
	/**/
	}
/**/
/**/
/**/
/**/
/* Relax EROSION surface */
if (0==0)
for (m1=0;m1<xnumx-1;m1++)
	{
	/* Calc x derivative of y position */
	dydx=(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]);
	/* Relax surface for critical slope */
	if(dydx>slopemax)
		{
		dy=((ep[m1+1]-ep[m1])-slopemax*(gx[m1+1]-gx[m1]))/2.0;
		ep[m1]  +=dy;
		ep[m1+1]-=dy;
/*
dydx=(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]);
printf("AAA %ld %e %e",m1,slopemax,dydx);getchar();
*/
		}
	if(dydx<-slopemax)
		{
		dy=((ep[m1+1]-ep[m1])+slopemax*(gx[m1+1]-gx[m1]))/2.0;
		ep[m1]  +=dy;
		ep[m1+1]-=dy;
/*
dydx=(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]);
printf("BBB %ld %e %e",m1,slopemax,dydx);getchar();
*/
		}
	}
	/**/
	/**/
	/**/
/* Add Erosion step */
ertimesum+=timestep;
/**/
/**/
/**/
/* Print Results */
if (printmod)
	{
	printf("\n EROSION STEP = %e YEARS    EROSION TIME = %e YEARS \n",timestep/3.15576e+7,ertimesum/3.15576e+7);
/*
getchar();
*/
	}
}
while(ertimesum<ertimesum0);
/* Restore timestep */
timestep=ertimesum0;
}
/* Errosion Surface progress */



/* Thermodynamic database use for ro, Cp */
void tdbasecalc(double mtk, double mpb, int mm2, long int mm1)
{
/* TD Database variables,  dTK,dPB - TK, PB step for tabulation in TD database */
double H0,H1,H2,H3,R0,R1,R2,R3,G0,G1,G2,G3,W0,W1,W2,W3,n,e;
/* Val Buffers */
int n1,n2,mm3,ynpb;
double mhh0,mhh1,mdhh,maa,mwa,dmwa,wro,mro,mcp,mbb,mgg,mkt,mkt1,pbmax,xold,kr01,kr1,kr10,xkr,krad;
long int m1=wn[0];
double sy1,e1;
/**/
/* Maximal pressure for the shallow database */
pbmax=pbmin+pbstp*(double)(pbnum-1);
/* Adiabate computing */
ynpb=0; if(1==0 && timesum<3.15576e+7*1e+3) {mpb*=timesum/(3.15576e+7*1e+3); ynpb=1;}
/**/
/* Reset TD variables */
eps[40]=eps[41]=eps[42]=eps[43]=eps[44]=eps[45]=0;
/**/
/* Thermal conductivity */
/* m895 Dry peridotite Fe=12 */
/* Olivine: Hoffmeister, 1999; Hoffmeister & Yuen, 2005 */
if(mpb<235000.0)
	{
	/* Lattice k */
	mkt1=(1.878+770.9/MINV(mtk,1200.0))*(1.0+4.26e-6*mpb);
	/* Radiative k 0.1 mm */
	kr01=pow(mtk/4000.0,3.0);
	/* Radiative k 1 mm */
	kr1=pow(mtk/1774.0,3.0);
	/* Radiative k 10 mm */
	xkr=pow(mtk/1636.0,10.0);
	xkr/=xkr+1.0; kr10=pow((mtk-1000.0*xkr)/1011.0,3.0)-0.7713*xkr;
	}
/* Perovskite: Hoffmeister, 1999; Hoffmeister & Yuen, 2005 */
else
	{
	/* Lattice k */
	mkt1=(1.291+1157.0/MINV(mtk,2100.0))*(1.0+2.50e-6*mpb);
	/* Radiative k 0.1 mm */
	kr01=pow(mtk/3591.0,3.0);
	/* Radiative k 1 mm */
	kr1=pow(mtk/2117.0,3.0);
	/* Radiative k 10 mm */
	xkr=pow(mtk/1500.0,4.0); xkr/=xkr+1.0;
	kr10=pow((mtk+4000.0*xkr)/5776.0,3.0)+2.822*xkr;
	}
krad=kr1;
/**/
/* Shallow TD base type */
if(mpb<pbmax && ynpb==0)
{
/* TD base type */
switch (mm2)
	{
	/* Dry Upper crust */
	case 5: mm3=11; break;
	/* Wet Upper crust */
	case 17: mm3=12; break;
	/* Dry Lower crust */
	case 6: mm3=13; break;
	/* Wet Lower crust */
	case 18: mm3=14; break;
	/* Sediments */
	case 2:
	case 3:
	case 4: mm3=5; break;
	/* Molten Sediments */
	case 37:
	case 25:
	case 22:
	case 23:
	case 24: mm3=6; break;
	/* Basalt */
	case 16:
	case 7: mm3=7; break;
	/* Molten Basalt */
	case 36:
	case 27: mm3=8; break;
	/* Gabbro */
	case 38:
	case 26:
	case 8: mm3=3; break;
	/* Molten Gabbro */
	case 28: mm3=4; break;
	/* Dry peridotite */
	case 9:
	case 12:
	case 14:
	case 10: mm3=0; break;
	/* Wet peridotite */
	case 13:
	case 11: mm3=1; break;
	/* Molten peridotite */
	case 34: mm3=2; break;
	/* Unknown type */
	default: {printf("Unknown rock type for TD database %d",mm2); exit(0);}
	}
/* ABCD-4Cell Number */
e=(mtk-tkmin)/tkstp;
if(e<0) e=0;
if(e>(double)(tknum-1)) e=(double)(tknum-1);
n=(mpb-pbmin)/pbstp;
if(n<0) n=0;
if(n>(double)(pbnum-1)) n=(double)(pbnum-1);
n1=(int)(e);
if(n1>tknum-2) n1=tknum-2;
n2=(int)(n);
if(n2>pbnum-2) n2=pbnum-2;
/* e,n Calc */
e=(e-(double)(n1));
n=(n-(double)(n2));
/* Ro H values */
/* 0 2 */
/* 1 3 */
R0=td[n1  ][n2  ][mm3][0]*1000.0;
R1=td[n1  ][n2+1][mm3][0]*1000.0;
R2=td[n1+1][n2  ][mm3][0]*1000.0;
R3=td[n1+1][n2+1][mm3][0]*1000.0;
H0=td[n1  ][n2  ][mm3][1]*1000.0*4.1837;
H1=td[n1  ][n2+1][mm3][1]*1000.0*4.1837;
H2=td[n1+1][n2  ][mm3][1]*1000.0*4.1837;
H3=td[n1+1][n2+1][mm3][1]*1000.0*4.1837;
W0=td[n1  ][n2  ][mm3][4];
W1=td[n1  ][n2+1][mm3][4];
W2=td[n1+1][n2  ][mm3][4];
W3=td[n1+1][n2+1][mm3][4];
G0=td[n1  ][n2  ][mm3][3]*1000.0;G0*=G0*R0;
G1=td[n1  ][n2+1][mm3][3]*1000.0;G1*=G1*R1;
G2=td[n1+1][n2  ][mm3][3]*1000.0;G2*=G2*R2;
G3=td[n1+1][n2+1][mm3][3]*1000.0;G3*=G3*R3;
/* Shear modulus calc by interpolation */
mgg=((G0*(1.0-n)+G1*n)*(1.0-e)+(G2*(1.0-n)+G3*n)*e);
/* Ro calc by interpolation */
mro=((R0*(1.0-n)+R1*n)*(1.0-e)+(R2*(1.0-n)+R3*n)*e);
/* Water wt% calc by interpolation */
mwa=((W0*(1.0-n)+W1*n)*(1.0-e)+(W2*(1.0-n)+W3*n)*e);
/* Add porocity fluid */
/* Erosion surface */
e1=(markx[mm1]-gx[m1])/(gx[m1+1]-gx[m1]);
sy1=marky[mm1]-(e1*ep[m1+1]+(1.0-e1)*ep[m1]);
if(marks0[mm2]>0 && sy1>0 && sy1<zmpor && mtk<tkpor) 
	{
	dmwa=marks0[mm2]*(tkpor-mtk)/(tkpor-273.15)*(zmpor-sy1)/zmpor;
	mwa+=dmwa;
	wro=1050.0;
	mro=mro/(1.0+dmwa*1e-2*(mro/wro-1.0));
/*
if(sy1>10000.0){printf("TD1 %d %d %e %e   %e %e   %e %e",mm2,mm3,mtk-273.15,sy1,mwa,mro,dmwa,zmpor);getchar();}
if(sy1>10000.0){printf("TD2 %d %d %e %e   %e %e   %e %e",mm2,mm3,mtk-273.15,sy1,mwa,mro,dmwa,zmpor);getchar();}
*/
	}
/* Cp calc by interpolation */
mcp=((H2-H0)*(1.0-n)+(H3-H1)*n)/tkstp;
if(mcp<1e+2) mcp=1e+2; else if(mcp>5e+4) mcp=5e+4;
/* Effective adiabatic betta=1/V*dV/dT=ro/T*[-dH/dP+V] calc by interpolation */
mbb=(2.0/(R1+R0)-(H1-H0)/pbstp/1e+5)*(1.0-e)+(2.0/(R3+R2)-(H3-H2)/pbstp/1e+5)*e;
mbb*=mro/mtk;
if(mbb<-1e-2) mbb=-1e-2; else if(mbb>1e-2) mbb=1e-2;
/* Effective compressibility term alpha=1/ro*d(ro)/dP calc by interpolation */
maa=(2.0/(R1+R0)*(R1-R0)*(1.0-e)+2.0/(R3+R2)*(R3-R2)*e)/pbstp/1e+5;
if(maa<0) maa=0;
/* Activation enthalpy recalc using enthalpy changes */
/* Current Enthalpy */
mhh1=((H0*(1.0-n)+H1*n)*(1.0-e)+(H2*(1.0-n)+H3*n)*e);
/* Pmin Enthalpy */
mhh0=(td[n1][0 ][mm3][1]*(1.0-e) + td[n1+1][0 ][mm3][1]*e)*1000.0*4.1837;
/* Enthalpy Difference calc */
mdhh=(mhh1-mhh0);
/*
{printf("TD1 %d %d %e %e   %e %e %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mgg,mro,mwa,mcp,mbb);getchar();}
{printf("TD1 %d %d %e %e   %e %e %e %e %e %e   %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mgg,mro,mwa,mcp,mbb,maa,mdhh,mhh1,mhh0);getchar();}
eps[47]=mhh1;
eps[48]=mhh0;
*/
/* Save TD variables */
eps[40]=mgg;
eps[41]=mro;
eps[42]=mwa;
eps[43]=mcp;
eps[44]=mbb;
eps[45]=maa;
eps[46]=mdhh;
eps[47]+=krad;
}
/**/
/**/
/* Deep TD base type */
if(1==0 || mpb>0.75*pbmax || ynpb==1)
{
switch (mm2)
	{
	/* MORB DATABASE */
	/* UPPER, LOWER Crust */
	case 5:
	case 6:
	case 17:
	case 18:
	case 25:
	case 26:
	case 37:
	case 38:
	/* Sediments */
	case 2:
	case 3:
	case 4:
	/* Molten Sediments */
	case 22:
	case 23:
	case 24:
	/* Basalt */
	case 16:
	case 7:
	/* Molten Basalt */
	case 36:
	case 27:
	/* Gabbro */
	case 8:
	/* Molten Gabbro */
	case 28: mm3=10; break;
	/**/
	/* PIROLITE DATABASE */
	/* Dry peridotite */
	case 9:
	case 12:
	case 14:
	case 10:
	/* Wet peridotite */
	case 13:
	case 11:
	/* Molten peridotite */
	case 34: mm3=9; break;
	/* Unknown type */
	default: {printf("Unknown rock type for TD database %d",mm2); exit(0);}
	}
/* ABCD-4Cell Number */
e=(mtk-tkmin1)/tkstp1;
if(e<0) e=0;
if(e>(double)(tknum1-1)) e=(double)(tknum1-1);
n=(mpb-pbmin1)/pbstp1;
if(n<0) n=0;
if(n>(double)(pbnum1-1)) n=(double)(pbnum1-1);
n1=(int)(e);
if(n1>tknum1-2) n1=tknum1-2;
n2=(int)(n);
if(n2>pbnum1-2) n2=pbnum1-2;
/* e,n Calc */
e=(e-(double)(n1));
n=(n-(double)(n2));
/* Ro H values */
/* 0 2 */
/* 1 3 */
R0=td[n1  ][n2  ][mm3][0]*1000.0;
R1=td[n1  ][n2+1][mm3][0]*1000.0;
R2=td[n1+1][n2  ][mm3][0]*1000.0;
R3=td[n1+1][n2+1][mm3][0]*1000.0;
H0=td[n1  ][n2  ][mm3][1]*1000.0*4.1837;
H1=td[n1  ][n2+1][mm3][1]*1000.0*4.1837;
H2=td[n1+1][n2  ][mm3][1]*1000.0*4.1837;
H3=td[n1+1][n2+1][mm3][1]*1000.0*4.1837;
W0=td[n1  ][n2  ][mm3][4];
W1=td[n1  ][n2+1][mm3][4];
W2=td[n1+1][n2  ][mm3][4];
W3=td[n1+1][n2+1][mm3][4];
G0=td[n1  ][n2  ][mm3][3]*1000.0;G0*=G0*R0;
G1=td[n1  ][n2+1][mm3][3]*1000.0;G1*=G1*R1;
G2=td[n1+1][n2  ][mm3][3]*1000.0;G2*=G2*R2;
G3=td[n1+1][n2+1][mm3][3]*1000.0;G3*=G3*R3;
/* Shear modulus calc by interpolation */
mgg=((G0*(1.0-n)+G1*n)*(1.0-e)+(G2*(1.0-n)+G3*n)*e);
/* Ro calc by interpolation */
mro=((R0*(1.0-n)+R1*n)*(1.0-e)+(R2*(1.0-n)+R3*n)*e);
/* Water wt% calc by interpolation */
mwa=0;
/* Water in crystals */
if(mm2!=9 && mm2!=10 && mm2!=14 && mpb<235000.0) 
	{
	dmwa=0.1;
	mwa+=dmwa;
	wro=1050.0;
	mro=100.0/((100.0-dmwa)/mro+dmwa/wro);
/*
{printf("TD1 %d %d %e %e   %e %e %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mgg,mro,mwa,mcp,mbb);getchar();}
*/
	}
/* Cp calc by interpolation */
mcp=((H2-H0)*(1.0-n)+(H3-H1)*n)/tkstp1;
if(mcp<1e+2) mcp=1e+2; else if(mcp>5e+4) mcp=5e+4;
/* Effective adiabatic betta=1/V*dV/dT=ro/T*[-dH/dP+V] calc by interpolation */
mbb=(2.0/(R1+R0)-(H1-H0)/pbstp1/1e+5)*(1.0-e)+(2.0/(R3+R2)-(H3-H2)/pbstp1/1e+5)*e;
mbb*=mro/mtk;
if(mbb<-1e-2) mbb=-1e-2; else if(mbb>1e-2) mbb=1e-2;
/* Effective compressibility term alpha=1/ro*d(ro)/dP calc by interpolation */
maa=(2.0/(R1+R0)*(R1-R0)*(1.0-e)+2.0/(R3+R2)*(R3-R2)*e)/pbstp1/1e+5;
if(maa<0) maa=0;
/* Activation enthalpy recalc using enthalpy changes */
/* Current Enthalpy */
mhh1=((H0*(1.0-n)+H1*n)*(1.0-e)+(H2*(1.0-n)+H3*n)*e);
/* Pmin Enthalpy */
mhh0=(td[n1][0 ][mm3][1]*(1.0-e) + td[n1+1][0 ][mm3][1]*e)*1000.0*4.1837;
/* Enthalpy Difference calc */
mdhh=(mhh1-mhh0);
/* Thermal conductivity */
mkt=mkt1+krad;
/**/
/*
if(timesum>0){printf("TD4  %d %d   %e %e     %e   %e",mm2,mm3,mtk-273.15,mpb/1000.0,mkt,pbmax);getchar();}
if(timesum>0){printf("TD4 %d %d %e %e %e  %e %e %e %e %e %e %e   %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,pbmax,mgg,mro,mwa,mcp,mbb,maa,mkt,mdhh,mhh1,mhh0);getchar();}
*/
/* Computing transitional parameters */
if(1==0 || mpb>pbmax || ynpb==1)
	{
	/* Save TD variables */
	eps[40]=mgg;
	eps[41]=mro;
	eps[42]=mwa;
	eps[43]=mcp;
	eps[44]=mbb;
	eps[45]=maa;
	eps[46]=mdhh;
	eps[47]=mkt;
/*
if(timesum>0){printf("TD2 %d %d %e %e %e  %e %e %e %e %e %e %e   %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,pbmax,mgg,mro,mwa,mcp,mbb,maa,mkt,mdhh,mhh1,mhh0);getchar();}
*/
	}
else
	{
	xold=(pbmax-mpb)/(0.25*pbmax);
	/* Save TD variables */
	mgg=mgg*(1.0-xold)+eps[40]*xold;
	mro=mro*(1.0-xold)+eps[41]*xold;
	mwa=mwa*(1.0-xold)+eps[42]*xold;
	mcp=mcp*(1.0-xold)+eps[43]*xold;
	mbb=mbb*(1.0-xold)+eps[44]*xold;
	maa=maa*(1.0-xold)+eps[45]*xold;
	mdhh=mdhh*(1.0-xold)+eps[46]*xold;
	mkt=mkt*(1.0-xold)+eps[47]*xold;
	eps[40]=mgg;
	eps[41]=mro;
	eps[42]=mwa;
	eps[43]=mcp;
	eps[44]=mbb;
	eps[45]=maa;
	eps[46]=mdhh;
	eps[47]=mkt;
/*
if(timesum>0){printf("TD3 %d %d %e %e %e  %e %e %e %e %e %e %e   %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,pbmax,mgg,mro,mwa,mcp,mbb,maa,mkt,mdhh,mhh1,mhh0);getchar();}
*/
	}
}
/*
{printf("TD2 %d %d %e %e   %e %e %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mgg,mro,mwa,mcp,mbb);getchar();}
*/
}
/* Thermodynamic database use for ro, Cp */




/* Check plastic marker yielding */
void plastic()
{
/* Counters buffers */
long int mm1;
int mm2;
double mpb,mtk,mnu;
for (mm1=0;mm1<marknum;mm1+=gridmod)
	{
	/* Marker type */
	mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
	/* Check markers out of grid */
	/**/
	if(markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && markk[mm1]>0 && markt[mm1]<50)
		{
		/**/
		/* P, T parameters calc */
		allintere((double)(markx[mm1]),(double)(marky[mm1]));
		mpb=eps[10]*1e-5;
		mtk=(double)(markk[mm1]);
		/**/
		/* Thermodynamic database use for ro, Cp */
		if (densimod==2)
		if(mm2>1)
			{
			/* Compute TD variables */
			tdbasecalc(mtk,mpb,mm2,mm1);
			}
		/**/
		/* Marker Rheology */
		mnu=viscalc(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1,mm2,1);
		}
	}
}
/* Check plastic marker yielding */
