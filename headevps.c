/* MACROS --------------------------------------------- */
/* Min, Max, Abs */
#define MINV(a,b) ((a)<=(b)? (a):(b))
#define MAXV(a,b) ((a)>=(b)? (a):(b))
#define ABSV(a)   ((a)>=0? (a):-(a))
/**/
/* Main Sizes */
#define MAXMAT  20000000 /* General Matrix Size */
#define MAXVAL  20000000 /* General Matrix Size */
#define MAXXLN  2201  /* Max vertical line num in Grid */
#define MAXYLN  301 /* Max horizontal line num in Grid */
#define MAXXMR   5 /* Max marker in cell for X direction */
#define MAXYMR   5 /* Max marker in cell for Y direction */
#define MAXTMR   200 /* Max markers types */
#define MAXPOS   200 /* Max Pos on Line num for wi[], wn[] buffers */
#define MAXFLN   3000 /* Max Output file names Num */
#define MAXNOD    MAXXLN*MAXYLN  /* Max Total nodes num for all grid */
#define MAXCEL    (MAXXLN-1)*(MAXYLN-1)  /* Max Total nodes num for all grid */
#define MAXPAR    MAXXLN*MAXYLN*4   /* Max Total par num for all grid */
#define MAXMRK   50000000 /* Max total marker Num */
#define MAXBON    MAXPAR /* Max Total bondary condition Equations for all grid */
/* End MACROS --------------------------------------------- */



/* ARRAYS --------------------------------------------- */
/* Processing+Service Arrays for solutions */
/* val0[] - matrix contents */
/* fre0[] - free member for lines */
/* bufv[] - buffer for matrix organisation */
/* lin0[] - line numbers for matrix contents */
/* num0[] - pos numbers for line in  val0[], pos0[] */
/* pos0[] - first pos numbers for line in  val0[], pos0[] */
/* sol0[],sol1[] - solution buffers */
/* wn[],wi[] - Line Num,Koef buffer */
double sol0[MAXPAR],sol1[MAXPAR];
double val0[MAXMAT],fre0[MAXPAR],bufv[MAXPAR];
int lin0[MAXMAT],num0[MAXPAR],cur0[MAXPAR];
int pos0[MAXPAR];
int bufn[MAXPOS]; 
int wn[MAXPOS];
double wi[MAXPOS];
/**/
/**/
/**/
/* PARDISO SOLVER ARRAYS */
/* ia[] - first coloumn of each row */
/* ja[] - coloumn position of each koef in the respective row from left to right */
/* a[]  - value of each koef */
/* b[] - Right hand side */
/* x[] - solutions */
int ia[MAXPAR];             
int ja[MAXVAL];                  
double b[MAXPAR];           
double x[MAXPAR];           
double a[MAXVAL];                 
/**/
/**/
/**/
/* Nodes information */
/* gx[], gy[] - coordinates of gridlines */
/* nd[] - Viscosity for SIGxx, SIGyy in Stokes equation, Pa*sek */
/* nu[] - Viscosity for SIGxy in Stokes equation, Pa*sek */
/* mu[] - Standard viskosity for node */
/* ep[] - Surface trace */
/* et[] - adiabatic modulus */
/* gg[], gd[] - shear modulus */
/* ro[] - densyty, kg/m^3 */
/* tk[] - themperature K */
/* cp[] - Heat capacity, J/kg */
/* kt[] - Thermal conductivity koef, Wt/m/K */
/* ht[] - Heat Sources, Wt/m3 */
/* vx[],vy[] - Cur circle Vx,Vy,  m/sek */
/* tk0[] - Last circle TK */
/* wa0[], wa1[] - Water content */
/* td[] - Thermodynamic data base */
/* bondm[] - bondary position in bondn[],bondv[] for cur par (0) - no bondary */
/* bondn[] - PAR1+1,PAR2+1,PAR3+1 num in bond equat: CURPAR=CONST+KOEF1*PAR1 */
/* bondv[] - CONST,KOEF1,KOEF2,KOEF3 val in bond equat: CURPAR=CONST+KOEF1*PAR1 */
/* pr[] - pressure, Pa */
/* pr0[] - last cycle pressure, Pa */
/* pkf[] - koef for pressure aproximation */
/* exx[], exy[] - strain rates for cells (exx[]) and nodes (exy[]) */
/* esp[] - rotation rate for nodes */
/* sxx[], sxy[] - stresses for cells (sxx[]) and nodes (sxy[]) */
/* sxxe[], sppe[], sxye[] - old elastic stresses, pressures for cells (sxxe[], sppe[]) and nodes (sxye[]) */
/* exxe[], exye[] - old strain rates from markers cells (exxe[] and nodes (exye[]) */
/* dro[] - Dro/Dt density changes with time */
/* drp[] - Dro/DP density changes */
double gx[MAXXLN],gy[MAXYLN];
double nu[MAXNOD],nd[MAXNOD],sxxe[MAXNOD],sppe[MAXNOD],sxye[MAXNOD],mu[MAXNOD],ep[MAXNOD],et[MAXNOD],gg[MAXNOD],gd[MAXNOD],ro[MAXNOD];
double nu0[MAXNOD],nd0[MAXNOD],sxxe0[MAXNOD],sppe0[MAXNOD],sxye0[MAXNOD],exxe[MAXNOD],exxe0[MAXNOD],exye[MAXNOD],exye0[MAXNOD],ep0[MAXNOD],et0[MAXNOD],gg0[MAXNOD],gd0[MAXNOD],ro0[MAXNOD];
double tk[MAXNOD],cp[MAXNOD],kt[MAXNOD],ht[MAXNOD];
double cp0[MAXNOD],kt0[MAXNOD],ht0[MAXNOD];
double vx[MAXNOD],vy[MAXNOD];
double tk0[MAXNOD], tk1[MAXNOD], tk2[MAXNOD], tk3[MAXNOD];
double wa0[MAXNOD], wa1[MAXNOD];
double bondv[MAXBON][4];
double pr[MAXNOD],pkf[20];
long int bondm[MAXPAR],bondn[MAXBON][3];
double td[360][360][15][5];
double exx[MAXNOD],dro[MAXNOD],exy[MAXNOD],esp[MAXNOD];
double sxx[MAXNOD],drp[MAXNOD],sxy[MAXNOD];
double sxx0[MAXNOD],sxy0[MAXNOD];
double dro0[MAXNOD],drp0[MAXNOD];
/**/
/**/
/**/
/* Markers and Rock information */
/* Markers information */
/* markx[], marky[] - X,Y of markers */
/* markk[] - temperature in K for marker */
/* markd[] - marker density, kg/m3 */
/* marke[] - finite strain for marker */
/* markp[] - pressure for marker */
/* markexx[], markexy[] - strain rate for marker */
/* markxx[], markxy[] - deviatoric stress components for marker */
/* markv[] - marker advected viscosity */
/* markt[] - rock type of markers */
/* Information for different rock types */
/* markim[] -  Immobility of marker type  Y(1)/No(0) */
/* marknu[], markdh[], markdv[], markss[] markmm[] -  Koef in ductile rheology Eq */
/* markll[], marka0[], marka1[], markb0[] markb1[], marke0[], marke1[], markf0[], markf1[] - Koef in brittle rheology Eq */
/* markgg[] - elastic shear modulus, Pa */
/* markn0[], markn1[], marks0[], marks1[] - viscosity and stress limits for individual rock types */
/* markro[], markaa[], markbb[], markcp[], markkt[], markkf[], markkv[], markht[] - ro,aro,bro,cp,kt, Tkt, Pkt, ht */
float markx[MAXMRK],marky[MAXMRK];
float markk[MAXMRK],markd[MAXMRK],marke[MAXMRK],markw[MAXMRK];
float markxx[MAXMRK],markv[MAXMRK],markxy[MAXMRK];
float markp[MAXMRK],markexx[MAXMRK],markexy[MAXMRK];
char markt[MAXMRK];
double marknu[MAXTMR],markdh[MAXTMR],markdv[MAXTMR],markss[MAXTMR],markmm[MAXTMR];
double markgg[MAXTMR],markll[MAXTMR],marka0[MAXTMR],marka1[MAXTMR],markb0[MAXTMR],markb1[MAXTMR],marke0[MAXTMR],marke1[MAXTMR],markf0[MAXTMR],markf1[MAXTMR];
double markro[MAXTMR],markaa[MAXTMR],markbb[MAXTMR],markcp[MAXTMR],markkt[MAXTMR],markht[MAXTMR],markkf[MAXTMR],markkp[MAXTMR];
double markn0[MAXTMR],markn1[MAXTMR],marks0[MAXTMR],marks1[MAXTMR];
int markim[MAXTMR];
/**/
/**/
/**/
/* Service, Connection Buffers */
/* vxy[] - Cur Vx,Vy */
/* eps[] - Cur Eps tenzors */
/* errbuf[] - Error buffer */
/* nunu[] - NUik buffer */
/* sa[] - input string buffer */
/* fl1in[] - input data file name */
/* fl1itp - input data file type */
/* fl1otp - output data file type */
/* fl0out[] - output data All file names */
/* fl0otp[] - output data files types */
/* fl0stp[],fl0cyc[] - output data file X,T,time steps cyc0max for Out files */
/* fl1out[] - output data Cur file name */
/* topo0out[] - output data ALL file names with only topography - ARMEL UPDATE */
/* topo1out[] - output data Cur file name with only topography - ARMEL UPDATE */
/* data0out[] - output data ALL file names with selected 2D data - ARMEL UPDATE */
/* data1out[] - output data Cur file name with selected 2D data - ARMEL UPDATE */
/* xn[], cn[] - high order FD cordinates and coefficients */
/* *fl - Streem for File Input/Output */
double vxy[20],errbuf[20],nunu[9][5],eps[100];
char sa[2500];
char fl1in[50],fl1out[50],topo1out[50],data1out[50];
int fl1itp,fl1otp;
char fl0out[MAXFLN][50],topo0out[MAXFLN][50],data0out[MAXFLN][50];
int fl0otp[MAXFLN];
double fl0stp[MAXFLN][10];
int fl0cyc[MAXFLN];
double xn[1000],cn[1000][10];
FILE *fl,*fl1;
/**/
/* End ARRAYS --------------------------------------------- */


/* VARIABLES -------------------------------------------------- */
/**/
/* Grid Parameters */
/* xnumx,ynumy - num of lines in grid in X,Y directions */
/* xnumx1,ynumy1 - num of cells in grid for X,Y directions */
/* mnumx,mnumy - num of markers in one cell for X,Y directions */
/* xsize,ysize - size of grid in X,Y directions, m */
/* GXKOEF,GYKOEF - Gravitation g for X,Y directions, m/sek^2 */
/* pinit - pressure on the upper boundary, Pa */
/* xstpx,ystpy - X,Y steps for grid */
/* kfx,kfy,kfxx,kfxy,kfyy,numdx,numdy - Koef for numeric differentiation */
/* cellnum - total Cell Num */
/* nodenum - total Node Num */
/* marknum - total Marker Num */
/* rocknum - rock types num */
/* bondnum - bondary condition equation num */
/* nodenum2 - nodenum*2 */
/* nodenum3 - nodenum*3 */
long int xnumx,ynumy,mnumx,mnumy;
long int xnumx1,ynumy1;
double xsize,ysize;
double GXKOEF,GYKOEF;
double pinit;
double xstpx,ystpy;
double kfx,kfy,kfxx,kfxy,kfyy,numdx,numdy,mardx,mardy;
long int cellnum,nodenum;
long int marknum,bondnum;
int rocknum=0;
long int nodenum2, nodenum3;
/**/
/**/
/**/
/* Service parameters */
/* printmod - print information on the monitor Y(1)/N(0) */
/* intermod - order of interpolation from nodes to markers 0-linear >=1-Non linear */
/* intermod1 - order of interpolation from markers to nodes 0-linear >=1-Non linear */
/* densimod - mode of  density calculation: 0-constant, 1-PT-dependent */
/* outgrid - marker move out of grid Y(0)/N(1) Orthogonal Only (2) */
/* fl0num - number of otput file Names */
/* pos0cur - Cur pos Counter in val0[] lin0[] */
long int pos0cur=0,printmod;
int fl0num,intermod,intermod1,outgrid=0,densimod;
/**/
/**/
/**/
/* Termodynamic database parameters <loadjxxx.c> */
/* coolzone,coolkoef - parameters for cooling zone of oceanic geotherm - ARMEL'S UPDATE */
double tkmin,pbmin,tkstp,pbstp;
int tknum,pbnum;
double tkmin1,pbmin1,tkstp1,pbstp1;
int tknum1,pbnum1;
double tkpor,zmpor,vyfluid,vymelt,dmwamin,tdeep,zdeep,dxwater,dywater,deserp,dyserp;
double lambfld,coolzone,coolkoef;
/**/
/**/
/**/
/* Errosion/Sedimentation parameters */
/* erosmod - errosion/sedimentation Y(1)/N(0)  */
/* eroslev - errosion level, m */
/* eroscon - constant erosion rate, m/s */
/* eroskoe - increment of erosion rate with elevation, 1/s */
/* sedilev - sedimentation level, m */
/* sedicon - constant sedimentation rate, m/s */
/* sedikoe - increment of sedimentation rate with burial, 1/s */
/* sedimnum - Num of cycles of sedimentation */
/* sedimcyc - Num of cycles of sedimentation for One Layer */
/* waterlev - water/air boundary, m */
/* basalty - basalt melting depth, m */
/* dehydrmin, dehydrmax - serpentine dehydration min, max depth, m */
int erosmod;
double eroslev,eroscon,eroskoe;
double sedilev,sedicon,sedikoe;
double waterlev,basalty,dehydrmin,dehydrmax,slopemax;
int sedimnum=0,sedimcyc=3;
/**/
/**/
/**/
/* Motion parameters */
/* cyc0max - Num of circles for cur calculation */
/* maxxystep - max Distance change for one time step, m */
/* xelvismin - min viscoelasticity factor */
/* maxtkstep - max Themperature for one time step, K */
/* maxtmstep - max Time for one time step, sek */
/* timebond - time limit for internal boundary conditions from start */
/* timestep - time step, sek */
/* timestepe - time step for viscoelasticity, sek */
/* timesum - time from start, sek */
/* timedir - direction of motion in time:1 - forward, -1 - backward */
/* movemod - solve Stoks+Continuity equations Y(1)/N(0) */
/* tempmod - solve Heat Transfer equation Y(1)/N(0) */
/* markmod - move markers Y(1)/N(0) */
/* ratemod - reset velosity and pressure Y(1)/N(0) */
/* gridmod - recalc grid parameters Y(1)/N(0) */
int cyc0max,movemod,tempmod,markmod,ratemod,gridmod;
double maxxystep,xelvismin,maxtkstep,maxtmstep,timestep,timestepe,timesum,timedir;
double plastmark,plastnew,plastnew1,plastold,plastmax,timebond;
/**/
/**/
/**/
/* General iteration parameters in <gaus.c> */
/* ckoef - Cur Koef for Zeydel method of iteration */
double ckoef;
/**/
/**/
/**/
/* V parameters in <move.c> */
/* DIVVMIN,STOKSMIN - Min Valid absolut Err val for Contin,Stokes Eq */
/* stoksmod - dNu/dT, dNu/dP Y(0)/N(1) in Stokes calc */
/* viscmod - Viscoelasticity Y(1,2)/N(0) */
/* stoksfd - Order of FD for EPSxx,EPSyy,EPSxy calculation */
/* nukoef - No function */
/* nubeg,nuend,nucontr - Min, Max, Max/Min limits of viscozity */
/* hidry - max depth with hidrostatic pressure of pore fluid */
/* hidrl - brittle weackening factor for hidrostatic pressure of pore fluid */
/* strmin, strmax - min max values of stress Pa allowed */
/* stredif - numerical stress relaxation coefficient */
double DIVVMIN,STOKSMIN;
int stoksmod,stoksfd;
double viscmod;
double nubeg,nuend,nucontr,hidry,hidrl,strmin,strmax,stredif;
/**/
/**/
/**/
/* T parameters in <heat.c> */
/* HEATMIN - Min Valid absolut Err val for Heat Equation */
/* heatmod - dK/dX account Y(1)/N(0) */
/* heatfd - Order of FD for Qx, Qy calculation */
/* heatdif - Numerical diffusion coefficient */
/* frictyn - Viscouse friction heat Y(1)/N(0) */
/* adiabyn - adiabatic heat calculation: N(0)/Y(1) */
double HEATMIN;
double heatdif;
int heatmod,heatfd,frictyn,adiabyn;
/**/
/* Standard Volume in <mark2d.c> */
double rostand=0;
/* End VARIABLES -------------------------------------------------- */



/* FUNCTONS PROTOTYPES ------------------------------------------- */
/**/
/* <load.c> FILE */
/* ffscanf() - Load single word from file fl without empty lines */
/* ffscanf1() - Load single word from file fl1 without empty lines */
/* loadconf() - Load configuration from mode.t3c */
/* loader() - Load information from data file */
/* saver() - Save information to data file */
/* gridcheck() - Calc,Check parameters of Grid */
void ffscanf();
void ffscanf1();
int loadconf();
void loader();
void saver(int, int);
void gridcheck();
/**/
/* <gaus.c> FILE */
/* gausmat3() - Solve system of linear equation by economic frontal Gauss method */
/* gausmat4() - Solve system of linear equation by Pardiso solver */
int gausmat3(int, long int, long int);
int gausmat4(int, long int, long int);
/**/
/* <move.c> FILE */
/* viterate() - General soubr for Vx,Vy,P Calc by Iterativ method */
/* xstokserr() -  Right part or  Err in X Stokes Equat for cur node calc */
/* ystokserr() -  Right part or  Err in Y Stokes Equat for cur node calc */
/* conterr() -  Right part or  Err in Contin Equat for cur node calc */
/* xbonderr() -  Right part or  Err in Boundary vX Equat for cur node calc */
/* ybonderr() -  Right part or  Err in Boundary vY Equat for cur node calc */
/* pbonderr() -  Right part or  Err in Boundary P Equat for cur cell calc */
/* maxvelstep() - Max time step for markers definition */
/* sxxcalc() etc. - Value or add EPS and SIG equations */
/* fdweight() - weight for FD calculation after Fornberg (1996) */
void viterate(int);
double xstokserr(long int, long int, int);
double ystokserr(long int, long int, int);
double conterr(long int, long int, int);
double xbonderr(long int, int);
double ybonderr(long int, int);
double pbonderr(long int, int);
void maxvelstep();
double sxxcalc(long int, long int, double);
double sxycalc(long int, long int, double);
void fdweight(int, int, double);
/**/
/* <heat.c> FILE */
/* titerate() -  Themperature recalc after time step */
/* heatserr() -  Right+Right part or  Err in Heat Equat for cur node calc */
/* tbonderr() -  Right part or  Err in Boundary T Equat for cur node calc */
/* tkrecalc() -  Calc average temperature after new marker position */
/* qxcalc(),qycalc() - Coefficients or value for Qx,Qy  Equations */ 
void titerate(int);
double heaterr(long int, long int, int);
double tbonderr(long int, int);
void tkrecalc();
double qxcalc(long int, long int, double);
double qycalc(long int, long int, double);
/**/
/* <mark.c> FILE */
/* movemark() - move markers by Runge-Kutta method */
/* erosmark() - Errosion/Sedimentation Function for markers */
/* erosion() - Errosion Surface Tracing */
/* allinter() - Vx,Vy,P,T,EPS,SIG calc for marker by interpolation */
/* allinterv() - Vx,Vy calc for marker by interpolation */
/* allinteri() - Vx,Vy, EPS calc for marker by interpolation */
/* allintert() - T calc for marker by interpolation */
/* allintere() - EPS calc for marker by interpolation */
/* allinters() - Vx,Vy,EPS*SIG calc for marker by interpolation */
/* allinterp() - P calc for marker by interpolation */
/* allinterd() - Sxx,Syy,Sxy calc for markers */
/* depthp() - P calc for marker as function of depth below the surface */
/* nodewt() - Weights for horisontal and vertical nodes calculation for marker interpolation */ 
/* ronurecalc() - recalc ro[],nu[] etc after new marker position */
/* dencalc() -  Calc density for given P,T after ro equation */
/* viscalc() -  Calc viscosity for given P,T after rheol equation */
/* antigor() - Mantle transformation to Antigorite */
/* basalt() - Basalt wet solidus depth along subduction plate definition */
/* hydration() - Hydration front progress recalc */
/* dehydration() - Serpentine dehydration min,max depth */
/* m1serch(), m2serch() - serch of nearest upper-left node of cell for  current location */
/* m1serch1(), m2serch1() - serch of nearest node for current location */
/* melting(), meltpart(), meltpart1() - melting of rocks account */
/* tdbasecalc() - computing TD variables from the database */
/* plastic() - check plastic yelding conditions */
void movemark();
void ronurecalc();
void allinter(double, double);
void allinterv(double, double);
void allintert(double, double);
void allintere(double, double);
void allinters(double, double);
void allinteri(double, double);
void allinterp(double, double);
void allinterd(double, double);
void depthp(double, double);
void nodewt(long int, long int, long int, long int, double, double, int, int);
double dencalc(double, double, double, double, int);
double viscalc(double, double, double, double, long int, int, int);
void antigor(double, double, double, double, long int);
double basalt();
double hydration();
double hydration2();
double dehydration(int);
void erosmark(long int, int);
void erosion();
long int m1serch(double);
long int m2serch(double);
double meltpart(double, double, double, double, long int, int);
void meltpart1(double, double, int);
void melting(double, double, long int);
void tdbasecalc(double, double, int, long int);
void plastic();
/**/
/* End FUNCTONS PROTOTYPES ------------------------------------------- */


