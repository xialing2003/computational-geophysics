/* 
class_tsunami.c
Program to simulate tsunami waves on a 2D Cartesian grid. The program uses
a 4th-order finite difference solution of the equation
Ptt = div * GH grad P
where P is the height of the Tsunami wave above sea level, G is the
acceleration due to gravity (a constant), and H is the ocean depth.
The speed of the wave is sqrt(GH).

R. Clayton, Caltech, Jan 2005

COMPILE:
  gcc -o class_tsunami class_tsunami.c -lm

RUN:
  ./class_tsunami
*/
#include	<stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h> 

float *v, *p1, *p2;
float *f1, *f2;
float *vx, *vy;
int	ord	= 2;
int	nx	=1000;
int	ny	=800;
int	nt	=3000;
float	h	=10.0;
float	dt	=10.0;
char	vmodel[]	="bathymetry.in";
char	output[]	="slices_homo_2nd.out";
int	itprint	=10;  /* time steps between print messages */
int	itslice	=100;  /* time steps between slice outputs */
float	latref	=-40.0;
float	lonref	=35;
float	slat	=3.30;
float	slon	=95.87;

#define V(ix,iy)	v[(ix) +nx*(iy)]
#define P1(ix,iy)	p1[(ix) +nx*(iy)]
#define P2(ix,iy)	p2[(ix) +nx*(iy)]
#define Vx(ix,iy)   vx[(ix) +nx*(iy)]
#define Vy(ix,iy)   vy[(ix) +nx*(iy)]
#define gauss(ix,iy) gauss[(ix) + 3*(iy)]
/* 2nd order second-derviative coefficients */
#define B1	-2.0	
#define B2	 1.0	
/* 4th order second-derviative coefficients */
#define C1	-2.5
#define C2	 1.3333333333
#define C3	-0.0833333333
/* 4th order first-derviative coefficients */
#define D1	-0.0833333333
#define D2	 0.6666666667
/* to decide if this is the homogeneous case or the heterogeneous case*/
#define flag 0
/* define smoothing values here*/
#define exp0 0.3319106653
#define exp1 0.1221031099
#define exp2 0.0449192238

float	mindepth	= 10.0; /* minimum depth to consider still ocean (m) */
#define G 0.00985	/* acc of gravity in km/sec**2 */
int ixref	=0;
int iyref	=0;

int main(int ac, char **av){
	int it, ix, iy, fd;
	int ixs, iys;
	float *tmp, vel, f, val, velmax, Vx_max=0, Vy_max=0, V_max=0;
	float Px,Py,V2;
	double norm(), sqrt();
	float gauss[9] = {exp2, exp1, exp2, exp1, exp0, exp1, exp2, exp1, exp2};

	fprintf(stderr,"order= %d\n",ord);

	v= (float *)(malloc(4*nx*ny));
	f1= (float *)(malloc(4*nx*ny));
	f2= (float *)(malloc(4*nx*ny));
    vx = (float *)(malloc(4*nx*ny));
    vy = (float *)(malloc(4*nx*ny));

	if(v == NULL || f1 == NULL || f2 == NULL)
	   {
		fprintf(stderr,"cannot alloc memory\n");
		exit(-1);
	   }
	if( (fd= open(vmodel,0)) < 0)
	   {
		fprintf(stderr,"cannot open velocity file=%s\n",vmodel);
		exit(-1);
	   }
	if( read(fd,v,4*nx*ny) != 4*nx*ny )
	   {
		fprintf(stderr,"read error in velocity file=%s\n",vmodel);
		exit(-1);
	   }
	close(fd);
	output_slice(v,nx,ny,-1.0);

	/* convert depth to velocity v= sqrt(g*depth).
	 * set values for land (pos. depths) to negative as flag
	 */
	for(iy=1; iy < ny-1; iy++)
	for(ix=1; ix < nx-1; ix++)
	{
		V(ix,iy) = V(ix-1,iy-1)*gauss(0,0) + V(ix-1,iy)*gauss(0,1) + V(ix-1,iy+1)*gauss(0,2) + 
					V(ix,iy-1)  *gauss(1,0) + V(ix,iy)  *gauss(1,1) + V(ix,iy+1)  *gauss(1,2) + 
					V(ix+1,iy-1)*gauss(2,0) + V(ix+1,iy)*gauss(2,1) + V(ix+1,iy+1)*gauss(2,2);
	}
	
	velmax= 0.0;
	for(iy=0; iy<ny; iy++)
	for(ix=0; ix<nx; ix++)
	   {
		val= -V(ix,iy); /* make depth positive */
		/* note 0.001 to convert depth (m) to km */
		if(val > mindepth) vel= sqrt(G*val*0.001);
		 else		   vel= 0.00001;
		if(vel > velmax) velmax= vel;
		if(vel > 0.0) V(ix,iy)= vel;//*vel*dt*dt/(h*h);
		 else	      V(ix,iy)= 0.00001;  
	   }
	fprintf(stdout,"maximum velocity = %8.4f (km/s)\n",velmax);
	fprintf(stdout,"nx= %d ny=%d nt=%d h=%8.4f dt=%8.4f\n",nx,ny,nt,h,dt);
	
    /*compute the first derivatives of velocity Vx = 2vdvdx, Vy = 2vdvdy*/
	for(iy=0; iy<ny; iy++)
	{
		Vx(0,iy) = 0;
		Vx(nx-1,iy) = 0;
	}
    for(ix=0; ix<nx; ix++)
	{
		Vx(ix,0) = 0;
		Vy(ix,ny-1) = 0;
	}
    for(iy=1; iy<ny-1; iy++)
    for(ix=1; ix<nx-1; ix++)
        {
			if(V(ix,iy)<0){
				Vx(ix,iy) = 0;
				Vy(ix,iy) = 0;
				continue;
			}
            if(ix == 1 || ix == nx-2 || iy ==1 || iy == ny-2)
			{
				Vx(ix,iy) = V(ix,iy)*(V(ix+1,iy)-V(ix-1,iy))/h;
                Vy(ix,iy) = V(ix,iy)*(V(ix,iy+1)-V(ix,iy-1))/h;
			}else{
				Vx(ix,iy) = 2*V(ix,iy)*(D1 * (V(ix+2,iy)-V(ix-2,iy)) + D2 * (V(ix+1,iy)-V(ix-1,iy)))/h;
                Vy(ix,iy) = 2*V(ix,iy)*(D1 * (V(ix,iy+2)-V(ix,iy-2)) + D2 * (V(ix,iy+1)-V(ix,iy-1)))/h;
			}
        }

	/* point the memory planes to real memory and zero it */
	p1= f1;
	p2= f2;
	zap(p1,nx*ny);
	zap(p2,nx*ny);

	/* add source */
	xcoord_convert(slat,slon,&ixs,&iys);
	fprintf(stderr,"source %8.3f %9.3f %4d %4d\n",slat,slon,ixs,iys);
	/*  source is placed on a grid:
	 *    1/4  1/2  1/4
	 *    1/2   1   1/2
	 *    1/4  1/2  1/4
	 */
	P2(ixs,iys)= 1.0;
	P2(ixs +1,iys  )= P2(ixs -1,iys  )= P2(ixs   ,iys+1)= P2(ixs   ,iys-1)= 0.5;
	P2(ixs+1,iys+1)= P2(ixs-1,iys+1)= P2(ixs+1,iys-1)= P2(ixs-1,iys-1)= 0.25;
	fprintf(stderr, "Vx(ixs,iys) %8.4f Vy(ixs,iys) %8.4f \n", Vx(ixs,iys), Vy(ixs,iys));
	
	float maxP;
	/* loop over time steps */
	for(it= 0; it<nt; it++)
	   {
		maxP = 0;
		/* loop over x-y plane */
		for(iy=1; iy < ny-1; iy++)
		for(ix=1; ix < nx-1; ix++)
		   {
			/* ignore points on land */
			if(V(ix,iy) < 0.0)
			   {
				P1(ix,iy)= 0.0;
				continue;
			   }
			if(ord == 2 || ix == 1 || ix == nx-2 || iy ==1 || iy == ny-2)
			   {
				/* 2nd-order */
				P1(ix,iy)= 2.0*(1.0+B1*V(ix,iy)*V(ix,iy)*dt*dt/h/h)*P2(ix,iy) - P1(ix,iy)
					 + B2*V(ix,iy)*V(ix,iy)*dt*dt/h/h*(P2(ix+1,iy)+P2(ix-1,iy)+P2(ix,iy+1)+P2(ix,iy-1));
				if(flag)
				{
					P1(ix,iy) += 0.5/h * (Vx(ix,iy) * (P2(ix+1,iy)-P2(ix-1,iy))) + 
										Vy(ix,iy) * (P2(ix,iy+1)-P2(ix,iy-1));
				}
			   }
			 else
			   {
				/* 4th-order interior */
				// fprintf(stdout,"Add 4th order code here\n");
				V2 = V(ix,iy) * V(ix,iy) * dt*dt/h/h;
                P1(ix,iy)= 2.0*(1.0+C1*V2)*P2(ix,iy) - P1(ix,iy)
					 + C2*V2*(P2(ix+1,iy)+P2(ix-1,iy)+P2(ix,iy+1)+P2(ix,iy-1))
                     + C3*V2*(P2(ix+2,iy)+P2(ix-2,iy)+P2(ix,iy+2)+P2(ix,iy-2));
				// exit(-1);
				if(flag)
				{
					Px = D1/h*(P2(ix+2,iy)-P2(ix-2,iy)) + D2/h*(P2(ix+1,iy)-P2(ix-1,iy));
					Py = D1/h*(P2(ix,iy+2)-P2(ix,iy-2)) + D2/h*(P2(ix,iy+1)-P2(ix,iy-1));
					P1(ix,iy) += Vx(ix,iy) * Px	+ Vy(ix,iy) * Py;
				}
			   }
			if(P1(ix,iy) > maxP) maxP = P1(ix,iy);
		   }
	
		/* Dirichlet boundary conditions */
		for(ix=0,    iy=0;    ix<nx; ix++) P1(ix,iy)= 0.0;
		for(ix=0,    iy=ny-1; ix<nx; ix++) P1(ix,iy)= 0.0;
		for(ix=0,    iy=0;    iy<ny; iy++) P1(ix,iy)= 0.0;
		for(ix=nx-1, iy=0;    iy<ny; iy++) P1(ix,iy)= 0.0;

		if(it%itprint == 0)
			fprintf(stderr,"done it=%3d norm=%14.3e maxP = %8.4f\n",it,
				norm(&P1(0,0),nx*ny),maxP);
		if(it%itslice == 0)
			output_slice(p1,nx,ny,(double)(it*dt));
		/* rotate the memory pointers */
		tmp= p1; p1= p2; p2= tmp;
	   }
   }

#define DEG2KM	111.195
#define DEG2R	  0.0174532
void xcoord_convert(double slat,double slon,int *ixs,int *iys)
   {
	/* convert (lat,lon) to grid coords (ixs,iys) */
	double cos();
	double x, y;
	y= (slat-latref)*DEG2KM;
	x= (slon-lonref)*DEG2KM*cos(slat*DEG2R);
	*ixs = x/h + ixref;
	*iys = y/h + iyref;
   }

double norm(float *x, int n)
   {
	/* compute the norm of a vector or plane */
	double sum, val, sqrt();
	int i;

	sum= 0.0;
	for(i=0; i<n; i++) sum += x[i]*x[i];
	val= sqrt( sum/(double)(n) );
	return(val);
   }

void zap(float *x, int n)
   {
	/* zero out a field */
	int i;
	for(i=0; i<n; i++) x[i]= 0.0;
   }

double
getmax(float *x, int n)
   {
	/* find the absolute max of a plane */
	int i;
	double fabs(), max;
	max= fabs(x[0]);
	for(i=1; i<n; i++) if(fabs(x[i]) > max) max= fabs(x[i]);
	return(max);
   }

/* output a slice of the field.
 * Note: the field is output at every istep samples
 *       the field is reversed in y
 */
float line[1000];
int outfd	= -1;
int	istep	= 2;
void output_slice(float *x,int nx,int ny,double t)
   {
	double max, getmax(), val;
	int ix, iy, i, ival;

	if(outfd < 0) outfd= creat(output,0664);
	if(outfd < 0)
	   {
		fprintf(stderr,"cannot create plot file= %s\n",output);
		exit(-1);
	   }
	for(iy=ny-1; iy >= 0; iy -= istep)
	   {
		for(ix=0, i=0; ix <nx; ix += istep, i++)
			line[i]= x[ix+iy*nx];
		write(outfd,line,4*i);
	   }
   }
