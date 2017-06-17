/************************************************************
   fleet.c by Travis Burkitt, November 29, 1988
   -- extensively modified by John Barron, 1991
   -- This program is the result of "sticking" three
      separate programs together. Hence, it uses less
      computation time but requires more space. If you
      need more space, use limit to increase the datasize.
      limit -h gives the hard limits of your system.
************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <fcntl.h>
#include <time.h>
#include <sys/types.h>


#define PI M_PI
#define PIC_X  325
#define PIC_Y  325
#define PIC_T  22
#define PMODE 0644
#define TRUE 1
#define FALSE 0
#define BETA 0.001 
#define KERNEL_X 20
#define KERNEL_Y 20
#define KERNEL_T 20
#define NUMBER 5
#define HEAD 32
#define THRESHOLDED 0
#define NOT_THRESHOLDED 1
#define TORAD (PI/180.0)
#define NRADIUS 4
#define NO_NORM_VEL (NRADIUS*2+1)*(NRADIUS*2+1)*22 
#define NO_UNKNOWNS 6
#define UNDEFINED 100.0
#define NO_VALUE 100.0
#define END_OF_VALUES 1000000.0
#define REMOVE 0 /* Remove Previous Filtered Files? */
#define PRINT 0


typedef unsigned char cimages[PIC_T][PIC_X][PIC_Y];
typedef float images[NUMBER][PIC_X][PIC_Y];
typedef float imageplane[PIC_X][PIC_Y]; 
typedef struct { float Re, Im;} Imag,*ImagP;
typedef int int3[3];
typedef float float3[3];
typedef Imag Imag3[3];





cimages pic;
images IfRe,IfIm,savet,savex,savett,DC;
int NUMFILES,OFFSET_X,OFFSET_Y,OFFSET_T,STARTFILE,OLD_WAY;
float maxamp;
long time1,time2,time3,time4,temporal_offset;
FILE *fp,*fpa;
int pic_x,pic_y,pic_t;
unsigned char thresholded[22][PIC_X][PIC_Y];
float E[PIC_X][PIC_Y];
float dot(),L2norm(),PsiER(),PsiEN();
double cal_velocity();
float normal_velocities[22][PIC_X][PIC_Y][2];

float full_velocities[PIC_X][PIC_Y][2];
float correct_velocities[PIC_X][PIC_Y][2],tauF,tauA;
int normal_ct[PIC_X][PIC_Y];
double cal_inverse();
int fdn,fdf,RUN1,RUN2,BINARY;




/******************************************************************/
main(argc,argv)
int argc;
char **argv;
{
char command[100],correct_filename[100];
char filename[100],path1[100],path2[100],path3[100];
float sigma,tau,percent_maxamp;
int i,central_image,size;

if(argc < 6 || argc > 17)
	{
	printf("Usage: %s <file stem> <sigma> <central image> <percent maxamp> \n",argv[0]);
	printf("       <input path> <filter path> <output path> [-C <full correct filename>] [-NF or -F] [-V or -NV] [-B <cols> <rows>] [-O <tauA> <tauF>]\n");
	printf("<filename> - the filename stem\n");
	printf("<sigma> - spatio-temporal sigma for the filters\n");
	printf("<central image> - the image the velocities are computed for\n");
	printf("<percent_amp> - percentage of maximum amplitude\n");
	printf("<input path> - where the input images are\n");
	printf("<filter path> - where the filter results are put or retrieved from\n");
	printf("<output path> - where the output goes\n");
	printf("-C <correct filename> - the correct velocity file name\n");
	printf("   default is unknown - no error analysis is done in this case\n");
	printf("-F  - deflault if absent - perform gabor filtering\n");
	printf("-NF - don't perform filtering\n");
	printf("-V  - deflault if absent - perform thresholding and velocity computation\n");
	printf("-NV - don't perform thresholding and velocity computation\n");
	printf("-B <cols> <rows> - input the size of binary files: assumes the files\n");
	printf("   are not black and white SUN rasterfiles\n");
	printf("-O <tauA> <tauF> - perform thresholding the old way - supply\n");
	printf(" amplitude and frequency thresholds\n");
	printf("-N <tau> - perform thresholding the new way using a single threshold\n");
	printf("           Default:  -N 1.25\n");
	printf("NOTE: non-minus arguments must be given in the specified order\n");
	printf("%d arguments specified\n",argc);
	exit(1);
	}
else
{
printf("Command line:");
for(i=0;i<argc;i++) printf("%s ",argv[i]);
printf("\n");
}

/* Compute size of the input data from sigma */
sscanf(argv[2],"%f",&sigma);
sscanf(argv[3],"%d",&central_image);
sscanf(argv[4],"%f",&percent_maxamp);
strcpy(path1,".");
strcpy(path2,".");
strcpy(path3,".");
printf("sigma=%f\n",sigma);
printf("Central image: %d\n",central_image);
printf("Amplitude threshold: %f percent\n",percent_maxamp*100.0);
strcpy(path1,argv[5]); 
strcpy(path2,argv[6]); 
strcpy(path3,argv[7]); 
printf("<input path>: %s\n<filter path>: %s\n<output path>: %s\n",
	path1,path2,path3);


size = 6*sigma+1;
if(size % 2 == 0) size += 1;
OFFSET_X = OFFSET_Y = OFFSET_T = size/2.0;
NUMFILES = size+4; /* Filter 5 adjacent images */
pic_t = 5;
if(NUMFILES > PIC_T || NUMBER > NUMFILES-2*OFFSET_T || NUMBER < pic_t) 
	{ 
	printf("Fatal error: not enough room for image data\n");
	printf("NUMFILES: %d PIC_T: %d OFFSET_T: %d\n",
		NUMFILES,PIC_T,OFFSET_T);
	exit(1);
	}
	
STARTFILE = central_image-size/2-2;
printf("Offsets in x, y and t: %d\n",size/2);
printf("NUMFILES: %d STARTFILE: %d\n",NUMFILES,STARTFILE);
fflush(stdout);
RUN1 = TRUE; /* Perform filtering */
RUN2 = TRUE; /* Perform thresholding and compute the image velocities */
BINARY = FALSE;
OLD_WAY = UNDEFINED;
strcpy(correct_filename,"unknown");
i = 9;
while(i <= argc)
	{
	if(strcmp(argv[i-1],"-NF")==0) RUN1=FALSE;
	else
	if(strcmp(argv[i-1],"-NV")==0) RUN2=FALSE;
	else
	if(strcmp(argv[i-1],"-C")==0)
		{
		strcpy(correct_filename,argv[i]);
		i++;
		}
	else
	if(strcmp(argv[i-1],"-B")==0)
		{
		sscanf(argv[i],"%d",&pic_y);
		sscanf(argv[i+1],"%d",&pic_x);
		i += 2;
		BINARY = TRUE;
		}
	else
	if(strcmp(argv[i-1],"-O")==0)
		{
		sscanf(argv[i],"%f",&tauA);
		sscanf(argv[i+1],"%f",&tauF);
		i += 2;
		OLD_WAY = TRUE;
		}
	else
	if(strcmp(argv[i-1],"-N")==0)
		{
		sscanf(argv[i],"%f",&tau);
		i++;
		OLD_WAY = FALSE;
		}
	i++;
	}
if(RUN1==FALSE && RUN2==FALSE)
	{
	printf("Nothing is done - program quits\n");
	exit(1);
	}
if(OLD_WAY == UNDEFINED)
	{
	printf("Default threshold: tau=1.25 used\n");
	tau = 1.25;
	OLD_WAY = FALSE;
	exit(1);
	}
if(RUN1) printf("Gabor filtering specified\n");
else printf("No Gabor filtering - using previous results\n");
if(RUN2) printf("Thresholding and Velocity Computation specified\n");
printf("Correct velocity file name: %s\n",correct_filename);
fflush(stdout);


/* Perform timing measurements */
temporal_offset = 100; /* Each program requires at 100 seconds */
sprintf(filename,"fleet.%stime",argv[1]);
fp = fopen(filename,"w");
fprintf(fp,"The time data is in file: %s\n",filename);
time1 = time(NULL);
fprintf(fp,"Start time before convolution: %d\n",time1);
fflush(fp);
time2 = time1;
time3 = time1;
time4 = time1;

if(RUN1)
{
if(REMOVE)
{
sprintf(command,"rm %s/%sfilter*",path2,argv[1]);
printf("%s\n",command);
system(command);
}
fflush(fp); 
/* Filter the images */
filter_images(path1,path2,argv[1],sigma);

time2 = time(NULL);
if(time2 > time1+temporal_offset)
{
fprintf(fp,"\nEnd time for the convolution: %d\n",time2);
fprintf(fp,"Time in seconds  needed for convolution: %d\n",
	(time2-time1));
fprintf(fp,"Time in minutes needed for convolution: %f\n",
	(time2-time1)*1.0/60.0);
fprintf(fp,"Time in hours needed for convolution: %f\n",
	(time2-time1)*1.0/3600.0);
fflush(fp);
}
}
/* have to set pic_x and pic_y values anyway for threholding
   and velocity calculation using the central image if
   black and white rasterfile images. */
else if(BINARY==FALSE)
    {
    set_dimensions(path1,argv[1],central_image);
    fprintf(fp,"\nConvolution not performed\n"); 
    fflush(fp);
    }


if(RUN2)
{
/* Threshold filtered results and compute velocities */
thresh_and_compute(argv[1],path2,path3,correct_filename,sigma,
	tau,tauA,tauF,percent_maxamp);

time3 = time(NULL);
if(time3 > time2+temporal_offset)
{
fprintf(fp,"\nEnd time for the thesholding: %d\n",time3);
fprintf(fp,"Time in seconds  needed for thresholding: %d\n",
	(time3-time2));
fprintf(fp,"Time in minutes needed for thresholding: %f\n",
	(time3-time2)*1.0/60.0);
fprintf(fp,"Time in hours needed for thresholding: %f\n",
	(time3-time2)*1.0/3600.0);
fflush(fp);
}

time4 = time(NULL);
fprintf(fp,"\n\nStart time: %d\n",time1);
fprintf(fp,"Convolution end time: %d\n",time2);
fprintf(fp,"Thresholding and Image Velocity Calculation end time: %d\n",time3);
fprintf(fp,"End time: %d\n",time4);
fprintf(fp,"\n\nElasped time in seconds: %d\n",time4-time1);
fprintf(fp,"Elasped time in minutes: %f\n",(time4-time1)*1.0/60.0);
fprintf(fp,"Elasped time in hours: %f\n",(time4-time1)*1.0/3600.0);
fflush(fp);
}
else { 
     fprintf(fp,"\nThresholding and Velocity Computation not performed\n");
     fflush(fp); 
     }
return(1);
}

/*********************************************************************/

filter_images(path1,path2,s,sigma)
char s[100],path1[100],path2[100];
float sigma;
{
float maxamp=0.0,tuning_info[22][2];
float filter[23][3];
char name[100];



sprintf(name,"%s/%s",path1,s);
printf("\nReading Files...\n");
read_image_files(name);
initfilters(filter,tuning_info,sigma);
printf("Files read and filters initialized\n");
fflush(stdout);

all_filters(s,path2,&maxamp,filter,sigma);
printf("Max amplitude == %15.10f\n",maxamp);
} /* End of main */


/************************************************************
-Convolves images with each filter and computes amplitudes. 
-writes files.
************************************************************/
all_filters(name,path,maxamp,filter,sigma)
float *maxamp,filter[23][3],sigma;
char name[100],path[100];
{
char name1[100];
int n;
float  Gsx[2*KERNEL_X+1], Gsy[2*KERNEL_Y+1], Gst[2*KERNEL_T+1],
       Gcx[2*KERNEL_X+1], Gcy[2*KERNEL_Y+1], Gct[2*KERNEL_T+1],
       Gx[2*KERNEL_X+1], Gy[2*KERNEL_Y+1], Gt[2*KERNEL_T+1];


/* Convolve image with a Gaussian  --- Subtract from cosine filter output */
/* We do not use the flicker channel */

(*maxamp) = 0.0;
for(n=0;n<22;n++)  /* Start of n loop */             
	{ 
	printf("\n\nCalculating filter %d...\n",n);    
	calc_gabors(Gsx,Gsy,Gst,Gcx,Gcy,Gct,
		   filter[n][0],filter[n][1],filter[n][2],sigma);
	printf("k1=%f k2=%f k3=%f\n",filter[n][0],filter[n][1],filter[n][2]);
	if(PRINT) coefficients_and_sums(Gsx,Gsy,Gst,Gcx,Gcy,Gct,sigma);	
	printf("Convolving Gabors...\n");      
	fflush(stdout);
	convolve_gabors(Gsx,Gsy,Gst,Gcx,Gcy,Gct);

	/***** subtract DC component, 0.001 G from cos output  *****/
	calc_gaussian(Gx,Gy,Gt,sigma);
	calc_DC(Gx,Gy,Gt);
	subpic(IfRe,DC); 

	calc_maxamp(maxamp);
	sprintf(name1,"%s%s",name,"filter");
	write_If_data(name1,path,n,IfRe,IfIm);
	} /* End of n loop */
/* Save the maximum amplitude in a file in case the program is run without
   filtering */
sprintf(name1,"%samp",name);
fpa = fopen(name1,"w");
fprintf(fpa,"%f",(*maxamp));
fclose(fpa);
} /* End of all_filters */


/************************************************************
  Computes amplitudes of each filter output and saves in 
  intermediate files. Saves a running maximum amplitude.
  USES GLOBALS images IfRe,IfIm;
************************************************************/
calc_maxamp(maxamp)
float *maxamp;
{
int i,j,k;
float amplitude;
float filtmax = 0.0;	
	
i = 2; /* Use middle image */
for(j=0;j<pic_x;j++) 
	{
       	for (k=0;k<pic_y;k++)
		{
		amplitude = sqrt((float)(IfRe[i][j][k]*IfRe[i][j][k] +
					 IfIm[i][j][k]*IfIm[i][j][k]));
		if (amplitude > (*maxamp))
			(*maxamp) = amplitude;
		if (amplitude > filtmax)
			filtmax = amplitude;
		}
	}
printf("Current Max amplitude == %12.8f\n",(*maxamp));
printf("Filter Max amplitude == %12.8f\n",filtmax);
} /* End of calc_amp */


/*******************************************************************/
/* Compute the Gaussian Coefficients 				   */
/*******************************************************************/
calc_gaussian(Gx,Gy,Gt,sigma)
float Gx[],Gy[],Gt[],sigma;
{
int i;
float val,twopi;

twopi = sqrt((float)(2.0*PI));

for(i=(-OFFSET_T);i<=OFFSET_T;i++) 
	{
        val = i;
        Gt[i+OFFSET_T] = 1.0/(twopi*sigma) * exp(-(val*val/(2.0*sigma*sigma)));
        }
 
for(i=(-OFFSET_X);i<=OFFSET_X;i++) 
	{
        val = i;
        Gx[i+OFFSET_X] = 1.0/(twopi*sigma) * exp(-(val*val/(2.0*sigma*sigma)));
        }
 
for(i=(-OFFSET_Y);i<=OFFSET_Y;i++) 
	{
        val = i;
        Gy[i+OFFSET_Y] = 1.0/(twopi*sigma) * exp(-(val*val/(2.0*sigma*sigma))); 
	}
} /* End of calc_gaussian */


/************************************************************
  Compute values for the 1-D sine and cosine gabors. 
************************************************************/
calc_gabors(Gsx,Gsy,Gst,Gcx,Gcy,Gct,wx,wy,wt,sigma)
float Gsx[],Gsy[],Gst[],Gcx[],Gcy[],Gct[];
float wx,wy,wt,sigma;
{
int i;
float val,twopi,constant;
int SWITCH;

SWITCH = 2; /* 0 or 2 to reverse the order of the coefficients */ 

twopi = sqrt((float)(2.0*PI));
constant = 1.0;

printf("\n");
for(i=(-OFFSET_T);i<=OFFSET_T;i++) 
	{
	val = i;
	Gst[i+OFFSET_T-SWITCH*i] = 1.0/(twopi*sigma)* sin(constant*val*wt) *
				        exp(-(val*val/(2.0*sigma*sigma)));
	}
for(i=(-OFFSET_X);i<=OFFSET_X;i++) 
	{
	val = i;
        Gsx[i+OFFSET_X-SWITCH*i] = 1.0/(twopi*sigma)* sin(constant*val*wx) *
				        exp(-(val*val/(2.0*sigma*sigma)));
	}
for(i=(-OFFSET_Y);i<=OFFSET_Y;i++) 
	{
	val = i;
        Gsy[i+OFFSET_Y-SWITCH*i] = 1.0/(twopi*sigma)* sin(constant*val*wy) *
				        exp(-(val*val/(2.0*sigma*sigma))); 
	} 

/* Cosine gabors */ 
for(i=(-OFFSET_T);i<=OFFSET_T;i++) 
	{
	val = i;
        Gct[i+OFFSET_T-SWITCH*i] = 1.0/(twopi*sigma)* (cos(constant*val*wt)) *
				        exp(-(val*val/(2.0*sigma*sigma)));
	}
for(i=(-OFFSET_X);i<=OFFSET_X;i++) 
	{
	val = i;
        Gcx[i+OFFSET_X-SWITCH*i] = 1.0/(twopi*sigma)* (cos(constant*val*wx)) *
				        exp(-(val*val/(2.0*sigma*sigma)));
	}
for(i=(-OFFSET_Y);i<=OFFSET_Y;i++) 
	{
	val = i;
        Gcy[i+OFFSET_Y-SWITCH*i] = 1.0/(twopi*sigma)* (cos(constant*val*wy)) *
				        exp(-(val*val/(2.0*sigma*sigma)));
	}
} /* End of calc_gabors */


/***********************************************************************/
/* Remove the DC component, i.e. 0.001 of the Gaussian		       */
/***********************************************************************/
calc_DC(Gx,Gy,Gt)
float Gx[],Gy[],Gt[]; 
{
printf("\n\nComputing DC component as 0.001G ...\n");
conv_t(Gt,DC,pic);
conv_x(Gx,savet,DC);
conv_y(Gy,DC,savet);
mult(DC,BETA);
} /* End of calc_DC */



/************************************************************
 Convolve images with several 1-D sine/cosine gabor
 filters.  This is equivalent to convoluting with a 3-D
 filter.    ** see HEEGER 87 -- appendix B **
  USES GLOBALS  cimages pic
		images  IfRe,IfIm;
		images 	savet, savex;
************************************************************/
convolve_gabors(Gsx,Gsy,Gst,Gcx,Gcy,Gct)
float Gsx[],Gsy[],Gst[],Gcx[],Gcy[],Gct[];
{
/* Do 3D convolution as sums/products of 1D convolutions */


conv_t(Gct,savet,pic);
conv_x(Gcx,savex,savet);
conv_y(Gsy,IfIm,savex);
conv_y(Gcy,IfRe,savex); 

conv_x(Gsx,savex,savet);
conv_y(Gcy,savet,savex);
addpic(IfIm,savet); 
conv_y(Gsy,savet,savex);
subpic(IfRe,savet); 

conv_t(Gst,savett,pic);
conv_x(Gsx,savex,savett);
conv_y(Gsy,savet,savex);
subpic(IfIm,savet); 
conv_y(Gcy,savet,savex);
subpic(IfRe,savet);

conv_x(Gcx,savex,savett);
conv_y(Gcy,savet,savex);
addpic(IfIm,savet); 
conv_y(Gsy,savet,savex);
subpic(IfRe,savet); 

fflush(stdout);
} /* End of convolve_gabors */



/************************************************************
  Convolve 1d gabor in time
************************************************************/
conv_t(G,save,pic)
float G[];
images save;
cimages pic;
{
register int i,j,k,m;

for(i=0;i<NUMFILES-2*OFFSET_T;i++)
for(j=OFFSET_X;j<pic_x-OFFSET_X;j++) 
for(k=OFFSET_Y;k<pic_y-OFFSET_Y;k++) 
	{
	save[i][j][k] = 0;
	for(m=i;m<=i+2*OFFSET_T;m++) 
		{
	        save[i][j][k] += (pic[m][j][k] * G[m-i]);
		}
	}
} /* End of conv_t */


/************************************************************
  Convolve 1d gabor in x
************************************************************/
conv_x(G,save,pic)
float G[];
images save,pic;
{
register int i,j,k,m,start;

for(i=0;i<NUMFILES-2*OFFSET_T;i++)
for(j=OFFSET_X;j<pic_x-OFFSET_X;j++) 
for(k=OFFSET_Y;k<pic_y-OFFSET_Y;k++)  
	{ 
	save[i][j][k] = 0;
	start = j-OFFSET_X;
	for (m=start;m<=j+OFFSET_X;m++) 
		{ 
		save[i][j][k] += (pic[i][m][k] * G[m-start]);
                } 
	}
} /* End of conv_x */


/************************************************************
  Convolve 1d gabor in y
************************************************************/
conv_y(G,save,pic) 
float G[]; 
images save,pic;
{
register int i,j,k,m,start;

for(i=0;i<NUMFILES-2*OFFSET_T;i++)	
for(j=OFFSET_X;j<pic_x-OFFSET_X;j++) 
for(k=OFFSET_Y;k<pic_y-OFFSET_Y;k++)  
	{ 
	save[i][j][k] = 0;
	start = k-OFFSET_Y;
	for(m=start;m<=k+OFFSET_Y;m++)
		{
		save[i][j][k] += (pic[i][j][m] * G[m-start]);
                }
	}
} /* End of conv_y */


/************************************************************
  calculate res = res + data
      Used to add results of two convolutions.
************************************************************/
addpic(res,data)
images res,data; 
{
register int i,j,k;

for(i=0;i<NUMFILES-2*OFFSET_T;i++)  /*******/
for(j=0;j<pic_x;j++)
for(k=0;k<pic_y;k++)  
	{
	res[i][j][k] += data[i][j][k];
        }
}


/************************************************************
  calculate res = res - data
      Used to subtract results of two convolutions.
************************************************************/
subpic(res,data)
images res,data; 
{
register int i,j,k;

for(i=0;i<NUMFILES-2*OFFSET_T;i++)     
for(j=0;j<pic_x;j++)
for(k=0;k<pic_y;k++) 
	{
	res[i][j][k] -= data[i][j][k];
	}
}



/************************************************************
   Multiply elements of 3D array by a constant
************************************************************/
mult(pic,val)
images pic;
float val;
{
register int i,j,k;
for(i=0;i<NUMFILES-2*OFFSET_T;i++)    
for(j=0;j<pic_x;j++)
for(k=0;k<pic_y;k++)  
	{
	pic[i][j][k] = pic[i][j][k] * val; 
	}
}




/************************************************************
   Read raster files into internal 3-D array.
   USES GLOBAL: cimages pic
   Also sets pic_x and pic_y (the actual picture dimensions)
   which are global
************************************************************/
read_image_files(s)
char *s;
{
unsigned char header[HEAD];
char fname[100];
int i,j,fd,ints[8],ONCE,bytes;

ONCE = TRUE;
bytes = 0;
for(i=0;i<NUMFILES;i++) 
	{
	sprintf(fname,"%s%d",s,i+STARTFILE);
 	if((fd = open(fname,O_RDONLY)) >0)
		{
		if(!BINARY)
		{
		if(ONCE)
			{
			read(fd,ints,HEAD);
			pic_y = ints[1];
			pic_x = ints[2];
			ONCE = FALSE;
			if(pic_y > PIC_Y || pic_x > PIC_X)
				{
				printf("Fatal error - not enough room\n");
				printf("Required size: %d times %d\n",pic_y,pic_x);
				exit(1);
				}
			}
		else read(fd,header,HEAD);
		bytes = 32;
		}

		/* Read row by row */
		for(j=0;j<pic_x;j++)
			bytes += read(fd,&pic[i][j][0],pic_y);
		printf("File %s read -- %d bytes\n",fname,bytes);
		fflush(stdout);
		}
	      else 
		{
		printf("File %s does not exist in read_image_files.\n",fname);
		exit(1);
		}
	}
printf("Size of input data: %d columns, %d rows\n",pic_y,pic_x);
} /* End of read_image_files */

/***************************************************************/
/* Given image n, read the width and height dimensions         */
/***************************************************************/
set_dimensions(path,name,n)
char path[100],name[100];
{
char fname[100];
int fd,ints[8];

sprintf(fname,"%s/%s%d",path,name,n);
if((fd = open(fname,O_RDONLY)) >0)
	{
	read(fd,ints,HEAD);
	pic_y = ints[1];
	pic_x = ints[2];
	}

printf("Size of input data: %d columns, %d rows\n",pic_y,pic_x);
printf("Size data obtained from %s\n",fname);
} /* End of set_dimensions  */


 

/************************************************************************/
/* Print out the coefficients and their sums for each sine and cosine   */
/* 1d gabor. Useful for debugging.					*/
/************************************************************************/
coefficients_and_sums(Gsx,Gsy,Gst,Gcx,Gcy,Gct,sigma)
float Gsx[],Gsy[],Gst[],Gcx[],Gcy[],Gct[],sigma;
{
int i;
float sum;

printf("Value of sigma: %f\n",sigma);
/* Sine Gabors */
printf("\nCoefficients of Gsx:\n");
for(i=0;i<2*OFFSET_X+1;i++) 
	{ 
	if(i%8==0) printf("\n"); 
	printf("%f ",Gsx[i]); 
	}
sum=0.0;
for(i=0;i<2*OFFSET_X+1;i++) sum += Gsx[i];
printf("\nSum=%f\n",sum);

printf("\nCoefficients of Gsy:\n");
for(i=0;i<2*OFFSET_Y+1;i++) 
	{ 
	if(i%8==0) printf("\n"); 
	printf("%f ",Gsy[i]); 
	}
sum=0.0;
for(i=0;i<2*OFFSET_Y+1;i++) sum += Gsy[i];
printf("\nSum=%f\n",sum);

printf("\nCoefficients of Gst:\n");
for(i=0;i<2*OFFSET_T+1;i++) 
	{ 
	if(i%8==0) printf("\n"); 
	printf("%f ",Gst[i]); 
	}
sum=0.0;
for(i=0;i<2*OFFSET_T+1;i++) sum += Gst[i];
printf("\nSum=%f\n",sum);


/* Cosine Gabors */
printf("\nCoefficients of Gcx:\n");
for(i=0;i<2*OFFSET_X+1;i++) 
	{ 
	if(i%8==0) printf("\n"); 
	printf("%f ",Gcx[i]); 
	}
sum=0.0;
for(i=0;i<2*OFFSET_X+1;i++) sum += Gcx[i];
printf("\nSum=%f\n",sum);

printf("\nCoefficients of Gcy:\n");
for(i=0;i<2*OFFSET_Y+1;i++) 
	{ 
	if(i%8==0) printf("\n"); 
	printf("%f ",Gcy[i]); 
	}
sum=0.0;
for(i=0;i<2*OFFSET_Y+1;i++) sum += Gcy[i];
printf("\nSum=%f\n",sum);

printf("\nCoefficients of Gct:\n");
for(i=0;i<2*OFFSET_T+1;i++) 
	{ 
	if(i%8==0) printf("\n"); 
	printf("%f ",Gct[i]); 
	}
sum=0.0;
for(i=0;i<2*OFFSET_T+1;i++) sum += Gct[i];
printf("\nSum=%f\n\n\n",sum);

fflush(stdout);
} /* End of coefficients and suns */


/*********************************************************
  Threshold the filter output and then compute image 
  velocities using the thresholded results
*********************************************************/
thresh_and_compute(s,path2,path3,correct_name,sigma,
		   tau,tauA,tauF,percent_maxamp)
char s[100],path2[100],path3[100],correct_name[100];
float tau,percent_maxamp,sigma,tauA,tauF;
{
char name[100],fullname[100],normalname[100];
float filter[23][3],tuning_info[22][2];
float min_angle,max_angle,ave_error,st_dev,density,residual;
int num_thresh[22],n,int_size_x,int_size_y,no_bytes,fd_correct,i;
float size_x,size_y,actual_x,actual_y,offset_x,offset_y;
float norm_err,norm_st_dev,min_norm_angle,max_norm_angle;
float phi[PIC_X][PIC_Y][3];

printf("\nStarting thresholding\n");
fflush(stdout);


/* Read maxamp from file in current directory in case filtering is
   not done during this run (and hence maxamp has not been calculated) */
sprintf(name,"%samp",s);
if((fpa = fopen(name,"r"))<0)
	{
	printf("Fatal error: %s does't exist\n",name);
	exit(1);
	}
fscanf(fpa,"%f",&maxamp);
fclose(fpa);
printf("Maximum amplitude: %f\n",maxamp);
fflush(stdout);

initialize_thresholds(thresholded);

printf("Maximum amplitude threshold percent: %f\n",percent_maxamp*100.0);
fflush(stdout);

initfilters(filter,tuning_info,sigma);
sprintf(name,"%s/%s%s",path2,s,"filter");

if(!OLD_WAY)
threshold1(name,thresholded,maxamp,num_thresh,filter,percent_maxamp,tau); 
else
threshold2(name,thresholded,maxamp,num_thresh,filter,percent_maxamp,tauA,tauF);

produce_thesh_report(num_thresh,filter,tuning_info);
fflush(stdout);

if(!OLD_WAY)
sprintf(fullname,"%s/fleet.%s%dF-%3.2f",path3,s,(int)(sigma*10),tau);
else
sprintf(fullname,"%s/fleet.%s%dF-%3.2f-%3.2f",path3,s,(int)(sigma*10),tauA,tauF);
if((fdf = creat(fullname,0644))<0)
	{
        printf("Error creating file %s.\n",fullname);
        exit(1); 
	}
printf("\nFile %s opened\n",fullname);

if(!OLD_WAY)
sprintf(normalname,"%s/fleet.%s%dN-%3.2f",path3,s,(int)(sigma*10),tau);
else
sprintf(normalname,"%s/fleet.%s%dN-%3.2f-%3.2f",path3,s,(int)(sigma*10),tauA,tauF);
if((fdn = creat(normalname,0644))<0)
	{
        printf("Error creating file %s.\n",normalname);
        exit(1); 
	}
printf("\nFile %s opened\n",normalname);
fflush(stdout);

for(n=0;n<22;n++)
{
read_If_data(name,n,IfRe,IfIm);
printf("Real and Imaginary data read for filter %s%d\n",name,n);
fflush(stdout);
Dphi(IfRe,IfIm,phi,&filter[n][0],n);
printf("DELTA phi computed for filter %d\n",n);
fflush(stdout);
calc_normal_velocities(phi,normal_velocities,n,thresholded);
printf("Normal velocities computed for filter %d\n",n);
fflush(stdout);
}

output_normal_velocities(fdn,normal_velocities);
fflush(stdout);
printf("Starting full image velocity calculation\n");
fflush(stdout);
calc_full_velocities(normal_velocities,full_velocities,E); 
printf("Full image computed\n");
fflush(stdout);
output_full_velocities(fdf,full_velocities);
n = OFFSET_X;
if(OFFSET_Y > n) n = OFFSET_Y;

if(strcmp("unknown",correct_name)!=0)
{
/* Read the correct velocity data */
fd_correct = open(correct_name,O_RDONLY);
no_bytes = 0;
no_bytes += read(fd_correct,&actual_y,4);
no_bytes += read(fd_correct,&actual_x,4);
no_bytes += read(fd_correct,&size_y,4);
no_bytes += read(fd_correct,&size_x,4);
no_bytes += read(fd_correct,&offset_y,4);
no_bytes += read(fd_correct,&offset_x,4);
if(offset_x != 0.0 || offset_y != 0.0 || actual_x != size_x || actual_y != size_y)
	{
	printf("Fatal error: something wrong with correct velocity data\n");
	printf("Actual y: %f Actual x: %f\n",actual_y,actual_x);
	printf("Size y: %f Size x: %f\n",size_y,size_x);
	printf("Offset y: %f Offset x: %f\n",offset_y,offset_x);
	exit(1);
	}
int_size_y = size_y;
int_size_x = size_x;
for(i=0;i<int_size_x;i++)
	no_bytes += read(fd_correct,&correct_velocities[i][0][0],int_size_y*8);
printf("\nFile %s opened and read\n",correct_name);
printf("Size of correct velocity data: %d %d\n",int_size_y,int_size_x);
printf("%d bytes read\n",no_bytes);
fflush(stdout);

calc_statistics(correct_velocities,int_size_x,int_size_y,full_velocities,
		E,pic_x,pic_y,n,&ave_error,&st_dev,&density,&residual,
		&min_angle,&max_angle,normal_velocities,&norm_err,&norm_st_dev,
		&min_norm_angle,&max_norm_angle);
printf("Average angle error: %f Standard deviation: %f\n",ave_error,st_dev);
printf("Density: %f Residual: %f\n",density,residual);
printf("Minimum angle error: %f Maximum angle error: %f\n",min_angle,max_angle);
printf("Average normal angle error: %f Normal Standard Deviation: %f\n",
	norm_err,norm_st_dev);
printf("Minimum normal angle error: %f Maximum normal angle error: %f\n",
	min_norm_angle,max_norm_angle);
fflush(stdout);
}
printf("\nEnd of Program\n");
fflush(stdout);
}



/************************************************************
   Write FLOATS from 3-D internal array into files.
************************************************************/
write_If_data(name,path,n,Re,Im)   
char name[100],path[100];
int n;
images Re,Im;
{
char fname[100];
int i,j,fd;


sprintf(fname,"%s/%s%d",path,name,n);
printf("Writing file: %s in write_If_data\n",fname);
if((fd = creat(fname,0755)) >=0)
	{
        for (i=0;i<NUMFILES-2*OFFSET_T;i++) 
		{
		for(j=0;j<pic_x;j++)
		{
                if(write(fd,&Re[i][j][0],pic_y*sizeof(float))
			!=pic_y*sizeof(float))
			{
                        printf("Write error in write_If_data for Re: i=%d j=%d\n",i,j);
			exit(1);
			}
		}
		for(j=0;j<pic_x;j++)
		{
                if(write(fd,&Im[i][j][0],pic_y*sizeof(float))
			!=pic_y*sizeof(float))
			{
                        printf("Write error in write_If_data for Im: i=%d j=%d\n",i,j);
			exit(1);
			}
		}
		}
	}
     else
	{
	printf("Error creating file %s.\n",fname);
	return EXIT_FAILURE;
	}
close(fd);
} /* End of write_If_data */


/************************************************************
   Reads FLOATS from 3-D internal array into files.
************************************************************/
read_If_data(s,n,Re,Im)
char *s;
int n;
images Re,Im;
{
char fname[100];
int i,j,fd,read_re,read_im;

sprintf(fname,"%s%d",s,n);
printf("\nReading file: %s in read_If_data\n",fname);
fflush(stdout);
if((fd = open(fname,O_RDONLY,0644))>=0)
	{
        for(i=0;i<NUMFILES-2*OFFSET_T;i++) 
		{
		for(j=0;j<pic_x;j++)
		{
                if(read_re=read(fd,&Re[i][j][0],pic_y*sizeof(float))
			!=pic_y*sizeof(float))
			{
                        printf("Read error in read_If_data for Re: i=%d j=%d\n",i,j);
			printf("%d bytes read - %d bytes required\n",read_re,
				pic_y*sizeof(float));
			exit(1);
			}
		}
		for(j=0;j<pic_x;j++)
		{
                if(read_im=read(fd,&Im[i][j][0],pic_y*sizeof(float))
			!=pic_y*sizeof(float))
			{
                        printf("Read error in read_If_data for Im: i=%d j=%d\n",i,j);
			printf("%d bytes read - %d bytes required\n",read_im,
				pic_y*sizeof(float));
			exit(1);
			}
		}
		}
	}
	else
	{
	printf("Error opening file %s.\n",fname);
	exit(1);
	}
close(fd);
printf("File %s read in read_If_data\n",fname);
fflush(stdout);
} /* End of read_If_data */



/************************************************************
   Compute the complex 4 point demodulation difference kernel
************************************************************/
calc_complex_kernel(DEL_If_Re,DEL_If_Im,k0)
float DEL_If_Re[5],DEL_If_Im[5],k0;
{
DEL_If_Re[0] = -cos(2*k0)/12.0;
DEL_If_Re[1] = 8.0*cos(k0)/12.0;
DEL_If_Re[2] = 0.0;
DEL_If_Re[3] = -8.0*cos(k0)/12.0;
DEL_If_Re[4] = cos(2.0*k0)/12.0;

DEL_If_Im[0] = sin(2*k0)/12.0;
DEL_If_Im[1] = -8.0*sin(k0)/12.0;
DEL_If_Im[2] = 0.0;
DEL_If_Im[3] = -8.0*sin(k0)/12.0;
DEL_If_Im[4] = sin(2.0*k0)/12.0;
}


/************************************************************
   Apply 1D complex kernel in the x direction
************************************************************/
apply_kernel_x(Re,Im,Re_kernel,Im_kernel,x,y,t,Re_sum,Im_sum,kx)
images Re,Im;
float Re_kernel[5],Im_kernel[5],kx;
float *Re_sum,*Im_sum;
int x,y,t;
{
int i;

(*Re_sum) = 0.0;
(*Im_sum) = 0.0;
for(i=(-2);i<=2;i++)
	{
	(*Re_sum) += Re_kernel[i+2]*Re[t][x-i][y]-Im_kernel[i+2]*Im[t][x-i][y];
	(*Im_sum) += Im_kernel[i+2]*Re[t][x-i][y]+Re_kernel[i+2]*Im[t][x-i][y];
	}
(*Re_sum) = (*Re_sum) - kx*Im[t][x][y];
(*Im_sum) = (*Im_sum) + kx*Re[t][x][y];
}

/************************************************************
   Apply 1D complex kernel in the y direction
************************************************************/
apply_kernel_y(Re,Im,Re_kernel,Im_kernel,x,y,t,Re_sum,Im_sum,ky)
images Re,Im;
float Re_kernel[5],Im_kernel[5];
float *Re_sum,*Im_sum,ky;
int x,y,t;
{
int i;

(*Re_sum) = 0.0;
(*Im_sum) = 0.0;
for(i=(-2);i<=2;i++)
	{
	(*Re_sum) += Re_kernel[i+2]*Re[t][x][y-i]-Im_kernel[i+2]*Im[t][x][y-i];
	(*Im_sum) += Im_kernel[i+2]*Re[t][x][y-i]+Re_kernel[i+2]*Im[t][x][y-i];
	}
(*Re_sum) = (*Re_sum) - ky*Im[t][x][y];
(*Im_sum) = (*Im_sum) + ky*Re[t][x][y];
}


/************************************************************
   Apply 1D complex kernel in the t direction.
************************************************************/
apply_kernel_t(Re,Im,Re_kernel,Im_kernel,x,y,t,Re_sum,Im_sum,kt)
images Re,Im;
float Re_kernel[5],Im_kernel[5];
float *Re_sum,*Im_sum,kt;
int x,y,t;
{
int i;

(*Re_sum) = 0.0;
(*Im_sum) = 0.0;
for(i=(-2);i<=2;i++)
	{
	(*Re_sum) += Re_kernel[i+2]*Re[t-i][x][y]-Im_kernel[i+2]*Im[t-i][x][y];
	(*Im_sum) += Im_kernel[i+2]*Re[t-i][x][y]+Re_kernel[i+2]*Im[t-i][x][y];
	}
(*Re_sum) = (*Re_sum) - kt*Im[t][x][y];
(*Im_sum) = (*Im_sum) + kt*Re[t][x][y];
}


/************************************************************
   Computes Del Phi 
************************************************************/
Dphi(Re,Im,phi,kf,n)
images Re,Im;
float kf[3];
float phi[PIC_X][PIC_Y][3];
int n;
{
int i,j,k;
float DEL_If_Re_t[5],DEL_If_Im_t[5];
float DEL_If_Re_x[5],DEL_If_Im_x[5];
float DEL_If_Re_y[5],DEL_If_Im_y[5];
float Re_sum,Im_sum,amplitude2;
printf("In Dphi k0=%f k1=%f k2=%f\n",kf[0],kf[1],kf[2]);
fflush(stdout);

for(i=0;i<PIC_X;i++)
for(j=0;j<PIC_Y;j++)
for(k=0;k<3;k++)
	phi[i][j][k] = UNDEFINED;

calc_complex_kernel(DEL_If_Re_x,DEL_If_Im_x,kf[0]);
calc_complex_kernel(DEL_If_Re_y,DEL_If_Im_y,kf[1]);
calc_complex_kernel(DEL_If_Re_t,DEL_If_Im_t,kf[2]);

for(j=OFFSET_X;j<pic_x-OFFSET_X;j++)
for(k=OFFSET_Y;k<pic_y-OFFSET_Y;k++)
	{
	amplitude2=(Re[2][j][k]*Re[2][j][k]+Im[2][j][k]*Im[2][j][k]);
	if(amplitude2 != 0.0)
	{
	apply_kernel_x(Re,Im,DEL_If_Re_x,DEL_If_Im_x,j,k,2,&Re_sum,&Im_sum,kf[0]);
	phi[j][k][0] = (Im_sum*Re[2][j][k]-Re_sum*Im[2][j][k])/amplitude2;
	apply_kernel_y(Re,Im,DEL_If_Re_y,DEL_If_Im_y,j,k,2,&Re_sum,&Im_sum,kf[1]);
	phi[j][k][1] = (Im_sum*Re[2][j][k]-Re_sum*Im[2][j][k])/amplitude2;
	apply_kernel_t(Re,Im,DEL_If_Re_t,DEL_If_Im_t,j,k,2,&Re_sum,&Im_sum,kf[2]);
	phi[j][k][2] = (Im_sum*Re[2][j][k]-Re_sum*Im[2][j][k])/amplitude2;
	}
	else { phi[j][k][0] = phi[j][k][1] = phi[j][k][2] = UNDEFINED; }
	}
printf("Dphi has been computed for filter %d\n\n",n);
fflush(stdout);
}

/**********************************************************/
/* Initialize the threshold array			  */
/**********************************************************/
initialize_thresholds(thresholded)
unsigned char thresholded[22][PIC_X][PIC_Y];
{
int i,j,k;

for(i=0;i<22;i++)
{
for(j=0;j<PIC_X;j++)
for(k=0;k<PIC_Y;k++)
	thresholded[i][j][k] = THRESHOLDED;

for(j=OFFSET_X;j<pic_x-OFFSET_X;j++)
for(k=OFFSET_Y;k<pic_y-OFFSET_Y;k++)
	thresholded[i][j][k] = NOT_THRESHOLDED;
}
printf("Thresholds initialized to NOT_THRESHOLDED except for borders\n");
fflush(stdout);
}

/************************************************************
  Threshold the filter results based on a single 
  amplitude/frequency constraint
************************************************************/
threshold1(name,thresholded,maxamp,num_thresh,filter,percent_maxamp,tau)
char name[100];
unsigned char thresholded[22][PIC_X][PIC_Y];
float maxamp,filter[23][3],percent_maxamp,tau;
int num_thresh[22];
{
images Re,Im;
float ampthresh,Re_sum,Im_sum,dA[3],sigma,sigmak_tau,sigmak_tau_2;
float beta,base,f0,b,amplitude,amp2,maxA,minA,sum3;
int n,i,j,count[22];
float DEL_If_Re_t[5],DEL_If_Im_t[5];
float DEL_If_Re_x[5],DEL_If_Im_x[5];
float DEL_If_Re_y[5],DEL_If_Im_y[5];
float diff[3],phi[PIC_X][PIC_Y][3];

ampthresh = percent_maxamp*maxamp;
printf("Amplitude threshold: %f (%f percent of %f)\n",
	ampthresh,percent_maxamp*100.0,maxamp);

beta = 0.8;
base = 2.0;
b = pow(base,beta); /* 1.7411 */
maxA = -HUGE;
minA =  HUGE;

for(i=0;i<22;i++) num_thresh[i] = 0;
for(i=0;i<22;i++) count[i] = 0;

for(n=0;n<22;n++)
{
f0 = L2norm(&filter[n][0],3);
sigma = (b+1)/((b-1)*f0);
sigmak_tau = tau/sigma;
sigmak_tau_2 = sigmak_tau*sigmak_tau;

read_If_data(name,n,Re,Im);
Dphi(Re,Im,phi,&filter[n][0],n);

calc_complex_kernel(DEL_If_Re_x,DEL_If_Im_x,filter[n][0]);
calc_complex_kernel(DEL_If_Re_y,DEL_If_Im_y,filter[n][1]);
calc_complex_kernel(DEL_If_Re_t,DEL_If_Im_t,filter[n][2]);


for(i=OFFSET_X;i<pic_x-OFFSET_X;i++)
for(j=OFFSET_Y;j<pic_y-OFFSET_Y;j++)
{
if(thresholded[n][i][j] == NOT_THRESHOLDED)
	{
	amplitude=sqrt(Re[2][i][j]*Re[2][i][j]+Im[2][i][j]*Im[2][i][j]);
	amp2 = amplitude*amplitude;
	if(amplitude != 0.0)
		{
		apply_kernel_x(Re,Im,DEL_If_Re_x,DEL_If_Im_x,i,j,2,
			      &Re_sum,&Im_sum,filter[n][0]);
		dA[0] =  (Re[2][i][j]*Re_sum+Im[2][i][j]*Im_sum)/amplitude;
		apply_kernel_y(Re,Im,DEL_If_Re_y,DEL_If_Im_y,i,j,2,
		 	      &Re_sum,&Im_sum,filter[n][1]);
		dA[1] =  (Re[2][i][j]*Re_sum+Im[2][i][j]*Im_sum)/amplitude;
		apply_kernel_t(Re,Im,DEL_If_Re_t,DEL_If_Im_t,i,j,2,
		 	      &Re_sum,&Im_sum,filter[n][2]);
		dA[2] =  (Re[2][i][j]*Re_sum+Im[2][i][j]*Im_sum)/amplitude;

		diff[0] = (phi[i][j][0]-filter[n][0]);
		diff[0] = dA[0]*dA[0]/amp2 + diff[0]*diff[0];
		diff[1] = (phi[i][j][1]-filter[n][1]);
		diff[1] = dA[1]*dA[1]/amp2 + diff[1]*diff[1];
		diff[2] = (phi[i][j][2]-filter[n][2]);
		diff[2] = dA[2]*dA[2]/amp2 + diff[2]*diff[2];
		if(FALSE)
		{
		printf("\ni: %d j: %d sum3: %f\n",i,j,diff[0]+diff[1]+diff[2]);
		printf("diff: %f %f %f\n",diff[0],diff[1],diff[2]);
		printf("phi - k: %f %f %f\n",phi[i][j][0]-filter[n][0],
			phi[i][j][1]-filter[n][1],phi[i][j][2]-filter[n][2]);
		printf("dA: %f %f %f\n",dA[0],dA[1],dA[2]);
		printf("phi: %f %f %f\n",phi[i][j][0],phi[i][j][1],phi[i][j][2]);
		}
		}
	else 
		{
		diff[0] = diff[1] = diff[2] = HUGE;
		}

	sum3 = diff[0]+diff[1]+diff[2];
	if(sum3 > maxA) maxA = sum3;
	if(sum3 < minA) minA = sum3;
	if(amplitude < ampthresh || sum3 > sigmak_tau_2)
		{
		/* Threshold this filter response */
		thresholded[n][i][j] = THRESHOLDED;
		num_thresh[n]++;
		}
	else 
		{
		count[n]++;
		}
	}
}
}
printf("\nFilter responses thresholded\n");
printf("maxA: %f minA: %f\n",maxA,minA);
printf("Thresholding results:\n");
for(n=0;n<22;n++)
printf("Filter: %3d Unthresholded: %d\n",n,count[n]);
fflush(stdout);
}


/************************************************************
  Threshold the filter results based on amplitude and
  frequency constaints
  The old way of doing things, i.e. use two thresholds
************************************************************/
threshold2(name,thresholded,maxamp,num_thresh,filter,percent_maxamp,tauA,tauF)
char name[100];
unsigned char thresholded[22][PIC_X][PIC_Y];
float maxamp,filter[23][3],percent_maxamp,tauA,tauF;
int num_thresh[22];
{
images Re,Im;
float ampthresh,Re_sum,Im_sum,dA[3],sigma;
float beta,base,f0,b,amplitude,sigmak_tauF;
int n,i,j,threshA[22],threshF1[22],threshF2[22],already;
float DEL_If_Re_t[5],DEL_If_Im_t[5];
float DEL_If_Re_x[5],DEL_If_Im_x[5];
float DEL_If_Re_y[5],DEL_If_Im_y[5];
float diff[3],phi[PIC_X][PIC_Y][3];

ampthresh = percent_maxamp*maxamp;
printf("Amplitude threshold: %f (%f percent of %f)\n",
	ampthresh,percent_maxamp*100.0,maxamp);


beta = 0.8;
base = 2.0;
b = pow(base,beta); /* 1.7411 */

for(n=0;n<22;n++) 
	{
	num_thresh[n] = 0;
	threshA[n] = 0;
	threshF1[n] = threshF2[n] = 0;
	}

for(n=0;n<22;n++)
{
f0 = L2norm(&filter[n][0],3);
sigma = (b+1)/((b-1)*f0);
sigmak_tauF = tauF/sigma;

read_If_data(name,n,Re,Im);
Dphi(Re,Im,phi,&filter[n][0],n);

calc_complex_kernel(DEL_If_Re_x,DEL_If_Im_x,filter[n][0]);
calc_complex_kernel(DEL_If_Re_y,DEL_If_Im_y,filter[n][1]);
calc_complex_kernel(DEL_If_Re_t,DEL_If_Im_t,filter[n][2]);


for(i=OFFSET_X;i<pic_x-OFFSET_X;i++)
for(j=OFFSET_Y;j<pic_y-OFFSET_Y;j++)
{
if(thresholded[n][i][j] == NOT_THRESHOLDED)
	{
	amplitude=sqrt(Re[2][i][j]*Re[2][i][j]+Im[2][i][j]*Im[2][i][j]);
	if(amplitude != 0.0)
		{
		apply_kernel_x(Re,Im,DEL_If_Re_x,DEL_If_Im_x,i,j,2,&Re_sum,&Im_sum,filter[n][0]);
		dA[0] =  (Re[2][i][j]*Re_sum+Im[2][i][j]*Im_sum)/amplitude;
		apply_kernel_y(Re,Im,DEL_If_Re_y,DEL_If_Im_y,i,j,2,&Re_sum,&Im_sum,filter[n][1]);
		dA[1] =  (Re[2][i][j]*Re_sum+Im[2][i][j]*Im_sum)/amplitude;
		apply_kernel_t(Re,Im,DEL_If_Re_t,DEL_If_Im_t,i,j,2,&Re_sum,&Im_sum,filter[n][2]);
		dA[2] =  (Re[2][i][j]*Re_sum+Im[2][i][j]*Im_sum)/amplitude;

		diff[0] = (phi[i][j][0]-filter[n][0]);
		diff[1] = (phi[i][j][1]-filter[n][1]);
		diff[2] = (phi[i][j][2]-filter[n][2]);
		}
	else 
		{
		dA[0] = dA[1] = dA[2] = HUGE;
		diff[0] = diff[1] = diff[2] = HUGE;
		}

	/* Threshold the old way */
	if((amplitude < ampthresh) || (L2norm(&diff[0],3) >= sigmak_tauF) ||
	    ((amplitude != 0.0) && (sigma*L2norm(&dA[0],3))/amplitude > tauA))
		{
		/* Threshold this filter response */
		thresholded[n][i][j] = THRESHOLDED;
		num_thresh[n]++;
		}
	/* Accumulate thresholding statistics */
	already = FALSE;
	if((amplitude < ampthresh) ||
	   (sigma*L2norm(&dA[0],3))/amplitude > tauA) 
		{
		threshA[n]++;
		already = TRUE;
		}
	if(L2norm(diff,3) >= sigmak_tauF)
		{
		if(already==FALSE) threshF1[n]++;
		else threshF2[n]++;
		}
	}
}
}

if(FALSE)
{
printf("\nFilter responses thresholded\n");
printf("\n\nAmplitude Threshold Report\n");
printf("  Filter   Number Thresholded\n");
for(n=0;n<22;n++) printf("%10d   %10d\n",n,threshA[n]);
printf("\n\nFrequency Threshold Report\n");
printf("\nAdditional filter responses thresholded\n");
printf("  Filter   Number Thresholded\n");
for(n=0;n<22;n++) printf("%10d   %10d\n",n,threshF1[n]);
printf("\n\nFrequency Threshold Report\n");
printf("\nWould have been frequency thresholded if not already amplitude thresholded\n");
printf("  Filter   Number Thresholded\n");
for(n=0;n<22;n++) printf("%10d   %10d\n",n,threshF2[n]);
printf("\n\nOverall Threshold Report\n");
printf("  Filter   Number Thresholded\n");
for(n=0;n<22;n++) printf("%10d   %10d\n",n,num_thresh[n]);
fflush(stdout);
}
}


/*****************************************************************/
/* Initialize the filters and tuning information for each filter */
/*****************************************************************/
initfilters(filter,tuning_info,sigma)
float3 filter[];
float tuning_info[22][2],sigma;
{

float constant1,constant2,beta,base;
float wavelength,b,mu,speed,lambdas,lambdat,k1,k2,k3,rho,theta,two_pi;
int ii;

two_pi = 2.0*M_PI;
mu = 1.0;
beta = 0.8;
base = 2.0;
b = pow(base,beta); /* 1.7411 */

constant2 = sqrt(3.0);
constant1 = 1.0/constant2;


tuning_info[0][0] = 0.0; tuning_info[0][1] = 0.0;
tuning_info[1][0] = 30.0; tuning_info[1][1] = 0.0;
tuning_info[2][0] = 60.0; tuning_info[2][1] = 0.0;
tuning_info[3][0] = 90.0; tuning_info[3][1] = 0.0;
tuning_info[4][0] = 120.0; tuning_info[4][1] = 0.0;
tuning_info[5][0] = 150.0; tuning_info[5][1] = 0.0;
tuning_info[6][0] = 0.0; tuning_info[6][1] = constant1;
tuning_info[7][0] = 36.0; tuning_info[7][1] = constant1;
tuning_info[8][0] = 72.0; tuning_info[8][1] = constant1;
tuning_info[9][0] = 108.0; tuning_info[9][1] = constant1;
tuning_info[10][0] = 144.0; tuning_info[10][1] = constant1;
tuning_info[11][0] = 180.0; tuning_info[11][1] = constant1;
tuning_info[12][0] = 216.0; tuning_info[12][1] = constant1;
tuning_info[13][0] = 252.0; tuning_info[13][1] = constant1;
tuning_info[14][0] = 288.0; tuning_info[14][1] = constant1;
tuning_info[15][0] = 324.0; tuning_info[15][1] = constant1;
tuning_info[16][0] = 0.0; tuning_info[16][1] = constant2;
tuning_info[17][0] = 60.0; tuning_info[17][1] = constant2;
tuning_info[18][0] = 120.0; tuning_info[18][1] = constant2;
tuning_info[19][0] = 180.0; tuning_info[19][1] = constant2;
tuning_info[20][0] = 240.0; tuning_info[20][1] = constant2;
tuning_info[21][0] = 300.0; tuning_info[21][1] = constant2;
 
 
printf("\nFilter frequency information:\n");
for(ii=0;ii<22;ii++)
	{
        theta = tuning_info[ii][0];
        speed = tuning_info[ii][1];
 
        wavelength = (two_pi*(b-1)*sigma)/(mu*(b+1));
        lambdas = wavelength*(sqrt(speed*speed+1.0));
        if(fabs(speed) > 0.000001) lambdat = -lambdas/speed;
        else lambdat = HUGE;
 
        k3 = two_pi/lambdat;
        rho = two_pi/lambdas;
        k1 = rho*sin(theta*two_pi/360.0);
        k2 = -rho*cos(theta*two_pi/360.0);
        filter[ii][0] = -k2;
        filter[ii][1] = k1;
        filter[ii][2] = k3;
	printf("filter %2d   %10.6f  %10.6f  %10.6f\n",ii,filter[ii][0],filter[ii][1],filter[ii][2]);
	}

printf("Filter and tuning information initialized\n");
if(FALSE)
{
theta = 108.0;
speed = 1.0/sqrt(3.0);
wavelength = (two_pi*(b-1)*sigma)/(mu*(b+1));
lambdas = wavelength*(sqrt(speed*speed+1.0));
if(fabs(speed) > 0.000001) lambdat = -lambdas/speed;
else lambdat = HUGE;
k3 = two_pi/lambdat;
rho = two_pi/lambdas;
k2 = rho*sin(theta*two_pi/360.0);
k1 = rho*cos(theta*two_pi/360.0);
printf("\nFrequencies of sinusoidal moving at 108 degrees, speed 1/sqrt(3):\n");
printf("\nSpeed: %f Direction: %f k1=%f k2=%f k3=%f\n",speed,theta,k1,k2,k3);
fflush(stdout);
}
} /* End of initfilters */

/******************************************************************/
/* Produce thresholding report 					  */
/******************************************************************/
produce_thesh_report(num_thresh,filter,tuning_info)
int num_thresh[22];
float filter[23][3],tuning_info[22][2];
{
int i;

printf("\n\nThreshold Report\n");
printf("\nNumber of filter responses amplitude thresholded:\n");
for(i=0;i<22;i++)
	{
	printf("Filter:%d Frequencies; %f %f %f\n",i,
		filter[i][0],filter[i][1],filter[i][2]);
	printf(" Angle: %f Speed: %f Number thresholded: %d Percentage thresholded: %f\n",
		tuning_info[i][0],tuning_info[i][1],num_thresh[i],
		((num_thresh[i]*1.0)/
		((pic_x-2*OFFSET_X)*(pic_y-2*OFFSET_Y)*1.0))*100.0);
	}
printf("\n\n\n");
fflush(stdout);
}




/************************************************************
   Returns the norm of a vector v of length n.
************************************************************/
float L2norm(v,n)
float v[];
int n;
{
int i;
float sum = 0.0;
for (i=0;i<n; i++) 
	sum += (float)(v[i]*v[i]);
sum = sqrt(sum);
return(sum);
}



/************************************************************
   Compute normal velocities for filter n output using
   Del Phi data.
************************************************************/
calc_normal_velocities(phi,normal_velocities,n,thresholded)
float phi[PIC_X][PIC_Y][3];
float normal_velocities[22][PIC_X][PIC_Y][2];
unsigned char thresholded[22][PIC_X][PIC_Y];
int n;
{
int i,j,count;
float denom;

for(i=0;i<PIC_X;i++)
for(j=0;j<PIC_Y;j++)
	{
	normal_velocities[n][i][j][0] = normal_velocities[n][i][j][1] = UNDEFINED;
	}

count = 0;
for(i=OFFSET_X;i<pic_x-OFFSET_X;i++)
for(j=OFFSET_Y;j<pic_y-OFFSET_Y;j++)
	if(thresholded[n][i][j] == NOT_THRESHOLDED)
		{
		if(phi[i][j][0] == UNDEFINED || phi[i][j][1] == UNDEFINED ||
		   phi[i][j][2] == UNDEFINED)
			{
			printf("\nFatal error: phi undefined but not thresholded\n");
			printf("phi: %f %f %f\n",phi[i][j][0],phi[i][j][1],phi[i][j][2]);
			printf("Re: %f Im: %f\n",IfRe[2][i][j],IfIm[2][i][j]);
			printf("n=%d i=%d j=%d\n",n,i,j);
			exit(1);
			}
		denom=(phi[i][j][0]*phi[i][j][0]+phi[i][j][1]*phi[i][j][1]);
		normal_velocities[n][i][j][0] = -(phi[i][j][0]*phi[i][j][2])/denom;
		normal_velocities[n][i][j][1] = -(phi[i][j][1]*phi[i][j][2])/denom;
		count ++;
		}
printf("%d normal velocities computed\n",count);
}

/************************************************************
   Compute all full velocities from 5*5 neighbourhoods of
   normal velocities.
************************************************************/
calc_full_velocities(normal_velocities,full_velocities,E)
float normal_velocities[22][PIC_X][PIC_Y][2],E[PIC_X][PIC_Y];
float full_velocities[PIC_X][PIC_Y][2];
{
int i,j,k,l,m,ii,count,ct_cond_yes,ct_cond_no,no_count;
float size,n1,n2,residual_error;
float difference[NO_NORM_VEL],fvn[NO_NORM_VEL];
double v[3][3][2],vn[NO_NORM_VEL];
double J[NO_NORM_VEL][NO_UNKNOWNS],product[NO_NORM_VEL];
double JI[NO_UNKNOWNS][NO_NORM_VEL],condnum;

for(i=0;i<PIC_X;i++)
for(j=0;j<PIC_Y;j++)
	{
	E[i][j] = full_velocities[i][j][0] = full_velocities[i][j][1] = NO_VALUE;
	normal_ct[i][j] = 0;
	}

ct_cond_yes = ct_cond_no = no_count = 0;;

for(i=(2*OFFSET_X+NRADIUS);i<pic_x-(2*OFFSET_X+NRADIUS);i++)
for(j=(2*OFFSET_Y+NRADIUS);j<pic_y-(2*OFFSET_Y+NRADIUS);j++)
	{
	count = 0;
	for(k=0;k<22;k++)
	for(l=i-2;l<=i+2;l++)
	for(m=j-2;m<=j+2;m++)
		{
		if(thresholded[k][l][m] == NOT_THRESHOLDED)
			{
			size = L2norm(&normal_velocities[k][l][m][0],2);
			if(size>0.000001)
			{
			n1 = normal_velocities[k][l][m][0]/size;
			n2 = normal_velocities[k][l][m][1]/size;
			J[count][0] = n1;
			J[count][1] = n1*(l-i);
			J[count][2] = n1*(m-j);
			J[count][3] = n2;
			J[count][4] = n2*(l-i);
			J[count][5] = n2*(m-j);
			vn[count] = size;
			count++;
			}
			}
		}
	if(i % 10 == 0 && j==(pic_y-(2*OFFSET_Y+NRADIUS)-1))
		{
		printf("\nFull velocities for row %d have been processed\n",i);
		}

	normal_ct[i][j] = count;
	if(count >= 6) /* At least 6 normal velocities are needed */
		condnum = cal_velocity(count,J,JI,vn,v,product);
	else condnum = HUGE;

	for(ii=0;ii<count;ii++) 
		{
		difference[ii] = (float) (vn[ii] - product[ii]);
		difference[ii] = difference[ii]*difference[ii];
		fvn[ii] = (float) vn[ii];
		}

	if(condnum < 10.0 && 
	  (residual_error=L2norm(difference,count)/L2norm(fvn,count)) < 0.5)
		{
		E[i][j] = residual_error;
		full_velocities[i][j][0] = (float) v[1][1][0];
		full_velocities[i][j][1] = (float) v[1][1][1];
		ct_cond_yes++;
		}
	else if(count >= 6)
		{
		full_velocities[i][j][0] = NO_VALUE;
		full_velocities[i][j][1] = NO_VALUE;
		ct_cond_no++;
		}
	else no_count++;
	}
printf("\nNumber of image locations where full velocity is computed: %d Percentage: %f\n", ct_cond_yes,ct_cond_yes*100.0/(ct_cond_yes+ct_cond_no+no_count));
printf("Number of image locations where full velocity is not computed: %d Percentage: %f\n", (ct_cond_no+no_count),(ct_cond_no+no_count)*100.0/(ct_cond_yes+ct_cond_no+no_count));
printf("Number of image positions where image velocity calculation was attempted but failed: %d\n",ct_cond_no);
printf("Number of image positions where image velocity calculation was not attempted: %d\n",no_count);
}


/************************************************************
   Output full velocities using Burkitt format
************************************************************/
output_full_velocities(fdf,full_velocities)
float full_velocities[PIC_X][PIC_Y][2];
int fdf;
{
float x,y;
int i,j,bytes,no_novals,no_vals;

if(fdf<0)
	{
	printf("\nFatal error: full velocity file not opened\n");
	exit(1);
	}
/* original size */
y = pic_x;
x = pic_y;
write(fdf,&x,4);
write(fdf,&y,4);

/* size of result data */
y = pic_x-(4*OFFSET_X+2*NRADIUS);
x = pic_y-(4*OFFSET_Y+2*NRADIUS);
write(fdf,&x,4);
write(fdf,&y,4);

/* offset to start of data */
y = 2*OFFSET_X+NRADIUS;
x = 2*OFFSET_Y+NRADIUS;
write(fdf,&x,4);
write(fdf,&y,4);
bytes = 24;

no_novals = no_vals = 0;
/* Prepare velocities for output, i.e. rotate by 90 degrees */
for(i=(2*OFFSET_X+NRADIUS);i<pic_x-(2*OFFSET_X+NRADIUS);i++)
for(j=(2*OFFSET_Y+NRADIUS);j<pic_y-(2*OFFSET_Y+NRADIUS);j++)
	{
	x = -full_velocities[i][j][0];
	y = full_velocities[i][j][1];
	if(full_velocities[i][j][0] != NO_VALUE && 
	   full_velocities[i][j][1] != NO_VALUE)
		{
		full_velocities[i][j][0] = y;
		full_velocities[i][j][1] = x;
		no_vals++;
		}
	else
		{
		no_novals++;
		}
	}
for(i=(2*OFFSET_X+NRADIUS);i<pic_x-(2*OFFSET_X+NRADIUS);i++)
	{
	bytes += write(fdf,&full_velocities[i][2*OFFSET_Y+NRADIUS][0],(pic_y-(4*OFFSET_Y+2*NRADIUS))*8);
	}
close(fdf);
printf("\nFull velocities output: %d bytes\n",bytes);
printf("Number of positions with velocity: %d\n",no_vals);
printf("Number of positions without velocity: %d\n",no_novals);
printf("Percentage of full velocities: %f\n",
	no_vals/(1.0*(no_vals+no_novals))*100.0);
fflush(stdout);
}


/************************************************************
   Output normal velocities
************************************************************/
output_normal_velocities(fdn,normal_velocities)
float normal_velocities[22][PIC_X][PIC_Y][2];
int fdn;
{
int i,j,k,num,no_bytes,no_points,no_normal_vels;
float x,y,last,data[5000];

no_points = 0;
no_normal_vels = 0;
last = END_OF_VALUES;
if(fdn<0)
	{
	printf("\nFatal error: normal velocity file not opened\n");
	exit(1);
	}
/* original size */
y = pic_x;
x = pic_y;
write(fdn,&x,4);
write(fdn,&y,4);

/* size of result data */
y = pic_x-(4*OFFSET_X+2*NRADIUS);
x = pic_y-(4*OFFSET_Y+2*NRADIUS);
write(fdn,&x,4);
write(fdn,&y,4);

/* offset to start of data */
y = 2*OFFSET_X+NRADIUS;
x = 2*OFFSET_Y+NRADIUS;
write(fdn,&x,4);
write(fdn,&y,4);

no_bytes = 24;

for(i=(2*OFFSET_X+NRADIUS);i<pic_x-(2*OFFSET_X+NRADIUS);i++)
for(j=(2*OFFSET_Y+NRADIUS);j<pic_y-(2*OFFSET_Y+NRADIUS);j++)
	{
	x = j;
	y = i;
	num = 2;
	/* write(fdn,&x,4);
	write(fdn,&y,4); */
	data[0] = x; data[1] = y;
	no_points++;
	for(k=0;k<22;k++)
		{
		if(thresholded[k][i][j] == NOT_THRESHOLDED)
			{
			x = -normal_velocities[k][i][j][0];
			y = normal_velocities[k][i][j][1]; 
			/* write(fdn,&y,4);
			write(fdn,&x,4); */
			data[num] = y;
			num++;
			data[num] = x;
			num++;
			no_normal_vels++;
			}
		}
	/* write(fdn,&last,4); */
	data[num] = last;
	num++;
	if(num > 5000) 
		{
		printf("Fatal error: data array not big enough for normal velocity data\n");
		exit(1);
		}
	no_bytes += write(fdn,data,num*4);
	}
close(fdn);
printf("\nNormal velocities output: %d bytes\n",no_bytes);
printf("%d normal velocities at %d points\n",no_normal_vels,no_points);
fflush(stdout);
}



/******************************************************************/
/* This function accepts a collection of r normal image velocity  */
/* measurements and a r*6 matrix computed using normal directions */
/* and image locations and computes the full image velocity for   */
/* some local neighbourhood. A linear approximation is used to    */
/* relate normal image velocity in some small neighbourhood to    */
/* full image velocity at one point.				  */
/******************************************************************/
double cal_velocity(r,J,JI,vn,v,product)
double v[3][3][2],vn[NO_NORM_VEL];
double J[NO_NORM_VEL][NO_UNKNOWNS],product[NO_NORM_VEL];
double JI[NO_UNKNOWNS][NO_NORM_VEL];
int r;
{
int i,j;
double condnum,alphabeta[NO_UNKNOWNS];

condnum=cal_inverse(r,J,JI);
if(condnum != HUGE)
{
for(i=0;i<NO_UNKNOWNS;i++)
	{
	alphabeta[i] = 0.0;
	for(j=0;j<r;j++)
		alphabeta[i] = alphabeta[i] + JI[i][j]*vn[j];
	}

for(i=0;i<r;i++)
	{
	product[i] = 0.0;
	for(j=0;j<NO_UNKNOWNS;j++)
		product[i] = product[i] + J[i][j]*alphabeta[j];
	}
/* Compute full velocity field in 3*3 neighbourhood */
for(i=(-1);i<=1;i++)
for(j=(-1);j<=1;j++)
	{
	v[i+1][j+1][0] = alphabeta[0]+alphabeta[1]*(i)+alphabeta[2]*(j);
	v[i+1][j+1][1] = alphabeta[3]+alphabeta[4]*(i)+alphabeta[5]*(j);
	}

}
return(condnum);
}

/*************************************************************************/
/* Compute the pseudo-inverse of J using its SVD                         */
/*************************************************************************/
double cal_inverse(r,J,JI)
double J[NO_NORM_VEL][NO_UNKNOWNS];
double JI[NO_UNKNOWNS][NO_NORM_VEL];
int r;
{
int size;
double 	VT[NO_UNKNOWNS][NO_UNKNOWNS],
	VV[NO_UNKNOWNS][NO_UNKNOWNS],
	U[NO_UNKNOWNS][NO_NORM_VEL],
	DI[NO_UNKNOWNS][NO_UNKNOWNS],
	UT[NO_NORM_VEL][NO_UNKNOWNS],
	V[NO_UNKNOWNS][NO_UNKNOWNS],
	D[NO_UNKNOWNS][NO_UNKNOWNS],
	JT[NO_UNKNOWNS][NO_NORM_VEL],
	UU[NO_UNKNOWNS][NO_NORM_VEL],
	I[NO_UNKNOWNS][NO_UNKNOWNS],
	W[NO_UNKNOWNS+1],temp[NO_NORM_VEL],
	zero[NO_UNKNOWNS],condnum,min,max;

int i,j,k,m,n,mdim,error,Ierr,job;

size = r;
n = NO_UNKNOWNS; m = r; mdim = NO_NORM_VEL;
for(i=0;i<NO_UNKNOWNS;i++)
for(j=0;j<size;j++)
	JT[i][j] = J[j][i];

job = 21;

/* Call limpack double percision SVD routine */
dsvdc(JT,&mdim,&m,&n,W,zero,UU,&mdim,VV,&n,temp,&job,&Ierr);  

/* Undo "Fortran" damage, i.e. transpose the data and convert to
   single percision.  */
for(i=0;i<size;i++)
for(j=0;j<NO_UNKNOWNS;j++)
	{
	U[j][i] = UU[j][i];
	UT[i][j] = UU[j][i];
	}
for(i=0;i<NO_UNKNOWNS;i++)
for(j=0;j<NO_UNKNOWNS;j++)
	{
	V[i][j] = VV[i][j];
	D[i][j] = DI[i][j] = 0.0;
	}
for(i=0;i<NO_UNKNOWNS;i++)
for(j=0;j<NO_UNKNOWNS;j++)
	VT[i][j] = V[j][i];
/* Compute condition number as the ratio of maximum to
   minimum diagonal element of the SVD diagonal matrix */
min = HUGE;
max = 0.0;
for(i=0;i<NO_UNKNOWNS;i++)
	{
	D[i][i] = W[i];
	if(fabs(W[i]) > max) max = fabs(W[i]);
	if(fabs(W[i]) < min) min = fabs(W[i]);
	}
if(min != 0.0) condnum = max/min;
else condnum=HUGE;
/* Compute the inverse of the diagonal matrix */
for(i=0;i<NO_UNKNOWNS;i++) if(D[i][i] !=0.0) DI[i][i] = 1.0/D[i][i];
	else DI[i][i] = HUGE;

/* Check correctness of SVD, compute I = JI*J */
for(i=0;i<NO_UNKNOWNS;i++)
for(j=0;j<NO_UNKNOWNS;j++)
	{
	VV[i][j] = 0.0;
	for(k=0;k<NO_UNKNOWNS;k++) VV[i][j] = VV[i][j] + VT[i][k]*DI[k][j];
	}
for(i=0;i<NO_UNKNOWNS;i++)
for(j=0;j<size;j++)
	{
	JI[i][j] = 0.0;
	for(k=0;k<NO_UNKNOWNS;k++) JI[i][j] = JI[i][j] + VV[i][k]*U[k][j];
	}

for(i=0;i<NO_UNKNOWNS;i++)
for(j=0;j<NO_UNKNOWNS;j++)
	{
	I[i][j] = 0.0;
	for(k=0;k<size;k++)
		I[i][j] = I[i][j] + JI[i][k]*J[k][j];
	}

/* Check I */
if (condnum < 100000.0) 
{
error = FALSE;
for(i=0;i<NO_UNKNOWNS;i++)
for(j=0;j<NO_UNKNOWNS;j++)
	{
	if(i==j) 
		{
		if(fabs(I[i][i])<0.9999 || fabs(I[i][i])>1.0001) error=TRUE;
		}
	else if(fabs(I[i][j]) > 0.0001) error=TRUE;
	}
if(error==TRUE)
	{
	printf("\nFatal error: Pseudo-inverse of Jacobian wrong\n");
	printf("The identity matrix:\n");
	for(i=0;i<NO_UNKNOWNS;i++)
		{
		for(j=0;j<NO_UNKNOWNS;j++) printf("%8.3f ",I[i][j]);
		printf("\n");
		}
	fflush(stdout);
	/*exit(1); */
	}
}
return(condnum);
}



/************************************************************************/
/* Compute average angle, standard deviation, density (as a percentage) */
/* and average residual error density as a percentage	     		*/
/************************************************************************/
calc_statistics(correct_vels,int_size_x,int_size_y,full_vels,E,pic_x,pic_y,n,
		ave_error,st_dev,density,residual,min_angle,max_angle,
		norm_vels,norm_err,norm_st_dev,min_norm_angle,max_norm_angle)
float full_vels[PIC_X][PIC_Y][2],*ave_error,*density,*st_dev,*residual;
float correct_vels[PIC_X][PIC_Y][2],E[PIC_X][PIC_Y],*min_angle,*max_angle;
int n,pic_x,pic_y,int_size_x,int_size_y;
float norm_vels[22][PIC_X][PIC_Y][2];
float *norm_err,*norm_st_dev,*min_norm_angle,*max_norm_angle;
{
int count,no_count,i,j,a,b,nn,norm_count;
float sumX2,temp,uva[2],uve[2],sum2;

count = no_count = norm_count = 0;
sum2 = sumX2 = 0.0;
(*min_angle) = HUGE;
(*max_angle) = 0.0;
(*ave_error) = (*st_dev) = 0.0;
(*density) = (*residual) = 0.0;
(*min_norm_angle) = HUGE;
(*max_norm_angle) = -HUGE;
(*norm_err) = (*norm_st_dev) = 0.0;

for(i=(2*OFFSET_X+NRADIUS);i<pic_x-(2*OFFSET_X+NRADIUS);i++)
{
for(j=(2*OFFSET_Y+NRADIUS);j<pic_y-(2*OFFSET_Y+NRADIUS);j++)
	{
	if(full_vels[i][j][0] == NO_VALUE && full_vels[i][j][1] == NO_VALUE)
		no_count++;
	else 
	  {
	  count++;
	  uve[0] = full_vels[i][j][0]; uve[1] = full_vels[i][j][1];
	  uva[0] = correct_vels[i][j][0]; uva[1] = correct_vels[i][j][1];
	  temp = PsiER(uve,uva);
	  (*ave_error) += temp;
	  if(E[i][j] == NO_VALUE) 
		{ printf("Fatal error: E has no value\n"); exit(1); }
	  (*residual) += E[i][j];
	  sumX2 += temp*temp;
	  if(temp < (*min_angle)) (*min_angle) = temp;
	  if(temp > (*max_angle)) (*max_angle) = temp;
	  }
	for(nn=0;nn<22;nn++)
	{
	if(norm_vels[nn][i][j][0]!=NO_VALUE && norm_vels[nn][i][j][1]!=NO_VALUE)
	  {
	  norm_count++;
	  uve[0] = norm_vels[nn][i][j][1]; uve[1] = -norm_vels[nn][i][j][0];
	  uva[0] = correct_vels[i][j][0]; uva[1] = correct_vels[i][j][1];
	  temp = PsiEN(uve,uva);
	  (*norm_err) += temp;
	  sum2 += temp*temp;
	  if(temp < (*min_norm_angle)) (*min_norm_angle) = temp;
	  if(temp > (*max_norm_angle)) (*max_norm_angle) = temp;
	  }
	}
	}
}
(*density) = (count*100.0)/(count+no_count);
if(count != 0) (*ave_error) = (*ave_error)/(count*1.0);
if(count > 1) (*st_dev) = sqrt((sumX2 - count*(*ave_error)*(*ave_error))/(count-1));
else (*st_dev) = 0.0;
if(norm_count != 0) (*norm_err) = (*norm_err)/(norm_count*1.0);
if(count > 1) (*norm_st_dev) = sqrt((sum2 - norm_count*(*norm_err)*(*norm_err))/(norm_count-1));
else (*norm_st_dev) = 0.0;
if(count != 0) (*residual) = (*residual)/count;
if((*ave_error) == 0.0) { (*min_angle) = (*max_angle) = 0.0; }
if((*norm_err) == 0.0) { (*min_norm_angle) = (*max_norm_angle) = 0.0; }
fflush(stdout);
}

/************************************************************
 Full Image Velocity Angle Error
************************************************************/
float PsiER(ve,va)
float ve[2],va[2];
{
float nva;
float nve;
float v,r,temp;
float VE[3],VA[3];

VE[0] = ve[0];
VE[1] = ve[1];
VE[2] = 1.0;

VA[0] = va[0];
VA[1] = va[1];
VA[2] = 1.0;

nva = L2norm(VA,3);
nve = L2norm(VE,3);
v =  (VE[0]*VA[0]+VE[1]*VA[1]+1.0)/(nva*nve);

/**  sometimes roundoff error causes problems **/
if(v>1.0 && v < 1.0001) v = 1.0;

r = acos(v)*180.0/PI;

if (!(r>=0.0 && r< 180.0))
{
printf("ERROR in PSIER()...\n r=%8.4f  v=%8.4f  nva=%8.4f nve= %8.4f\n",r,v,nva,nve);
printf("va=(%f,%f) ve=(%f,%f)\n",va[0],va[1],ve[0],ve[1]);
}

return r;
}

/************************************************************
 Normal Image Velocity Angle Error
************************************************************/
float PsiEN(ve,va)
float ve[2],va[2];
{
float nva,nve;
float v1,v2,term1,term2;
float n[2];

nva = L2norm(va,2), nve = L2norm(ve,2);
if(nve > 0.00000001)
{
n[0] = ve[0]/nve;
n[1] = ve[1]/nve;
v1 = (va[0]*n[0] + va[1]*n[1]-nve) ;
v2 = v1/(sqrt((1+nva*nva))*sqrt((1+nve*nve)));
v1 =  asin(v2)*180.0/PI;

if(!(v1>=-90.0 && v1<=90.0))
	{
       	printf("ERROR in PSIEN()  v1: %f ve: %f\n",v1,v2);
	printf("nve: %f nva: %f\n",nve,nva);
	printf("n: %f %f\n",n[0],n[1]);
	printf(" ve: %f %f va: %f %f\n",ve[0],ve[1],va[0],va[1]);
	fflush(stdout);
	}
}
else v1 = NO_VALUE;
	
return(v1);
}

