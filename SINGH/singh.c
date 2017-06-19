/**************************************************************
   singh.c 
   -- Written by John Barron, 1992
   -- Implementation of Ajit Singh, ICCV, 1990, pp168-177.
   -- See also "Optic Flow Computation: A Unified Perspective"
      IEEE Computer Society Press 1992.
**************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/types.h>


#define PIC_X  325
#define PIC_Y  325
#define PMODE 0644
#define TRUE 1
#define FALSE 0
#define HEAD 32
#define TORAD (M_PI/180.)
#define TODEG (180./M_PI)
#define NO_VALUE 100.0
#define ARRSIZE 9
#define MAX_ARRSIZE 17
#define WEIGHTSIZE 21
#define NUMFILES 3
#define FULL 69
#define NORMAL 96
#define NOTHING 0
#define CELL_RASTER 0
#define WHITE 255
#define BLACK 0
#define PRINT 0
#define PRINT2 0
#define NO_ROWS 2
#define NO_COLS 2
#define DIM 2
#define DIFF_THRESH 1e-5
#define SVD_PRINT 0
#define TOTAL 64
#define PT_X1 75
#define PT_Y1 60
#define PT_X2 60
#define PT_Y2 75
#define SMALL 0.1
int EXTRA_OFFSET_X=0;
int EXTRA_OFFSET_Y=0;
int CUT_X=0,CUT_Y=0;
int PRINT_ITERATIONS=0;
int PRINT_LAPLACIAN=0;
int PRINT1=0;
int PRINT11=0;
int OFFSET_X=10;
int OFFSET_Y=10;



typedef unsigned char cimages[NUMFILES][PIC_Y][PIC_X];
typedef float fimages[NUMFILES][PIC_Y][PIC_X];

float cal_normal_angle_error(),cal_full_angle_error();
float dot(),PsiER(),PsiEN(),L2norm(),bigL2norm();
int check_eigen_calc();


float calc[ARRSIZE][ARRSIZE];
int fd_correct,SAMPLE,PREVIOUS1,PREVIOUS2,pic_x,pic_y,no_bytes,fd_vels;
float offset_x,offset_y,size_x,size_y,actual_x,actual_y;
int int_size_x,int_size_y,inverse22(),fdcov1,fdcov2;
int COUNT_SINGULAR,STEP,PRINT_FLOWS1,PRINT_FLOWS2;
int MAX_NUMBER_OF_ITERATIONS,SUBSET;
long time1,time2,time3,time4,temporal_offset;
float normal_velocities[PIC_Y][PIC_X][2];
float full_velocities[PIC_Y][PIC_X][2],threshold_velocity[PIC_Y][PIC_X][2];
float Ucc[PIC_Y][PIC_X][2],Scc[PIC_Y][PIC_X][2][2],S[PIC_Y][PIC_X][2][2];
float correct_velocities[PIC_Y][PIC_X][2],TAU,SsumI[PIC_Y][PIC_X][2][2];
float Un[2][PIC_Y][PIC_X][2],Sn[2][PIC_Y][PIC_X][2][2];
float SccI[PIC_Y][PIC_X][2][2],Ssum[PIC_Y][PIC_X][2][2];
float SccI_Ucc[PIC_Y][PIC_X][2],Ua[PIC_Y][PIC_X][2];
int fdf1,fdf2,BINARY,CORRECT_VELOCITIES,header_ints[8],STEP2,FILE_NUMBER;
int LAPLACIAN;
FILE *fp;

cimages inpic;
fimages fpic;
char inname[100],filename[100],path1[100],path2[100],outname[100];


/******************************************************************/
int main(int argc, char **argv)
{
char command[100],correct_filename[100],name[100];
float tau1,tau2,sigma;
int n,N,i,central_image,old_pic_x,old_pic_y,steps1,steps2;
int threshold1,threshold2,w,offset;
float ave_error,st_dev,density,min_angle,max_angle;
float low_tau1,high_tau1,inc_tau1;
float low_tau2,high_tau2,inc_tau2;

FILE_NUMBER = -1;
COUNT_SINGULAR = 0;
SUBSET = FALSE;

if(argc<7 || argc>36)
	{
        printf("Usage: %s <filename stem> <central_image> <input data path> <output data path> [-C <filename> -B <cols> <rows> -T1 <low> <high> <steps> -T2 <low> <high> <steps>  -P1 -P2 -PR1 -PR2 -PRI -PRL -n <int> -N <int> -w <int> -i <int> -NL -s <int> -SUB <int> <int> <int> <int>]\n",argv[0]);
	printf("Use 2 images starting with stem.center_image\n");
	printf("<input data path> - directory where images are\n");
	printf("<output data path> - directory where the flow fields are stored\n");
	printf("-C <file> - correct velocity data provided and error analysis\n");
	printf("            will be performed on it\n");
	printf("-B <cols> <rows> - Input binary images-specify dimensions \n");
	printf("                   not necessary for rasterfiles\n");
	printf("-T1 <float> <float> <int> - threshold the step 1 velocities using the smallest\n");
	printf("              eigenvalue of the step 1covariance matrix at each position\n");
	printf("-T2 <float> <float> <int> - threshold the step 2 velocities using the smallest\n");
	printf("              eigenvalue of the step 2 covariance matrix at each position\n");
	printf("-P1 - do not perform step 1 of the compuatation, instead read the\n");
	printf("      previously computed step 1 velocities and covariance matrices\n");
	printf("-P2 - do not perform step 2 of the compuatation, instead read the\n");
	printf("      previously computed step 2 velocities and covariance matrices\n");
	printf("    If -P2 is used but -P1 is not the program terminates after step 1\n");
	printf("-PR1 - print thresholded flow fields for step 1\n");
	printf("-PR2 - print thresholded flow fields for step 2\n");
	printf("-PRI - print flow fields after each iteration\n");
	printf("     - These flow fields occupy a LOT of ddisk space\n");
	printf("-PRL - print the laplacian images\n");
	printf("-n - neighbourhood size (2n+1) in step1 computation (default 2)\n");
	printf("-N - maximum displacement in pixels (-u,-v <= N <= u,v) (default and maximum value is/can be 4)\n");
	printf("-w - window size (2w+1) for step 2 computation [currently must be 1 or 2] (default is 1)\n");
	printf("-i - number of iterations for step 2 (default 10)\n");
	printf("-s <int> - specify the distance of the adjacent left and right frames\n");
	printf("           from the central frame - default is 1\n");
	printf("-NL - use the original images as input\n");
	printf("-SUB <x> <y> <size_x> <size_y> - compute velocity for a subset of the images\n");
	printf("     compute flow for a subarea starting at (x,y) and of size (size_x,size_y)\n");
	printf("     Note: you must take offsets into account!\n");
	printf("-C, -B, -T1, -T2, -P1, -P2, -n, -N, -w, -i, -PR1, -PR2, PRI, -PRL, -NL and -SUB can be in any order and are optional\n");
	exit(1);
        }
else
	{
	printf("%d arguments\n",argc);
	printf("Command line:");
	for(i=0;i<argc;i++) printf("%s ",argv[i]);
	printf("\n");
	}


sscanf(argv[2],"%d",&central_image);
strcpy(path1,argv[3]); 
strcpy(path2,argv[4]); 
printf("<input path>: %s\n<output path>: %s\n",path1,path2);

printf("SMALL: %f\n",SMALL);

BINARY = FALSE;
CORRECT_VELOCITIES = FALSE;
SAMPLE = 1;
STEP = 1;
PREVIOUS1 = FALSE;
PREVIOUS2 = FALSE;
LAPLACIAN = TRUE;
PRINT_FLOWS1 = FALSE;
PRINT_FLOWS2 = FALSE;
PRINT_ITERATIONS = FALSE;
PRINT_LAPLACIAN = FALSE;
strcpy(correct_filename,"unknown"); 
i = 6; 
threshold1 = threshold2 = FALSE;
STEP2 = TRUE;
w = 1;
n = 2;
N = 4;
MAX_NUMBER_OF_ITERATIONS = 10;
while(i <= argc) 
	{
	if(strcmp(argv[i-1],"-S1")==0) 
		{
		sscanf(argv[i],"%d",&SAMPLE);
		STEP2 = FALSE;
		i++;
		}
	else
	if(strcmp(argv[i-1],"-T1")==0) 
		{
		sscanf(argv[i],"%f",&low_tau1);
		sscanf(argv[i+1],"%f",&high_tau1);
		sscanf(argv[i+2],"%d",&steps1);
		threshold1 = TRUE;
		if((low_tau1 >= high_tau1) || steps1 <= 0)
			{
			printf("Error: threshold ranges for step1 are wrong - no thresholding done\n");
			threshold1 = FALSE;
			}
		else
			{
			inc_tau1 = (high_tau1-low_tau1)/steps1;
			}
		i+=3;
		}
	else
	if(strcmp(argv[i-1],"-T2")==0) 
		{
		sscanf(argv[i],"%f",&low_tau2);
		sscanf(argv[i+1],"%f",&high_tau2);
		sscanf(argv[i+2],"%d",&steps2);
		threshold2 = TRUE;
		if((low_tau2 >= high_tau2) || steps2 <= 0)
			{
			printf("Error: threshold ranges for step2 are wrong - no thresholding done\n");
			threshold2 = FALSE;
			}
		else
			{
			inc_tau2 = (high_tau2-low_tau2)/steps2;
			}
		i+=3;
		}
	else if(strcmp(argv[i-1],"-C")==0)
		{
		strcpy(correct_filename,argv[i]);
		CORRECT_VELOCITIES = TRUE;
		i++;
		}
	else if(strcmp(argv[i-1],"-B")==0)
		{
		sscanf(argv[i],"%d",&pic_y);
		sscanf(argv[i+1],"%d",&pic_x);
		BINARY = TRUE;
		i += 2;
		}
	else if(strcmp(argv[i-1],"-s")==0)
		{
		sscanf(argv[i],"%d",&STEP);
		i += 1;
		}
	else if(strcmp(argv[i-1],"-P1")==0)
		{
		PREVIOUS1 = TRUE;
		}
	else if(strcmp(argv[i-1],"-NL")==0)
		{
		LAPLACIAN = FALSE;
		}
	else if(strcmp(argv[i-1],"-PR1")==0)
		{
		PRINT_FLOWS1 = TRUE;
		}
	else if(strcmp(argv[i-1],"-PR2")==0)
		{
		PRINT_FLOWS2 = TRUE;
		}
	else if(strcmp(argv[i-1],"-PRI")==0)
		{
		PRINT_ITERATIONS = TRUE;
		}
	else if(strcmp(argv[i-1],"-PRL")==0)
		{
		PRINT_LAPLACIAN = TRUE;
		}
	else if(strcmp(argv[i-1],"-w")==0)
		{
		sscanf(argv[i],"%d",&w);
		if(w!=1 && w!=2)
			{
			printf("Fatal error: w must be 1 or 2\n");
			exit(1);
			}
		i+= 1;
		}
	else if(strcmp(argv[i-1],"-n")==0)
		{
		sscanf(argv[i],"%d",&n);
		i+= 1;
		}
	else if(strcmp(argv[i-1],"-N")==0)
		{
		sscanf(argv[i],"%d",&N);
		i+= 1;
		}
	else if(strcmp(argv[i-1],"-i")==0)
		{
		sscanf(argv[i],"%d",&MAX_NUMBER_OF_ITERATIONS);
		i+= 1;
		}
	else if(strcmp(argv[i-1],"-P2")==0)
		{
		PREVIOUS2 = TRUE;
		}
	else if(strcmp(argv[i-1],"-SUB")==0)
		{
		sscanf(argv[i],"%d",&EXTRA_OFFSET_X);
		sscanf(argv[i+1],"%d",&EXTRA_OFFSET_Y);
		sscanf(argv[i+2],"%d",&CUT_X);
		sscanf(argv[i+3],"%d",&CUT_Y);
		SUBSET = TRUE;
		i+=4;
		}
	else    {
		printf("Fatal error invalid %dth argument: %s\n",i,argv[i-1]);
		exit(1);
		}
	i++;
	}
if(!PREVIOUS1 && PREVIOUS2) STEP2 = FALSE;
if(threshold1)
	printf("Step 1 velocities are thresholded for tau=%f to %f in increments of %f\n",
		low_tau1,high_tau1,inc_tau1);
if(threshold2)
	printf("Step 2 velocities are thresholded for tau=%f to %f in increments of %f\n",
		low_tau2,high_tau2,inc_tau2);
printf("Velocity Range: -%d <= u,v <= %d\n",N,N);
printf("n=%d neighbourhood size for SSD surface calculation\n",n);
printf("w=%d window size for step 2\n",w);
printf("Left image: %d\n",central_image-STEP);
printf("Central image: %d\n",central_image);
printf("Right image: %d\n",central_image+STEP);

/* Create appropriate file names for step 1 results */
if(!PREVIOUS1)
{
if(STEP==1)
sprintf(outname,"%s/singh.step1.%sF-n-%d-w-%d-N-%d",path2,argv[1],n,w,N);
else
sprintf(outname,"%s/singh.step1.%sF-n-%d-w-%d-N-%d-s-%d",path2,argv[1],n,w,N,STEP);
if((fdf1 = creat(outname,0755))<0)
	{
	printf("Error creating file %s.\n",outname);
	exit(1);
	}
printf("File %s opened for output\n",outname);
if(STEP==1)
sprintf(outname,"%s/singh.step1.%sC-n-%d-w-%d-N-%d",path2,argv[1],n,w,N);
else
sprintf(outname,"%s/singh.step1.%sC-n-%d-w-%d-N-%d-s-%d",path2,argv[1],n,w,N,STEP);
if((fdcov1 = creat(outname,0755))<0)
	{
	printf("Error creating file %s.\n",outname);
	exit(1);
	}
printf("File %s opened for output\n",outname);
}
else
{
if(STEP==1)
sprintf(outname,"%s/singh.step1.%sF-n-%d-w-%d-N-%d",path2,argv[1],n,w,N);
else
sprintf(outname,"%s/singh.step1.%sF-n-%d-w-%d-N-%d-s-%d",path2,argv[1],n,w,N,STEP);
if((fdf1 = open(outname,O_RDONLY))<0)
	{
	printf("Fatal error: file %s does not exit\n",outname);
	exit(1);
	}
printf("File %s opened for input\n",outname);
if(STEP==1)
sprintf(outname,"%s/singh.step1.%sC-n-%d-w-%d-N-%d",path2,argv[1],n,w,N);
else
sprintf(outname,"%s/singh.step1.%sC-n-%d-w-%d-N-%d-s-%d",path2,argv[1],n,w,N,STEP);
if((fdcov1 = open(outname,O_RDONLY))<0)
	{
	printf("Fatal error: file %s does not exist\n",outname);
	exit(1);
	}
printf("File %s opened for input\n",outname);
}

/* Create appropriate file names for step 2 results */
if(!PREVIOUS2)
{
if(STEP==1)
sprintf(outname,"%s/singh.step2.%sF-n-%d-w-%d-N-%d",path2,argv[1],n,w,N);
else
sprintf(outname,"%s/singh.step2.%sF-n-%d-w-%d-N-%d-s-%d",path2,argv[1],n,w,N,STEP);
if((fdf2 = creat(outname,0755))<0)
	{
	printf("Error creating file %s.\n",outname);
	exit(1);
	}
printf("File %s opened for output\n",outname);
if(STEP==1)
sprintf(outname,"%s/singh.step2.%sC-n-%d-w-%d-N-%d",path2,argv[1],n,w,N);
else
sprintf(outname,"%s/singh.step2.%sC-n-%d-w-%d-N-%d-s-%d",path2,argv[1],n,w,N,STEP);
if((fdcov2 = creat(outname,0755))<0)
	{
	printf("Error creating file %s.\n",outname);
	exit(1);
	}
printf("File %s opened for output\n",outname);
fflush(stdout);
}
else
{
if(STEP==1)
sprintf(outname,"%s/singh.step2.%sF-n-%d-w-%d-N-%d",path2,argv[1],n,w,N);
else
sprintf(outname,"%s/singh.step2.%sF-n-%d-w-%d-N-%d-s-%d",path2,argv[1],n,w,N,STEP);
if((fdf2 = open(outname,O_RDONLY))<0)
	{
	printf("Fatal error: file %s does not exit\n",outname);
	exit(1);
	}
printf("File %s opened for input\n",outname);
if(STEP==1)
sprintf(outname,"%s/singh.step2.%sC-n-%d-w-%d-N-%d",path2,argv[1],n,w,N);
else
sprintf(outname,"%s/singh.step2.%sC-n-%d-w-%d-N-%d-s-%d",path2,argv[1],n,w,N,STEP);
if((fdcov2 = open(outname,O_RDONLY))<0)
	{
	printf("Fatal error: file %s does not exist\n",outname);
	exit(1);
	}
printf("File %s opened for input\n",outname);
}

/* Process correct velocity data */
if(strcmp("unknown",correct_filename)!=0)
{
/* Read the correct velocity data */
fd_correct = open(correct_filename,O_RDONLY);
no_bytes = 0;
no_bytes += read(fd_correct,&actual_x,4);
no_bytes += read(fd_correct,&actual_y,4);
no_bytes += read(fd_correct,&size_x,4);
no_bytes += read(fd_correct,&size_y,4);
no_bytes += read(fd_correct,&offset_x,4);
no_bytes += read(fd_correct,&offset_y,4);
if(offset_x != 0.0 || offset_y != 0.0 || actual_x != size_x || actual_y != size_y)
	{
	printf("Fatal error: something wrong with correct velocity data\n");
	printf("Actual y: %f Actual x: %f\n",actual_y,actual_x);
	printf("Size y: %f Size x: %f\n",size_y,size_x);
	printf("Offset y: %f Offset x: %f\n",offset_y,offset_x);
	exit(1);
	}
int_size_x = size_x;
int_size_y = size_y;
for(i=0;i<int_size_y;i++)
	no_bytes += read(fd_correct,&correct_velocities[i][0][0],int_size_x*8);
printf("\nFile %s opened and read\n",correct_filename);
printf("Size of correct velocity data: %d %d --- ",int_size_x,int_size_y);
printf("%d bytes read\n",no_bytes);
fflush(stdout);
}

/* Perform timing measurements */
temporal_offset = 100; /* Each program requires at 100 seconds */
if(STEP!=1)
sprintf(filename,"%s.%stime-n-%d-w-%d-s-%d",argv[0],argv[1],n,w,STEP);
else
sprintf(filename,"%s.%stime-n-%d-w-%d",argv[0],argv[1],n,w);
fp = fopen(filename,"w");
fprintf(fp,"The time data is in file: %s\n",filename);
time1 = time(NULL);
fprintf(fp,"Start time: %d\n",time1);
fflush(fp);
time2 = time1;

/* Read data */
sprintf(name,"%s/%s",path1,argv[1]);
read_image_files(name,inpic,&pic_x,&pic_y,header_ints,central_image,STEP);

/*****************************************************************/
/* Set parameters for a subset of the image sequence             */
/*****************************************************************/
old_pic_x = pic_x; old_pic_y = pic_y;
if(SUBSET)
{
printf("\n");
pic_x = EXTRA_OFFSET_X+CUT_X+OFFSET_X; 
pic_y = EXTRA_OFFSET_Y+CUT_Y+OFFSET_Y;
EXTRA_OFFSET_X -= OFFSET_X;
EXTRA_OFFSET_Y -= OFFSET_Y;
if(CUT_X <= 5 || CUT_Y <=5)
	{
	printf("Fatal error: subarea too small - must be at least 6*6\n");
	printf("Subarea specified: %d*%d\n",CUT_X,CUT_Y);
	exit(1);
	}
if((CUT_X+OFFSET_X+EXTRA_OFFSET_X > old_pic_x-OFFSET_X) ||
   (CUT_Y+OFFSET_Y+EXTRA_OFFSET_Y > old_pic_y-OFFSET_Y))
	{
	printf("Fatal error: subarea parameters incorrect\n");
	printf("Starting Coordinates: %d,%d  ",OFFSET_X+EXTRA_OFFSET_X,OFFSET_Y+EXTRA_OFFSET_Y);
	printf("Size: %d * %d\n",CUT_X,CUT_Y);
	printf("Too big by %d,%d\n",CUT_X+EXTRA_OFFSET_X+2*OFFSET_X-old_pic_x,
			    	    CUT_Y+EXTRA_OFFSET_Y+2*OFFSET_Y-old_pic_y);
	exit(1);
	}
printf("\n\nSubset of image sequence used\n");
printf("Starting Coordinates: %d,%d  ",OFFSET_X+EXTRA_OFFSET_X,OFFSET_Y+EXTRA_OFFSET_Y);
printf("Size: %d * %d\n",CUT_X,CUT_Y);
}

/* Compute laplacian images */
if(LAPLACIAN && !PREVIOUS1) 
	{
	laplacian(inpic,fpic,1.0,pic_x,pic_y); 
	if(PRINT_LAPLACIAN) write_image_files(argv[1],inpic,pic_x,pic_y,header_ints,central_image,STEP);
	}
else make_float(inpic,fpic,pic_x,pic_y);

/**********************************************/
/* Step 1: Conservation Information recovery  */
/**********************************************/
if(!PREVIOUS1)
{
comp1(fpic,pic_x,pic_y,Ucc,Scc,n,N,STEP,old_pic_x,old_pic_y);

/* Output velocities and rearrange them */
printf("\nOutputing flow data\n");
output_velocities(fdf1,Ucc,old_pic_x,old_pic_y,"Full",SAMPLE);
output_covariances(fdcov1,Scc,old_pic_x,old_pic_y);

}
else
{
input_velocities(fdf1,Ucc,old_pic_x,old_pic_y);
input_covariances(fdcov1,Scc,old_pic_x,old_pic_y);
}

if(strcmp("unknown",correct_filename)!=0)
{
printf("\nStep 1:\n");
calc_statistics(correct_velocities,Ucc,old_pic_x,old_pic_y,
	&ave_error,&st_dev,&density,&min_angle,&max_angle,SAMPLE);
printf("\nError Analysis Performed\n");
printf("Average Error: %f Standard Deviation: %f Density: %f\n",
	ave_error,st_dev,density);
printf("Minimum Angle Error: %f Maximum Angle Error: %f\n",min_angle,max_angle);
fflush(stdout);

/* Threshold step 1 velocities */
if(threshold1)
	{
	tau1 = low_tau1;
	printf("\n-----------------------------------------------------\n");
	printf("Step 1 Thresholding");
	printf("\n-----------------------------------------------------\n");
	while(tau1 <= high_tau1)
	{
	threshold_velocities(Ucc,Scc,threshold_velocity,pic_x,pic_y,tau1);
	printf("\nThreshold: %f\n",tau1);
	calc_statistics(correct_velocities,threshold_velocity,old_pic_x,old_pic_y,
		&ave_error,&st_dev,&density,&min_angle,&max_angle,SAMPLE);
	printf("\nError Analysis Performed\n");
	printf("Average Error: %f Standard Deviation: %f Density: %f\n",
		ave_error,st_dev,density);
	printf("Minimum Angle Error: %f Maximum Angle Error: %f\n",min_angle,max_angle);
	fflush(stdout);
	if(PRINT_FLOWS1)
	{
	sprintf(outname,"%s/singh.step1.%sF-tau-%4.2f",path2,argv[1],tau1);
	if((fd_vels = creat(outname,0755))<0)
	{
	printf("Error creating file %s.\n",outname);
	exit(1);
	}
	printf("File %s opened for output\n",outname);
	output_velocities(fd_vels,threshold_velocity,old_pic_x,old_pic_y,"Full",SAMPLE);
	}
	tau1 += inc_tau1;
	}
	printf("\n-----------------------------------------------------\n");
	}
}

if(!STEP2) exit(1);

/**********************************************/
/* Step 2: Neighbourhood Information recovery */
/**********************************************/
if(!PREVIOUS2)
{
comp2(path2,argv[1],Ucc,Scc,S,pic_x,pic_y,full_velocities,correct_filename,w,old_pic_x,old_pic_y); 


/* Output velocities for step 2 computation  */
printf("\nOutputing flow data for step 2\n");
output_velocities(fdf2,full_velocities,old_pic_x,old_pic_y,"Full",1);
output_covariances(fdcov2,S,old_pic_x,old_pic_y);
}
else
{
offset = (2*w+1)/2;
OFFSET_X += offset;
OFFSET_Y += offset; 
input_velocities(fdf2,full_velocities,old_pic_x,old_pic_y);
input_covariances(fdcov2,S,old_pic_x,old_pic_y);
}

if(strcmp("unknown",correct_filename)!=0)
{
printf("\nStep 2:\n");
calc_statistics(correct_velocities,full_velocities,old_pic_x,old_pic_y,
	&ave_error,&st_dev,&density,&min_angle,&max_angle,1);
printf("\nError Analysis Performed\n");
printf("Average Error: %f Standard Deviation: %f Density: %f\n",
	ave_error,st_dev,density);
printf("Minimum Angle Error: %f Maximum Angle Error: %f\n",min_angle,max_angle);

/* Threshold step 2 velocities */
if(threshold2)
	{
	tau2 = low_tau2;
	printf("\n-----------------------------------------------------\n");
	printf("Step 2 Thresholding");
	printf("\n-----------------------------------------------------\n");
	while(tau2 <= high_tau2)
	{
	threshold_velocities(full_velocities,S,threshold_velocity,pic_x,pic_y,tau2);
	printf("\nThreshold: %f\n",tau2);
	calc_statistics(correct_velocities,threshold_velocity,old_pic_x,old_pic_y,
		&ave_error,&st_dev,&density,&min_angle,&max_angle,SAMPLE);
	printf("\nError Analysis Performed\n");
	printf("Average Error: %f Standard Deviation: %f Density: %f\n",
		ave_error,st_dev,density);
	printf("Minimum Angle Error: %f Maximum Angle Error: %f\n",min_angle,max_angle);
	if(PRINT_FLOWS2)
	{
	sprintf(outname,"%s/singh.step2.%sF-tau-%4.2f",path2,argv[1],tau2);
	if((fd_vels = creat(outname,0755))<0)
	{
	printf("Error creating file %s.\n",outname);
	exit(1);
	}
	printf("File %s opened for output\n",outname);
	output_velocities(fd_vels,threshold_velocity,pic_x,pic_y,"Full",SAMPLE);
	}
	tau2 += inc_tau2;
	}
	printf("\n-----------------------------------------------------\n");
	}
}

printf("\nProcessing Finished\n");
fflush(stdout);

time2 = time(NULL);
if(time2 > time1+temporal_offset)
{
fprintf(fp,"\nEnd time: %d\n",time2);
fprintf(fp,"End Time in seconds: %d\n",(time2-time1));
fprintf(fp,"End Time in minutes: %f\n",(time2-time1)*1.0/60.0);
fprintf(fp,"End Time in hours: %f\n",(time2-time1)*1.0/3600.0);
fflush(fp);
}
}



/************************************************************
   Read raster files into internal 3-D array.
   Also sets pic_x and pic_y (the actual picture dimensions)
   which are global
************************************************************/
read_image_files(s,pic,pic_x,pic_y,header_ints,central_image,STEP)
char *s;
cimages pic;
int header_ints[8],*pic_x,*pic_y,central_image,STEP;
{
char fname[100];
int i,j,fd,ints[8],ONCE,bytes;
unsigned char header[HEAD];

ONCE = TRUE;
for(i=0;i<NUMFILES;i++) 
	{
	sprintf(fname,"%s%d",s,central_image+(i-1)*STEP);
 	if((fd = open(fname,O_RDONLY)) >0)
		{
		if(!BINARY)
		{
		if(ONCE)
			{
			read(fd,ints,HEAD);
			for(j=0;j<8;j++) header_ints[j] = ints[j];
			(*pic_x) = ints[1];
			(*pic_y) = ints[2];
			ONCE = FALSE;
			if((*pic_y) > PIC_Y || (*pic_x) > PIC_X)
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
		for(j=0;j<(*pic_y);j++)
			bytes += read(fd,&pic[i][j][0],(*pic_x));
		printf("File %s read -- %d bytes\n",fname,bytes);
		bytes = 0;
		fflush(stdout);
		}
	      else 
		{
		printf("File %s does not exist in read_image_files.\n",fname);
		exit(1);
		}
	}
} /* End of read_image_files */



/************************************************************
   write raster files from internal 3-D array.
   Also sets pic_x and pic_y (the actual picture dimensions)
   which are global
************************************************************/
write_image_files(s,pic,pic_x,pic_y,header_ints,central_image,STEP)
char *s;
cimages pic;
int header_ints[8],pic_x,pic_y,central_image,STEP;
{
char fname[100];
int i,j,fd,ints[8],ONCE,bytes;
unsigned char header[HEAD];

ONCE = TRUE;
for(i=0;i<NUMFILES;i++) 
	{
	sprintf(fname,"laplacian.%s%d",s,central_image+(i-1)*STEP);
 	if((fd = creat(fname,0755)) >= 0)
		{
		write(fd,header_ints,HEAD);
		bytes = 32;

		/* Read row by row */
		for(j=0;j<(pic_y);j++)
			bytes += write(fd,&pic[i][j][0],(pic_x));
		printf("File %s written -- %d bytes\n",fname,bytes);
		fflush(stdout);
		}
	      else 
		{
		printf("File %s does not exist in write_image_files.\n",fname);
		exit(1);
		}
	}
} /* End of write_image_files */




/*****************************************************************/
/* Compute velocities						 */
/*****************************************************************/
comp1(fpic,pic_x,pic_y,Ucc,Scc,n,N,STEP,old_pic_x,old_pic_y)
float Ucc[PIC_Y][PIC_X][2],Scc[PIC_Y][PIC_X][2][2];
fimages fpic;
int pic_x,pic_y,n,N,STEP,old_pic_x,old_pic_y;
{
int i,j,vel_type,SKIP,ct,max_a,max_b;
float vx,vy,step_size;
float uve[2],uva[2],N2;
float SSDdata[ARRSIZE][ARRSIZE],SI[2][2],condition_number;
float covariance[2][2],mean[2],eigenvalues[2],eigenvector1[2],eigenvector2[2];
	
printf("\nComputing Velocity Information via step 1...\n");

for(i=0;i<PIC_Y;i++)
for(j=0;j<PIC_X;j++)
	{
	Ucc[i][j][0] = Ucc[i][j][1] = NO_VALUE;
	Scc[i][j][0][0] = Scc[i][j][0][1] = NO_VALUE;
	Scc[i][j][1][0] = Scc[i][j][1][1] = NO_VALUE;
	}

N2 = 2.0*N;
step_size = 1.0*STEP;

/* Compute velocity for every location */
printf("Processing rows for step 1:\n");
ct = 0;
for(i=OFFSET_Y+EXTRA_OFFSET_Y;i<pic_y-OFFSET_Y;i++) 
{
printf("%3d ",i);
ct++;
if((ct % 15) == 0) printf("\n");
fflush(stdout);
for(j=OFFSET_X+EXTRA_OFFSET_X;j<pic_x-OFFSET_X;j++) 
{
if(i%SAMPLE==0 && j%SAMPLE==0)
	{
	fflush(stdout);
	vel_type = NOTHING;
	compute_SSD_surface(fpic,i,j,n,N,SSDdata,&max_a,&max_b);
	calc_mean_and_covariance1(SSDdata,mean,covariance,
		-4.0+max_b,4.0+max_b,-4.0+max_a,4.0+max_a,
		eigenvalues,eigenvector1,eigenvector2,&vel_type,N);
	vx =  mean[0];
	vy = -mean[1];
	Scc[i][j][0][0] = covariance[0][0];
	Scc[i][j][0][1] = covariance[0][1];
	Scc[i][j][1][0] = covariance[1][0];
	Scc[i][j][1][1] = covariance[1][1];
	Ucc[i][j][0] = vx/step_size;
	Ucc[i][j][1] = vy/step_size;
	}
}
}
}
 
/*************************************************************/
/* Smooth velocity field using velocities computed in step 1 */
/* Step 2 Neighbourhood Information			     */
/*************************************************************/
comp2(path,name,Ucc,Scc,S,pic_x,pic_y,full_velocity,correct_filename,w,old_pic_x,old_pic_y)
float Ucc[PIC_Y][PIC_X][2],Scc[PIC_Y][PIC_X][2][2];
int pic_x,pic_y,w,old_pic_x,old_pic_y;
float full_velocity[PIC_Y][PIC_X][2],S[PIC_Y][PIC_X][2][2];
unsigned char correct_filename[100],name[100],path[100];
{
int fd,i,j,k,l,FIRST,SECOND,weighted_size,RELAXATION,temp,ct;
float neighbour_velocities[WEIGHTSIZE][WEIGHTSIZE][2];
float weights[WEIGHTSIZE][WEIGHTSIZE],diff[2],mean[2];
float SnI[2][2],vec1[2],vec2[2],bigdiff[PIC_Y][PIC_X][2];
float eigenvector1[2],eigenvector2[2],eigenvalues[2],max_diff;
float ave_error,st_dev,min_angle,max_angle,density,size;
int iteration_number,offset,SINGULAR,no_vels;

FIRST = 0;
SECOND = 1;

max_diff = 0.0;
for(i=0;i<PIC_Y;i++)
for(j=0;j<PIC_X;j++)
	{
	bigdiff[i][j][0] = bigdiff[i][j][1] = NO_VALUE;
	Un[FIRST][i][j][0] = Un[FIRST][i][j][1] = NO_VALUE;
	Un[SECOND][i][j][0] = Un[SECOND][i][j][1] = NO_VALUE;
	full_velocity[i][j][0] = full_velocity[i][j][1] = NO_VALUE;
	}
for(i=0;i<old_pic_y;i++)
for(j=0;j<old_pic_x;j++)
	{
	Un[FIRST][i][j][0] = Un[SECOND][i][j][0] = Ucc[i][j][0];
	Un[FIRST][i][j][1] = Un[SECOND][i][j][1] = Ucc[i][j][1];
	}
if(strcmp("unknown",correct_filename)!=0)
{
calc_statistics(correct_velocities,&Un[FIRST][0][0][0],old_pic_x,old_pic_y,
	&ave_error,&st_dev,&density,&min_angle,&max_angle,1);
printf("Error Analysis Performed at start of step 2\n");
printf("Average Error: %f Standard Deviation: %f Density: %f\n",
	ave_error,st_dev,density);
printf("Minimum Angle Error: %f Maximum Angle Error: %f\n\n",min_angle,max_angle);
fflush(stdout);
} 

compute_weights(weights,w);
weighted_size = 2*w+1;
offset = weighted_size/2;
OFFSET_X += offset;
OFFSET_Y += offset;

/* Initialization: compute SccI and multiply by Ucc, compute the
   weighted_size*weighted_size velocity neighbourhood and its
   mean and covariance matrix */

printf("Computing velocity information via step 2\n");
printf("Performing Initialization...\n");
for(i=OFFSET_Y+EXTRA_OFFSET_Y;i<pic_y-OFFSET_Y;i++) 
{
for(j=OFFSET_X+EXTRA_OFFSET_X;j<pic_x-OFFSET_X;j++) 
	{
	/* Compute the inverse of Scc and multiply by Ucc */
	SINGULAR=inverse22(&Scc[i][j][0][0],&SccI[i][j][0][0],SVD_PRINT);
	mult21(&SccI[i][j][0][0],&Ucc[i][j][0],&SccI_Ucc[i][j][0]);

	/* Compute neighbourhood velocities for initialization */
	no_vels = 0;
	for(k=(-offset);k<=offset;k++)
	for(l=(-offset);l<=offset;l++)
		{
		if(Un[FIRST][i+k][j+l][0] != NO_VALUE &&
		   Un[FIRST][i+k][j+l][1] != NO_VALUE)
		   {
		   no_vels++;
		   neighbour_velocities[k+offset][l+offset][0] = Un[FIRST][i+k][j+l][0];
		   neighbour_velocities[k+offset][l+offset][1] = Un[FIRST][i+k][j+l][1];
		   }
		else
		   {
		   printf("FIRST i=%d j=%d k=%d l=%d\n",i,j,k,l);
		   for(k=(-offset);k<=offset;k++)
		   for(l=(-offset);l<=offset;l++)
		   printf("vel for k=%d l=%d: %f,%f\n",k,l,Un[FIRST][i+k][j+l][0],Un[FIRST][i+k][j+l][1]);
		   printf("Fatal error: computed velocity undefined during initialization\n");
		   exit(1);
		   }
		}
	calc_mean_and_covariance2(weights,neighbour_velocities,
		weighted_size,&Ua[i][j][0],&Sn[FIRST][i][j][0][0]);
	}
}
printf("\nInitialization complete for step 2\n");

/* Relaxation Computation using neighbourhood velocities */
RELAXATION = TRUE;
iteration_number = 0;
printf("\n");
while(RELAXATION && iteration_number < MAX_NUMBER_OF_ITERATIONS)
{
printf("Iteration %d\n",iteration_number);
iteration_number++;
RELAXATION = FALSE;
size = 0.0;
for(i=OFFSET_Y+EXTRA_OFFSET_Y;i<pic_y-OFFSET_Y;i++) 
{
for(j=OFFSET_X+EXTRA_OFFSET_X;j<pic_x-OFFSET_X;j++) 
	/* Perform an iteration */
	{
	SINGULAR = inverse22(&Sn[FIRST][i][j][0][0],SnI,SVD_PRINT);
	add22(&SccI[i][j][0][0],SnI,&Ssum[i][j][0][0]);
	SINGULAR += inverse22(&Ssum[i][j][0][0],&SsumI[i][j][0][0],SVD_PRINT);
	mult21(SnI,&Ua[i][j][0],vec1);
	add21(&SccI_Ucc[i][j][0],vec1,vec2);
	mult21(&SsumI[i][j][0][0],vec2,&Un[SECOND][i][j][0]);
	/*
	if(i==32 && j==74)
	{
	printf("SsumI[32][74]:\n"); 
	printf(" %f %f\n",SsumI[32][74][0][0],SsumI[32][74][0][1]);
	printf(" %f %f\n",SsumI[32][74][1][0],SsumI[32][74][1][1]);
	printf("SccI[32][74]:\n"); 
	printf(" %f %f\n",SccI[32][74][0][0],SccI[32][74][0][1]);
	printf(" %f %f\n",SccI[32][74][1][0],SccI[32][74][1][1]);
	printf("SnI:\n"); 
	printf(" %f %f\n",SnI[0][0],SnI[0][1]);
	printf(" %f %f\n",SnI[1][0],SnI[1][1]);
	printf("vec1: %f %f\n",vec1[0],vec1[1]);
	printf("vec2: %f %f\n",vec2[0],vec2[1]);
	printf("SccI_Ucc[32][74]: %f %f\n",SccI_Ucc[32][74][0],SccI_Ucc[32][74][1]);
	printf("Ua[32][74]: %f %f\n",Ua[32][74][0],Ua[32][74][1]);
	printf("U[FIRST][32][74]: %f %f\n",Un[FIRST][32][74][0],Un[FIRST][32][74][1]);
	printf("U[SECOND][32][74]: %f %f\n",Un[SECOND][32][74][0],Un[SECOND][32][74][1]);
	}
	*/
	diff[0] = Un[FIRST][i][j][0]-Un[SECOND][i][j][0];
	diff[1] = Un[FIRST][i][j][1]-Un[SECOND][i][j][1];
	bigdiff[i][j][0] =  diff[0];
	bigdiff[i][j][1] =  diff[1];
	size  = L2norm(diff,2);
	if(size > max_diff) max_diff = size;
	if(size > DIFF_THRESH) RELAXATION = TRUE;
	}
}

if(!RELAXATION) printf("\nConvergence detected - iterative calculations are stopped\n");

/* Prepare for next iteration */
if(RELAXATION)
{
for(i=OFFSET_Y+EXTRA_OFFSET_Y;i<pic_y-OFFSET_Y;i++) 
{
for(j=OFFSET_X+EXTRA_OFFSET_X;j<pic_x-OFFSET_X;j++) 
{
/* Compute neighbourhood velocities for initialization */
for(k=(-offset);k<=offset;k++)
for(l=(-offset);l<=offset;l++)
	{
	if(Un[SECOND][i+k][j+l][0] != NO_VALUE && 
	   Un[SECOND][i+k][j+l][1] != NO_VALUE)
	   {
	   neighbour_velocities[k+offset][l+offset][0] = Un[SECOND][i+k][j+l][0];
	   neighbour_velocities[k+offset][l+offset][1] = Un[SECOND][i+k][j+l][1];
	   }
	else
	   {
	   printf("SECOND i=%d j=%d k=%d l=%d\n",i,j,k,l);
	   for(k=(-offset);k<=offset;k++)
	   for(l=(-offset);l<=offset;l++)
		printf("vel for k=%d l=%d: %f,%f\n",k,l,Un[SECOND][i+k][j+l][0],Un[SECOND][i+k][j+l][1]);
	   printf("\nFatal error: computed velocity undefined during iteration %d\n",iteration_number);
	   exit(1); 
	   }
	calc_mean_and_covariance2(weights,neighbour_velocities,weighted_size,
			  &Ua[i][j][0],&Sn[SECOND][i][j][0][0]);
	}
}
}

if(PRINT_ITERATIONS)
{
sprintf(outname,"%s/singh.iteration.step2.%sF-%d",path,name,iteration_number-1);
if((fd = creat(outname,0755))<0)
	{
	printf("Error creating file %s.\n",outname);
	exit(1);
	}
printf("\nFile %s opened\n",outname);
output_velocities(fd,&Un[SECOND][0][0][0],old_pic_x,old_pic_y,"Full",1);
}

/* Compute error statistics */
if(strcmp("unknown",correct_filename)!=0)
{
calc_statistics(correct_velocities,&Un[SECOND][0][0][0],old_pic_x,old_pic_y,
	&ave_error,&st_dev,&density,&min_angle,&max_angle,1);
printf("Error Analysis Performed\n");
printf("Average Error: %f Standard Deviation: %f Density: %f\n",
	ave_error,st_dev,density);
printf("Minimum Angle Error: %f Maximum Angle Error: %f\n",min_angle,max_angle);
printf("L2norm of difference: %14.9f\n",bigL2norm(bigdiff,pic_x));
printf("Maximum individual velocity difference: %f\n\n",max_diff);
fflush(stdout);
} 
} 
/* Switch the values of FIRST and SECOND before next iteration */
temp = FIRST;
FIRST = SECOND;
SECOND = temp;
}


no_vels = 0;
for(i=OFFSET_Y+EXTRA_OFFSET_Y;i<pic_y-OFFSET_Y;i++)
for(j=OFFSET_X+EXTRA_OFFSET_X;j<pic_x-OFFSET_X;j++)
	{
	full_velocity[i][j][0] = Un[FIRST][i][j][0];
	full_velocity[i][j][1] = Un[FIRST][i][j][1];
	S[i][j][0][0] = SsumI[i][j][0][0];
	S[i][j][1][0] = SsumI[i][j][1][0];
	S[i][j][0][1] = SsumI[i][j][0][1];
	S[i][j][1][1] = SsumI[i][j][1][1];
	if(!(full_velocity[i][j][0]==NO_VALUE && full_velocity[i][j][1]==NO_VALUE)) no_vels++;
	else { printf("Velocity at %d,%d undefined\n",i,j); }
	}


}

/************************************************************************/
/* Threshold the velocities based of eigenvalues of covariance matrices */
/************************************************************************/
threshold_velocities(input_velocity,covariances,full_velocity,pic_x,pic_y,tau)
float input_velocity[PIC_Y][PIC_X][2],covariances[PIC_Y][PIC_X][2][2];
float full_velocity[PIC_Y][PIC_X][2],tau;
int pic_x,pic_y;
{
int i,j;
float eigenvalues[2],eigenvector1[2],eigenvector2[2];

for(i=0;i<PIC_Y;i++)
for(j=0;j<PIC_X;j++)
	full_velocity[i][j][0] = full_velocity[i][j][1] = NO_VALUE;

/* Threshold the final result */
for(i=OFFSET_Y+EXTRA_OFFSET_Y;i<pic_y-OFFSET_Y;i++)
for(j=OFFSET_X+EXTRA_OFFSET_X;j<pic_x-OFFSET_X;j++)
	{
	if(!(input_velocity[i][j][0] == NO_VALUE && 
	     input_velocity[i][j][1] == NO_VALUE))
	{
	comp_eigen(&covariances[i][j][0][0],eigenvalues,eigenvector1,eigenvector2);
	if(fabs(eigenvalues[0]) < tau) 
		{
		full_velocity[i][j][0] = input_velocity[i][j][0];
		full_velocity[i][j][1] = input_velocity[i][j][1];
		}
	else
		{
		full_velocity[i][j][0] = NO_VALUE;
		full_velocity[i][j][1] = NO_VALUE;
		}
	}
	else
	{
	full_velocity[i][j][0] = NO_VALUE;
	full_velocity[i][j][1] = NO_VALUE;
	}
	}
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
   Output full velocities using Burkitt format
************************************************************/
output_velocities(fdf,velocities,pic_x,pic_y,type,SAMPLE)
float velocities[PIC_Y][PIC_X][2];
char type[100];
int fdf,pic_x,pic_y,SAMPLE;
{
float x,y;
int i,j,bytes,no_novals,no_vals;

if(fdf<0)
	{
	printf("\nFatal error: full velocity file not opened\n");
	exit(1);
	}
/* Original size */
x = pic_x;
y = pic_y;
write(fdf,&x,4);
write(fdf,&y,4);

/* Size of result data */
x = (pic_x-2*OFFSET_X);
y = (pic_y-2*OFFSET_Y);
write(fdf,&x,4);
write(fdf,&y,4);

/* Offset to start of data */
x = OFFSET_Y; 
y = OFFSET_X; 
write(fdf,&x,4);
write(fdf,&y,4);
bytes = 24;

no_novals = no_vals = 0;
/* Count velocity positions */
for(i=OFFSET_Y;i<pic_y-OFFSET_Y;i++)
for(j=OFFSET_X;j<pic_x-OFFSET_X;j++)
if(i%SAMPLE==0 && j%SAMPLE==0)
	{
	if(velocities[i][j][0] != NO_VALUE && 
	   velocities[i][j][1] != NO_VALUE)
		{
		no_vals++;
		}
	else
		{
		no_novals++;
		}
	}

for(i=OFFSET_Y;i<pic_y-OFFSET_Y;i++)
	{
	bytes += write(fdf,&velocities[i][OFFSET_X][0],(pic_x-(2*OFFSET_X))*8);
	}
close(fdf);
printf("\n%s velocities output: %d bytes\n",type,bytes);
printf("Number of positions with velocity: %d\n",no_vals);
printf("Number of positions without velocity: %d\n",no_novals);
printf("Percentage of %s velocities: %f\n",type,
	no_vals/(1.0*(no_vals+no_novals))*100.0);
fflush(stdout);
}


/************************************************************
   Output full velocities using Burkitt format
************************************************************/
output_covariances(fdf,covariances,pic_x,pic_y)
float covariances[PIC_Y][PIC_X][2][2];
int fdf,pic_x,pic_y;
{
float x,y;
int i,j,bytes,no_novals,no_vals;

if(fdf<0)
	{
	printf("\nFatal error: covariance file not opened\n");
	exit(1);
	}
/* Original size */
x = pic_x;
y = pic_y;
write(fdf,&x,4);
write(fdf,&y,4);

/* Size of result data */
x = (pic_x-2*OFFSET_X);
y = (pic_y-2*OFFSET_Y);
write(fdf,&x,4);
write(fdf,&y,4);

/* Offset to start of data */
x = OFFSET_Y; 
y = OFFSET_X; 
write(fdf,&x,4);
write(fdf,&y,4);
bytes = 24;


for(i=OFFSET_Y;i<pic_y-OFFSET_Y;i++)
	{
	bytes += write(fdf,&covariances[i][OFFSET_X][0][0],(pic_x-(2*OFFSET_X))*16);
	}
close(fdf);
printf("\nFull velocity covariance matrices output: %d bytes\n",bytes);
fflush(stdout);
}



/************************************************************
   Input full velocities using Burkitt format
************************************************************/
input_velocities(fdf,velocities,pic_x,pic_y)
float velocities[PIC_Y][PIC_X][2];
int fdf,pic_x,pic_y;
{
float x,y;
int i,j,bytes;

for(i=0;i<PIC_Y;i++)
for(j=0;j<PIC_X;j++)
	velocities[i][j][0] = velocities[i][j][1] = NO_VALUE; 

if(fdf<0)
	{
	printf("\nFatal error: full velocity file not opened\n");
	exit(1);
	}
/* Original size */
read(fdf,&x,4);
read(fdf,&y,4);

/* Size of result data */
read(fdf,&x,4);
read(fdf,&y,4);

/* Offset to start of data */
read(fdf,&x,4);
read(fdf,&y,4);
OFFSET_Y = x;
OFFSET_X = y;
bytes = 24;

for(i=OFFSET_Y;i<pic_y-OFFSET_Y;i++)
	{
	bytes += read(fdf,&velocities[i][OFFSET_X][0],(pic_x-(2*OFFSET_X))*8);
	}
close(fdf);
printf("\nVelocities input: %d bytes\n",bytes);
fflush(stdout);
}


/************************************************************
   Output full velocities using Burkitt format
************************************************************/
input_covariances(fdf,covariances,pic_x,pic_y)
float covariances[PIC_Y][PIC_X][2][2];
int fdf,pic_x,pic_y;
{
float x,y;
int i,j,bytes;

if(fdf<0)
	{
	printf("\nFatal error: covariance file not opened\n");
	exit(1);
	}
/* Original size */
read(fdf,&x,4);
read(fdf,&y,4);

/* Size of result data */
read(fdf,&x,4);
read(fdf,&y,4);

/* Offset to start of data */
read(fdf,&x,4);
read(fdf,&y,4);
bytes = 24;


for(i=OFFSET_Y;i<pic_y-OFFSET_Y;i++)
	{
	bytes += read(fdf,&covariances[i][OFFSET_X][0][0],(pic_x-(2*OFFSET_X))*16);
	}
close(fdf);
printf("\nFull velocity covariance matrices input: %d bytes\n",bytes);
fflush(stdout);
}





/***************************************************************************/
/* Compute average angle, standard deviation and density (as a percentage) */
/***************************************************************************/
calc_statistics(correct_vels,full_vels,pic_x,pic_y,
		ave_error,st_dev,density,min_angle,max_angle,SAMPLE)
float full_vels[PIC_Y][PIC_X][2],*ave_error,*density,*st_dev;
float correct_vels[PIC_Y][PIC_X][2],*min_angle,*max_angle;
int pic_x,pic_y,SAMPLE;
{
int count,no_count,i,j;
float sumX2,temp,uva[2],uve[2],sum2,temp1,sum1;
float error[PIC_Y][PIC_X];

count = no_count = 0;
sum2 = sumX2 = 0.0;
(*min_angle) = HUGE;
(*max_angle) = 0.0;
(*ave_error) = (*st_dev) = 0.0;
(*density) = 0.0;
for(i=OFFSET_Y;i<pic_y-OFFSET_Y;i++)
{
for(j=OFFSET_X;j<pic_x-OFFSET_X;j++)
if(i%SAMPLE==0 && j%SAMPLE==0)
	{
	if(full_vels[i][j][0] == NO_VALUE && full_vels[i][j][1] == NO_VALUE)
		{
		no_count++;
		error[i][j] = NO_VALUE;
		}
	else 
	  {
	  count++;
	  uve[0] = full_vels[i][j][0]; uve[1] = full_vels[i][j][1];
	  uva[0] = correct_vels[i][j][0]; uva[1] = correct_vels[i][j][1];
	  temp = PsiER(uve,uva,i,j);
	  error[i][j] = temp;
	  (*ave_error) += temp;
	  sumX2 += temp*temp;
	  if(temp < (*min_angle)) (*min_angle) = temp;
	  if(temp > (*max_angle)) (*max_angle) = temp;
	  }
	}
}
if(count != 0) (*ave_error) = (*ave_error)/(count*1.0);
if(count > 1) 
	{
	sum1 = 0.0;
	for(i=OFFSET_Y;i<pic_y-OFFSET_Y;i++)
	{
	for(j=OFFSET_X;j<pic_x-OFFSET_X;j++)
	if(i%SAMPLE==0 && j%SAMPLE==0)
	if(error[i][j] != NO_VALUE)
			{
			temp1 = error[i][j] - (*ave_error);
			temp1 = temp1*temp1;
			sum1 += temp1;
			}
	}
	(*st_dev) = sqrt(sum1/(count-1));
	}
else (*st_dev) = 0.0;

(*density) = (count*100.0)/(count+no_count);

if((*ave_error) == 0.0) { (*min_angle) = (*max_angle) = 0.0; }
fflush(stdout);
}

/************************************************************
 Full Image Velocity Angle Error
************************************************************/
float PsiER(ve,va,i,j)
float ve[2],va[2];
int i,j;
{
float nva;
float nve;
float v,r;
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

r = acos(v)*TODEG;

if(!(r>=0.0 && r< 180.0))
{
printf("ERROR in PsiER()...\n r=%8.4f  v=%8.4f  nva=%8.4f nve= %8.4f\n",
	r,v,nva,nve);
printf("va=(%f,%f) ve=(%f,%f)\n",va[0],va[1],ve[0],ve[1]);
printf("i=%d j=%d\n",i,j);
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
float v1,v2;
float n[2];

nva = L2norm(va,2), nve = L2norm(ve,2);
if(nve > 0.00000001)
{
n[0] = ve[0]/nve;
n[1] = ve[1]/nve;
v1 = (va[0]*n[0] + va[1]*n[1]-nve) ;
v2 = v1/(sqrt((1+nva*nva))*sqrt((1+nve*nve)));
v1 =  asin(v2)*TODEG;

if(!(v1>=-90.0 && v1<=90.0))
	{
       	printf("ERROR in PsiEN()  v1: %f ve: %f\n",v1,v2);
	printf("nve: %f nva: %f\n",nve,nva);
	printf("n: %f %f\n",n[0],n[1]);
	printf(" ve: %f %f va: %f %f\n",ve[0],ve[1],va[0],va[1]);
	fflush(stdout);
	}
}
else v1 = NO_VALUE;
	
return fabs(v1);
}




/**************************************************************************/
/* Compute all eigenvalues and eigenvectors of a real symmetric matrix    */
/* a[DIM][DIM]. On output elements of a above the disgonal are destroyed. */
/* d[DIM] returns the eigenvalues of a. v[DIM][DIM] is a matrix whose     */
/* columns  contain, on output, the normalized eigenvectors of a. nrot    */
/* returns the number of Jacobi rotations that were required.		  */
/**************************************************************************/
void jacobi(aa,n,d,v,nrot)
float aa[DIM][DIM],d[DIM],v[DIM][DIM];
int n,*nrot;
{
int j,iq,ip,i;
float thresh,theta,tau,t,sm,s,h,g,c;
float b[DIM],z[DIM],a[DIM][DIM];

if(n!=DIM) 
	{ 
	fprintf(stderr,"\nFatal error: n not DIM %d in jacobi\n",DIM); 
	exit(1); 
	}
for(ip=0;ip<n;ip++) /* Initialize to the identity matrix */
	{
	for(iq=0;iq<n;iq++) v[ip][iq] = 0.0;
	for(iq=0;iq<n;iq++) a[ip][iq] = aa[ip][iq]; /* Don't destroy aa */
	v[ip][ip] = 1.0;
	}
/* Initialize b and d to the diagonals of a */
for(ip=0;ip<n;ip++)
	{
	b[ip] = d[ip] = a[ip][ip];
	z[ip] = 0.0;
	}
*nrot = 0;
for(i=0;i<100;i++)
	{
	sm = 0.0;
	for(ip=0;ip<(n-1);ip++)
		{
		for(iq=ip+1;iq<n;iq++)
			sm += fabs(a[ip][iq]);
		}

	/* Normal return, which relies on quadratic convergence to
	   machine underflow */
	if(sm == 0.0) return;

	if(i<3) thresh=0.2*sm/(n*n); /* on the first three sweeps */
	else thresh = 0.0; /* the rest of the sweeps */

	for(ip=0;ip<(n-1);ip++)
		{
		for(iq=ip+1;iq<n;iq++)
			{
			g = 100.0*fabs(a[ip][iq]);
			/* After 4 sweeps skip the rotation if the
			   off diagonal element is small */
			if(i>3 && fabs(d[ip])+g == fabs(d[ip])
			       && fabs(d[iq])+g == fabs(d[iq])) a[ip][iq] = 0.0;
			else if(fabs(a[ip][iq]) > thresh)
				{
				h = d[iq]-d[ip];
				if(fabs(h)+g==fabs(h)) t=(a[ip][iq])/h;
				else 
				  {
				  theta = 0.5*h/(a[ip][iq]);
				  t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
				  if(theta < 0.0) t = -t;
				  }
				c = 1.0/sqrt(1.0+t*t);
				s = t*c;
				tau = s/(1.0+c);
				h = t*a[ip][iq];
				z[ip] -= h;
				z[iq] += h;
				d[ip] -= h;
				d[iq] += h;
				a[ip][iq] = 0.0;
				for(j=0;j<ip-1;j++)
					rotate(a,j,ip,j,iq,&h,&g,s,tau);
				for(j=ip+1;j<iq-1;j++)
					rotate(a,ip,j,j,iq,&h,&g,s,tau);
				for(j=iq+1;j<n;j++)
					rotate(a,ip,j,iq,j,&h,&g,s,tau);
				for(j=0;j<n;j++)
					rotate(v,j,ip,j,iq,&h,&g,s,tau);
				++(*nrot);
				}
			}
		}
	for(ip=0;ip<n;ip++)
		{
		b[ip] += z[ip];
		d[ip] = b[ip];
		z[ip] = 0.0;
		}
	}
/* fprintf(stderr,"\nFatal error: too many iterations in jacobi\n"); */
}

/*********************************************************************/
/* Do rotations required by Jacobi Transformation		     */
/*********************************************************************/
rotate(a,i,j,k,l,h,g,s,tau)
float a[DIM][DIM],s,tau;
int i,j,k,l;
float *h,*g;
{
(*g) = a[i][j];
(*h) = a[k][l];
a[i][j] = (*g)-s*((*h)+(*g)*tau);
a[k][l] = (*h)+s*((*g)-(*h)*tau);
}


/*********************************************************************/
/* Compute the normalized mean and convariance matrix		     */
/* Step1: Recovery of Conservation Information.			     */
/*********************************************************************/
calc_mean_and_covariance1(calc,mean,covariance,u_low,u_high,v_low,v_high,
	eigenvalues,eigenvector1,eigenvector2,vel_type,N)
float calc[ARRSIZE][ARRSIZE]; 
float mean[2],covariance[2][2],u_low,v_low,u_high,v_high,eigenvalues[2];
float eigenvector1[2],eigenvector2[2];
int *vel_type,N;
{
int i,j,k,Q,fdr,a,b,c,d,size;
float inc[2],f_sum,uv[2],min,max,vx,vy,af,bf,est[2],temp[2][2],u,v,av,bv;
unsigned char raster[ARRSIZE][ARRSIZE],name[100],path[100],command[100];

/* Compute normalized weighted mean */ 
(*vel_type) = FULL;
max = -HUGE;
min = HUGE;
f_sum = 0.0;
size = 2*N+1;
for(i=0;i<size;i++)
for(j=0;j<size;j++)
	{
	f_sum += calc[i][j];
	if(calc[i][j] > max) { max = calc[i][j]; a=j; b=i; }
	if(calc[i][j] < min) { min = calc[i][j]; c=j; d=i; }
	}

inc[0] = (u_high-u_low)/((size-1)*1.0);
inc[1] = (v_high-v_low)/((size-1)*1.0);
uv[0] = u_low;
uv[1] = v_low;

est[0] = uv[0]+(inc[0])*a;
est[1] = uv[1]+(inc[1])*b;

/* Compute weighted means */
af = bf = 0.0;
for(i=0;i<size;i++)
for(j=0;j<size;j++)
	{
	bv = uv[1]+inc[1]*i;
	av = uv[0]+inc[1]*j;
	bf += bv*calc[i][j];
	af += av*calc[i][j];
	}
af = af/f_sum;
bf = bf/f_sum;
mean[0] = af;
mean[1] = bf;

/* Compute normalized weighted covariance matrix */
for(i=0;i<2;i++)
for(j=0;j<2;j++)
	temp[i][j] = covariance[i][j] = 0.0;

for(i=0;i<size;i++)
for(j=0;j<size;j++)
	{
	bv = uv[1]+inc[1]*i;
	av = uv[0]+inc[0]*j;
	temp[0][0] = covariance[0][0] += (bv-bf)*(bv-bf)*calc[i][j];
	temp[1][1] = covariance[1][1] += (av-af)*(av-af)*calc[i][j];
	temp[0][1] = covariance[0][1] += (bv-bf)*(av-af)*calc[i][j];
	}
temp[1][0] = covariance[1][0] = covariance[0][1];
for(i=0;i<2;i++)
for(j=0;j<2;j++)
	{
	temp[i][j] = (covariance[i][j] /= f_sum);
	}
comp_eigen(temp,eigenvalues,eigenvector1,eigenvector2);
(*vel_type) = FULL;
}

/************************************************************************/
/* Compute eigenvalues and eigenvectors of the covariance matrix	*/
/************************************************************************/
comp_eigen(A,value,vector1,vector2)
float A[DIM][DIM],value[DIM],vector1[DIM],vector2[DIM];
{
float v[2][2],temp,diff1[2],diff2[2],length1,length2,angle;
float eigenvalues2[2],eigenvectors2[2][2];
int nrot,SWAP;

SWAP = FALSE;
jacobi(A,2,value,v,&nrot);
/* if(value[0] < 0.0) 
	{ value[0] = -value[0]; v[0][0] = -v[0][0]; v[1][0] = -v[1][0]; }
if(value[1] < 0.0) 
	{ value[1] = -value[1]; v[0][1] = -v[0][1]; v[1][1] = -v[1][1]; } */
/* The largest eigenvalue should be the first */
if(fabs(value[0]) > fabs(value[1]))
	{
	vector1[0] = v[0][0]; vector1[1] = v[1][0];
	vector2[0] = v[0][1]; vector2[1] = v[1][1];
	}
else
	{
	temp = value[0];
	value[0] = value[1];
	value[1] = temp;
	vector2[0] = v[0][0]; vector2[1] = v[1][0];
	vector1[0] = v[0][1]; vector1[1] = v[1][1];
	temp = v[0][0]; v[0][0] = v[0][1]; v[0][1] = temp;
	temp = v[1][0]; v[1][0] = v[1][1]; v[1][1] = temp;
	SWAP = TRUE;
	}

/* Check eigenvalue/eigenvector calculation */
if(check_eigen_calc(A,value,v,nrot,diff1,diff2,&length1,&length2,&angle)==FALSE)
	{
	printf("\n********************************************\n");
	printf("Fatal error: eigenvalue/eigenvector error\n");
	printf("eigenvalues: %f %f\n",value[0],value[1]);
	printf("eigenvector1: %f %f\n",v[0][0],v[1][0]);
	printf("eigenvector2: %f %f\n",v[0][1],v[1][1]);
	if(SWAP) printf("Eigenvalues/eigenvectors are swapped from original order\n");
	printf("\n      A:\n");
	printf("%12.6f %12.6f\n",A[0][0],A[0][1]);
	printf("%12.6f %12.6f\n",A[1][0],A[1][1]);
	printf("Determinant of A: %f\n",A[0][0]*A[1][1]-A[0][1]*A[1][0]);
	printf("nrot: %d\n",nrot);
	printf("Angle between two eigenvectors: %f degrees\n",angle);
	printf("Difference length for eigenvector1\n");
	printf("Difference: %f %f Length: %f\n",diff1[0],diff1[1],length1);
	printf("Difference length for eigenvector2\n");
	printf("Difference: %f %f Length: %f\n",diff2[0],diff2[1],length2);

	/* Check eigenvalues/eigenvectors for 2*2 matrix using
	   Anandan's method outlined in his thesis. */
	eigenvalues2[0] = 0.5*((A[0][0]+A[1][1]) - 
	  sqrt((A[0][0]-A[1][1])*(A[0][0]-A[1][1])+4.0*A[0][1]*A[1][0]));
	eigenvalues2[1] = 0.5*((A[0][0]+A[1][1]) +
	  sqrt((A[0][0]-A[1][1])*(A[0][0]-A[1][1])+4.0*A[0][1]*A[1][0]));
	angle = atan2(eigenvalues2[0]-A[0][0],A[0][1]);
	eigenvectors2[0][0] = -cos(angle);
	eigenvectors2[1][0] = -sin(angle);
	eigenvectors2[0][1] = -sin(angle);
	eigenvectors2[1][1] = cos(angle);
	printf("\nUsing Anandan's calculation:\n");
	printf("Angle of rotation: %f degrees\n",angle*180/M_PI);
	printf("eigenvalues: %f %f\n",eigenvalues2[0],eigenvalues2[1]);
	printf("eigenvector1: %f %f\n",eigenvectors2[0][0],eigenvectors2[1][0]);
	printf("eigenvector2: %f %f\n",eigenvectors2[0][1],eigenvectors2[1][1]);
	check_eigen_calc(A,eigenvalues2,eigenvectors2,nrot,diff1,diff2,
			 &length1,&length2,&angle);
	printf("Angle between two eigenvectors: %f degrees\n",angle);
	printf("Difference length for eigenvector1\n");
	printf("Difference: %f %f Length: %f\n",diff1[0],diff1[1],length1);
	printf("Difference length for eigenvector2\n");
	printf("Difference: %f %f Length: %f\n",diff2[0],diff2[1],length2);
	printf("Original eigenvalues set to infinity\n");
	value[0] = value[1] = HUGE;
	printf("********************************************\n\n");
	fflush(stdout);
	if(FALSE) exit(1);
	}
}



/*********************************************************************/
/* Check eigenvector and eigenvalue computation for 2*2 matrix	     */
/*********************************************************************/
int check_eigen_calc(mm,d,v,n,diff1,diff2,length1,length2,angle)
float mm[2][2],d[2],v[2][2],diff1[2],diff2[2],*angle,*length1,*length2;
int n;
{
int status;

status = TRUE;
/* Compute angle between two eigenvectors - should be orthogonal */
(*angle)=acos((v[0][0]*v[0][1]+v[1][0]*v[1][1])/
	 (sqrt(v[0][0]*v[0][0]+v[1][0]*v[1][0])*
	  sqrt(v[0][1]*v[0][1]+v[1][1]*v[1][1])))*180.0/M_PI;
if((*angle) < 89.5 && (*angle) > 90.5) 
	{
	status = FALSE;
	}

/* Eigenvector test */
diff1[0] = mm[0][0]*v[0][0]+mm[0][1]*v[1][0];
diff1[1] = mm[1][0]*v[0][0]+mm[1][1]*v[1][0];
diff1[0] = diff1[0] - d[0]*v[0][0];
diff1[1] = diff1[1] - d[0]*v[1][0];
if(((*length1)=sqrt(diff1[0]*diff1[0]+diff1[1]*diff1[1])) > 0.1) 
	{
	status = FALSE;
	}
diff2[0] = mm[0][0]*v[0][1]+mm[0][1]*v[1][1];
diff2[1] = mm[1][0]*v[0][1]+mm[1][1]*v[1][1];
diff2[0] = diff2[0] - d[1]*v[0][1];
diff2[1] = diff2[1] - d[1]*v[1][1];
if(((*length2)=sqrt(diff2[0]*diff2[0]+diff2[1]*diff2[1])) > 0.1) 
	{
	status = FALSE;
	}
if(n > 50) 
	{
	status = FALSE;
	}
return(status);
}

/*****************************************************************/
/* Compute the SSD surface for point (a,b)			 */
/*****************************************************************/
compute_SSD_surface(fpic,x,y,n,N,SSDdata,max_a,max_b)
fimages fpic;
int x,y,n,N,*max_a,*max_b;
float SSDdata[ARRSIZE][ARRSIZE];
{
int u,v,u_index,v_index,i,j,a,b,fdr,pixel,size,c,d;
float term1,term2,k,sum1,sum2,constant,sum_min,sum_max,max,min;
float SSDvalues[MAX_ARRSIZE][MAX_ARRSIZE];
float SSDvalues1[MAX_ARRSIZE][MAX_ARRSIZE];
float SSDvalues2[MAX_ARRSIZE][MAX_ARRSIZE];
float SSDmin,SSDmax,SSDmag;
unsigned char raster[MAX_ARRSIZE][MAX_ARRSIZE];
unsigned path[100],name[100],command[100];

constant = 0.95;
size = 2*N;

/****************************************************************/
/* Compute the SSD surfaces for left and right images about x,y */
/* in a 4N+1 bu 4N+1 area					*/
/****************************************************************/
for(u=(-size);u<=size;u++)
for(v=(-size);v<=size;v++)
	{
	u_index = u+size;
	v_index = v+size;
	sum1 = 0.0;
	sum2 = 0.0;
	for(i=(-n);i<=n;i++)
	for(j=(-n);j<=n;j++)
		{
		term1 = (fpic[1][x+i][y+j]-fpic[2][x+i+u][y+j+v]);
		term2 = (fpic[1][x+i][y+j]-fpic[0][x+i+u][y+j+v]); 
		sum1 += term1*term1;
		sum2 += term2*term2; 
		}
	SSDvalues1[u_index][v_index] = sum1;
	SSDvalues2[2*size-u_index][2*size-v_index] = sum2; 
	}


/****************************************************************************/
/* Add the left and right SSD surfaces and compute minimum SSD coordinates  */
/* in 2N+1 by 2N+1 area centered about the origin of 4N+1 by 4N+1 SSD data  */
/****************************************************************************/
sum_min =  HUGE;
sum_max = -HUGE;
SSDmin =  HUGE;
SSDmax = -HUGE;
a = b = 0;
SSDmag = 0.0;
for(u=(-size);u<=size;u++)
for(v=(-size);v<=size;v++)
	{
	u_index = u+size;
	v_index = v+size;
	SSDvalues[u_index][v_index] = SSDvalues1[u_index][v_index]+
				      SSDvalues2[u_index][v_index];
	/* Only look for peak in 2N+1 by 2N+1 region centered at origin */
	if(u >= -N && u <= N && v >= -N && v <= N)
	{
	if((fabs(SSDvalues[u_index][v_index]-sum_min) <= SMALL && (u*u+v*v) < SSDmag) ||
	   ((sum_min-SSDvalues[u_index][v_index]) > SMALL))
		{ 
		sum_min=SSDvalues[u_index][v_index]; 
		a=u; b=v; 
		SSDmag = (u*u+v*v);
		}
	if(SSDvalues[u_index][v_index] > sum_max) 
		{ sum_max=SSDvalues[u_index][v_index]; }
	}
	if(SSDvalues[u_index][v_index] > SSDmax) SSDmax = SSDvalues[u_index][v_index];
	if(SSDvalues[u_index][v_index] < SSDmin) SSDmin = SSDvalues[u_index][v_index];
	}
(*max_a) = a;
(*max_b) = b;

/* A test: compute the locations of n=10 minimum SSD responses */
if(FALSE)
{
if((x>=PT_X1-2 && x<= PT_X1+2&& y>=PT_Y1-2 && y<=PT_Y1+2) ||
   (x>=PT_X2-2 && x<= PT_X2+2&& y>=PT_Y2-2 && y<=PT_Y2+2))
{
printf("\nBig Minimum data at image point: %d,%d\n",x,y);
printf("max_a: %d max_b: %d\n",*max_a,*max_b);
compute_big_n_minimums(SSDvalues,size,20); 
}
}

/*******************************************/
/* Center SSD computation by extracting    */
/* 2N+1 by 2N+1 SSD sub-surface about peak */
/*******************************************/
sum_min = HUGE;
c = d = 0;
for(u=(-N);u<=N;u++)
for(v=(-N);v<=N;v++)
	{
	u_index = u+N;
	v_index = v+N;
	SSDdata[u_index][v_index] = SSDvalues[u_index+a+N][v_index+b+N];
	/* Compute minimum non-zero SSD value */
	if(SSDdata[u_index][v_index] < sum_min && 
	   SSDdata[u_index][v_index] != 0.0) 
		{ sum_min=SSDdata[u_index][v_index]; c=u_index; d=v_index;}
	}

if(((x>=PT_X1-2 && x<= PT_X1+2&& y>=PT_Y1-2 && y<=PT_Y1+2) ||
    (x>=PT_X2-2 && x<= PT_X2+2&& y>=PT_Y2-2 && y<=PT_Y2+2)) && FALSE)
	printf("c=%d d=%d u=%d v=%d\n",c,d,c-N,d-N);
/* A test: compute the locations of n=10 minimum SSD responses */
if(FALSE)
{
if((x>=PT_X1-2 && x<= PT_X1+2&& y>=PT_Y1-2 && y<=PT_Y1+2) ||
   (x>=PT_X2-2 && x<= PT_X2+2&& y>=PT_Y2-2 && y<=PT_Y2+2))
{
printf("\nMinimum data at image point: %d,%d about: %d %d\n",x,y,a,b);
compute_n_minimums(SSDdata,N,20); 
}
}
/**********************************************************************/
/* Compute a k value and recompute the SSD surface using it.	      */
/**********************************************************************/
if(sum_min < 1.0e-6) k = 0.0085;
else k = -log(constant)/sum_min;

min =  HUGE;
max = -HUGE;
for(u=(-N);u<=N;u++)
for(v=(-N);v<=N;v++)
	{
	u_index = u+N;
	v_index = v+N;
	SSDdata[u_index][v_index] = exp(-k*SSDdata[u_index][v_index]);
	if(SSDdata[u_index][v_index] > max) max = SSDdata[u_index][v_index];
	if(SSDdata[u_index][v_index] < min) min = SSDdata[u_index][v_index];
	}

} 

/*********************************************************************/ 
/* Compute the normalized mean and convariance matrix		     */ 
/* Step2: Recovery of Neighbourhood Information.		     */ 
/*********************************************************************/
calc_mean_and_covariance2(weights,velocities,weighted_size,mean,covariance)
float weights[WEIGHTSIZE][WEIGHTSIZE],velocities[WEIGHTSIZE][WEIGHTSIZE][2]; 
float mean[2],covariance[2][2];
int weighted_size;
{
int i,j,k,Q,fdr,a,b;
float vx,vy,w_sum;

/* Compute weighted means - sum of weights is 1.0 */
w_sum = vx = vy = 0.0;
for(i=0;i<weighted_size;i++)
for(j=0;j<weighted_size;j++)
	{
	if(velocities[i][j][0] != NO_VALUE && velocities[i][j][1] != NO_VALUE)
		{
		w_sum += weights[i][j];
		vx += weights[i][j]*velocities[i][j][0];
		vy += weights[i][j]*velocities[i][j][1];
		}
	}
mean[0] = vx;
mean[1] = vy;

/* Compute normalized weighted covariance matrix */
for(i=0;i<2;i++)
for(j=0;j<2;j++)
	covariance[i][j] = 0.0;

/* Sum of the weights is 1.0 so no need to normalize */
for(i=0;i<weighted_size;i++)
for(j=0;j<weighted_size;j++)
	{
	if(velocities[i][j][0] != NO_VALUE && velocities[i][j][1] != NO_VALUE)
	{
	covariance[0][0] += (velocities[i][j][0]-vx)*(velocities[i][j][0]-vx)*weights[i][j];
	covariance[1][1] += (velocities[i][j][1]-vy)*(velocities[i][j][1]-vy)*weights[i][j];
	covariance[0][1] += (velocities[i][j][0]-vx)*(velocities[i][j][1]-vy)*weights[i][j];
	}
	}
covariance[0][0] /= w_sum;
covariance[0][1] /= w_sum;
covariance[1][1] /= w_sum;
covariance[1][0] = covariance[0][1];
}


/**********************************************************************/
/* Invert 2*2 matrix if possible				      */
/**********************************************************************/
int old_inverse22(M,MI)
float M[2][2],MI[2][2];
{
float denominator;

/* Invert 2*2 matrix */
denominator = M[0][0]*M[1][1]-M[1][0]*M[0][1]; /* The determinant of M */
if(fabs(denominator) >= 0.000000001)
{
MI[0][0] = M[1][1]/denominator;
MI[0][1] = -M[0][1]/denominator;
MI[1][0] = -M[1][0]/denominator;
MI[1][1] = M[0][0]/denominator;
printf("            M                          MI\n");
printf("%12.8f %12.8f    %12.8f %12.8f\n",M[0][0],M[0][1],MI[0][0],MI[0][1]);
printf("%12.8f %12.8f    %12.8f %12.8f\n",M[1][0],M[1][1],MI[1][0],MI[1][1]);
printf("Determinant: %f\n",denominator);
return(FALSE); /* Not singular */
}
else
{
MI[0][0] = 1.0;
MI[0][1] = 0.0;
MI[1][0] = 0.0;
MI[1][1] = HUGE;
printf("            M                          MI\n");
printf("%12.8f %12.8f    %12.8f %12.8f\n",M[0][0],M[0][1],MI[0][0],MI[0][1]);
printf("%12.8f %12.8f    %12.8f %12.8f\n",M[1][0],M[1][1],MI[1][0],MI[1][1]);
printf("Determinant: %f\n",denominator);
return(TRUE); /* Singular */
}
}

/***********************************************************************/
/* Add two vectors, a and b, and place result in c		       */
/***********************************************************************/
add21(a,b,c)
float a[2],b[2],c[2];
{
c[0] = a[0] + b[0];
c[1] = a[1] + b[1];
}

/***********************************************************************/
/* Add two 2*2 matrices S1 and S2 and save the result in 2*2 matrix s3 */
/***********************************************************************/
add22(S1,S2,S3)
float S1[2][2],S2[2][2],S3[2][2];
{
S3[0][0] = S1[0][0] + S2[0][0];
S3[1][0] = S1[1][0] + S2[1][0];
S3[0][1] = S1[0][1] + S2[0][1];
S3[1][1] = S1[1][1] + S2[1][1];
}

/***********************************************************************/
/* Multiple a 2*2 matrix and a 2*1 vector			       */
/***********************************************************************/
mult21(A,b,x)
float A[2][2],b[2],x[2];
{
x[0] = A[0][0]*b[0]+A[0][1]*b[1];
x[1] = A[1][0]*b[0]+A[1][1]*b[1];
}

/**********************************************************************/
/* Compute 2D Gaussian weights					      */
/**********************************************************************/
compute_weights(weights,w)
float weights[WEIGHTSIZE][WEIGHTSIZE];
int w;
{
int i,j,offset,weighted_size;
float kernel[WEIGHTSIZE],term1,term2,sum;

weighted_size = 2*w+1;

if(w==1)
{
kernel[0] = 0.25; kernel[1] = 0.5; kernel[2] = 0.25;
sum = 0.0;
for(i=0;i<(weighted_size);i++)
for(j=0;j<(weighted_size);j++)
	{
	weights[i][j] = kernel[i]*kernel[j];
	sum += weights[i][j];
	}
if(sum!=1.0)
	{
	printf("Fatal error: sum of Gaussian coefficients is zero\n");
	exit(1);
	}
}
else
if(w==2)
{
kernel[0] = 0.0625; kernel[1] = 0.25; kernel[2] = 0.375;
kernel[3] = 0.25; kernel[4] = 0.0625;
sum = 0.0;
for(i=0;i<(weighted_size);i++)
for(j=0;j<(weighted_size);j++)
	{
	weights[i][j] = kernel[i]*kernel[j];
	sum += weights[i][j];
	}
if(sum!=1.0)
	{
	printf("Fatal error: sum of Gaussian coefficients is zero\n");
	exit(1);
	}
}
else
{
printf("\nFatal error: window size not 1 or 2: w=%d\n",w);
exit(1);
}

printf("\nWeights:\n");
for(i=0;i<(weighted_size);i++)
{
for(j=0;j<(weighted_size);j++)
	printf("%8.5f ",weights[i][j]);
printf("\n");
}
}


/*************************************************************************/
/* Compute the pseudo-inverse of J using its SVD                         */
/*************************************************************************/
int inverse22(J,JI, print)
float J[NO_ROWS][NO_COLS], JI[NO_COLS][NO_ROWS];
int print;
{
extern void dsvdc(double [], int* , int*, int*, double [], double [], double [], int*, double [], int* , double [], int*, int*); 
int size;
float VT[NO_COLS][NO_COLS];
float U[NO_COLS][NO_ROWS],DI[NO_COLS][NO_COLS];
float UT[NO_ROWS][NO_COLS],V[NO_COLS][NO_COLS];
float D[NO_COLS][NO_COLS],I[NO_COLS][NO_COLS];
double JT[NO_COLS][NO_ROWS], W[NO_ROWS+1], zero[NO_COLS];
double UU[NO_COLS][NO_ROWS],temp[NO_ROWS],VV[NO_COLS][NO_COLS];
float min,max,cond,best,worst;
int i,j,k,m,n,mdim,error,Ierr,job;

size = NO_ROWS;
n = NO_COLS; m = NO_ROWS; mdim = NO_ROWS;
for(i=0;i<NO_COLS;i++)
for(j=0;j<size;j++)
	JT[i][j] = J[j][i];

job = 21;
dsvdc(JT,&mdim,&m,&n,W,zero,UU,&mdim,VV,&n,temp,&job,&Ierr); 

for(i=0;i<size;i++)
for(j=0;j<NO_COLS;j++)
	{
	U[j][i] = UU[j][i];
	UT[i][j] = UU[j][i];
	}
for(i=0;i<NO_COLS;i++)
for(j=0;j<NO_COLS;j++)
	{
	V[i][j] = VV[i][j];
	D[i][j] = DI[i][j] = 0.0;
	}
for(i=0;i<NO_COLS;i++)
for(j=0;j<NO_COLS;j++)
	VT[i][j] = V[j][i];

for(i=0;i<NO_COLS;i++) D[i][i] = W[i];
min = W[0];
max = W[0];
for(i=1;i<NO_COLS;i++)
	{
	if(W[i] > max) max=W[i];
	if(W[i] < min) min=W[i];
	}
if(min != 0.0) { cond = max/min; worst = 1.0/min; }
else { cond = HUGE; worst = HUGE; }
best = 1.0/max;

for(i=0;i<NO_COLS;i++) if(D[i][i] !=0.0) DI[i][i] = 1.0/D[i][i];
	else DI[i][i] = HUGE;
if(print) printf("\nIerr=%d DI[0][0]: %f DI[1][1]: %f\n",Ierr,DI[0][0],DI[1][1]);

/* Check correctness of SVD, compute I = JI*J */
if(DI[0][0] <= 10.0 && DI[1][1] <= 10.0)
{
for(i=0;i<NO_COLS;i++)
for(j=0;j<NO_COLS;j++)
	{
	V[i][j] = 0.0;
	for(k=0;k<NO_COLS;k++) V[i][j] = V[i][j] + VT[i][k]*DI[k][j];
	}
for(i=0;i<NO_COLS;i++)
for(j=0;j<size;j++)
	{
	JI[i][j] = 0.0;
	for(k=0;k<NO_COLS;k++) JI[i][j] = JI[i][j] + V[i][k]*U[k][j];
	}

for(i=0;i<NO_COLS;i++)
for(j=0;j<NO_COLS;j++)
	{
	I[i][j] = 0.0;
	for(k=0;k<size;k++)
		I[i][j] = I[i][j] + JI[i][k]*J[k][j];
	}
/* Check I */
error = FALSE;
for(i=0;i<NO_COLS;i++)
for(j=0;j<NO_COLS;j++)
	{
	if(i==j) 
		{
		if(fabs(I[i][i])<0.9999 || fabs(I[i][i])>1.0001) error=TRUE;
		}
	else if(fabs(I[i][j]) > 0.0001) error=TRUE;
	}
if(error==TRUE)
	{
	printf("\nFatal error: Pseudo-inverse of J wrong\n");
	printf("The identity matrix:\n");
	for(i=0;i<NO_COLS;i++)
		{
		for(j=0;j<NO_COLS;j++) printf("%6.3f ",I[i][j]);
		printf("\n");
		}
	fflush(stdout);
	exit(1);
	}
}
else /* Set inverse elements over 10 to 10 and compute inverse */
{
if(DI[0][0] > 10.0) DI[0][0] = 10.0;
if(DI[1][1] > 10.0) DI[1][1] = 10.0; 
COUNT_SINGULAR++;
/*
JI[0][0] = 10.0; JI[0][1] = 0.0;
JI[1][0] = 0.0; JI[1][1] = 10.0;
*/

if(TRUE)
{
for(i=0;i<NO_COLS;i++)
for(j=0;j<NO_COLS;j++)
	{
	VV[i][j] = 0.0;
	for(k=0;k<NO_COLS;k++) VV[i][j] = VV[i][j] + VT[i][k]*DI[k][j];
	}
for(i=0;i<NO_COLS;i++)
for(j=0;j<size;j++)
	{
	JI[i][j] = 0.0;
	for(k=0;k<NO_COLS;k++) JI[i][j] = JI[i][j] + VV[i][k]*U[k][j];
	}
}
}
if(print)
	{
	printf("Diagonal(s) greater than 10\n");
	printf("       J                      JI\n");
	printf("%12.8f %12.8f     %12.8f %12.8f\n",J[0][0],J[0][1],JI[0][0],JI[0][1]);
	printf("%12.8f %12.8f     %12.8f %12.8f\n",J[1][0],J[1][1],JI[1][0],JI[1][1]);
	}

return(FALSE); /* The inverse matrix is guaranteed to be non-singular */
}

/************************************************************
  Convolve images with center surround filter, i.e. a DOG
  (Difference of Gaussians) 
************************************************************/
laplacian(cpic,fpic,sigma,pic_x,pic_y)
int pic_x,pic_y;
cimages cpic;
fimages fpic;
float sigma;
{
int i,j,k,l,m,OFFSET,count,no_count;
float sigma2,begin,end;
float val2,pic2[PIC_Y][PIC_X];
float h2[TOTAL],const2,x1,x2;
float pic22[PIC_Y][PIC_X],min,max;

/* Use the original image as first image and subtract a second
   Gaussian blurred image (with sigma2 from it */
sigma2 = sigma;
const2 = 1.0/(sqrt(2.0*M_PI)*sigma2);
OFFSET = (int) (6*sigma2+1)/2.0;
for(i=0,x1=(-OFFSET),x2=(-OFFSET);i<=2*OFFSET;i++,x1++,x2++) 
	{
	h2[i] = const2 * exp((-x2*x2)/(2.0*sigma2*sigma2));
	}

/* Do convolutions using separable 2D filters */
printf("\nComputing Laplacian as a center-surround filter...\n");
printf("sigma1=%f sigma2=%f\n",0.0,sigma2);
for(i=0;i<NUMFILES; i++)
	{
	printf("\n\n");
        for(j=0;j<pic_y;j++)
		{
                for(k=OFFSET;k<pic_x-OFFSET;k++)
			{
			val2 = 0.0;
	        	for(l=(-OFFSET);l<=OFFSET;l++) 
				{
              			val2 += cpic[i][j][k+l]*h2[l+OFFSET];
				} 	 
       			pic2[j][k] = val2;
			}
		}  
	printf("Image #%d - column convolution completed\n",i);
	fflush(stdout);
        for(j=OFFSET;j<pic_y-OFFSET;j++)
		{
                for(k=OFFSET;k<pic_x-OFFSET;k++)
			{
			val2 = 0.0;
	        	for(m=(-OFFSET);m<=OFFSET;m++) 
				{
              			val2 += pic2[j+m][k]*h2[m+OFFSET];
				} 	 
       			pic22[j][k] = val2;
			}
		}  
	printf("Image #%d - row convolution completed\n",i);
	fflush(stdout);

	/* Subtract blurred images to obtain Laplacian of Gaussian */
	min =  HUGE;
	max = -HUGE;
	for(j=OFFSET;j<pic_y-OFFSET;j++)
	for(k=OFFSET;k<pic_x-OFFSET;k++)
		{
		fpic[i][j][k] = (((float) cpic[i][j][k]) - pic22[j][k]);
		if(fpic[i][j][k] < min) min=fpic[i][j][k];
		if(fpic[i][j][k] > max) max=fpic[i][j][k];
		}
	printf("Gaussian images subtracted\n");
	printf("Minimum value: %f Maximun value: %f\n",min,max);
	fflush(stdout);

	for(j=OFFSET;j<pic_y-OFFSET;j++)
	for(k=OFFSET;k<pic_x-OFFSET;k++)
		cpic[i][j][k] = ((int) (fpic[i][j][k]-min)/(max-min)*255);
	}
}

/*************************************************************************/
/* Make the image grayvalues stored in cpic as unsigned characters into  */
/* float numbers stored in fpic.					 */
/*************************************************************************/
make_float(cpic,fpic,pic_x,pic_y)
int pic_x,pic_y;
cimages cpic;
fimages fpic;
{
int i,j,k;
for(k=0;k<3;k++)
for(i=0;i<pic_y;i++)
for(j=0;j<pic_x;j++)
	fpic[k][i][j] = (float) cpic[k][i][j];
}

/*************************************************************************/
/* Compute the n maximum data values and their locations in the array.   */
/*************************************************************************/
compute_n_minimums(data,N,n)
int N,n;
float data[ARRSIZE][ARRSIZE];
{
int i,j,u,v,u_index,v_index,u_min,v_min;
float min_value;
float SSD[ARRSIZE][ARRSIZE],SSDmag;

/* Copy SSD data so we can modify it without affecting the SSD data */
for(i=0;i<2*N+1;i++)
{
for(j=0;j<2*N+1;j++)
	{
	SSD[i][j] = data[i][j];
	}
}
for(i=0;i<n;i++)
{
min_value = HUGE;
SSDmag = HUGE;
u_min = v_min = 0;
for(u=(-N);u<=N;u++)
for(v=(-N);v<=N;v++)
	{
	u_index = u+N;
	v_index = v+N;
	if((fabs(SSD[u_index][v_index]-min_value) <= SMALL && (u*u+v*v < SSDmag)) || SSD[u_index][v_index] < min_value)
		{
		min_value = SSD[u_index][v_index];
		SSDmag = u*u+v*v;
		u_min = u;
		v_min = v;
		}
	}
SSD[u_min+N][v_min+N] = HUGE;
printf("%dth minimum: %f at %d,%d\n",i,min_value,u_min,v_min);
}
}

/*************************************************************************/
/* Compute the n maximum data values and their locations in the array.   */
/*************************************************************************/
compute_big_n_minimums(data,N,n)
int N,n;
float data[MAX_ARRSIZE][MAX_ARRSIZE];
{
int i,j,u,v,u_index,v_index,u_min,v_min;
float min_value;
float SSD[MAX_ARRSIZE][MAX_ARRSIZE],SSDmag;

/* Copy SSD data so we can modify it without affecting the SSD data */
for(i=0;i<2*N+1;i++)
{
for(j=0;j<2*N+1;j++)
	{
	SSD[i][j] = data[i][j];
	}
}
for(i=0;i<n;i++)
{
min_value = HUGE;
SSDmag = 0.0;
u_min = v_min = 0;
for(u=(-N);u<=N;u++)
for(v=(-N);v<=N;v++)
	{
	u_index = u+N;
	v_index = v+N;
	if((fabs(SSD[u_index][v_index]-min_value) <= SMALL && (u*u+v*v < SSDmag)) || 
	   (min_value-SSD[u_index][v_index] > SMALL))
		{
		min_value = SSD[u_index][v_index];
		SSDmag = u*u+v*v;
		u_min = u;
		v_min = v;
		}
	}
SSD[u_min+N][v_min+N] = HUGE;
printf("%dth minimum: %f at %d,%d SSDmag: %d\n",i,min_value,u_min,v_min,u_min*u_min+v_min*v_min);
}
}

float bigL2norm(bigdiff,n)
float bigdiff[PIC_X][PIC_Y][2];
int n;
{
int i,j,ct;
float sum;

sum = 0.0;
ct = 0;
for(i=0;i<n;i++)
for(j=0;j<n;j++)
	{
	if(bigdiff[i][j][0] != NO_VALUE && bigdiff[i][j][1] != NO_VALUE)
		{
		sum += bigdiff[i][j][1]*bigdiff[i][j][1]+bigdiff[i][j][0]*bigdiff[i][j][0];
		ct++;
		}
	}
return(sqrt(sum));
}
