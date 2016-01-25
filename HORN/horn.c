/************************************************************
   horn.c
   Travis Burkitt, modified by John Barron, 1992
   May  1988  -- written for MASSCOMP
    -- May 8, 1989  modified for SUN 3/60
    --- Implementation of Horn and Schunck's algorithm
	for calculating an optical flow field.
************************************************************/
#include  <fcntl.h>
#include  <stdio.h>
#include  <math.h>
#include "rasterfile.h"


#define BORDER 2
#define HEAD 32
#define PI  M_PI
#define PIC_X 675
#define PIC_Y 675
#define PIC_T 15
#define FIVE 5
#define PMODE 0644
#define TRUE 1
#define FALSE 0
#define HEAD 32
#define NO_VALUE 100.0
#define OUTPUT_SMOOTH 1
#define BIG_MAG 1000000.0
#define HUGE 999999999.9


unsigned char inpic[PIC_T][PIC_X][PIC_Y];
unsigned char pic[FIVE][PIC_X][PIC_Y];
float floatpic[FIVE][PIC_X][PIC_Y];
unsigned header[HEAD];
float Ix[PIC_X][PIC_Y],Iy[PIC_X][PIC_Y];
float It[PIC_X][PIC_Y],full_vels[PIC_X][PIC_Y][2];
float correct_vels[PIC_X][PIC_Y][2],full_vels1[PIC_X][PIC_Y][2]; 
float diff_x(),diff_y(),diff_t(),difference(),temp_vels[PIC_X][PIC_Y][2];
float PsiER(),norm(),alpha,fmin();
int pic_x,pic_y,pic_t,THRESHOLD,STANDARD,BINARY,int_size_x,int_size_y;
float actual_x,actual_y,size_x,size_y,offset_x,offset_y;
int startx,starty,endx,endy,step,WRITE_SMOOTH;


/*********************************************************************/
/*   Main program 						     */
/*********************************************************************/
main(argc,argv)
int argc;
char **argv;
{
int fd,fdf,offset,size,start,end,middle,i,j,num;
int fd_correct,no_bytes,numpass,time;
float ave_error,st_dev,density,min_angle,max_angle,sigma,tau;
unsigned char header[HEAD],path1[100],path2[100],path3[100];
char full_name[100],norm_name[100],correct_filename[100];

if(argc < 7 || argc > 19)
	{
        printf("Usage: %s <filename stem> <alpha> <sigma> <central file number> <number of iterations> <input path> <output path> [-S <smooth path> -C <full correct filename> -B <cols> <rows> -T <tau> -H or -MH]\n\n",argv[0]);
	printf("<filename stem> - stem of input filenames\n");
	printf("<alpha> - Lagrange multiplier value\n");
	printf("<sigma> - the standard deviation used in smoothing the files\n");
	printf("   sigma==0.0 means no smoothing: use 2 or 5 unsmoothed images\n");
	printf("<central file number> - the image file number for which\n");
	printf("   image velocity is to be computed\n");
	printf("<number of iterations> - number of iterations in the velocity calculation\n");
	printf("<input path> - directory where input data resides\n");
	printf("<output path> - directory where computed flow fields put\n");
	printf("-S <smooth path> - directory where smoothed data is to be put\n");
	printf("   if not present smoothed files are not written\n"); 
	printf("-B <cols> <rows> - use binary file of size <cols>*<rows> characters\n");
	printf("   instead of black and white rasterfiles as image input\n");
	printf("   image size read from rasterfile header if used\n");
	printf("-C <correct filename> - perform error analysis using correct velocity field\n");
	printf("-T <tau> - threshold the computed velocities on spatial intensity gradient\n");
	printf("-H Standard Horn and Schunck\n");
	printf("   perform H and S differencing on 2 input images\n");
	printf(" sigma is ignored for -H\n");
	printf("-MH Non-standard Horn and Schunck - default\n");
	printf("    perform 4-point central differences on 5 input images\n");
	exit(1);
        }

printf("\n------------------------------\n");
printf("Command line: ");
for(i=0;i<argc;i++) printf("%s ",argv[i]);
printf("\n%d arguments\n",argc-1);

sscanf(argv[2],"%f",&alpha);
sscanf(argv[3],"%f",&sigma);
sscanf(argv[4],"%d",&middle);
sscanf(argv[5],"%d",&numpass);
printf("alpha=%f\n",alpha);
printf("sigma=%f\n",sigma);
printf("Central image: %d\n",middle);
printf("Number of iterations: %d\n",numpass);
strcpy(path1,argv[6]); 
strcpy(path2,"."); 
strcpy(path3,argv[7]); 

printf("Input directory: %s\n",path1);
printf("Output directory: %s\n",path3);

BINARY = FALSE;
STANDARD = FALSE;
THRESHOLD = FALSE;
WRITE_SMOOTH = FALSE;
i=8;
strcpy(correct_filename,"unknown");
while(i<argc)
{
if(strcmp("-H",argv[i])==0) 
	{ 
	STANDARD = TRUE; 
	i++;
	}
else
if(strcmp("-MH",argv[i])==0) 
	{ 
	STANDARD = FALSE; 
	i++;
	}
else
if(strcmp("-C",argv[i])==0) 
	{ 
	strcpy(correct_filename,argv[i+1]); 
	i+=2;
	}
else
if(strcmp("-T",argv[i])==0) 
	{ 
	sscanf(argv[i+1],"%f",&tau); 
	i+=2;
	THRESHOLD = TRUE;
	}
else
if(strcmp("-S",argv[i])==0)
	{
	strcpy(path2,argv[i+1]);
	WRITE_SMOOTH = TRUE;
	i+=2;
	}
else
if(strcmp("-B",argv[i])==0) 
	{
	sscanf(argv[i+1],"%d",&pic_y);
	sscanf(argv[i+2],"%d",&pic_x);
	BINARY = TRUE;
	i += 3;
	}
else
	{
	printf("Invalid option %s specified as argument %d -  program terminates\n",argv[i],i);
	exit(1);
	}
}

if(WRITE_SMOOTH) printf("Smoothed data directory: %s\n",path2);
else printf("Smoothed images not written\n");
fflush(stdout);

printf("Correct velocity file: %s\n",correct_filename);
fflush(stdout);
if(strcmp(correct_filename,"unknown")!=0)
{
/* Read the correct velocity data */
if((fd_correct = open(correct_filename,O_RDONLY))==NULL)
	{
	printf("Fatal error in opening file %s\n",correct_filename);
	printf("fd_correct: %d\n",fd_correct);
	exit(1);
	}
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
	no_bytes += read(fd_correct,&correct_vels[i][0][0],int_size_y*8);
printf("\nFile %s opened and read\n",correct_filename);
printf("Size of correct velocity data: %d %d\n",int_size_y,int_size_x);
printf("%d bytes read\n",no_bytes);
fflush(stdout);
}

if(!STANDARD)
{
size = 6*sigma+1;
if(size%2==0) size = size+1;
offset = size/2+2; /* Add 2 as neighbourhood size offset */
start = middle-offset;
end = middle+offset;
printf("Size: %d Offset: %d Start: %d End: %d\n",size,offset,start,end);
printf("%d images required\n",end-start+1);
if((end-start+1) > PIC_T)
	{
	printf("Fatal error: not enough room for the images\n");
	exit(1);
	}
if(end < start) 
	{ 
	printf("Fatal error: specify images in ascending order\n"); 
	exit(1);
	}

read_and_smooth3D(path1,argv[1],sigma,floatpic,pic,inpic,start,middle,end,header);
if(sigma!=0.0 && WRITE_SMOOTH)
writefiles(path2,argv[1],pic,sigma,pic_t,pic_x,pic_y,middle-2,middle+2,header);
}
else 
{
start = middle;
end = middle+1;
offset = 2;
readfiles(path1,argv[1],inpic,pic_t,&pic_x,&pic_y,start,end,header);
/* Copy the input images unchanged */
for(i=0;i<pic_x;i++)
for(j=0;j<pic_y;j++)
	{
	pic[0][i][j] = inpic[0][i][j];
	pic[1][i][j] = inpic[1][i][j];
	pic[2][i][j] = 0;
	pic[3][i][j] = 0;
	pic[4][i][j] = 0;
	floatpic[0][i][j] = (int) inpic[0][i][j];
	floatpic[1][i][j] = (int) inpic[1][i][j];
	floatpic[2][i][j] = 0.0;
	floatpic[3][i][j] = 0.0;
	floatpic[4][i][j] = 0.0;
	}
}

printf("Number of Columns: %d Number of Rows: %d\n",pic_y,pic_x);
if(pic_x > PIC_X || pic_y > PIC_Y)
	{ 
	printf("Fatal error: images are too big\n");
	exit(1);
	}

printf("\n");
fflush(stdout);
if(STANDARD==TRUE)
{
if(THRESHOLD) sprintf(full_name,"%s/horn.original.%sF-%4.2f",path3,argv[1],tau);
else sprintf(full_name,"%s/horn.original.%sF",path3,argv[1]);
}
else
{
if(THRESHOLD) sprintf(full_name,"%s/horn.modified.%sF-%4.2f",path3,argv[1],tau);
else sprintf(full_name,"%s/horn.modified.%sF",path3,argv[1]);
}

printf("Output full velocities go to file %s\n",full_name);
printf("\n");
fflush(stdout);

startx = 0;
starty = 0;
endx = pic_x-1;
endy = pic_y-1;

startx += BORDER;
starty += BORDER;
endx -= BORDER;
endy -= BORDER;

num = 2;
step = 1;

time = 0;
if(STANDARD)
	{
	calcIx(Ix,floatpic,time);
	calcIy(Iy,floatpic,time);
	calcIt(It,floatpic,time);
	}
else compute_ders(Ix,Iy,It,floatpic,pic_t,pic_x,pic_y,offset);
if(FALSE)
{
if(STANDARD) printf("Standard Horn and Schunck Derivatives\n");
else printf("4 point central difference derivatives\n");
print_ders(Ix,Iy,It);
}

for(i=0;i<PIC_X;i++)
for(j=0;j<PIC_Y;j++)
	{
	full_vels[i][j][0] = full_vels[i][j][1] = 0.0;
	full_vels1[i][j][0] = full_vels1[i][j][1] = 0.0;
	}

/* Perform iterations */
for(i=0;i<numpass;i+=2) 
	{
	printf("%3dth iteration\n",i);
	fflush(stdout);
	calc_vels(full_vels,full_vels1,Ix,Iy,It);
	printf("The improvement: %f\n",difference(full_vels,full_vels1,pic_x,pic_y));
	fflush(stdout);
if(strcmp(correct_filename,"unknown")!=0)
	{
	rearrange(full_vels1,temp_vels);
	calc_statistics(correct_vels,int_size_x,int_size_y,temp_vels,
     		pic_x,pic_y,2*offset,&ave_error,&st_dev,&density,&min_angle,&max_angle);
	printf("Error: %f St Dev: %f Density: %f\n",ave_error,st_dev,density);
	fflush(stdout);
	}
	printf("%3dth iteration\n",i+1);
	calc_vels(full_vels1,full_vels,Ix,Iy,It);
	printf("The improvement: %f\n",difference(full_vels,full_vels1,pic_x,pic_y));
if(strcmp(correct_filename,"unknown")!=0)
	{
	rearrange(full_vels,temp_vels);
	calc_statistics(correct_vels,int_size_x,int_size_y,temp_vels,
     		pic_x,pic_y,2*offset,&ave_error,&st_dev,&density,&min_angle,&max_angle);
	printf("Error: %f St Dev: %f Density: %f\n",ave_error,st_dev,density);
	fflush(stdout);
	}
	}
if(THRESHOLD) threshold(full_vels,Ix,Iy,tau,pic_x,pic_y);

if((fdf=creat(full_name,0644))!=NULL)
	output_velocities(fdf,"Full",full_vels,pic_x,pic_y,2*offset);
else printf("Error in opening %s file\n\n",full_name);
printf("\nVelocities computed and output\n");
fflush(stdout);

if(strcmp(correct_filename,"unknown")!=0)
{
calc_statistics(correct_vels,int_size_x,int_size_y,full_vels,
     pic_x,pic_y,2*offset,&ave_error,&st_dev,&density,&min_angle,&max_angle);
printf("\n\nComputed Statistics\n");
printf("Error: %f St Dev: %f\n",ave_error,st_dev);
printf("Density: %f\n",density);
printf("Minimum angle error: %f Maximum angle error: %f\n",min_angle,max_angle);
}
fflush(stdout);
}



/****************************************************************/
/* Compute x derivatives					*/
/****************************************************************/
calcIx(Ex,floatpic,t)
float Ex[PIC_X][PIC_Y];
float floatpic[FIVE][PIC_X][PIC_Y];
int t;
{
int i,j;

printf("****** calculating Ex ******\n"); 
for(i=startx;i<=endx;i+=step)
for(j=starty;j<=endy;j+=step) 
	{
	Ex[i][j] = (floatpic[t][i+step][j] + floatpic[t][i+step][j+step] +
	    	    floatpic[t+1][i+step][j] + floatpic[t+1][i+step][j+step])/4.0
		  -(floatpic[t][i][j] + floatpic[t][i][j+step] +
		    floatpic[t+1][i][j] + floatpic[t+1][i][j+step])/4.0;
	if(Ex[i][j] > BIG_MAG)
		{
		printf("Ex too large at i=%d j=%d Ex=%f\n",i,j,Ex[i][j]);
		exit(1);
		}
	}
}


/****************************************************************/
/* Compute y derivatives					*/
/****************************************************************/
calcIy(Ey,floatpic,t)
float Ey[PIC_X][PIC_Y];
float floatpic[FIVE][PIC_X][PIC_Y];
int t;
{
int i,j;

printf("****** calculating Ey ******\n");
for(i=startx;i<=endx;i+=step)
for(j=starty;j<=endy;j+=step) 
	{
	Ey[i][j] = (floatpic[t][i][j+step] + floatpic[t][i+step][j+step] + 
		    floatpic[t+1][i][j+step] + floatpic[t+1][i+step][j+step])/4.0
		  -(floatpic[t][i][j] + floatpic[t][i+step][j] +
		    floatpic[t+1][i][j] + floatpic[t+1][i+step][j])/4.0;
	if(Ey[i][j] > BIG_MAG)
		{
		printf("Ey too large at i=%d j=%d Ey=%f\n",i,j,Ey[i][j]);
		exit(1);
		}
	}
}

/****************************************************************/
/* Compute t derivatives					*/
/****************************************************************/
calcIt(Et,floatpic,t)
float Et[PIC_X][PIC_Y];
float floatpic[FIVE][PIC_X][PIC_Y];
int t;
{
int i,j;

printf("****** calculating Et ******\n"); 
for(i=startx;i<=endx;i+=step)
for(j=starty;j<=endy;j+=step)
	{
	Et[i][j] = (floatpic[t+1][i][j] + floatpic[t+1][i+step][j] +
		    floatpic[t+1][i][j+step] + floatpic[t+1][i+step][j+step])/4.0
		  -(floatpic[t][i][j] + floatpic[t][i+step][j] +
	      	    floatpic[t][i][j+step] + floatpic[t][i+step][j+step])/4.0;
	if(Et[i][j] > BIG_MAG)
		{
		printf("Et too large at i=%d j=%d Ex=%f\n",i,j,Et[i][j]);
		exit(1);
		}
	}
}


/****************************************************************/
/* Compute average u value in neighbour about i,j		*/
/****************************************************************/
vels_avg(vels,ave)
float vels[PIC_X][PIC_Y][2],ave[PIC_X][PIC_Y][2];
{
int i,j;
for(i=startx+1;i<endx;i++)
for(j=starty+1;j<endy;j++)
	{
	ave[i][j][0] = (vels[i-1][j][0]+vels[i][j+1][0]+
			vels[i+1][j][0]+vels[i][j-1][0])/6.0 +
       		       (vels[i-1][j-1][0]+vels[i-1][j+1][0]+
		        vels[i+1][j+1][0]+vels[i+1][j-1][0])/12.0;
	ave[i][j][1] = (vels[i-1][j][1]+vels[i][j+1][1]+
			vels[i+1][j][1]+vels[i][j-1][1])/6.0 +
       		       (vels[i-1][j-1][1]+vels[i-1][j+1][1]+
			vels[i+1][j+1][1]+vels[i+1][j-1][1])/12.0;
	}

/* Copy average values of neighbourhoods for boundaries */
for(i=startx;i<=endx;i++)
	{
	ave[i][starty][0] = ave[i][starty+1][0];
	ave[i][endy][0] = ave[i][endy-1][0];
	ave[i][starty][1] = ave[i][starty+1][1];
	ave[i][endy][1] = ave[i][endy-1][1];
	}
for(j=starty+1;j<=endy;j++)
	{
	ave[startx][j][0] = ave[startx+1][j][0];
	ave[endx][j][0] = ave[endx-1][j][0];
	ave[startx][j][1] = ave[startx+1][j][1];
	ave[endx][j][1] = ave[endx-1][j][1];
	}
/* Corner Points */
ave[startx][0][0] = ave[startx+1][1][0];
ave[startx][0][1] = ave[startx+1][1][1];
ave[startx][endy][0] = ave[startx+1][endy-1][0];
ave[startx][endy][1] = ave[startx+1][endy-1][1];
ave[endx][0][0] = ave[endx-1][1][0];
ave[endx][0][1] = ave[endx-1][1][1];
ave[endx][endy][0] = ave[endx-1][endy-1][0];
ave[endx][endy][1] = ave[endx-1][endy-1][1];
}



/****************************************************************/
/* Compute u,v values						*/
/****************************************************************/
calc_vels(vels,vels1,Ex,Ey,Et)
float vels[PIC_X][PIC_Y][2],vels1[PIC_X][PIC_Y][2];
float Ex[PIC_X][PIC_Y],Ey[PIC_X][PIC_Y],Et[PIC_X][PIC_Y];
{
int i,j,k;
float mag,ave[PIC_X][PIC_Y][2];

printf("****** Computing Velocity ******\n"); 
fflush(stdout);
vels_avg(vels1,ave);
for(i=startx;i<=endx;i+=step)
for(j=starty;j<=endy;j+=step) 
	{
	vels[i][j][0] = ave[i][j][0]-Ex[i][j]*
		(Ex[i][j]*ave[i][j][0]+Ey[i][j]*ave[i][j][1]+Et[i][j])
		/(alpha*alpha+Ex[i][j]*Ex[i][j]+Ey[i][j]*Ey[i][j]);
	vels[i][j][1] = ave[i][j][1]-Ey[i][j]*
		(Ex[i][j]*ave[i][j][0]+Ey[i][j]*ave[i][j][1]+Et[i][j])
		/(alpha*alpha+Ex[i][j]*Ex[i][j]+Ey[i][j]*Ey[i][j]);
	mag = sqrt(vels[i][j][0]*vels[i][j][0]+vels[i][j][1]*vels[i][j][1]);
	if(mag > 5.0 && FALSE) 
		{
		printf("Velocity magnitude of %f at %d %d is over 5.0\n",mag,i,j) ;
                }
	}
}



/****************************************************************/
/* Print derivative information					*/
/****************************************************************/
print_ders(Ex,Ey,Et)
float Ex[PIC_X][PIC_Y],Ey[PIC_X][PIC_Y],Et[PIC_X][PIC_Y];
{
int i,j;

printf("********************  Ex  ********************\n");
for(i=50;i<=60;i+=step) 
	{
	for(j=50;j<=60;j+=step) printf("%5.1f ",Ex[i][j]);
	printf ("\n");
    	}

printf("********************  Ey  ********************\n");
for(i=50;i<=60;i+=step) 
	{
	for(j=50;j<=60;j+=step) printf("%5.1f ",Ey[i][j]);
	printf("\n");
    	}

printf("********************  Et  ********************\n");
for(i=50;i<=60;i+=step) 
	{
	for(j=50;j<=60;j+=step) printf("%5.1f ",Et[i][j]);
	printf("\n");
    	}
}




/*********************************************************************/
/* Compute spatio-temporal derivatives 				     */
/*********************************************************************/
compute_ders(Ix,Iy,It,floatpic,pic_t,pic_x,pic_y,n)
float Ix[PIC_X][PIC_Y];
float Iy[PIC_X][PIC_Y];
float It[PIC_X][PIC_Y];
float floatpic[PIC_T][PIC_X][PIC_Y];
int n,pic_t,pic_x,pic_y;
{
int i,j;
float kernel[5];

for(i=0;i<PIC_X;i++)
for(j=0;j<PIC_X;j++)
	{
	Ix[i][j] = Iy[i][j] = It[i][j] = 0.0;
	}

calc_diff_kernel(kernel);
for(i=n;i<pic_x-n;i++)
for(j=n;j<pic_y-n;j++)
	{
	Ix[i][j] = diff_x(floatpic,kernel,i,j,2);
	Iy[i][j] = diff_y(floatpic,kernel,i,j,2);
	It[i][j] = diff_t(floatpic,kernel,i,j,2);
	/* printf("i:%d j:%d Ix:%f Iy:%f It:%f\n",i,j,Ix[i][j],Iy[i][j],It[i][j]); */
	}
printf("Spatio-Temporal Intensity Derivatives Computed\n");
fflush(stdout);
}

/*********************************************************************/
/* Compute a 4 point central difference kernel			     */
/*********************************************************************/
calc_diff_kernel(diff_kernel)
float diff_kernel[5];
{
diff_kernel[0] = -1.0/12.0;
diff_kernel[1] = 8.0/12.0;
diff_kernel[2] = 0.0;
diff_kernel[3] = -8.0/12.0;
diff_kernel[4] = 1.0/12.0;
}


/************************************************************
   Apply 1D real kernels in the x direction
************************************************************/
float diff_x(floatpic,kernel,x,y,n)
float floatpic[PIC_T][PIC_X][PIC_Y];
float kernel[5];
int x,y,n;
{
int i;
float sum;

sum = 0.0;
for(i=(-n);i<=n;i++)
	{
	sum += kernel[i+2]*floatpic[2][x+i][y];
	}
return(sum);
}

/************************************************************
   Apply 1D real kernels in the y direction
************************************************************/
float diff_y(floatpic,kernel,x,y,n)
float floatpic[PIC_T][PIC_X][PIC_Y];
float kernel[5];
int x,y,n;
{
int i;
float sum;

sum = 0.0;
for(i=(-n);i<=n;i++)
	{
	sum += kernel[i+2]*floatpic[2][x][y+i];
	}
return(sum);
}


/************************************************************
   Apply 1D real kernels in the t direction.
************************************************************/
float diff_t(floatpic,kernel,x,y,n)
float floatpic[PIC_T][PIC_X][PIC_Y];
float kernel[5];
int x,y,n;
{
int i;
float sum;

sum = 0.0;
for(i=(-n);i<=n;i++)
	{
	sum += kernel[i+2]*floatpic[2+i][x][y];
	}
return(sum);
}


/*********************************************************************/
/* Read input images and perform Gaussian smoothing.		     */
/*********************************************************************/
readfiles(path,s,pic,pic_t,pic_x,pic_y,start,end,header)
char s[100],path[100];
int pic_t,*pic_x,*pic_y,start,end;
unsigned char header[HEAD];
unsigned char pic[PIC_T][PIC_X][PIC_Y];
{
char fname[100];
int i,j,k,fp,fd,time,no_bytes;
unsigned char temp_save[PIC_X][PIC_Y];
int ONCE;
int ints[8];


printf("Reading Files...\n");
time = -1;
if((end-start) < 0)
	{
	printf("\nSpecified time for writing file incorrect\n");
	exit(1);
	}
ONCE = TRUE;
for(i=start;i<=end;i++) 
	{
	time++;
	no_bytes = 0;
	sprintf(fname,"%s/%s%d",path,s,i);
 	if((fp=open(fname,O_RDONLY)) >0)
		{
		if(!BINARY)
		{
		if(ONCE)
		{
		no_bytes += read(fp,ints,HEAD);
		(*pic_y) = ints[1];
		(*pic_x) = ints[2];
		ONCE = FALSE;
		}
		else no_bytes += read(fp,header,HEAD);
		}

		for(j=0;j<(*pic_x);j++)
			no_bytes += read(fp,&pic[time][j][0],(*pic_y));
		printf("File %s read (%d bytes)\n",fname,no_bytes);
		no_bytes = 0;
		}
	else 
		{
		printf("File %s does not exist in readfiles.\n",fname);
		exit(1);
		}
	fflush(stdout);
	}
} /* End of readfiles */


/*********************************************************************/
/*   Write smoothed files					     */
/*********************************************************************/
writefiles(path,s,result,sigma,pic_t,pic_x,pic_y,start,end,header)
char s[100],path[100];
float sigma;
int pic_t,pic_x,pic_y,start,end;
unsigned char result[FIVE][PIC_X][PIC_Y];
unsigned char header[HEAD];
{
char fname[100];
int i,j,k,fp,time,no_bytes;

printf("\nWriting smoothed files...\n");

time = -1;
for (i=start;i<=end;i++) 
	{
	no_bytes = 0;
	time++;
	sprintf(fname,"%s/smoothed.%s%d-%3.1f",path,s,i,sigma);
 	if((fp=creat(fname,0644))!=NULL)
		{
		no_bytes += write(fp,&header[0],HEAD); /* Write 32 byte raster header */
		for(j=0;j<pic_x;j++)
			no_bytes += write(fp,&result[time][j][0],pic_y);
		printf("File %s written (%d bytes)\n",fname,no_bytes);
		}
	else 
		{
		printf("File %s cannot be written\n",fname);
		exit(1);
		}
	fflush(stdout);
	close(fp);
	}
} /* End of writefiles */


/*********************************************************************/
/* Smooth the image sequence using a 3D separable Gaussian filter.   */
/*********************************************************************/
read_and_smooth3D(path,stem,sigma,floatpic,pic,inpic,start,middle,end,header)
char stem[100],path[100];
float sigma;
unsigned char pic[FIVE][PIC_X][PIC_Y];
unsigned char inpic[PIC_T][PIC_X][PIC_Y];
unsigned char header[HEAD];
float floatpic[FIVE][PIC_X][PIC_Y];
int start,end,middle;
{
int fd,n,size,i,j,k,time,frame;
char name[100];

for(k=0;k<FIVE;k++)
for(i=0;i<PIC_X;i++)
for(j=0;j<PIC_Y;j++)
	floatpic[k][i][j] = 0.0;

pic_t = end-start+1;


readfiles(path,stem,inpic,pic_t,&pic_x,&pic_y,start,end,header);

printf("Number of columns: %d Number of rows: %d\n",pic_y,pic_x);
if(pic_x > PIC_X || pic_y > PIC_Y)
	{
	printf("Fatal error: images are too big\n");
	exit(1);
	}
fflush(stdout);

time = -1;
frame = middle-2;
if(OUTPUT_SMOOTH) printf("Start: %d End: %d Middle: %d\n",start,end,middle);
for(frame=middle-2;frame<=middle+2;frame++)
	{
	time++;
	convolve_Gaussian(inpic,floatpic,pic,sigma,pic_t,pic_x,pic_y,
		start,middle-2+time,time);
	}
if(sigma!=0.0) printf("Input files smoothed\n");
}

/************************************************************************/
/* Perform 3D Gaussian smoothing by separable convolution 		*/
/************************************************************************/
convolve_Gaussian(inpic,floatpic,pic,sigma,pic_t,pic_x,pic_y,start,frame,time)
unsigned char inpic[PIC_T][PIC_X][PIC_Y];
unsigned char pic[FIVE][PIC_X][PIC_Y];
float floatpic[FIVE][PIC_X][PIC_Y];
int pic_x,pic_y,pic_t,frame,time,start;
float sigma;
{
float mask[100],term,product,sum;
int size,i,j,k,offset,a,b;
float pic0[PIC_X][PIC_Y],pic1[PIC_X][PIC_Y];


if(OUTPUT_SMOOTH)
{
printf("\nStart of 3D convolution\n");
printf("Time: %d Frame: %d\n",time,frame);
}

fflush(stdout);
size = (int) 6*sigma+1;
if(size%2==0) size = size+1;
offset = size/2;
if(pic_t < size)
	{
	printf("\nFatal error: not enough images\n");
	exit(1);
	}
sum = 0.0;

if(sigma != 0.0)
for(i=0;i<size;i++)
	{
	mask[i] = (1.0/(sqrt(2.0*3.141592654)*sigma))*
		exp(-(i-offset)*(i-offset)/(2.0*sigma*sigma));
	sum = sum+mask[i];
	}
else { mask[0] = 1.0; sum = 1.0; }

if(OUTPUT_SMOOTH && sigma!=0.0)
{
printf("Size: %d Offset: %d\nMask values: ",size,offset); 
for(i=0;i<size;i++)
	printf("%f ",mask[i]);
printf("\nSum of mask values: %f\n",sum); 
}
for(i=0;i<pic_x;i++)
for(j=0;j<pic_y;j++)
	pic0[i][j] = pic1[i][j] = 0.0;


if(sigma != 0.0)
{
for(i=0;i<pic_x;i++)
for(j=0;j<pic_y;j++)
	{
	term = 0.0;
	for(a=-offset;a<=offset;a++)
		{
		term = term +(inpic[a+frame-start][i][j]*mask[a+offset]);
		}
	pic0[i][j] = term;
	}
if(OUTPUT_SMOOTH) printf("Convolution in t direction completed\n");
for(i=offset;i<pic_x-offset;i++)
for(j=offset;j<pic_y-offset;j++)
	{
	term = 0.0;
	for(a=-offset;a<=offset;a++)
		{
		term = term + (pic0[i+a][j]*mask[a+offset]);
		}
	pic1[i][j] = term;
	}
if(OUTPUT_SMOOTH) printf("Convolution in x direction completed\n");

for(i=offset;i<pic_x-offset;i++)
for(j=offset;j<pic_y-offset;j++)
	{
	term = 0.0;
	for(b=-offset;b<=offset;b++)
		{
		term = term + (pic1[i][j+b])*mask[b+offset];
		}
	floatpic[time][i][j] = term;
	if(term > 255.0) term = 255.0;
	if(term < 0.0) term = 0.0;
	pic[time][i][j] = (int) (term+0.5);
	}
if(OUTPUT_SMOOTH) printf("Convolution in y direction completed\n");
}
else /* No smoothing */
{
printf("No smoothing: frame=%d frame-start=%d time=%d\n",frame,frame-start,time);
fflush(stdout);
for(i=0;i<pic_x;i++)
for(j=0;j<pic_y;j++)
	{
	pic[time][i][j] = inpic[frame-start][i][j];
	floatpic[time][i][j] = inpic[frame-start][i][j];
	}
}

if(OUTPUT_SMOOTH) printf("End of Convolution\n");
fflush(stdout);
}


/************************************************************
   Output full velocities using old Burkitt format
************************************************************/
output_velocities(fdf,s,full_velocities,pic_x,pic_y,n)
float full_velocities[PIC_X][PIC_Y][2];
char s[100];
int fdf,n;
{
float x,y;
int i,j,bytes,no_novals,no_vals,NORMAL;

NORMAL = FALSE;
if(strcmp(s,"Normal")==0) NORMAL = TRUE;
if(fdf==NULL)
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
y = pic_x-2*n;
x = pic_y-2*n;
write(fdf,&x,4);
write(fdf,&y,4);

/* offset to start of data */
y = n;
x = n;
write(fdf,&x,4);
write(fdf,&y,4);
bytes = 24;

no_novals = no_vals = 0;
/* Prepare velocities for output, i.e. rotate by 90 degrees */
for(i=n;i<pic_x-n;i++)
for(j=n;j<pic_y-n;j++)
	{
	if(full_velocities[i][j][0] != NO_VALUE && 
	   full_velocities[i][j][1] != NO_VALUE)
		{
		no_vals++;
		x = full_velocities[i][j][0];
		y = full_velocities[i][j][1];
		full_velocities[i][j][0] = y;
		full_velocities[i][j][1] = -x;
		}
	else
		{
		no_novals++;
		}
	}
for(i=n;i<pic_x-n;i++)
	{
	bytes += write(fdf,&full_velocities[i][n][0],(pic_y-2*n)*8);
	}
close(fdf);
printf("\n%s velocities output from output_velocities: %d bytes\n",s,bytes);
printf("Number of positions with velocity: %d\n",no_vals);
printf("Number of positions without velocity: %d\n",no_novals);
printf("Percentage of %s velocities: %f\n",s,
	no_vals/(1.0*(no_vals+no_novals))*100.0);
fflush(stdout);
}


/***************************************************************/
/*  Compute error statistics				       */
/***************************************************************/
calc_statistics(correct_vels,int_size_x,int_size_y,full_vels,
	pic_x,pic_y,n,ave_error,st_dev,density,min_angle,max_angle)
float full_vels[PIC_X][PIC_Y][2],*ave_error,*density,*st_dev;
float correct_vels[PIC_X][PIC_Y][2],*min_angle,*max_angle;
int n,pic_x,pic_y,int_size_x,int_size_y;
{
int full_count,no_full_count,i,j,a,b,total_count;
float sumX2,temp,uva[2],uve[2];

full_count = no_full_count = total_count = 0;
sumX2 = 0.0;
(*min_angle) = HUGE;
(*max_angle) = -HUGE;
(*ave_error) = (*st_dev) = (*density) = 0.0;
for(i=n;i<pic_x-n;i++)
{
for(j=n;j<pic_y-n;j++)
	{
	if(full_vels[i][j][0] != NO_VALUE && full_vels[i][j][1] != NO_VALUE)
	  {
	  full_count++;
	  uve[0] = full_vels[i][j][0]; uve[1] = full_vels[i][j][1];
	  uva[0] = correct_vels[i][j][0]; uva[1] = correct_vels[i][j][1];
	  temp = PsiER(uve,uva);
	  (*ave_error) += temp;
	  sumX2 += temp*temp;
	  if(temp < (*min_angle)) (*min_angle) = temp;
	  if(temp > (*max_angle)) (*max_angle) = temp;
	  }
	else no_full_count++;
	total_count++;
	}
}

if(full_count != 0) (*ave_error) = (*ave_error)/full_count;
else (*ave_error) = 0.0;
if(full_count > 1) 
{
temp = fabs((sumX2 - full_count*(*ave_error)*(*ave_error))/(full_count-1));
(*st_dev) = sqrt(temp);
}
else (*st_dev) = 0.0;
(*density) = full_count*100.0/(total_count*1.0);

if((*ave_error) == 0.0) { (*min_angle) = (*max_angle) = 0.0; }

if(FALSE)
{
printf("\nIn calc_statistics\n");
printf("%d full velocities\n",full_count);
printf("%d positons without full velocity\n",no_full_count);
printf("%d positions in total\n",total_count);
fflush(stdout);
}
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

nva = norm(VA,3);
nve = norm(VE,3);
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
   returns the norm of a vector v of length n.
************************************************************/
float norm(v,n)
float v[];
int n;
{
int i;
float sum = 0.0;

for (i=0;i<n; i++)
        sum += (v[i]*v[i]);
sum = sqrt(sum);
return sum;
}

/**********************************************************/
/* Compute the difference between two flow fields	  */
/**********************************************************/
float difference(v1,v2,pic_x,pic_y)
float v1[PIC_X][PIC_Y][2],v2[PIC_X][PIC_Y][2];
int pic_x,pic_y;
{
int i,j,n;
float sum,t1,t2;

sum = 0.0;
n = 0;
for(i=0;i<pic_x;i++)
for(j=0;j<pic_y;j++)
	{
	t1 = v1[i][j][0]-v2[i][j][0];
	t2 = v1[i][j][1]-v2[i][j][1];
	sum += sqrt(t1*t1+t2*t2);
	n++;
	}
sum = sum/n;
return(sum);
}


/******************************************************************/
/* Threshold the computed velocities on the basis of the spatial  */
/* intensity gradient.						  */
/******************************************************************/
threshold(full_vels,Ix,Iy,tau,pic_x,pic_y)
float full_vels[PIC_X][PIC_Y][2],Ix[PIC_X][PIC_Y],Iy[PIC_X][PIC_Y];
float tau;
int pic_x,pic_y;
{
int i,j,count;
float gradient2,tau2;

count = 0;
tau2 = tau*tau;
for(i=0;i<pic_x;i++)
for(j=0;j<pic_y;j++)
	{
	gradient2 = Ix[i][j]*Ix[i][j]+Iy[i][j]*Iy[i][j]; 
	/* gradient2 = fmin(Ix[i][j]*Ix[i][j],Iy[i][j]*Iy[i][j]); */
	if(gradient2 < tau2)
		{
		count++;
		full_vels[i][j][0] = NO_VALUE;
		full_vels[i][j][1] = NO_VALUE;
		}
	}
printf("Threshold: %f\n",tau);
printf("%d velocities thresholded\n",count);
}

float fmin(x,y)
float x,y;
{
if(x<y) return(x);
else return(y);
}

/*************************************************************/
/* Rearrange velocity field v1 (in output format) and place  */
/* the result in v2 so that error analysis can be performed. */
/*************************************************************/
rearrange(v1,v2)
float v1[PIC_X][PIC_Y][2],v2[PIC_X][PIC_Y][2];
{
int i,j;
for(i=0;i<PIC_X;i++)
for(j=0;j<PIC_Y;j++)
	{
	if(v1[i][j][0] != NO_VALUE && v1[i][j][1] != NO_VALUE)
		{
		v2[i][j][0] =  v1[i][j][1];
		v2[i][j][1] = -v1[i][j][0];
		}
	else
		{
		v2[i][j][0] = NO_VALUE;
		v2[i][j][1] = NO_VALUE;
		}
	}
}
