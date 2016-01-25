/*********************************************************************
	Lucas and Kanade, 1981
	With non-hierarchical modifications using Simoncelli, 
		Adelson and Heeger, CVPR 91
	cc -g lucas.c -lm -o lucas
**********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "rasterfile.h"

#define N 2
#define HALF_C 3.14159 
#define PIC_X 512
#define PIC_Y 512
#define PIC_T 17
#define FIVE 5
#define PMODE 0644
#define TRUE 1
#define FALSE 0
#define HEAD 32
#define NO_VALUE 100.0
#define MAX_I 25
#define DEBUG 0


unsigned char inpic[PIC_T][PIC_X][PIC_Y];
unsigned char pic[FIVE][PIC_X][PIC_Y];
float floatpic[FIVE][PIC_X][PIC_Y];
unsigned char header[HEAD];
float Ix[PIC_X][PIC_Y],Iy[PIC_X][PIC_Y],E[PIC_X][PIC_Y];
float It[PIC_X][PIC_Y],full_vels[PIC_X][PIC_Y][2],Imag[PIC_X][PIC_Y];
float norm_vels1[PIC_X][PIC_Y][2],correct_vels[PIC_X][PIC_Y][2]; 
float norm_vels2[PIC_X][PIC_Y][2];
float diff1[2],diff2[2];
int pic_x,pic_y,pic_t;
int BINARY,int_size_x,int_size_y,RAW_STATISTICS,OUTPUT_SMOOTH;
float actual_x,actual_y,size_x,size_y,offset_x,offset_y,RAW_MAG;
float MAX_RAW_MAG;
FILE *fd_temp[8];
float data[8][150][150];
float pic0[PIC_X][PIC_Y],pic1[PIC_X][PIC_Y];

void calc_statistics(float correct_vels[PIC_X][PIC_Y][2],
				float norm_vels1[PIC_X][PIC_Y][2],
				float norm_vels2[PIC_X][PIC_Y][2],
				int int_size_x, int int_size_y,
				float full_vels[PIC_X][PIC_Y][2],
				float E[PIC_X][PIC_Y], int pic_x,
				int pic_y, int n, float *ave_error,
				float *st_dev,float *density, float *residual,
				float *min_angle, float *max_angle,
				float *norm_ave_error1,float *norm_st_dev1,
				float *norm_density1, float *norm_min_angle1,
				float *norm_max_angle1, float *norm_ave_error2,
				float *norm_st_dev2, float *norm_density2,
				float *norm_min_angle2, float *norm_max_angle2);
int check_eigen_calc(float mm[N][N],
float d[N], float v[N][N], int n, 
float diff1[N], float diff2[N], float *length1, float *length2,
float *angle);
void compute_ders(
float Ix[PIC_X][PIC_Y],
float Iy[PIC_X][PIC_Y],
float It[PIC_X][PIC_Y],
float floatpic[PIC_T][PIC_X][PIC_Y],
int pic_t, int pic_x, int pic_y, int n, int middle,
char *filename);

void read_and_smooth3D(char path[], char stem[],
float sigma, float floatpic[][PIC_X][PIC_Y],
 unsigned char pic[][PIC_X][PIC_Y],
unsigned char inpic[][PIC_X][PIC_Y],
int start, int middle, int end,
unsigned char header[]);
void writefiles(char path[], char s[],
unsigned char result[][PIC_X][PIC_Y],
float sigma,
int pic_t, int pic_x, int pic_y, int start, int end,
unsigned char header[HEAD]);
void compute_vels(
float Ix[PIC_X][PIC_Y],
float Iy[PIC_X][PIC_Y],
float It[PIC_X][PIC_Y],
float full_vels[PIC_X][PIC_Y][2],
float norm_vels1[PIC_X][PIC_Y][2],
float norm_vels2[PIC_X][PIC_Y][2],
int pic_x, int pic_y, int n,
float tau_D,
int flag, float E[PIC_X][PIC_Y]);
void calc_diff_kernel(float diff_kernel[]);
void convolve_Gaussian(
unsigned char inpic[PIC_T][PIC_X][PIC_Y],
float floatpic[FIVE][PIC_X][PIC_Y],
unsigned char pic[FIVE][PIC_X][PIC_Y],
float sigma, int pic_t, int pic_x, int pic_y,
int start, int frame, int time);
float diff_t(
float floatpic[PIC_T][PIC_X][PIC_Y],
float kernel[5],
int x, int y, int n);
float diff_x(
float floatpic[PIC_T][PIC_X][PIC_Y],
float kernel[5],
int x, int y, int n);
float diff_y(
float floatpic[PIC_T][PIC_X][PIC_Y],
float kernel[5],
int x, int y, int n);
void jacobi(float aa[N][N], int n, float d[N], float v[N][N],
int *nrot);
float norm(float v[],int n);
void output_velocities(FILE *fdf, char s[],
float full_velocities[PIC_X][PIC_Y][2],
int pic_x, int pic_y,
int n);
float PsiEN(float ve[],float va[]);
float PsiER(float ve[], float va[]);
void readfiles(char path[], char s[],
unsigned char pic[PIC_T][PIC_X][PIC_Y],
float sigma, int pic_t, int *pic_x, int *pic_y, int start, int end,
unsigned char header[HEAD]);
void rotate(float a[N][N], int i, int j, int k, int l,
float *h, float *g, float s, float tau);


/*********************************************************************/
/*   Main program 						     */
/*********************************************************************/
int main(argc,argv)
int argc;
char **argv;
{

FILE * fdf, *fdn,*fdr;
int offset,size,start,end,middle,flag,i;
FILE *fd_correct;
int no_bytes;
float ave_error,st_dev,residual,density,min_angle,max_angle;
float norm_ave_error1,norm_st_dev1,norm_density1,norm_min_angle1,norm_max_angle1;
float norm_ave_error2,norm_st_dev2,norm_density2,norm_min_angle2,norm_max_angle2;
float sigma,tau_D;
unsigned char header[HEAD],path1[100],path2[100],path3[100];
char name[100],full_name[100],norm_name[100],raw_name[100],correct_filename[100];


if(argc < 6 || argc > 17)
	{
        printf("Usage: %s <filename stem> <sigma> <central file number> <threshold> <input path> <output path> [-S <smoothed path> -C <full correct filename> -B <cols> <rows>  -L -M]\n",argv[0]);
	printf("Read enough images so that smoothing with a 3D Gaussian\n");
	printf(" with <sigma> leaves 5 smoothed images centered at ");
	printf("<central file number>\n");
	printf("<input path> - directory where input data resides\n");
	printf("<output path> - directory where computed flow fields put\n");
	printf("sigma is the standard deviation used in smoothing the files\n");
	printf("sigma==0.0 means no smoothing: 5 unsmoothed files are used\n");
	printf("-L - standard Lucas and Kanade with presmoothing and thresholding\n");
	printf("-M - modified Lucas and Kanade (Simoncelli, Adelson and Heeger)\n");
	printf("default -L if any flag other than -M specified or no flag specified\n");
	printf("-B <cols> <rows> - use binary instead of black and white rasterfiles as image input\n");
	printf("-S <smoothed path> - creates appropriate files for the smoothed input data\n");
	printf("   and writes to the specifed directory\n");
	printf("-C <correct filename> - perform error analysis using correct velocity field\n");
	printf("-T <float> - specify a magnitude threshold for spatial derivatives\n");
	printf("             for the raw normal velocity data\n");
	printf("If -T 0.0 is specified histogram and cumulative statistics are printed\n");
	printf("If -T is not specified default value is 0.0 with no statistics printed\n");
	exit(1);
        }
printf("\n--------------------------------\nCommand line: ");
for(i=0;i<argc;i++) printf("%s ",argv[i]);
printf("\n%d arguments\n",argc-1);
sprintf(name,"testdata/%s",argv[1]);
sscanf(argv[2],"%f",&sigma);
sscanf(argv[3],"%d",&middle);
sscanf(argv[4],"%f",&tau_D);
strcpy(path1,".");
strcpy(path2,".");
strcpy(path3,".");
printf("sigma=%f\n",sigma);
printf("Central image: %d\n",middle);
printf("Tau threshold: %f\n",tau_D);
strcpy(path1,argv[5]); 
strcpy(path3,argv[6]); 
printf("Input directory: %s\n",path1);
printf("Output directory: %s\n",path3);
fflush(stdout);

if(DEBUG)
{
for(i=0;i<8;i++)
	{
	sprintf(name,"data.%d",i);
	fd_temp[i] = fopen(name,"wb");
	printf("File %s created\n",name);
	}
fflush(stdout);
}

BINARY = FALSE;
OUTPUT_SMOOTH = FALSE;
RAW_STATISTICS = FALSE;
RAW_MAG = 0.0;
flag = FALSE;
i= 7;
strcpy(correct_filename,"unknown");
while(i<argc)
{
if(strcmp("-M",argv[i])==0) { flag = TRUE; i++; }
else
if(strcmp("-L",argv[i])==0) i++;
else
if(strcmp("-C",argv[i])==0) { strcpy(correct_filename,argv[i+1]); i+=2;}
else
if(strcmp("-T",argv[i])==0)
	{
	sscanf(argv[i+1],"%f",&RAW_MAG);
	if(RAW_MAG==0.0) RAW_STATISTICS = TRUE;
	i += 2;
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
if(strcmp("-S",argv[i])==0)
	{
	strcpy(path2,argv[i+1]); 
	printf("Smoothed data directory: %s\n",path2);
	OUTPUT_SMOOTH = TRUE;
	i += 2;
	}
else
	{
	printf("Invalid option %s - program terminates\n",argv[i]);
	exit(1);
	}
}


printf("Correct velocity file: %s\n",correct_filename);
fflush(stdout);
if(strcmp(correct_filename,"unknown")!=0)
{
/* Read the correct velocity data */
if((fd_correct = fopen(correct_filename,"rb"))==NULL)
	{
	printf("Fatal error in opening file %s\n",correct_filename);
	printf("fd_correct: %d\n",fd_correct);
	exit(1);
	}
no_bytes = 0;
no_bytes += (4*fread(&actual_y,4, 1, fd_correct));
no_bytes += (4*fread(&actual_x,4, 1, fd_correct));
no_bytes += (4*fread(&size_y,4, 1, fd_correct));
no_bytes += (4*fread(&size_x,4,1, fd_correct));
no_bytes += (4*fread(&offset_y,4, 1, fd_correct));
no_bytes += (4*fread(&offset_x,4, 1, fd_correct));
if(offset_x != 0.0 || offset_y != 0.0 || actual_x != size_x || actual_y != size_y)
	{
	printf("Fatal error: something wrong with correct velocity data\n");
	printf("Actual y: %f Actual x: %f\n",actual_y,actual_x);
	printf("Size y: %f Size x: %f\n",size_y,size_x);
	printf("Offset y: %f Offset x: %f\n",offset_y,offset_x);
	exit(1);
	}
int_size_y = (int)size_y;
int_size_x = (int)size_x;
for(i=0;i<int_size_x;i++)
	no_bytes += (4 * fread(&correct_vels[i][0][0],4, int_size_y*2, fd_correct));
printf("\nFile %s opened and read\n",correct_filename);
printf("Size of correct velocity data: %d %d\n",int_size_y,int_size_x);
printf("%d bytes read\n",no_bytes);
fflush(stdout);
}

size = (int) (6*sigma+1);
if(size%2==0) size = size+1;
offset = size/2+2; /* Add 2 as neighbourhood size offset */
start = middle-offset;
end = middle+offset;
printf("Size: %d Offset: %d Start: %d End: %d\n",size,offset,start,end);
printf("%d images required\n",end-start+1);
if(end < start) 
	{ 
	printf("Specify images in ascending order\n"); 
	exit(1);
	}

read_and_smooth3D(path1,argv[1],sigma,floatpic,pic,inpic,start,middle,end,header);
if(OUTPUT_SMOOTH && sigma!= 0.0)
writefiles(path2,argv[1],pic,sigma,pic_t,pic_x,pic_y,middle-2,middle+2,header);
else if(sigma != 0.0) printf("\nSmoothed images not output\n");

printf("Number of Columns: %d Number of Rows: %d\n",pic_y,pic_x);
if(pic_x > PIC_X || pic_y > PIC_Y)
	{ 
	printf("Fatal error: images are too big\n");
	exit(1);
	}

printf("\n");
if(flag==FALSE)
{
sprintf(full_name,"%s/lucas.%s%dF-%4.2f-%3.1f",path3,argv[1],middle,tau_D,sigma);
sprintf(norm_name,"%s/lucas.%s%dN-%4.2f-%3.1f",path3,argv[1],middle,tau_D,sigma);
sprintf(raw_name,"%s/lucas.%s%dR-%4.2f-%3.1f",path3,argv[1],middle,tau_D,sigma);
}
else
{
sprintf(full_name,"%s/Mlucas.%s%dF-%4.2f-%3.1f",path3,argv[1],middle,tau_D,sigma);
sprintf(norm_name,"%s/Mlucas.%s%dN-%4.2f-%3.1f",path3,argv[1],middle,tau_D,sigma);
sprintf(raw_name,"%s/Mlucas.%s%dR-%4.2f-%3.1f",path3,argv[1],middle,tau_D,sigma);
}
printf("Output full velocities go to file %s\n",full_name);
printf("Output least squares normal velocities go to file %s\n",norm_name);
printf("Output raw normal velocities go to file %s\n",raw_name);
printf("\n");
fflush(stdout);

compute_ders(Ix,Iy,It,floatpic,pic_t,pic_x,pic_y,offset,middle,argv[1]);
compute_vels(Ix,Iy,It,full_vels,norm_vels1,norm_vels2,pic_x,pic_y,2*offset,tau_D,flag,E);
if((fdf=fopen(full_name,"wb"))!=NULL){
	
	output_velocities(fdf,"Full",full_vels,pic_x,pic_y,2*offset);
}
else printf("Error in opening %s file\n\n",full_name);
if((fdn=fopen(norm_name,"wb"))!=NULL){
	
	output_velocities(fdn,"Normal",norm_vels1,pic_x,pic_y,2*offset);
}
else printf("Error in opening %s file\n\n",norm_name);
if((fdr=fopen(raw_name,"wb"))!=NULL){
	
	output_velocities(fdr,"Normal",norm_vels2,pic_x,pic_y,2*offset);
}
else printf("Error in opening %s file\n\n",raw_name);
printf("\n");
fflush(stdout);

if(strcmp(correct_filename,"unknown")!=0)
{
calc_statistics(correct_vels,norm_vels1,norm_vels2,int_size_x,int_size_y,full_vels,E,
	pic_x,pic_y,2*offset,&ave_error,&st_dev,&density,&residual,&min_angle,&max_angle,
	&norm_ave_error1,&norm_st_dev1,&norm_density1,&norm_min_angle1,&norm_max_angle1,
	&norm_ave_error2,&norm_st_dev2,&norm_density2,&norm_min_angle2,&norm_max_angle2);
printf("\n\nComputed Statistics\n");
printf("Lambda_2: %f  Error: %f St Dev: %f\n",tau_D,ave_error,st_dev);
printf("Lambda_2: %f  Density: %f\n",tau_D,density);
printf("Lambda_2: %f  Residual: %f\n",tau_D,residual);
printf("Minimum angle error: %f Maximum angle error: %f\n",min_angle,max_angle);
printf("\nLeast Squares Normal Velocity Results\n");
printf("Normal Error: %f Normal Standard Deviation: %f\n",
	norm_ave_error1,norm_st_dev1);
printf("Normal Density: %f\nMinimum Normal Angle: %f Maximum Normal Angle: %f\n",
	norm_density1,norm_min_angle1,norm_max_angle1);
printf("\nRaw Normal Velocity Results\n");
printf("Normal Error: %f Normal Standard Deviation: %f\n",
	norm_ave_error2,norm_st_dev2);
printf("Normal Density: %f\nMinimum Normal Angle: %f Maximum Normal Angle: %f\n",
	norm_density2,norm_min_angle2,norm_max_angle2);
}
fflush(stdout);
}


/*********************************************************************/
/* Compute spatio-temporal derivatives 				     */
/*********************************************************************/
void compute_ders(
float Ix[PIC_X][PIC_Y],
float Iy[PIC_X][PIC_Y],
float It[PIC_X][PIC_Y],
float floatpic[PIC_T][PIC_X][PIC_Y],
int pic_t, int pic_x, int pic_y, int n, int middle,
char *filename)
{
int i,j;
float kernel[5];
FILE *fp;
char realname[100];

sprintf(realname, "%sIxyt%d", filename, middle);
fp = fopen(realname, "wb");
fwrite(&pic_y, sizeof(int), 1, fp);
fwrite(&pic_x, sizeof(int), 1, fp);
fwrite(&n, sizeof(int), 1, fp);
fwrite(&n, sizeof(int), 1, fp);


for(i=0;i<PIC_X;i++)
for(j=0;j<PIC_Y;j++)
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
	fwrite(&Ix[i][j], sizeof(float), 1, fp);
	fwrite(&Iy[i][j], sizeof(float), 1, fp);
	fwrite(&It[i][j], sizeof(float), 1, fp);
	/* printf("i:%d j:%d Ix:%f Iy:%f It:%f\n",i,j,Ix[i][j],Iy[i][j],It[i][j]); */
	}
printf("Spatio-Temporal Intensity Derivatives Computed\n");
fflush(stdout);
fclose(fp);
}

/*********************************************************************/
/* Compute a 4 point central difference kernel			     */
/*********************************************************************/
void calc_diff_kernel(float diff_kernel[])
{
diff_kernel[0] = -1.0f/12.0f;
diff_kernel[1] = 8.0f/12.0f;
diff_kernel[2] = 0.0f;
diff_kernel[3] = -8.0f/12.0f;
diff_kernel[4] = 1.0f/12.0f;
}


/************************************************************
   Apply 1D real kernels in the x direction
************************************************************/
float diff_x(
float floatpic[PIC_T][PIC_X][PIC_Y],
float kernel[5],
int x, int y, int n)
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
float diff_y(
float floatpic[PIC_T][PIC_X][PIC_Y],
float kernel[5],
int x, int y, int n)
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
float diff_t(
float floatpic[PIC_T][PIC_X][PIC_Y],
float kernel[5],
int x, int y, int n)
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

/************************************************************
   Compute velocities
************************************************************/
void compute_vels(
float Ix[PIC_X][PIC_Y],
float Iy[PIC_X][PIC_Y],
float It[PIC_X][PIC_Y],
float full_vels[PIC_X][PIC_Y][2],
float norm_vels1[PIC_X][PIC_Y][2],
float norm_vels2[PIC_X][PIC_Y][2],
int pic_x, int pic_y, int n,
float tau_D,
int flag, float E[PIC_X][PIC_Y])
{
float mag,M[2][2],MI[2][2],B[2],denominator;
float eigenvalues[2],eigenvectors[2][2],length1,length2;
float angle,temp1,v1,v2;
float sigma1,sigma2,sigmap,temp,diff1[N],diff2[N];
int i,j,k,l,full_count,norm_count1,norm_count2,no_count,nrot,eigen_count,no_swaps;
float eigenvalues2[2],eigenvectors2[2][2],weight[5][5],sum,coeff[5],v[2];
int mag_zero;
int kk;


printf("Eigenvalue Threshold: %f\n",tau_D);
printf("Threshold on Raw Normal Velocities: %f\n",RAW_MAG);
MAX_RAW_MAG = (float) -HUGE_VAL;
mag_zero = 0;
fflush(stdout);
/* Parameter values as specified in Simoncelli, Adelson and Heeger, page 313  */
sigma1 = 0.08f;
sigma2 = 1.0f;
sigmap = 2.0f;

/* Compute weights */
sum = 0.0;
coeff[0] = coeff[4] = 0.0625;
coeff[1] = coeff[3] = 0.25;
coeff[2] = 0.375;
printf("Coefficients\n");
for(i=0;i<5;i++)
{
for(j=0;j<5;j++)
	{
	weight[i][j] = coeff[i]*coeff[j];
	printf("%12.8f \n",weight[i][j]); 
	sum += weight[i][j];
	}
printf("\n"); 
}
printf("Sum=%f\n",sum);
full_count = 0;
norm_count1 = norm_count2 = 0;
no_count = 0;
eigen_count = 0;
no_swaps = 0;
for(i=0;i<PIC_X;i++)
for(j=0;j<PIC_Y;j++)
	{
	full_vels[i][j][0] = full_vels[i][j][1] = NO_VALUE;
	norm_vels1[i][j][0] = norm_vels1[i][j][1] = NO_VALUE;
	norm_vels2[i][j][0] = norm_vels2[i][j][1] = NO_VALUE;
	if(DEBUG) for(kk=0;kk<8;kk++) data[kk][i][j] = 0.0;
	Imag[i][j] = E[i][j] = 0.0;
	}

for(i=n;i<pic_x-n;i++)
for(j=n;j<pic_y-n;j++)
	{
	M[0][0] = M[1][1] = M[0][1] = M[1][0] = 0.0;
	B[0] = B[1] = 0.0;
	mag = 1.0;
	/* Compute on 5*5 neighbourhood */
	for(k=(-2);k<=2;k++)
	for(l=(-2);l<=2;l++)
		{
		if(flag==TRUE) mag = sigma1*(Ix[i+k][j+l]*Ix[i+k][j+l]+
		                             Iy[i+k][j+l]*Iy[i+k][j+l])+sigma2;
		M[0][0] = M[0][0] + weight[k+2][l+2]*(Ix[i+k][j+l]*Ix[i+k][j+l])/mag;
		M[1][1] = M[1][1] + weight[k+2][l+2]*(Iy[i+k][j+l]*Iy[i+k][j+l])/mag;
		M[0][1] = M[0][1] + weight[k+2][l+2]*(Ix[i+k][j+l]*Iy[i+k][j+l])/mag;
		B[0] = B[0] + weight[k+2][l+2]*(Ix[i+k][j+l]*It[i+k][j+l])/mag;
		B[1] = B[1] + weight[k+2][l+2]*(Iy[i+k][j+l]*It[i+k][j+l])/mag;
		}
	if(DEBUG)
	{
	data[0][i][j] = M[0][0];
	data[1][i][j] = M[0][1];
	data[2][i][j] = M[1][1];
	data[3][i][j] = B[0];
	data[4][i][j] = B[1];
	data[5][i][j] = Ix[i][j];
	data[6][i][j] = Iy[i][j];
	data[7][i][j] = It[i][j];
	}
	M[1][0] = M[0][1]; /* The M array is symmetric */
	if(flag==TRUE)
		{
		M[0][0] = M[0][0] + 1.0f/sigmap;
		M[1][1] = M[1][1] + 1.0f/sigmap;
		}
	/* Invert 2*2 matrix */
	denominator = M[0][0]*M[1][1]-M[1][0]*M[0][1]; /* The determinant of M */
	MI[0][0] = M[1][1]/denominator;
	MI[0][1] = -M[0][1]/denominator;
	MI[1][0] = -M[1][0]/denominator;
	MI[1][1] = M[0][0]/denominator;

	jacobi(M,2,eigenvalues,eigenvectors,&nrot);
	if(check_eigen_calc(M,eigenvalues,eigenvectors,nrot,diff1,diff2,
		&length1,&length2,&angle)==FALSE)
		{
		if(FALSE)
		{
		printf("\n********************************************\n");
		printf("Fatal error: eigenvalue/eigenvector error\n");
		printf("i=%d j=%d\n",i,j);
		printf("eigenvalues: %f %f\n",eigenvalues[0],eigenvalues[1]);
		printf("eigenvector1: %f %f\n",eigenvectors[0][0],eigenvectors[1][0]);
		printf("eigenvector2: %f %f\n",eigenvectors[0][1],eigenvectors[1][1]);
		printf("\n      M:                        MI\n");
		printf("%12.6f %12.6f %12.6f %12.6f\n",M[0][0],M[0][1],MI[0][0],MI[0][1]);
		printf("%12.6f %12.6f %12.6f %12.6f\n",M[1][0],M[1][1],MI[1][0],MI[1][1]);
		printf("B: %f %f\n",B[0],B[1]);
		printf("Determinant of M: %f\n",denominator);
		printf("nrot: %d\n",nrot);
		printf("Angle between two eigenvectors: %f degrees\n",angle);
		printf("Difference length for eigenvector1\n");
		printf("Difference: %f %f Length: %f\n",diff1[0],diff1[1],length1);
		printf("Difference length for eigenvector2\n");
		printf("Difference: %f %f Length: %f\n",diff2[0],diff2[1],length2);
		/* Check eigenvalues/eigenvectors for 2*2 matrix using
		   closed form method as in Anandan's thesis */
		eigenvalues2[0] = 0.5*((M[0][0]+M[1][1]) - 
		  sqrt((M[0][0]-M[1][1])*(M[0][0]-M[1][1])+4.0*M[0][1]*M[1][0]));
		eigenvalues2[1] = 0.5*((M[0][0]+M[1][1]) +
		  sqrt((M[0][0]-M[1][1])*(M[0][0]-M[1][1])+4.0*M[0][1]*M[1][0]));
		angle = (float) atan2(eigenvalues2[0]-M[0][0],M[0][1]);
		eigenvectors2[0][0] = (float) -cos(angle);
		eigenvectors2[1][0] = (float) -sin(angle);
		eigenvectors2[0][1] = (float) -sin(angle);
		eigenvectors2[1][1] = (float) cos(angle);
		printf("\nUsing Anandan's calculation:\n");
		printf("Angle of rotation: %f degrees\n",angle*180/HALF_C);
		printf("eigenvalues: %f %f\n",eigenvalues2[0],eigenvalues2[1]);
		printf("eigenvector1: %f %f\n",eigenvectors2[0][0],eigenvectors2[1][0]);
		printf("eigenvector2: %f %f\n",eigenvectors2[0][1],eigenvectors2[1][1]);
		check_eigen_calc(M,eigenvalues2,eigenvectors2,nrot,diff1,diff2,
				 &length1,&length2,&angle);
		printf("Angle between two eigenvectors: %f degrees\n",angle);
		printf("Difference length for eigenvector1\n");
		printf("Difference: %f %f Length: %f\n",diff1[0],diff1[1],length1);
		printf("Difference length for eigenvector2\n");
		printf("Difference: %f %f Length: %f\n",diff2[0],diff2[1],length2);
		printf("********************************************\n\n");
		fflush(stdout);
		if(FALSE) exit(1);
		}
		eigen_count++;
		}
	else 
	{
	/* Sort the eigenvalues and the corresponding eigenvectors */
	/* Most likely, already ordered				   */
	if(eigenvalues[0] < eigenvalues[1]) /* Largest eigenvalue first */
		{
		/* swap eigenvalues */
		temp = eigenvalues[0];
		eigenvalues[0] = eigenvalues[1];
		eigenvalues[1] = temp;
		/* swap eigenvector components */
		temp = eigenvectors[0][0];
		eigenvectors[0][0] = eigenvectors[0][1];
		eigenvectors[0][1] = temp;
		temp = eigenvectors[1][0];
		eigenvectors[1][0] = eigenvectors[1][1];
		eigenvectors[1][1] = temp;
		no_swaps++;
		}
	
	/* Full velocity if spread of M is small */
	if(eigenvalues[0] >= tau_D && eigenvalues[1] >= tau_D)
	{
	if(denominator > 0.0)
		{
		full_vels[i][j][1] = -(v1= -(MI[0][0]*B[0]+MI[0][1]*B[1]));
		full_vels[i][j][0] =  (v2= -(MI[1][0]*B[0]+MI[1][1]*B[1]));
		/* Compute the residual */
		for(k=(-2);k<=2;k++)
		for(l=(-2);l<=2;l++)
			{
			temp1 = weight[k+2][l+2]*
			      (Ix[i+k][j+l]*v1+
			       Iy[i+k][j+l]*v2+It[i+k][j+l]);
			/* temp2 = weight[k+2][l+2]*It[i+k][j+l]; */
			E[i][j] += (temp1*temp1);
			}
		full_count++;
		}
	else { full_vels[i][j][0] = full_vels[i][j][1] = NO_VALUE; }
	}
	/* Normal velocity if spread of M is small in one direction only */
	else if(eigenvalues[0] > tau_D && fabs(denominator) > 0.00000001)
	/* Normal velocity if spread of MI is small in one direction only */
		{
		/* Project v onto that direction */
		v[0] = -(MI[0][0]*B[0]+MI[0][1]*B[1]);
		v[1] = -(MI[1][0]*B[0]+MI[1][1]*B[1]);
		norm_vels1[i][j][1] = -(v[0]*eigenvectors[0][0]+
				      v[1]*eigenvectors[1][0])*eigenvectors[0][0];
		norm_vels1[i][j][0] = (v[0]*eigenvectors[0][0]+
				      v[1]*eigenvectors[1][0])*eigenvectors[1][0];
		norm_count1++;
		/* if(i >= 50 && i <= 60 && j >= 50 && j <= 60)
		printf("Normal velocity at i:%d j:%d: %f %f\n",i,j,norm_vels1[i][j][0],norm_vels1[i][j][1]); */
		}
	else 
		{
		/* printf("No velocity\n"); */
		no_count++;
		}
	}

	/* Compute type 2 normal velocity */
	mag = (Ix[i][j]*Ix[i][j] + Iy[i][j]*Iy[i][j]);
	Imag[i][j] = (float) sqrt(mag);
	if(Imag[i][j] > MAX_RAW_MAG) MAX_RAW_MAG = Imag[i][j];
	if(Imag[i][j] > RAW_MAG)
	{
	norm_vels2[i][j][1] =  It[i][j]*Ix[i][j]/mag;
	norm_vels2[i][j][0] = -It[i][j]*Iy[i][j]/mag;
	norm_count2++;
	}
	else mag_zero++;
	}

printf("%d full velocities computed\n",full_count);
printf("%d least squares normal velocities computed\n",norm_count1);
printf("%d raw normal velocities computed\n",norm_count2);
printf("%d locations where velocity information thresholded\n",no_count);
printf("%d locations where eigenvalue/eigenvector calculation failed\n",eigen_count);
printf("%d locations with spatial gradient is zero\n",mag_zero);
fflush(stdout);
if(DEBUG) for(kk=0;kk<8;kk++) fwrite(&data[kk][0][0], sizeof(float), 150*150, fd_temp[kk]);
}

/************************************************************
   Output full velocities using old Burkitt format
************************************************************/
void output_velocities(FILE *fdf, char s[],
float full_velocities[PIC_X][PIC_Y][2],
int pic_x, int pic_y,
int n)
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
y = (float) pic_x;
x = (float) pic_y;
fwrite(&x,4, 1, fdf);
fwrite(&y,4, 1, fdf);

/* size of result data */
y = (float) (pic_x-2*n);
x = (float) (pic_y-2*n);
fwrite(&x,4, 1, fdf);
fwrite(&y,4, 1, fdf);

/* offset to start of data */
y = (float) n;
x = (float) n;
fwrite(&x,4, 1, fdf);
fwrite(&y,4, 1, fdf);
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
		}
	else
		{
		no_novals++;
		}
	}
for(i=n;i<pic_x-n;i++)
	{
	bytes += fwrite(&full_velocities[i][n][0], 4, (pic_y-2*n)*2, fdf);
	}
fclose(fdf);
printf("\n%s velocities output from output_velocities: %d velocities\n",s,bytes);
printf("Number of positions with velocity: %d\n",no_vals);
printf("Number of positions without velocity: %d\n",no_novals);
printf("Percentage of %s velocities: %f\n",s,
	no_vals/(1.0*(no_vals+no_novals))*100.0);
fflush(stdout);
}


/************************************************************************/
/* Compute all eigenvalues and eigenvectors of a real symmetric matrix  */
/* a[N][N]. On output elements of a above the disgonal are destroyed.   */
/* d[N] returns the eigenvalues of a. v[N][N] is a matrix whose columns */
/* contain, on output, the normalized eigenvectors of a. nrot returns   */
/* the number of Jacobi rotations that were required.			*/
/************************************************************************/
void jacobi(float aa[N][N], int n, float d[N], float v[N][N],
int *nrot)
{
int j,iq,ip,i;
float thresh,theta,tau,t,sm,s,h,g,c;
float b[N],z[N],a[N][N];

if(n!=N) { fprintf(stderr,"\nFatal error: n not N in jacobi\n",N); exit(1); }
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
			sm += (float) fabs(a[ip][iq]);
		}

	/* Normal return, which relies on quadratic convergence to
	   machine underflow */
	if(sm == 0.0) return;

	if(i<3) thresh=(float) (0.2*sm/(n*n)); /* on the first three sweeps */
	else thresh = 0.0; /* the rest of the sweeps */

	for(ip=0;ip<(n-1);ip++)
		{
		for(iq=ip+1;iq<n;iq++)
			{
			g = (float)(100.0*fabs(a[ip][iq]));
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
				  theta = 0.5f*h/(a[ip][iq]);
				  t = (float)(1.0/(fabs(theta)+sqrt(1.0+theta*theta)));
				  if(theta < 0.0) t = -t;
				  }
				c = (float)(1.0/sqrt(1.0+t*t));
				s = t*c;
				tau = s/(1.0f+c);
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
void rotate(float a[N][N], int i, int j, int k, int l,
float *h, float *g, float s, float tau)
{
(*g) = a[i][j];
(*h) = a[k][l];
a[i][j] = (*g)-s*((*h)+(*g)*tau);
a[k][l] = (*h)+s*((*g)-(*h)*tau);
}

/*********************************************************************/
/* Check eigenvector and eigenvalue computation for 2*2 matrix	     */
/*********************************************************************/
int check_eigen_calc(float mm[N][N],
float d[N], float v[N][N], int n, 
float diff1[N], float diff2[N], float *length1, float *length2,
float *angle)
{
int status;

status = TRUE;
/* Compute angle between two eigenvectors - should be orthogonal */
(*angle)=acos((v[0][0]*v[0][1]+v[1][0]*v[1][1])/
	 (sqrt(v[0][0]*v[0][0]+v[1][0]*v[1][0])*
	  sqrt(v[0][1]*v[0][1]+v[1][1]*v[1][1])))*180.0/3.1415926535;
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

/***************************************************************/
/*  Compute error statistics				       */
/***************************************************************/
void calc_statistics(float correct_vels[PIC_X][PIC_Y][2],
				float norm_vels1[PIC_X][PIC_Y][2],
				float norm_vels2[PIC_X][PIC_Y][2],
				int int_size_x, int int_size_y,
				float full_vels[PIC_X][PIC_Y][2],
				float E[PIC_X][PIC_Y], int pic_x,
				int pic_y, int n, float *ave_error,
				float *st_dev,float *density, float *residual,
				float *min_angle, float *max_angle,
				float *norm_ave_error1,float *norm_st_dev1,
				float *norm_density1, float *norm_min_angle1,
				float *norm_max_angle1, float *norm_ave_error2,
				float *norm_st_dev2, float *norm_density2,
				float *norm_min_angle2, float *norm_max_angle2)
{
FILE *fp;
int full_count,norm_count1,no_full_count,no_norm_count1,i,j,total_count;
int norm_count2,no_norm_count2;
float normal_sumX2_1,normal_sumX2_2,sumX2,temp,uva[2],uve[2],bin_err2;
float histogram_error[100],histogram_error2[100],bin_sum,bin_ave,bin_st_dev;
float bin_density;
int histogram_count[100],int_mag,bin_ct;

full_count = norm_count1 = norm_count2 = no_full_count = no_norm_count1 = no_norm_count2 = total_count = 0;
normal_sumX2_1 = normal_sumX2_2 = sumX2 = 0.0;
(*norm_min_angle1) = (*min_angle) = (float) HUGE_VAL;
(*norm_max_angle1) = (*max_angle) = 0.0f;
(*norm_min_angle2) = (float) HUGE_VAL;
(*norm_max_angle2) = 0.0f;
(*ave_error) = (*st_dev) = (*density) = (*residual) = 0.0f;
(*norm_ave_error1) = (*norm_st_dev1) = (*norm_density1) = 0.0f;
(*norm_ave_error2) = (*norm_st_dev2) = (*norm_density2) = 0.0f;

for(i=0;i<100;i++) { histogram_error[i]=histogram_error2[i]=0.0; histogram_count[i]=0; }


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
	  if(E[i][j] == NO_VALUE) 
		{ printf("Fatal error: E has no value\n"); exit(1); }
	  (*residual) += E[i][j];
	  sumX2 += temp*temp;
	  if(temp < (*min_angle)) (*min_angle) = temp;
	  if(temp > (*max_angle)) (*max_angle) = temp;
	  }
	else no_full_count++;

	if(norm_vels1[i][j][0] != NO_VALUE && norm_vels1[i][j][1] != NO_VALUE)
	  {
	  norm_count1++;
	  uve[0] = norm_vels1[i][j][0]; uve[1] = norm_vels1[i][j][1];
	  uva[0] = correct_vels[i][j][0]; uva[1] = correct_vels[i][j][1];
	  temp = PsiEN(uve,uva);
	  (*norm_ave_error1) += temp;
	  normal_sumX2_1 += temp*temp;
	  if(temp < (*norm_min_angle1)) (*norm_min_angle1) = temp;
	  if(temp > (*norm_max_angle1)) (*norm_max_angle1) = temp;
	  }
	else no_norm_count1++;

	if(norm_vels2[i][j][0] != NO_VALUE && norm_vels2[i][j][1] != NO_VALUE)
	  {
	  norm_count2++;
	  uve[0] = norm_vels2[i][j][0]; uve[1] = norm_vels2[i][j][1];
	  uva[0] = correct_vels[i][j][0]; uva[1] = correct_vels[i][j][1];
	  temp = PsiEN(uve,uva);
	  (*norm_ave_error2) += temp;
	  normal_sumX2_2 += temp*temp;
	  if(temp < (*norm_min_angle2)) (*norm_min_angle2) = temp;
	  if(temp > (*norm_max_angle2)) (*norm_max_angle2) = temp;
	  int_mag = (int) (Imag[i][j]);
	  if(int_mag >= MAX_I) int_mag = MAX_I;
	  histogram_error[int_mag] += temp;
	  histogram_error2[int_mag] += temp*temp;
	  histogram_count[int_mag]++;
	  }
	else no_norm_count2++;

	total_count++;
	}
}
(*density) = (full_count*100.0f)/(total_count);
(*norm_density1) = (norm_count1*100.0f)/(total_count);
(*norm_density2) = (norm_count2*100.0f)/(total_count);

if(full_count != 0) (*ave_error) = (*ave_error)/full_count;
else (*ave_error) = 0.0;
if(norm_count1 != 0) (*norm_ave_error1) = (*norm_ave_error1)/norm_count1;
else (*norm_ave_error1) = 0.0;
if(norm_count2 != 0) (*norm_ave_error2) = (*norm_ave_error2)/norm_count2;
else (*norm_ave_error2) = 0.0;

if(full_count > 1) 
{
temp =(float) (fabs((sumX2 - full_count*(*ave_error)*(*ave_error))/(full_count-1)));
(*st_dev) = (float) sqrt(temp);
}
else (*st_dev) = 0.0f;
if(norm_count1 > 1) 
(*norm_st_dev1) = (float) (sqrt((normal_sumX2_1 - norm_count1*(*norm_ave_error1)*(*norm_ave_error1))/(norm_count1-1)));
else (*norm_st_dev1) = 0.0;
if(norm_count2 > 1) 
(*norm_st_dev2) = (float) (sqrt((normal_sumX2_2 - norm_count2*(*norm_ave_error2)*(*norm_ave_error2))/(norm_count2-1)));
else (*norm_st_dev2) = 0.0;

if(full_count != 0) (*residual) = (*residual)/full_count;
if((*ave_error) == 0.0) { (*min_angle) = (*max_angle) = 0.0; }
if((*norm_ave_error1) == 0.0) { (*norm_min_angle1) = (*norm_max_angle1) = 0.0; }
if((*norm_ave_error2) == 0.0) { (*norm_min_angle2) = (*norm_max_angle2) = 0.0; }

printf("\nIn calc_statistics\n");
printf("%d full velocities\n",full_count);
printf("%d positons without full velocity\n",no_full_count);
printf("%d positions in total\n",total_count);
printf("%d least squares normal velocities\n",norm_count1);
printf("%d positons without least squares normal velocity\n",no_norm_count1);
printf("%d positions with full or least squares normal velocity\n",
	total_count-full_count-norm_count1);
printf("%d raw normal velocities\n",norm_count2);
printf("%d positons without raw normal velocity\n",no_norm_count2);

fp =fopen("lucas.plot.data","w");
fprintf(fp,"%d\n",MAX_I);
if(RAW_STATISTICS)
{
printf("\n                             Raw Normal Velocity Analysis\n");
printf("\n                                  Histogram Error\n\n");
for(i=0;i<MAX_I;i++) 
	{
     	if(histogram_count[i] != 0.0) bin_ave = histogram_error[i]/(histogram_count[i]*1.0f);
     	else bin_ave = 0.0;
	if(histogram_count[i] > 1) temp = ((histogram_error2[i]-histogram_count[i]*bin_ave*bin_ave)/(histogram_count[i]-1.0f));
	else temp = 0.0;
	bin_density = histogram_count[i]*100.0f/total_count;
	bin_st_dev = (float) sqrt(fabs(temp));
     	printf("Bin:%5.2f Average Error:%8.5f St. Dev:%8.5f Count:%5d Density:%6.3f\n",i+0.5,bin_ave,bin_st_dev,histogram_count[i],bin_density);
     	fprintf(fp,"%f %f %f %f\n",i+0.5,bin_ave,bin_st_dev,bin_density);
	}
fprintf(fp,"\n\n\n");
fprintf(fp,"%d\n",MAX_I);
printf("\n                                 Cumulative Error\n\n");
bin_sum = 0.0; bin_ct=0; bin_err2 = 0.0;
for(i=0;i<MAX_I;i++)
     {
     bin_sum += histogram_error[i];
     bin_err2 += histogram_error2[i];
     bin_ct += histogram_count[i];
     }
for(i=0;i<MAX_I;i++)
     {
     if(bin_ct != 0) bin_ave = bin_sum/bin_ct; else bin_ave = 0.0;
     bin_density = bin_ct*100.0f/total_count;
     if(bin_ct > 1)  temp = ((bin_err2-bin_ct*bin_ave*bin_ave)/(bin_ct-1));
     else temp = 0.0;
     bin_st_dev = (float) sqrt(fabs(temp));
     printf("Bin:>%5.2f Average Error:%8.5f St. Dev:%8.5f Count:%5d Density:%6.3f\n",i*1.0,bin_ave,bin_st_dev,bin_ct,bin_density);
     fprintf(fp,"%f %f %f %f\n",i*1.0,bin_ave,bin_st_dev,bin_density);
     bin_sum -= histogram_error[i];
     bin_err2 -= histogram_error2[i];
     bin_ct -= histogram_count[i];
     }
fprintf(fp,"\n");
fclose(fp);
}
fflush(stdout);
}

/************************************************************
 Full Image Velocity Angle Error
************************************************************/
float PsiER(float ve[], float va[])
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

nva = norm(VA,3);
nve = norm(VE,3);
v =  (VE[0]*VA[0]+VE[1]*VA[1]+1.0f)/(nva*nve);

/**  sometimes roundoff error causes problems **/
if(v>1.0 && v < 1.0001) v = 1.0;

r = (float) (acos(v)*180.0/HALF_C);

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
float norm(float v[],int n)
{
int i;
float sum = 0.0;

for (i=0;i<n; i++)
        sum += (v[i]*v[i]);
sum = (float) sqrt(sum);
return sum;
}


/************************************************************
 Normal Image Velocity Angle Error
************************************************************/
float PsiEN(float ve[],float va[])
{
float nva,nve;
float v1,v2;
float n[2];

nva = norm(va,2), nve = norm(ve,2);
if(nve > 0.00000001)
{
n[0] = ve[0]/nve;
n[1] = ve[1]/nve;
v1 = (va[0]*n[0] + va[1]*n[1]-nve) ;
v2 = (float) (v1/(sqrt((1.0+nva*nva))*sqrt((1.0+nve*nve))));
v1 =  (float) (asin(v2)*180.0/HALF_C);

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
	
return (float) fabs(v1);
}

/*********************************************************************/
/* Smooth the image sequence using a 3D separable Gaussian filter.   */
/*********************************************************************/
void read_and_smooth3D(char path[], char stem[],
float sigma, float floatpic[][PIC_X][PIC_Y],
 unsigned char pic[][PIC_X][PIC_Y],
unsigned char inpic[][PIC_X][PIC_Y],
int start, int middle, int end,
unsigned char header[])
{
int i,j,k,time,frame;

for(k=0;k<FIVE;k++)
for(i=0;i<PIC_X;i++)
for(j=0;j<PIC_Y;j++)
	floatpic[k][i][j] = 0.0;

pic_t = end-start+1;


readfiles(path,stem,inpic,sigma,pic_t,&pic_x,&pic_y,start,end,header);
if(sigma==0.0)
	{
	for(k=0;k<FIVE;k++)
	for(i=0;i<pic_x;i++)
	for(j=0;j<pic_y;j++)
		floatpic[k][i][j] = inpic[k][i][j];
	return;
	}

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
printf("Input files smoothed\n");
}

/************************************************************************/
/* Perform 3D Gaussian smoothing by separable convolution 		*/
/************************************************************************/
void convolve_Gaussian(
unsigned char inpic[PIC_T][PIC_X][PIC_Y],
float floatpic[FIVE][PIC_X][PIC_Y],
unsigned char pic[FIVE][PIC_X][PIC_Y],
float sigma, int pic_t, int pic_x, int pic_y,
int start, int frame, int time)
{
float mask[100],term,sum;
int size,i,j,offset,a,b;


if(OUTPUT_SMOOTH)
{
printf("\nStart of 3D convolution\n");
printf("Time: %d Frame: %d\n",time,frame);
}

fflush(stdout);
size = (int) (6*sigma+1);
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

if(OUTPUT_SMOOTH)
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
for(j=0;j<pic_y;j++)
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
	if(term > 255.0) term = 255.0;
	if(term < 0.0) term = 0.0;
	floatpic[time][i][j] = term;
	pic[time][i][j] = (int) (term+0.5);
	}
if(OUTPUT_SMOOTH) printf("Convolution in y direction completed\n");
}
else /* No smoothing */
for(i=0;i<pic_x;i++)
for(j=0;j<pic_y;j++)
	pic[time][i][j] = inpic[frame-start][i][j];

if(OUTPUT_SMOOTH) printf("End of Convolution\n");
fflush(stdout);
}


/*********************************************************************/
/* Read input images and perform Gaussian smoothing.		     */
/*********************************************************************/
void readfiles(char path[], char s[],
unsigned char pic[][PIC_X][PIC_Y],
float sigma, int pic_t, int *pic_x, int *pic_y, int start, int end,
unsigned char header[HEAD])

{
char fname[200];
int i,j,time,no_bytes;
int ONCE;
int ints[8];
FILE *fp;


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
	sprintf(fname,"%s/%s%d.ras",path,s,i);
 	if((fp=fopen(fname,"rb")) != NULL)
		{
		if(!BINARY)
		{
		if(ONCE)
		{
		no_bytes += (4*fread(ints, sizeof(int), 8, fp));
		(*pic_y) = ints[1];
		(*pic_x) = ints[2];
		ONCE = FALSE;
		}
		else no_bytes += fread(header, sizeof(unsigned char), HEAD, fp);
		}

		for(j=0;j<(*pic_x);j++)
			no_bytes += fread(&pic[time][j][0],sizeof(unsigned char), (*pic_y), fp);
		printf("File %s read (%d bytes)\n",fname,no_bytes);
		no_bytes = 0;
		fclose(fp);
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
void writefiles(char path[], char s[],
unsigned char result[FIVE][PIC_X][PIC_Y],
float sigma,
int pic_t, int pic_x, int pic_y, int start, int end,
unsigned char header[HEAD])
{
char fname[100];
int i,j,time,no_bytes;
FILE *fd;

if(sigma==0.0) return;
printf("\nWriting smoothed files...\n");

time = -1;
for (i=start;i<=end;i++) 
	{
	no_bytes = 0;
	time++;
	sprintf(fname,"%s/smoothed.%s%d-%3.1f",path,s,i,sigma);
 	if((fd=fopen(fname,"wb"))!=NULL)
		{
		no_bytes += fwrite(&header[0], sizeof(unsigned char), HEAD, fd); /* Write 32 byte raster header */
		for(j=0;j<pic_x;j++)
			no_bytes += fwrite(&result[time][j][0],sizeof(unsigned char), pic_y, fd);
		printf("File %s written (%d bytes)\n",fname,no_bytes);
		}
	else 
		{
		printf("File %s cannot be written\n",fname);
		exit(1);
		}
	fflush(stdout);
	fclose(fd);
	}
} /* End of writefiles */
