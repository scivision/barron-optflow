#include <math.h>
#include <stdio.h>
#include <fcntl.h>

/* 
           NAME : const.h 
   PARAMETER(S) : none
 
        PURPOSE : Definition of the constants for the 
                  implementation of Uras' motion 
                  detection approach.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : September 13 1990

*/

#define TRUE        1
#define FALSE       0

#define NO_REG      0
#define TR1         1
#define TR2         2
#define PI          3.14159

#define H           32
#define SMSK        5
#define REGION_SIZE 8
#define X           316
#define Y           316
#define Z           20

#define DEF_S1      3.0 
#define DEF_S2      1.5
#define DEF_S3      3.0

#define N_HISTO     4
#define N_BINS      100

#define MAX_COND    100000.0

#define NRADIUS     4
#define SKIP        1

#define MAXFLOW      20.0
#define NO_ERROR     0
#define SA_OVERFLOW  1
#define NO_FLOW      2

/* 
           NAME : type.h 
   PARAMETER(S) : none
 
        PURPOSE : Type definitions for images
                  and related data structures.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : September 13 1990

*/

typedef struct t_raster {
                          int m, width, height, depth, length,
                              type, maptype, maplength ;
                        } raster_t ;

typedef struct t_beaudet {
                           float g[SMSK][SMSK] ;
                           float f ;
                           int m ;
                         } beaudet_t ;

typedef struct t_kernel {
                          float k[X], f ;
                          int m ;
                        } kernel_t ;

typedef struct t_disp_vect {
                             float x, y ;  
                           } disp_vect_t ;

typedef struct t_param {
                         float fxx, fyy, fxy, ft, fxt, fyt, fx, fy, 
                               discr, gauss, cond ;
                         int err ;
                       } param_t ;

typedef struct t_pos { int i, j ; } pos_t ;

typedef disp_vect_t disp_field512_t[X*Y] ;

typedef param_t param512_t[X*Y] ;

typedef float image512_t[X*Y] ;

typedef struct t_qnode {
                         int             res, sizx, sizy, sizz, ofst, level ;
                         image512_t      *gauss_ptr[Z] ;
                         param512_t      *param_ptr[Z] ;
                         disp_field512_t *flow_ptr ;
                         struct t_qnode  *forth, *back ;
                       } qnode_t, *qnode_ptr_t ;

typedef struct t_histo {
                         float std, avg ;
                         int   freq ;
                       } histo_t ;

/* 
           NAME : str.h
   PARAMETER(S) : none
 
        PURPOSE : Definitions of constants and types for 
                  string manipulation.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : June 25 1990
*/

#define STRLENGTH 80

typedef char string[STRLENGTH] ;

/* 
           NAME : extvar.h
   PARAMETER(S) : none
 
        PURPOSE : declaration of external variables; all
                  masks are external.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : September 13 1990
*/

beaudet_t Ix, Iy, Ixx, Ixy ;
kernel_t ker1, ker2, ker3, C ;
int KERNEL_X, KERNEL_Y ;

/* 
           NAME : concat(s1,s2,s3) ;
   PARAMETER(S) : s1, s2 : strings to concat;
                      s3 : output string.
 
        PURPOSE : Concats s1 and s2 into s3.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : June 25 1990
*/

concat(s1,s2,s3)
string s1, s2, s3 ;

{ string t ;
  int i, j ;

  for (i = 0 ; (t[i] = s1[i]) != '\0' ; i++) ;
  for (j = i ; (t[j] = s2[j - i]) != '\0' ; j++) ;
  t[j] = '\0' ;
  for (i = 0 ; i <= j ; i++) {
    s3[i] = t[i] ;
  }
}

/* 
           NAME : condition(loc,f)
   PARAMETER(S) : loc : location on flow;
                    f : flow pointer.
 
        PURPOSE : returns condition number of estimate at loc in f

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : November 7 1990
*/

float condition(loc,f)
pos_t loc ;
qnode_ptr_t f ;

{ float cond ;

  if ((loc.i != -1) && (loc.j != -1)) {
    cond = (*f->param_ptr[f->sizz/2])[f->res*loc.i + loc.j].cond ;
  }
  else {
    cond = 0.0 ;
  }
  return(cond) ;
}

/* 
           NAME : convolve(c1,c2,ker1,ker2,n)
   PARAMETER(S) :     c1,c2 : pointers on image cube nodes;
                  ker1,ker2 : kernels  used for 3D convolution;
                          n : number of image frames in cubes.
 
        PURPOSE : Performs a 3D convolution on images of c1 in c2
                  using 2 1D kernels.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 20 1990

*/

convolve(c1,c2,ker1,ker2,n)
qnode_ptr_t c1, c2 ;
kernel_t ker1, ker2 ;
int n ;

{ image512_t *t ;
  float s ;
  int i, j, k, l, h1, h2 ;

  h2 = ker2.m/2 ;
  
  for (k = h2 ; k < n - h2 ; k++) {
    for (i = 0 ; i < c1->sizy ; i++) {
      for (j = 0 ; j < c1->sizx ; j++) {
        s = 0.0 ;
        for (l = 0 ; l < ker2.m ; l++) {
          s += (*c1->gauss_ptr[l+k-h2])[c1->res*i + j]*ker2.k[l] ;
        }
        (*c2->gauss_ptr[k-h2])[c2->res*i + j] = s/ker2.f ;
      }
    }
  }

  h1 = ker1.m/2 ;

  for (k = h2 ; k < n - h2 ; k++) {

    t = (image512_t *)malloc(sizeof(image512_t)) ;

    for (i = 0 ; i < c2->sizy ; i++) {
      for (j = h1 ; j < c2->sizx - h1 ; j++) {
        s = 0.0 ;
        for (l = 0 ; l < ker1.m ; l++) {
          s += (*c2->gauss_ptr[k-h2])[c2->res*i + j+l-h1]*ker1.k[l] ;
        }
        (*t)[c2->res*i + j] = s/ker1.f ;
      }
    }
  
    for (i = h1 ; i < c2->sizy - h1 ; i++) {
      for (j = 0 ; j < c2->sizx ; j++) {
        s = 0.0 ;
        for (l = 0 ; l < ker1.m ; l++) {
          s += (*t)[c2->res*(i+l-h1) + j]*ker1.k[l] ;
        }
        (*c2->gauss_ptr[k-h2])[c2->res*i + j] = s/ker1.f ;
      }
    }
    free((image512_t *)t) ;
  }
  c2->ofst += h1 ;
}

/* 
           NAME : qnode_ptr_t create_node(l,r,sizx,sizy,sizz,ofst) 
   PARAMETER(S) :   l          : level in the pyramid;
                    r          : resolution at level l;
                    sizx, sizy : image size;
                    sizz       : sequence length;
                    ofst       : offset from image boundaries.
 
        PURPOSE : Creation of a node with 
                  level l.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 17 1990

*/

qnode_ptr_t create_node(l,r,sizx,sizy,sizz,ofst)
int l, r, sizx, sizy, sizz, ofst ;

{ qnode_ptr_t p ; 
  
  p = (qnode_ptr_t)malloc(sizeof(qnode_t)) ;
  p->res = r ;
  p->sizx = sizx ;
  p->sizy = sizy ;
  p->sizz = sizz ;
  p->ofst = ofst ;
  p->level = l ;
  return(p) ;
}

/* 
           NAME : delete_node(p,h,q)
   PARAMETER(S) : p : pointer on the node to be deleted;
                  h : head list pointer;
                  q : tail list pointer. 

        PURPOSE : Deletion of the node pointed by p.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 17 1990

*/

delete_node(p,h,q)
qnode_ptr_t p, *h, *q ;

{ qnode_ptr_t s, t ;

  s = p->back ;
  t = p->forth ;
  if (t != (qnode_ptr_t)NULL) {
    t->back = s ;
  }
  else {
    *q = s ;
  }
  if (s != (qnode_ptr_t)NULL) {
    s->forth = t ;
  }
  else {
    *h = t ;
  }
  free(p) ;
}

/* 
           NAME : determinant(loc,f)
   PARAMETER(S) : loc : location on flow;
                    f : flow pointer.
 
        PURPOSE : returns determinant of estimate at loc in f

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : November 7 1990
*/

float determinant(loc,f)
pos_t loc ;
qnode_ptr_t f ;

{ float det ;

  if ((loc.i != -1) && (loc.j != -1)) {
    det = (*f->param_ptr[f->sizz/2])[f->res*loc.i + loc.j].gauss ;
  }
  else {
    det = 0.0 ;
  }
  return(det) ;
}

/* 
           NAME : dt(c,q,x,y)
   PARAMETER(S) :      c : image cube node;
                       q : image parameter node ;
                    x, y : image location for computation.
 
        PURPOSE : computes 1st order derivative with respect
                  to time at x,y for middle image of cube c.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 1 1990
*/

float dt(c,q,x,y)
qnode_ptr_t c, q ;
int x, y ;

{ extern kernel_t C ;
  float d ;
  int i, k, h ;

  h = C.m/2 ;
  k = c->sizz/2 ;

  d = 0.0 ;
 for (i = -h ; i <= h ; i++) {
    d += (*c->gauss_ptr[i+h])[c->res*x + y]*C.k[i+h] ;
  }
  d /= C.f ;
  return(d) ;
}

/* 
           NAME : overflow(x,y,sizx,sizy,ofst,h,l)
   PARAMETER(S) : x,y       : image coordinates ;
                  sizx,sizy : image size;
                  ofst      : offset from image boundaries;
                    h       : even half size of odd masks applied;
                    l       : number of mask applications.
 
        PURPOSE : Returns TRUE if x,y in image, 
                  FALSE otherwise.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : July 23 1990
*/

int overflow(x,y,sizx,sizy,ofst,h,l) 
int x, y, sizx, sizy, ofst, h, l ;

{
  return((x < ofst + h*l) || 
         (x >= sizx - (ofst + h*l)) || 
         (y < ofst + h*l) || 
         (y >= sizy - (ofst + h*l))) ;
}

/* 
           NAME : dux(p,x,y)
   PARAMETER(S) : p : flow node;
                x,y : flow location for computation.
 
        PURPOSE : computes 1st order partial derviative with respect
                  to x at x,y for component u of flow p.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 4 1990
*/

float dux(p,x,y)
qnode_ptr_t p ;
int x, y ;

{ extern beaudet_t Ix ;
  float d ;
  int i, j, h ;

  h = Ix.m/2 ;
  d = 0.0 ;

  if (!overflow(x,y,p->sizy,p->sizx,p->ofst,h,3)) {
    for (i = -h ; i <= h ; i++) {
      for (j = -h ; j <= h ; j++) {
        d = d + (*p->flow_ptr)[p->res*(x+i) + y+j].x*Ix.g[i+h][j+h] ;
      }
    }
  }
  return(d/Ix.f) ;
}

/* 
           NAME : duy(p,x,y)
   PARAMETER(S) : p : flow node;
                x,y : flow location for computation.
 
        PURPOSE : computes 1st order partial derviative with respect
                  to y at x,y for component u of flow p.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 4 1990
*/

float duy(p,x,y)
qnode_ptr_t p ;
int x, y ;

{ extern beaudet_t Ix ;
  float d ;
  int i, j, h ;
  
  h = Ix.m/2 ;
  d = 0.0 ;
  
  if (!overflow(x,y,p->sizy,p->sizx,p->ofst,h,3)) {
    for (i = -h ; i <= h ; i++) {
      for (j = -h ; j <= h ; j++) {
        d = d + (*p->flow_ptr)[p->res*(x+i) + y+j].x*Ix.g[j+h][i+h] ;
      }
    }
  }
  return(d/Ix.f) ;
}

/* 
           NAME : dvx(p,x,y)
   PARAMETER(S) : p : flow node;
                x,y : flow location for computation.
 
        PURPOSE : computes 1st order partial derviative with respect
                  to x at x,y for component v of flow p.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 4 1990
*/

float dvx(p,x,y)
qnode_ptr_t p ;
int x, y ;

{ extern beaudet_t Ix ;
  float d ;
  int i, j, h ;

  h = Ix.m/2 ;
  d = 0.0 ; 

  if (!overflow(x,y,p->sizy,p->sizx,p->ofst,h,3)) {
    for (i = -h ; i <= h ; i++) {
      for (j = -h ; j <= h ; j++) {
        d = d + (*p->flow_ptr)[p->res*(x+i) + y+j].y*Ix.g[i+h][j+h] ;
      }
    }
  }
  return(d/Ix.f) ;
}

/* 
           NAME : dvy(p,x,y)
   PARAMETER(S) : p : flow node;
                x,y : flow location for computation.
 
        PURPOSE : computes 1st order partial derviative with respect
                  to y at x,y for component v of flow p.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 4 1990
*/

float dvy(p,x,y)
qnode_ptr_t p ;
int x, y ;

{ extern beaudet_t Ix ;
  float d ;
  int i, j, h ;

  h = Ix.m/2 ;
  d = 0.0 ;

  if (!overflow(x,y,p->sizy,p->sizx,p->ofst,h,3)) {
    for (i = -h ; i <= h ; i++) {
      for (j = -h ; j <= h ; j++) {
        d = d + (*p->flow_ptr)[p->res*(x+i) + y+j].y*Ix.g[j+h][i+h] ;
      }
    }
  }
  return(d/Ix.f) ;
}

/* 
           NAME : dx(p,q,x,y,z)
   PARAMETER(S) : p : image cube node;
                  q : parameters of flow (image derivatives);
              x,y,z : image location for computation.
 
        PURPOSE : computes 1st order partial derviative with respect
                  to x at x,y for image p->gauss_ptr[z].

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 1 1990
*/

float dx(p,q,x,y,z)
qnode_ptr_t p, q ;
int x, y, z ;

{ extern kernel_t C ;
  float d ;
  int i, h ;

  h = C.m/2 ;
  d = 0.0 ;

  if (overflow(x,y,p->sizy,p->sizx,p->ofst,h,1)) {
    (*q->param_ptr[z])[q->res*x + y].err = SA_OVERFLOW ;
  }
  else {
    (*q->param_ptr[z])[p->res*x + y].err = NO_ERROR ;
    for (i = -h ; i <= h ; i++) {
      d += (*p->gauss_ptr[z])[p->res*(x+i) + y]*C.k[i+h] ;
    }
    d /= C.f ;
  }
  return(d) ;
}

/* 
           NAME : dxt(p,x,y)
   PARAMETER(S) : p : image parameter node;
                x,y : image location for computation. 

        PURPOSE : computes partial derivative with respect to
                  x,t at x,y.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 1 1990
*/

float dxt(p,x,y)
qnode_ptr_t p ;
int x, y ;

{ extern kernel_t C ;
  float d ;
  int i, h ;

  h = C.m/2 ;
  d = 0.0 ;

  if  ((*p->param_ptr[0])[p->res*x + y].err == NO_ERROR) {
    for (i = -h ; i <= h ; i++) {
      d += (*p->param_ptr[i+h])[p->res*x + y].fx*C.k[i+h] ;
    }
    d /= C.f ;
  }
  return(d) ;
}

/* 
           NAME : dxx(q,x,y,z)
   PARAMETER(S) : q : parameters of flow (image derivatives);
              x,y,z : image location for computation.
 
        PURPOSE : computes 2nd order partial derviative with respect
                  to x at x,y. 

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 1 1990
*/

float dxx(q,x,y,z)
qnode_ptr_t q ;
int x, y, z ;

{ extern kernel_t C ;
  float d ;
  int i, h ;

  h = C.m/2 ;
  d = 0.0 ;

  if (overflow(x,y,q->sizy,q->sizx,q->ofst,h,2)) {
    (*q->param_ptr[z])[q->res*x + y].err = SA_OVERFLOW ;
  }
  else {
    (*q->param_ptr[z])[q->res*x + y].err = NO_ERROR ;
    for (i = -h ; i <= h ; i++) {
      d += (*q->param_ptr[z])[q->res*(x+i) + y].fx*C.k[i+h] ;
    }
    d /= C.f ;
  }
  return(d) ;
}

/* 
           NAME : dxy(q,x,y,z)
   PARAMETER(S) : q : parameters of flow (image derivatives);
              x,y,z : image location for computation.
 
        PURPOSE : computes 2nd order partial derviative with respect
                  to x,y at x,y. 

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 1 1990
*/

float dxy(q,x,y,z)
qnode_ptr_t q ;
int x, y, z ;

{ extern kernel_t C ;
  float d ;
  int i, h ;

  h = C.m/2 ;
  d = 0.0 ;

  if (overflow(x,y,q->sizy,q->sizx,q->ofst,h,2)) {
    (*q->param_ptr[z])[q->res*x + y].err = SA_OVERFLOW ;
  }
  else {
    (*q->param_ptr[z])[q->res*x + y].err = NO_ERROR ;
    for (i = -h ; i <= h ; i++) {
      d += (*q->param_ptr[z])[q->res*(x+i) + y].fy*C.k[i+h] ;
    }
    d /= C.f ;
  }
  return(d) ;
}

/* 
           NAME : dy(p,q,x,y,z)
   PARAMETER(S) : p : image cube node;
                  q : parameters of flow (image derivatives);
              x,y,z : image location for computation.
 
        PURPOSE : computes 1st order partial derviative with respect
                  to y at x,y for image p->gauss_ptr[z].

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 1 1990
*/

float dy(p,q,x,y,z)
qnode_ptr_t p, q ;
int x, y, z ;

{ extern kernel_t C ;
  float d ;
  int i, h ;

  h = C.m/2 ;
  d = 0.0 ;

  if (overflow(x,y,p->sizy,p->sizx,p->ofst,h,1)) {
    (*q->param_ptr[z])[q->res*x + y].err = SA_OVERFLOW ;
  }
  else {
    (*q->param_ptr[z])[p->res*x + y].err = NO_ERROR ;
    for (i = -h ; i <= h ; i++) {
      d += (*p->gauss_ptr[z])[p->res*(x) + y + i]*C.k[i+h] ;
    }
    d /= C.f ;
  }
  return(d) ;
}

/* 
           NAME : dyt(p,x,y)
   PARAMETER(S) : p : image parameter node;
                x,y : image location for computation. 

        PURPOSE : computes partial derivative with respect to
                  y,t at x,y.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 1 1990
*/

float dyt(p,x,y)
qnode_ptr_t p ;
int x, y ;

{ extern kernel_t C ;
  float d ;
  int i, h ;

  h = C.m/2 ;
  d = 0.0 ;

  if  ((*p->param_ptr[0])[p->res*x + y].err == NO_ERROR) {
    for (i = -h ; i <= h ; i++) {
      d += (*p->param_ptr[i+h])[p->res*x + y].fy*C.k[i+h] ;
    }
    d /= C.f ;
  }
  return(d) ;
}

/* 
           NAME : dyy(q,x,y,z)
   PARAMETER(S) : q : parameters of flow (image derivatives);
              x,y,z : image location for computation.
 
        PURPOSE : computes 2nd order partial derviative with respect
                  to y at x,y. 

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 1 1990
*/

float dyy(q,x,y,z)
qnode_ptr_t q ;
int x, y, z ;

{ extern kernel_t C ;
  float d ;
  int i, h ;

  h = C.m/2 ;
  d = 0.0 ;

  if (overflow(x,y,q->sizy,q->sizx,q->ofst,h,2)) {
    (*q->param_ptr[z])[q->res*x + y].err = SA_OVERFLOW ;
  }
  else {
    (*q->param_ptr[z])[q->res*x + y].err = NO_ERROR ;
    for (i = -h ; i <= h ; i++) {
      d += (*q->param_ptr[z])[q->res*x + y+i].fy*C.k[i+h] ;
    }
    d /= C.f ;
  }
  return(d) ;
}

/* 
           NAME : usage()
   PARAMETER(S) : None.
 
        PURPOSE : Prints the appropriate usage and 
                  stops the run.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : June 25 1990
*/

usage()

{ 
  fprintf(stderr,"Usage: uras input-path output-path seq-name. [-SG n.n] [-TG n.n] [-M n] [-R 1 | -R 2 n.n] [-B rows cols] [-F n.n] [-C corr-vel-file nbins1 increment1 nbins2 increment2]\n") ;
  fprintf(stderr,"[-SG n.n]         : sigma value for spatial gaussian\n") ;
  fprintf(stderr,"                    default value is 3.0\n") ;
  fprintf(stderr,"[-TG n.n]         : sigma value for temporal gaussian\n") ; 
  fprintf(stderr,"                    default value is 1.5\n") ;
  fprintf(stderr,"[-M n]            : middle frame of image sequence\n") ;
  fprintf(stderr,"                    default value is 0\n") ;
  fprintf(stderr,"[-R 1 | -R 2 n.n] : n = 1 : applies TR1 regularization\n") ; 
  fprintf(stderr,"                    n = 2 : applies TR2 regularization\n") ;
  fprintf(stderr,"                    with sigma = n.n\n") ;
  fprintf(stderr,"[-B cols rows]    : for binary input files without header\n") ;
  fprintf(stderr,"[-F n.n]          : filters out unreliable estimates\n") ;
  fprintf(stderr,"                    using the Gaussian curv. (det(H))\n") ;
  fprintf(stderr,"[-C corr-vel-file nbins1 increment1 nbins2 increment2]\n") ;
  fprintf(stderr,"                  : error histogram\n") ;  
  fprintf(stderr,"                    corr-vel-file : file of correct velocities\n") ;
  fprintf(stderr,"                    nbins1, nbins2 : number of bins for histograms\n") ;
  fprintf(stderr,"                    increment1, increment2 : value of gap between bins\n") ;
  exit(-1) ;
}

/* 
           NAME : error(n)
   PARAMETER(S) : n : error number.
 
        PURPOSE : Prints the appropriate error message and 
                  stops the run.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : June 25 1990
*/

error(n)
int n ;

{ 
  switch(n) {
    
     case 1 : fprintf(stderr,"error %d : wrong number of arguments\n",n) ;
              break ;
     case 2 : fprintf(stderr,"error %d : invalid option\n",n) ;
              break ;
     case 3 : fprintf(stderr,"error %d : maximum number of frames is %d\n",n,Z) ;
              break ;
     case 4 : fprintf(stderr,"error %d : read error\n",n) ;
              exit(-1) ;
              break ;
     case 5 : fprintf(stderr,"error %d : header read error\n",n) ;
              exit(-1) ;
              break ;
     case 6 : fprintf(stderr,"error %d : file open error\n",n) ;
              exit(-1) ;
              break ;
     case 7 : fprintf(stderr,"error %d : write error\n",n) ;
              exit(-1) ;
              break ;
     case 8 : fprintf(stderr,"error %d : header write error\n",n) ;
              exit(-1) ;
              break ;
     case 9 : fprintf(stderr,"error %d : file create error\n",n) ;
              exit(-1) ;
              break ;
     case 10: fprintf(stderr,"error %d : maximum image size is %d by %d\n",n,X,Y) ;
              break ;
     case 11: fprintf(stderr,"error %d : option -R n mandatory\n",n) ;
              break ;
     case 12: fprintf(stderr,"error %d : no filtering with -R 2\n",n) ;
              break ;
     case 13: fprintf(stderr,"error %d : no histogram with -R 2\n",n) ;
              break ;
     case 14: fprintf(stderr,"error %d : undefined error\n",n) ; 
              break ;
     case 15: fprintf(stderr,"error %d : maximum number of bins is %d\n",
              n,N_BINS) ;
              break ;
    default : fprintf(stderr,"error %d : undefined error\n",n) ;
              exit(-1) ;
              break ;
  }
  usage() ;
}

/* 
           NAME : frames(s)
   PARAMETER(S) : s : Length of temporal Gaussian.
 
        PURPOSE : returns the number of frames necessary for
                  the smoothing and the computation of
                  4-point central differences.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 25 1992
*/

int frames(s)
int s ;

{ 
  if (s - 1 + SMSK > Z) {
    error(3) ;
  }
  else {
   return(s - 1  + SMSK) ;
  }
}

/* 
           NAME : start_number(middle,s)
   PARAMETER(S) : middle : number of middle frame;
                       s : length of temporal mask for derivatives;
 
        PURPOSE : computes the number of the start frame given
                  the middle frame and the length of the temporal
                  masks for derivatives. 
                 

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 25 1992
*/

int start_number(middle,s)
int middle, s ;

{ return(middle - (SMSK + s - 1)/2) ; }

/* 
           NAME : free_cube(p,n) ;
   PARAMETER(S) : p : pointer on an image cube node;
                  n : number of frames.
 
        PURPOSE : deallocates n images to p->gauss_ptr[i].

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : September 15 1990
*/

free_cube(p,n)
qnode_ptr_t p ;

{ int i ;

  for (i = 0 ; i < n ; i++) {
    free((image512_t *)p->gauss_ptr[i]) ;
  }
}

/* 
           NAME : free_flow(p,n) ;
   PARAMETER(S) : p : pointer on an flow and parameter node.
 
        PURPOSE : deallocates a flow field to p->flow_ptr
                  and n arrays of flow parameters to 
                  p->param_ptr[i].

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : September 15 1990
*/

free_flow(p,n)
qnode_ptr_t p ;
int n ;

{ int i ;

  free((disp_field512_t *)p->flow_ptr);
  for (i = 0 ; i < n ; i++) {
    free((param512_t *)p->param_ptr[i]) ;
  }
}

/* 
           NAME : generate_gauss(ker,sig)
   PARAMETER(S) : ker : 1D  kernel;
                  sig : sigma.
 
        PURPOSE : Inits the kernel values with G(x,sig).

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 6 1990

*/

generate_gauss(ker,sig)
kernel_t *ker ;
float sig ;

{ int i, h ;

  h = (int)(sig*6.0 + 1.0) ;
  if (h%2 == 0) {
    h++ ;
  }
  (*ker).f = 0.0 ;
  if (sig != 0.0) {
    for (i = -(int)(h/2.0) ; i <= (int)(h/2.0) ; i++) {
      (*ker).k[i+(int)(h/2.0)] = exp(-pow((float)i,2.0)/(2.0*pow(sig,2.0)))/
      (sqrt(2.0*PI)*sig) ;
      (*ker).f += (*ker).k[i+(int)(h/2.0)] ;
    }
  }
  else {
    (*ker).k[0] = 1.0 ;
    (*ker).f = 1.0 ;
  }
  (*ker).m = h ; 
}

/* 
           NAME : init_beaudet(Ix,Ixx,Ixy)
   PARAMETER(S) :  Ix : operator for 1st order x derivative;
                  Ixx : operator for 2nd order x derivative;
                  Ixy : operator for partial x and y derivatives.
 
        PURPOSE : Initialize the image operators defined by Beaudet
                  to obtain numerical derivatives.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : September 13 1990
*/

init_beaudet(Ix,Ixx,Ixy)
beaudet_t *Ix, *Ixx, *Ixy ;

{
  (*Ix).g[0][0] = -2 ; (*Ix).g[1][0] = -1 ; (*Ix).g[2][0] = 0 ; 
  (*Ix).g[3][0] = 1 ;  (*Ix).g[4][0] = 2 ;
  (*Ix).g[0][1] = -2 ; (*Ix).g[1][1] = -1 ; (*Ix).g[2][1] = 0 ; 
  (*Ix).g[3][1] = 1 ;  (*Ix).g[4][1] = 2 ;
  (*Ix).g[0][2] = -2 ; (*Ix).g[1][2] = -1 ; (*Ix).g[2][2] = 0 ; 
  (*Ix).g[3][2] = 1 ;  (*Ix).g[4][2] = 2 ;
  (*Ix).g[0][3] = -2 ; (*Ix).g[1][3] = -1 ; (*Ix).g[2][3] = 0 ; 
  (*Ix).g[3][3] = 1 ;  (*Ix).g[4][3] = 2 ;
  (*Ix).g[0][4] = -2 ; (*Ix).g[1][4] = -1 ; (*Ix).g[2][4] = 0 ; 
  (*Ix).g[3][4] = 1 ;  (*Ix).g[4][4] = 2 ;
  (*Ix).f = 50.0 ;
  (*Ix).m = 5 ;
      
  (*Ixx).g[0][0] = 2 ;  (*Ixx).g[1][0] = -1 ;  (*Ixx).g[2][0] = -2 ; 
  (*Ixx).g[3][0] = -1 ;  (*Ixx).g[4][0] = 2 ;
  (*Ixx).g[0][1] = 2 ; (*Ixx).g[1][1] = -1 ; (*Ixx).g[2][1] = -2; 
  (*Ixx).g[3][1] = -1 ; (*Ixx).g[4][1] = 2 ;
  (*Ixx).g[0][2] = 2 ; (*Ixx).g[1][2] = -1 ; (*Ixx).g[2][2] = -2 ; 
  (*Ixx).g[3][2] = -1 ; (*Ixx).g[4][2] = 2 ;
  (*Ixx).g[0][3] = 2 ; (*Ixx).g[1][3] = -1 ; (*Ixx).g[2][3] = -2 ; 
  (*Ixx).g[3][3] = -1 ; (*Ixx).g[4][3] = 2 ;
  (*Ixx).g[0][4] = 2 ;  (*Ixx).g[1][4] = -1 ;  (*Ixx).g[2][4] = -2 ; 
  (*Ixx).g[3][4] = -1 ;  (*Ixx).g[4][4] = 2 ;
  (*Ixx).f = 35.0 ;
  (*Ixx).m = 5 ;

  (*Ixy).g[0][0] = 4 ;  (*Ixy).g[1][0] = 2 ;  (*Ixy).g[2][0] = 0 ; 
  (*Ixy).g[3][0] = -2 ; (*Ixy).g[4][0] = -4 ;
  (*Ixy).g[0][1] = 2 ;  (*Ixy).g[1][1] = 1 ;  (*Ixy).g[2][1] = 0 ; 
  (*Ixy).g[3][1] = -1 ; (*Ixy).g[4][1] = -2 ;
  (*Ixy).g[0][2] = 0 ;  (*Ixy).g[1][2] = 0 ;  (*Ixy).g[2][2] = 0 ; 
  (*Ixy).g[3][2] = 0 ;  (*Ixy).g[4][2] = 0 ;
  (*Ixy).g[0][3] = -2 ; (*Ixy).g[1][3] = -1 ; (*Ixy).g[2][3] = 0 ; 
  (*Ixy).g[3][3] = 1 ;  (*Ixy).g[4][3] = 2 ;
  (*Ixy).g[0][4] = -4 ; (*Ixy).g[1][4] = -2 ; (*Ixy).g[2][4] = 0 ; 
  (*Ixy).g[3][4] = 2 ;  (*Ixy).g[4][4] = 4 ;
  (*Ixy).f = 100.0 ;
  (*Ixy).m = 5 ;
}

/* 
           NAME : init_central(C)
   PARAMETER(S) : C : 1-D mask contraining 4-point central difference
                      factors.
 
        PURPOSE : Inits a 4-point central difference mask. 

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 25 1992
*/

init_central(C)
kernel_t *C ;

{
  (*C).k[0] = -1.0 ; (*C).k[1] = 8.0 ; 
  (*C).k[2] = 0.0 ; (*C).k[3] = -8.0 ; (*C).k[4] = 1.0 ;
  (*C).m = 5 ;
  (*C).f = 12.0 ;
}

/* 
           NAME : init_cube(p,n) ;
   PARAMETER(S) : p : pointer on an image cube node;
                  n : number of image frames to allocate.
 
        PURPOSE : allocates n images to p->gauss_ptr[i].

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : September 15 1990
*/

init_cube(p,n)
qnode_ptr_t p ;
int n ;

{ int i ;

  for (i = 0 ; i < n ; i++) { 
    p->gauss_ptr[i] = (image512_t *)malloc(sizeof(image512_t)) ;
  }
}

/* 
           NAME : init_flow(p,n) ;
   PARAMETER(S) : p : pointer on an parameter cube node;
                  n : number of parameter frames to allocate.
 
        PURPOSE : allocates n parameter frames to p->param_ptr[i]
                  and a flow field to p->flow_ptr.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : September 15 1990
*/

init_flow(p,n)
qnode_ptr_t p ;
int n ;

{ int i ;

  p->flow_ptr = (disp_field512_t *)malloc(sizeof(disp_field512_t)) ;
  for (i = 0 ; i < n ; i++) { 
    p->param_ptr[i] = (param512_t *)malloc(sizeof(param512_t)) ;
  }
}

/* 
           NAME : init_list(h,q) 
   PARAMETER(S) : h : head list pointer;
                  q : tail list pointer.
 
        PURPOSE : Initialization of the head
                  and tail pointers to NULL.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 17 1990

*/

init_list(h,q)
qnode_ptr_t *h, *q ;

{ 
  *h = (qnode_ptr_t)NULL ;
  *q = (qnode_ptr_t)NULL ;
}

/* 
           NAME : insert_node(p,h,q) 
   PARAMETER(S) : p : pointer on the node to insert;
                  h : head pointer of the list;
                  q : tail pointer of the list.
 
        PURPOSE : Insertion of the node p in a doubly 
                  linked list in order of the level
                  numbers (p->level).

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 17 1990

*/

insert_node(p,h,q) 
qnode_ptr_t p, *h, *q ;

{ qnode_ptr_t s, t ;

  t = *h ;
  s = (qnode_ptr_t)NULL ;
  while ((t != (qnode_ptr_t)NULL) && (p->level < t->level)) {
    s = t ;
    t = t->forth ;
  }
  if (s != (qnode_ptr_t)NULL) {
    s->forth = p ;
  }
  else {
    *h = p ;
  }
  if (t != (qnode_ptr_t)NULL) {
    t->back = p ;
  }
  else {
    *q = p ;
  }
  p->back = s ;
  p->forth = t ;
}


/* 
           NAME : itoa(n,s)
   PARAMETER(S) : n : integer value ;
                  s : output string.
 
        PURPOSE : Converts integer into string.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : June 25 1990
*/

itoa(n,s)
int n ;
string s ; 

{ int i, j, digit, neg ;
  string t ;

  if (n < 0) { 
    neg = TRUE ;
  }
  else {
    neg = FALSE ;
  }
  if (n == 0) {
    s[0] = '0' ;
    s[1] = '\0' ;
  }
  else {
    i = 0 ;
    while (abs(n) > 0) {
      digit = abs(n) % 10 ;
      t[i] = (char)((int)'0' + digit) ;
      n = n / 10 ;
      i++ ;
    }
    if (neg) {
      t[i] = '-' ; 
      i++ ;
    }
    for (j = 0 ; j < i ; j++) {
      s[j] = t[i - j - 1] ;
    }
    s[i] = '\0' ;
  }
}

/* 
           NAME : mag(u)
   PARAMETER(S) : u : 2D flow vector.
 
        PURPOSE : returns the magnitude of u.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 1 1990
*/

float mag(u)
disp_vect_t u ;

{ 
  return(sqrt(pow(u.x,2.0) + pow(u.y,2.0))) ;
}

/* 
           NAME : min_cond(s,f,n,m,trsh,sf)
   PARAMETER(S) : s : 1D array containing at most n image flow locations;
                  f : image parameter node;
                  n : maximum number of image locations in s;
                  m : actual number of image locations in s;
               trsh : threshold for condition number;
                 sf : thresholding boolean value.
 
        PURPOSE : determines among the locations in s, which one has a
                  minimal condition number.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 4 1992
*/

pos_t min_cond(s,f,n,m,trsh,sf)
pos_t s[X] ;
qnode_ptr_t f ;
int n, m ;
float trsh ;
int sf ;

{ pos_t loc, err_loc ;
  float min_c, gauss ;
  int l ;

  min_c = MAX_COND + 1 ;
  loc.i = -1 ;
  loc.j = -1 ;
  err_loc.i = -1 ;
  err_loc.j = -1 ;
  if (m < n) {
    n = m ;
  }
  for (l = 0 ; l < n ; l++) {
    if ((*f->param_ptr[f->sizz/2])[f->res*s[l].i + s[l].j].cond < min_c) {
      min_c = (*f->param_ptr[f->sizz/2])[f->res*s[l].i + s[l].j].cond ;
      gauss = (*f->param_ptr[f->sizz/2])[f->res*s[l].i + s[l].j].gauss ;
      loc = s[l] ;
    }
  }
  if (sf) {
    if (gauss > trsh) {
      return(loc) ;
    }
    else {
      return(err_loc) ;
    }
  }
  else {
    return(loc) ;
  }
}

/* 
           NAME : sort(discr_val,sample,num,n)
   PARAMETER(S) : discr_val : discriminant values for image positions
                              contained in sample;
                  sample    : image positions to be sorted with respect
                              to their discriminant values;
                  num       : number of elements in sample;
                  n         : image region size.
 
        PURPOSE : sorts the n first elements of sample[1..m] with respect
                  to the discriminant value.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : September 15 1990
*/

sort(discr_val,sample,num,n) 
float discr_val[X] ;
pos_t sample[X] ;
int num, n ;

{ pos_t ts ; 
  float td ;
  int i, j, index ;

  if (num < n) {
    n = num ;
  }
  for (i = 0 ; i < n - 1 ; i++) {
    index = i ;
    for (j = index + 1 ; j < num ; j++) {
      if (discr_val[index] > discr_val[j]) {
        index = j ;
      }
    }
    td = discr_val[index] ;
    discr_val[index] = discr_val[i] ;
    discr_val[i] = td ;
   
    ts = sample[index] ;
    sample[index] = sample[i] ;
    sample[i] = ts ;
  }
}

/* 
           NAME : min_discr(f,i,j,n,sample,num)
   PARAMETER(S) :   f      : image parameter node;
                  i,j      : upper left corner location of a n*n image area;
                    n      : image area size;
                    sample : valid image locations of image area;
                    num    : number of image locations in sample.
 
        PURPOSE : extracts the n image locations from the n*n image area
                  for which the discriminant value is minimal.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : November 7 1990
*/

min_discr(f,i,j,n,sample,num)
qnode_ptr_t f ;
int i, j, n ;
pos_t sample[X] ;
int *num ;

{ float discr_val[X] ;
  int k, l ;

  *num = 0 ; 
  for (k = 0 ; k < n ; k++) {
    for (l = 0 ; l < n ; l++) {
      if ((*f->param_ptr[f->sizz/2])[f->res*(i+k) + j+l].err == NO_ERROR) {
        sample[*num].i = i + k ;
        sample[*num].j = j + l ;
        discr_val[*num] = (*f->param_ptr[f->sizz/2])[f->res*(i+k) + j+l].discr ;
        (*num)++ ;
      }
    }
  }
  sort(discr_val,sample,*num,n) ;
}

/* 
           NAME : norm(V,n)
   PARAMETER(S) : V : vector;
                  n : length of V.
 
        PURPOSE : returns L2 norm of V

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : November 7 1990
*/

float norm(V,n) 
float V[] ;
int n ;

{ int i ;
  float sum ; 

  sum = 0.0 ;
  for (i = 0 ; i < n ; i++) {
    sum += (V[i]*V[i]) ;
  }
  return(sqrt(sum)) ;
}

/* 
           NAME : pgetrast(fn,hd,bf,sx,sy,r) 
   PARAMETER(S) : fn    : filename;
                  hd    : image header (raster file);
                  bf    : Pointer on a 2D array containing the image;
                  sx,sy : image size to read;
                   r    : resolution of the array.

        PURPOSE : Reads rasterfile specified
                  by fn into bf.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 27 1989, modified March 20 1990.

*/

pgetrast(fn,hd,bf,sx,sy,r) 
char *fn ;
unsigned char hd[H] ;
float bf[X*Y] ;
int sx, sy, r ;

{ int fd, i, j ;
  unsigned char c ;

  if ((fd = open(fn,O_RDONLY)) > 0) { ;
    if (read(fd,hd,H) == H) { 
      for (i = 0 ; i < sy ; i++) {
        for (j = 0 ; j < sx ; j++) {
          if (read(fd,&c,sizeof(c)) != 1) {
            error(4) ;
          }
          bf[r*i + j] = (float)c ;
        }
      }
    }
    else {
      error(5) ;
    }
  }
  else {
    error(6) ;
  }
  close(fd) ;
}

/* 
           NAME : Bpgetrast(fn,bf,sx,sy,r) 
   PARAMETER(S) : fn    : filename;
                  bf    : Pointer on a 2D array containing the image;
                  sx,sy : image size to read;
                   r    : resolution of the array.

        PURPOSE : Reads rasterfile without header specified
                  by fn into bf.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 27 1989, modified March 20 1990.

*/

Bpgetrast(fn,bf,sx,sy,r) 
char *fn ;
float bf[X*Y] ;
int sx, sy, r ;

{ int fd, i, j ;
  unsigned char c ;

  if ((fd = open(fn,O_RDONLY)) > 0) { ;
    for (i = 0 ; i < sy ; i++) {
      for (j = 0 ; j < sx ; j++) {
        if (read(fd,&c,sizeof(c)) != 1) {
          error(4) ;
        }
        bf[r*i + j] = (float)c ;
      }
    }
  }
  else {
    error(6) ;
  }
  close(fd) ;
}

/* 
           NAME : pputrast(fn,hd,bf,sx,sy,r) 
   PARAMETER(S) : fn    : filename;
                  hd    : image header (raster file);
                  bf    : Pointer on a 2D array containing the image;
                  sx,sy : image size to write;
                  r     : resolution of the array.
 
        PURPOSE : Writes the contents of bf into a 
                  rasterfile specified by fn.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 27 1989, modified March 20 1990.

*/

pputrast(fn,hd,bf,sx,sy,r) 
char *fn ;
unsigned char hd[H] ;
float bf[X*Y] ;
int sx, sy, r ;

{ int fd, i, j ;
  unsigned char c ;

  if ((fd = creat(fn,0600)) > 0) { ;
    if (write(fd,hd,H) == H) { 
      for (i = 0 ; i < sy ; i++) { 
        for (j = 0 ; j < sx ; j++) { 
          c = (unsigned char)bf[r*i + j] ;
          if (write(fd,&c,sizeof(c)) != 1) {
            error(7) ;
          }
        }
      }
    }
    else {
      error(8) ; 
    }
  }
  else {
    error(9) ;
  }
  close(fd) ;
}

/* 
           NAME : propagate(u,det,cn,q,i,j,n)
   PARAMETER(S) :   u : flow vector to be propagated;
                  det : determinant;
                   cn : condition number;
                    q : flow field for propagation;
                  i,j : coordinates of upper left corner of image area;
                    n : image area size.
 
        PURPOSE : propagates u, det, cn in a n*n square image area with upper 
                  left corner coordinates (i,j).

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : November 7 1990
*/

propagate(u,det,cn,q,i,j,n)
disp_vect_t u ;
float det, cn ;
qnode_ptr_t q ;
int i, j, n ;

{ int k, l ;

  for (k = 0 ; k < n ; k++) {
    for (l = 0 ; l < n ; l++) {
      (*q->flow_ptr)[q->res*(i+k) + j+l] = u ;
      (*q->param_ptr[q->sizz/2])[q->res*(i+k) + j+l].gauss = det ;
      (*q->param_ptr[q->sizz/2])[q->res*(i+k) + j+l].cond = cn ;
    }
  }
}

/* 
           NAME : psi_error(ve,va)
   PARAMETER(S) : ve : estimated flow vector;
                  va : accurate flow vector.
 
        PURPOSE : computes angular error of ve with respect to va
                  using Fleet [90] angular error metric

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : November 7 1990
*/

float psi_error(ve,va)
disp_vect_t ve, va ;

{ float norm(), v ;
  float VE[3], VA[3] ;

  VE[0] = ve.x ;
  VE[1] = ve.y ;
  VE[2] = 1.0 ;

  VA[0] = va.x ;
  VA[1] = va.y ;
  VA[2] = 1.0 ;

  v = (VE[0]*VA[0] + VE[1]*VA[1] + 1.0)/(norm(VA,3)*norm(VE,3)) ;
  if ((v > 1.0) && (v < 1.0001)) {
    v = 1.0 ;
  }
  return((float)(acos(v))*180.0/PI) ;
}

/* 
           NAME : raster_size(fn,num,x,y)
   PARAMETER(S) : fn : raster file name; 
                 num : frame number;
                   x : raster width;
                   y : raster height.
 
        PURPOSE : Reads the header of raterfile fn to get the
                  size. 

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 22 1992
*/

raster_size(fn,num,x,y)
char *fn ;
int num, *x, *y ;

{ raster_t hd ;
  string s1, s2 ;
  int fd ;

  itoa(num,s1) ;
  concat(fn,s1,s2) ;
  if ((fd = open(s2,O_RDONLY)) > 0) {
    if (read(fd,&hd,sizeof(hd)) == sizeof(hd)) {
      *x = hd.width ;
      *y = hd.height ;
      if (*x > X || *y > Y) {
        error(10) ;
      }
    }
  }
  else {
    error(6) ;
  }
  close(fd) ;
}

/* 
           NAME : binary(row,col,x,y)
   PARAMETER(S) : row : number of rows in input image;
                  col : number of cols in input image;
                    x : raster width;
                    y : raster height.
 
        PURPOSE : sets the size of images when option -B used.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 22 1992
*/

binary_size(row,col,x,y)
int row, col, *x, *y ;

{ *x = col ;
  *y = row ;
  if ((*x > X) || (*y > Y)) {
    error(10) ;
  }
}

/* 
           NAME : vector(loc,f)
   PARAMETER(S) : loc : flow field location;
                    f : node pointer.
 
        PURPOSE : returns the vector located at loc in f.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : November 7 1990
*/

disp_vect_t vector(loc,f)
pos_t loc ;
qnode_ptr_t f ;

{ disp_vect_t u ;
  
  if ((loc.i != -1) && (loc.j != -1)) {
    u = (*f->flow_ptr)[f->res*loc.i + loc.j] ;
  }
  else {
    u.x = -100.0 ;
    u.y = 100.0 ;
  }
  return(u) ;
}

/* 
           NAME : regul1(f,n,sf,trsh)
   PARAMETER(S) : f : node pointer;
                  n : image area size;
                 sf : thresholding boolean value;
               trsh : threshold for Gaussian curvature.
 
        PURPOSE : operates procedure 1 (TR1 described in Uras et al.[88])
                  for regularization of flow field f.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : November 7 1990
*/

regul1(f,n,sf,trsh)
qnode_ptr_t f ;
int n, sf ;
float trsh ;

{ pos_t min_cond(), sample[X], rep_loc ;
  disp_vect_t vector(), u ;
  float condition(), determinant(), cond_num, det ;
  int num, xf, yf, i, j, k ;

  xf = (int)(((f->sizx-f->ofst-NRADIUS-(f->ofst+NRADIUS))%n)/2.0 + 0.5) ;
  yf = (int)(((f->sizy-f->ofst-NRADIUS-(f->ofst+NRADIUS))%n)/2.0 + 0.5) ;
  KERNEL_X += xf ;
  KERNEL_Y += yf ;

  for (i = f->ofst+NRADIUS+yf ; i < f->sizy-f->ofst-NRADIUS-yf ; i += n) {
    for (j = f->ofst+NRADIUS+xf ; j < f->sizx-f->ofst-NRADIUS-xf ; j += n) {
      min_discr(f,i,j,n,sample,&num) ;
      rep_loc = min_cond(sample,f,n,num,trsh,sf) ;
      u = vector(rep_loc,f) ;
      cond_num = condition(rep_loc,f) ;
      det = determinant(rep_loc,f) ;
      propagate(u,det,cond_num,f,i,j,n) ;
    }
  }
}

/* 
           NAME : regul2(fl,ker)
   PARAMETER(S) :  fl : node pointer;
                  ker : 1D kernel for convolution of flow.
   
 
        PURPOSE : operates procedure 2 (TR2 described in Uras et al.[88])
                  on flow field fl.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : November 7 1990
*/

regul2(fl,ker)
qnode_ptr_t fl ;
kernel_t ker ;

{ disp_field512_t *t ;
  disp_vect_t s ;
  int h, i, j, k ;

  h = ker.m/2 ;
  t = (disp_field512_t *)malloc(sizeof(disp_field512_t)) ;

  for (i = fl->ofst+NRADIUS ; i < fl->sizy-fl->ofst-NRADIUS ; i++) {
    for (j = fl->ofst+NRADIUS+h ; j < fl->sizx-fl->ofst-NRADIUS-h ; j++) {
      s.x = 0.0 ;
      s.y = 0.0 ;
      for (k = 0 ; k < ker.m ; k++) {
        s.x = s.x + (*fl->flow_ptr)[fl->res*i + j+k-h].x*ker.k[k] ;
        s.y = s.y + (*fl->flow_ptr)[fl->res*i + j+k-h].y*ker.k[k] ;
      }
      (*t)[fl->res*i + j].x = s.x/ker.f ;
      (*t)[fl->res*i + j].y = s.y/ker.f ;
    }
  }

     
  for (i = fl->ofst+NRADIUS+h ; i < fl->sizy-fl->ofst-NRADIUS-h ; i++) {
    for (j = fl->ofst+NRADIUS ; j < fl->sizx-fl->ofst-NRADIUS ; j++) {
      s.x = 0.0 ;
      s.y = 0.0 ;
      for (k = 0 ; k < ker.m ; k++) {
        s.x = s.x + (*t)[fl->res*(i+k-h) + j].x*ker.k[k] ;
        s.y = s.y + (*t)[fl->res*(i+k-h) + j].y*ker.k[k] ;
      }
      (*fl->flow_ptr)[fl->res*i + j].x = s.x/ker.f ;
      (*fl->flow_ptr)[fl->res*i + j].y = s.y/ker.f ;
    }
  }
  free((disp_field512_t *)t) ;
}

/*
           NAME : write_velocity(fn,flow_field)
   PARAMETER(S) : flow_field : the 2D array containing velocities 
                               or displacements.

        PURPOSE : Output field using Travis Burkitt's format.

         AUTHOR : Travis Burkitt, updated by Steven Beauchemin
             AT : University of Western Ontario
           DATE : May 7 1990
*/

write_velocity(fn,p)
qnode_ptr_t p ;
char *fn ;

{ extern int KERNEL_X, KERNEL_Y ;
  float x, y ; 
  int i, j, fdf, bytes ;

  if ((fdf=creat(fn,0600)) < 1) {
    error(6) ;
  }
  
  x = p->sizx ;
  y = p->sizy ;
  write(fdf,&x,4) ;
  write(fdf,&y,4) ;
  
  x = ((p->sizx-KERNEL_X-NRADIUS)-(KERNEL_X+NRADIUS)+SKIP-1)/SKIP ;
  y = ((p->sizy-KERNEL_Y-NRADIUS)-(KERNEL_Y+NRADIUS)+SKIP-1)/SKIP ;
  write(fdf,&x,4);
  write(fdf,&y,4);
  
  x = (KERNEL_X+NRADIUS+SKIP-1)/SKIP ;
  y = (KERNEL_Y+NRADIUS+SKIP-1)/SKIP ;
  write(fdf,&x,4) ;
  write(fdf,&y,4) ;
  bytes = 24 ;
  
  for(i = KERNEL_Y + NRADIUS ; i < p->sizy - KERNEL_Y - NRADIUS ; i++) {
    for(j = KERNEL_X + NRADIUS ; j < p->sizx - KERNEL_X - NRADIUS ; j++) {
      x = (*p->flow_ptr)[p->res*i + j].y ;
      y = -(*p->flow_ptr)[p->res*i + j].x ;
      write(fdf,&x,4) ;
      write(fdf,&y,4) ;
      bytes += 8 ;
    }
  }
  close(fdf) ;
}

/* 
           NAME : valid_option(argc,argv,in_path,out_path,sigma1,sigma2,
                  sigma3,histo,re,sf,trsh,i_fname,v_fname,c_fname,h_fname,
                  nbin_1,incr_1,nbin_2,incr_2,n_frame,binary,row,col) 

   PARAMETER(S) : argc : argument count;
                  argv : argument values;
               in_path : path name for input data;
              out_path : path name for output data;
                sigma1 : spatial sigma;
                sigma2 : temporal sigma;
                sigma3 : flow field sigma (TR2);
                 histo : error histogram option;
                    re : regularization option;
                    sf : filtering option;
                  trsh : threshold value for filtering;
               i_fname : input filename;
               v_fname : velocity filename;
               c_fname : correct velocities file name;
               h_fname : hisogram data file name ;
                nbin_1 : number of bins in histo 1;
                incr_1 : increment in histo 1 ;
                nbin_2 : number of bins in histo 2 ;
                incr_2 : increment in histo2;
               n_frame : number of required frames;
                binary : input files without header;
                   row : number of rows in input files;
                   col : number of cols in input files.
 
        PURPOSE : Validates line arguments.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : June 25 1990
*/

valid_option(argc,argv,in_path,out_path,sigma1,sigma2,sigma3,histo,re,sf,trsh,
             i_fname,v_fname,c_fname,h_fname,nbin_1,incr_1,nbin_2,incr_2,
             n_frame,binary,row,col) 
int argc ;
char *argv[] ;
int *nbin_1, *nbin_2, *n_frame, *histo, *re, *sf, *binary, *row, *col ;
float *incr_1, *incr_2, *sigma1, *sigma2, *sigma3, *trsh ;
string in_path, out_path, i_fname, v_fname, c_fname, h_fname ;

{ int i ;
  string ext, t0, str_sg, str_tg, str_r2, str_trsh ;

  if (argc == 1) {
    usage() ;
  }
  if ((argc <= 23) && (argc >= 4)) {
    *binary = FALSE ;
    *re = TR1 ;
    *sf = FALSE ;
    *histo = FALSE ;
    *trsh = 0.0 ;
    *sigma1 = DEF_S1 ;
    strcpy(str_sg,"3.0") ;
    *sigma2 = DEF_S2 ;
    strcpy(str_tg,"1.5") ;
    *sigma3 = DEF_S3 ; 
    strcpy(str_r2,"3.0") ;
    *n_frame = 0 ;
    *nbin_1 = 0 ;
    *nbin_2 = 0 ;
    strcpy(i_fname,"Undefined") ;
    strcpy(v_fname,"Undefined") ;
    strcpy(c_fname,"Undefined") ;
    strcpy(h_fname,"Undefined") ;
  }
  else {
    error(1) ;
  }
  i = 4 ;
  while (i < argc) {
    if (strcmp("-R",argv[i]) == 0) {
      if (i + 1 < argc) {
        sscanf(argv[i+1],"%d",re) ;
        i += 2 ;
        if (*re == TR2) {
          if (i < argc) {
            sscanf(argv[i],"%f",sigma3) ;
            strcpy(str_r2,argv[i]) ;
            i += 1 ;
          }
          else {
            error(1) ;
          }
        }
      }
      else {
        error(1) ;
      }
    }
    else {
      if (strcmp("-F",argv[i]) == 0) {
        if (i + 1 < argc) {
          sscanf(argv[i+1],"%f",trsh) ;
          strcpy(str_trsh,argv[i+1]) ;
          *sf = TRUE ;
          i += 2 ;
        }
        else {
          error(1) ;
        }
      }
      else {
        if (strcmp("-C",argv[i]) == 0) {
          if (i + 5 < argc) {
            strcpy(c_fname,argv[i+1]) ;
            sscanf(argv[i+2],"%d",nbin_1) ;
            sscanf(argv[i+3],"%f",incr_1) ;
            sscanf(argv[i+4],"%d",nbin_2) ;
            sscanf(argv[i+5],"%f",incr_2) ;
            *histo = TRUE ;
            i += 6 ;
          }
          else {
            error(1) ;
          }
        } 
        else {
          if (strcmp("-SG", argv[i]) == 0) {
            if (i + 1 < argc) {
              sscanf(argv[i+1],"%f",sigma1) ;
              strcpy(str_sg,argv[i+1]) ;
              i += 2 ;
            }
            else {
              error(1) ;
            }
          }
          else {
            if (strcmp("-TG",argv[i]) == 0) {
              if (i + 1 < argc) {
                sscanf(argv[i+1],"%f",sigma2) ;
                strcpy(str_tg,argv[i+1]) ;
                i += 2 ;
              }
              else {
                error(1) ;
              } 
            }
            else {
              if (strcmp("-M",argv[i]) == 0) {
                if (i + 1 < argc) {
                  sscanf(argv[i+1],"%d",n_frame) ;
                  i += 2 ;
                }
                else {
                  error(1) ;
                }
              }
              else {
                if (strcmp("-B",argv[i]) == 0) { 
                  if (i + 2 < argc) {
                    sscanf(argv[i+1],"%d",col) ;
                    sscanf(argv[i+2],"%d",row) ;
                    *binary = TRUE ;
                    i += 3 ;
                  }
                  else {
                    error(1) ;
                  }
                }
                else {
                  error(2) ;
                }  
              }
            }
          }
        }
      }
    }
  }
  if ((*sf) && (*re == TR2)) {
    error(12) ;
  }
  if ((*nbin_1 > N_BINS) || (*nbin_2 > N_BINS)) {
    error(15) ;
  }
  if ((*re == TR2) && (*histo)) {
    error(13) ;
  }
  if ((*re != TR1) && (*re != TR2) && (*re != NO_REG)) {
    error(2) ;
  }
  concat(argv[1],"/",in_path) ;
  concat(argv[2],"/",out_path) ;
  concat(in_path,argv[3],i_fname) ;
  strcpy(ext,"-sg") ;
  concat(ext,str_sg,ext) ;
  concat(ext,"-tg",ext) ;
  concat(ext,str_tg,ext) ;
  concat(ext,"-m",ext) ;
  itoa(*n_frame,t0) ;
  concat(ext,t0,ext) ;
  concat(ext,"-r",ext) ;
  if (*re == TR1) {
    concat(ext,"1",ext) ;
  }
  else {
    concat(ext,"2-",ext) ;
    concat(ext,str_r2,ext) ;
  }
  if (*sf) {
    concat(ext,"-f",ext) ;
    concat(ext,str_trsh,ext) ;
  }
  concat(out_path,argv[0],v_fname) ;
  concat(v_fname,".",v_fname) ;
  concat(v_fname,argv[3],v_fname) ;
  concat(v_fname,"F",v_fname) ;
  concat(v_fname,ext,v_fname) ;
  if (*histo) {
    concat(out_path,argv[0],h_fname) ;
    concat(h_fname,".",h_fname) ;
    concat(h_fname,argv[3],h_fname) ;
    concat(h_fname,"H",h_fname) ;
    concat(h_fname,ext,h_fname) ;
  }
  printf(" Generated velocity filename  : %s\n",v_fname) ;
  printf(" Generated histogram filename : %s\n",h_fname) ;
  printf(" Input stem name              : %s\n",i_fname) ;
  printf(" Correct velocity filename    : %s\n",c_fname) ;
  fflush(stdout) ;
}

/* 
           NAME : prod_histo(f,c_fn,h_fn,nbin_1,incr_1,nbin_2,incr_2,
                             ttl_err,ttl_std,dens)
   PARAMETER(S) :    f : flow field with parameters;
                  c_fn : correct velocities file name;
                  h_fn : histogram file name;
                nbin_1 : number of bins in histo 1;
                incr_1 : increment per bin in histo 1;
                nbin_2 : number of bins in histo 2;
                incr_2 : increment per bin in histo 2;
               ttl_err : error in flow field;
               ttl_std : standard deviation;
                  dens : density of flow field.
 
        PURPOSE : produces an error histogram
                  h1: error vs determinant;
                  h2: error vs condition number;
                  h3: h1 cumulated;
                  h4: h2 cumulated.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : February 20 1992
*/

prod_histo(f,c_fn,h_fn,nbin_1,incr_1,nbin_2,incr_2,ttl_err,ttl_std,dens)
qnode_ptr_t f ;
string c_fn, h_fn ;
float incr_1, incr_2, *ttl_err, *ttl_std, *dens ;
int nbin_1, nbin_2 ;

{ qnode_ptr_t create_node(), q ;
  disp_vect_t u ;
  histo_t histo1[N_BINS], histo2[N_BINS], c_histo1[N_BINS], c_histo2[N_BINS] ;
  float max_det, max_cond, psi_error(), actual_x, actual_y, size_x, size_y, 
        offset_x, offset_y, det, cond, density, err, avg_err, t, x, y ;
  int fdf, nbytes, i, j, k, index1, index2, ttl_freq, abs_freq ;
  FILE *fdp ;
  extern int KERNEL_X, KERNEL_Y ;

  q = create_node(0,f->res,f->sizx,f->sizy,f->sizz,0) ;
  q->flow_ptr = (disp_field512_t *)malloc(sizeof(disp_field512_t)) ;
  if ((fdf = open(c_fn,O_RDONLY)) <= 0) {
    error(6) ;
  }

  nbytes = 0 ;
  nbytes += read(fdf,&actual_x,sizeof(float)) ;
  nbytes += read(fdf,&actual_y,sizeof(float)) ;
  nbytes += read(fdf,&size_x,sizeof(float)) ;
  nbytes += read(fdf,&size_y,sizeof(float)) ;
  nbytes += read(fdf,&offset_x,sizeof(float)) ;
  nbytes += read(fdf,&offset_y,sizeof(float)) ;
  
  for (i = (int)offset_y ; i < (int)actual_y ; i++) {
    for (j = (int)offset_x ; j < (int)actual_x ; j++) {
      nbytes += read(fdf,&y,sizeof(float)) ;
      nbytes += read(fdf,&x,sizeof(float)) ;
      (*q->flow_ptr)[q->res*i + j].y = y ;
      (*q->flow_ptr)[q->res*i + j].x = -x ;
    }
  }
  close(fdf) ;

  max_det = incr_1*(float)nbin_1 ;
  max_cond = incr_2*(float)nbin_2 ;
  for (i = 0 ; i < N_BINS ; i++) {
    histo1[i].avg = 0.0 ;
    histo1[i].std = 0.0 ;
    histo1[i].freq = 0 ;
    histo2[i].avg = 0.0 ;
    histo2[i].std = 0.0 ;
    histo2[i].freq = 0 ;
    c_histo1[i].avg = 0.0 ;
    c_histo1[i].std = 0.0 ;
    c_histo1[i].freq = 0 ;
    c_histo2[i].avg = 0.0 ;
    c_histo2[i].std = 0.0 ;
    c_histo2[i].freq = 0 ;
    abs_freq = 0 ;
    ttl_freq = 0 ;
    avg_err = 0.0 ;
    *ttl_err = 0.0 ;
    *ttl_std = 0.0 ;
  }
  for (i = KERNEL_Y + NRADIUS ; i < f->sizy - KERNEL_Y - NRADIUS ; i++) {
    for (j = KERNEL_X + NRADIUS ; j < f->sizx -  KERNEL_X - NRADIUS ; j++) {
      det = (*f->param_ptr[f->sizz/2])[f->res*i + j].gauss ;
      cond = (*f->param_ptr[f->sizz/2])[f->res*i + j].cond ;
      index1 = (int)((det/max_det)*(float)nbin_1) ;
      index2 = (int)((cond/max_cond)*(float)nbin_2) ;
      u = (*f->flow_ptr)[f->res*i+j] ;
      abs_freq++ ;
      if (u.x != -100.0 && u.y != 100.0) {
        err = psi_error((*f->flow_ptr)[f->res*i+j],(*q->flow_ptr)[q->res*i+j]) ;
        avg_err += err ;
        ttl_freq++ ;
        if (index1 >= nbin_1) {
          index1 = nbin_1 - 1 ;
        }
        histo1[index1].avg += err ;
        histo1[index1].freq++ ;
        if (index2 >= nbin_2) {
          index2 = nbin_2 - 1 ;
        }
        histo2[index2].avg += err ;
        histo2[index2].freq++ ;
        for (k = 0 ; k <= index1 ; k++) {
          c_histo1[k].avg += err ;
          c_histo1[k].freq++ ;
        }
        for (k = index2 ; k < nbin_2 ; k++) {
          c_histo2[k].avg += err ;
          c_histo2[k].freq++ ;
        }
      } 
    }
  }
  if (abs_freq != 0) {
    *dens = (float)ttl_freq/(float)abs_freq ;
  }
  else {
    *dens = 0.0 ;
  }
  *ttl_err = avg_err/(float)ttl_freq ;
  avg_err = *ttl_err ;
  for (i = 0 ; i < nbin_1 ; i++) {
    if (histo1[i].freq != 0) {
      histo1[i].avg /= (float)histo1[i].freq ;
    }
    if (c_histo1[i].freq != 0) {
      c_histo1[i].avg /= (float)c_histo1[i].freq ;
    } 
  }
  for (i = 0 ; i < nbin_2 ; i++) {
    if (histo2[i].freq != 0) {
      histo2[i].avg /= (float)histo2[i].freq ;
    }
    if (c_histo2[i].freq != 0) {
      c_histo2[i].avg /= (float)c_histo2[i].freq ;
    }
  }
  for (i = KERNEL_Y + NRADIUS ; i < f->sizy - KERNEL_Y - NRADIUS ; i++) {
    for (j = KERNEL_X + NRADIUS ; j < f->sizx - KERNEL_X - NRADIUS ; j++) {
      det = (*f->param_ptr[f->sizz/2])[f->res*i + j].gauss ;
      cond = (*f->param_ptr[f->sizz/2])[f->res*i + j].cond ;
      index1 = (int)((det/max_det)*(float)nbin_1) ;
      index2 = (int)((cond/max_cond)*(float)nbin_2) ;
      u = (*f->flow_ptr)[f->res*i+j] ;
      if (u.x != -100.0 && u.y != 100.0) {
        err = psi_error((*f->flow_ptr)[f->res*i+j],(*q->flow_ptr)[q->res*i+j]) ;
        *ttl_std += pow(avg_err - err,2.0) ;
        if (index1 >= nbin_1) {
          index1 = nbin_1 - 1 ;
        }
        histo1[index1].std += pow(histo1[index1].avg - err,2.0) ;
        if (index2 >= nbin_2) {
          index2 = nbin_2 - 1 ;
        }
        histo2[index2].std += pow(histo2[index2].avg - err,2.0) ;
        for (k = 0 ; k <= index1 ; k++) {
          c_histo1[k].std += pow(c_histo1[k].avg - err,2.0) ;
        }
        for (k = index2 ; k < nbin_2 ; k++) {
          c_histo2[k].std += pow(c_histo2[k].avg - err,2.0) ;
        }
      }
    }
  }
  *ttl_std = sqrt(*ttl_std/(float)ttl_freq) ;
  for (i = 0 ; i < nbin_1 ; i++) {
    if (histo1[i].freq != 0) {
      histo1[i].std = sqrt(histo1[i].std/(float)histo1[i].freq) ;
    }
    if (c_histo1[i].freq != 0) {
      c_histo1[i].std = sqrt(c_histo1[i].std/(float)c_histo1[i].freq) ;
    } 
  }
  for (i = 0 ; i < nbin_2 ; i++) {
    if (histo2[i].freq != 0) { 
      histo2[i].std = sqrt(histo2[i].std/(float)histo2[i].freq) ;
    }
    if (c_histo2[i].freq != 0) { 
      c_histo2[i].std = sqrt(c_histo2[i].std/(float)c_histo2[i].freq) ;
    }
  }
  if ((fdp = fopen(h_fn,"w")) == NULL) {
    error(6) ;
  }

  density = c_histo1[0].freq ;
  fprintf(fdp,"%2d\n",N_HISTO) ;
  fprintf(fdp,"%3d\n",nbin_1) ;
  fprintf(fdp,"%5.7f\n",max_det - incr_1/2.0) ;
  for (i = 0 ; i < nbin_1 ; i++) {
    t = max_det/(float)nbin_1 ;
    fprintf(fdp,"%5.7f %10.7f %10.7f %3.7f\n", t*(float)(i)+t/2.0, 
    histo1[i].avg, histo1[i].std, (float)histo1[i].freq/density) ;
  }
  fprintf(fdp,"\n\n\n") ;

  fprintf(fdp,"%3d\n",nbin_1) ;
  fprintf(fdp,"%5.7f\n",max_det - incr_1/2.0) ;
  for (i = 0 ; i < nbin_1 ; i++) {
    t = max_det/(float)nbin_1 ; 
    fprintf(fdp,"%5.7f %10.7f %10.7f %3.7f\n", t*(float)(i)+t/2.0, 
    c_histo1[i].avg, c_histo1[i].std, (float)c_histo1[i].freq/density) ;
  }
  fprintf(fdp,"\n\n\n") ;

  fprintf(fdp,"%3d\n",nbin_2) ;
  fprintf(fdp,"%5.7f\n",max_cond - incr_2/2.0) ;
  for (i = 0 ; i < nbin_2 ; i++) {
    t = max_cond/(float)nbin_2 ;
    fprintf(fdp,"%5.7f %10.7f %10.7f %3.7f\n", t*(float)(i)+t/2.0, 
    histo2[i].avg, histo2[i].std, (float)histo2[i].freq/density) ;
  }
  fprintf(fdp,"\n\n\n") ;

  fprintf(fdp,"%3d\n",nbin_2) ;
  fprintf(fdp,"%5.7f\n",max_cond - incr_2/2.0) ;
  for (i = 0 ; i < nbin_2 ; i++) {
    t = max_cond/(float)nbin_2 ; 
    fprintf(fdp,"%5.7f %10.7f %10.7f %3.7f\n", t*(float)(i)+t/2.0, 
    c_histo2[i].avg, c_histo2[i].std, (float)c_histo2[i].freq/density) ;
  }
  fclose(fdp) ;
}

/* 
           NAME : load_frames(fn,start,n_frame,x,y,cub,bin)
   PARAMETER(S) : fn : filename if the image sequence;
               start : start number of frames;
             n_frame : number of frames needed;
                 x,y : image size; 
                 cub : pointer on an image cube;
                 bin : binary input files.
 
        PURPOSE : Loads the image sequence from frame start to
                  frame start + n_frame into image cube cub1.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 25 1992
*/

load_frames(fn,start,n_frame,x,y,cub,bin)
string fn ;
int start, n_frame, x, y, bin ;
qnode_ptr_t cub ;

{ unsigned char hd[H] ;
  string s1, s2 ;
  int i ;

  for (i = 0 ; i < n_frame ; i++) {
    itoa(start,s1) ;
    concat(fn,s1,s2) ;
    if (!bin) {
      pgetrast(s2,hd,cub->gauss_ptr[i],x,y,X) ;
    }
    else {
      Bpgetrast(s2,cub->gauss_ptr[i],x,y,X) ;
    }
    start++ ;
  }
}

/* 
           NAME : compute_deriv(c,f)
   PARAMETER(S) :     c : pointer on an image cube node;
                      f : pointer on a flow field and
                          its parameters.
 
        PURPOSE : Computes image derivatives needed for
                  the second-order method of Uras et al [88].

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : September 16 1990
*/

compute_deriv(c,f)
qnode_ptr_t c, f ;

{ float dxx(), dyy(), dxy(), dx(), dy(), dt(), dxt(), dyt() ;
  int i, j, k ;
param_t P ;
  for (i = f->ofst ; i < f->sizy - f->ofst ; i++) {
    for (j = f->ofst ; j < f->sizx - f->ofst ; j++) {
      for (k = 0 ; k < f->sizz ; k++) {
        (*f->param_ptr[k])[f->res*i + j].fx = dx(c,f,i,j,k) ;
        (*f->param_ptr[k])[f->res*i + j].fy = dy(c,f,i,j,k) ;
      }
    }
  }
  k = c->sizz/2.0 ; 
  for (i = f->ofst ; i < f->sizy - f->ofst ; i++) {
    for (j = f->ofst ; j < f->sizx - f->ofst ; j++) { 
      if ((*f->param_ptr[k])[f->res*i + j].err == NO_ERROR) { 
        (*f->param_ptr[k])[f->res*i + j].ft = dt(c,f,i,j) ;
        (*f->param_ptr[k])[f->res*i + j].fxt = dxt(f,i,j) ;
        (*f->param_ptr[k])[f->res*i + j].fyt = dyt(f,i,j) ;
      }
    }
  } 
  for (i = f->ofst ; i < f->sizy - f->ofst ; i++) {
    for (j = f->ofst ; j < f->sizx - f->ofst ; j++) { 
      if ((*f->param_ptr[k])[f->res*i + j].err == NO_ERROR) { 
        (*f->param_ptr[k])[f->res*i + j].fxx = dxx(f,i,j,k) ;
        (*f->param_ptr[k])[f->res*i + j].fyy = dyy(f,i,j,k) ;
        (*f->param_ptr[k])[f->res*i + j].fxy = dxy(f,i,j,k) ;
      }
    }
  }
}

/* 
           NAME : compute_discr(fl)
   PARAMETER(S) : flow field node.
 
        PURPOSE : computes ||Mt*gradient(U)||/||gradient(It)||
                  and gaussian curvature for all non-error points.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : November 7 1990
*/

compute_discr(fl)
qnode_ptr_t fl ;

{ param_t P ;
  float dux(), duy(), dvx(), dvy(), ux, uy, vx, vy, MtU, It, cmax, cmin ;
  int i, j ;

  for (i = fl->ofst ; i < fl->sizy - fl->ofst ; i++) {
    for (j = fl->ofst ; j < fl->sizx - fl->ofst ; j++) {
      P = (*fl->param_ptr[fl->sizz/2])[fl->res*i + j] ;
      if (P.err == NO_ERROR) {
        ux = dux(fl,i,j) ;
        uy = duy(fl,i,j) ;
        vx = dvx(fl,i,j) ;
        vy = dvy(fl,i,j) ;
        MtU = pow(P.fx*ux + P.fy*vx,2.0) + pow(P.fx*uy + P.fy*vy,2.0) ;
        It = pow(P.fxt,2.0) + pow(P.fyt,2.0) ;
        P.discr = sqrt(MtU/It) ;
        P.gauss = fabs(P.fxx*P.fyy - pow(P.fxy,2.0)) ;
        cmax = 0.5*((P.fxx+P.fyy)+sqrt(pow(P.fxx-P.fyy,2.0)+4.0*P.fxy*P.fxy)) ;
        cmin = 0.5*((P.fxx+P.fyy)-sqrt(pow(P.fxx-P.fyy,2.0)+4.0*P.fxy*P.fxy)) ;
        if (fabs(cmin) < fabs(cmax)) {
          if (cmin != 0.0) {
            P.cond = fabs(cmax)/fabs(cmin) ;
          }
          else {
            P.cond = MAX_COND ;
          }
        }
        else {
          if (cmax != 0.0) {
            P.cond = fabs(cmin)/fabs(cmax) ;
          }
          else {
            P.cond = MAX_COND ;
          }
        }
        (*fl->param_ptr[fl->sizz/2])[fl->res*i + j] = P ;
      }
    }
  }
}

/* 
           NAME : compute_flow(fl)
   PARAMETER(S) : fl : pointer on a flow field and its parameters.
 
        PURPOSE : Computes image displacements with
                  the second-order method of Nagel[87].

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : September 16 1990
*/

compute_flow(fl)
qnode_ptr_t fl ;

{ disp_vect_t u, v ;
  param_t P ;
  float mag() ;
  int i, j ;

  for (i = fl->ofst ; i < fl->sizy - fl->ofst ; i++) {
    for (j = fl->ofst ; j < fl->sizx - fl->ofst ; j++) { 
      P = (*fl->param_ptr[fl->sizz/2])[fl->res*i + j] ;
      if (P.err == NO_ERROR) {      
        if (P.fxx*P.fyy - pow(P.fxy,2.0) != 0.0) {
          u.x = (P.fyt*P.fxy - P.fxt*P.fyy)/(P.fxx*P.fyy - pow(P.fxy,2.0)) ;
          u.y = (P.fxt*P.fxy - P.fyt*P.fxx)/(P.fxx*P.fyy - pow(P.fxy,2.0)) ;
          if (mag(u) > MAXFLOW) { 
            v.x = (u.x/mag(u))*MAXFLOW ;
            v.y = (u.y/mag(u))*MAXFLOW ;
            u = v ;
          }
        }
        else {  
          u.x = 0.0 ;
          u.y = 0.0 ;
          (*fl->param_ptr[fl->sizz/2])[fl->res*i + j].err = NO_FLOW ;
        }
      }
      else {
        u.x = 0.0 ;
        u.y = 0.0 ; 
      }
      (*fl->flow_ptr)[fl->res*i + j] = u ;
    }
  }
}

/* 
           NAME : uras.c
   PARAMETER(S) : line parameters.

        PURPOSE : Computes optic displacement
                  between two image frames.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : September 16 1990
*/

main(argc,argv)
int argc ;
char *argv[] ;

{ qnode_ptr_t create_node(), cub1h, cub1q, cub2h, cub2q, flh, flq,
              cub1, cub2, fl ; 
  unsigned char hd[H] ;
  float sigma1, sigma2, sigma3, trsh, incr_1, incr_2, avg_err, std, density ;
  int histo, sm, sf, mid_num, start_num, x, y, z, nbin_1, nbin_2,
  binary, row, col, i ;
  string in_path, out_path, i_fname, v_fname, c_fname, h_fname ;

  valid_option(argc,argv,in_path,out_path,&sigma1,&sigma2,&sigma3,&histo,&sm,
               &sf,&trsh,i_fname,v_fname,c_fname,h_fname,&nbin_1,&incr_1,
               &nbin_2,&incr_2,&mid_num,&binary,&row,&col) ;
  printf("STAGE 1: Options Validated\n") ; fflush(stdout) ;
  generate_gauss(&ker1,sigma1) ;
  generate_gauss(&ker2,sigma2) ;
  generate_gauss(&ker3,sigma3) ;
  z = frames(ker2.m) ;
  start_num = start_number(mid_num,ker2.m) ;
  if (!binary) {
    raster_size(i_fname,start_num,&x,&y) ;
  }
  else {
    binary_size(row,col,&x,&y) ;
  }
  init_list(&cub1h,&cub1q) ;
  init_list(&cub2h,&cub2q) ;
  init_list(&flh,&flq) ;
  cub1 = create_node(0,X,x,y,z,0) ;
  cub2 = create_node(0,X,x,y,z - ker2.m + 1,0) ;
  insert_node(cub1,&cub1h,&cub1q) ;
  insert_node(cub2,&cub2h,&cub2q) ;
  init_cube(cub1,z) ;
  init_cube(cub2,z - ker2.m + 1) ;
  printf("STAGE 2: Initializations Completed\n") ; fflush(stdout) ;
  load_frames(i_fname,start_num,z,x,y,cub1,binary) ;
  printf("STAGE 3: Input Data Read\n") ; fflush(stdout) ;
  convolve(cub1,cub2,ker1,ker2,z) ;
  printf("STAGE 4: Convolutions Completed\n") ; fflush(stdout) ;
  free_cube(cub1,z) ;
  delete_node(cub1,&cub1h,&cub1q) ;
  fl = create_node(0,X,x,y,z - ker2.m + 1,cub2->ofst) ;
  insert_node(fl,&flh,&flq) ;
  init_flow(fl,z - ker2.m + 1) ;
  init_central(&C) ;
  init_beaudet(&Ix,&Ixx,&Ixy) ;
  compute_deriv(cub2,fl) ;
  printf("STAGE 5: Derivatives Computed\n") ; fflush(stdout) ;
  compute_flow(fl) ;
  printf("STAGE 6: Flow Estimated\n") ; fflush(stdout) ;
  if (sm != NO_REG) {
    if (sm == TR1) {
      KERNEL_X = (int)(3.0*sigma1 + 1.0) ;
      KERNEL_Y = (int)(3.0*sigma1 + 1.0) ;
      compute_discr(fl) ;
      regul1(fl,REGION_SIZE,sf,trsh) ;
    }
    else {
      KERNEL_X = (int)(3.0*sigma1 + 1.0) + (int)(3.0*sigma3 + 1.0) ;
      KERNEL_Y = (int)(3.0*sigma1 + 1.0) + (int)(3.0*sigma3 + 1.0) ;
      regul2(fl,ker3) ;
    }
    printf("STAGE 7: Regularization Completed\n") ; fflush(stdout) ;
  }
  if (histo) {
    prod_histo(fl,c_fname,h_fname,nbin_1,incr_1,nbin_2,incr_2,&avg_err,&std,
               &density) ;
    printf("STAGE 8: Histograms Produced\n") ; fflush(stdout) ;
  }
  write_velocity(v_fname,fl) ;
  printf("STAGE 9: Flow Written to File\n") ; fflush(stdout) ;
  free_cube(cub2,z - ker2.m + 1) ;
  free_flow(fl,z - ker2.m + 1) ;
  delete_node(cub2,&cub2h,&cub2q) ;
  delete_node(fl,&flh,&flq) ; 
  printf("STAGE 10: End of Program\n\n") ;
  if (histo) {
    printf("         Average Angular Error: %10.5f\n", avg_err) ;
    printf("            Standard Deviation: %10.5f\n", std) ;
    printf("                       Density: %10.5f\n", density*100.0) ;
    fflush(stdout) ;
  }
}
