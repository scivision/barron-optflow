#include <math.h>
#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
/* 
           NAME : const.h 
   PARAMETER(S) : none
 
        PURPOSE : Definition of the constants for the 
                  implementation of Nagel[87]'s motion 
                  detection approach.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : September 13 1990

*/

#define TRUE        1
#define FALSE       0

#define H           32
#define SMSK        5
#define X           316
#define Y           316
#define Z           20

#define DEF_S1      3.0 
#define DEF_S2      1.5
#define ALPHA       5.0
#define DELTA       1.0

#define N_HISTO     2
#define N_BINS      100

#define NRADIUS     9
#define SKIP        1

#define MAXFLOW      20.0
#define NO_ERROR     0
#define SA_OVERFLOW  1
#define NO_FLOW      2

#define RELAXITER    10

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

typedef struct t_kernel {
                          float k[X], f ;
                          int m ;
                        } kernel_t ;

typedef struct t_disp_vect {
                             float x, y ;  
                           } disp_vect_t ;

typedef struct t_param {
                         float fxx, fyy, fxy, ft, fx, fy, mag ;
                         int err ;
                       } param_t ;


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

kernel_t  ker1, ker2, C ;
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

{ float d ;

  d = (*p->flow_ptr)[p->res*(x-1) + y].x - (*p->flow_ptr)[p->res*(x+1) + y].x ;
  return(0.5*d) ;
}

/* 
           NAME : wlapu(p,x,y,a,c)
   PARAMETER(S) : p : flow node;
                x,y : flow location for computation;
                a,c : scaling factors.
 
        PURPOSE : computes a weighted version of uxx + uyy.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 4 1990
*/

float wlapu(p,x,y,a,c)
qnode_ptr_t p ;
int x, y ;
float a, c ;

{ float d ;

  d = a*((*p->flow_ptr)[p->res*(x-1) + y].x +
  (*p->flow_ptr)[p->res*(x+1) + y].x) +
  c*((*p->flow_ptr)[p->res*x + (y-1)].x +
  (*p->flow_ptr)[p->res*x + (y+1)].x) ;
   return(d) ;
}

/* 
           NAME : wlapv(p,x,y,a,c)
   PARAMETER(S) : p : flow node;
                x,y : flow location for computation;
                a,c : scaling factors.
 
        PURPOSE : computes a weighted version of vxx + vyy.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 4 1990
*/

float wlapv(p,x,y,a,c)
qnode_ptr_t p ;
int x, y ;
float a, c ;

{ float d ;

  d = a*((*p->flow_ptr)[p->res*(x-1) + y].y +
  (*p->flow_ptr)[p->res*(x+1) + y].y) +
  c*((*p->flow_ptr)[p->res*x + (y-1)].y +
  (*p->flow_ptr)[p->res*x + (y+1)].y) ;
  return(d) ;
}

/* 
           NAME : lapu(p,x,y)
   PARAMETER(S) : p : flow node;
                x,y : flow location for computation.
 
        PURPOSE : computes uxx + uyy as in Horn and Schunck's [81].

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 4 1990
*/

float lapu(p,x,y)
qnode_ptr_t p ;
int x, y ;

{ float d ;

  d = (1.0/6.0)*((*p->flow_ptr)[p->res*(x-1) + y].x +
  (*p->flow_ptr)[p->res*(x+1) + y].x +
  (*p->flow_ptr)[p->res*x + (y-1)].x +
  (*p->flow_ptr)[p->res*x + (y+1)].x) + 
  (1.0/12.0)*((*p->flow_ptr)[p->res*(x-1) + y-1].x +
  (*p->flow_ptr)[p->res*(x-1) + y+1].x +
  (*p->flow_ptr)[p->res*(x+1) + y-1].x +
  (*p->flow_ptr)[p->res*(x+1) + y+1].x) ;
  return(d) ;
}

/* 
           NAME : lapv(p,x,y)
   PARAMETER(S) : p : flow node;
                x,y : flow location for computation.
 
        PURPOSE : computes vxx + vyy as in Horn and Schunck's [81].

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 4 1990
*/

float lapv(p,x,y)
qnode_ptr_t p ;
int x, y ;

{ float d ;

  d = (1.0/6.0)*((*p->flow_ptr)[p->res*(x-1) + y].y +
  (*p->flow_ptr)[p->res*(x+1) + y].y +
  (*p->flow_ptr)[p->res*x + (y-1)].y +
  (*p->flow_ptr)[p->res*x + (y+1)].y) + 
  (1.0/12.0)*((*p->flow_ptr)[p->res*(x-1) + y-1].y +
  (*p->flow_ptr)[p->res*(x-1) + y+1].y +
  (*p->flow_ptr)[p->res*(x+1) + y-1].y +
  (*p->flow_ptr)[p->res*(x+1) + y+1].y) ;
  return(d) ;
}

/* 
           NAME : duxy(p,x,y)
   PARAMETER(S) : p : flow node;
                x,y : flow location for computation.
 
        PURPOSE : computes 1st order partial derviative with respect
                  to x y at x,y for component u of flow p.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 4 1990
*/

float duxy(p,x,y)
qnode_ptr_t p ;
int x, y ;

{ float dux(), d ;
  d = dux(p,x,y-1) - dux(p,x,y+1) ;
  return(0.5*d) ; 
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

{ float d ;
  
  d = (*p->flow_ptr)[p->res*x + (y-1)].x - (*p->flow_ptr)[p->res*x + (y+1)].x ;
  return(0.5*d) ;
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

{ float d ;

  d = (*p->flow_ptr)[p->res*(x-1) + y].y - (*p->flow_ptr)[p->res*(x+1) + y].y ;
  return(0.5*d) ;
}

/* 
           NAME : dvxy(p,x,y)
   PARAMETER(S) : p : flow node;
                x,y : flow location for computation.
 
        PURPOSE : computes 1st order partial derviative with respect
                  to x y at x,y for component v of flow p.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 4 1990
*/

float dvxy(p,x,y)
qnode_ptr_t p ;
int x, y ;

{ float dvx(), d ;

  d = dvx(p,x,y-1) - dvx(p,x,y+1) ;
  return(0.5*d) ;
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

{ float d ;

  d = (*p->flow_ptr)[p->res*x + (y-1)].y - (*p->flow_ptr)[p->res*x + (y+1)].y ;
  return(0.5*d) ;
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
  fprintf(stderr,"Usage: nagel input-path output-path seq-name. [-SG n.n] [-TG n.n] [-M n] [-I n] [-A n.n] [-D n.n] [-B cols rows] [-F n.n] [-C corr-vel-file nbins increment]\n") ;
  fprintf(stderr,"[-SG n.n]        : sigma value for spatial gaussian\n") ;
  fprintf(stderr,"                   default value is 3.0\n") ;
  fprintf(stderr,"[-TG n.n]        : sigma value for temporal gaussian\n") ; 
  fprintf(stderr,"                   default value is 1.5\n") ;
  fprintf(stderr,"[-M n]           : middle frame of image sequence\n") ;
  fprintf(stderr,"                   default value is 0\n") ;
  fprintf(stderr,"[-I n]           : number of relaxation iterations\n") ;
  fprintf(stderr,"                   default value is 10\n") ;
  fprintf(stderr,"[-A n.n]         : smoothing parameter alpha\n") ;
  fprintf(stderr,"                   default value is 5.0\n") ;
  fprintf(stderr,"[-D n.n]         : delta parameter\n") ;
  fprintf(stderr,"                   default value is 1.0\n") ;
  fprintf(stderr,"[-B cols rows]   : for binary input files without header\n") ;
  fprintf(stderr,"[-F n.n]         : filters out unreliable estimates\n") ;
  fprintf(stderr,"                   using the magnitude of local gradient\n") ;
  fprintf(stderr,"[-C corr-vel-file nbins increment]\n") ;
  fprintf(stderr,"                 : error histogram\n") ;  
  fprintf(stderr,"                   corr-vel-file : file of correct velocities\n") ;
  fprintf(stderr,"                   nbins : number of bins for histogram\n") ;
  fprintf(stderr,"                   increment : value of gap between bins\n") ;
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
     case 11: fprintf(stderr,"error %d : undefined error\n",n) ;
              break ;
     case 12: fprintf(stderr,"error %d : undefined error\n",n) ;
              break ;
     case 13: fprintf(stderr,"error %d : undefined error\n",n) ;
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
           NAME : filter(p,q,tresh)
   PARAMETER(S) : p : pointer on a flow node ;
                  q : pointer on a parameter node ;
                  tresh : confidence measure thresholding value.
 
        PURPOSE : Cancels out displacement estimates with conf.
                  measure Gaussian curvature lower than tresh.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : January 21 1991 ;
*/

filter(p,q,trsh)
qnode_ptr_t p, q ;
float trsh ;

{ param_t P ; 
  disp_vect_t Uerr ;
  int i, j ;

  Uerr.x = -100.0 ;
  Uerr.y = 100.0 ;

  for (i = 0 ; i < p->sizy ; i++) {
    for (j = 0 ; j < p->sizx ; j++) {
      P = (*q->param_ptr[q->sizz/2])[q->res*i + j] ;
      if (P.mag <= trsh || P.err != NO_ERROR) {
        (*p->flow_ptr)[p->res*i + j] = Uerr ;
      }
    }
  }
}

/* 
           NAME : screen(p,maxv)
   PARAMETER(S) : p : pointer on a flow node ;
               maxv : maximum magnitude on velocities.
 
        PURPOSE : deletes velocity estimates having magnitudes
                  greater than maxv.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : January 21 1991 ;
*/

screen(p,maxv)
qnode_ptr_t p ;
float maxv ;

{ disp_vect_t U, Uerr ;
  int i, j ;

  Uerr.x = -100.0 ;
  Uerr.y = 100.0 ;

  for (i = 0 ; i < p->sizy ; i++) {
    for (j = 0 ; j < p->sizx ; j++) {
      U = (*p->flow_ptr)[p->res*i + j] ;
      if (sqrt(pow(U.x,2.0)+pow(U.y,2.0)) > maxv) {
        (*p->flow_ptr)[p->res*i + j] = Uerr ;
      }
    }
  }
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
                       s : length of temporal mask for derivatives.
 
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

void free_cube(qnode_ptr_t p, int n)
{ 
  for (int i = 0 ; i < n ; i++) {
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
      (sqrt(2.0*M_PI)*sig) ;
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
           NAME : Bpputrast(fn,hd,bf,sx,sy,r) 
   PARAMETER(S) : fn    : filename;
                  bf    : Pointer on a 2D array containing the image;
                  sx,sy : image size to write;
                  r     : resolution of the array.
 
        PURPOSE : Writes the contents of bf into a 
                  binary file specified by fn.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 27 1989, modified March 20 1990.

*/

Bpputrast(fn,bf,sx,sy,r) 
char *fn ;
float bf[X*Y] ;
int sx, sy, r ;

{ int fd, i, j ;
  unsigned char c ;

  if ((fd = creat(fn,0600)) > 0) { ;
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
    error(9) ;
  }
  close(fd) ;
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
  return((float)(acos(v))*180.0/M_PI) ;
}

/* 
           NAME : load_frames(fn,hd,start,n_frame,x,y,cub,bin)
   PARAMETER(S) : fn : filename if the image sequence;
                  hd : image header (raster file);
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

load_frames(fn,hd,start,n_frame,x,y,cub,bin)
string fn ;
unsigned char hd[H] ;
int start, n_frame, x, y, bin ;
qnode_ptr_t cub ;

{ string s1, s2 ;
  int i ;

  for (i = 0 ; i < n_frame ; i++) {
    sprintf(s1,"%d",start);
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
           NAME : dump_frames(fn,hd,start,n_frame,x,y,cub,bin)
   PARAMETER(S) : fn : filename if the image sequence;
                  hd : image header (raster file);
               start : start number of frames;
             n_frame : number of frames needed;
                 x,y : image size; 
                 cub : pointer on an image cube;
                 bin : binary input files.
 
        PURPOSE : dumps the image sequence from frame start to
                  frame start + n_frame from image cube cub1.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 25 1992
*/

dump_frames(fn,hd,start,n_frame,x,y,cub,bin)
string fn ;
unsigned char hd[H] ;
int start, n_frame, x, y, bin ;
qnode_ptr_t cub ;

{ string s1, s2 ;
  int i ;

  for (i = 0 ; i < n_frame ; i++) {
    sprintf(s1,"%d",start);
    concat(fn,s1,s2) ;
    if (!bin) {
      pputrast(s2,hd,cub->gauss_ptr[i],x,y,X) ;
    }
    else {
      Bpputrast(s2,cub->gauss_ptr[i],x,y,X) ;
    }
    start++ ;
  }
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

  sprintf(s1,"%d",num);
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
    error(7) ;
  }
  close(fd) ;
}

/* 
           NAME : binary_size(row,col,x,y)
   PARAMETER(S) : row : number of rows in input image;
                  col : number of cols in input image;
                   x : raster width;
                   y : raster height.
 
        PURPOSE : sets the size of images when option -B used.
                  size. 

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 22 1992
*/

binary_size(row,col,x,y)
int row, col, *x, *y ;

{ *x = col ;
  *y = row ;
  if (*x > X || *y > Y) {
    error(10) ;
  }
}

/* 
           NAME : flow_error(f,q,avg,std)
   PARAMETER(S) : f,q : flow field pointers;
             avg, std : error statistics.
 
        PURPOSE : computes the average angular error and 
                  standard deviation of flow f.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 22 1992
*/

flow_error(f,q,avg,std) 
qnode_ptr_t f, q ;
float *avg, *std ;

{ disp_vect_t u ;
  float psi_error(), err ;
  int freq, i, j ;

  freq = 0 ;
  *avg = 0.0 ;
  *std = 0.0 ;
  for (i = KERNEL_Y + NRADIUS ; i < f->sizy - KERNEL_Y - NRADIUS ; i++) {
    for (j = KERNEL_X + NRADIUS ; j < f->sizx - KERNEL_X - NRADIUS ; j++) {
      u = (*f->flow_ptr)[f->res*i + j] ;
      if (u.x != -100.0 && u.y != 100.0) {
        *avg += psi_error((*f->flow_ptr)[f->res*i+j],
        (*q->flow_ptr)[q->res*i+j]) ;
        freq++ ;
      }
    }
  }
  *avg = *avg/(float)freq ;
  for (i = KERNEL_Y + NRADIUS ; i < f->sizy - KERNEL_Y - NRADIUS ; i++) {
    for (j = KERNEL_X + NRADIUS ; j < f->sizx - KERNEL_X - NRADIUS ; j++) {
      u = (*f->flow_ptr)[f->res*i + j] ;
      if (u.x != -100.0 && u.y != 100.0) {
        err = psi_error((*f->flow_ptr)[f->res*i+j],(*q->flow_ptr)[q->res*i+j]) ;
        *std += pow(*avg - err,2.0) ;
      }
    }
  }
  *std = sqrt(*std/(float)freq) ;
}

/* 
           NAME : flow_l2_norm(q,p)
   PARAMETER(S) : q,p : flow field pointers.
 
        PURPOSE : computes the average l2 norm of the difference
                  between flow q and p.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 22 1992
*/

float flow_l2_norm(q,p) 
qnode_ptr_t q, p ;

{ disp_vect_t Uold, Unew, T ;
  float sum ;
  int i, j, n ;

  n = 0 ;
  sum = 0.0 ;

  for (i = KERNEL_Y + NRADIUS ; i < p->sizy - KERNEL_Y - NRADIUS ; i++) {
    for (j = KERNEL_X + NRADIUS ; j < p->sizx - KERNEL_X - NRADIUS ; j++) {
      Uold = (*q->flow_ptr)[q->res*i + j] ;
      Unew = (*p->flow_ptr)[p->res*i + j] ;
      T.x = Uold.x - Unew.x ;
      T.y = Uold.y - Unew.y ;
      sum += sqrt(pow(T.x,2.0) + pow(T.y,2.0)) ;
      n++ ;
    }
  }
  return(sum/(float)n);
}

/* 
           NAME : load_velocity(q,fn)
   PARAMETER(S) :  q : flow field pointer ;
                  fn : correct velocity filename.
 
        PURPOSE : loads a correct velocity file fn into q.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 22 1992
*/

load_velocity(q,fn) 
qnode_ptr_t q ;
string fn ;

{ qnode_ptr_t create_node() ;
  float actual_x, actual_y, size_x, size_y, offset_x, offset_y, x, y ;
  int fdf, nbytes, i, j ;

  if ((fdf = open(fn,O_RDONLY)) <= 0) {
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
} 

/* 
           NAME : nagel_relax(p,alpha,delta,n,fn,prod_err)
   PARAMETER(S) : p : displacement field to relax and image derivatives;
              alpha : smoothing paramenter;
              delta : delta parameter ;
                  n : number of iterations;
                 fn : correct velocity filename;
           prod_err : produce error statistics for each iteration.
 
        PURPOSE : applies the smoothness constraint defined in Nagel[87].

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 4 1990
*/

nagel_relax(p,alpha,delta,n,fn,prod_err)
qnode_ptr_t p ;
float alpha, delta ;
int n, prod_err ;
string fn ;

{ qnode_ptr_t create_node(), q, r ;
  param_t f ;
  float dux(), duy(), dvx(), dvy(), duxy(), dvxy(), wlapu(), wlapv(),
        flow_l2_norm(), a, b, c, d, A, B, C, D, 
        u_bar, v_bar, ux, uy, vx, vy, uxy, vxy, Ku, Kv, avg, std ;
  int   i, j, k ;

  q = create_node(0,p->res,p->sizx,p->sizy,p->sizz,p->ofst) ;
  q->flow_ptr = (disp_field512_t *)malloc(sizeof(disp_field512_t)) ;

  if (prod_err) {
    r = create_node(0,p->res,p->sizx,p->sizy,p->sizz,p->ofst) ;
    r->flow_ptr = (disp_field512_t *)malloc(sizeof(disp_field512_t)) ;
    load_velocity(r,fn) ;
  }

  for (i = p->ofst ; i < p->sizy ; i++) {
    for (j = p->ofst ; j < p->sizx ; j++) {
      (*p->flow_ptr)[p->res*i + j].x = 0.0 ;
      (*p->flow_ptr)[p->res*i + j].y = 0.0 ;
    }
  }
  for (k = 0 ; k < n ; k++) {
    for (i = p->ofst ; i < p->sizy ; i++) {
      for (j = p->ofst ; j < p->sizx ; j++) {
        (*q->flow_ptr)[q->res*i + j] = (*p->flow_ptr)[p->res*i + j] ;
      }
    }
    for (i = p->ofst ; i < p->sizy - p->ofst ; i++) {
      for (j = p->ofst ; j < p->sizx - p->ofst ; j++) {
        if ((*p->param_ptr[p->sizz/2])[p->res*i + j].err == NO_ERROR) {
          f = (*p->param_ptr[p->sizz/2])[p->res*i + j] ;
          d = pow(f.fx,2.0) + pow(f.fy,2.0) + 2.0*delta ;
          a = (pow(f.fy,2.0) + delta)/d ;
          b = (-f.fx*f.fy)/d ;
          c = (pow(f.fx,2.0) + delta)/d ;
          A = f.fyy + 2.0*(f.fxx*a + f.fxy*b) ;
          B = -f.fxy + 2.0*(f.fxx*b + f.fxy*c) ;
          C = -f.fxy + 2.0*(f.fxy*a + f.fyy*b) ;
          D = f.fxx + 2.0*(f.fxy*b + f.fyy*c) ;

          ux = dux(q,i,j) ;
          uy = duy(q,i,j) ;
          vx = dvx(q,i,j) ;
          vy = dvy(q,i,j) ;
          uxy = duxy(q,i,j) ;
          vxy = dvxy(q,i,j) ;

          u_bar = wlapu(q,i,j,a,c) - 2.0*b*uxy ;
          v_bar = wlapv(q,i,j,a,c) - 2.0*b*vxy ;

          Ku = 0.5*(-(f.fx*(ux*A + uy*B) + f.fy*(ux*C + uy*D))/d + u_bar) ;
          Kv = 0.5*(-(f.fx*(vx*A + vy*B) + f.fy*(vx*C + vy*D))/d + v_bar) ;
      
          (*p->flow_ptr)[p->res*i + j].x = 
          Ku - f.fx*(f.fx*Ku + f.fy*Kv + f.ft)
          /(pow(f.fx,2.0) + pow(f.fy,2.0) + 2.0*pow(alpha,2.0)) ;

          (*p->flow_ptr)[p->res*i + j].y = 
          Kv - f.fy*(f.fx*Ku + f.fy*Kv + f.ft)
          /(pow(f.fx,2.0) + pow(f.fy,2.0) + 2.0*pow(alpha,2.0)) ;
        }
      }
    }
    printf(" Iteration %3d :\n",k+1) ;
    if (prod_err) {
      flow_error(p,r,&avg,&std) ;
      printf("    Average Angular Error : %10.5f\n",avg) ; 
      printf("       Standard Deviation : %10.5f\n",std) ;
    }
    printf("       L2 Norm Difference : %10.5f\n",flow_l2_norm(q,p)) ;
    fflush(stdout) ;
  }
  free((disp_field512_t *)q->flow_ptr) ;
  free(q) ;
  if (prod_err) { 
    free((disp_field512_t *)r->flow_ptr) ;
    free(r) ;
  } 
}

/* 
           NAME : horn_relax(p,alpha,n,fn,prod_err)
   PARAMETER(S) : p : displacement field to relax and image derivatives;
              alpha : smoothing paramenter;
                  n : number of iterations;
                 fn : correct velocity filename;
           prod_err : produce error statistics for each iteration.
 
        PURPOSE : applies the smoothness constraint defined in 
                  Horn and Schunck's [81].

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 25 1993
*/

horn_relax(p,alpha,n,fn,prod_err)
qnode_ptr_t p ;
float alpha ;
int n, prod_err ;
string fn ;

{ qnode_ptr_t create_node(), q, r ;
  param_t f ;
  float lapu(), lapv(), flow_l2_norm(), u_bar, v_bar, avg, std ;
  int i, j, k ;

  q = create_node(0,p->res,p->sizx,p->sizy,p->sizz,p->ofst) ;
  q->flow_ptr = (disp_field512_t *)malloc(sizeof(disp_field512_t)) ;

  if (prod_err) {
    r = create_node(0,p->res,p->sizx,p->sizy,p->sizz,p->ofst) ;
    r->flow_ptr = (disp_field512_t *)malloc(sizeof(disp_field512_t)) ;
    load_velocity(r,fn) ;
  }

  for (i = p->ofst ; i < p->sizy ; i++) {
    for (j = p->ofst ; j < p->sizx ; j++) {
      (*p->flow_ptr)[p->res*i + j].x = 0.0 ;
      (*p->flow_ptr)[p->res*i + j].y = 0.0 ;
    }
  }
  for (k = 0 ; k < n ; k++) {
    for (i = p->ofst ; i < p->sizy ; i++) {
      for (j = p->ofst ; j < p->sizx ; j++) {
        (*q->flow_ptr)[q->res*i + j] = (*p->flow_ptr)[p->res*i + j] ;
      }
    }
    for (i = p->ofst ; i < p->sizy - p->ofst ; i++) {
      for (j = p->ofst ; j < p->sizx - p->ofst ; j++) {
        if ((*p->param_ptr[p->sizz/2])[p->res*i + j].err == NO_ERROR) {
          f = (*p->param_ptr[p->sizz/2])[p->res*i + j] ;

          u_bar = lapu(q,i,j) ;
          v_bar = lapv(q,i,j) ;

          (*p->flow_ptr)[p->res*i + j].x = 
          u_bar - f.fx*(f.fx*u_bar + f.fy*v_bar + f.ft)
          /(pow(f.fx,2.0) + pow(f.fy,2.0) + pow(alpha,2.0)) ;

          (*p->flow_ptr)[p->res*i + j].y = 
          v_bar - f.fy*(f.fx*u_bar + f.fy*v_bar + f.ft)
          /(pow(f.fx,2.0) + pow(f.fy,2.0) + pow(alpha,2.0)) ;
        }
      }
    }
    printf(" Iteration %3d :\n",k+1) ;
    if (prod_err) {
      flow_error(p,r,&avg,&std) ;
      printf("    Average Angular Error : %10.5f\n",avg) ; 
      printf("       Standard Deviation : %10.5f\n",std) ;
    }
    printf("       L2 Norm Difference : %10.5f\n",flow_l2_norm(q,p)) ;
    fflush(stdout) ;
  }
  free((disp_field512_t *)q->flow_ptr) ;
  free(q) ;
  if (prod_err) { 
    free((disp_field512_t *)r->flow_ptr) ;
    free(r) ;
  } 
}

/* 
           NAME : valid_option(argc,argv,in_path,out_path,sigma1,sigma2,histo,
                  sf,trsh,i_fname,v_fname,c_fname,h_fname,nbin,incr,
                  n_frame,alpha,iter,binary,row,col,horn,vmag,maxv)

   PARAMETER(S) : argc : argument count;
                  argv : argument values;
               in_path : path name for input data;
              out_path : path name for output data;
                sigma1 : spatial sigma;
                sigma2 : temporal sigma;
                 histo : error histogram option;
                    sf : filtering option;
                  trsh : threshold value for filtering;
               i_fname : input filename;
               v_fname : velocity filename;
               c_fname : correct velocities file name;
               h_fname : hisogram data file name ;
                  nbin : number of bins in histo 1;
                  incr : increment in histo 1 ;
               n_frame : number of required frames;
                 alpha : smoothing parameter (r;elaxation);
                 delta : delta parameter ;
                  iter : number iterations for relaxation;
                binary : input files without header;
                   row : number of rows in input files;
                   col : number of cols in input files.
                  horn : Horn and Schunck regularization option;
                  vmag : removes outliers;
                  maxv : maximum velocity magnitudes tolerated.
 
        PURPOSE : Validates line arguments.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : June 25 1990
*/

valid_option(argc,argv,in_path,out_path,sigma1,sigma2,histo,sf,trsh,
             i_fname,v_fname,c_fname,h_fname,nbin,incr,n_frame,alpha,
             delta,iter,binary,row,col,horn,vmag,maxv) 
int argc ;
char *argv[] ;
int *nbin, *iter, *n_frame, *histo, *sf, *binary, *row, *col, *horn, *vmag ;
float *incr, *sigma1, *sigma2, *trsh, *alpha, *delta, *maxv ;
string in_path, out_path, i_fname, v_fname, c_fname, h_fname ;

{ int i ;
  string ext, t0, str_sg, str_tg, str_trsh, str_alpha, str_delta ;

  if (argc == 1) {
    usage() ;
  }
  if ((argc <= 28) && (argc >= 4)) {
    *binary = FALSE ;
    *horn = FALSE ;
    *vmag = FALSE ;
    *sf = FALSE ;
    *histo = FALSE ;
    *trsh = 0.0 ;
    *maxv = 0.0 ;
    *sigma1 = DEF_S1 ;
    strcpy(str_sg,"3.0") ;
    *sigma2 = DEF_S2 ;
    strcpy(str_tg,"1.5") ;
    *iter = RELAXITER ;
    *alpha = ALPHA ;
    strcpy(str_alpha,"5.0") ;
    *delta = DELTA ;
    strcpy(str_delta,"1.0") ;
    *n_frame = 0 ;
    *nbin = 0 ;
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
    if (strcmp("-A",argv[i]) == 0) {
      if (i + 1 < argc) {
        sscanf(argv[i+1],"%f",alpha) ;
        strcpy(str_alpha,argv[i+1]) ;
        i += 2 ;
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
          if (i + 3 < argc) {
            strcpy(c_fname,argv[i+1]) ;
            sscanf(argv[i+2],"%d",nbin) ;
            sscanf(argv[i+3],"%f",incr) ;
            *histo = TRUE ;
            i += 4 ;
          }
          else {
            error(1) ;
          }
        } 
        else {
          if (strcmp("-I",argv[i]) == 0) {
            if (i + 1 < argc) {
              sscanf(argv[i+1],"%d",iter) ;
              i += 2 ;
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
                    if (strcmp("-H",argv[i]) == 0) {
                      *horn = TRUE ;
                      i++ ;
                    }
                    else {
                      if (strcmp("-D",argv[i]) == 0) {
                        if (i + 1 < argc) {
                          sscanf(argv[i+1],"%f",delta) ;
                          strcpy(str_delta,argv[i+1]) ;
                          i += 2 ;
                        }
                        else {
                          error(1) ;
                        }
                      }
                      else {
                        if (strcmp("-V",argv[i]) == 0) {
                          if (i + 1 < argc) {
                            *vmag = TRUE ;
                            sscanf(argv[i+1],"%f",maxv) ;
                            i += 2 ;
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
        }
      }
    }
  }
  if (*nbin > N_BINS) {
    error(15) ;
  }
  concat(argv[1],"/",in_path) ;
  concat(argv[2],"/",out_path) ;
  concat(in_path,argv[3],i_fname) ;
  strcpy(ext,"-sg") ;
  concat(ext,str_sg,ext) ;
  concat(ext,"-tg",ext) ;
  concat(ext,str_tg,ext) ;
  concat(ext,"-m",ext) ;
  sprintf(t0,"%d",*n_frame);
  concat(ext,t0,ext) ;
  concat(ext,"-i",ext) ;
  sprintf(t0,"%d",*iter);
  concat(ext,t0,ext) ;
  concat(ext,"-a",ext) ;
  concat(ext,str_alpha,ext) ;
  if (!(*horn)) {
    concat(ext,"-d",ext) ;
    concat(ext,str_delta,ext) ;
  }
  if (*horn) {
    concat(ext,"-h",ext) ;
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
           NAME : prod_histo(f,c_fn,h_fn,nbin,incr,ttl_err,ttl_std,dens)
   PARAMETER(S) :    f : flow field with parameters;
                  c_fn : correct velocities file name;
                  h_fn : histogram file name;
                  nbin : number of bins in histo 1;
                  incr : increment per bin in histo 1;
               ttl_err : error in flow field;
               ttl_std : standard deviation;
                  dens : denstiy of flow field.
 
        PURPOSE : produces an error histogram
                  h1: error vs confidence measure Gaussian curvature;
                  h2: h1 cumulated.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : February 20 1992
*/

prod_histo(f,c_fn,h_fn,nbin,incr,ttl_err,ttl_std,dens)
qnode_ptr_t f ;
string c_fn, h_fn ;
float incr, *ttl_err, *ttl_std, *dens ;
int nbin ;

{ qnode_ptr_t create_node(), q ;
  disp_vect_t u ;
  histo_t histo1[N_BINS], c_histo1[N_BINS] ;
  float psi_error(), max_mag, mag, err, avg_err, density, t, x, y ;
  int fdf, nbytes, i, j, k, index1, ttl_freq, abs_freq ;
  FILE *fdp ;
  extern int KERNEL_X, KERNEL_Y ;

  q = create_node(0,f->res,f->sizx,f->sizy,f->sizz,0) ;
  q->flow_ptr = (disp_field512_t *)malloc(sizeof(disp_field512_t)) ;
  load_velocity(q,c_fn) ;

  max_mag = incr*(float)nbin ;
  for (i = 0 ; i < N_BINS ; i++) {
    histo1[i].avg = 0.0 ;
    histo1[i].std = 0.0 ;
    histo1[i].freq = 0 ;
    c_histo1[i].avg = 0.0 ;
    c_histo1[i].std = 0.0 ;
    c_histo1[i].freq = 0 ;
  }
  abs_freq = 0 ;
  ttl_freq = 0 ;
  avg_err = 0.0 ;
  *ttl_err = 0.0 ;
  *ttl_std = 0.0 ;

  for (i = KERNEL_Y + NRADIUS ; i < f->sizy - KERNEL_Y - NRADIUS ; i++) {
    for (j = KERNEL_X + NRADIUS ; j < f->sizx -  KERNEL_X - NRADIUS ; j++) {
      mag = (*f->param_ptr[f->sizz/2])[f->res*i + j].mag ;
      index1 = (int)((mag/max_mag)*(float)nbin) ;
      u = (*f->flow_ptr)[f->res*i+j] ;
      abs_freq++ ;
      if (u.x != -100.0 && u.y != 100.0) {
        err = psi_error((*f->flow_ptr)[f->res*i+j],(*q->flow_ptr)[q->res*i+j]) ;
        avg_err += err ;
        ttl_freq++ ;
        if (index1 >= nbin) {
          index1 = nbin - 1 ;
        }
        histo1[index1].avg += err ;
        histo1[index1].freq++ ;
        for (k = 0 ; k <= index1 ; k++) {
          c_histo1[k].avg += err ;
          c_histo1[k].freq++ ;
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
  for (i = 0 ; i < nbin ; i++) {
    if (histo1[i].freq != 0) {
      histo1[i].avg /= (float)histo1[i].freq ;
    }
    if (c_histo1[i].freq != 0) {
      c_histo1[i].avg /= (float)c_histo1[i].freq ;
    } 
  }
  for (i = KERNEL_Y + NRADIUS ; i < f->sizy - KERNEL_Y - NRADIUS ; i++) {
    for (j = KERNEL_X + NRADIUS ; j < f->sizx - KERNEL_X - NRADIUS ; j++) {
      mag = (*f->param_ptr[f->sizz/2])[f->res*i + j].mag ;
      index1 = (int)((mag/max_mag)*(float)nbin) ;
      u = (*f->flow_ptr)[f->res*i+j] ;
      if (u.x != -100.0 && u.y != 100.0) {
        err = psi_error((*f->flow_ptr)[f->res*i+j],(*q->flow_ptr)[q->res*i+j]) ;
        *ttl_std += pow(avg_err - err,2.0) ;
        if (index1 >= nbin) {
          index1 = nbin - 1 ;
        }
        histo1[index1].std += pow(histo1[index1].avg - err,2.0) ;
        for (k = 0 ; k <= index1 ; k++) {
          c_histo1[k].std += pow(c_histo1[k].avg - err,2.0) ;
        }
      }
    }
  }
  *ttl_std = sqrt(*ttl_std/(float)ttl_freq) ;
  for (i = 0 ; i < nbin ; i++) {
    if (histo1[i].freq != 0) {
      histo1[i].std = sqrt(histo1[i].std/(float)histo1[i].freq) ;
    }
    if (c_histo1[i].freq != 0) {
      c_histo1[i].std = sqrt(c_histo1[i].std/(float)c_histo1[i].freq) ;
    } 
  }
  if ((fdp = fopen(h_fn,"w")) == NULL) {
    error(7) ;
  }

  density = c_histo1[0].freq ;
  fprintf(fdp,"%2d\n",N_HISTO) ;
  fprintf(fdp,"%3d\n",nbin) ;
  fprintf(fdp,"%5.7f\n",max_mag - incr/2.0) ;
  for (i = 0 ; i < nbin ; i++) {
    t = max_mag/(float)nbin ;
    fprintf(fdp,"%5.7f %10.7f %10.7f %3.7f\n", t*(float)(i)+t/2.0, 
    histo1[i].avg, histo1[i].std, (float)histo1[i].freq/density) ;
  }
  fprintf(fdp,"\n\n\n") ;

  fprintf(fdp,"%3d\n",nbin) ;
  fprintf(fdp,"%5.7f\n",max_mag - incr/2.0) ;
  for (i = 0 ; i < nbin ; i++) {
    t = max_mag/(float)nbin ; 
    fprintf(fdp,"%5.7f %10.7f %10.7f %3.7f\n", t*(float)(i)+t/2.0, 
    c_histo1[i].avg, c_histo1[i].std, (float)c_histo1[i].freq/density) ;
  }
  fclose(fdp) ;
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
           NAME : compute_deriv(c,f,horn)
   PARAMETER(S) :     c : pointer on an image cube node;
                      f : pointer on a flow field and
                          its parameters;
                   horn : Horn and Schunck [81] regularization option.
 
        PURPOSE : Computes image derivatives needed for
                  the second-order method of Nagel [87].

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : September 16 1990
*/

compute_deriv(c,f,horn)
qnode_ptr_t c, f ;
int horn ;

{ float dxx(), dyy(), dxy(), dx(), dy(), dt() ;
  param_t P ;
  int i, j, k ;

  for (i = f->ofst ; i < f->sizy - f->ofst ; i++) {
    for (j = f->ofst ; j < f->sizx - f->ofst ; j++) {
      for (k = 0 ; k < f->sizz ; k++) {
        (*f->param_ptr[k])[f->res*i + j].fx = dx(c,f,i,j,k) ;
        (*f->param_ptr[k])[f->res*i + j].fy = dy(c,f,i,j,k) ;
      }
    }
  }
  k = c->sizz/2 ; 
  for (i = f->ofst ; i < f->sizy - f->ofst ; i++) {
    for (j = f->ofst ; j < f->sizx - f->ofst ; j++) { 
      if ((*f->param_ptr[k])[f->res*i + j].err == NO_ERROR) { 
        (*f->param_ptr[k])[f->res*i + j].ft = dt(c,f,i,j) ;
      }
    }
  }
  if (!horn) { 
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
  for (i = f->ofst ; i < f->sizy - f->ofst ; i++) {
    for (j = f->ofst ; j < f->sizx - f->ofst ; j++) {
      if ((*f->param_ptr[k])[f->res*i + j].err == NO_ERROR) { 
        P = (*f->param_ptr[k])[f->res*i + j] ;
        (*f->param_ptr[k])[f->res*i + j].mag =
        sqrt(pow(P.fx,2.0) + pow(P.fy,2.0)) ;
      }
    }
  }
}

/* 
           NAME : nagel.c
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
  float sigma1, sigma2, trsh, incr, alpha, delta, avg_err, std, density,
        maxv ;
  int histo, sf, start_num, mid_num, x, y, z, nbin, iter, binary, horn,
      row, col, vmag ;
  string in_path,out_path,i_fname, v_fname, c_fname, h_fname ;

  valid_option(argc,argv,in_path,out_path,&sigma1,&sigma2,&histo,&sf,&trsh,
  i_fname,v_fname,c_fname,h_fname,&nbin,&incr,&mid_num,&alpha,&delta,&iter,
  &binary,&row,&col,&horn,&vmag,&maxv) ;
  printf("STAGE 1: Options Validated\n") ; fflush(stdout) ;
  generate_gauss(&ker1,sigma1) ;
  generate_gauss(&ker2,sigma2) ;
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
  load_frames(i_fname,hd,start_num,z,x,y,cub1,binary) ;
  printf("STAGE 3: Input Data Read\n") ; fflush(stdout) ;
  convolve(cub1,cub2,ker1,ker2,z) ;
  printf("STAGE 4: Convolutions Completed\n") ; fflush(stdout) ;
  free_cube(cub1,z) ;
  delete_node(cub1,&cub1h,&cub1q) ;
  fl = create_node(0,X,x,y,z - ker2.m + 1,cub2->ofst) ;
  insert_node(fl,&flh,&flq) ;
  init_flow(fl,z - ker2.m + 1) ;
  init_central(&C) ;
  compute_deriv(cub2,fl,horn) ;
  printf("STAGE 5: Derivatives Computed\n") ; fflush(stdout) ;
  KERNEL_X = (int)(3.0*sigma1 + 1.0) ;
  KERNEL_Y = (int)(3.0*sigma1 + 1.0) ;
  if (horn) {
    horn_relax(fl,alpha,iter,c_fname,histo) ;
  }
  else {
    nagel_relax(fl,alpha,delta,iter,c_fname,histo) ;
  }
  printf("STAGE 6: Regularization Completed\n") ; fflush(stdout) ;
  if (sf) {
    filter(fl,fl,trsh) ;
    printf("STAGE 7: Filtering Out Completed\n") ; fflush(stdout) ;
  }
  if (histo) {
    prod_histo(fl,c_fname,h_fname,nbin,incr,&avg_err,&std,&density) ;
    printf("STAGE 8: Histograms Produced\n") ; fflush(stdout) ;
  }
  if (vmag) {
    screen(fl,maxv) ;
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
