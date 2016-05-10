#include <string.h>
#include <math.h>
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>

/* 
           NAME : const.h 
   PARAMETER(S) : none
 
        PURPOSE : Definition of the constants for the 
                  implementation of Anandan's hierarchical 
                  motion detection approach.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 13 1990

*/

#define TRUE      1
#define FALSE     0

#define K1        150.0
#define K2        1.0
#define K3        0.0

#define CMIN      0
#define CMAX      1

#define S         7
#define H         32
#define N_FRAME   2
#define SMALL     0.001

#define MAXGRAY   255
#define RELAXITER 10
#define MAXLEVEL  4
#define DEFLEVEL  1
#define GAUSS     1 
#define LAP       2
#define MAXSSD    2340901.0
#define SSDVAL    S*S

#define N1        512
#define N2        256
#define N3        128
#define N4        64
#define N5        32
#define N6        16
#define N7        8 

#define KERNEL_X  8
#define KERNEL_Y  8
#define NRADIUS   4
#define SKIP      1

#define N_HISTO    2
#define N_BINS     100


void error(int);

/* 
           NAME : type.h 
   PARAMETER(S) : none
 
        PURPOSE : Type definitions for pyramids
                  and related data structures.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 13 1990

*/

typedef struct t_raster {
                          int m, width, height, depth, length,
                              type, maptype, maplength ;
                        } raster_t ;

typedef struct t_kernel {
                          float k[N1], f ;
                          int m ;
                        } kernel_t ;

typedef struct t_beaudet {
                           float g[S][S] ;
                           float f ;
                           int m ;
                         } beaudet_t ;

typedef struct t_disp_vect {
                             float x, y ;  
                           } disp_vect_t ;

typedef struct t_param {
                         float th, cmin, cmax ;
                       } param_t ;

typedef struct t_search_area {
                               float k[S][S] ;
                             } search_area_t ;

typedef disp_vect_t disp_field512_t[N1*N1] ;
typedef disp_vect_t disp_field256_t[N2*N2] ;
typedef disp_vect_t disp_field128_t[N3*N3] ;
typedef disp_vect_t disp_field64_t[N4*N4] ;
typedef disp_vect_t disp_field32_t[N5*N5] ;
typedef disp_vect_t disp_field16_t[N6*N6] ;
typedef disp_vect_t disp_field8_t[N7*N7] ;

typedef param_t param512_t[N1*N1] ;
typedef param_t param256_t[N2*N2] ;
typedef param_t param128_t[N3*N3] ;
typedef param_t param64_t[N4*N4] ;
typedef param_t param32_t[N5*N5] ;
typedef param_t param16_t[N6*N6] ;
typedef param_t param8_t[N7*N7] ;

typedef float image512_t[N1*N1] ;
typedef float image256_t[N2*N2] ;
typedef float image128_t[N3*N3] ;
typedef float image64_t[N4*N4] ;
typedef float image32_t[N5*N5] ;
typedef float image16_t[N6*N6] ;
typedef float image8_t[N7*N7] ; 

typedef struct t_qnode {
                         int             res, sizx, sizy, level ;
                         image512_t      *gauss_ptr, *lap_ptr ; 
                         param512_t      *param_ptr ;
                         disp_field512_t *flow_ptr ;
                         struct t_qnode  *forth, *back ;
                       } qnode_t, *qnode_ptr_t ;

typedef struct t_histo {
                         float std, avg ;
                         int freq ;
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
           DATE : April 10 1990
*/

beaudet_t Ix, Ixx, Ixy ;
kernel_t kera, kerb ;

/* 
           NAME : alloc_flow(p)
   PARAMETER(S) : p : pointer on a node.
 
        PURPOSE : Allocates a vector field of size p->res to p->flow_ptr.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 9 1990

*/

void alloc_flow(p) 
qnode_ptr_t p ;

{
  switch(p->res) {
    case N1 : p->flow_ptr = (disp_field512_t *)malloc(sizeof(disp_field512_t)) ;
              break ;
    case N2 : p->flow_ptr = (disp_field512_t *)malloc(sizeof(disp_field256_t)) ;
              break ;
    case N3 : p->flow_ptr = (disp_field512_t *)malloc(sizeof(disp_field128_t)) ;
              break ;
    case N4 : p->flow_ptr = (disp_field512_t *)malloc(sizeof(disp_field64_t)) ;
              break ;
    case N5 : p->flow_ptr = (disp_field512_t *)malloc(sizeof(disp_field32_t)) ;
              break ;
    case N6 : p->flow_ptr = (disp_field512_t *)malloc(sizeof(disp_field16_t)) ;
              break ;
    case N7 : p->flow_ptr = (disp_field512_t *)malloc(sizeof(disp_field8_t)) ;
              break ;
    default : error(11) ;
              break ;
  }
}

/* 
           NAME : alloc_gauss(p)
   PARAMETER(S) : p : pointer on a node.
 
        PURPOSE : Allocates an image of size p->res by p->res
                  and attaches it to p->gauss_ptr.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 13 1990

*/

void alloc_gauss(p)
qnode_ptr_t p ;

{
  switch(p->res) {
    case N1 : p->gauss_ptr = (image512_t *)malloc(sizeof(image512_t)) ;
              break ; 
    case N2 : p->gauss_ptr = (image512_t *)malloc(sizeof(image256_t)) ;
              break ; 
    case N3 : p->gauss_ptr = (image512_t *)malloc(sizeof(image128_t)) ;
              break ;
    case N4 : p->gauss_ptr = (image512_t *)malloc(sizeof(image64_t)) ;
              break ;
    case N5 : p->gauss_ptr = (image512_t *)malloc(sizeof(image32_t)) ;
              break ;
    case N6 : p->gauss_ptr = (image512_t *)malloc(sizeof(image16_t)) ;
              break ;
    case N7 : p->gauss_ptr = (image512_t *)malloc(sizeof(image8_t)) ;
              break ;
    default:  error(11) ;
              break ;
  }     
} 

/* 
           NAME : alloc_image(p)
   PARAMETER(S) : p : pointer on a node.
 
        PURPOSE : Allocates two images of size p->res by p->res
                  and attaches it to p->gauss_ptr and
                  p->lap_ptr.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 13 1990

*/

void alloc_image(p)
qnode_ptr_t p ;

{
  switch(p->res) {
    case N1 : p->gauss_ptr = (image512_t *)malloc(sizeof(image512_t)) ;
              p->lap_ptr = (image512_t *)malloc(sizeof(image512_t)) ;
              break ; 
    case N2 : p->gauss_ptr = (image512_t *)malloc(sizeof(image256_t)) ;
              p->lap_ptr = (image512_t *)malloc(sizeof(image256_t)) ;
              break ; 
    case N3 : p->gauss_ptr = (image512_t *)malloc(sizeof(image128_t)) ;
              p->lap_ptr = (image512_t *)malloc(sizeof(image128_t)) ;
              break ;
    case N4 : p->gauss_ptr = (image512_t *)malloc(sizeof(image64_t)) ;
              p->lap_ptr = (image512_t *)malloc(sizeof(image64_t)) ;
              break ;
    case N5 : p->gauss_ptr = (image512_t *)malloc(sizeof(image32_t)) ;
              p->lap_ptr = (image512_t *)malloc(sizeof(image32_t)) ;
              break ;
    case N6 : p->gauss_ptr = (image512_t *)malloc(sizeof(image16_t)) ;
              p->lap_ptr = (image512_t *)malloc(sizeof(image16_t)) ;
              break ;
    case N7 : p->gauss_ptr = (image512_t *)malloc(sizeof(image8_t)) ;
              p->lap_ptr = (image512_t *)malloc(sizeof(image8_t)) ;
              break ;
    default:  error(11) ;
              break ;
  }     
} 


/* 
           NAME : alloc_lap(p)
   PARAMETER(S) : p : pointer on a node.
 
        PURPOSE : Allocates an image of size p->res by p->res
                  and attaches it to p->lap_ptr.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 13 1990

*/

void alloc_lap(p)
qnode_ptr_t p ;

{
  switch(p->res) {
    case N1 : p->lap_ptr = (image512_t *)malloc(sizeof(image512_t)) ;
              break ; 
    case N2 : p->lap_ptr = (image512_t *)malloc(sizeof(image256_t)) ;
              break ; 
    case N3 : p->lap_ptr = (image512_t *)malloc(sizeof(image128_t)) ;
              break ;
    case N4 : p->lap_ptr = (image512_t *)malloc(sizeof(image64_t)) ;
              break ;
    case N5 : p->lap_ptr = (image512_t *)malloc(sizeof(image32_t)) ;
              break ;
    case N6 : p->lap_ptr = (image512_t *)malloc(sizeof(image16_t)) ;
              break ;
    case N7 : p->lap_ptr = (image512_t *)malloc(sizeof(image8_t)) ;
              break ;
    default:  error(11) ;
              break ;
  }     
} 

/* 
           NAME : alloc_param(p)
   PARAMETER(S) : p : pointer on a node.
 
        PURPOSE : Allocates a parameter field of size p->res to p->param_ptr.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : May 2 1990

*/

void alloc_param(p) 
qnode_ptr_t p ;

{
  switch(p->res) {
    case N1 : p->param_ptr = (param512_t *)malloc(sizeof(param512_t)) ;
              break ;
    case N2 : p->param_ptr = (param512_t *)malloc(sizeof(param256_t)) ;
              break ;
    case N3 : p->param_ptr = (param512_t *)malloc(sizeof(param128_t)) ;
              break ;
    case N4 : p->param_ptr = (param512_t *)malloc(sizeof(param64_t)) ;
              break ;
    case N5 : p->param_ptr = (param512_t *)malloc(sizeof(param32_t)) ;
              break ;
    case N6 : p->param_ptr = (param512_t *)malloc(sizeof(param16_t)) ;
              break ;
    case N7 : p->param_ptr = (param512_t *)malloc(sizeof(param8_t)) ;
              break ;
    default : error(11) ;
              break ;
  }
}

/* 
           NAME : concat(s1,s2,s3) ;
   PARAMETER(S) : s1, s2 : strings to concat;
                      s3 : output string.
 
        PURPOSE : Concats s1 and s2 into s3.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : June 25 1990
*/

void concat(s1,s2,s3)
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
           NAME : qnode_ptr_t create_node(l,r,sx,sy) 
   PARAMETER(S) : l      : level in the pyramid;
                  r      : resolution at level l;
                  sx, sy : size of the image at level l.
 
        PURPOSE : Creation of a node with 
                  level l.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 17 1990

*/

qnode_ptr_t create_node(l,r,sx,sy)
int l, r, sx, sy ;

{ qnode_ptr_t p ; 
  
  p = (qnode_ptr_t)malloc(sizeof(qnode_t)) ;
  p->res = r ;
  p->sizx = sx ;
  p->sizy = sy ;
  p->level = l ;
  return(p) ;
}

/* 
           NAME : dealloc_flow(p)
   PARAMETER(S) : p : pointer on a node.
 
        PURPOSE : deallocates the vector field attached to p.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 9 1990

*/

void dealloc_flow(p) 
qnode_ptr_t p ;

{
  switch(p->res) {
    case N1 : free((disp_field512_t *)p->flow_ptr) ;
              break ;
    case N2 : free((disp_field256_t *)p->flow_ptr) ;
              break ;
    case N3 : free((disp_field128_t *)p->flow_ptr) ;
              break ;
    case N4 : free((disp_field64_t *)p->flow_ptr) ;
              break ;
    case N5 : free((disp_field32_t *)p->flow_ptr) ;
              break ;
    case N6 : free((disp_field16_t *)p->flow_ptr) ;
              break ;
    case N7 : free((disp_field8_t *)p->flow_ptr) ;
              break ;
    default : error(11) ;
              break ;
  }
}

/*
           NAME : dealloc_gauss(p)
   PARAMETER(S) : p : pointer on a node.
 
        PURPOSE : Deallocates the Gaussian image attached to p.
                  
         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 30 1990

*/

void dealloc_gauss(p)
qnode_ptr_t p ;

{ 
  switch(p->res) {
    case N1 : free((image512_t *)p->gauss_ptr) ;  
              break ;
    case N2 : free((image256_t *)p->gauss_ptr) ;  
              break ;
    case N3 : free((image128_t *)p->gauss_ptr) ;   
              break ;
    case N4 : free((image64_t *)p->gauss_ptr) ;
              break ;
    case N5 : free((image32_t *)p->gauss_ptr) ; 
              break ;
    case N6 : free((image16_t *)p->gauss_ptr) ;
              break ;
    case N7 : free((image8_t *)p->gauss_ptr) ;
              break ;
    default : error(11) ;
              break ;
  }
}

/* 
           NAME : dealloc_image(p)
   PARAMETER(S) : p : pointer on a node.
 
        PURPOSE : Deallocates the images attached to p.
                  
         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 13 1990

*/

void dealloc_image(p)
qnode_ptr_t p ;

{ 
  switch(p->res) {
    case N1 : free((image512_t *)p->gauss_ptr) ;
              free((image512_t *)p->lap_ptr) ; 
              break ;
    case N2 : free((image256_t *)p->gauss_ptr) ;
              free((image256_t *)p->lap_ptr) ; 
              break ;
    case N3 : free((image128_t *)p->gauss_ptr) ; 
              free((image128_t *)p->lap_ptr) ;
              break ;
    case N4 : free((image64_t *)p->gauss_ptr) ;
              free((image64_t *)p->lap_ptr) ;
              break ;
    case N5 : free((image32_t *)p->gauss_ptr) ; 
              free((image32_t *)p->lap_ptr) ;
              break ;
    case N6 : free((image16_t *)p->gauss_ptr) ;
              free((image16_t *)p->lap_ptr) ;
              break ;
    case N7 : free((image8_t *)p->gauss_ptr) ;
              free((image8_t *)p->lap_ptr) ;
              break ;
    default : error(11) ;
              break ;
  }
}

/*
           NAME : dealloc_lap(p)
   PARAMETER(S) : p : pointer on a node.
 
        PURPOSE : Deallocates the Laplacian image attached to p.
                  
         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 30 1990

*/

void dealloc_lap(p)
qnode_ptr_t p ;

{ 
  switch(p->res) {
    case N1 : free((image512_t *)p->lap_ptr) ;  
              break ;
    case N2 : free((image256_t *)p->lap_ptr) ;  
              break ;
    case N3 : free((image128_t *)p->lap_ptr) ;   
              break ;
    case N4 : free((image64_t *)p->lap_ptr) ;
              break ;
    case N5 : free((image32_t *)p->lap_ptr) ; 
              break ;
    case N6 : free((image16_t *)p->lap_ptr) ;
              break ;
    case N7 : free((image8_t *)p->lap_ptr) ;
              break ;
    default : error(11) ;
              break ;
  }
}

/* 
           NAME : dealloc_param(p)
   PARAMETER(S) : p : pointer on a node.
 
        PURPOSE : Deallocates the parameter field of size p->res 
                  at p->param_ptr.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : May 2 1990

*/

void dealloc_param(p) 
qnode_ptr_t p ;

{
  switch(p->res) {
    case N1 : free((param512_t *)p->param_ptr) ;
              break ;
    case N2 : free((param256_t *)p->param_ptr) ;
              break ;
    case N3 : free((param128_t *)p->param_ptr) ;
              break ;
    case N4 : free((param64_t *)p->param_ptr) ;
              break ;
    case N5 : free((param32_t *)p->param_ptr) ;
              break ;
    case N6 : free((param16_t *)p->param_ptr) ;
              break ;
    case N7 : free((param8_t *)p->param_ptr) ;
              break ;
    default : error(11) ;
              break ;
  }
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

void delete_node(p,h,q)
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
           NAME : dx(s)
   PARAMETER(S) : s : ssd surface.
 
        PURPOSE : Computes 1st order partial derivative with
                  respect to x, locally located at the center
                  of the ssd surface.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 10 1990
*/

float dx(s) 
search_area_t *s ;

{ extern beaudet_t Ix ;
  float d ;
  int i, j ;

  d = 0.0 ;
  for(i = 0 ; i < Ix.m ; i++) {
    for(j = 0 ; j < Ix.m ; j++) {
      d += Ix.g[i][j]*(*s).k[i][j] ;
    }
  }
  return(d/Ix.f) ;
}

/* 
           NAME : dxx(s)
   PARAMETER(S) : s : ssd surface.
 
        PURPOSE : Computes 2nd order partial derivative with
                  respect to x, locally located at the center
                  of the ssd surface.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 10 1990
*/

float dxx(s) 
search_area_t *s ;

{ extern beaudet_t Ixx ;
  float d ;
  int i, j ;

  d = 0.0 ;
  for(i = 0 ; i < Ixx.m ; i++) {
    for(j = 0 ; j < Ixx.m ; j++) {
      d += Ixx.g[i][j]*(*s).k[i][j] ;
    }
  }
  return(d/Ixx.f) ;
}

/* 
           NAME : dxy(s)
   PARAMETER(S) : s : ssd surface.
 
        PURPOSE : Computes 2nd order partial derivative with
                  respect to x and y, locally located at the center
                  of the ssd surface.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 10 1990
*/

float dxy(s) 
search_area_t *s ;

{ extern beaudet_t Ixy ;
  float d ;
  int i, j ;

  d = 0.0 ;
  for(i = 0 ; i < Ixy.m ; i++) {
    for(j = 0 ; j < Ixy.m ; j++) {
      d += Ixy.g[i][j]*(*s).k[i][j] ;
    }
  }
  return(d/Ixy.f) ;
}

/* 
           NAME : dy(s)
   PARAMETER(S) : s : ssd surface.
 
        PURPOSE : Computes 1st order partial derivative with
                  respect to y, locally located at the center
                  of the ssd surface (uses the transpose of Ix).

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 10 1990
*/

float dy(s) 
search_area_t *s ;

{ extern beaudet_t Ix ;
  float d ;
  int i, j ;

  d = 0.0 ;
  for(i = 0 ; i < Ix.m ; i++) {
    for(j= 0 ; j < Ix.m ; j++) {
      d += Ix.g[j][i]*(*s).k[i][j] ;
    }
  }
  return(d/Ix.f) ;
}

/* 
           NAME : dyy(s)
   PARAMETER(S) : s : ssd surface.
 
        PURPOSE : Computes 2nd order partial derivative with
                  respect to y, locally located at the center
                  of the ssd surface (uses the transpose of Ixx).

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 10 1990
*/

float dyy(s) 
search_area_t *s ;

{ extern beaudet_t Ixx ;
  float d ;
  int i, j ;

  d = 0.0 ;
  for(i = 0 ; i < Ixx.m ; i++) {
    for(j = 0 ; j < Ixx.m ; j++) {
      d += Ixx.g[j][i]*(*s).k[i][j] ;
    }
  }
  return(d/Ixx.f) ;
}

/* 
           NAME : filter(p,q,f_type,tresh)
   PARAMETER(S) : p : pointer on a flow node ;
                  q : pointer on a parameter node ;
                  f_type : type of filtering (on Cmax or Cmin) ;
                  tresh : confidence measure thresholding value.
 
        PURPOSE : Cancels out displacement estimates with conf.
                  measure magnitude lower than tresh.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : January 21 1991 ;
*/

void filter(p,q,f_type,trsh)
qnode_ptr_t p, q ;
int f_type ;
float trsh ;

{ param_t P ; 
  disp_vect_t Uerr ;
  int i, j ;

  Uerr.x = -100.0 ;
  Uerr.y = 100.0 ;

  for (i = 0 ; i < p->sizy ; i++) {
    for (j = 0 ; j < p->sizx ; j++) {
      P = (*q->param_ptr)[q->res*i + j] ;
      if (f_type == CMIN) {
        if (P.cmin <= trsh) {
          (*p->flow_ptr)[p->res*i + j] = Uerr ;
        }
      }
      else {
        if (P.cmax <= trsh) {
          (*p->flow_ptr)[p->res*i + j] = Uerr ;
        }
      }
    }
  }
}


/* 
           NAME : init_beaudet3(Ix,Ixx,Ixy)
   PARAMETER(S) :  Ix : 3 by 3 operator for 1st order x derivative;
                  Ixx : 3 by 3 operator for 2nd order x derivative;
                  Ixy : 3 by 3 operator for partial x and y derivatives.
 
        PURPOSE : Initialize the image operators defined by Beaudet
                  to obtain numerical derivatives.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 10 1990
*/

void init_beaudet3(Ix,Ixx,Ixy)
beaudet_t *Ix, *Ixx, *Ixy ;

{
  (*Ix).g[0][0] = -1.0 ; (*Ix).g[1][0] = 0.0 ; (*Ix).g[2][0] = 1.0 ;
  (*Ix).g[0][1] = -1.0 ; (*Ix).g[1][1] = 0.0 ; (*Ix).g[2][1] = 1.0 ;
  (*Ix).g[0][2] = -1.0 ; (*Ix).g[1][2] = 0.0 ; (*Ix).g[2][2] = 1.0 ;
  (*Ix).f = 6.0 ;
  (*Ix).m = 3 ;

  (*Ixx).g[0][0] = 1.0 ; (*Ixx).g[1][0] = -2.0 ; (*Ixx).g[2][0] = 1.0 ;
  (*Ixx).g[0][1] = 1.0 ; (*Ixx).g[1][1] = -2.0 ; (*Ixx).g[2][1] = 1.0 ;
  (*Ixx).g[0][2] = 1.0 ; (*Ixx).g[1][2] = -2.0 ; (*Ixx).g[2][2] = 1.0 ;
  (*Ixx).f = 3.0 ;
  (*Ixx).m = 3 ;

  (*Ixy).g[0][0] = 1.0 ;  (*Ixy).g[1][0] = 0.0 ; (*Ixy).g[2][0] = -1.0 ;
  (*Ixy).g[0][1] = 0.0 ;  (*Ixy).g[1][1] = 0.0 ; (*Ixy).g[2][1] = 0.0 ;
  (*Ixy).g[0][2] = -1.0 ; (*Ixy).g[1][2] = 0.0 ; (*Ixy).g[2][2] = 1.0 ;
  (*Ixy).f = 4.0 ;
  (*Ixy).m = 3 ;
}

/* 
           NAME : init_beaudet5(Ix,Ixx,Ixy)
   PARAMETER(S) :  Ix : operator for 1st order x derivative;
                  Ixx : operator for 2nd order x derivative;
                  Ixy : operator for partial x and y derivatives.
 
        PURPOSE : Initialize the image operators defined by Beaudet
                  to obtain numerical derivatives.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : September 13 1990
*/

void init_beaudet5(Ix,Ixx,Ixy)
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
           NAME : init_beaudet7(Ix,Ixx,Ixy)
   PARAMETER(S) :  Ix : operator for 1st order x derivative;
                  Ixx : operator for 2nd order x derivative;
                  Ixy : operator for partial x and y derivatives.
 
        PURPOSE : Initialize the image operators defined by Beaudet
                  to obtain numerical derivatives.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : September 13 1990
*/

void init_beaudet7(Ix,Ixx,Ixy)
beaudet_t *Ix, *Ixx, *Ixy ;

{
  (*Ix).g[0][0] = -3 ; (*Ix).g[1][0] = -2 ; (*Ix).g[2][0] = -1 ; 
  (*Ix).g[3][0] = 0 ;  (*Ix).g[4][0] = 1 ;  (*Ix).g[5][0] = 2 ;
  (*Ix).g[6][0] = 3 ;
  (*Ix).g[0][1] = -3 ; (*Ix).g[1][1] = -2 ; (*Ix).g[2][1] = -1 ; 
  (*Ix).g[3][1] = 0 ;  (*Ix).g[4][1] = 1 ;  (*Ix).g[5][1] = 2 ;
  (*Ix).g[6][1] = 3 ;
  (*Ix).g[0][2] = -3 ; (*Ix).g[1][2] = -2 ; (*Ix).g[2][2] = -1 ; 
  (*Ix).g[3][2] = 0 ;  (*Ix).g[4][2] = 1 ;  (*Ix).g[5][2] = 2 ;
  (*Ix).g[6][2] = 3 ;
  (*Ix).g[0][3] = -3 ; (*Ix).g[1][3] = -2 ; (*Ix).g[2][3] = -1 ; 
  (*Ix).g[3][3] = 0 ;  (*Ix).g[4][3] = 1 ;  (*Ix).g[5][3] = 2 ;
  (*Ix).g[6][3] = 3 ;
  (*Ix).g[0][4] = -3 ; (*Ix).g[1][4] = -2 ; (*Ix).g[2][4] = -1 ; 
  (*Ix).g[3][4] = 0 ;  (*Ix).g[4][4] = 1 ;  (*Ix).g[5][4] = 2 ;
  (*Ix).g[6][4] = 3 ;
  (*Ix).g[0][5] = -3 ; (*Ix).g[1][5] = -2 ; (*Ix).g[2][5] = -1 ;
  (*Ix).g[3][5] = 0 ;  (*Ix).g[4][5] = 1 ;  (*Ix).g[5][5] = 2 ;
  (*Ix).g[6][5] = 3 ;
  (*Ix).g[0][6] = -3 ; (*Ix).g[1][6] = -2 ; (*Ix).g[2][6] = -1 ;
  (*Ix).g[3][6] = 0 ;  (*Ix).g[4][6] = 1 ;  (*Ix).g[5][6] = 2 ;
  (*Ix).g[6][6] = 3 ;
  (*Ix).f = 196.0 ;
  (*Ix).m = 7 ;
      
  (*Ixx).g[0][0] = 5 ;  (*Ixx).g[1][0] = 0 ;  (*Ixx).g[2][0] = -3 ; 
  (*Ixx).g[3][0] = -4 ; (*Ixx).g[4][0] = -3 ; (*Ixx).g[5][0] = 0 ;
  (*Ixx).g[6][0] = 5 ;
  (*Ixx).g[0][1] = 5 ;  (*Ixx).g[1][1] = 0 ;  (*Ixx).g[2][1] = -3 ; 
  (*Ixx).g[3][1] = -4 ; (*Ixx).g[4][1] = -3 ; (*Ixx).g[5][1] = 0 ;
  (*Ixx).g[6][1] = 5 ;
  (*Ixx).g[0][2] = 5 ;  (*Ixx).g[1][2] = 0 ;  (*Ixx).g[2][2] = -3 ; 
  (*Ixx).g[3][2] = -4 ; (*Ixx).g[4][2] = -3 ; (*Ixx).g[5][2] = 0 ;
  (*Ixx).g[6][2] = 5 ;
  (*Ixx).g[0][3] = 5 ;  (*Ixx).g[1][3] = 0 ;  (*Ixx).g[2][3] = -3 ; 
  (*Ixx).g[3][3] = -4 ; (*Ixx).g[4][3] = -3 ; (*Ixx).g[5][3] = 0 ;
  (*Ixx).g[6][3] = 5 ;
  (*Ixx).g[0][4] = 5 ;  (*Ixx).g[1][4] = 0 ;  (*Ixx).g[2][4] = -3 ; 
  (*Ixx).g[3][4] = -4 ; (*Ixx).g[4][4] = -3 ; (*Ixx).g[5][4] = 0 ;
  (*Ixx).g[6][4] = 5 ;
  (*Ixx).g[0][5] = 5 ;  (*Ixx).g[1][5] = 0 ;  (*Ixx).g[2][5] = -3 ; 
  (*Ixx).g[3][5] = -4 ; (*Ixx).g[4][5] = -3 ; (*Ixx).g[5][5] = 0 ;
  (*Ixx).g[6][5] = 5 ;
  (*Ixx).g[0][6] = 5 ;  (*Ixx).g[1][6] = 0 ;  (*Ixx).g[2][6] = -3 ; 
  (*Ixx).g[3][6] = -4 ; (*Ixx).g[4][6] = -3 ; (*Ixx).g[5][6] = 0 ;
  (*Ixx).g[6][6] = 5 ;
  (*Ixx).f = 294.0 ;
  (*Ixx).m = 7 ;

  (*Ixy).g[0][0] = 9 ;  (*Ixy).g[1][0] = 6 ;  (*Ixy).g[2][0] = 3 ; 
  (*Ixy).g[3][0] = 0 ;  (*Ixy).g[4][0] = -3 ; (*Ixy).g[5][0] = -6 ;
  (*Ixy).g[6][0] = -9 ;
  (*Ixy).g[0][1] = 6 ;  (*Ixy).g[1][1] = 4 ;  (*Ixy).g[2][1] = 2 ; 
  (*Ixy).g[3][1] = 0 ;  (*Ixy).g[4][1] = -2 ; (*Ixy).g[5][1] = -4 ;
  (*Ixy).g[6][1] = -6 ;
  (*Ixy).g[0][2] = 3 ;  (*Ixy).g[1][2] = 2 ;  (*Ixy).g[2][2] = 1 ; 
  (*Ixy).g[3][2] = 0 ;  (*Ixy).g[4][2] = -1 ; (*Ixy).g[5][2] = -2 ;
  (*Ixy).g[6][2] = -3 ;
  (*Ixy).g[0][3] = 0 ;  (*Ixy).g[1][3] = 0 ;  (*Ixy).g[2][3] = 0 ; 
  (*Ixy).g[3][3] = 0 ;  (*Ixy).g[4][3] = 0 ;  (*Ixy).g[5][3] = 0 ;
  (*Ixy).g[6][3] = 0 ;
  (*Ixy).g[0][4] = -3 ; (*Ixy).g[1][4] = -2 ; (*Ixy).g[2][4] = -1 ; 
  (*Ixy).g[3][4] = 0 ;  (*Ixy).g[4][4] = 1 ;  (*Ixy).g[5][4] = 2 ;
  (*Ixy).g[6][4] = 3 ;
  (*Ixy).g[0][5] = -6 ; (*Ixy).g[1][5] = -4 ; (*Ixy).g[2][5] = -2 ; 
  (*Ixy).g[3][5] = 0 ;  (*Ixy).g[4][5] = 2 ;  (*Ixy).g[5][5] = 4 ;
  (*Ixy).g[6][5] = 6 ;
  (*Ixy).g[0][6] = -9 ; (*Ixy).g[1][6] = -6 ; (*Ixy).g[2][6] = -3 ; 
  (*Ixy).g[3][6] = 0 ;  (*Ixy).g[4][6] = 3 ;  (*Ixy).g[5][6] = 6 ;
  (*Ixy).g[6][6] = 9 ;
  (*Ixy).f = 784.0 ;
  (*Ixy).m = 7 ;
}

/* 
           NAME : init_kernel_a(ker)
   PARAMETER(S) : ker : the kernel.
 
        PURPOSE : Inits the kernel values.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 6 1990

*/

void init_kernel_a(ker) 
kernel_t *ker ;

{  
  (*ker).k[0] = 1.0 ;
  (*ker).k[1] = 5.0;
  (*ker).k[2] = 8.0 ;
  (*ker).k[3] = 5.0 ;
  (*ker).k[4] = 1.0 ;
  (*ker).m = 5 ;
  (*ker).f = 20.0 ;
}

/* 
           NAME : init_kernel_b(ker)
   PARAMETER(S) : ker : the kernel.
 
        PURPOSE : Inits the kernel values.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 6 1990

*/

void init_kernel_b(ker) 
kernel_t *ker ;

{  
  (*ker).k[0] = 1.0 ;
  (*ker).k[1] = 5.0 ;
  (*ker).k[2] = 8.0 ;
  (*ker).k[3] = 5.0 ;
  (*ker).k[4] = 1.0 ;
  (*ker).m = 5 ;
  (*ker).f = 10.0 ;
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

void init_list(h,q)
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

void insert_node(p,h,q) 
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
           NAME : inside(x,y,rx,ry)
   PARAMETER(S) : x,y   : image coordinates;
                  rx,ry : image size.
 
        PURPOSE : Returns the truth value of (x,y) in image 
                  of size rx*ry.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 9 1990
*/

int inside(x,y,rx,ry)
int x, y, rx, ry ;

{
  return(((x >= 0) && (x < rx) && (y >= 0) && (y < ry))) ;
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


void itoa(n,s)
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
                  sx,sy : image size;
                   r    : resolution of the array.

        PURPOSE : Reads rasterfile specified
                  by fn into bf.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 27 1989, modified March 20 1990.

*/

void pgetrast(fn,hd,bf,sx,sy,r) 
char *fn ;
unsigned char hd[H] ;
float bf[N1*N1] ;
int sx, sy, r ;

{ 
  int fd;
  int i, j;
  unsigned char c;

  if ((fd = open(fn,O_RDONLY)) > 0) { 
    if (read(fd,hd,H) == H) { 
      for (i = 0 ; i < sy ; i++) {
        for (j = 0 ; j < sx ; j++) {
          if (read(fd,&c,sizeof(c)) != 1) {
            error(5) ;
          }
          bf[r*i + j] = (float)c ;
        }
      }
    }
    else {
      error(6) ;
    }
  }
  else {
    error(7) ;
  }
  close(fd) ;
}

/* 
           NAME : Bpgetrast(fn,bf,sx,sy,r) 
   PARAMETER(S) : fn    : filename;
                  bf    : Pointer on a 2D array containing the image;
                  sx,sy : image size;
                   r    : resolution of the array.

        PURPOSE : Reads rasterfile without header specified
                  by fn into bf.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 27 1989, modified March 20 1990.

*/

void Bpgetrast(fn,bf,sx,sy,r) 
char *fn ;
float bf[N1*N1] ;
int sx, sy, r ;

{ int fd, i, j ;
  unsigned char c ;

  if ((fd = open(fn,O_RDONLY)) > 0) { 
    for (i = 0 ; i < sy ; i++) {
      for (j = 0 ; j < sx ; j++) {
        if (read(fd,&c,sizeof(c)) != 1) {
          error(5) ;
        }
        bf[r*i + j] = (float)c ;
      }
    }
  }
  else {
    error(7) ;
  }
  close(fd) ;
}

/* 
           NAME : pputrast(fn,hd,bf,sx,sy,r) 
   PARAMETER(S) : fn     : filename;
                  hd     : image header (raster file);
                  bf     : Pointer on a 2D array containing the image;
                  sx, sy : image size;
                   r     : resolution of the array.
 
        PURPOSE : Writes the contents of bf into a 
                  rasterfile specified by fn.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : October 27 1989, modified March 20 1990.

*/

void pputrast(fn,hd,bf,sx,sy,r) 
char *fn ;
unsigned char hd[H] ;
float bf[N1*N1] ;
int sx, sy, r ;

{ int fd, i, j ;
  unsigned char c ;

  if ((fd = creat(fn,0600)) > 0) { ;
    if (write(fd,hd,H) == H) { 
      for (i = 0 ; i < sy ; i++) { 
        for (j = 0 ; j < sx ; j++) { 
          c = (unsigned char)bf[r*i + j] ;
          if (write(fd,&c,sizeof(c)) != 1) {
            error(8) ;
          }
        }
      }
    }
    else {
      error(9) ; 
    }
  }
  else {
    error(10) ;
  }
  close(fd) ;
}

/* 
           NAME : project(U,r)
   PARAMETER(s) : U : vector from next coarser level;
                  r : resolution reduction ratio. 

        PURPOSE : Projects down the 2D vector U to the next finer 
                  level as an initial estimate for the match criterion.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : June 25 1990
*/

disp_vect_t project(U,r)
disp_vect_t U ;
int r ;

{  
  U.x = U.x*r ;
  U.y = U.y*r ;
  if (U.x > 0.0) {
    U.x = (int)(U.x + 0.5) ;
  }
  else {
    U.x = (int)(U.x - 0.5) ;
  }
  if (U.y > 0.0) {
    U.y = (int)(U.y + 0.5) ;
  }
  else {
    U.y = (int)(U.y - 0.5) ;
  }
  return(U) ;
}

/* 
           NAME : rotate(x,y,t)
   PARAMETER(S) : x,y : vector points;
                    t : rotational angle (radians).
 
        PURPOSE : Rotates the vector (x,y) of t radians.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 11 1990
*/

void rotate(x,y,t) 
float *x, *y, t ;

{ float x1, y1 ;

  x1 = *x ; y1 = *y ;
  *x =  x1*cos(t) + y1*sin(t) ;
  *y = -x1*sin(t) + y1*cos(t) ;
}

/* 
           NAME : sample(p1,p2) 
   PARAMETER(S) : p1 : pointer on an arbitrary node of the pyramid;
                  p2 : pointer on the next coarser level to p1.
 
        PURPOSE : Subsamples from finer to coarser resolution.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 26 1990

*/

void sample(p1,p2)
qnode_ptr_t p1, p2 ;

{ int f, i, j ;
  
  f = p1->res/p2->res ;
  for (i = 0 ; i < p2->sizy ; i++) {
    for (j = 0 ; j < p2->sizx ; j++) {
      (*p2->gauss_ptr)[p2->res*i + j] = 
      (*p1->lap_ptr)[p1->res*(i*f + 1) + j*f + 1] ;
    }
  }
}

/* 
           NAME : ssd(f,g,x1,y1,x2,y2)
   PARAMETER(S) : f : laplacian image containing the pixel (x1,y1);
                  g : laplacian image containing the pixel (x2,y2);
                  x1,y1 : pixel coordinates in f; 
                  x2,y2 : pixel coordinates in g.

        PURPOSE : Computes the Gaussian weighted sum of differences
                  between (x1,y1) and (x2,y2).

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 9 1990

*/

float ssd(f,g,x1,y1,x2,y2) 
qnode_ptr_t f, g ;
int x1, y1, x2, y2 ;

{ extern kernel_t kera ;
  int k1, k2, xr1, yr1, xr2, yr2 ;
  float s ;

  s = 0.0 ;
  for (k1 = -kera.m/2 ; k1 <= kera.m/2; k1++) {
    if (x1 + k1 < 0) {
      xr1 = 0 ;
    }
    else {
      if (x1 + k1 >= f->sizy) {
        xr1 = f->sizy - 1 ;
      }
      else {
        xr1 = x1 + k1 ;
      }
    }
    if (x2 + k1 < 0) {
      xr2 = 0 ;
    }
    else {
      if (x2 + k1 >= g->sizy) {
        xr2 = g->sizy - 1 ;
      }
      else {
        xr2 = x2 + k1 ;
      }
    }
    for (k2 = -kera.m/2 ; k2 <= kera.m/2 ; k2++) {
      if (y1 + k2 < 0) {
        yr1 = 0 ; 
      }
      else {
        if (y1 + k2 >= f->sizx) {
          yr1 = f->sizx - 1 ;
        }
        else {
          yr1 = y1 + k2 ;
        }
      }
      if (y2 + k2 < 0) {
        yr2 = 0 ; 
      }
      else {
        if (y2 + k2 >= g->sizx) {
          yr2 = g->sizx - 1 ;
        }
        else {
          yr2 = y2 + k2 ;
        }
      }
      s = s + (float)kera.k[k1+kera.m/2]*(float)kera.k[k2+kera.m/2]*
      pow(((*f->lap_ptr)[f->res*xr1 + yr1] - 
           (*g->lap_ptr)[g->res*xr2 + yr2]),2.0) ;
    }
  }
  return(s/(kera.f*kera.f)) ;
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

void usage()

{ 
  fprintf(stderr,"Usage: anandan input-path output-path seq-name. -N n1 n2 [-W n] [-L n] [-I n] [-NR] [-P] [-F n.n] [-B cols rows] [-C corr-vel-file nbins increment]\n") ;
  fprintf(stderr,"-N n1 n2 : frame numbers n1 and n2\n") ;
  fprintf(stderr,"[-W n]   : correlation window size (3, 5 or 7)\n") ;
  fprintf(stderr,"           default value is 3\n") ;
  fprintf(stderr,"[-L n]   : number of hierarchical levels\n") ;
  fprintf(stderr,"           default value is 1\n") ;
  fprintf(stderr,"[-I n]   : number of relaxation iterations\n") ;
  fprintf(stderr,"           default value is 10\n") ;
  fprintf(stderr,"[-NR]    : no regularization\n") ;
  fprintf(stderr,"         : default is regulariztion\n") ;
  fprintf(stderr,"[-P]     : pixel accuracy only\n") ;
  fprintf(stderr,"         : default is sub-pixel accuracy\n") ;
  fprintf(stderr,"[-F n.n] : filters out unreliable estimates\n") ;
  fprintf(stderr,"           using confidence measure Cmin\n") ;
  fprintf(stderr,"[-B cols rows]\n") ;
  fprintf(stderr,"         : for binary input files without header\n") ;
  fprintf(stderr,"[-C corr-vel-file nbins increment]\n") ;
  fprintf(stderr,"         : error histogram\n") ;
  fprintf(stderr,"           corr-vel-file : file of correct velocities\n") ;
  fprintf(stderr,"           nbins : number of bins for histogram\n") ;
  fprintf(stderr,"           increment : value of gap between bins\n") ;
  exit(-1) ;
}

/* 
           NAME : vec_rotate(U,t)
   PARAMETER(S) : U : 2D vector;
                  t : rotational angle (radians).
 
        PURPOSE : Rotates the vector U(x,y) of t radians.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : May 3 1990
*/

disp_vect_t vec_rotate(U,t) 
disp_vect_t U ;
float t ;

{ disp_vect_t Up ;

  Up.x = U.x ; Up.y = U.y ;
  U.x  =  Up.x*cos(t) + Up.y*sin(t) ;
  U.y  = -Up.x*sin(t) + Up.y*cos(t) ;
  return(U) ;
}

/*
           NAME : write_velocity(fn,p,div)
   PARAMETER(S) : fn  : file name for flow field;
                  p   : node pointer on a flow field;
                  div : flow estimates divisor.

        PURPOSE : Output field using Travis Burkitt's format.

         AUTHOR : Travis Burkitt, updated by Steven Beauchemin
             AT : University of Western Ontario
           DATE : May 7 1990
*/

void write_velocity(fn,p,div)
qnode_ptr_t p ;
char *fn ;
float div ;

{ float x, y ; 
  int i, j, fdf, bytes ;
  if ((fdf=creat(fn,0600)) < 1) {
    error(7) ;
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
      if (x != 100.0 && y != 100.0) {
        x /= div ;
        y /= div ;
      }
      write(fdf,&x,4) ;
      write(fdf,&y,4) ;
      bytes += 8 ;
    }
  }
  close(fdf) ;
}

/* 
           NAME : best_match(im1,im2,ssd_v,x,y,cx,cy)
   PARAMETER(S) : im1,im2 : Contiguous time-varying laplacian images;
                    ssd_v : ssd value to be centered at the
                            best match computed from the 3 by 3 area 
                            around (cx,cy);
                      x,y : image coordinates of the pixel in im1 to
                            be matched;
                    cx,cy : image coordinates in im2 of the center 
                            of the n by n search area of the best
                            match for pixel (x,y) from im1.
 
        PURPOSE : Finds the best match (using the minimization of
                  the ssd surface) of pixel (x,y) in laplacian image 1
                  in the n by n search area of image 2 centered 
                  at image coordinates (cx,cy) in image 2.
                  Also returns the ssd value at the best match.
                  Note: displacement is computed to pixel accuracy.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 9 1990

*/

disp_vect_t best_match(im1,im2,ssd_v,x,y,cx,cy)
qnode_ptr_t im1, im2 ;
float *ssd_v ;
int x, y, cx, cy ;

{ float ssd(), v[SSDVAL] ;
  disp_vect_t U[SSDVAL], V ;
  int pos, dx, dy, h, i, j, k ;
  extern beaudet_t Ix ;

  k = 0 ; 
  h = Ix.m/2 ;
  *ssd_v = MAXSSD ;
  for(dx = -h ; dx <= h ; dx++) {                /* search in n by n area */
    for (dy = -h ; dy <= h ; dy++) {                    /* around (cx,cy) */
      if (inside(cx + dx,cy + dy,im2->sizy,im2->sizx)) {
        v[k] = ssd(im1,im2,x,y,cx + dx,cy + dy) ;  /* cumulate ssd values */
        U[k].x = (cx + dx) - x ;
        U[k].y = (cy + dy) - y ;
        if (v[k] < *ssd_v) {                      /* run SSD minimization */
          *ssd_v = v[k] ;
          pos = k ;
        }
        k++ ;
      }
    }
  }
  j = 1 ; 
  V = U[pos] ;  
  for (i = 0 ; i < k ; i++) {          /* compute variation of ssd values */
    if ((i != pos) && (fabs(v[pos] - v[i]) == 0.0)){
      V.x += U[i].x ;            /* average disp. estimates with same SSD */
      V.y += U[i].y ;                /* NOTE: this case is extremely rare */
      j++ ;
    }
  }
  V.x = V.x/(float)j ;
  V.y = V.y/(float)j ; 
  return(V) ;
}

/* 
           NAME : conf_measure(im1,im2,fl,prm,x,y,sp)
   PARAMETER(S) : im1,im2 : time-varying Laplacian images;
                       fl : pointer on a displacement-field node;
                      prm : parameter-field of the disp. field 
                            (theta,cmin,cmax);
                      x,y : disp.-field coordinates for which the 
                            computation is achieved;
                       sp : adds sub-pixel accurate disp. when true.
 
        PURPOSE : Computes the principal curvatures of the ssd
                  surface, the rotational angle of the eigen
                  basis from the (x,y) space, sub-pixel
                  accuracy displacement in the eigenbasis
                  and the confidence measures defined by
                  Anandan.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 9 1990
*/

void conf_measure(im1,im2,fl,prm,x,y,sp) 
qnode_ptr_t im1, im2, fl, prm ;
int x, y, sp ;

{ disp_vect_t vec_rotate(), U, Sd ;
  param_t P ;
  search_area_t ssd_s ;
  float ssd(), dx(), dy(), dxx(), dyy(), dxy(), Sx, Sy, Sxx, Syy, Sxy,
  Cmax, Cmin ;
  int abort, h, i, j ;
  extern beaudet_t Ix ;

  abort = FALSE ;
  h = Ix.m/2 ;
  U = (*fl->flow_ptr)[fl->res*x + y] ;
  for(i = -h ; i <= h ; i++) {       /* computation of 3 by 3 ssd surface */
    for(j = -h ; j <= h ; j++) {
      if (inside(x + (int)U.x + i,y + (int)U.y + j,fl->sizy,fl->sizx)) {
        ssd_s.k[i+h][j+h] = ssd(im1,im2,x,y,x + (int)U.x + i,y + (int)U.y +j) ;
      }
      else {
        ssd_s.k[i+h][j+h] = -1.0 ;
        abort = TRUE ;
      }
    }
  }

  if (!abort) { 
    Sx = dx(&ssd_s) ;       /* derivatives from Beaudet's image operators */
    Sy = dy(&ssd_s) ;
    Sxx = dxx(&ssd_s) ;
    Syy = dyy(&ssd_s) ;
    Sxy = dxy(&ssd_s) ;

             /* computation of the principal curvatures along (Emin,Emax) */

    if ((Sxx - Syy)*(Sxx - Syy) + 4.0*Sxy*Sxy >= 0.0) {
      Cmax = 0.5*((Syy + Sxx) + sqrt((Sxx - Syy)*(Sxx - Syy) + 4.0*Sxy*Sxy)) ;
      Cmin = 0.5*((Syy + Sxx) - sqrt((Sxx - Syy)*(Sxx - Syy) + 4.0*Sxy*Sxy)) ;
    }
    else {
      Cmin = 0.0 ;
      Cmax = 0.0 ;
    }

    P.cmax = Cmax/(K1 + K2*ssd_s.k[h][h] + K3*Cmax) ;
    P.cmin = Cmin/(K1 + K2*ssd_s.k[h][h] + K3*Cmin) ;
    P.th = atan2(Cmin - Sxx,Sxy) ;           /* angle from Emin to x axis */
    rotate(&Sx,&Sy,P.th) ;               /* derivatives along (Emin,Emax) */
  
    if (Cmin != 0.0) {
      Sd.x = Sx/Cmin ;           /* determination of the extremum along A */
    }
    else {
      Sd.x = 0.0 ;
    }
    if (Cmax != 0.0) {
      Sd.y = Sy/Cmax ;           /* determination of the extremum along B */
    }
    else {
      Sd.y = 0.0 ;
    }
              /* use the second order derivatives to find if the extremum */
                                      /* is a minimum of the ssd function */
  
    if ((fabs(Sd.x) > 1.0) || (Cmin <= 0.0)) {
      Sd.x = 0.0 ;
      P.cmin = 0.0 ;
    }
    if ((fabs(Sd.y) > 1.0) || (Cmax <= 0.0)) {
      Sd.y = 0.0 ;
      P.cmax = 0.0 ;
    }
    if (sp) {
      U = vec_rotate(U,-P.th) ;
      U.x = U.x + Sd.x ;                  /* sub pixel displacement added */
      U.y = U.y + Sd.y ;      
      U = vec_rotate(U,P.th) ;
    } 
  
    (*fl->flow_ptr)[fl->res*x + y] = U ;
    (*prm->param_ptr)[prm->res*x + y] = P ;
  }
  else {
    U.x = 0.0 ; U.y = 0.0 ;
    (*fl->flow_ptr)[fl->res*x + y] = U ;
    P.cmin = 0.0 ; P.cmax = 0.0 ;
    (*prm->param_ptr)[prm->res*x + y] = P ;
  }
}

/* 
           NAME : convolve(p,ker)
   PARAMETER(S) :   p : pointer on the current node of
                        the image to be convolved;
                  ker : kernel (mask) used for convolution.
 
        PURPOSE : Performs a 2D convolution on (*p->gauss_ptr)[i][j]
                  using a 1D kernel (ker) and puts the result in
                  (*p->lap_ptr)[i][j]. 

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 20 1990

*/

void convolve(p,ker) 
qnode_ptr_t p ;
kernel_t ker ;


{ qnode_ptr_t create_node(), t ;
  int i, j, k, m, e ;
  float s ;

  t = create_node(0,p->res,p->sizx,p->sizy) ;
  alloc_gauss(t) ;
  e = (float)ker.m/2.0 ;
  
  for (i = 0 ; i < p->sizy ; i++) {
    for (j = 0 ; j < p->sizx ; j++) {
      s = 0.0 ;
      for (k = 0 ; k < ker.m ; k++) {
        m = j + (k - e) ;
        if (m >= p->sizx) {
          m = p->sizx - (m%(p->sizx - 1) + 1) ;
        }
        else {
          if (m < 0) {
            m = -m ;
          }
        }
        s = s + (*p->gauss_ptr)[p->res*i + m]*ker.k[k] ;
      }
      (*t->gauss_ptr)[t->res*i + j] = s/ker.f ;
    }
  }

  for (i = 0 ; i < p->sizy ; i++) {
    for (j = 0 ; j < p->sizx ; j++) {
      s = 0.0 ;
      for (k = 0 ; k < ker.m ; k++) {
        m = i + (k - e) ;
        if (m >= p->sizy) {
          m = p->sizy - (m%(p->sizy - 1) + 1) ;
        }
        else {
          if (m < 0) {
            m = -m ;
          }
        }
        s = s + (*t->gauss_ptr)[t->res*m + j]*ker.k[k] ;
      }
      (*p->lap_ptr)[p->res*i + j] = s/ker.f ;
    }
  }
  dealloc_gauss(t) ;
  free(t) ;
}

/* 
           NAME : create_pyramid(h,q,n,r,sx,sy,f) 
   PARAMETER(S) : h      : head list pointer;
                  q      : tail list pointer;
                  n      : number of levels;
                  r      : resolution at level 0;
                  sx, sy : size of the image at level 0;
                  f      : resolution reduction rate (2 for 2:1).
 
        PURPOSE : Creation of a pyramid with n levels of
                  resolution r  at level 0 and reduced by
                  f at every level.  Each level is 
                  attached to a doubly linked list node.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 18 1990

*/

void create_pyramid(h,q,n,r,sx,sy,f)
qnode_ptr_t *h, *q ;
int n, r, sx, sy, f ;

{ qnode_ptr_t create_node(), p ;
  int i ;

 for (i = 0 ; i < n ; i++) {
    p = create_node(i,r,sx,sy) ;
    insert_node(p,h,q) ;
    alloc_image(p) ;
    sx = sx/f ;
    sy = sy/f ;
    r = r/f ;
  }
}

/* 
           NAME : delete_pyramid(h,q)
   PARAMETER(S) : h : head list pointer;
                  q : tail list pointer.
 
        PURPOSE : Deletion of all nodes in the
                  list and the corresponding images.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 17 1990

*/

void delete_pyramid(h,q)
qnode_ptr_t *h, *q ;

{ qnode_ptr_t t ;

  t = *h ;
  while (*h != (qnode_ptr_t)NULL) {
    dealloc_image(t) ;
    delete_node(t,h,q) ; 
    t = *h ;
  }
}

/* 
           NAME : dump_flow(h,q,f,div)
   PARAMETER(S) : h, q : Head and queue pointers on
                         an optic-field pyramid;
                     f : file name;
                   div : flow estimates divisor.
 
        PURPOSE : Dumps all levels of a pyramidal optic-field
                  under the names f"0", ... , f"n-1".

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : June 25 1990
*/

void dump_flow(h,q,f,div)
qnode_ptr_t *h, *q ;
char *f ;
float div ;
{ 
    qnode_ptr_t t ;

  /* while (*q != (qnode_ptr_t)NULL) {
    t = *q ;
    itoa(t->level,s) ;
    concat(f,s,fname) ;
    write_velocity(fname,t,div) ;
    dealloc_flow(t) ;
    dealloc_param(t) ;
    delete_node(t,h,q) ;
  } */
  t = *q ;
  write_velocity(f,t,div) ;
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

void error(n)
int n ;

{ 
  switch(n) {
    
     case 1 : fprintf(stderr,"error %d : wrong number of arguments\n",n) ;
              break ;
     case 2 : fprintf(stderr,"error %d : invalid option\n",n) ;
              break ;
     case 3 : fprintf(stderr,"error %d : incompatible options\n",n) ;
              break ;
     case 4 : fprintf(stderr,"error %d : maximum number of levels is %d\n",n,MAXLEVEL) ;
              break ;
     case 5 : fprintf(stderr,"error %d : read error\n",n) ;
              exit(-1) ;
              break ;
     case 6 : fprintf(stderr,"error %d : header read error\n",n) ;
              exit(-1) ;
              break ;
     case 7 : fprintf(stderr,"error %d : file open error\n",n) ;
              exit(-1) ;
              break ;
     case 8 : fprintf(stderr,"error %d : write error\n",n) ;
              exit(-1) ;
              break ;
     case 9 : fprintf(stderr,"error %d : header write error\n",n) ;
              exit(-1) ;
              break ;
     case 10: fprintf(stderr,"error %d : file create error\n",n) ;
              exit(-1) ;
              break ;
     case 11: fprintf(stderr,"error %d : undefined error\n",n) ; 
              break ;
     case 12: fprintf(stderr,"error %d : undefined error\n",n) ;
              break ;
     case 13: fprintf(stderr,"error %d : maximum number of bins is %d\n",
              n,N_BINS) ;
              exit(-1) ;
              break ;
     case 14: fprintf(stderr,"error %d : image size must be between %d and %d\n", n,N1,N3) ;
              exit(-1) ;
              break ;
     case 15: fprintf(stderr,"error %d : incompatible image sizes\n",n) ;
              exit(-1) ;
              break ; 
     case 16: fprintf(stderr,"error %d : window sizes are 3, 5 or 7\n",n) ;
              break ;
     case 17: fprintf(stderr,"error %d : -N required\n",n) ;
              break ;
     case 18: fprintf(stderr,"error %d : frame numbers must be different\n",n) ;
              exit(-1) ;
              break ;
    default : fprintf(stderr,"error %d : undefined error\n",n) ;
              exit(-1) ;
              break ;
  }
  usage() ;
}

/* 
           NAME : filenames(fn,n1,n2,f1,f2,div)
   PARAMETER(S) : fn  : image sequence filename;
                  n1  : number of first frame;
                  n2  : number of second frame;
                  f1  : complete filename of image 1;
                  f2  : complete filename of image 2;
                  div : divisor of flow estimates.
 
        PURPOSE : Determines the correct filenames given
                  the stem name and frame numbers.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : June 23 1992
*/

void filenames(fn,n1,n2,f1,f2,div)
string fn, f1, f2 ;
int n1, n2 ;
float *div ;

{ string tf1, tf2 ;

  itoa(n1,tf1) ;
  itoa(n2,tf2) ;
  concat(fn,tf1,f1) ;
  concat(fn,tf2,f2) ;
  *div = n2 - n1 ;
}

/* 
           NAME : laplacian(p)
   PARAMETER(S) : p : node pointer.
 
        PURPOSE : Computes the laplacian of an image using
                  a DOG approximation.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : March 26 1990

*/

void laplacian(p)
qnode_ptr_t p ;

{ extern kernel_t kerb ;
  qnode_ptr_t create_node(), q ;
  int i, j ;

  q = create_node(0,p->res,p->sizx,p->sizy) ; 
  alloc_lap(q) ;
  q->gauss_ptr = p->lap_ptr ;
  for (i = 0 ; i < q->sizy ; i++) {
    for (j = 0 ; j < q->sizx ; j++) {
      if ((i%2 == 0) || (j%2 == 0)) {
        (*q->gauss_ptr)[q->res*i + j] = 0.0 ;
      }
    }
  }
  convolve(q,kerb) ;
  p->lap_ptr = q->lap_ptr ;
  dealloc_gauss(q) ;
  free(q) ;
  for (i = 0 ; i < p->sizy ; i++) {
    for (j = 0 ; j < p->sizx ; j++) {
      (*p->lap_ptr)[p->res*i + j] = 
      (*p->gauss_ptr)[p->res*i + j] - (*p->lap_ptr)[p->res*i + j] ;
    }
  }
}

/* 
           NAME : psi_error(ve,va,div)
   PARAMETER(S) : ve : estimated flow vector;
                  va : accurate flow vector.
 
        PURPOSE : computes angular error of ve with respect to va
                  using Fleet [90] angular error metric

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : November 7 1990
*/

float psi_error(ve,va,div)
disp_vect_t ve, va ;
float div ;

{ float norm(), v ;
  float VE[3], VA[3] ;

  VE[0] = ve.x/div ;
  VE[1] = ve.y/div ;
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
           NAME : raster_size(fn,x,y,a)
   PARAMETER(S) : fn : raster file name; 
                   x : raster width;
                   y : raster height;
                   a : suitable array size (power of 2).
 
        PURPOSE : Reads the header of raterfile fn to get the
                  size. Computes the suitable array size for
                  the image (array size is a power of 2).

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 22 1992
*/

void raster_size(fn,x,y,a)
char *fn ;
int *x, *y, *a ;

{ raster_t hd ;
  int fd, max ;

  if ((fd = open(fn,O_RDONLY)) > 0) {
    if (read(fd,&hd,sizeof(hd)) == sizeof(hd)) {
      *x = hd.width ;
      *y = hd.height ;
      if (*x > *y) {
        max = *x ;
      }
      else {
        max = *y ;
      }
      if (max <= N1 && max > N2) {
        *a = N1 ;
      }
      else if (max <= N2 && max > N3) {
        *a = N2 ;
      }
      else if (max <= N3 && max > N4) {
        *a = N3 ;
      }
      else {
        error(14) ;
      }
    }
    else {
      error(6) ;
    }
  }
  else {
    error(7) ;
  }
  close(fd) ;
}

/* 
           NAME : binary_size(x,y,a)
   PARAMETER(S) : x : raster width;
                  y : raster height;
                  a : suitable array size (power of 2).
 
        PURPOSE : Computes the suitable array size for
                  the image (array size is a power of 2).

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 22 1992
*/

void binary_size(x,y,a)
int x, y, *a ;

{ int max ;

  if (x > y) {
    max = x ;
  }
  else {
    max = y ;
  }
  if (max <= N1 && max > N2) {
    *a = N1 ;
  }
  else if (max <= N2 && max > N3) {
    *a = N2 ;
  }
  else if (max <= N3 && max > N4) {
    *a = N3 ;
  }
  else {
    error(14) ;
  }
}

/* 
           NAME : valid_option(argc,argv,in_path,out_path,n1,n2,level,h_type,
                  histo,sm,sp, sf,pyr,ssd,f_type,trsh,v_fname,c_fname,h_fname,
                  nbin,incr,iter,binary,row,col) 

   PARAMETER(S) : argc : argument count;
                  argv : argument values;
               in_path : path name for input data;         
              out_path : path name for output data;
                 n1,n2 : frame numbers;
                 level : number of hierarchical levels;
                h_type : type of histogram (on Cmin or Cmax);
                 histo : error histogram option;
                    sm : regularization option;
                    sp : sub-pixel accuracy option;
                    sf : filtering option;
                   pyr : pyramid option (lap or gauss) ;
                   ssd : size of ssd window (3,5, or 7) ;
                f_type : type of filtering (on Cmin or Cmax);
                  trsh : threshold value for filtering;
               i_fname : input file name;
               v_fname : velocity file name;
               c_fname : correct velocities file name;
               h_fname : hisogram data file name ;
                  nbin : number of bins in histo 1;
                  incr : increment in histo 1;
                  iter : number of iterations for relaxation;
                binary : input files without header;
                   row : number of rows in input files;
                   col : number of cols in input files. 
 
        PURPOSE : Validates line arguments.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : June 25 1990
*/

void valid_option(argc,argv,in_path,out_path,n1,n2,level,h_type,histo,sm,sp,sf,pyr,
ssd,f_type,trsh,i_fname,v_fname,c_fname,h_fname,nbin,incr,iter,binary,row,col)
char *argv[] ;
int argc, *n1, *n2, *iter, *level, *h_type, *histo, *sm, *sp, *sf, *pyr, 
    *f_type, *nbin, *ssd, *binary, *row, *col ;
float *trsh, *incr ;
string in_path, out_path, i_fname, v_fname, c_fname, h_fname ;

{ string ext, t0, str_trsh ;
  int i, fn ;

  if (argc == 1) {
    usage() ;
  }
  if ((argc <= 24) && (argc >= 4)) {
    *binary = FALSE ;
    *sm = TRUE ;
    *sp = TRUE ;
    *sf = FALSE ;
    *pyr = LAP ;
    *level = DEFLEVEL ;
    *histo = FALSE ;
    fn = FALSE ;
    *ssd = 3 ;
    *trsh = 0.0 ;
    *nbin = 0 ;
    *iter = RELAXITER ;
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
    if (strcmp("-NR",argv[i]) == 0) {
      *sm = FALSE ;
      i += 1 ;
    }
    else {
      if (strcmp("-F",argv[i]) == 0) {
        if (i + 1 < argc) {
          sscanf(argv[i+1],"%f",trsh) ;
          strcpy(str_trsh,argv[i+1]) ; 
          if (*sf == TRUE) { 
            error(3) ;
          }
          *sf = TRUE ;
          *f_type = CMIN ;
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
            *h_type = CMIN ;
            *histo = TRUE ;
            i += 4 ;
          }
          else {
            error(1) ;
          }
        } 
        else {
          if (strcmp("-P",argv[i]) == 0) {
            *sp = FALSE ; 
            i += 1 ;
          }
          else {
            if (strcmp("-L",argv[i]) == 0) {
              if ( i + 1 < argc) {
                sscanf(argv[i+1],"%d",level) ;
                i += 2 ;
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
                if (strcmp("-W",argv[i]) == 0) {
                  if (i + 1 < argc) {
                    sscanf(argv[i+1],"%d",ssd) ;
                    i += 2 ;
                  }
                  else {
                    error(1) ;
                  }
                }
                else {
                  if (strcmp("-N",argv[i]) == 0) {
                    if (i + 2 < argc) {
                      sscanf(argv[i+1],"%d",n1) ;
                      sscanf(argv[i+2],"%d",n2) ;
                      fn = TRUE ;
                      i += 3 ;
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
    }
  }
  if (*nbin > N_BINS) {
    error(13) ;
  }
  if (*level > MAXLEVEL) {
    error(4) ;
  }
  if (*ssd != 3 && *ssd != 5 && *ssd != 7) {
    error(16) ;
  }
  if (!fn) {
    error(17) ;
  }
  if (*n1 == *n2) {
    error(18) ;
  }
  concat(argv[1],"/",in_path) ;
  concat(argv[2],"/",out_path) ;
  
  concat(in_path,argv[3],i_fname) ;
  
  strcpy(ext,"-n") ;
  itoa(*n1,t0) ;
  concat(ext,t0,ext) ;
  concat(ext,"-",ext) ;
  itoa(*n2,t0) ;
  concat(ext,t0,ext) ;
  concat(ext,"-w",ext) ;
  itoa(*ssd,t0) ;
  concat(ext,t0,ext) ;
  concat(ext,"-l",ext) ;
  itoa(*level,t0) ; 
  concat(ext,t0,ext) ;
  concat(ext,"-i",ext) ;
  itoa(*iter,t0) ;
  concat(ext,t0,ext) ;
  if (!(*sm)) {
    concat(ext,"-r",ext) ;
  }
  if (!(*sp)) { 
    concat(ext,"-s",ext) ;
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
           NAME : valid_size(x,y)
   PARAMETER(S) : x : image sizes in X for frames 1 and 2;
                  t : image sizes in Y for frames 1 and 2.
 
        PURPOSE : Validates the compatibility of image sizes

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 24 1992 ;
*/

void valid_size(x,y)
int x[N_FRAME], y[N_FRAME] ;

{ if (x[0] != x[1] || y[0] != y[1]) {
    error(15) ;
  } 
}

/* 
           NAME : relax(f,g,iter)
   PARAMETER(S) : f : pointer on a vector field node;
                  g : pointer on a parameter field node;
               iter : number of iterations for relaxation.
 
        PURPOSE : Applies a smoothness constraint to the vector field f
                  by using the Gauss-Seidel relaxation algorithm.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 15 1990
*/

void relax(f,g,iter) 
qnode_ptr_t f, g ;
int iter ;

{ qnode_ptr_t create_node(), p, q ;
  disp_vect_t vec_rotate(), Ukp, D, Emax, Emin ;
  float u, v, pCmax, pCmin ;
  int x, y, k, w ;

  p = create_node(0,f->res,f->sizx,f->sizy) ;
  q = create_node(0,f->res,f->sizx,f->sizy) ;
  alloc_flow(p) ;
  alloc_flow(q) ;

                                            /* copy f into q so that Uo = D */

  for(x = 0 ; x < f->sizy ; x++) {
    for (y = 0 ; y < f->sizx ; y++) {
      (*q->flow_ptr)[q->res*x + y] = (*f->flow_ptr)[f->res*x + y] ;
    }
  }

                     /* for k = 0 to convergence, compute Uk prime and Uk+1 */

  for(k = 0 ; k < iter ; k++) {
    for (x = 0 ; x < f->sizy ; x++) {
      for (y = 0 ; y < f->sizx ; y++) {

        Ukp.x = 0.0 ;                            /* computation of Uk prime */
        Ukp.y = 0.0 ;
        w = 0.0 ;

        if(inside(x + 1,y,f->sizy,f->sizx)) {
          u = (*q->flow_ptr)[q->res*(x + 1) + y].x ;
          v = (*q->flow_ptr)[q->res*(x + 1) + y].y ;
          Ukp.x = Ukp.x + u ;
          Ukp.y = Ukp.y + v ;
          w++ ;
        }
        if(inside(x - 1,y,f->sizy,f->sizx)) {
          u = (*q->flow_ptr)[q->res*(x - 1) + y].x ; 
          v = (*q->flow_ptr)[q->res*(x - 1) + y].y ;
          Ukp.x = Ukp.x + u ;
          Ukp.y = Ukp.y + v ;
          w++ ;
        }
        if(inside(x,y + 1,f->sizy,f->sizx)) {
          u = (*q->flow_ptr)[q->res*x + y + 1].x ; 
          v = (*q->flow_ptr)[q->res*x + y + 1].y ;
          Ukp.x = Ukp.x + u ;
          Ukp.y = Ukp.y + v ;
          w++ ;
        }
        if(inside(x,y - 1,f->sizy,f->sizx)) {
          u = (*q->flow_ptr)[q->res*x + y - 1].x ; 
          v = (*q->flow_ptr)[q->res*x + y - 1].y ;
          Ukp.x = Ukp.x + u ;
          Ukp.y = Ukp.y + v ;
          w++ ;
        }

                          /* average of the neighborhood of (u(x,y),v(x,y)) */

        Ukp.x = Ukp.x/(float)w ;
        Ukp.y = Ukp.y/(float)w ;
        (*p->flow_ptr)[p->res*x + y] = Ukp ;
      }
    }

    for(x = 0 ; x < f->sizy ; x++) {                 /* computation of Uk+1 */
      for(y = 0 ; y < f->sizx ; y++) {

        Ukp = (*p->flow_ptr)[p->res*x + y] ;
        pCmax = (*g->param_ptr)[g->res*x + y].cmax/((*g->param_ptr)[g->res*x + y].cmax + 1) ;
        pCmin = (*g->param_ptr)[g->res*x + y].cmin/((*g->param_ptr)[g->res*x + y].cmin + 1) ;
        D = (*f->flow_ptr)[f->res*x + y] ;
        Emin.x = 1.0 ; Emin.y = 0.0 ;
        Emax.x = 0.0 ; Emax.y = 1.0 ;
        Emin = vec_rotate(Emin,-(*g->param_ptr)[g->res*x + y].th) ;
        Emax = vec_rotate(Emax,-(*g->param_ptr)[g->res*x + y].th) ;

        (*q->flow_ptr)[q->res*x + y].x = Ukp.x + 
        pCmax*((D.x - Ukp.x)*Emax.x + (D.y - Ukp.y)*Emax.y)*Emax.x +
        pCmin*((D.x - Ukp.x)*Emin.x + (D.y - Ukp.y)*Emin.y)*Emin.x ;

        (*q->flow_ptr)[q->res*x + y].y = Ukp.y + 
        pCmax*((D.x - Ukp.x)*Emax.x + (D.y - Ukp.y)*Emax.y)*Emax.y +
        pCmin*((D.x - Ukp.x)*Emin.x + (D.y - Ukp.y)*Emin.y)*Emin.y ;
      }
    }
  }

  for(x = 0 ; x < f->sizy ; x++) {
    for(y = 0 ; y < f->sizx ; y++) {
      (*f->flow_ptr)[f->res*x + y] = (*q->flow_ptr)[q->res*x + y] ;
    }
  }
  dealloc_flow(p) ;
  dealloc_flow(q) ;
  free(p) ;
  free(q) ;
}

/* 
           NAME : compute_flow(im1h,im1q,im2h,im2q,flh,flq,sm,sp,iter)
   PARAMETER(S) : im1h,im1q : head and queue pointers on a 
                              pyramid of Laplacian images (image1);
                  im2h,im2q : head and queue pointers on a 
                              pyramid of Laplacian images (image2);
                    flh,flq : head and queue pointers on a list
                              of displacement-fields;
                         sm : relaxation option;
                         sp : sub-pixel accuracy option;
                       iter : number of iterations for relaxation.
 
        PURPOSE : Computes the sub-pixel accurate displacement 
                  field using overlapped projection.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 9 1990
*/

void compute_flow(im1h,im1q,im2h,im2q,flh,flq,sm,sp,iter)
qnode_ptr_t *im1h, *im1q, *im2h, *im2q, *flh, *flq ;
int sm, sp, iter ;

{ disp_vect_t best_match(), project(), U[S], V ; 
  qnode_ptr_t fl ; 
  int i, j, k, m, r, ip, jp, min ;
  float ssd[S], min_ssd ;

  if ((*flq)->back == (qnode_ptr_t)NULL) {   /* single level, no projection */
    alloc_param(*flq) ;
    for(i = 0 ; i < (*flq)->sizy ; i++) {
      for(j = 0 ; j < (*flq)->sizx ; j++) {
        (*(*flq)->flow_ptr)[(*flq)->res*i + j] = 
        best_match(*im1h,*im2h,&min_ssd,i,j,i,j) ;
        conf_measure(*im1h,*im2h,*flq,*flq,i,j,sp) ;
      }
    }
  }
  else {              /* displacement-vector projection from coarse to fine */
    fl = (*flq)->back ;
    r = (*flq)->res/fl->res ;

    for(i = 0 ; i < (*flq)->sizy ; i++) {
      ip = i - (i + 1) % r ;
      for(j = 0 ; j < (*flq)->sizx ; j++) {
        jp = j - (j + 1) % r ;
        k = 0 ;
        if (inside(ip/r,jp/r,fl->sizy,fl->sizx)) { 
          V = project((*fl->flow_ptr)[fl->res*(ip/r) + jp/r],r) ;
          U[k] = best_match(*im1h,*im2h,&min_ssd,i,j,(int)(V.x+i),(int)(V.y+j)) ;
          ssd[k] = min_ssd ;
          k++ ;
        }
        if (inside(ip/r + 1,jp/r,fl->sizy,fl->sizx)) { 
          V = project((*fl->flow_ptr)[fl->res*(ip/r + 1) + jp/r],r) ;
          U[k] = best_match(*im1h,*im2h,&min_ssd,i,j,(int)(V.x+i),(int)(V.y+j)) ;
          ssd[k] = min_ssd ;
          k++ ;
        }
        if (inside(ip/r,jp/r + 1,fl->sizy,fl->sizx)) {
          V = project((*fl->flow_ptr)[fl->res*(ip/r) + jp/r + 1],r) ;
          U[k] = best_match(*im1h,*im2h,&min_ssd,i,j,(int)(V.x+i),(int)(V.y+j)) ;
          ssd[k] = min_ssd ;
          k++ ;
        }
        if (inside(ip/r + 1,jp/r + 1,fl->sizy,fl->sizx)) { 
          V = project((*fl->flow_ptr)[fl->res*(ip/r + 1) + jp/r + 1],r) ;
          U[k] = best_match(*im1h,*im2h,&min_ssd,i,j,(int)(V.x+i),(int)(V.y+j)) ;
          ssd[k] = min_ssd ;
          k++ ;
        }
        min_ssd = MAXSSD ;    /* run minimum selection of the ssd measures */
        for(m = 0 ; m < k ; m++) {
          if (ssd[m] < min_ssd) {
            min_ssd = ssd[m] ;
            min = m ;
          }
        }
        (*(*flq)->flow_ptr)[(*flq)->res*i + j] = U[min] ;
      }
    }
    alloc_param(*flq) ;
    for(i = 0 ; i < (*flq)->sizy ; i++) {
      for(j = 0 ; j < (*flq)->sizx ; j++) {
        conf_measure(*im1h,*im2h,*flq,*flq,i,j,sp) ;
      }
    }
  }
  dealloc_lap(*im1h) ;                      /* freeing processed levels of */
  dealloc_lap(*im2h) ;             /* pyramid 1 and 2 and  parameter field */
  delete_node(*im1h,im1h,im1q) ;
  delete_node(*im2h,im2h,im2q) ;
  if (sm) { 
    relax(*flq,*flq,iter) ; 
  }
}

/* 
           NAME : cons_gauss(q)
   PARAMETER(S) : q : pointer on a node;
 
        PURPOSE : Constructs the Gaussian images
                  in the pyramid pointed by q.
                  The pointer q must point on the
                  node corresponding to the highest
                  resolution.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 6 1990

*/

void cons_gauss(q)
qnode_ptr_t q ;

{ extern kernel_t kera ;

  while (q != (qnode_ptr_t)NULL) {
    convolve(q,kera) ;
    if (q->back != (qnode_ptr_t)NULL) {
      sample(q,q->back) ;
    }
    dealloc_lap(q) ;
    q->lap_ptr = q->gauss_ptr ;
    q = q->back ;
  }
}

/* 
           NAME : cons_lap(q)
   PARAMETER(S) : q : pointer on a node;
 
        PURPOSE : Constructs the laplacian images
                  in the pyramid pointed by q.
                  The pointer q must point on the
                  node corresponding to the highest
                  resolution.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 6 1990

*/

void cons_lap(q)
qnode_ptr_t q ;

{ extern kernel_t kera ;

  while (q != (qnode_ptr_t)NULL) {
    convolve(q,kera) ;
    if (q->back != (qnode_ptr_t)NULL) {
      sample(q,q->back) ;
    }
    laplacian(q) ;
    dealloc_gauss(q) ;
    q->gauss_ptr = q->lap_ptr ;
    q = q->back ;
  }
}

/*
           NAME : prod_histo(f,c_fn,h_fn,nbin,incr,h_type,ttl_err,ttl_std,
                             dens,div)
   PARAMETER(S) :    f : flow field with parameters;
                  c_fn : correct velocities file name;
                  h_fn : histogram file name;
                  nbin : number of bins in histo 1;
                  incr : increment per bin in histo 1;
                h_type : type of histogram (cmin or cmax);
               ttl_err : error in flow field;
               ttl_std : standard deviation;
                  dens : dens.
		  div  : step size (usually 1)
 
        PURPOSE : produces an error histogram
                  h1: error vs confidence measure Cmin;
                  h2: h1 cumulated.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : February 20 1992
*/

void prod_histo(f,c_fn,h_fn,nbin,incr,h_type,ttl_err,ttl_std,dens,div)
qnode_ptr_t f ;
string c_fn, h_fn ;
float incr, *ttl_err, *ttl_std, *dens, div ;
int nbin, h_type ;

{ qnode_ptr_t create_node(), q ;
  disp_vect_t u ;
  histo_t histo1[N_BINS], c_histo1[N_BINS] ;
  float psi_error(), actual_x, actual_y, size_x, size_y, 
        offset_x, offset_y, max_qty, qty, err, avg_err, density, t, x, y ;
  int fdf, nbytes, i, j, k, index1, ttl_freq, abs_freq ;
  FILE *fdp ;

  q = create_node(0,f->res,f->sizx,f->sizy) ;
  alloc_flow(q) ;
  if ((fdf = open(c_fn,O_RDONLY)) <= 0) {
    error(7) ;
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
      if(x!=100.0 && y!=100.0) {
      (*q->flow_ptr)[q->res*i + j].y = y ;
      (*q->flow_ptr)[q->res*i + j].x = -x ;
      }
      else {
      (*q->flow_ptr)[q->res*i + j].y = y;
      (*q->flow_ptr)[q->res*i + j].x = -x;
      }
    }
  }
  close(fdf) ;

  max_qty = incr*(float)nbin ;
  for (i = 0 ; i < N_BINS ; i++) {
    histo1[i].avg = 0.0 ;
    histo1[i].std = 0.0 ;
    histo1[i].freq = 0 ;
    c_histo1[i].avg = 0.0 ;
    c_histo1[i].std = 0.0 ;
    c_histo1[i].freq = 0 ;
    abs_freq = 0 ;
    ttl_freq = 0 ;
    avg_err = 0.0 ;
    *ttl_err = 0.0 ;
    *ttl_std = 0.0 ;
  }
  for (i = KERNEL_Y + NRADIUS ; i < f->sizy - KERNEL_Y - NRADIUS ; i++) {
    for (j = KERNEL_X + NRADIUS ; j < f->sizx -  KERNEL_X - NRADIUS ; j++) {
      if (h_type == CMIN) {
        qty = (*f->param_ptr)[f->res*i + j].cmin ;
      }
      else {
        qty = (*f->param_ptr)[f->res*i + j].cmax ;
      }
      index1 = (int)((qty/max_qty)*(float)nbin) ;
      u = (*f->flow_ptr)[f->res*i+j] ;
      abs_freq++ ;
      if (u.x != -100.0 && u.y != 100.0) {
        err = psi_error((*f->flow_ptr)[f->res*i+j],
                        (*q->flow_ptr)[q->res*i+j],div) ;
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
  *ttl_err = avg_err/(float)ttl_freq ;
  avg_err = *ttl_err ;
  if (abs_freq != 0) {
    *dens = (float)ttl_freq/(float)abs_freq ;
  }
  else {
    *dens = 0.0 ;
  }
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
      if (h_type == CMIN) {
        qty = (*f->param_ptr)[f->res*i + j].cmin ;
      }
      else {
        qty = (*f->param_ptr)[f->res*i + j].cmax ;
      }
      index1 = (int)((qty/max_qty)*(float)nbin) ;
      u = (*f->flow_ptr)[f->res*i+j] ;
      if (u.x != -100.0 && u.y != 100.0) {
        err = psi_error((*f->flow_ptr)[f->res*i+j],(*q->flow_ptr)[q->res*i+j],div) ;
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
  fprintf(fdp,"%5.7f\n",max_qty - incr/2.0) ;
  for (i = 0 ; i < nbin ; i++) {
    t = max_qty/(float)nbin ;
    fprintf(fdp,"%5.7f %10.7f %10.7f %3.7f\n", t*(float)(i)+t/2.0, 
    histo1[i].avg, histo1[i].std, (float)histo1[i].freq/density) ;
  }
  fprintf(fdp,"\n\n\n") ;

  fprintf(fdp,"%3d\n",nbin) ;
  fprintf(fdp,"%5.7f\n",max_qty - incr/2.0) ;
  for (i = 0 ; i < nbin ; i++) {
    t = max_qty/(float)nbin ; 
    fprintf(fdp,"%5.7f %10.7f %10.7f %3.7f\n", t*(float)(i)+t/2.0, 
    c_histo1[i].avg, c_histo1[i].std, (float)c_histo1[i].freq/density) ;
  }
  fclose(fdp) ;
}

/* 
           NAME : anandan.c
   PARAMETER(S) : line parameters.
 
        PURPOSE : Computes optic flow between two image frames
                  using Anandan[87]'s algorithm.

         AUTHOR : Steven Beauchemin
             AT : University of Western Ontario
           DATE : April 9 1990
*/

int main(argc,argv)
int argc ;
char *argv[] ;

{ qnode_ptr_t create_node(), im1h, im1q, im2h, im2q, flh, flq, fl ;
  unsigned char hd[H] ;
  float trsh, incr, div, avg_err, std, density ;
  int n1, n2, level, sm, sp, sf, pyr, f_type, h_type, histo, nbin, iter, ssd,
  x[N_FRAME], y[N_FRAME], a, binary, row, col ;
  string in_path, out_path, v_fname, i_fname, c_fname, h_fname, f1, f2 ;

  valid_option(argc,argv,in_path,out_path,&n1,&n2,&level,&h_type,&histo,&sm,
               &sp,&sf,&pyr,&ssd,&f_type,&trsh,i_fname,v_fname,c_fname,h_fname,
               &nbin,&incr,&iter,&binary,&row,&col) ;
  printf("STAGE 1: Options Validated\n") ; fflush(stdout) ;
  init_list(&im1h,&im1q) ;
  init_list(&im2h,&im2q) ;
  init_list(&flh,&flq) ;
  init_kernel_a(&kera) ;
  init_kernel_b(&kerb) ;
  filenames(i_fname,n1,n2,f1,f2,&div) ;
  if (!binary) {
    raster_size(f1,&(x[0]),&(y[0]),&a) ;
    raster_size(f2,&(x[1]),&(y[1]),&a) ;
    valid_size(x,y) ;
  }
  else {
    x[0] = col ; x[1] = col ;
    y[0] = row ; y[1] = row ;
    binary_size(col,row,&a) ;
  }
  printf("STAGE 2: Initializations Completed\n") ; fflush(stdout) ;
  create_pyramid(&im1h,&im1q,level,a,x[0],y[0],2) ;
  create_pyramid(&im2h,&im2q,level,a,x[1],y[1],2) ;
  if (!binary) {
    pgetrast(f1,hd,im1q->gauss_ptr,x[0],y[0],a) ;
    pgetrast(f2,hd,im2q->gauss_ptr,x[1],y[1],a) ;
  }
  else {
    Bpgetrast(f1,im1q->gauss_ptr,col,row,a) ;
    Bpgetrast(f2,im2q->gauss_ptr,col,row,a) ;
  }
  printf("STAGE 3: Input Data Read\n") ; fflush(stdout) ;
  if (pyr == LAP) {
    cons_lap(im1q) ;
    cons_lap(im2q) ;
  }
  else {
    cons_gauss(im1q) ;
    cons_gauss(im2q) ;
  }
  printf("STAGE 4: Pyramids Created\n") ; fflush(stdout) ;
  if (ssd == 3) {
    init_beaudet3(&Ix,&Ixx,&Ixy) ;
  }
  else {
    if (ssd == 5) {
      init_beaudet5(&Ix,&Ixx,&Ixy) ;
    }
    else {
      if (ssd == 7) {
        init_beaudet7(&Ix,&Ixx,&Ixy) ;
      }
    }
  }
  level -= 1 ;
  while ((im1h != (qnode_ptr_t)NULL) && (im2h != (qnode_ptr_t)NULL)) {
    fl = create_node(level,im1h->res,im1h->sizx,im1h->sizy) ;
    insert_node(fl,&flh,&flq) ;
    alloc_flow(fl) ;
    compute_flow(&im1h,&im1q,&im2h,&im2q,&flh,&flq,sm,sp,iter) ;
    level-- ;
  }
  printf("STAGE 5: Flow Estimated\n") ; fflush(stdout) ;
  if (sf) {
    filter(flq,flq,f_type,trsh) ;
    printf("STAGE 6: Filtering Out Completed\n") ; fflush(stdout) ;
  }
  if (histo) {
    prod_histo(flq,c_fname,h_fname,nbin,incr,h_type,&avg_err,&std,&density,div) ;
    printf("STAGE 7: Histograms Produced\n") ; fflush(stdout) ;
  }
  dump_flow(&flh,&flq,v_fname,div) ;
  printf("STAGE 8: Flow Field Written to File\n") ; fflush(stdout) ;
  printf("STAGE 9: End of Program\n\n") ;
  if (histo) {
    printf("         Average Angular Error: %10.5f\n", avg_err) ;
    printf("            Standard Deviation: %10.5f\n", std) ;
    printf("                       Density: %10.5f\n", density*100.0) ;
    fflush(stdout) ;
  }
}

