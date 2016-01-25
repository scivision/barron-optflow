#include <fcntl.h>
#include <stdio.h>
#include <math.h>

#define SPACE 10

FILE *fo ;

int fi, cut=0, insubsamp=0, insub=1, subsamp=0, sub=1, caption =1,
	latex = 0, suboffset = 0 ;

int numx, numy, totx, toty, offx, offy, x1, x2, yl, y2 ;
float width, getval() ;
float multiply= 0.5, numinch = 7.5, pixwidth = 2250.0 ;
float scale ;
float offset ; 


main(argc,argv)
int argc;
char **argv;

{ int argcount = 1;

  while (argcount < argc && argv[argcount][0] == '-') {
    switch (argv[argcount][1]) {
      case 'c': cut = 1;           
                break;
      case 'm': sscanf(argv[argcount]+2,"%f",&multiply);  
                break;
      case 'l': pixwidth = 810.0;
                numinch = 2.7;
                caption = 0;
                latex = 1;
                break;
      case 'f': pixwidth = 1800.0;
                numinch = 6.0;
                caption = 0;
                latex = 1;
                break;
      case 'S': insubsamp = 1;
                sscanf(argv[argcount]+2,"%d",&insub);  
                break;
      case 's': subsamp = 1;
                sscanf(argv[argcount]+2,"%d",&sub);  
                break;
      case 'o': sscanf(argv[argcount]+2,"%d",&suboffset);  
                break;
      case 'n': caption = 0;
                break;
      default : printf("Invalid option\n");
	        usage(argv[0]);
    }		
    argcount++;
  }
  if ((argc-argcount)!=2) usage(argv[0]);
  initialize(argv[argc-2],argv[argc-1]);
  execute();
  close(fo);
}

usage(s)
char *s;

{

  printf("Usage:  %s [-c] [-l] [-f] [-sN] [-mN.N] [-n] <data-filename>  <postscript-outfile>\n",s);
  printf("\t-c    : cut  - prompts for corners of box to limit flow field.\n");
  printf("\t-l    : output for latex 2.7in.\n"); 
  printf("\t-f    : output for latex 6.0in.\n"); 
  printf("\t-sN   : subsample - N is integer - uses every Nth point in x & y.\n");
  printf("\t-mN.N : multiply (scale) - multiplies all vectors by N.N\n");
  printf("\t-n    : no caption printed.\n");
  printf("\t-oN   : choose which pixel to use when subsampling\n"); 
  exit(-1);
}

initialize(in,out)
char *in, *out ;
/* Read raster files into internal 3-D array. */

{ char fname[100] ;
  int  i, j, k, fp, ok ;

  sub = sub*insub;
  if ((fo = fopen(out,"w")) ==NULL ) printf("Error creating file %s.\n",out);

  if ( (fi = open(in,O_RDONLY)) <=0 ) {
    printf("File %s does not exist.\n",in);
    exit(0);
  }
  totx = (int)getval(-1,-1);
  toty = (int)getval(-1,-2);
  numx = (int)getval(-1,-1);
  numy = (int)getval(-1,-2);
  offx = (int)getval(-1,-1);
  offy = (int)getval(-1,-2);

  if (cut) {
    printf("Offset : in X: %d  in Y: %d\n",offx,offy);
    ok = 1;
    while (ok) { 
      printf("Low X Value (%d - %d) :",offx,totx-offx);
      scanf("%d",&x1);
      if (x1 >= offx && x1 <= totx - offx) ok = 0;
    }
    ok = 1;
    while (ok) { 
      printf("High X Value (%d - %d):",x1,totx-offx);
      scanf("%d",&x2);
      if (x2 > x1 && x2 <= totx - offx) ok = 0;	
    }
    ok = 1;
    while (ok) { 
      printf("Low Y Value (%d - %d) :",offy,toty-offy);
      scanf("%d",&yl);
      if (yl >= offy && yl <= toty - offy) 
      ok = 0;	
    }
    ok = 1;
    while (ok) { 
      printf("High Y Value (%d - %d):",yl,toty-offy);
      scanf("%d",&y2);
      if (y2 > yl && y2 <= toty - offy) ok = 0;	
    }
  } 
  scale = numinch*72.0/((float)(totx+2*sub+1)*SPACE);
  width = (float)2.0*(totx+2*sub+1)*SPACE/pixwidth;

  fprintf(fo,"%%! Input file name: %s\n",in) ;
  fprintf(fo,"%%! subsampled by %d scaled by %f\n",sub,multiply); 
  fprintf(fo,"%%! image size is %d in X and %d in Y\n",totx,toty); 
  fprintf(fo,"%%! offset in X is %d and %d in Y\n",offx,offy) ;
  if (cut) {
    fprintf(fo,"%%! flow printed is from %d to %d in X and from %d to %d in Y\n", x1,x2,yl,y2) ;
  }
  else {
    fprintf(fo,"%%! flow printed is from %d to %d in X and from %d to %d in Y\n", offx,totx-offx,offy,toty-offy) ;
  }
  fprintf(fo,"%%!\n") ;
  fprintf(fo,"%5.3f %5.3f scale\n", scale,scale);		
  if (latex) offset = 0.0; else offset = 30.0/scale;
}


execute()

{ int i, j, ii, jj, x, y, iang, leftx, lefty ;
  float U, V, angle, arrow, mag ;
	
  leftx = totx - numx - offx ;
  lefty = toty - numy - offy;

  fprintf(fo,"newpath\n");
  fprintf(fo,"  %d %d  moveto\n",(int)offset,(int)offset);
  fprintf(fo,"  %d %d rlineto\n",(totx+2*sub+1)*SPACE,0);
  fprintf(fo,"  %d %d rlineto\n",0,(toty+2*sub+1)*SPACE);
  fprintf(fo,"  %d %d rlineto\n",-(totx+2*sub+1)*SPACE,0);
  fprintf(fo,"  closepath\n");
  fprintf(fo,"  %5.2f setlinewidth\n",width);
  fprintf(fo,"stroke\n\n");
  if (caption) {
    char name[100],exp[100];
    printf("Enter Name for caption:\n");
    gets(name);
    printf("Enter Experiment for caption:\n");
    gets(exp);	
    fprintf(fo,"/Times-Roman findfont\n"); fprintf(fo,"  %6.2f scalefont\n  setfont\n",1.0/scale*12.0);
    fprintf(fo,"  %6.2f %6.2f  moveto\n",1.0/scale*72,1.0/scale*72*10.5);
    fprintf(fo,"  (%s --- %s) show\n",name,exp);
    fprintf(fo,"  %6.2f %6.2f  moveto\n",1.0/scale*72,1.0/scale*72*10.2);
    fprintf(fo,"  (Image size %d x %d    Subsampled by %d    scaled by %7.3f) show\n\n",totx,toty,sub,multiply);
  }
  arrow = 2*width;
  for (ii=0;ii<numy;ii++)
    for (jj=0;jj<numx;jj++) {
      if (insubsamp) {
        i = ii*insub;
        j = jj*insub;
      } 
      else {
        i = ii,j = jj;
      }
      U = getval(i,j);
      V = getval(i,j);

      if (U == 100.0 && V == 100.0) continue; 
      if (cut && (i+offy<yl || i+yl>y2 || j+offx<x1 || j+x1>x2)) continue;
      if (subsamp && ((i+offx-suboffset)%sub != 0 || (j+offy-suboffset)%sub != 0)) continue;
      if (cut) {
        x = (j+x1+sub+1)*SPACE;
        y = (toty-(i+yl)+sub)*SPACE ;
      } 
      else {
        x = (j+offx+sub+1)*SPACE;
        y = (toty-(i+offy)+sub)*SPACE;  
      }
      if (U == 0.0) {
        if (V>0.0) {
          angle = 90.0; 
        }
        else {
         angle = -90.0;
        }
      }
      else {
        angle = atan(V/U);  
        if (U<0.0) {
          angle += M_PI;
        }
        angle = angle*360/(2*M_PI);
      }
      mag = (sqrt(U*U + V*V)*(SPACE*multiply));
      x = (int)offset+x;
      y = (int)offset+y;

      fprintf(fo,"newpath\n  %d  %d translate\n",x,y);
      if (mag>0.001) fprintf(fo,"  %5.2f rotate\n",angle);
      fprintf(fo,"  0 0  moveto\n");
      if (mag<=0.001) fprintf(fo,"  %5.2f 0 rlineto\n",width);
      else {
        fprintf(fo,"  %5.2f 0 rlineto\n",mag); 
        fprintf(fo,"  %5.2f %5.2f rlineto\n",-arrow,arrow);
        fprintf(fo,"  %5.2f %5.2f rlineto\n",arrow,-arrow);
        fprintf(fo,"  %5.2f %5.2f rlineto\n",-arrow,-arrow);
      }

      fprintf(fo,"  %5.2f  setlinewidth\nstroke\n",width);
      if (mag>0.001) fprintf(fo,"  %5.2f rotate\n",-angle);
      fprintf(fo,"  %d  %d translate\n\n",-x,-y);
    }
    if (!latex) fprintf(fo,"showpage\n",-x,-y);
	
}


float getval(i,j)

{ float v ;

  if (read(fi,&v,sizeof(float))!= sizeof(float)) {
     printf("Read Error: %d %d\n",i,j);
  }
  return(v);
}


flength()

{ return(lseek(fi,0,3)); }
