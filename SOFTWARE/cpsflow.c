/************************************************************
   flow.c
   Travis Burkitt
   May 9 , 1989
	-- reads file containing velocities,
		creates postscript file to draw flow field.
************************************************************/
#include        <fcntl.h>
#include        <stdio.h>
#include        <math.h>

#define SPACE 10

FILE *fo;
int 	fi,cut=0,insubsamp=0,insub=1,  subsamp=0,sub=1,caption =1,
	showpage = 1,suboffset = 0,border=1;
int numx,numy,totx,toty,offx,offy,x1,x2,yl,y2;
float width,getval();
float multiply= 0.5,numinch = 7.5,pixwidth = 2250.0;
float scale;
float offset; 
float borderx,bordery;

main(argc,argv)
int argc;
char **argv;
{
	int argcount = 1;

	while (argcount < argc && argv[argcount][0] == '-') {
		switch (argv[argcount][1]) {
			case 'c': cut = 1;           
			break;
			case 'm': sscanf(argv[argcount]+2,"%f",&multiply);  
			break;
			case 'l': pixwidth = 1800.0;   /* latex */
				  numinch = 6.0;
				  caption = 0;
				  showpage = 0;
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
			case 'b': border = 0;
			break;
			case 'w': pixwidth = 1800.0;
					  numinch = 6.0;
					  caption = -1; 
			break;
			case 'f': caption = -1;     /* use filename as caption */
			break;
			default: fprintf(stderr,"Invalid option\n");
				 usage(argv[0]);
		}		
		argcount++;
	}
	if ((argc-argcount)!=1) 
		usage(argv[0]);
	initialize(argc,argv);
	execute(argv[argc-1]); close(fo);
}



usage(s)
char *s;
{
fprintf(stderr,"Usage:  %s [-c] [-l] [-sN] [-mN.N] [-n] <data-filename>\n",s);
fprintf(stderr,"\t-c    == cut  - prompts for corners of box to limit flow field.\n");
fprintf(stderr,"\t-sN   == subsample - N is integer - uses every Nth point in x & y.\n");
fprintf(stderr,"\t-mN.N == multiply (scale) - multiplies all vectors by N.N\n");
fprintf(stderr,"\t-n    == no caption printed.\n");
fprintf(stderr,"\t-f    == use filename as caption.\n");
fprintf(stderr,"\t-oN   == choose which pixel to use when subsampling\n"); 
fprintf(stderr,"\t-l    == output for latex.\n"); 
fprintf(stderr,"\t\t - no caption\n");
fprintf(stderr,"\t\t - no showpage command inserted\n");
fprintf(stderr,"\t\t - flow field smaller to fit text margins\n");
fprintf(stderr,"\t-w    == like -l option but with show page (test latex)\n");
fprintf(stderr,"\tPostscript goes to standard output.\n");
exit(1);
}



initialize(argc,argv)
int argc;
char **argv;
/************************************************************
   Read raster files into internal 3-D array.
************************************************************/
{

        char fname[100];
        int  i,j,k,fp,ok;


	sub = sub*insub;
       if ( (fi = open(argv[argc-1],O_RDONLY)) <=0 ) {
        	fprintf(stderr,"File %s does not exist.\n",argv[argc-1]);
        	exit(0);
		}

		totx = (int)getval(-1,-1);
		toty = (int)getval(-1,-2);
		numx = (int)getval(-1,-1);
		numy = (int)getval(-1,-2);
		offx = (int)getval(-1,-1);
		offy = (int)getval(-1,-2);
		fflush(stdout);
		
		if (!border) {
			totx = numx;
			toty = numy;
			offx = offy = 0;
		}

	if (cut) {
		ok = 1;
		while (ok) { 
			fprintf(stderr,"Low X Value (0 - %d) [%d - %d]:",
						totx,offx,totx-offx);
		  	scanf("%d",&x1);
			if (x1<=totx) 
				ok = 0;	
		}
		ok = 1;
		while (ok) { 
			fprintf(stderr,"High X Value (%d - %d):",x1,totx);
		  	scanf("%d",&x2);
			if (x2>x1 && x2<=numx) 
				ok = 0;	
		}
		ok = 1;
		while (ok) { 
			fprintf(stderr,"Low Y Value (0 - %d) [%d - %d]:",
								toty,offy,toty-offy);
		  	scanf("%d",&yl);
			if (yl<=toty) 
				ok = 0;	
		}
		ok = 1;
		while (ok) { 
			fprintf(stderr,"High Y Value (%d - %d):",yl,toty);
		  	scanf("%d",&y2);
			if (y2>yl && y2<=numy) 
				ok = 0;	
		}
		scale = numinch*72.0/((float)(x2-x1+2*sub+1)*SPACE);
   		width = (float)2.0*(x2-x1+2*sub+1)*SPACE/pixwidth;
	} else {
	        scale = numinch*72.0/((float)(totx+2*sub+1)*SPACE);
   	        width = (float)2.0*(totx+2*sub+1)*SPACE/pixwidth;
	}
	printf("%%!\n"); 

	printf("%%%%Created by %s -- written by Travis Burkitt\n",argv[0]);
	printf("%%%%Command used: ");
	for (i = 0; i<argc; i++)
		printf("%s ",argv[i]);
	printf("\n");

	printf("%5.3f %5.3f scale\n", scale,scale);		
	printf("2 setmiterlimit\n\n");

	borderx = (totx+2*sub+1)*SPACE;   /* corner of border */
	bordery = (toty+2*sub+1)*SPACE;

	if ( !showpage ) {        /* latex output */
		offset = 0;       /* latex will position, so leave no margin */
/**		borderx -= 30.0/scale;
		bordery -= 30.0/scale;  **/
	} else
		offset = 30.0/scale;
fflush(stdout);
}


execute(infile)
char *infile;
{
	int i,j,ii,jj,x,y,iang,leftx,lefty;
	float U,V,angle,arrow,mag,temp;
	
	leftx = totx - numx - offx;
	lefty = toty - numy - offy;

	printf("newpath\n");
	printf("  %d %d  moveto\n",(int)offset,(int)offset);
	printf("  %d %d rlineto\n",(int)borderx,0);
	printf("  %d %d rlineto\n",0,(int)bordery);
	printf("  %d %d rlineto\n",-(int)borderx,0);
	printf("  closepath\n");
	printf("  %5.2f setlinewidth\n",width);
	printf("stroke\n\n");

	printf("/line\n{ /transx exch def\n  /transy exch def\n  /mag exch def\n");
	printf("  /angle exch def\n  /arrow exch def\n  /width exch def\n");
	printf("\n  transx transy translate\n  angle rotate\n  0 0 moveto\n");
	printf("  mag 0 rlineto\n  arrow neg arrow rlineto\n");
	printf("  arrow arrow neg rlineto\n  arrow neg arrow neg rlineto\n");
	printf("  width setlinewidth\n  stroke\n  angle neg rotate\n");
	printf("  transx neg transy neg translate\n} def\n\n\n");

	if (caption) {
		char name[100],exp[100];
		if (caption == -1) {
			printf("/Times-Roman findfont\n"); 
			printf("  %6.2f scalefont\n  setfont\n",1.0/scale*12.0);
			printf("  %6.2f %6.2f  moveto\n",1.0/scale*72,1.0/scale*72*10.5);
			printf("  (Filename : %s) show\n",infile);
			printf("  %6.2f %6.2f  moveto\n",1.0/scale*72,1.0/scale*72*10.2);
			printf("  (Image size %d x %d    Subsampled by %d    scaled by %7.3f) show\n\n",totx,toty,sub,multiply);
		} else {
			fprintf(stderr,"Enter Name for caption:\n");
			gets(name);
			fprintf(stderr,"Enter Experiment for caption:\n");
			gets(exp);	
			printf("/Times-Roman findfont\n");
			printf("  %6.2f scalefont\n  setfont\n",1.0/scale*12.0);
			printf("  %6.2f %6.2f  moveto\n",1.0/scale*72,1.0/scale*72*10.5);
			printf("  (%s --- %s) show\n",name,exp);
			printf("  %6.2f %6.2f  moveto\n",1.0/scale*72,1.0/scale*72*10.2);
			printf("  (Image size %d x %d    Subsampled by %d    scaled by %7.3f) show\n\n",totx,toty,sub,multiply);
		}
	}
	fflush(stdout);
	arrow = 2*width;
	jj = (int) getval(1,1);
	while(jj != 1000000)
		{
		ii = (int) getval(1,1);
		if (insubsamp) 
			{
			i = ii*insub;
			j = jj*insub;
			} else 
			{
			i = ii;
			j = jj;
			}

		while((V=getval(1,1))!=1000000.0)
		{
		U = getval(ii,jj);
		temp = U;
		U = V;
		V = temp;

		if (cut && (i+offy<yl || i+offy>y2 || j+offx<x1 || j+offx> x2))
			continue;
		if (subsamp && ((i+offx-suboffset)%sub != 0 || 
			       (j+offy-suboffset)%sub != 0))
			continue;
		else printf("\%\% U=%f V=%f i=%d j=%d\n",U,V,i,j);
		if (cut) 
			{
		    	x = (j-x1+offx+sub+1)*SPACE;
		    	y = (i+offy-yl+sub);
			} 
			else 
			{
		    	x = (j+sub+1)*SPACE;
		    	y = (toty-(i)+sub) * SPACE;  
			}
		if (U == 0.0)
		if (V>0.0)
			angle = 90.0;
		else   angle = -90.0;
	        else {
		     angle = atan(V/U);  
		     if (U<0.0) 	
			angle += M_PI;
		     angle = angle*360/(2*M_PI);
		     }
		mag = (sqrt(U*U + V*V)*(SPACE*multiply));

		x = (int)offset+x;
		y = (int)offset+y;


		if (mag>0.01) 
		   printf("%5.2f %5.2f %5.2f %5.2f %d %d line\n",
					width,arrow,angle,mag,y,x);
		else {
		     if (mag < width) 
			mag = width;
		     printf("%5.2f 0 0 %5.2f %d %d line\n",width,mag,y,x);
		     fflush(stdout);
		     }
		}
	jj = (int) getval(1,1);
	}
	if(showpage) printf("showpage\n");
		
}



float getval(i,j)
{
	
	float v;
	if (read(fi,&v,sizeof(float))!= sizeof(float))
        	fprintf(stderr,"At end of data\n"),v=1000000.0;
	return(v);
}



flength()
{
	return(lseek(fi,0,3));
}
