/*
  Header File For Operating Xlib(TM)
  Original Program By Ken Tanaka
  Re-write Program By Yakuza Isono
  Re-Re-writed(multi-window version) Program By Ken Tanaka
  1/27/1990 Version 0.5
  program added around g_pset function By H.Kimura
  2/10/1993 Version 0.6
  */

/* Include Files--------------------------------------------- */

#include                <X11/Xlib.h>
#include                <math.h>
#include		<stdio.h>
#include		<stdlib.h>
#include		<string.h>
/*-----bitmap file--------------------*/ 



#define Max_Window_Number 20

#define SCREEN_X(W, X)  ((int)(((X) - WX1[(W)]) * FACTX[(W)] + VX1[(W)]))
#define SCREEN_Y(W, Y)  ((int)(((Y) - WY1[(W)]) * FACTY[(W)] + VY1[(W)]))

#define NUMBER_OF_LINESTYLES    5

#define LineDotted		1
#define LineDashed		2
#define LineDashDot		3
#define LineDashDotDotted	4

#define gray1_width 16
#define gray1_height 16
static char gray1_bits[] = {
  (char)0xff, (char)0xff, (char)0x01, (char)0x00, (char)0x01, (char)0x80, (char)0x39, (char)0x5e, (char)0x11, (char)0x20, (char)0x01, (char)0x10,
   (char)0x23, (char)0x7e, (char)0x01, (char)0x01, (char)0xa1, (char)0x1f, (char)0xfb, (char)0x5e, (char)0x05, (char)0x5c, (char)0x0f, (char)0x78,
  (char)0x19, (char)0x19, (char)0x39, (char)0x1d, (char)0x05, (char)0x00, (char)0xfd, (char)0x04}; 

static char dotted_list[10] = {5, 5};
static char dashed_list[10] = {15, 5};
static char dashdot_list[10] = {15, 3, 5, 3};
static char dashdotdotted_list[10] = {15, 3, 3, 3, 3, 3};

/* isono yori */
extern char **environ;

/* Define Constants--------------------------------------------- */

static	Window                  windowId[Max_Window_Number],root;
static	Display                 *display;
static	XSetWindowAttributes    xswa[Max_Window_Number];
static	GC                      gc[Max_Window_Number];
static  GC                      white_gc[Max_Window_Number];
static  GC                      black_gc[Max_Window_Number];
static	XEvent                  xevent[Max_Window_Number];
static	XGCValues               gcvalues[Max_Window_Number];
static  Font                    font[Max_Window_Number];
static	int                     screen;
static	int                     xx,yy;
static	double                  VX1[Max_Window_Number],VY1[Max_Window_Number];
static  double                  VX2[Max_Window_Number],VY2[Max_Window_Number]; 
static	double                  WX1[Max_Window_Number],WY1[Max_Window_Number];
static  double                  WX2[Max_Window_Number],WY2[Max_Window_Number];
static	double                  FACTX[Max_Window_Number],FACTY[Max_Window_Number];
static	double                  xx1[Max_Window_Number],yy1[Max_Window_Number];
static  double                  xx2[Max_Window_Number],yy2[Max_Window_Number];

static  Pixmap                  gray1;
static Colormap cmap;
static XColor c0,black,white,red,blue,green,yellow,magenta,cyan,orange,brown,pink,gray,salmon,violet,aquamarine,khaki,plum,sienna,thistle,palegreen,seagreen,yellowgreen,navy,wheat,gold,darkturquoise,deepskyblue,lightskyblue,powderblue,paleturquoise;
static int ColorNum[40];

/* Library Definition-------------------------------------------- */
void g_open(char dname[])
{
  display = XOpenDisplay(dname);
  screen  = DefaultScreen(display);
  root    = RootWindow(display,screen);
}
/*------------------------------------------------------------------*/
void g_open_window(int window_number,double x_width,double y_width,char* window_name)
{
  XEvent		ev;
  
  xswa[window_number].background_pixel = WhitePixel(display,screen);
  xswa[window_number].border_pixel     = BlackPixel(display,screen);
  xswa[window_number].backing_store    = Always;
  
  VX1[window_number] = 0.0;
  VY1[window_number] = 0.0;
  VX2[window_number] = x_width;
  VY2[window_number] = y_width;


 
  windowId[window_number] =
    XCreateWindow(display,root,0,0,(int)x_width,(int)y_width,5,CopyFromParent,
		  InputOutput,CopyFromParent,
		  CWBackPixel | CWBorderPixel | CWBackingStore,
		  &xswa[window_number]); 
  XStoreName(display,windowId[window_number] ,window_name); 
  gcvalues[window_number].function   = GXcopy;
  gcvalues[window_number].foreground = BlackPixel(display,screen);
  gcvalues[window_number].background = WhitePixel(display,screen);
  gc[window_number] = XCreateGC(display,windowId[window_number],
				GCFunction | GCForeground | GCBackground,
				&gcvalues[window_number]);
  white_gc[window_number] = XCreateGC(display,windowId[window_number],
                                      0,NULL);


/*  XMapWindow(display,windowId[window_number]); */
  XMapRaised(display, windowId[window_number]);
  XSync(display, 0);
  XSelectInput(display, windowId[window_number], ExposureMask);

}

/*---------------------------------------------------------------*/
void g_close(int window_number)
{
  XDestroyWindow(display, windowId[window_number]);
}
/*----------------------------------------------------------------*/
void g_window(int window_number,double x1,double y1,double x2,double y2)
{
  if(x1 < x2) {
    WX1[window_number] = x1;
    WX2[window_number] = x2;
  } else {
    WX1[window_number] = x2;
    WX2[window_number] = x1;
  }
  if(y1 < y2) {
    WY1[window_number] = y1;
    WY2[window_number] = y2;
  } else {
    WY1[window_number] = y2;
    WY2[window_number] = y1;
  }
  FACTX[window_number] = (VX2[window_number] - VX1[window_number])
    / (WX2[window_number] - WX1[window_number]);
  FACTY[window_number] = (VY2[window_number] - VY1[window_number])
    / (WY2[window_number] - WY1[window_number]);
}




/*-----------------------------------------------------------*/
void g_setFont(int window_number,char *fontname)
{
  font[window_number] = XLoadFont(display, fontname);
  XSetFont(display, gc[window_number], font[window_number]);
}
/*-----------------------------------------------------------*/
void g_string(int window_number,double x,double y,char *str)
{
  int xx, yy;

  xx = SCREEN_X(window_number, x);
  yy = SCREEN_Y(window_number, y);  
  XDrawString(display, windowId[window_number], gc[window_number],
	      xx, yy, str, strlen(str));
}

/*-----------------------------------------------------------*/
void g_setColor(int window_number,int color_num)
{
  XSetForeground(display,gc[window_number],ColorNum[color_num]);
}

void g_set_ColorPixcel()
{
  cmap=DefaultColormap(display,0);
  XAllocNamedColor(display,cmap,"white",&white,&c0);     ColorNum[0]=white.pixel;
  XAllocNamedColor(display,cmap,"black",&black,&c0);    ColorNum[1]=black.pixel;
  XAllocNamedColor(display,cmap,"red",&red,&c0);     ColorNum[2]=red.pixel;
  XAllocNamedColor(display,cmap,"blue",&blue,&c0);    ColorNum[3]=blue.pixel;
  XAllocNamedColor(display,cmap,"green",&green,&c0);   ColorNum[4]=green.pixel;
  XAllocNamedColor(display,cmap,"yellow",&yellow,&c0);  ColorNum[5]=yellow.pixel;
  XAllocNamedColor(display,cmap,"magenta",&magenta,&c0); ColorNum[6]=magenta.pixel;
  XAllocNamedColor(display,cmap,"cyan",&cyan,&c0);    ColorNum[7]=cyan.pixel;
  XAllocNamedColor(display,cmap,"orange",&orange,&c0);  ColorNum[8]=orange.pixel;
  XAllocNamedColor(display,cmap,"brown",&brown,&c0); ColorNum[9]=brown.pixel;
  XAllocNamedColor(display,cmap,"pink",&pink,&c0);    ColorNum[10]=pink.pixel;
  XAllocNamedColor(display,cmap,"gray",&gray,&c0);    ColorNum[11]=gray.pixel;
  XAllocNamedColor(display,cmap,"salmon",&salmon,&c0);    ColorNum[12]=salmon.pixel;
  XAllocNamedColor(display,cmap,"violet",&violet,&c0);    ColorNum[13]=violet.pixel;
  XAllocNamedColor(display,cmap,"aquamarine",&aquamarine,&c0);    ColorNum[14]=aquamarine.pixel;
  XAllocNamedColor(display,cmap,"khaki",&khaki,&c0);    ColorNum[15]=khaki.pixel;
 XAllocNamedColor(display,cmap,"plum",&plum,&c0);    ColorNum[16]=plum.pixel;
 XAllocNamedColor(display,cmap,"sienna",&sienna,&c0);    ColorNum[17]=sienna.pixel;
 XAllocNamedColor(display,cmap,"thistle",&thistle,&c0);    ColorNum[18]=thistle.pixel;
  XAllocNamedColor(display,cmap,"pale green",&palegreen,&c0);    ColorNum[19]=palegreen.pixel;
  XAllocNamedColor(display,cmap,"sea green",&seagreen,&c0);    ColorNum[20]=seagreen.pixel;
  XAllocNamedColor(display,cmap,"yellow green",&yellowgreen,&c0);    ColorNum[21]=yellowgreen.pixel;
  XAllocNamedColor(display,cmap,"navy",&navy,&c0);    ColorNum[22]=navy.pixel;
  XAllocNamedColor(display,cmap,"wheat",&wheat,&c0);    ColorNum[23]=wheat.pixel;
  XAllocNamedColor(display,cmap,"gold",&gold,&c0);    ColorNum[24]=gold.pixel;
  XAllocNamedColor(display,cmap,"darkturquoise",&darkturquoise,&c0);    ColorNum[25]=darkturquoise.pixel;
  XAllocNamedColor(display,cmap,"deepskyblue",&deepskyblue,&c0);    ColorNum[26]=deepskyblue.pixel;
  XAllocNamedColor(display,cmap,"lightskyblue",&lightskyblue,&c0);    ColorNum[27]=lightskyblue.pixel;
  XAllocNamedColor(display,cmap,"powderblue",&powderblue,&c0);    ColorNum[28]=powderblue.pixel;
  XAllocNamedColor(display,cmap,"paleturquoise",&paleturquoise,&c0);    ColorNum[29]=paleturquoise.pixel;
}


/*-----------------------------------------------------------*/
double g_getXLength(int window_number,int pixels)
{
  return((double)pixels / FACTX[window_number]);
}
 
/*-----------------------------------------------------------*/
double g_getYLength(int window_number,int pixels)
{
  return((double)pixels / FACTY[window_number]);
}
 
/*-----------------------------------------------------------*/
void g_pset(int window_number,double x,double y)
{
  int xx,yy;
  
  xx = SCREEN_X(window_number, x);
  yy = SCREEN_Y(window_number, y);
  XDrawPoint(display,windowId[window_number],gc[window_number],xx,yy);
  /* program added by H.Kimura 1993.2.10 */
}

/*-------------------------------------------------------------*/
void g_line(int window_number,double x1,double y1,double x2,double y2,int wide)
{
  int	xx1,yy1,xx2,yy2;
 
  xx1 = SCREEN_X(window_number, x1);
  yy1 = SCREEN_Y(window_number, y1);
  xx2 = SCREEN_X(window_number, x2);
  yy2 = SCREEN_Y(window_number, y2);
  
  XSetLineAttributes(display, gc[window_number], wide, LineSolid, CapButt, JoinMiter); 
  XDrawLine(display,windowId[window_number],gc[window_number],xx1,yy1,xx2,yy2); 
}

/*------------------------------------------------------------------*/
void g_styledline(int window_number,double x1,double y1,double x2,double y2,int style,int wide)
{
  int	xx1, yy1, xx2, yy2, tmp;

  xx1 = SCREEN_X(window_number, x1);
  yy1 = SCREEN_Y(window_number, y1);
  xx2 = SCREEN_X(window_number, x2);
  yy2 = SCREEN_Y(window_number, y2);

  tmp =(style == LineSolid) ? LineSolid : LineOnOffDash;

  XSetLineAttributes(display, gc[window_number], wide, tmp, CapButt, JoinMiter);

  switch (style) {
    case LineDotted:
    XSetDashes(display, gc[window_number], 0, dotted_list, 2);
    break;
    case LineDashed:
    XSetDashes(display, gc[window_number], 0, dashed_list, 2);
    break;
    case LineDashDot:
    XSetDashes(display, gc[window_number], 0, dashdot_list, 4);
    break;
    case LineDashDotDotted:
    XSetDashes(display, gc[window_number], 0, dashdotdotted_list, 6);
    break;
    default:
    break;
  }

  XDrawLine(display,windowId[window_number],gc[window_number],xx1,yy1,xx2,yy2);

}

/*----------------------------------------------------------*/
void g_styledlines(int window_number,double *xpts,double *ypts,int n,int style)
{
  XPoint *pts;
  int	 tmp;
  
  pts = (XPoint *)malloc(sizeof(XPoint) * n);
  for (tmp = 0; tmp < n; tmp++) {
    pts[tmp].x = (short)SCREEN_X(window_number, xpts[tmp]);
    pts[tmp].y = (short)SCREEN_Y(window_number, ypts[tmp]);
  }

  tmp =(style == LineSolid) ? LineSolid : LineOnOffDash;
  XSetLineAttributes(display, gc[window_number], 1, tmp, CapButt, JoinMiter);

  switch (style) {
    case LineDotted:
    XSetDashes(display, gc[window_number], 0, dotted_list, 2);
    break;
    case LineDashed:
    XSetDashes(display, gc[window_number], 0, dashed_list, 2);
    break;
    case LineDashDot:
    XSetDashes(display, gc[window_number], 0, dashdot_list, 4);
    break;
    case LineDashDotDotted:
    XSetDashes(display, gc[window_number], 0, dashdotdotted_list, 6);
    break;
    default:
    break;
  }

  XDrawLines(display,windowId[window_number],gc[window_number],
	     pts, n, CoordModeOrigin);
  
  free(pts);
}

/*-----------------------------------------------------*/
void g_box(int window_number,double x1,double y1,double	x2,double y2,int wide)
{
  int 	xx1,yy1,xx2,yy2;
  unsigned int dxx,dyy;
  double X1,X2,Y1,Y2;

  X1=x1; Y1=y1; X2=x2; Y2=y2;
  
  xx1 = SCREEN_X(window_number, X1);
  yy1 = SCREEN_Y(window_number, Y1);
  xx2 = SCREEN_X(window_number, X2);
  yy2 = SCREEN_Y(window_number, Y2);

  dxx = (unsigned int)(abs(xx2 - xx1));
  dyy = (unsigned int)(abs(yy2 - yy1));

    
  XSetLineAttributes(display, gc[window_number], wide, LineSolid, CapButt, JoinMiter); /* wide */

  XDrawRectangle(display,windowId[window_number],gc[window_number],xx1,yy1,dxx,dyy);
}

/*--------------------------------------------------------*/
void g_boxfill(int window_number,double	x1,double y1,double x2,double y2)
{
  int 	xx1,yy1,xx2,yy2;
  unsigned int dxx,dyy;
  double X1,X2,Y1,Y2;

  X1=x1; Y1=y1; X2=x2; Y2=y2;
  
  xx1 = SCREEN_X(window_number, X1);
  yy1 = SCREEN_Y(window_number, Y1);
  xx2 = SCREEN_X(window_number, X2);
  yy2 = SCREEN_Y(window_number, Y2);

  dxx = (unsigned int)(abs(xx2 - xx1));
  dyy = (unsigned int)(abs(yy2 - yy1));

 
  XFillRectangle(display,windowId[window_number],gc[window_number],xx1,yy1,dxx,dyy);
}



/*-----------------------------------------------------------------*/
void g_boxfill_w(int window_number,double x1,double y1,double x2,double y2)
{
  int 	xx1,yy1,xx2,yy2;
  unsigned int dxx,dyy;
  double X1,X2,Y1,Y2;

  X1=x1; Y1=y1; X2=x2; Y2=y2;
  
  xx1 = SCREEN_X(window_number, X1);
  yy1 = SCREEN_Y(window_number, Y1);
  xx2 = SCREEN_X(window_number, X2);
  yy2 = SCREEN_Y(window_number, Y2);

  dxx = (unsigned int)(abs(xx2 - xx1));
  dyy = (unsigned int)(abs(yy2 - yy1));

  XFillRectangle(display,windowId[window_number],white_gc[window_number],xx1,yy1,dxx,dyy);
}


/*-----------------------------------------------------------------*/
void g_boxfill2(int window_number,double x1,double y1,double x2,double y2)
{
  int 	xx1,yy1,xx2,yy2;
  unsigned int dxx,dyy;
  double X1,X2,Y1,Y2;

  X1=x1; Y1=y1; X2=x2; Y2=y2;
  
  xx1 = SCREEN_X(window_number, X1);
  yy1 = SCREEN_Y(window_number, Y1);
  xx2 = SCREEN_X(window_number, X2);
  yy2 = SCREEN_Y(window_number, Y2);

  dxx = (unsigned int)(abs(xx2 - xx1));
  dyy = (unsigned int)(abs(yy2 - yy1));

  XSetTile(display,gc[window_number],gray1);

  XFillRectangle(display,windowId[window_number],gc[window_number],xx1,yy1,dxx,dyy);
}


/*------------------------------------------------------------------*/
void g_circle(int window_number,double x,double y,double r)
{
  int		xx,yy;
  unsigned	rx,ry;
  
  xx = SCREEN_X(window_number, x);
  yy = SCREEN_Y(window_number, y);
  rx = (unsigned)fabs(2. * r * FACTX[window_number]);
  ry = (unsigned)fabs(2. * r * FACTY[window_number]);
  XDrawArc(display,windowId[window_number],gc[window_number],xx,yy,rx,ry,0,23040);
}
/*-----------------------------------------------------------------*/
void g_circlefill(int window_number,double x,double y,double r)
{
  int		xx,yy;
  unsigned	rx,ry;
  
  xx = SCREEN_X(window_number, x);
  yy = SCREEN_Y(window_number, y);
  rx = (unsigned)fabs(2. * r * FACTX[window_number]);
  ry = (unsigned)fabs(2. * r * FACTY[window_number]);
  XFillArc(display,windowId[window_number],gc[window_number],xx,yy,rx,ry,0,23040);
  
}
/*-----------------------------------------------------------------*/
void g_clearWindow(int window_number)
{
  XClearWindow(display,windowId[window_number]);
}


/*-----------------------------------------*/ 
void g_make_pixpat(int window_number,char *bits,int width,int height )
{
  gray1 = XCreatePixmapFromBitmapData(display,windowId[window_number],bits,
  width,height,BlackPixel(display,screen),WhitePixel(display,screen),
  DefaultDepth(display,screen));

}
/*-----------------------------------------*/
void g_set_pixpat(int window_number)
{
  g_make_pixpat(window_number,gray1_bits,gray1_width,gray1_height );

}
/*------------------------------------------*/
char *get_host_name(char *host)
{
	int  i,flag;
	char *ptr;

	flag = 1;
	i = 0;
	while((flag == 1) || (*environ[i] == '\0')){
		if(strncmp("DISPRAY=",environ[i],8)){
			flag = 0;
			ptr = environ[i];
		}
		i++;
	}
        for(i = 0;i < 8;i++){
		ptr++;
	}
	strcpy(host,ptr);
	return(host);
}

void g_flush( void )
{
  XFlush(display);
}

/*-----------------------------------------------------*/





void g_polygonfill(int window_number, double pointsD[][2], int npoints )
{
  XPoint points[npoints];
  for( int i = 0; i < npoints; ++i ){
    points[i].x = SCREEN_X(window_number, pointsD[i][0] );
    points[i].y = SCREEN_Y(window_number, pointsD[i][1] );
  }
 
  XFillPolygon(display,windowId[window_number],gc[window_number],points,npoints,Nonconvex,CoordModeOrigin);
}




void g_polygon(int window_number, double pointsD[][2], int npoints, int wide )
{
  int	xx1,yy1,xx2,yy2;

  for( int i = 0; i < npoints-1; ++i ){
    xx1 = SCREEN_X(window_number, pointsD[i][0]);
    yy1 = SCREEN_Y(window_number, pointsD[i][1]);
    xx2 = SCREEN_X(window_number, pointsD[i+1][0]);
    yy2 = SCREEN_Y(window_number, pointsD[i+1][1]);
  
    XSetLineAttributes(display, gc[window_number], wide, LineSolid, CapButt, JoinMiter); 
    XDrawLine(display,windowId[window_number],gc[window_number],xx1,yy1,xx2,yy2); 
  }

  xx1 = SCREEN_X(window_number, pointsD[npoints-1][0]);
  yy1 = SCREEN_Y(window_number, pointsD[npoints-1][1]);
  xx2 = SCREEN_X(window_number, pointsD[0][0]);
  yy2 = SCREEN_Y(window_number, pointsD[0][1]);
  
  XSetLineAttributes(display, gc[window_number], wide, LineSolid, CapButt, JoinMiter); 
  XDrawLine(display,windowId[window_number],gc[window_number],xx1,yy1,xx2,yy2); 
}









