/*
	Header File For Operating Xlib(TM)
	Original Program By Ken Tanaka
	Re-write Program By Yakuza Isono
	1989/12/18 Version 0.5
*/

void g_init();
void g_open_window(int,double,double,char*);
void g_open(char []);
void g_close(int);
void g_window(int,double,double,double,double);
void g_setFont(int,char *);
void g_string(int,double,double,char *);
void g_setColor(int,int);
void g_set_ColorPixcel(); 
double g_getXLength(int,int);
double g_getYLength(int,int);
void g_pset(int,double,double);
void g_line(int,double,double,double,double,int); 
void g_styledline(int,double,double,double,double,int,int);
void g_styledlines(int,double *,double *,int,int);
void g_box(int,double,double,double,double,int);
void g_boxfill(int,double,double,double,double);
void g_boxfill_w(int,double,double,double,double);
void g_boxfill2(int,double,double,double,double);
void g_circle(int,double,double,double);
void g_circlefill(int,double,double,double);
void g_clearWindow(int);
void g_make_pixpat(int,char *,int,int);
void g_set_pixpat(int);
char *get_host_name(char *);
void g_flush( void );
void g_polygonfill(int, double [][2], int );
void g_polygon(int, double [][2], int, int );
