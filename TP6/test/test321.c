
//------------------------------------------------------
// module  : Tp4-IFT2425-2.c
// author  : 
// date    : 
// version : 1.0
// language: C++
// note    :
//------------------------------------------------------
//  

//------------------------------------------------
// FICHIERS INCLUS -------------------------------
//------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <new>

/************************************************************************/
/* WINDOWS						          	*/
/************************************************************************/
#include <X11/Xutil.h>

Display   *display;
int	  screen_num;
int 	  depth;
Window	  root;
Visual*	  visual;
GC	  gc;

//------------------------------------------------
// DEFINITIONS -----------------------------------                       
//------------------------------------------------
#define CARRE(X) ((X)*(X))
#define CUBE(X) ((X)*(X)*(X))
#define OUTPUT_FILE "Tp4-Img-II"
#define VIEW_PGM    "xv" 

#define DEBUG 1
#define TROIS 3

//-Cst-Modele
#define X_1 0.0
#define Y_1 1.0
#define X_2 -1.0/sqrt(2.0)
#define Y_2 -1.0/2.0
#define X_3 +1.0/2*sqrt(2.0)
#define Y_3 -1.0/2.0
#define C 0.25  
#define R 0.1  
#define D 0.3

 
//-Cst-Runge-Kutta
#define H            0.1       
#define T_0          0.0                 
#define T_F          20.0      
#define NB_INTERV (T_F-T_0)/H
   
 //-Cst-Image                             
#define WIDTH  128     
#define HEIGHT 128                
#define MAX_X  4.0                
#define MAX_Y  4.0  
#define EVOL_GRAPH 3000
              
#define WHITE     255
#define GREYWHITE 230
#define GREY      200
#define GREYDARK  120
#define BLACK       0   

struct pos { 
  float x;
  float y; 
};

struct speed {
  float x;
  float y;
};

//#define float double
//------------------------------------------------
// GLOBAL CST ------------------------------------                       
//------------------------------------------------
float Xmin=0.0;
float Xmax=0.0;
float Ymin=0.0;
float Ymax=0.0; 

float xx_1=((WIDTH/MAX_X)*X_1)+(WIDTH/2);
float yy_1=(-(HEIGHT/MAX_Y)*Y_1)+(HEIGHT/2);
float xx_2=((WIDTH/MAX_X)*X_2)+(WIDTH/2);
float yy_2=(-(HEIGHT/MAX_Y)*Y_2)+(HEIGHT/2);
float xx_3=((WIDTH/MAX_X)*X_3)+(WIDTH/2);
float yy_3=(-(HEIGHT/MAX_Y)*Y_3)+(HEIGHT/2);

float X_1_INI;
float X_2_INI;
float X_3_INI;
float X_4_INI;

/************************************************************************/
/* OPEN_DISPLAY()							*/
/************************************************************************/
int open_display()
{
  if ((display=XOpenDisplay(NULL))==NULL)
   { printf("Connection impossible\n");
     return(-1); }

  else
   { screen_num=DefaultScreen(display);
     visual=DefaultVisual(display,screen_num);
     depth=DefaultDepth(display,screen_num);
     root=RootWindow(display,screen_num);
     return 0; }
}

/************************************************************************/
/* FABRIQUE_WINDOW()							*/
/* Cette fonction crée une fenetre X et l'affiche à l'écran.	        */
/************************************************************************/
Window fabrique_window(char *nom_fen,int x,int y,int width,int height,int zoom)
{
  Window                 win;
  XSizeHints      size_hints;
  XWMHints          wm_hints;
  XClassHint     class_hints;
  XTextProperty  windowName, iconName;

  char *name=nom_fen;

  if(zoom<0) { width/=-zoom; height/=-zoom; }
  if(zoom>0) { width*=zoom;  height*=zoom;  }

  win=XCreateSimpleWindow(display,root,x,y,width,height,1,0,255);

  size_hints.flags=PPosition|PSize|PMinSize;
  size_hints.min_width=width;
  size_hints.min_height=height;

  XStringListToTextProperty(&name,1,&windowName);
  XStringListToTextProperty(&name,1,&iconName);
  wm_hints.initial_state=NormalState;
  wm_hints.input=True;
  wm_hints.flags=StateHint|InputHint;
  class_hints.res_name=nom_fen;
  class_hints.res_class=nom_fen;

  XSetWMProperties(display,win,&windowName,&iconName,
                   NULL,0,&size_hints,&wm_hints,&class_hints);

  gc=XCreateGC(display,win,0,NULL);

  XSelectInput(display,win,ExposureMask|KeyPressMask|ButtonPressMask| 
               ButtonReleaseMask|ButtonMotionMask|PointerMotionHintMask| 
               StructureNotifyMask);

  XMapWindow(display,win);
  return(win);
}

/****************************************************************************/
/* CREE_XIMAGE()							    */
/* Crée une XImage à partir d'un tableau de float                          */
/* L'image peut subir un zoom.						    */
/****************************************************************************/
XImage* cree_Ximage(float** mat,int z,int length,int width)
{
  int lgth,wdth,lig,col,zoom_col,zoom_lig;
  float somme;
  unsigned char	 pix;
  unsigned char* dat;
  XImage* imageX;

  /*Zoom positiv*/
  /*------------*/
  if (z>0)
  {
   lgth=length*z;
   wdth=width*z;

   dat=(unsigned char*)malloc(lgth*(wdth*4)*sizeof(unsigned char));
   if (dat==NULL)
      { printf("Impossible d'allouer de la memoire.");
        exit(-1); }

  for(lig=0;lig<lgth;lig=lig+z) for(col=0;col<wdth;col=col+z)
   { 
    pix=(unsigned char)mat[lig/z][col/z];
    for(zoom_lig=0;zoom_lig<z;zoom_lig++) for(zoom_col=0;zoom_col<z;zoom_col++)
      { 
       dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+0)]=pix;
       dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+1)]=pix;
       dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+2)]=pix;
       dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+3)]=pix; 
       }
    }
  } /*--------------------------------------------------------*/

  /*Zoom negatifv*/
  /*------------*/
  else
  {
   z=-z;
   lgth=(length/z);
   wdth=(width/z);

   dat=(unsigned char*)malloc(lgth*(wdth*4)*sizeof(unsigned char));
   if (dat==NULL)
      { printf("Impossible d'allouer de la memoire.");
        exit(-1); }

  for(lig=0;lig<(lgth*z);lig=lig+z) for(col=0;col<(wdth*z);col=col+z)
   {  
    somme=0.0;
    for(zoom_lig=0;zoom_lig<z;zoom_lig++) for(zoom_col=0;zoom_col<z;zoom_col++)
     somme+=mat[lig+zoom_lig][col+zoom_col];
           
     somme/=(z*z);    
     dat[((lig/z)*wdth*4)+((4*(col/z))+0)]=(unsigned char)somme;
     dat[((lig/z)*wdth*4)+((4*(col/z))+1)]=(unsigned char)somme;
     dat[((lig/z)*wdth*4)+((4*(col/z))+2)]=(unsigned char)somme;
     dat[((lig/z)*wdth*4)+((4*(col/z))+3)]=(unsigned char)somme; 
   }
  } /*--------------------------------------------------------*/

  imageX=XCreateImage(display,visual,depth,ZPixmap,0,(char*)dat,wdth,lgth,16,wdth*4);
  return (imageX);
}

/****************************************************************************/
/* CREE_XIMAGECOUL()							    */
/* Crée une XImage à partir d'un tableau 3 d de float                       */
/* L'image peut subir un zoom.						    */
/****************************************************************************/
XImage* cree_XimageCoul(float*** matRVB,int z,int length,int width)
{
  int i;
  int lgth,wdth,lig,col,zoom_col,zoom_lig;
  float somme;
  float sum[3];
  unsigned char	 pixR,pixV,pixB,pixN;
  unsigned char* dat;
  XImage* imageX;

  /*Zoom positif*/
  /*------------*/
  if (z>0)
  {
   lgth=length*z;
   wdth=width*z;

   dat=(unsigned char*)malloc(lgth*(wdth*4)*sizeof(unsigned char));
   if (dat==NULL)
      { printf("Impossible d'allouer de la memoire.");
        exit(-1); }

  for(lig=0;lig<lgth;lig=lig+z) for(col=0;col<wdth;col=col+z)
   { 
    pixR=(unsigned char)matRVB[0][lig/z][col/z];
    pixV=(unsigned char)matRVB[1][lig/z][col/z];
    pixB=(unsigned char)matRVB[2][lig/z][col/z];
    somme=(1.0/3.0)*(pixR+pixV+pixB);
    pixN=(unsigned char)somme;

    for(zoom_lig=0;zoom_lig<z;zoom_lig++) for(zoom_col=0;zoom_col<z;zoom_col++)
      { 
       dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+0)]=pixB; 
       dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+1)]=pixV; 
       dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+2)]=pixR; 
       dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+3)]=0; 
       }
    }
  } /*--------------------------------------------------------*/

  /*Zoom negatif*/
  /*------------*/
  else
  {
   z=-z;
   lgth=(length/z);
   wdth=(width/z);

   dat=(unsigned char*)malloc(lgth*(wdth*4)*sizeof(unsigned char));
   if (dat==NULL)
      { printf("Impossible d'allouer de la memoire.");
        exit(-1); }

  for(lig=0;lig<(lgth*z);lig=lig+z) for(col=0;col<(wdth*z);col=col+z)
   {  
    sum[0]=sum[1]=sum[2]=0.0;
    
    for(i=0;i<3;i++)
    for(zoom_lig=0;zoom_lig<z;zoom_lig++) for(zoom_col=0;zoom_col<z;zoom_col++)
     sum[i]+=matRVB[i][lig+zoom_lig][col+zoom_col];
       
    for(i=0;i<3;i++)  sum[i]/=(z*z); 

     dat[((lig/z)*wdth*4)+((4*(col/z))+0)]=(unsigned char)sum[1];
     dat[((lig/z)*wdth*4)+((4*(col/z))+1)]=(unsigned char)sum[1];
     dat[((lig/z)*wdth*4)+((4*(col/z))+2)]=(unsigned char)sum[1];
     dat[((lig/z)*wdth*4)+((4*(col/z))+3)]=(unsigned char)sum[1]; 
   }
  } /*--------------------------------------------------------*/

  imageX=XCreateImage(display,visual,depth,ZPixmap,0,(char*)dat,wdth,lgth,16,wdth*4);
  return (imageX);
}

//------------------------------------------------
// FUNCTIONS -------------------------------------                       
//------------------------------------------------
//-------------------------//
//-- Matrice de Double ----//
//-------------------------//
//---------------------------------------------------------
// Alloue de la memoire pour une matrice 1d de float
//----------------------------------------------------------
float* dmatrix_allocate_1d(int hsize)
 {
  float* matrix;
  matrix=new float[hsize]; return matrix; }

//----------------------------------------------------------
// Alloue de la memoire pour une matrice 2d de float
//----------------------------------------------------------
float** dmatrix_allocate_2d(int vsize,int hsize)
 {
  float** matrix;
  float *imptr;

  matrix=new float*[vsize];
  imptr=new float[(hsize)*(vsize)];
  for(int i=0;i<vsize;i++,imptr+=hsize) matrix[i]=imptr;
  return matrix;
 }

//----------------------------------------------------------
// alloue de la memoire pour une matrice 3d de float
//----------------------------------------------------------
float*** dmatrix_allocate_3d(int dsize,int vsize,int hsize)
 {
  float*** matrix;

  matrix=new float**[dsize];

  for(int i=0;i<dsize;i++)
    matrix[i]=dmatrix_allocate_2d(vsize,hsize);
  return matrix;
 }

//----------------------------------------------------------
// Libere la memoire de la matrice 1d de float
//----------------------------------------------------------
void free_dmatrix_1d(float* pmat)
{ delete[] pmat; }

//----------------------------------------------------------
// Libere la memoire de la matrice 2d de float
//----------------------------------------------------------
void free_dmatrix_2d(float** pmat)
{ delete[] (pmat[0]);
  delete[] pmat;}

//----------------------------------------------------------
// libere la memoire de la matrice 3d de float
//----------------------------------------------------------
void free_dmatrix_3d(float*** pmat,int dsize)
{
 for(int i=0;i<dsize;i++)
  {
   delete[] (pmat[i][0]);
   delete[] (pmat[i]);
   }
 delete[] (pmat);
}

//----------------------------------------------------------
// Sauvegarde de l'image de nom <name> au format ppm        
//----------------------------------------------------------
void SaveImagePpm(char* Name,float*** matrvb,int wdth,int lgth)
 {
  int i,j;
  char buff[200];
  FILE* fuser;

  //extension
  strcpy(buff,Name);
  strcat(buff,".ppm");

  //ouverture fichier
  fuser=fopen(buff,"w");
    if (fuser==NULL) 
        { printf(" probleme dans la sauvegarde de %s",buff); 
          exit(-1); }

  //affichage
  printf("\n  Sauvegarde de %s au format %s",buff,".ppm");
  fflush(stdout);

  //sauvegarde de l'entete
  fprintf(fuser,"P6");
  fprintf(fuser,"\n# IMG Module");
  fprintf(fuser,"\n%d %d",lgth,wdth);
  fprintf(fuser,"\n255\n");

  //enregistrement
  for(i=0;i<wdth;i++) for(j=0;j<lgth;j++) 
    {
     fprintf(fuser,"%c",(char)matrvb[0][i][j]);
     fprintf(fuser,"%c",(char)matrvb[1][i][j]);
     fprintf(fuser,"%c",(char)matrvb[2][i][j]);
    }
       
  //fermeture fichier
   fclose(fuser); 
 }

//------------------------------------------------------------------------
// plot_point
//
// Affiche entre x dans [-MAX_X/2  MAX_X/2]
//               y dans [-MAX_Y/2  MAX_Y/2]                
//------------------------------------------------------------------------
void plot_point(float** MatPts,float** MatPict,int NbPts)
{
 int x_co,y_co;
 int i,j,k;

 //Init
 for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++)  MatPict[i][j]=GREYWHITE;

 for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++) 
   { if ((fabs(i-yy_1)+fabs(j-xx_1))<10) MatPict[i][j]=GREYDARK;
     if ((fabs(i-yy_2)+fabs(j-xx_2))<10) MatPict[i][j]=GREYDARK;
     if ((fabs(i-yy_3)+fabs(j-xx_3))<10) MatPict[i][j]=GREYDARK; }

 //Loop
 for(k=0;k<NbPts;k++)
    { x_co=(int)((WIDTH/MAX_X)*MatPts[k][0]);
      y_co=-(int)((HEIGHT/MAX_Y)*MatPts[k][1]);
      y_co+=(HEIGHT/2);
      x_co+=(WIDTH/2);
      if (DEBUG) printf("[%d::%d]",x_co,y_co); 
      if ((x_co<WIDTH)&&(y_co<HEIGHT)&&(x_co>0)&&(y_co>0)) 
	 MatPict[y_co][x_co]=BLACK; 
    }
}

//------------------------------------------------------------------------
// Fill_Pict
//------------------------------------------------------------------------
void Fill_Pict(float** MatPts,float** MatPict,int PtsNumber,int NbPts)
{
 int i,j;
 int x_co,y_co;
 int k,k_Init,k_End;

 //Init
 for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++) 
   { if (MatPict[i][j]!=GREYWHITE) MatPict[i][j]=GREY;
     if ((fabs(i-yy_1)+fabs(j-xx_1))<10) MatPict[i][j]=GREYDARK;
     if ((fabs(i-yy_2)+fabs(j-xx_2))<10) MatPict[i][j]=GREYDARK;
     if ((fabs(i-yy_3)+fabs(j-xx_3))<10) MatPict[i][j]=GREYDARK; }

 //Loop
 k_Init=PtsNumber;
 k_End=(k_Init+EVOL_GRAPH)%NbPts;
 for(k=k_Init;k<k_End;k++)
    { k=(k%NbPts);
      x_co=(int)((WIDTH/MAX_X)*MatPts[k][0]);
      y_co=-(int)((HEIGHT/MAX_Y)*MatPts[k][1]);
      y_co+=(HEIGHT/2);
      x_co+=(WIDTH/2);
      if ((x_co<WIDTH)&&(y_co<HEIGHT)&&(x_co>0)&&(y_co>0)) 
         MatPict[y_co][x_co]=BLACK; }
}

//------------------------------------------------
// Prototype de Fonctions ------------------------
//------------------------------------------------

float sumxy(float,float,float*,float*);
pos x_n_plus1(speed,pos);
float acceleration(float,float,float,float*,float*);

//------------------------------------------------
// FONCTIONS TPs----------------------------------                      
//------------------------------------------------


/*
void acceleration(accel acc, speed vt, pos p) {
  acc.x = -(R * vt.x - sumxy(p.x, p.y) + C * p.x);
  acc.y = -(R * vt.y - sumxy(p.y, p.x) + C * p.y);
}
*/

float acceleration(float s, float p1, float p2, float* arrx, float* arry) {
  return -(R * s - sumxy(p1,p2,arrx,arry) + C * p1);
}

/* calculate the sum for accel on one coordinates.
   for values in respect to y, just send sum(y,x) */
float sumxy(float x, float y, float* arrx, float* arry) {
  float sum1 = (arrx[0] - x) / CUBE(sqrt(CARRE(arrx[0] - x) + CARRE(arry[0] - y) + CARRE(D)));
  float sum2 = (arrx[1] - x) / CUBE(sqrt(CARRE(arrx[1] - x) + CARRE(arry[1] - y) + CARRE(D)));
  float sum3 = (arrx[2] - x) / CUBE(sqrt(CARRE(arrx[2] - x) + CARRE(arry[2] - y) + CARRE(D)));
  return sum1+sum2+sum3;
}

pos x_n_plus1(speed s, pos p) {
  int i;
  float x, y, v, k;
  float arrx[3] = {X_1, X_2, X_3};
  float arry[3] = {Y_1, Y_2, Y_3};
  for (i=0;i<2;i++) {
    if (i == 0) 
    {
      x = p.x;  y = p.y;  v = s.x;
    }
    else 
    {
      x = p.y;  y = p.x;  v = s.y;
      float tmp[3] = {X_1, X_2, X_3};
      arrx[0] = arry[0]; arrx[1] = arry[1]; arrx[2] = arry[2];
      arry[0] = tmp[0]; arry[1] = tmp[1]; arry[2] = tmp[2];
    }
    float k1 = H * v;
    float k2 = H * (v + (H/4.0) * 
               acceleration(v, x + k1/4.0, y, arrx, arry));
    float k3 = H * (v + (3.0*H/8.0) * 
               acceleration(v, x + 3.0*k1/32.0 + 9.0*k2/32.0, y, arrx, arry));
    float k4 = H * (v + (12.0*H/13.0) * 
               acceleration(v, x + 1932.0*k1/2197.0 - 7200.0*k2/2197.0 + 7296.0*k3/2197.0, y, arrx, arry));
    float k5 = H * (v + H * 
               acceleration(v, x + 439.0*k1/216.0 - 8.0*k2 + 3680.0*k3/513.0 
                                 - 845.0*k4/4104.0, y, arrx, arry));
    float k6 = H * (v + (H/2.0) * 
               acceleration(v, x - 8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 
                                 + 1859.0*k4/4104.0 - 11.0*k5/40.0, y, arrx, arry));

    k = 16.0*k1/135.0 
      + 6656.0*k3/12825.0 
      + 28561.0*k4/56430.0
      - 9.0*k5/50.0 
      + 2.0*k6/55.0;
    if (i == 0)
      p.x = p.x + k;
    else
      p.y = p.y + k;
  }
  return p;
}

float distance(float x, float x0, float y, float y0) {
  return fabs(x - x0) + fabs(y - y0);
}

int nearestMagnet(float x, float y, float* xArr, float* yArr) {
  int i = 0;
  float d = distance(x,xArr[0],y,yArr[0]);
  float smallD = d;
  float smallDIndex = i;
  for (i=1; i < 3; i++) {
    d = distance(x,xArr[i],y,yArr[i]);
    if (d <= smallD) {
      smallD = d;
      smallDIndex = i;
    }
  }
  return smallDIndex;
}
//----------------------------------------------------------
//----------------------------------------------------------
// MAIN  
//----------------------------------------------------------
//----------------------------------------------------------
int main (int argc, char **argv)
{
  int i,j,k;
  int flag_graph;
  int zoom;

  XEvent ev;
  Window win_ppicture;
  XImage *x_ppicture;
  char   nomfen_ppicture[100]; 
  char BufSystVisu[100];

  //>AllocMemory
  float*** MatPict=dmatrix_allocate_3d(TROIS,HEIGHT,WIDTH);
  float** MatPts=dmatrix_allocate_2d((int)(NB_INTERV),2);
  
  //>Init
  for(k=0;k<TROIS;k++) for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++) MatPict[k][i][j]=0;
  for(i=0;i<2;i++) for(j=0;j<(int)(NB_INTERV);j++) MatPts[i][j]=0.0;
  flag_graph=1;
  zoom=1;


  //---------------------------------------------------------------------
  //>Question 2 
  //---------------------------------------------------------------------  

  //Il faut travailler ici ...et dans > // FONCTIONS TPs

  //Un exemple ou la matrice de points MatPict est remplie
  //par une image couleur donné par l'équation d'en bas... et non pas par 
  //les bassins d'attractions

  //for(k=0;k<TROIS;k++) for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++) 
  //   {  MatPict[k][i][j]=(i+j*k*i)%255; }

  //Un exemple ou la matrice de points MatPict est remplie
  //par une image en niveaux de gris  donné par l'équation d'en bas... et non pas par 
  //la vitesse de convergence
  pos p;
  speed v;
  v.x = 0.0; v.y = 0.0;
  float color;
  int nMagnetI;
  float nMagnetD;
  int inRange;
  int cMagnetI;
  int converging = 0;

  float xarr[3] = {X_1, X_2, X_3};
  float yarr[3] = {Y_1, Y_2, Y_3};

  for (i=0; i < HEIGHT; i++) {
    for (j=0; j < WIDTH; j++) {
      p.x = (j / (float)HEIGHT) * MAX_X - MAX_X / 2.0;
      p.y = (1.0 - i /  (float)WIDTH) * MAX_Y - MAX_Y / 2.0;

      nMagnetI = nearestMagnet(p.x,p.y,xarr,yarr);
      nMagnetD = distance(p.x,xarr[nMagnetI],p.y,yarr[nMagnetI]);
      inRange = nMagnetD < 0.5;

      for(k=1;k<(int)(NB_INTERV);k++)
        { 
          p = x_n_plus1(v,p);

          v.x = v.x + H * acceleration(v.x, p.x, p.y, xarr, yarr);
          v.y = v.y + H * acceleration(v.y, p.y, p.x, yarr, xarr);
          cMagnetI = nearestMagnet(p.x,p.y,xarr,yarr);
          if (cMagnetI != nMagnetI) {
            nMagnetI = cMagnetI;
            converging = 0;
          }
          nMagnetD = distance(p.x,xarr[cMagnetI],p.y,yarr[cMagnetI]);
          if (nMagnetD >= 0.5) {
            converging = 0;
          } else {
            converging++;
          }
        }
        color = converging > 20 ? 255 - converging : 255;
        MatPict[0][i][j] = MatPict[1][i][j] = MatPict[2][i][j] = color;
    }
  }




/*
  for(k=0;k<TROIS;k++) for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++) 
     {  MatPict[0][i][j]=(i+j*k*i)%255; 
        MatPict[1][i][j]=(i+j*k*i)%255;
        MatPict[2][i][j]=(i+j*k*i)%255;  }
*/
 
   //--Fin Question 2-----------------------------------------------------


  //>Save&Visu de MatPict
  SaveImagePpm((char*)OUTPUT_FILE,MatPict,HEIGHT,WIDTH);
  

  //--------------------------- 

  //>Affiche Statistique
  printf("\n\n Stat:  Xmin=[%.2f] Xmax=[%.2f] Ymin=[%.2f] Ymax=[%.2f]\n",Xmin,Xmax,Ymin,Ymax);

 //--------------------------------------------------------------------------------
 //-------------- visu sous XWINDOW -----------------------------------------------
 //--------------------------------------------------------------------------------
 if (flag_graph)
 {
 //>Uuverture Session Graphique
 if (open_display()<0) printf(" Impossible d'ouvrir une session graphique");
 sprintf(nomfen_ppicture,"Évolution du Graphe");
 win_ppicture=fabrique_window(nomfen_ppicture,10,10,HEIGHT,WIDTH,zoom);
 x_ppicture=cree_XimageCoul(MatPict,zoom,HEIGHT,WIDTH);

 printf("\n\n Pour quitter,appuyer sur la barre d'espace");
 fflush(stdout);

  //boucle d'evenements
  for(;;)
     {
      XNextEvent(display,&ev);
       switch(ev.type)
        {
	 case Expose:   

         XPutImage(display,win_ppicture,gc,x_ppicture,0,0,0,0,x_ppicture->width,x_ppicture->height);  
         break;

         case KeyPress: 
         XDestroyImage(x_ppicture);

         XFreeGC(display,gc);
         XCloseDisplay(display);
         flag_graph=0;
         break;
         }
   if (!flag_graph) break;
   }
 } 
       
 //>Retour  
 printf("\n Fini... \n\n\n");
 return 0;
}
