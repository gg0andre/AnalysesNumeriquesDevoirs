//------------------------------------------------------
// module  : Tp2-IFT2425-1.c
// author  : Andre-Viet tran (20051185)
// date    : 23 avril 2020
// version : 1.0
// language: C++
// note    :
//------------------------------------------------------
//  

//------------------------------------------------
// FICHIERS INCLUS -------------------------------
//------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/                              
/*------------------------------------------------*/
#define CARRE(X) ((X)*(X))
#define DEBUG 0
#define NBITER_GS 1000

//-------------------------
//--- Windows -------------
//-------------------------
#include <X11/Xutil.h>

Display   *display;
int	  screen_num;
int 	  depth;
Window	  root;
Visual*	  visual;
GC	  gc;

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
/* Crée une XImage à partir d'un tableau de float                           */
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

//-------------------------//
//-- Matrice de Flottant --//
//-------------------------//
//---------------------------------------------------------
//  alloue de la memoire pour une matrice 1d de float
//----------------------------------------------------------
float* fmatrix_allocate_1d(int hsize)
 {
  float* matrix;
  matrix=(float*)malloc(sizeof(float)*hsize); return matrix; }

//----------------------------------------------------------
//  alloue de la memoire pour une matrice 2d de float
//----------------------------------------------------------
float** fmatrix_allocate_2d(int vsize,int hsize)
 {
  int i;
  float** matrix;
  float* imptr;

  matrix=(float**)malloc(sizeof(float*)*vsize);
  imptr=(float*)malloc(sizeof(float)*(hsize)*(vsize));
  for(i=0;i<vsize;i++,imptr+=hsize) matrix[i]=imptr;
  return matrix;
 }

//----------------------------------------------------------
// Allocation matrix 3d float
//----------------------------------------------------------
float*** fmatrix_allocate_3d(int dsize,int vsize,int hsize)
 {
  int i;
  float*** matrix;

  matrix=(float***)malloc(sizeof(float**)*dsize); 

  for(i=0;i<dsize;i++)
    matrix[i]=fmatrix_allocate_2d(vsize,hsize);
  return matrix;
 }

//----------------------------------------------------------
// Sauvegarde de l'image de nom <name> au format pgm                        
//----------------------------------------------------------                
void SaveImagePgm(char* bruit,char* name,float** mat,int lgth,int wdth)
{
 int i,j;
 char buff[300];
 FILE* fic;

  //--extension--
  strcpy(buff,bruit);
  strcat(buff,name);
  strcat(buff,".pgm");

  //--ouverture fichier--
  fic=fopen(buff,"wb");
    if (fic==NULL) 
        { printf("Probleme dans la sauvegarde de %s",buff); 
          exit(-1); }
  printf("\n Sauvegarde de %s au format pgm\n",buff);

  //--sauvegarde de l'entete--
  fprintf(fic,"P5");
  fprintf(fic,"\n# IMG Module");
  fprintf(fic,"\n%d %d",wdth,lgth);
  fprintf(fic,"\n255\n");

  //--enregistrement--
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
	fprintf(fic,"%c",(char)mat[i][j]);
   
  //--fermeture fichier--
   fclose(fic); 
}

//-------------------------//
//---- Fonction Pour TP ---//
//-------------------------//
//----------------------------------------------------------
// Affiche Axes                       
//----------------------------------------------------------                
void DisplayAxe(float** Graph2D,int length,int width)
{
 int i;
 for(i=0;i<length;i++) Graph2D[length/2][i]=0.0;
}

//----------------------------------------------------------
// Affichage dans Tab2D/Graph2D                 
//----------------------------------------------------------                
void PrintInTab2D(float** Tab2D,int lgth,int wdth,float x,float y,float x_inf,float x_sup,float y_inf,float y_sup)
{
 float pix_x,pix_y;
 float x_scale=(x_sup-x_inf)/wdth;
 float y_scale=(y_sup-y_inf)/lgth;

 pix_x=(int)((x-x_inf)/x_scale);
 pix_y=(int)((y-y_inf)/y_scale);

 //printf("[%.0f:::%.0f]",pix_y,pix_x);

 if ((pix_x>0)&&(pix_x<lgth)&&(pix_y>0)&&(pix_y<wdth)) Tab2D[(int)pix_y][(int)pix_x]=0.0;
}

//----------------------------------------------------------//
// Écrire une matrice n*m fans un fichier                   //
//----------------------------------------------------------//
void writeMatrixF(FILE* info, float** mat, int n, int m)
{   
 int i,j;

 if (info) fprintf(info,"    ");

 //Loop
 for(i=0;i<n;i++)
   { for(j=0;j<m;j++)
     if (info) fprintf(info,"%6.2f", mat[i][j]);
     if (info) fprintf(info, "\n    "); }
}

//----------------------------------------------------------//
// Écrire une matrice n*m                                   //
//----------------------------------------------------------//
void writeMatrix(float** mat, int n, int m)
{   
 int i,j;

 //Loop
 for(i=0;i<n;i++)
   { for(j=0;j<m;j++)
     printf("%6.2f", mat[i][j]);
     printf("\n"); }
}

//-------------------------//
//--- Vos Fonctions ici ---//
//-------------------------//
void GaussSeidel(float** b,float** x,int lgth,float** A){
  // int i = 0;
  // x[i][i] = b[i][i] - (A[i][i] *  x[i][i] + A[i+1][i] * x[i+1][i]);
  // for (i=1; i<lgth-1; i++) 
  //   x[i][i] = b[i][i] - (A[i][i] *  x[i][i] + A[i+1][i] * x[i+1][i] + A[i-1][i] * x[i-1][i]);

  // x[i][i] = b[i][i] - (A[i][i] *  x[i][i] + A[i-1][i] * x[i-1][i]);

  int Nelem = 20;
  int i=0, j=0;
 
  for (i=0; i < lgth; i++) {
    x[j][i] = (1/A[i][i]) * (b[j][i] - x[j+1][i] * A[i][i+1]);

    for (j=1; j < Nelem - 2 ; j++) {
      x[j][i] = (1/A[i][i]) * (b[j][i] - (x[j+1][i] * A[i][i+1] + x[j-1][i] * A[i][i-1]));
    }

    x[j][i] = (1/A[i][i]) * (b[j][i] - x[j-1][i] * A[i-1][i]);
  }
}

//----------------------------------------------------------
//----------------------------------------------------------
// PROGRAMME PRINCIPAL -------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------
int main(int argc,char** argv)
{
 int   i,j,k;
 
 //Pour Xwindow
 XEvent ev;
 Window win_ppicture1,win_ppicture2;
 XImage *x_ppicture1,*x_ppicture2; 
 char   nomfen_ppicture1[100],nomfen_ppicture2[100]; 
 int    length,width;
 int    flag_graph;
 int    zoom;

 length=width=128;
 flag_graph=1;
 zoom=1;

 //Variables/Cst du problème
 float *x;
 float **A;
 float **Vct_fx;
 float **Vct_ux; 

 //#define NDIS 5000
 #define NDIS 100
 #define NTEM 20

 //Allocation Memoire Vecteur d'Images pour Windows
 float*** Vct_ImgIn=fmatrix_allocate_3d(NTEM,length,width);
 float*** Vct_ImgOut=fmatrix_allocate_3d(NTEM,length,width);

 //Allocation Memoire Vecteurs/Matrices
 x=fmatrix_allocate_1d(NDIS-1);
 Vct_fx=fmatrix_allocate_2d(NTEM,(NDIS-1));
 Vct_ux=fmatrix_allocate_2d(NTEM,(NDIS-1));
 A=fmatrix_allocate_2d((NDIS-1),(NDIS-1));

 //Initialisation Vecteur x[i]
 for(i=0;i<(NDIS-1); i++)  x[i]=(float)(i+1)/(float)NDIS;

 //Initialisation Vecteur fx[t][i]
 for(j=0;j<NTEM;j++) for(i=0;i<(NDIS-1);i++) Vct_fx[j][i]=x[i]*sin(1.5*M_PI*x[i]*x[i]*(2.0-(float)j/NTEM));

 //Initialisation Vecteur ux[t][i]
 for(j=0;j<NTEM;j++) for(i=0;i<(NDIS-1);i++) Vct_ux[j][i]=0.0;

 //Remplissage Matrices 
 for(i=0;i<(NDIS-1);i++) for(j=0;j<(NDIS-1);j++) 
    { A[i][j]=0.0; 
      if (i==j) { A[i][j]=2.0;  }
      else if (i == j+1 || i == j-1) A[i][j]=-1.0; }

 //Affichage Matrice A 
 if (DEBUG) { printf("\n\n\nMATRICE A :\n");
          writeMatrix(A,(NDIS-1),(NDIS-1)); }


//--------------------------------------------------------------------------------
// (*) La matrice carrée A de taille NDIS*NDIS est définie et initialisée 
//     avec des valeurs correctes
//
// (*)  Le vecteur x_i=[x_1,x_2,...x_n-1] est définie et initialisé correctement
//
// (*)  le vecteur f(t,x) est définie et initialisé correctement :: Vct_fx[][]
//
// A vous de faire le fonctions:
//      
//        Question 1.    void GaussSeidel(float* b,float* x,int lgth,float** A)
//                       pour le systeme A.x=b  (resultat dans le vecteur x)
//
//        Question 2.    Ecrire ensuite le code qui permettra de 
//                       trouver les différents vecteurs solutions
//                       Vct_ux[k][i] pour les différents membres de droites Vct_fx[k][i]
//                                      k=[0,...,NTEM=10]  i=[0,...,NDIS=100]
// 
//--------------------------------------------------------------------------------
// PROGRAMME ---------------------------------------------------------------------
//--------------------------------------------------------------------------------
 
 //>votre prog. principal ici
 //---------------------------

  GaussSeidel(Vct_fx, Vct_ux, NDIS, A);

  for(j=0;j<NTEM;j++) for(i=0;i<(NDIS-1);i++) Vct_ux[j][i] = Vct_ux[j][i] * (CARRE(NDIS));
  


//--------------------------------------------------------------------------------
// 
// Les fonctions ci-dessous vont representer les différents vecteurs 
//     Vct_fx[k] dans l'image Vct_ImgIn[k]
//  et Vct_ux[k] dans l'image Vct_ImgOut[k]        
//
//  et representer la déformation temporelle comme une suite de séquence d'images
//-------------------------------------------------------------------------------- 

 //Vecteur fx dans Vct_ImgIn qui sera ensuite visualisé
 for(k=0;k<NTEM;k++) for(i=0;i<length;i++) for(j=0;j<width;j++) 
   { Vct_ImgIn[k][i][j]=220.0;
     DisplayAxe(Vct_ImgIn[k],length,width); }

 for(j=0;j<NTEM;j++) for(i=0;i<(NDIS-1);i++) 
    PrintInTab2D(Vct_ImgIn[j],length,width,x[i],Vct_fx[j][i],0.0,1.0,-1.0,1.0); 

 //Vecteur ux dans Vct_ImgOut qui sera ensuite visualisé
 for(k=0;k<NTEM;k++) for(i=0;i<length;i++) for(j=0;j<width;j++) 
   { Vct_ImgOut[k][i][j]=220.0;
     DisplayAxe(Vct_ImgOut[k],length,width); }

 for(j=0;j<NTEM;j++) for(i=0;i<(NDIS-1);i++) 
    PrintInTab2D(Vct_ImgOut[j],length,width,x[i],Vct_ux[j][i],0.0,1.0,-0.05,0.05); 


//------------------------------
// Sauvegarde
//------------------------------ 
//SaveImagePgm("","In",Vct_ImgIn[0],length,width);
//SaveImagePgm("","Out",Vct_ImgOut[0],length,width);

//--------------------------------------------------------------------------------
//---------------- Visu sous XWINDOW ---------------------------------------------
//--------------------------------------------------------------------------------
//
// Pour visualiser les differentes images de la  séquence, 
//        appuyer successivement sur Return
// Pour quitter appuyer sur une autre touche
//--------------------------------------------------------------------------------
 if (flag_graph)
 {
 //ouverture session graphique
 if (open_display()<0) printf(" Impossible d'ouvrir une session graphique");
 sprintf(nomfen_ppicture1,"f(x)","");
 sprintf(nomfen_ppicture2,"u(x)","");
 win_ppicture1=fabrique_window(nomfen_ppicture1,10,10,width,length,zoom);
 win_ppicture2=fabrique_window(nomfen_ppicture2,10+width+10,10,width,length,zoom);
 
 printf("\n\n Pour quitter, appuyer sur la barre d'espace");
 printf("\n Pour visualiser les differentes images de la  séquence, appuyer successivement sur Return\n");
 fflush(stdout);

 //Event Loop
  k=0;
  for(;;)
     { x_ppicture1=cree_Ximage(Vct_ImgIn[k],zoom,length,width);
       x_ppicture2=cree_Ximage(Vct_ImgOut[k],zoom,length,width);

       XNextEvent(display,&ev);
       switch(ev.type)
        {
	 case Expose:   
         XPutImage(display,win_ppicture1,gc,x_ppicture1,0,0,0,0,x_ppicture1->width,x_ppicture1->height); 
         XPutImage(display,win_ppicture2,gc,x_ppicture2,0,0,0,0,x_ppicture2->width,x_ppicture2->height);
         break;

         case KeyPress:
         if (XLookupKeysym(&ev.xkey,0)==XK_Return)  
           { k++; k=k%NTEM;
             printf("[%d]",k); fflush(stdout);
             usleep(500000); 
             x_ppicture1=cree_Ximage(Vct_ImgIn[k],zoom,length,width);
             x_ppicture2=cree_Ximage(Vct_ImgOut[k],zoom,length,width);
             XPutImage(display,win_ppicture1,gc,x_ppicture1,0,0,0,0,x_ppicture1->width,x_ppicture1->height); 
             XPutImage(display,win_ppicture2,gc,x_ppicture2,0,0,0,0,x_ppicture2->width,x_ppicture2->height); }

         else 
           { XDestroyImage(x_ppicture1);
             XDestroyImage(x_ppicture2);
             XFreeGC(display,gc);
             XCloseDisplay(display);
             flag_graph=0;
             break; }
         }
   if (!flag_graph) break;
   }
 } 
       
 //retour sans probleme 
 printf("\n Fini... \n\n\n");
 return 0;
 }
 



