/*
 * mt_dice.c -
 *
 * $Id: mt_dice.c,v 1.0 2013/07/30 11:36:34 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2013/07/30
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <string.h>

#include <vtmalloc.h>

#include <mt_dice.h>
#include <vt_inrimage.h>

static int _verbose_ = 0;






const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}









int MT_ComputeDice3D( vt_image imageIn,
                             vt_image imageExt,
                             char *filename )
{
  char *proc = "MT_ComputeDice3D";
  int *inHist, *extHist;
  int **tab;
  double **Dice;

  int nonEmptyIn  = 0;
  int nonEmptyExt = 0;
  int *indicesIn  = NULL;
  int *indicesExt = NULL;

  int i,j;

  MT_ComputeTab3D(imageIn, imageExt, &tab, &inHist, &extHist, &indicesIn, &indicesExt, &nonEmptyIn, &nonEmptyExt);

  Dice = vtmalloc( nonEmptyExt*sizeof(double*), "Dice", proc );
  for(j=0;j<nonEmptyExt;j++)
  {
    Dice[j] = vtmalloc( nonEmptyIn*sizeof(double), "Dice[j]", proc );
    for(i=0;i<nonEmptyIn;i++)
    {
      Dice[j][i] = 2*((double)tab[j][i])/
          ((double)extHist[indicesExt[j]]+(double)inHist[indicesIn[i]]);
    }
  }

  /* Generation du fichier */

  const char *extension = get_filename_ext(filename);

  /* fprintf(stdout, "Filename extension is %s\n", extension); */

  if ( strcmp ( extension, "txt" ) == 0 || strcmp ( extension, "" ) == 0 ) {  /*  Output format is a text file */
      FILE *fichier = fopen (filename, "w" );
      if (fichier == NULL)
      {
        perror (filename);
      }
      else
      {
        /* Dice */
        fprintf(fichier, "%-10s", "Dice");
        for (j=0;j<nonEmptyExt;j++)
        {
          fprintf(fichier, "%7d", indicesExt[j]);
        }
        fprintf(fichier, "\n\n\n");
        for (i=0;i<nonEmptyIn;i++)
        {
          fprintf(fichier, "%-10d", indicesIn[i]);
          for (j=0;j<nonEmptyExt;j++)
          {
            fprintf(fichier, "%7.3f", Dice[j][i]);
          }
          fprintf(fichier, "\n\n");
        }

        /* Histogrammes */
        /*  HistIn */
        fprintf(fichier, "\n%-10s\n\n", "HistIn");
        for(i=0;i<nonEmptyIn;i++)
        {
          fprintf(fichier, "%10d\n", inHist[indicesIn[i]]);
        }
        /*  HistExt */
        fprintf(fichier, "\n%-10s\n\n", "HistExt");
        for(i=0;i<nonEmptyExt;i++)
        {
          fprintf(fichier, "%10d\n", extHist[indicesExt[i]]);
        }

      }
  }
  else {  /*  Output format is an image */
      vt_image imgDice;
      float ***array;
      VT_Image( &imgDice );
      VT_InitImage(&imgDice, filename, nonEmptyIn, nonEmptyExt, 1, FLOAT);
      VT_AllocImage(&imgDice);
      array = (float***) imgDice.array;
      for(j=0;j<nonEmptyExt;j++)
      for(i=0;i<nonEmptyIn;i++)
      {
          array[0][j][i]=(float)Dice[j][i];
      }
      VT_WriteInrimage(&imgDice);
      VT_FreeImage(&imgDice);
  }



  vtfree(inHist);
  vtfree(extHist);
  for(i=0;i<nonEmptyExt;i++) {
    vtfree(tab[i]);
    tab[i]=NULL;
    vtfree(Dice[i]);
    Dice[i]=NULL;
  }
  vtfree(tab);
  tab=NULL;
  vtfree(Dice);
  Dice=NULL;
  inHist=NULL;
  extHist=NULL;
  return (1);
}

int MT_ComputeTab3D(vt_image imageIn,
                    vt_image imageExt,
                    int ***Table,
                    int **inH,
                    int **extH,
                    int **indIn,
                    int **indExt,
                    int *sizeIn,
                    int *sizeExt
                    )
{
    char *proc = "MT_ComputeTab3D";
    int dimxIN=0,dimyIN=0,dimzIN=0;
    int dimxEXT=0,dimyEXT=0,dimzEXT=0;
    int x,y,z;
    int inMax=0, extMax=0;
    int *inHist=NULL, *extHist=NULL;
    int *indicesIn, *indicesExt;
    int **tab, **table;
    int nonEmptyIn  = 0;
    int nonEmptyExt = 0;

    unsigned char ***inu8;
    unsigned char ***extu8;

    unsigned short int ***in;
    unsigned short int ***ext;

    switch (imageIn.type)
    {
      case SCHAR:
      case UCHAR:
        inu8= (unsigned char ***)imageIn.array;
        dimxIN=imageIn.dim.x;
        dimyIN=imageIn.dim.y;
        dimzIN=imageIn.dim.z;
        in = vtmalloc( dimzIN*sizeof(unsigned short int **), "in", proc );
        for (z=0;z<dimzIN;z++)
        {
          in[z]=vtmalloc( dimyIN*sizeof(unsigned short int *), "in[z]", proc );
          for (y=0;y<dimyIN;y++)
          {
            in[z][y]=vtmalloc( dimxIN*sizeof(unsigned short int), "in[z][y]", proc );
            for(x=0;x<dimxIN;x++)
              in[z][y][x]=(unsigned short int)(inu8[z][y][x]);
          }
        }
        break;
      case SSHORT:
      case USHORT:
        in = (unsigned short int***)imageIn.array;
        break;
      case TYPE_UNKNOWN:
      default:
        fprintf(stderr,"Unknown type or not supported format of imageIn");
        return(0);
        break;
    }

    switch (imageExt.type)
    {
      case SCHAR:
      case UCHAR:
        extu8= (unsigned char ***)imageExt.array;
        dimxEXT=imageExt.dim.x;
        dimyEXT=imageExt.dim.y;
        dimzEXT=imageExt.dim.z;
        ext = vtmalloc( dimzEXT*sizeof(unsigned short int **), "ext", proc );
        for (z=0;z<dimzEXT;z++)
        {
          ext[z]=vtmalloc( dimyEXT*sizeof(unsigned short int *), "ext[z]", proc);
          for (y=0;y<dimyEXT;y++)
          {
            ext[z][y]=vtmalloc( dimxEXT*sizeof(unsigned short int), "ext[z][y]", proc );
            for(x=0;x<dimxEXT;x++)
              ext[z][y][x]=(unsigned short int)(extu8[z][y][x]);
          }
        }
        break;
      case SSHORT:
      case USHORT:
        ext = (unsigned short int***)imageExt.array;
        break;
      case TYPE_UNKNOWN:
      default:
        fprintf(stderr,"Unknown type or not supported format of imageIn");
        if (imageIn.type == UCHAR)
        {
          for (z=0;z<dimzIN;z++)
          {
            for (y=0;y<dimyIN;y++)
            {
              vtfree(in[z][y]);
              in[z][y]=NULL;
            }
            vtfree(in[z]);
            in[z]=NULL;
          }
          vtfree(in);
          in=NULL;
        }
        return(0);
        break;
    }

    int i,j;

    for (z=0;z<(int)imageIn.dim.z;z++)
    for (y=0;y<(int)imageIn.dim.y;y++)
    for (x=0;x<(int)imageIn.dim.x;x++)
    {
      if (inMax < (int)in[z][y][x])
        inMax = (int)in[z][y][x];
      if (extMax < (int)ext[z][y][x])
        extMax = (int)ext[z][y][x];
    }

    /*--allocation + initialisation tableaux qui serviront au calcul des intersections--*/
    tab =  vtmalloc((extMax+1)*sizeof(int*), "tab", proc );
    for (i=0;i<=extMax;i++)
      tab[i] =  vtmalloc((inMax+1)*sizeof(int), "tab[i]", proc );
    inHist =  vtmalloc((inMax+1)*sizeof(int), "inHist", proc );
    extHist =  vtmalloc((extMax+1)*sizeof(int), "extHist", proc );

    for (i=0;i<=inMax;i++) inHist[i]=0;
    for (i=0;i<=extMax;i++) extHist[i]=0;

    for (i=0;i<=inMax;i++)
    for (j=0;j<=extMax;j++)
      tab[j][i] = 0;

    for (z=0;z<(int)imageIn.dim.z;z++)
    for (y=0;y<(int)imageIn.dim.y;y++)
    for (x=0;x<(int)imageIn.dim.x;x++)
    {
      inHist[(int)in[z][y][x]]  += 1;
      extHist[(int)ext[z][y][x]]+= 1;

      tab[(int)ext[z][y][x]][(int)in[z][y][x]] += 1;
    }


    for (i=0;i<=inMax;i++)
      if (inHist[i]!=0)
        ++nonEmptyIn;

    indicesIn =  vtmalloc( nonEmptyIn*sizeof(int), "indicesIn", proc );
    nonEmptyIn=0;
    for (i=0;i<=inMax;i++)
      if (inHist[i]!=0)
        indicesIn[nonEmptyIn++]=i;

    for (i=0;i<=extMax;i++)
      if (extHist[i]!=0)
        ++nonEmptyExt;

    indicesExt =  vtmalloc( nonEmptyExt*sizeof(int), "indicesExt", proc );
    nonEmptyExt=0;
    for (i=0;i<=extMax;i++)
      if (extHist[i]!=0)
        indicesExt[nonEmptyExt++]=i;

    /*  tableau des intersections */
    table =  vtmalloc( nonEmptyExt*sizeof(int*), "table", proc );
    for(j=0;j<nonEmptyExt;j++)
    {
      table[j] =  vtmalloc( nonEmptyIn*sizeof(int), "table[j]", proc );
      for(i=0;i<nonEmptyIn;i++)
        table[j][i] = tab[indicesExt[j]][indicesIn[i]];
    }

    *Table = table;

    *indIn = indicesIn;
    *indExt = indicesExt;

    *inH = inHist;
    *extH = extHist;

    *sizeIn=nonEmptyIn;
    *sizeExt=nonEmptyExt;

    for (i=0;i<=extMax;i++) {
        vtfree(tab[i]);
        tab[i]=NULL;
    }
    vtfree(tab);
    tab=NULL;


    if (imageIn.type == UCHAR)
    {
      for (z=0;z<dimzIN;z++)
      {
        for (y=0;y<dimyIN;y++)
        {
          vtfree(in[z][y]);
          in[z][y]=NULL;
        }
        vtfree(in[z]);
        in[z]=NULL;
      }
      vtfree(in);
      in=NULL;
    }

    if (imageExt.type == UCHAR)
    {
      for (z=0;z<dimzEXT;z++)
      {
        for (y=0;y<dimyEXT;y++)
        {
          vtfree(ext[z][y]);
          ext[z][y]=NULL;
        }
        vtfree(ext[z]);
        ext[z]=NULL;
      }
      vtfree(ext);
      ext=NULL;
    }

    return(1);
}

int MT_ComputeIntersection3D( vt_image imageIn,
                             vt_image imageExt,
                             char *filename)
{
    int **tab=NULL;
    int *inHist=NULL, *extHist=NULL;
    int *indicesIn  = NULL;
    int *indicesExt = NULL;
    int nonEmptyIn, nonEmptyExt;
    int i,j;

    MT_ComputeTab3D(imageIn, imageExt, &tab, &inHist, &extHist, &indicesIn, &indicesExt, &nonEmptyIn, &nonEmptyExt);

    /* Generation du fichier */

    const char *extension = get_filename_ext(filename);

    /* fprintf(stdout, "Filename extension is %s\n", extension); */

    if ( strcmp ( extension, "txt" ) == 0 || strcmp ( extension, "" ) == 0 ) {  /*  Output format is a text file */
        FILE *fichier = fopen (filename, "w" );
        if (fichier == NULL)
        {
          perror (filename);
        }
        else
        {
          /* Intersection */
          fprintf(fichier, "%-10s", "Intersection");
          for (j=0;j<nonEmptyExt;j++)
          {
            fprintf(fichier, "%7d", indicesExt[j]);
          }
          fprintf(fichier, "\n\n\n");
          for (i=0;i<nonEmptyIn;i++)
          {
            fprintf(fichier, "%-10d", indicesIn[i]);
            for (j=0;j<nonEmptyExt;j++)
            {
              fprintf(fichier, "%10d", tab[j][i]);
            }
            fprintf(fichier, "\n\n");
          }

          /* Histogrammes */
          /*  HistIn */
          fprintf(fichier, "\n%-10s\n\n", "HistIn");
          for(i=0;i<nonEmptyIn;i++)
          {
            fprintf(fichier, "%10d\n", inHist[indicesIn[i]]);
          }
          /*  HistExt */
          fprintf(fichier, "\n%-10s\n\n", "HistExt");
          for(i=0;i<nonEmptyExt;i++)
          {
            fprintf(fichier, "%10d\n", extHist[indicesExt[i]]);
          }

        }
    }
    else {  /*  Output format is an image */
        vt_image imgTab;
        float ***array;
        VT_Image( &imgTab );
        VT_InitImage(&imgTab, filename, nonEmptyIn, nonEmptyExt, 1, FLOAT);
        VT_AllocImage(&imgTab);
        array = (float***) imgTab.array;
        for(j=0;j<nonEmptyExt;j++)
        for(i=0;i<nonEmptyIn;i++)
        {
            array[0][j][i]=(float)tab[j][i];
        }
        /* imgTab.buf = &tab; */
        VT_WriteInrimage(&imgTab);
        VT_FreeImage(&imgTab);
    }

    for(i=0;i<nonEmptyExt;i++) {
      vtfree(tab[i]);
      tab[i]=NULL;
    }
    vtfree(tab);
    tab=NULL;

    vtfree(inHist);
    inHist=NULL;
    vtfree(extHist);
    extHist=NULL;
    vtfree(indicesIn);
    indicesIn=NULL;
    vtfree(indicesExt);
    indicesExt=NULL;

    return(1);
}





int MT_ComputeDice3DTrsf( vt_image imageIn,
                             vt_image imageExt,
                          double *T_ext_in,
                             int **indicesInPtr,
                            int **indicesExtPtr,
                            int *nonEmptyInPtr,
                            int *nonEmptyExtPtr,
                            double ***DicePtr
                            )
{
  char *proc = "MT_ComputeDice3DTrsf";
  int dimxIN=0,dimyIN=0,dimzIN=0;
  int dimxEXT=0,dimyEXT=0,dimzEXT=0;
  int x,y,z;
  int inMax=0, extMax=0;
  int *inHist, *extHist;
  int **tab;
  double **Dice;

  int nonEmptyIn  = 0;
  int nonEmptyExt = 0;
  int *indicesIn  = NULL;
  int *indicesExt = NULL;

  unsigned char ***inu8;
  unsigned char ***extu8;

  unsigned short int ***in;
  unsigned short int ***ext;

  switch (imageIn.type)
  {
    case SCHAR:
    case UCHAR:
      inu8= (unsigned char ***)imageIn.array;
      dimxIN=imageIn.dim.x;
      dimyIN=imageIn.dim.y;
      dimzIN=imageIn.dim.z;
      in = vtmalloc(dimzIN*sizeof(unsigned short int **), "in", proc );
      for (z=0;z<dimzIN;z++)
      {
        in[z] = vtmalloc(dimyIN*sizeof(unsigned short int *), "in[z]", proc );
        for (y=0;y<dimyIN;y++)
        {
          in[z][y] = vtmalloc(dimxIN*sizeof(unsigned short int), "in[z][y]", proc );
          for(x=0;x<dimxIN;x++)
            in[z][y][x]=(unsigned short int)(inu8[z][y][x]);
        }
      }
      break;
    case SSHORT:
    case USHORT:
      in = (unsigned short int***)imageIn.array;
      break;
    case TYPE_UNKNOWN:
    default:
      fprintf(stderr,"Unknown type or not supported format of imageIn");
      return(0);
      break;
  }

  switch (imageExt.type)
  {
    case SCHAR:
    case UCHAR:
      extu8= (unsigned char ***)imageExt.array;
      dimxEXT=imageExt.dim.x;
      dimyEXT=imageExt.dim.y;
      dimzEXT=imageExt.dim.z;
      ext = vtmalloc(dimzEXT*sizeof(unsigned short int **), "ext", proc );
      for (z=0;z<dimzEXT;z++)
      {
        ext[z] = vtmalloc(dimyEXT*sizeof(unsigned short int *), "ext[z]", proc );
        for (y=0;y<dimyEXT;y++)
        {
          ext[z][y] = vtmalloc(dimxEXT*sizeof(unsigned short int), "ext[z][y]", proc );
          for(x=0;x<dimxEXT;x++)
            ext[z][y][x]=(unsigned short int)(extu8[z][y][x]);
        }
      }
      break;
    case SSHORT:
    case USHORT:
      ext = (unsigned short int***)imageExt.array;
      break;
    case TYPE_UNKNOWN:
    default:
      fprintf(stderr,"Unknown type or not supported format of imageIn");
      if (imageIn.type == UCHAR)
      {
        for (z=0;z<dimzIN;z++)
        {
          for (y=0;y<dimyIN;y++)
          {
            vtfree(in[z][y]);
            in[z][y]=NULL;
          }
          vtfree(in[z]);
          in[z]=NULL;
        }
        vtfree(in);
        in=NULL;
      }
      return(0);
      break;
  }

  int i,j;

  for (z=0;z<(int)imageIn.dim.z;z++)
    for (y=0;y<(int)imageIn.dim.y;y++)
      for (x=0;x<(int)imageIn.dim.x;x++)
      {
        if (inMax < (int)in[z][y][x])
          inMax = (int)in[z][y][x];
      }
  for (z=0;z<(int)imageExt.dim.z;z++)
    for (y=0;y<(int)imageExt.dim.y;y++)
      for (x=0;x<(int)imageExt.dim.x;x++)
      {
        if (extMax < (int)ext[z][y][x])
          extMax = (int)ext[z][y][x];
      }

  /*--allocation + initialisation tableaux qui serviront au dice--*/
  tab =  vtmalloc( (extMax+1)*sizeof(int*), "tab", proc );
  for (i=0;i<=extMax;i++)
    tab[i] =  vtmalloc( (inMax+1)*sizeof(int), "tab[i]", proc );
  inHist =  vtmalloc( (inMax+1)*sizeof(int), "inHist", proc );
  extHist = vtmalloc( (extMax+1)*sizeof(int), "extHist", proc );

  for (i=0;i<=inMax;i++) inHist[i]=0;
  for (i=0;i<=extMax;i++) extHist[i]=0;

  for (i=0;i<=inMax;i++)
    for (j=0;j<=extMax;j++)
      tab[j][i] = 0;

  int tx,ty,tz;

  for (z=0;z<(int)imageIn.dim.z;z++)
  for (y=0;y<(int)imageIn.dim.y;y++)
  for (x=0;x<(int)imageIn.dim.x;x++)
  {
        inHist[(int)in[z][y][x]]  += 1;
        tx=(int)(x*T_ext_in[0]+y*T_ext_in[1]+z*T_ext_in[2]+T_ext_in[3]+0.5);
        ty=(int)(x*T_ext_in[4]+y*T_ext_in[5]+z*T_ext_in[6]+T_ext_in[7]+0.5);
        tz=(int)(x*T_ext_in[8]+y*T_ext_in[9]+z*T_ext_in[10]+T_ext_in[11]+0.5);
        if (tx<0 || ty<0 || tz<0) continue;
        if (tx>=(int)imageExt.dim.x || ty>=(int)imageExt.dim.y || tz>=(int)imageExt.dim.z) continue;
        tab[(int)ext[tz][ty][tx]][(int)in[z][y][x]] += 1;
  }
  for (z=0;z<(int)imageExt.dim.z;z++)
  for (y=0;y<(int)imageExt.dim.y;y++)
  for (x=0;x<(int)imageExt.dim.x;x++)
  {
        extHist[(int)ext[z][y][x]]+= 1;
  }

  for (i=0;i<=inMax;i++)
  {
    if (inHist[i]!=0)
    {
      indicesIn = vtrealloc( indicesIn,(++nonEmptyIn)*sizeof(int), "indicesIn", proc );
      indicesIn[nonEmptyIn-1]=i;
    }
  }

  for (i=0;i<=extMax;i++)
  {
    if (extHist[i]!=0)
    {
      indicesExt = vtrealloc(indicesExt,(++nonEmptyExt)*sizeof(int), "indicesExt", proc );
      indicesExt[nonEmptyExt-1]=i;
    }
  }

  if (0)
  {
    for (j=0;j<nonEmptyExt;j++)
    {
      fprintf(stdout, "\t%d", indicesExt[j]);
    }
    fprintf(stdout, "\n");

    for (i=0;i<nonEmptyIn;i++)
    {
      fprintf(stdout, "%d", indicesIn[i]);
      for (j=0;j<nonEmptyExt;j++)
      {
        fprintf(stdout, "\t%d", tab[indicesExt[j]][indicesIn[i]]);
      }
      fprintf(stdout, "\n\n");
    }
  }

  Dice = vtmalloc(nonEmptyExt*sizeof(double*), "Dice", proc );
  for(j=0;j<nonEmptyExt;j++)
  {
    Dice[j] = vtmalloc(nonEmptyIn*sizeof(double), "Dice[j]", proc );
    for(i=0;i<nonEmptyIn;i++)
    {
      Dice[j][i] = 2*((double)tab[indicesExt[j]][indicesIn[i]])/
          ((double)extHist[indicesExt[j]]+(double)inHist[indicesIn[i]]);
    }
  }


  *DicePtr=Dice;
  *indicesInPtr=indicesIn;
  *indicesExtPtr=indicesExt;
  *nonEmptyInPtr=nonEmptyIn;
  *nonEmptyExtPtr=nonEmptyExt;

  /* Free memory */
  if (imageIn.type == UCHAR)
  {
    for (z=0;z<dimzIN;z++)
    {
      for (y=0;y<dimyIN;y++)
      {
        vtfree(in[z][y]);
        in[z][y]=NULL;
      }
      vtfree(in[z]);
      in[z]=NULL;
    }
    vtfree(in);
    in=NULL;
  }

  if (imageExt.type == UCHAR)
  {
    for (z=0;z<dimzEXT;z++)
    {
      for (y=0;y<dimyEXT;y++)
      {
        vtfree(ext[z][y]);
        ext[z][y]=NULL;
      }
      vtfree(ext[z]);
      ext[z]=NULL;
    }
    vtfree(ext);
    ext=NULL;
  }

  vtfree(inHist);
  vtfree(extHist);
  for(i=0;i<=extMax;i++) {
    vtfree(tab[i]);
    tab[i]=NULL;
  }
  vtfree(tab);
  tab=NULL;
  inHist=NULL;
  extHist=NULL;
  return (1);
}


int MT_ComputeDice2D( vt_image imageIn __attribute__ ((unused)),
                      vt_image imageExt __attribute__ ((unused)),
                      char *filename __attribute__ ((unused)) )
{       /* TODO */
  char *proc="MT_ComputeDice2D";
  fprintf(stderr,"%s: unimplemented function...\n",proc);




  return (1);
}











int MT_ComputeConfusionTab3D( vt_image labelIn,
                             vt_image imageExt,
                             char *filename )
{
  char *proc = "MT_ComputeConfusionTab3D";
  fprintf(stderr, "Warning: MT_ComputeConfusionTab3D method is obsolete, use rather MT_ComputeIntersection3D...");
  int dimxIN=0,dimyIN=0,dimzIN=0;
  int dimxEXT=0,dimyEXT=0,dimzEXT=0;
  int x,y,z;
  int inMax=0, extMax=0;
  int *inHist, *extHist;
  int **tab;
  int **confusionTab;

  int nonEmptyIn  = 0;
  int nonEmptyExt = 0;
  int *indicesIn  = NULL;
  int *indicesExt = NULL;

  unsigned char ***inu8;
  unsigned char ***extu8;
  unsigned short int ***in;
  unsigned short int ***ext;

  switch (labelIn.type)
  {
    case SCHAR:
    case UCHAR:
      inu8= (unsigned char ***)labelIn.array;
      dimxIN=labelIn.dim.x;
      dimyIN=labelIn.dim.y;
      dimzIN=labelIn.dim.z;
      in = vtmalloc(dimzIN*sizeof(unsigned short int **), "in", proc );
      for (z=0;z<dimzIN;z++)
      {
        in[z] = vtmalloc(dimyIN*sizeof(unsigned short int *), "in[z]", proc );
        for (y=0;y<dimyIN;y++)
        {
          in[z][y] = vtmalloc(dimxIN*sizeof(unsigned short int), "in[z][y]", proc );
          for(x=0;x<dimxIN;x++)
            in[z][y][x]=(unsigned short int)(inu8[z][y][x]);
        }
      }
      break;
    case SSHORT:
    case USHORT:
      in = (unsigned short int***)labelIn.array;
      break;
    case TYPE_UNKNOWN:
    default:
      fprintf(stderr,"Unknown type or not supported format of imageIn");
      return(0);
      break;
  }

  switch (imageExt.type)
  {
    case SCHAR:
    case UCHAR:
      extu8= (unsigned char ***)imageExt.array;
      dimxEXT=imageExt.dim.x;
      dimyEXT=imageExt.dim.y;
      dimzEXT=imageExt.dim.z;
      ext = vtmalloc(dimzEXT*sizeof(unsigned short int **), "ext", proc );
      for (z=0;z<dimzEXT;z++)
      {
        ext[z]=vtmalloc(dimyEXT*sizeof(unsigned short int *), "ext[z]", proc );
        for (y=0;y<dimyEXT;y++)
        {
          ext[z][y]=vtmalloc(dimxEXT*sizeof(unsigned short int), "ext[z][y]", proc );
          for(x=0;x<dimxEXT;x++)
            ext[z][y][x]=(unsigned short int)(extu8[z][y][x]);
        }
      }
      break;
    case SSHORT:
    case USHORT:
      ext = (unsigned short int***)imageExt.array;
      break;
    case TYPE_UNKNOWN:
    default:
      fprintf(stderr,"Unknown type or not supported format of imageIn");
      if (labelIn.type == UCHAR)
      {
        for (z=0;z<dimzIN;z++)
        {
          for (y=0;y<dimyIN;y++)
          {
            vtfree(in[z][y]);
            in[z][y]=NULL;
          }
          vtfree(in[z]);
          in[z]=NULL;
        }
        vtfree(in);
        in=NULL;
      }
      return(0);
      break;
  }

  int i,j;

  for (z=0;z<(int)labelIn.dim.z;z++)
    for (y=0;y<(int)labelIn.dim.y;y++)
      for (x=0;x<(int)labelIn.dim.x;x++)
      {
        if (inMax < (int)in[z][y][x])
          inMax = (int)in[z][y][x];
        if (extMax < (int)ext[z][y][x])
          extMax = (int)ext[z][y][x];
      }


  /*--allocation + initialisation tableaux qui serviront au dice--*/
  tab = vtmalloc((extMax+1)*sizeof(int*), "tab", proc );
  for (i=0;i<=extMax;i++)
    tab[i] = vtmalloc((inMax+1)*sizeof(int), "tab[i]", proc );
  inHist = vtmalloc((inMax+1)*sizeof(int), "inHist", proc );
  extHist =vtmalloc((extMax+1)*sizeof(int), "extHist", proc );

  for (i=0;i<=inMax;i++) inHist[i]=0;
  for (i=0;i<=extMax;i++) extHist[i]=0;

  for (i=0;i<=inMax;i++)
    for (j=0;j<=extMax;j++)
      tab[j][i] = 0;

  for (z=0;z<(int)labelIn.dim.z;z++)
    for (y=0;y<(int)labelIn.dim.y;y++)
      for (x=0;x<(int)labelIn.dim.x;x++)
      {
        inHist[(int)in[z][y][x]]  += 1;
        extHist[(int)ext[z][y][x]]+= 1;

        tab[(int)ext[z][y][x]][(int)in[z][y][x]] += 1;
      }

  if (_verbose_)
    fprintf(stdout, "inMax = %d\textMax = %d\n\n", inMax, extMax); /* ////////////////// */


  for (i=0;i<=inMax;i++)
  {
    if (inHist[i]!=0)
    {
      indicesIn = vtrealloc( indicesIn,(++nonEmptyIn)*sizeof(int), "indicesIn", proc );
      indicesIn[nonEmptyIn-1]=i;
    }
  }

  for (i=0;i<=extMax;i++)
  {
    if (extHist[i]!=0)
    {
      indicesExt = vtrealloc( indicesExt,(++nonEmptyExt)*sizeof(int), "indicesExt", proc );
      indicesExt[nonEmptyExt-1]=i;
    }
  }

  if (_verbose_)
  {
    for (j=0;j<nonEmptyExt;j++)
    {
      fprintf(stdout, "\t%d", indicesExt[j]);
    }
    fprintf(stdout, "\n");

    for (i=0;i<nonEmptyIn;i++)
    {
      fprintf(stdout, "%d", indicesIn[i]);
      for (j=0;j<nonEmptyExt;j++)
      {
        fprintf(stdout, "\t%d", tab[indicesExt[j]][indicesIn[i]]);
      }
      fprintf(stdout, "\n\n");
    }
  }

  confusionTab = vtmalloc( nonEmptyExt*sizeof(int*), "confusionTab", proc );
  for(j=0;j<nonEmptyExt;j++)
  {
    confusionTab[j] = vtmalloc(nonEmptyIn*sizeof(int), "confusionTab[j]", proc );
    for(i=0;i<nonEmptyIn;i++)
    {
      confusionTab[j][i] = tab[indicesExt[j]][indicesIn[i]]>0 ? 1 : 0;
    }
  }

  if (_verbose_)
  {
    for (j=0;j<nonEmptyExt;j++)
    {
      fprintf(stdout, "\t%d", indicesExt[j]);
    }
    fprintf(stdout, "\n");

    for (i=0;i<nonEmptyIn;i++)
    {
      fprintf(stdout, "%d", indicesIn[i]);
      for (j=0;j<nonEmptyExt;j++)
      {
        fprintf(stdout, "\t%d", confusionTab[j][i]);
      }
      fprintf(stdout, "\n\n");
    }
  }

  /* Generation du fichier */

  FILE *fichier = fopen (filename, "w" );
  if (fichier == NULL)
  {
    perror (filename);
  }
  else
  {
    /* tab */
    fprintf(fichier, "\n%-10s", "tab");
    for (j=0;j<nonEmptyExt;j++)
    {
      fprintf(fichier, "%10d", indicesExt[j]);
    }
    fprintf(fichier, "\n\n\n");
    for (i=0;i<nonEmptyIn;i++)
    {
      fprintf(fichier, "%-10d", indicesIn[i]);
      for (j=0;j<nonEmptyExt;j++)
      {
        fprintf(fichier, "%10d", tab[indicesExt[j]][indicesIn[i]]);
      }
      fprintf(fichier, "\n\n");
    }


    /* confusion */
    fprintf(fichier, "%-10s", "confusion");
    for (j=0;j<nonEmptyExt;j++)
    {
      fprintf(fichier, "%7d", indicesExt[j]);
    }
    fprintf(fichier, "\n\n\n");
    for (i=0;i<nonEmptyIn;i++)
    {
      fprintf(fichier, "%-10d", indicesIn[i]);
      for (j=0;j<nonEmptyExt;j++)
      {
        fprintf(fichier, "%7d", confusionTab[j][i]);
      }
      fprintf(fichier, "\n\n");
    }

    /* Histogrammes */
    /*  HistIn */
    fprintf(fichier, "\n%-10s\n\n", "HistIn");
    for(i=0;i<nonEmptyIn;i++)
    {
      fprintf(fichier, "%10d\n", inHist[indicesIn[i]]);
    }
    /*  HistExt */
    fprintf(fichier, "\n%-10s\n\n", "HistExt");
    for(i=0;i<nonEmptyExt;i++)
    {
      fprintf(fichier, "%10d\n", extHist[indicesExt[i]]);
    }

  }

  if (labelIn.type == UCHAR)
  {
    for (z=0;z<dimzIN;z++)
    {
      for (y=0;y<dimyIN;y++)
      {
        vtfree(in[z][y]);
        in[z][y]=NULL;
      }
      vtfree(in[z]);
      in[z]=NULL;
    }
    vtfree(in);
    in=NULL;
  }

  if (imageExt.type == UCHAR)
  {
    for (z=0;z<dimzEXT;z++)
    {
      for (y=0;y<dimyEXT;y++)
      {
        vtfree(ext[z][y]);
        ext[z][y]=NULL;
      }
      vtfree(ext[z]);
      ext[z]=NULL;
    }
    vtfree(ext);
    ext=NULL;
  }


  vtfree(inHist);
  vtfree(extHist);
  for(i=0;i<=extMax;i++) {
    vtfree(tab[i]);
    tab[i]=NULL;
  }
  for(i=0;i<nonEmptyExt;i++)
  {
    vtfree(confusionTab[i]);
    confusionTab[i]=NULL;
  }
  vtfree(tab);
  vtfree(confusionTab);
  tab=NULL;
  inHist=NULL;
  extHist=NULL;
  return (1);
}
















/* ////////////////////////////////////////////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////////////////////// */
/* //////////////////////////// SYMMETRY //////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////////////////////// */


int MT_ComputeSymmetryDice3D( vt_image* imageIn,
                             int x_pos,
                             char *filename )
{
  char *proc = "MT_ComputeSymmetryDice3D";
  int dimx,dimy,dimz;
  int imin, imax;
  int x,y,z;
  int gMax=0, dMax=0;
  int *gHist, *dHist;
  int **tab;
  double **Dice;

  int i,j;

  int nonEmptyG  = 0;
  int nonEmptyD = 0;
  int *indicesG  = NULL;
  int *indicesD = NULL;

  unsigned char ***inu8;

  unsigned short int ***in;

  dimx=imageIn->dim.x;
  dimy=imageIn->dim.y;
  dimz=imageIn->dim.z;

  switch (imageIn->type)
  {
    case SCHAR:
    case UCHAR:
      inu8= (unsigned char ***)imageIn->array;
      dimx=imageIn->dim.x;
      dimy=imageIn->dim.y;
      dimz=imageIn->dim.z;
      in = vtmalloc( dimz*sizeof(unsigned short int **), "in", proc );
      for (z=0;z<dimz;z++)
      {
        in[z] = vtmalloc( dimy*sizeof(unsigned short int *), "in[z]", proc );
        for (y=0;y<dimy;y++)
        {
          in[z][y] = vtmalloc( dimx*sizeof(unsigned short int), "in[z][y]", proc );
          for(x=0;x<dimx;x++)
            in[z][y][x]=(unsigned short int)(inu8[z][y][x]);
        }
      }
      break;
    case SSHORT:
    case USHORT:
      in = (unsigned short int***)imageIn->array;
      break;
    case TYPE_UNKNOWN:
    default:
      fprintf(stderr,"Unknown type or not supported format of image\n");
      return(0);
      break;
  }

  imax=(x_pos>dimx-x_pos-1) ? x_pos : dimx-x_pos-1;
  imin=(x_pos<dimx-x_pos-1) ? x_pos : dimx-x_pos-1;
  for (z=0;z<dimz;z++)
    for (y=0;y<dimy;y++)
      for (x=1;x<=imax;x++)
      {
        if (x_pos-x>=0 && gMax < (int)in[z][y][x_pos-x])
          gMax = (int)in[z][y][x_pos-x];
        if (x_pos+x<dimx && dMax < (int)in[z][y][x_pos+x])
          dMax = (int)in[z][y][x_pos+x];
      }

  /*--allocation + initialisation tableaux qui serviront au dice--*/
  tab = vtmalloc((dMax+1)*sizeof(int*), "tab", proc );
  for (i=0;i<=dMax;i++)
    tab[i] = vtmalloc((gMax+1)*sizeof(int), "tab[i]", proc );
  gHist = vtmalloc((gMax+1)*sizeof(int), "gHist", proc );
  dHist = vtmalloc((dMax+1)*sizeof(int), "dHist", proc );

  for (i=0;i<=gMax;i++) gHist[i]=0;
  for (i=0;i<=dMax;i++) dHist[i]=0;

  for (i=0;i<=gMax;i++)
  for (j=0;j<=dMax;j++)
    tab[j][i] = 0;

  for (z=0;z<dimz;z++)
    for (y=0;y<dimy;y++)
      for (x=1;x<=imax;x++)
      {
        if (x_pos-x >= 0) gHist[(int)in[z][y][x_pos-x]] += 1;
        if (x_pos+x<dimx) dHist[(int)in[z][y][x_pos+x]] += 1;

        if (x<=imin) tab[(int)in[z][y][x_pos+x]][(int)in[z][y][x_pos-x]] += 1;
      }

  for (i=0;i<=gMax;i++)
  {
    if (gHist[i]!=0)
    {
      indicesG = vtrealloc(indicesG,(++nonEmptyG)*sizeof(int), "indicesG", proc );
      indicesG[nonEmptyG-1]=i;
    }
  }

  for (i=0;i<=dMax;i++)
  {
    if (dHist[i]!=0)
    {
      indicesD = vtrealloc(indicesD,(++nonEmptyD)*sizeof(int), "indicesD", proc );
      indicesD[nonEmptyD-1]=i;
    }
  }

  if (_verbose_)
  {
    for (j=0;j<nonEmptyD;j++)
    {
      fprintf(stdout, "\t%d", indicesD[j]);
    }
    fprintf(stdout, "\n");

    for (i=0;i<nonEmptyG;i++)
    {
      fprintf(stdout, "%d", indicesG[i]);
      for (j=0;j<nonEmptyD;j++)
      {
        fprintf(stdout, "\t%d", tab[indicesD[j]][indicesG[i]]);
      }
      fprintf(stdout, "\n\n");
    }
  }

  Dice = vtmalloc( nonEmptyD*sizeof(double*), "Dice", proc );
  for(j=0;j<nonEmptyD;j++)
  {
    Dice[j] = vtmalloc( nonEmptyG*sizeof(double), "Dice[j]", proc );
    for(i=0;i<nonEmptyG;i++)
    {
      Dice[j][i] = 2*((double)tab[indicesD[j]][indicesG[i]])/
          ((double)dHist[indicesD[j]]+(double)gHist[indicesG[i]]);
    }
  }

  if (_verbose_)
  {
    for (j=0;j<=nonEmptyD;j++)
    {
      fprintf(stdout, "\t%d", indicesD[j]);
    }
    fprintf(stdout, "\n");

    for (i=0;i<nonEmptyG;i++)
    {
      fprintf(stdout, "%d", indicesG[i]);
      for (j=0;j<nonEmptyD;j++)
      {
        fprintf(stdout, "\t%.3f", Dice[j][i]);
      }
      fprintf(stdout, "\n\n");
    }
  }

  /* Generation du fichier */

  FILE *fichier = fopen (filename, "w" );
  if (fichier == NULL)
  {
    perror (filename);
  }
  else
  {
    /* Dice */
    fprintf(fichier, "%-10s", "Dice");
    for (j=0;j<nonEmptyD;j++)
    {
      fprintf(fichier, "%7d", indicesD[j]);
    }
    fprintf(fichier, "\n\n\n");
    for (i=0;i<nonEmptyG;i++)
    {
      fprintf(fichier, "%-10d", indicesG[i]);
      for (j=0;j<nonEmptyD;j++)
      {
        fprintf(fichier, "%7.3f", Dice[j][i]);
      }
      fprintf(fichier, "\n\n");
    }

    /* Histogrammes */
    /*  HistIn */
    fprintf(fichier, "\n%-10s\n\n", "HistGauche");
    for(i=0;i<nonEmptyG;i++)
    {
      fprintf(fichier, "%10d\n", gHist[indicesG[i]]);
    }
    /*  HistExt */
    fprintf(fichier, "\n%-10s\n\n", "HistDroite");
    for(i=0;i<nonEmptyD;i++)
    {
      fprintf(fichier, "%10d\n", dHist[indicesD[i]]);
    }

  }

  if (imageIn->type == UCHAR)
  {
    for (z=0;z<dimz;z++)
    {
      for (y=0;y<dimy;y++)
      {
        vtfree(in[z][y]);
        in[z][y]=NULL;
      }
      vtfree(in[z]);
      in[z]=NULL;
    }
    vtfree(in);
    in=NULL;
  }

  vtfree(gHist);
  vtfree(dHist);
  for(i=0;i<=dMax;i++) {
    vtfree(tab[i]);
    tab[i]=NULL;
  }
  vtfree(tab);
  tab=NULL;
  gHist=NULL;
  dHist=NULL;
  return (1);
}





















int MT_ComputeSymmetryDice3DBis( vt_image* imageIn,
                             int x_pos,
                             double *r )
{
  char *proc = "MT_ComputeSymmetryDice3DBis";
  int dimx,dimy,dimz;
  int imin, imax;
  int x,y,z;
  int gMax=0, dMax=0;
  int *gHist, *dHist;
  int **tab;
  double **Dice;

  int i,j;

  int nonEmptyG  = 0;
  int nonEmptyD = 0;
  int *indicesG  = NULL;
  int *indicesD = NULL;

  unsigned char ***inu8;

  unsigned short int ***in;

  dimx=imageIn->dim.x;
  dimy=imageIn->dim.y;
  dimz=imageIn->dim.z;

  switch (imageIn->type)
  {
    case SCHAR:
    case UCHAR:
      inu8= (unsigned char ***)imageIn->array;
      dimx=imageIn->dim.x;
      dimy=imageIn->dim.y;
      dimz=imageIn->dim.z;
      in = vtmalloc(dimz*sizeof(unsigned short int **), "in", proc );
      for (z=0;z<dimz;z++)
      {
        in[z] = vtmalloc(dimy*sizeof(unsigned short int *), "in[z]", proc );
        for (y=0;y<dimy;y++)
        {
          in[z][y] = vtmalloc(dimx*sizeof(unsigned short int), "in[z][y]", proc );
          for(x=0;x<dimx;x++)
            in[z][y][x]=(unsigned short int)(inu8[z][y][x]);
        }
      }
      break;
    case SSHORT:
    case USHORT:
      in = (unsigned short int***)imageIn->array;
      break;
    case TYPE_UNKNOWN:
    default:
      fprintf(stderr,"Unknown type or not supported format of imageIn");
      return(0);
      break;
  }

  imax=(x_pos>dimx-x_pos-1) ? x_pos : dimx-x_pos-1;
  imin=(x_pos<dimx-x_pos-1) ? x_pos : dimx-x_pos-1;
  for (z=0;z<dimz;z++)
    for (y=0;y<dimy;y++)
      for (x=1;x<=imax;x++)
      {
        if (x_pos-x>=0 && gMax < (int)in[z][y][x_pos-x])
          gMax = (int)in[z][y][x_pos-x];
        if (x_pos+x<dimx && dMax < (int)in[z][y][x_pos+x])
          dMax = (int)in[z][y][x_pos+x];
      }

  /*--allocation + initialisation tableaux qui serviront au dice--*/
  tab = vtmalloc((dMax+1)*sizeof(int*), "tab", proc );
  for (i=0;i<=dMax;i++)
    tab[i] = vtmalloc((gMax+1)*sizeof(int), "tab[i]", proc );
  gHist = vtmalloc((gMax+1)*sizeof(int), "gHist", proc );
  dHist = vtmalloc((dMax+1)*sizeof(int), "dHist", proc );

  for (i=0;i<=gMax;i++) gHist[i]=0;
  for (i=0;i<=dMax;i++) dHist[i]=0;

  for (i=0;i<=gMax;i++)
  for (j=0;j<=dMax;j++)
    tab[j][i] = 0;

  for (z=0;z<dimz;z++)
    for (y=0;y<dimy;y++)
      for (x=1;x<=imax;x++)
      {
        if (x_pos-x >= 0) gHist[(int)in[z][y][x_pos-x]] += 1;
        if (x_pos+x<dimx) dHist[(int)in[z][y][x_pos+x]] += 1;

        if (x<=imin) tab[(int)in[z][y][x_pos+x]][(int)in[z][y][x_pos-x]] += 1;
      }

  for (i=0;i<=gMax;i++)
  {
    if (gHist[i]!=0)
    {
      indicesG = vtrealloc(indicesG,(++nonEmptyG)*sizeof(int), "indicesG", proc );
      indicesG[nonEmptyG-1]=i;
    }
  }

  for (i=0;i<=dMax;i++)
  {
    if (dHist[i]!=0)
    {
      indicesD = vtrealloc(indicesD,(++nonEmptyD)*sizeof(int), "indicesD", proc );
      indicesD[nonEmptyD-1]=i;
    }
  }

  if (_verbose_)
  {
    for (j=0;j<nonEmptyD;j++)
    {
      fprintf(stdout, "\t%d", indicesD[j]);
    }
    fprintf(stdout, "\n");

    for (i=0;i<nonEmptyG;i++)
    {
      fprintf(stdout, "%d", indicesG[i]);
      for (j=0;j<nonEmptyD;j++)
      {
        fprintf(stdout, "\t%d", tab[indicesD[j]][indicesG[i]]);
      }
      fprintf(stdout, "\n\n");
    }
  }

  Dice = vtmalloc(nonEmptyD*sizeof(double*), "Dice", proc );
  for(j=0;j<nonEmptyD;j++)
  {
    Dice[j] = vtmalloc(nonEmptyG*sizeof(double), "Dice[j]", proc );
    for(i=0;i<nonEmptyG;i++)
    {
      Dice[j][i] = 2*((double)tab[indicesD[j]][indicesG[i]])/
          ((double)dHist[indicesD[j]]+(double)gHist[indicesG[i]]);
    }
  }

  *r=Dice[nonEmptyD-1][nonEmptyG-1];

  if (_verbose_)
  {
    for (j=0;j<=nonEmptyD;j++)
    {
      fprintf(stdout, "\t%d", indicesD[j]);
    }
    fprintf(stdout, "\n");

    for (i=0;i<nonEmptyG;i++)
    {
      fprintf(stdout, "%d", indicesG[i]);
      for (j=0;j<nonEmptyD;j++)
      {
        fprintf(stdout, "\t%.3f", Dice[j][i]);
      }
      fprintf(stdout, "\n\n");
    }
  }


  if (imageIn->type == UCHAR)
  {
    for (z=0;z<dimz;z++)
    {
      for (y=0;y<dimy;y++)
      {
        vtfree(in[z][y]);
        in[z][y]=NULL;
      }
      vtfree(in[z]);
      in[z]=NULL;
    }
    vtfree(in);
    in=NULL;
  }

  vtfree(gHist);
  vtfree(dHist);
  for(i=0;i<=dMax;i++) {
    vtfree(tab[i]);
    tab[i]=NULL;
  }
  vtfree(tab);
  tab=NULL;
  gHist=NULL;
  dHist=NULL;
  return (1);
}





















int MT_ComputeSymmetryDiceOblic3D( vt_image* imageIn,
                             double *n,
                             char *filename )
{
  char *proc = "MT_ComputeSymmetryDiceOblic3D";
  int dimx,dimy,dimz;
  int x,y,z;
  int symx, symy, symz;
  int gMax=0, dMax=0;
  int *gHist, *dHist;
  int **tab;
  double **Dice;
  double lambda;
  
  
  int i,j, tmp;

  int nonEmptyG  = 0;
  int nonEmptyD = 0;
  int *indicesG  = NULL;
  int *indicesD = NULL;

  unsigned char ***inu8;

  unsigned short int ***in;

  dimx=imageIn->dim.x;
  dimy=imageIn->dim.y;
  dimz=imageIn->dim.z;

  switch (imageIn->type)
  {
    case SCHAR:
    case UCHAR:
      inu8= (unsigned char ***)imageIn->array;
      dimx=imageIn->dim.x;
      dimy=imageIn->dim.y;
      dimz=imageIn->dim.z;
      in = vtmalloc(dimz*sizeof(unsigned short int **), "in", proc );
      for (z=0;z<dimz;z++)
      {
        in[z]=vtmalloc(dimy*sizeof(unsigned short int *), "in[z]", proc );
        for (y=0;y<dimy;y++)
        {
          in[z][y]=vtmalloc(dimx*sizeof(unsigned short int), "in[z][y]", proc );
          for(x=0;x<dimx;x++)
            in[z][y][x]=(unsigned short int)(inu8[z][y][x]);
        }
      }
      break;
    case SSHORT:
    case USHORT:
      in = (unsigned short int***)imageIn->array;
      break;
    case TYPE_UNKNOWN:
    default:
      fprintf(stderr,"Unknown type or not supported format of imageIn");
      return(0);
      break;
  }

  for (z=0;z<dimz;z++)
  for (y=0;y<dimy;y++)
  for (x=0;x<dimx;x++)
  {
    tmp=n[0]*x*imageIn->siz.x+n[1]*y*imageIn->siz.y+n[2]*z*imageIn->siz.z+n[3];
    if (tmp<0)
	{
	  /*  Gauche */
	  if(gMax < (int)in[z][y][x])
		gMax = (int)in[z][y][x];
	}
    if (tmp>0)
	{
	  /*  Droite */
	  if(dMax < (int)in[z][y][x])
		dMax = (int)in[z][y][x];
	}
  }

  /*--allocation + initialisation tableaux qui serviront au dice--*/
  tab = vtmalloc((dMax+1)*sizeof(int*), "tab", proc );
  for (i=0;i<=dMax;i++)
    tab[i] = vtmalloc((gMax+1)*sizeof(int), "tab[i]", proc );
  gHist = vtmalloc((gMax+1)*sizeof(int), "gHist", proc );
  dHist = vtmalloc((dMax+1)*sizeof(int), "dHist", proc );

  for (i=0;i<=gMax;i++) gHist[i]=0;
  for (i=0;i<=dMax;i++) dHist[i]=0;

  for (j=0;j<=dMax;j++)
  for (i=0;i<=gMax;i++)
    tab[j][i] = 0;


  for (z=0;z<dimz;z++)
  for (y=0;y<dimy;y++)
  for (x=0;x<dimx;x++)
  {
    lambda=n[0]*x*imageIn->siz.x+n[1]*y*imageIn->siz.y+n[2]*z*imageIn->siz.z+n[3];
    if (lambda>1) {
	  dHist[(int)in[z][y][x]] += 1;
	}
    if (lambda<-1) {
	  gHist[(int)in[z][y][x]] += 1;
      symx=(int)(x-2*lambda*n[0]/imageIn->siz.x+0.5);
      symy=(int)(y-2*lambda*n[1]/imageIn->siz.y+0.5);
      symz=(int)(z-2*lambda*n[2]/imageIn->siz.x+0.5);
	  if ( symx<0 || symx>=dimx ||symy<0 || symy>=dimy ||symz<0 || symz>=dimz )
		continue;
	  /* lambda=n[0]*symx+n[1]*symy+n[2]*symz+n[3]; */
	  /* if ((int)in[z][y][x] > 0 ) */
	  /* fprintf(stdout, "{%d %d %d}: %d \t[%d %d %d]: %d\n" , x,y,z, (int)in[z][y][x], symx, symy, symz,(int)in[symz][symy][symx]); */
	  tab[(int)in[symz][symy][symx]][(int)in[z][y][x]] += 1;
	}
  }

  for (i=0;i<=gMax;i++)
  {
    if (gHist[i]!=0)
    {
      indicesG = vtrealloc(indicesG,(++nonEmptyG)*sizeof(int), "indicesG", proc );
      indicesG[nonEmptyG-1]=i;
    }
  }

  for (i=0;i<=dMax;i++)
  {
    if (dHist[i]!=0)
    {
      indicesD = vtrealloc(indicesD,(++nonEmptyD)*sizeof(int), "indicesD", proc );
      indicesD[nonEmptyD-1]=i;
    }
  }

  /* fprintf(stdout, "gMax : %d \t dMax = %d \t nonEmptyG = %d \t nonEmptyD = %d\n", gMax, dMax, nonEmptyG, nonEmptyD); */
  
  if (_verbose_ )
  {
    for (j=0;j<nonEmptyD;j++)
    {
      fprintf(stdout, "\t%d", indicesD[j]);
    }
    fprintf(stdout, "\n");

    for (i=0;i<nonEmptyG;i++)
    {
      fprintf(stdout, "%d", indicesG[i]);
      for (j=0;j<nonEmptyD;j++)
      {
        fprintf(stdout, "\t%d", tab[indicesD[j]][indicesG[i]]);
      }
      fprintf(stdout, "\n\n");
    }
  }
  

  Dice = vtmalloc(nonEmptyD*sizeof(double*), "Dice", proc );
  for(j=0;j<nonEmptyD;j++)
  {
    Dice[j] = vtmalloc(nonEmptyG*sizeof(double), "Dice[j]", proc );
    for(i=0;i<nonEmptyG;i++)
    {
      Dice[j][i] = 2*((double)tab[indicesD[j]][indicesG[i]])/
          ((double)dHist[indicesD[j]]+(double)gHist[indicesG[i]]);
    }
  }

  if (_verbose_ )
  {
    for (j=0;j<nonEmptyD;j++)
    {
      fprintf(stdout, "\t%d", indicesD[j]);
    }
    fprintf(stdout, "\n");

    for (i=0;i<nonEmptyG;i++)
    {
      fprintf(stdout, "%d", indicesG[i]);
      for (j=0;j<nonEmptyD;j++)
      {
        fprintf(stdout, "\t%.3f", Dice[j][i]);
      }
      fprintf(stdout, "\n\n");
    }
  }

  /* Generation du fichier */

  FILE *fichier = fopen (filename, "w" );
  if (fichier == NULL)
  {
    perror (filename);
  }
  else
  {
    /* Dice */
    fprintf(fichier, "%-10s", "Dice");
    for (j=0;j<nonEmptyD;j++)
    {
      fprintf(fichier, "%7d", indicesD[j]);
    }
    fprintf(fichier, "\n\n\n");
    for (i=0;i<nonEmptyG;i++)
    {
      fprintf(fichier, "%-10d", indicesG[i]);
      for (j=0;j<nonEmptyD;j++)
      {
        fprintf(fichier, "%7.3f", Dice[j][i]);
      }
      fprintf(fichier, "\n\n");
    }

    /* Histogrammes */
    /*  HistIn */
    fprintf(fichier, "\n%-10s\n\n", "HistGauche");
    for(i=0;i<nonEmptyG;i++)
    {
      fprintf(fichier, "%10d\n", gHist[indicesG[i]]);
    }
    /*  HistExt */
    fprintf(fichier, "\n%-10s\n\n", "HistDroite");
    for(i=0;i<nonEmptyD;i++)
    {
      fprintf(fichier, "%10d\n", dHist[indicesD[i]]);
    }

  }

  if (imageIn->type == UCHAR)
  {
    for (z=0;z<dimz;z++)
    {
      for (y=0;y<dimy;y++)
      {
        vtfree(in[z][y]);
        in[z][y]=NULL;
      }
      vtfree(in[z]);
      in[z]=NULL;
    }
    vtfree(in);
    in=NULL;
  }

  vtfree(gHist);
  vtfree(dHist);
  for(i=0;i<=dMax;i++) {
    vtfree(tab[i]);
    tab[i]=NULL;
  }
  vtfree(tab);
  tab=NULL;
  gHist=NULL;
  dHist=NULL;
  return (1);
}











int MT_ComputeSymmetryDiceOblic3DBis( vt_image* imageIn,
                             double *n,
                             double *r )
{
  char *proc="MT_ComputeSymmetryDiceOblic3DBis";
  int dimx,dimy,dimz;
  int x,y,z;
  int symx, symy, symz;
  int gMax=0, dMax=0;
  int *gHist, *dHist;
  int **tab;
  double **Dice;
  double lambda;
  
  
  int i,j;

  int nonEmptyG  = 0;
  int nonEmptyD = 0;
  int *indicesG  = NULL;
  int *indicesD = NULL;

  unsigned char ***inu8;

  unsigned short int ***in;

  dimx=imageIn->dim.x;
  dimy=imageIn->dim.y;
  dimz=imageIn->dim.z;

  switch (imageIn->type)
  {
    case SCHAR:
    case UCHAR:
      inu8= (unsigned char ***)imageIn->array;
      dimx=imageIn->dim.x;
      dimy=imageIn->dim.y;
      dimz=imageIn->dim.z;
      in = vtmalloc(dimz*sizeof(unsigned short int **), "in", proc );
      for (z=0;z<dimz;z++)
      {
        in[z] = vtmalloc(dimy*sizeof(unsigned short int *), "in[z]", proc );
        for (y=0;y<dimy;y++)
        {
          in[z][y] = vtmalloc(dimx*sizeof(unsigned short int), "in[z][y]", proc );
          for(x=0;x<dimx;x++)
            in[z][y][x]=(unsigned short int)(inu8[z][y][x]);
        }
      }
      break;
    case SSHORT:
    case USHORT:
      in = (unsigned short int***)imageIn->array;
      break;
    case TYPE_UNKNOWN:
    default:
      fprintf(stderr,"%s: Unknown type or not supported format of imageIn", proc);
      return(0);
      break;
  }

  /* imax=(x_pos>dimx-x_pos-1) ? x_pos : dimx-x_pos-1; */
  /* imin=(x_pos<dimx-x_pos-1) ? x_pos : dimx-x_pos-1; */
  for (z=0;z<dimz;z++)
  for (y=0;y<dimy;y++)
  for (x=0;x<dimx;x++)
  {
	if (n[0]*x+n[1]*y+n[2]*z+n[3]<-1)
	{
	  /*  Gauche */
	  if(gMax < (int)in[z][y][x])
		gMax = (int)in[z][y][x];
	}
	if (n[0]*x+n[1]*y+n[2]*z+n[3]>1)
	{
	  /*  Droite */
	  if(dMax < (int)in[z][y][x])
		dMax = (int)in[z][y][x];
	}
  }

  /*--allocation + initialisation tableaux qui serviront au dice--*/
  tab = vtmalloc((dMax+1)*sizeof(int*), "tab", proc );
  if(tab == NULL) {
	if (imageIn->type == UCHAR) {
	  for (z=0;z<dimz;z++)
      {
		for (y=0;y<dimy;y++)
		{
          vtfree(in[z][y]);
          in[z][y]=NULL;
		}
		free(in[z]);
		in[z]=NULL;
      }
	  vtfree(in);
      in=NULL;
	}
	fprintf(stderr, "%s: unable to allocate tab\n", proc);
	return(0);
  }
  for (i=0;i<=dMax;i++) {
    tab[i] = vtmalloc((gMax+1)*sizeof(int), "tab[i]", proc );
    if (tab[i]==NULL) {
	  if (imageIn->type == UCHAR) {
		for (z=0;z<dimz;z++)
		{
		  for (y=0;y<dimy;y++)
		  {
			free(in[z][y]);
			in[z][y]=NULL;
		  }
		  vtfree(in[z]);
		  in[z]=NULL;
		}
		free(in);
		in=NULL;
	  }
	  for (j=0; j<i ; j++) vtfree(tab[j]);
	  vtfree(tab);
	  fprintf(stderr, "%s: unable to allocate tab[%d]\n", proc, i);
	  return(0);
	}
  }
  gHist = vtmalloc((gMax+1)*sizeof(int), "gHist", proc );
  if (gHist==NULL) {
	if (imageIn->type == UCHAR) {
	  for (z=0;z<dimz;z++)
	  {
		for (y=0;y<dimy;y++)
		{
		  vtfree(in[z][y]);
		  in[z][y]=NULL;
		}
		free(in[z]);
		in[z]=NULL;
	  }
	  vtfree(in);
	  in=NULL;
	}
	for (i=0; i<=dMax ; i++) vtfree(tab[i]);
	free(tab);
	fprintf(stderr, "%s: unable to allocate gHist\n", proc);
	return(0);
  }
  
  dHist = vtmalloc((dMax+1)*sizeof(int), "dHist", proc );
  if (dHist==NULL) {
	if (imageIn->type == UCHAR) {
	  for (z=0;z<dimz;z++)
	  {
		for (y=0;y<dimy;y++)
		{
		  vtfree(in[z][y]);
		  in[z][y]=NULL;
		}
		free(in[z]);
		in[z]=NULL;
	  }
	  vtfree(in);
	  in=NULL;
	}
	for (i=0; i<=dMax ; i++) vtfree(tab[i]);
	free(tab);
	free(gHist);
	fprintf(stderr, "%s: unable to allocate gHist\n", proc);
	return(0);
  }

  for (i=0;i<=gMax;i++) gHist[i]=0;
  for (i=0;i<=dMax;i++) dHist[i]=0;

  for (j=0;j<=dMax;j++)
  for (i=0;i<=gMax;i++)
    tab[j][i] = 0;

  for (z=0;z<dimz;z++)
  for (y=0;y<dimy;y++)
  for (x=0;x<dimx;x++)
  {
	lambda=n[0]*x+n[1]*y+n[2]*z+n[3];
    if (lambda>1) {
	  dHist[(int)in[z][y][x]] += 1;
	}
    if (lambda<-1) {
	  gHist[(int)in[z][y][x]] += 1;
	  symx=(int)(x-2*lambda*n[0]+0.5);
	  symy=(int)(y-2*lambda*n[1]+0.5);
	  symz=(int)(z-2*lambda*n[2]+0.5);
	  if ( symx<0 || symx>=dimx ||symy<0 || symy>=dimy ||symz<0 || symz>=dimz )
		continue;
	  tab[(int)in[symz][symy][symx]][(int)in[z][y][x]] += 1;
	}
  }

  for (i=0;i<=gMax;i++)
  {
    if (gHist[i]!=0)
    {
      indicesG = vtrealloc(indicesG,(++nonEmptyG)*sizeof(int), "indicesG", proc );
      indicesG[nonEmptyG-1]=i;
    }
  }

  for (i=0;i<=dMax;i++)
  {
    if (dHist[i]!=0)
    {
      indicesD = vtrealloc(indicesD,(++nonEmptyD)*sizeof(int), "indicesD", proc );
      indicesD[nonEmptyD-1]=i;
    }
  }

  if (_verbose_)
  {
    for (j=0;j<nonEmptyD;j++)
    {
      fprintf(stdout, "\t%d", indicesD[j]);
    }
    fprintf(stdout, "\n");

    for (i=0;i<nonEmptyG;i++)
    {
      fprintf(stdout, "%d", indicesG[i]);
      for (j=0;j<nonEmptyD;j++)
      {
        fprintf(stdout, "\t%d", tab[indicesD[j]][indicesG[i]]);
      }
      fprintf(stdout, "\n\n");
    }
  }

  Dice = vtmalloc( nonEmptyD*sizeof(double*), "Dice", proc );
  for(j=0;j<nonEmptyD;j++)
  {
    Dice[j] = vtmalloc( nonEmptyG*sizeof(double), "Dice[j]", proc );
    for(i=0;i<nonEmptyG;i++)
    {
      Dice[j][i] = 2*((double)tab[indicesD[j]][indicesG[i]])/
          ((double)dHist[indicesD[j]]+(double)gHist[indicesG[i]]);
    }
  }
  
  *r =  Dice[nonEmptyD-1][nonEmptyG-1];
  
  if (_verbose_)
  {
    for (j=0;j<=nonEmptyD;j++)
    {
      fprintf(stdout, "\t%d", indicesD[j]);
    }
    fprintf(stdout, "\n");

    for (i=0;i<nonEmptyG;i++)
    {
      fprintf(stdout, "%d", indicesG[i]);
      for (j=0;j<nonEmptyD;j++)
      {
        fprintf(stdout, "\t%.3f", Dice[j][i]);
      }
      fprintf(stdout, "\n\n");
    }
  }


  if (imageIn->type == UCHAR)
  {
    for (z=0;z<dimz;z++)
    {
      for (y=0;y<dimy;y++)
      {
        vtfree(in[z][y]);
        in[z][y]=NULL;
      }
      vtfree(in[z]);
      in[z]=NULL;
    }
    vtfree(in);
    in=NULL;
  }

  vtfree(gHist);
  vtfree(dHist);
  for(i=0;i<=dMax;i++) {
    vtfree(tab[i]);
    tab[i]=NULL;
  }
  vtfree(tab);
  tab=NULL;
  gHist=NULL;
  dHist=NULL;
  return (1);
}

static int findInVec(int v, int *l, int n)
{
    int i;
    if(v<l[0] || v>l[n-1])
        return(-1);

    for (i=0;i<n;i++)
        if(v==l[i]) return(i);
    return(-1);
}

int VT_ComputePairDices(vt_image *imageIn, vt_image *imageExt, double *T_ext_in, int *labelsIn, int *labelsExt, int npairs,  double *dices)
{
    char *proc="VT_ComputeDices";
    /* P_ext = T_ext_in * P_in */
    int i,in,ext;
    double **Dices;
    int *indicesIn;
    int *indicesExt;
    int nonEmptyIn;
    int nonEmptyExt;

    if ( MT_ComputeDice3DTrsf( *imageIn, *imageExt, T_ext_in, &indicesIn,&indicesExt,&nonEmptyIn,&nonEmptyExt,&Dices) != 1 ) {
        fprintf(stderr, "%s : Error while computing Dice values\n", proc);
        return(0);
    }


    for (i=0 ; i<npairs ; i++)
    {
        in=findInVec(labelsIn[i],indicesIn,nonEmptyIn);
        if (in<0) continue;
        ext=findInVec(labelsExt[i],indicesExt,nonEmptyExt);
        if (ext<0) continue;
        dices[i]=Dices[ext][in];
        /*
        if(VT_DiceAB(imageIn, imageExt, T_ext_in, labelsIn[i], labelsExt[i], &(dices[i])) != 1)
        {
            fprintf(sterr, "%s : Error while computing Dice value between labels %d and %d\n", proc, labelsIn[i], labelsExt[i]);
            return(0);
        }
        */
    }

    return(1);
}
/*
int VT_DiceAB(vt_image imageIn, vt_image imageExt, double *T_ext_in, int labelIn, int labelExt, double *dice)
{
    unsigned char *inU8;
    unsigned char *extU8;
    unsigned short int *inU16;
    unsigned short int *extU16;


    return(1);
}
*/
