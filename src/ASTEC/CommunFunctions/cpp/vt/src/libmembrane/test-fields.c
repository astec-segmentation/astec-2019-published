/*************************************************************************
 * test-fields.c -
 *
 * $Id: test-fields.c,v 2.0 2013/10/22 14:22:00 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 * 
 * CREATION DATE: 
 * 2013/06/26
 *
 * ADDITIONS, CHANGES
 *
 */

#include <transfo.h>
#include <vtmalloc.h>

#include <vt_common.h>

#include <mt_membrane3D.h>


typedef struct local_par {
  int niter; /*  == nangles dans TVmembrane.c */
  int nsticks; /*  nb stick fields par calcul de plate field */
  double scale;
  double zfact;
  vt_names names;
  enumTVmode mode;
  int cptStick;
  int cptPlate;
  int cptBall;
  int test;
  double r;
  double alpha;
} local_par;





/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );







static char *usage = "[filename] [-niter|-nangles %d] [-nsticks %d] [-cftv|-tvclassic] [-scale %lf] [-zfact %lf]\n\
\t [-stickonly | -ballonly | -plateonly] [-test [-rayon %lf -alpha %lf]] [-v] [-D] [-help]";

static char *detail = "\
\t si 'filename' est absent, on prendra stdout\n\
\t -[niter|nangles] %d : nombre d'iterations\n\t si cette option n'est pas presente, on fait une seule iteration\n\
\t -nsticks %d : nombre d'angles pour discretiser le plan normal au champ de tenseur plate\n\
\t -cftv : stick fields de type CFTV\n\
\t -tvclassic : stick fields de type tensor voting classique (medioni)\n\
\t -scale %lf : echelle du vote\n\
\t -zfact %lf : facteur de resolution selon l'axe des z\n\
\t -[stickonly|ballonly|plateonly] : ecriture des champs stick|ball|plate seulement\n\
\t -test : test fields\n\
\t -rayon : rayon fot test fields\n\
\t -alpha : angle for test fields\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
\t -help : help\n\
 $Revision: 1.1 $ $Date: 2013/10/22 11:28:56 $ $Author: gael $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;

  int i, n;

  mt_angles angles;
  int Nangles, Nsticks;
  double theta, phi, v[3];

  vt_3Dtensor *sfields, bfield, *pfields, *tfields;
  int hsf;
  int dimFields[3];
  double scale;
  double zfact;
  double theCoeff[3];
  enumTVmode mode;

  int cptStick = 0;
  int cptBall  = 0;
  int cptPlate = 0;

  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );


  /*--- debut du test ---*/
  
  cptStick = par.cptStick;
  cptPlate = par.cptPlate;
  cptBall = par.cptBall;
  fprintf(stdout,"Entering in MT_Compute3DAngles(&angles,%d)...\n", par.niter);

  MT_InitAngles( &angles );
  Nangles = MT_Compute3DAngles( &angles, par.niter );

  fprintf(stdout,"Number of computed spherical coordinates = %d\n", Nangles);

  for (i = 0; i<Nangles; i++)
  {
    theta=angles.angles[i][0];
    phi=angles.angles[i][1];
    SphericalAnglesToUnitVector( theta, phi, v);
    fprintf(stdout, "angles[%d][0]=%f\t angles[%d][1]=%f\tv[%d]=(%f ; %f ; %f)\n", i, theta, i, phi, i, v[0], v[1], v[2]);

  }
  
  /*  TEST FIELDS */
  if (par.test == 1)
  {
	dimFields[0]=dimFields[1]=dimFields[2]=(int)2*par.r+1;
	
    fprintf(stdout,"Entering in MT_Compute3DTestFields...\n");
    if (MT_Compute3DTestFields( &tfields, dimFields,
          &angles, par.r, par.alpha) != 1)
    {
      MT_FreeAngles( &angles );
      return (-1);
    }
    fprintf(stdout, "Fields computed.\n");

    /*--- ecriture des filtres obtenus ---*/
    for (n=0;n<Nangles;n++)
    {
      VT_Write3Dtensor(tfields+n);
    }
 
    /*--- liberation memoire ---*/
    for (n=0;n<Nangles;n++){
      VT_Free3Dtensor(tfields+n);
    }
    vtfree(tfields);
    tfields=NULL;
    MT_FreeAngles( &angles );
    return( 1 );
  }

  mode = par.mode;
  scale = par.scale;
  zfact = par.zfact;
  theCoeff[0] = theCoeff[1] = scale;
  theCoeff[2] = scale * zfact;
  hsf = ceil(sqrt(-pow(scale,2)*log(0.01)));
  dimFields[0] = dimFields[1] = dimFields[2] = 2*hsf+1;
  Nsticks = par.nsticks;

  fprintf(stdout, "Parameters: scale=%f\ttheCoeff=(%f ; %f ; %f)\thsf=%d\tdimFields=(%d, %d, %d)\n", scale, theCoeff[0], theCoeff[1], theCoeff[2], hsf, dimFields[0], dimFields[1], dimFields[2]);

  /*  STICK FIELDS */
  if(cptStick==1 || cptBall==1) {
  fprintf(stdout,"Entering in MT_Compute3DStickFields...\n");
  if (MT_Compute3DStickFields( &sfields, dimFields,
          &angles, theCoeff, mode) != 1)
  {
    MT_FreeAngles( &angles );
    return (-1);
  }
  }

  /*  BALL FIELD */
  if(cptBall==1) {
  fprintf(stdout,"Entering in MT_Compute3DBallFieldFromStick...\n");
  if (MT_Compute3DBallFieldFromStick( &bfield, sfields, dimFields,
          Nangles ) != 1)
  {
    for (n=0;n<Nangles;n++){
      VT_Free3Dtensor(&(sfields[n]));
    }
    vtfree(sfields);
    sfields=NULL;
    MT_FreeAngles( &angles );
    return (-1);
  }
  }
  /*  PLATE FIELDS */
  if(cptPlate==1) {
  fprintf(stdout,"Entering in MT_Compute3DPlateFieldsFromStick...\n");
  if (MT_Compute3DPlateFields(&pfields, dimFields,
      &angles, theCoeff, Nsticks, mode) != 1)
  {
    for (n=0;n<Nangles;n++){
      VT_Free3Dtensor(&(sfields[n]));
    }
    vtfree(sfields);
    sfields=NULL;
    VT_Free3Dtensor(&bfield);
    MT_FreeAngles( &angles );
    return (-1);
  }
  }

  fprintf(stdout, "Fields computed.\n");

  /*--- ecriture des filtres obtenus ---*/
  if(cptBall==1)
    VT_Write3Dtensor( &bfield );
  for (n=0;n<Nangles;n++)
  {
    if(cptStick==1)
      VT_Write3Dtensor(sfields+n);
    if(cptPlate==1)
      VT_Write3Dtensor(pfields+n);
  }

  /*--- liberation memoire ---*/
  for (n=0;n<Nangles;n++){
    if(cptStick==1 || cptBall==1)
    VT_Free3Dtensor(&(sfields[n]));
    if(cptPlate==1)
      VT_Free3Dtensor(pfields+n);
  }
  if(cptStick==1 || cptBall==1) {
    vtfree(sfields);
    sfields=NULL;
  }
  if(cptBall==1) {
    VT_Free3Dtensor(&bfield);
  }
  if(cptPlate==1) {
    vtfree(pfields);
    pfields=NULL;
  }
  MT_FreeAngles( &angles );



  return( 0 );
}








static void MT_Parse( int argc, 
		      char *argv[], 
		      local_par *par )
{
  int i, nb, status;
  char text[STRINGLENGTH];
  
  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) MT_ErrorParse("\n", 0 );
  
  /*--- lecture des parametres ---*/
  i = 1; nb = 0;
  while ( i < argc ) {
    if ( argv[i][0] == '-' ) {
      if ( argv[i][1] == '\0' ) {
	if ( nb == 0 ) {
	  /*--- standart output ---*/
	  strcpy( par->names.out, ">" );
	  nb += 1;
	}
      }
      /*--- arguments generaux ---*/
      else if ( strcmp ( argv[i], "-help" ) == 0 ) {
	MT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 ) {
	_VT_VERBOSE_ = 1;
	VT_IncrementVerboseInVtTube3D(  );
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
	_VT_DEBUG_ = 1;
	VT_IncrementDebugInVtTube3D(  );
      }

      /*--- arguments specifiques ---*/



      
      else if ( strcmp ( argv[i], "-niter" ) == 0 || strcmp ( argv[i], "-nangles" ) == 0  ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -niter...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->niter) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -niter...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-nsticks" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -nsticks...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->nsticks) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -nsticks...\n", 0 );
      }


      else if ( strcmp ( argv[i], "-scale" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -scale...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->scale) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -scale...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-zfact" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -zfact...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->zfact) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -zfact...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-cftv" ) == 0 ) {
        par->mode = CFTV;
      }

      else if ( strcmp ( argv[i], "-tvclassic" ) == 0 ) {
        par->mode = TVCLASSIC;
      }

      else if ( strcmp ( argv[i], "-stickonly" ) == 0 ) {
        par->cptBall = par->cptPlate = 0;
      }

      else if ( strcmp ( argv[i], "-ballonly" ) == 0 ) {
        par->cptStick = par->cptPlate = 0;
      }

      else if ( strcmp ( argv[i], "-plateonly" ) == 0 ) {
        par->cptStick = par->cptBall = 0;
      }

      else if ( strcmp ( argv[i], "-test" ) == 0 ) {
        par->test = 1;
      }

      else if ( strcmp ( argv[i], "-rayon" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -rayon...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->r) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -rayon...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-alpha" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -alpha...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->alpha) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -alpha...\n", 0 );
      }

      /*--- option inconnue ---*/
      else {
	sprintf(text,"unknown option %s\n",argv[i]);
	MT_ErrorParse(text, 0);
      }
    }
    /*--- saisie du nom de fichier ---*/
    else if ( argv[i][0] != 0 ) {
	strncpy( par->names.out, argv[i], STRINGLENGTH );  
	nb++;
      }
      else 
	MT_ErrorParse("too much file names when parsing\n", 0 );
    
    i += 1;
  }
  
  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb == 0) {
    strcpy( par->names.out, ">" );  /* standart output */
  }
 

}







static void MT_ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit( 1 );
}








static void MT_InitParam( local_par *par )
{
  VT_Names( &(par->names) );
  par->niter = 1;
  par->scale = 2.0;
  par->mode = TVCLASSIC;
  par->nsticks = 36;
  par->zfact = 1;
  par->cptStick=1;
  par->cptPlate=1;
  par->cptBall=1;
  par->test=0;
  par->r=100;
  par->alpha=0.1;
}
