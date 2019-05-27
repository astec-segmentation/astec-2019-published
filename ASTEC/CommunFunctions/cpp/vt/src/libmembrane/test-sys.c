#include <stdlib.h>
#include <stdio.h>
#include <string.h> 

#include <sys/types.h>
#include <sys/sysctl.h>
#include <sys/resource.h>
#include <sys/time.h>

#define   RUSAGE_SELF     0
#define   RUSAGE_CHILDREN     -1

char *intDescription[] = {
  "hw.activecpu",
  "hw.logicalcpu",
  "hw.logicalcpu_max",
  "hw.memsize",
  "hw.ncpu",
  "hw.pagesize",
  "hw.physicalcpu",
  "hw.physicalcpu_max",
  ""
};

int main( int argc, char *argv[] )
{
  int i;
  int val;
  size_t len = sizeof( int );
  struct rusage r;
  void *t;
  long int o;


  for ( i=0; strlen( intDescription[i] ) > 0; i++ ) {
    sysctlbyname( intDescription[i], &val, &len, NULL, 0);
    fprintf( stderr, "str[%d] ('%s') = %d \n", i, intDescription[i], val );
  }
  
  (void)getrusage( RUSAGE_SELF, &r );

  fprintf( stderr, "ru_maxrss = %ld\n", r.ru_maxrss );
  fprintf( stderr, "ru_ixrss = %ld\n", r.ru_ixrss );
  fprintf( stderr, "ru_idrss = %ld\n", r.ru_idrss );
  fprintf( stderr, "ru_isrss = %ld\n", r.ru_isrss );

  o = r.ru_maxrss;
  
  for ( i=0; i< 10 ; i++ ) {

   
    t = malloc( 1024 * 1024 );
  
    (void)getrusage( RUSAGE_SELF, &r );

    fprintf( stderr, "\n ... iteration %d : %ld - %ld = %ld \n", i, o, r.ru_maxrss, o-r.ru_maxrss );

    fprintf( stderr, "ru_maxrss = %ld\n", r.ru_maxrss );
    fprintf( stderr, "ru_ixrss = %ld\n", r.ru_ixrss );
    fprintf( stderr, "ru_idrss = %ld\n", r.ru_idrss );
    fprintf( stderr, "ru_isrss = %ld\n", r.ru_isrss );

    o = r.ru_maxrss;
  }
  
}
