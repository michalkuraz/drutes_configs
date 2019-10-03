#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/types.h>

#include "cmlfile.h"
#include "drutes.h"


#if defined ( __PARALLEL__ )
drutes::drutes ( my_mpi *oM ) : objective_function ()
#else
drutes::drutes ( void ) : objective_function ()
#endif
{
#if defined ( __PARALLEL__ )
  M=oM ;
#endif

    long i;
    FILE *g ;

    Dim = 8;

#if defined ( __MULTIOBJECTIVE__ )
    Obj = 2 ;
#else
    Obj = 1 ;
#endif

    Domain = new p_double [Dim];
    for ( i=0; i<Dim; i++ )
    {
        Domain[i] = new double [3];
        Domain[i][2] = 0.0;
    }

  //alpha   |    n   |    m       | theta_r    | theta_s  |    darcy's K_s |  alpha coupling term |  fracture ratio (w_f)

    Domain[0][0]=0.001;     //alpha
    Domain[0][1]=0.1;
    Domain[1][0]=1.5;     // n
    Domain[1][1]=7.0;
    Domain[2][0]=0.3;     // m
    Domain[2][1]=0.9;
    Domain[3][0]=0.0;     //theta_r
    Domain[3][1]=0.0001;
    Domain[4][0]=0.1; //theta_s
    Domain[4][1]=1.0;
    Domain[5][0]=1.0;   //darcy's
    Domain[5][1]=4850.0;
    Domain[6][0]=0.0;    //couple
    Domain[6][1]=0.5;
    Domain[7][0]=0.000001;     //fracture r
    Domain[7][1]=0.5;
  
    
    //optimum = new double[Obj] ;
    //for ( i=0; i<Obj; i++ )
    //  optimum[i] = 0.0;


    //precision = 1e-1;
    Return_to_domain=1;

//     cmlfile f( "parameters.cfg" ) ;
//     
//     f.set_labels( 1 ) ;
//     f.set_label_string(  0,"OsminaPrecision" ) ;
//      
//     f.get_value( 0,precision ) ;
//  
//     f.close() ;

    idx = 0 ;

    Grafy.load_exp_data ();

    // cleaning output file                                                                                                                                                                                           
    g = fopen( "drutes.sav", "wt" );
    fclose( g );
 

}

drutes::~drutes ( void ) 
{
   for ( int i=0; i<Dim; i++ )
    {
      if( Domain[i] ) delete [] Domain[i] ;
    }
   if( Domain ) delete [] Domain ;
   //if( optimum ) delete [] optimum ;
}


void drutes::value ( solution& S )
{
  double a, b, c ;
  objective_function::value ( S ) ;

  idx++ ;
  S.id = idx ;

  Grafy.create_oofem_input( S ) ;
  run_oofem( Grafy.in_file ) ;
  Grafy.load_data() ;

//
// Opravit MIGUEL
// Uperavit, co se bude optimalizovat (vahovani)
  a = 200.*Grafy.err_inter2(1) ;
  b = Grafy.err_start() ;
  c = Grafy.err_inter_horisontal(1) ;

  printf ( "@@@ %ld %lf %lf %lf %lf \n", idx, a, b, c, a + b + c ) ;

#if defined ( __MULTIOBJECTIVE__ )

  S.value[0] = a ;
  S.value[1] = b ;
  S.value[2] = c ;

#else

  S.value[0] = a + b + c ;

#endif

//   S.value[0] = Grafy.err_inter(1) ;
//   S.value[1] = Grafy.err_inter(2) ;
//   S.value[2] = Grafy.err_inter(3) ;
//   S.value[3] = Grafy.err_inter(4) ;

  Grafy.delete_data() ;

  FILE *f;
  long j ;
  double tmp=0.0 ;

  f = fopen( "drutes.sav", "at" );
  for ( j=0; j<Dim; j++ ) {
    fprintf( f, " %10g", S[j] );
  }
  for ( j=0; j<Obj; j++ ) {
    tmp += S.value[j] ;
    fprintf( f, " %10g", S.value[j] );
  }
  S.tmp = tmp ;
  fprintf( f, " %10g", tmp );
  fprintf( f, "\n" ) ;

  fclose( f );
}

void drutes::run_oofem ( char *fname )
{
  char** arguments = new char*[ 3 ];
    
  arguments[0] = new char[maxstring] ;
  arguments[1] = new char[maxstring] ;
  arguments[2] = NULL ;

//
// Opravit MIGUEL
// jmenu binarniho souboru (co se spousti)
  sprintf( arguments[0], "./startsim" ) ;
  sprintf( arguments[1], "%s", fname );
  
#if defined ( __DEBUG__ )
  printf( " I will execute: %s %s\n", arguments[0], arguments[1]  );
#endif

  pid_t cpid, w;
  int status;

  cpid = fork();
  if (cpid == -1) { perror("fork"); exit(EXIT_FAILURE); }

  if (cpid == 0) {            /* Code executed by child */
    //
    // Here, the outside program is run
    //
    execvp( arguments [0], arguments ) ; 
  } else {                    /* Code executed by parent */
    do {
      w = waitpid(cpid, &status, 0 );
      if (w == -1) { perror("waitpid"); exit(EXIT_FAILURE); }
      
      if (WIFEXITED(status)) {
#if defined ( __DEBUG__ )
	printf(" exited, status=%d\n", WEXITSTATUS(status));
#endif
      } else if (WIFSIGNALED(status)) {
	printf("!PROGRAM killed by signal %d\n", WTERMSIG(status));
      } else if (WIFSTOPPED(status)) {
	printf("!PROGRAM stopped by signal %d\n", WSTOPSIG(status));
      } 
    } while (!WIFEXITED(status) && !WIFSIGNALED(status));
  }

  delete arguments[0] ;
  delete arguments[1] ;
  delete [] arguments ;

    
#if defined ( __DEBUG__ )
  printf( " PROGRAM is done. O.K. \n\n" ) ;
#endif

}

