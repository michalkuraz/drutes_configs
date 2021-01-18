#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <time.h>

#include "cmlfile.h"
#include "opt1.h"


#if defined ( __PARALLEL__ )
opt1::opt1 ( my_mpi *oM ) : objective_function ()
#else
opt1::opt1 ( void ) : objective_function ()
#endif
{
#if defined ( __PARALLEL__ )
  M=oM ;
#endif

    long i;
    FILE *g ;

    Dim = 4;

#if defined ( __MULTIOBJECTIVE__ )
    Obj = 3 ;
#else
    Obj = 1 ;
#endif

    Domain = new p_double [Dim];
    for ( i=0; i<Dim; i++ )
    {
        Domain[i] = new double [3];
        Domain[i][2] = 0.0;
    }
/*
  p1 \ in ( 0.0001, 0.5 )
  p2 \in ( 1.050, 2.10 )
  p3 \in ( 0.250, 0.90 )
  p4 \in ( 0.0050, 10.00 )
  p5 \in ( 0.000, 0.10 )
*/

    Domain[0][0]=1e-16;
    Domain[0][1]=1e-1;
    
    Domain[1][0]=1e-16;
    Domain[1][1]=0.8;
    
    Domain[2][0] = 1e-16;
    Domain[2][1] = 0.8;

    Domain[3][0] = 1e-20;
    Domain[3][1] = 1e-5;
 
    
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


    // cleaning output file                                                                                                                                                                                           
    g = fopen( "opt1.sav", "wt" );
    fclose( g );
 

}

opt1::~opt1 ( void ) 
{
   for ( int i=0; i<Dim; i++ )
    {
      if( Domain[i] ) delete [] Domain[i] ;
    }
   if( Domain ) delete [] Domain ;
   //if( optimum ) delete [] optimum ;
}

void opt1::evaluate ( solution* S, long SelectedSize )
{

    FILE *f, *of;

    double a=0.0 ;
    long i, k, j ;
    double tmp=0.0 ;

    f = fopen( "pars.in", "wt" );
    
    for ( i=0; i<SelectedSize; i++ )
    {
      objective_function::value ( S[i] ) ;
 

      fprintf( f, "p  ");

#if defined ( __DEBUG__ )
      printf( "%ld : ", i);
#endif

      for ( j=0; j<Dim; j++ ) {
	fprintf( f, " %10g", S[i].CH[j] );

#if defined ( __DEBUG__ )
	printf( " %10g", S[i].CH[j] );
#endif

      }

      fprintf( f, "\n");

#if defined ( __DEBUG__ )
      printf( "\n");
#endif
    }

    fprintf( f, "end\n" ) ;

#if defined ( __DEBUG__ )
    printf( "\n" ) ;
#endif

    fclose( f );

    run_bash( ) ;

    char out_file[maxstring] ; 
  
    for ( i=0; i<SelectedSize; i++ )
      {
	idx++ ;
	S[i].id = idx ;

	sprintf( out_file, "%ld/objfnc.val", i+1 ) ; 
	of = fopen( out_file, "rt" );

	if ( of==NULL )
	  {
	    fprintf ( stderr,"\n\n Error during opening file %s", out_file ) ;
	    fprintf ( stderr,"\n (%s, line %d).\n",__FILE__,__LINE__);
	    exit ( 1 ) ;
	  }

	printf ( "@@@ ID %ld : ", S[i].id ) ; 

	/* READING FROM FILE */

#if defined ( __MULTIOBJECTIVE__ )
	for ( j=0; j<Obj; j++ ) {

	  fscanf( of, "%lf\n", &a );
	  S[i].value[j] = -a ;

	  printf ( "%lf ", S[i].value[j]  ) ;

	}

	printf ( "for "  ) ;
#else

	fscanf( of, "%lf\n", &a );
	S[i].value[0] = -a ;

	printf ( "%lf for ", S[i].value[0]  ) ;

#endif

	for ( j=0; j<Dim; j++ ) {
	  printf( " %10g", S[i].CH[j] );
	}
	printf( "\n");

	fclose( of ) ;

      }

    f = fopen( "opt1.sav", "at" );
    
    for ( i=0; i<SelectedSize; i++ )
      {
	tmp = 0.0 ;

	fprintf( f, "OBJs:" ) ;
	for ( j=0; j<Obj; j++ ) {
	  tmp += S[i].value[j] ;
	  fprintf( f, " %10g", S[i].value[j] );
	}
	S[i].tmp = tmp ;
	fprintf( f, " TMP: %10g", tmp );

	fprintf( f, "X[]:" ) ;
	for ( j=0; j<Dim; j++ ) {
	  fprintf( f, " %10g", S[i].CH[j] );
	}

	fprintf( f, "\n" ) ;
      }
 
   fclose( f );
}

void opt1::run_bash ( void )
{
  
  char** arguments = new char*[ 2 ];
    
  arguments[0] = new char[maxstring] ;
  arguments[1] = NULL ;

  sprintf( arguments[0], "./bash.sh" ) ;
  
#if defined ( __DEBUG__ )
  printf( " I will execute: %s \n", arguments[0] );
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

  delete [] arguments[0] ;
  delete [] arguments ;

    
#if defined ( __DEBUG__ )
  printf( " PROGRAM is done. O.K. \n\n" ) ;
#endif

}

void opt1::value ( solution& S )
{
    FILE *f;
    char out_file[maxstring] ; 
    double a=0.0 ;
    objective_function::value ( S ) ;

#if defined ( __PARALLEL__ )
    idx = S.id;
#else
    idx++ ;
    S.id = idx ;
#endif

    clock_t cas ;
    cas = clock() ;
    S.id = cas ;
    
    //run_oofem( S ) ;

    //sprintf( out_file, "objfnc_%ld_%ld.val", S.created, S.id ) ; 
    //f = fopen( out_file, "rt" );
    //fscanf( f, "%lf", &a );
    
    //a = rand() ;
    a = S.CH[0] ;
    //fclose( f ); 

    printf ( "@@@ %ld : %ld %lf \n", S.created, S.id, a ) ;

#if defined ( __MULTIOBJECTIVE__ )

  S.value[0] = -a ;
  //S.value[1] = b ;
  //S.value[2] = c ;

#else

  S.value[0] = -a ;

#endif

  long j ;
  double tmp=0.0 ;

  f = fopen( "opt1.sav", "at" );
  for ( j=0; j<Dim; j++ ) {
    fprintf( f, " %10g", S[j] );
  }
  for ( j=0; j<Obj; j++ ) {
    //tmp += S.value[j] ;
    fprintf( f, " %10g", S.value[j] );
  }
  //S.tmp = tmp ;
  //fprintf( f, " %10g", tmp );
  fprintf( f, "\n" ) ;

  fclose( f );
}

void opt1::run_oofem ( solution& S )
{
  
  char** arguments = new char*[ 8 ];
    
  
  
  arguments[0] = new char[maxstring] ;
  arguments[1] = new char[maxstring] ;
  arguments[2] = new char[maxstring] ;
  arguments[3] = new char[maxstring] ;
  arguments[4] = new char[maxstring] ;
  arguments[5] = new char[maxstring] ;
  arguments[6] = new char[maxstring] ;
  arguments[7] = NULL ;

  sprintf( arguments[0], "./bash.sh" ) ;
  sprintf( arguments[1], "%ld_%ld", S.created, S.id );
  sprintf( arguments[2], "%f", S[0] );
  sprintf( arguments[3], "%f", S[1] );
  sprintf( arguments[4], "%f", S[2] );
  sprintf( arguments[5], "%f", S[3] );
  sprintf( arguments[6], "%f", S[4] );
  
#if defined ( __DEBUG__ )
  printf( " I will execute: %s %s %s %s %s %s %s \n", arguments[0], arguments[1], arguments[2], arguments[3], arguments[4], arguments[5], arguments[6]  );
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

  delete [] arguments[0] ;
  delete [] arguments[1] ;
  delete [] arguments[2] ;
  delete [] arguments[3] ;
  delete [] arguments[4] ;
  delete [] arguments[5] ;
  delete [] arguments[6] ;
  delete [] arguments ;

    
#if defined ( __DEBUG__ )
  printf( " PROGRAM is done. O.K. \n\n" ) ;
#endif

}

