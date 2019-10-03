#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/types.h>

#include "cmlfile.h"
#include "test.h"

#if defined ( __PARALLEL__ )
test::test ( my_mpi *oM ) : objective_function ()
#else
test::test ( void ) : objective_function ()
#endif
{
#if defined ( __PARALLEL__ )
  M=oM ;
#endif

    long i;

    Dim = 8;
    Obj = 8;

    Domain = new p_double [Dim];
    for ( i=0; i<Dim; i++ )
    {
        Domain[i] = new double [3];
        Domain[i][2] = 0.0;
    }

    Domain[0][0]=20000.0; //e
    Domain[0][1]=50000.0;
    Domain[1][0]=0.1;     //nu
    Domain[1][1]=0.3;
    Domain[2][0]=3.0;     //k1
    Domain[2][1]=5.0;
    Domain[3][0]=100.0;   //k2
    Domain[3][1]=1000.0;
    Domain[4][0]=5.0;     //k3
    Domain[4][1]=15.0;
    Domain[5][0]=30.0;    //k4
    Domain[5][1]=200.0;
    Domain[6][0]=3.0;     //c3
    Domain[6][1]=5.0;
    //Domain[7][0]=0.2;     //c20
    //Domain[7][1]=5.0;
    Domain[7][0]=2.0;     //c20
    Domain[7][1]=50.0;
    
    optimum = new double[Obj] ;
    for ( i=0; i<Obj; i++ )
      optimum[i] = 0.0;

    cmlfile f( "parameters.cfg" ) ;
    
    f.set_labels( 1 ) ;
    f.set_label_string(  0,"TestPrecision" ) ;
     
    f.get_value( 0,precision ) ;
 
    f.close() ;

    Return_to_domain=1;

    idx = 0 ;
}

test::~test ( void ) 
{
}


void test::value ( solution& S )
{
  objective_function::value ( S ) ;

  double R[8], W[8] ;

  R[0] = S.CH[0];
  R[1] = S.CH[1];
  R[2] = pow( 10.0, -S.CH[2] );
  R[3] = S.CH[3];
  R[4] = S.CH[4];
  R[5] = S.CH[5];
  R[6] = S.CH[6];
  R[7] = S.CH[7]/10.0;

  W[0] = 32035.5 ;
  W[1] = 0.178759 ;
  W[2] = 0.0001069 ; //-log10(0.0001069) ;
  W[3] = 452.5 ;
  W[4] = 7.4167 ;
  W[5] = 48.417 ;
  W[6] = 4.4167 ;
  W[7] = 0.24 ;

  idx++ ;
  S.id = idx ;

  long i ;
  double tmp=0.0 ;

  //
  // For 1 objective
  //
  if ( Obj==1 ) {
    for ( i=0; i<8; i++ )
      {
	//	tmp += fabs( R[i]-W[i] ) ;
	tmp += ( R[i]-W[i] )*( R[i]-W[i] ) ;
      }
    tmp/= 8.0 ;
    
    S.value[0] =  - tmp ;
    S.computed = 1 ;
  }
  //
  // For 8 objective functions
  //
  else if ( Obj==8 ) {
    for ( i=0; i<8; i++ )
      {
	//	S.value[i]= -fabs( R[i]-W[i] )/W[i] ;
	S.value[i]= -( R[i]-W[i] )*( R[i]-W[i] )/(W[i]*W[i]) ;
      }
    S.computed = 1 ;
  }
  else {
      fprintf ( stderr,"\n\n No such number of objective functions" ) ;
      fprintf ( stderr,"\n in function test::false_enumerate (%s, line %d).\n", \
		__FILE__,__LINE__);
      exit ( 1 ) ;
  }    

}
