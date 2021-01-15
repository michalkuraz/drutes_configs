#if !defined ( __OBJECTIVE_FUNCTION__ )
#define        __OBJECTIVE_FUNCTION__

#if !defined ( NULL )
#define NULL 0
#endif

#include <vector>
#include <math.h>

#include "general.h"
#include "solution.h"

struct CDomainFormat
{
	CDomainFormat(double Min,double Max,double Step)
	{
		d[0] = Min;
		d[1] = Max;
		d[2] = Step;
	}
	double& operator [] (long i)
	{
	  return (d[i]) ;
	}
	double d[3];
};

typedef std::vector<CDomainFormat> CDomain;

class objective_function
{
    public:
        objective_function ( void )
        {
	  Dim=0;
	  Domain=NULL;
	  optimum=NULL;
	  outfile=NULL;
	  precision=0.0;
	  Return_to_domain=0;
        } 
        virtual ~objective_function ( void )
        {
	  long i;
	  for ( i=0; i<Dim; i++ ) if( Domain[i] ) delete [] Domain[i];
	  //if( Domain ) delete [] Domain;
	  //if( optimum ) delete [] optimum;
        } 
        virtual void value ( solution& oS ) 
	{
	  return_to_domain( oS ) ;
	}
        virtual void evaluate ( solution*, long ) {}
	virtual long update ( void )
	{
	  return 0 ;
	}
	virtual void alloc ( solution& oS )
	{
	  oS.Dim = Dim ;
	  oS.Obj = Obj ;
	  oS.value = new double[Obj] ;
	  if ( Obj > 1 ) {
	    oS.ranks = new double[Obj] ;
	  }
	  oS.CH = new double [Dim] ;
	}
	void return_to_domain ( solution& oS )
	{
	  long k ;

	  for ( k=0; k<Dim; k++ )
	    {
	      if ( Domain[k][2] != 0.0 ) 
		{
		  oS.CH[k] = floor ( oS.CH[k]/Domain[k][2]+0.5 )*Domain[k][2];
		}
	    }

	  if ( Return_to_domain )
	    {
	      for ( k=0; k<Dim; k++ )
		{
		  if ( oS.CH[k]<Domain[k][0] ) oS.CH[k]=Domain[k][0];
		  if ( oS.CH[k]>Domain[k][1] ) oS.CH[k]=Domain[k][1];
		}
	    }
	}

        long Dim ;                  ///< number of optimization variables
	long Obj ;                  ///< number of objective functions
	p_double* Domain;
	//std::vector<CDomainFormat> Domain;
        double *optimum;
        double precision;
        char *outfile;
        char name[256];
        int Return_to_domain;
};

#endif
