/**
    class solution

    File:         solution.h
  
    Description:  Class to store any solution 
                  within optimization projects (SIFEL, moot)
  
    Author:       Matej Leps, leps@cml.fsv.cvut.cz, 03/2007
  
*/  

#if !defined ( __SOLUTION__ )
#define        __SOLUTION__

#if !defined ( NULL )
#define NULL 0
#endif

#include <vector>
#include <string.h>

#include "general.h"

class solution
{
    public:
        solution ( void )
        {
      id = 0 ;
      Dim = 0 ;
      Obj = 1 ;
      CH = NULL ;
      value = NULL ;
      ranks = NULL ;
      author = 0 ;
      created = 0 ;
      computed = 0 ;
        }
        ~solution ( void )
    {
      if( value ) delete [] value ;
      if ( Obj > 1 ) {
        delete [] ranks ;
      }
      if ( CH ) delete [] CH ;
    }
    double& operator [] (long i)
    {
      return (CH[i]);
    }
    solution& operator = ( solution& b )
    {
      long i ;
      id = b.id ;
      author = b.author ;
      created = b.created ;
      computed = b.computed ;
      Dim = b.Dim ;
      //memcpy ( CH, b.CH, Dim*sizeof ( double ) );
      //memcpy ( value, b.value, Obj*sizeof ( double ) );
      for( i=0 ; i<b.Dim ; i++)
         CH[i] = b.CH[i] ;
      for( i=0 ; i<b.Obj ; i++)
         value[i] = b.value[i] ;
         
      if ( Obj > 1 ) {
        memcpy ( ranks, b.ranks, Obj*sizeof ( double ) );
      }
      tmp = b.tmp ;
      return *this;
    }

    long id ;        
    long author ;   
    long created ;  
    long computed ;          
        long Dim ;                  ///< number of optimization variables
    long Obj ;                  ///< number of objective functions
    double *CH ;                ///< vector of decision variables
    double *value ;             ///< objective functions values
    double *ranks ;             ///< objective functions ranking

    double tmp ;             ///< stores sums of objective functions (for BSF and CERAF)
};

typedef std::vector<solution> CSolutionSet;

#endif
