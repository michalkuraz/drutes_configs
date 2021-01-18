#if !defined ( __MY_MPI__ )
#define        __MY_MPI__

#if defined ( __PARALLEL__ )

#include <mpi.h>

#define tag_run     201
#define tag_end     202
#define tag_mode    203
#define tag_data    204
#define tag_id      205
#define tag_value   206

#define max_processes  30

class my_mpi
{
  public:
    my_mpi ( void ) { solved=0 ; mode=tag_end ; } ;
   ~my_mpi ( void ) {} ;
  public:
   int rank ;
   int size ;
   long solved ;
   long mode ;
   char uname[256];
   MPI_Status status ;
} ;

#endif

#endif
