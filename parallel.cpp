#include <mpi.h>

// -------------------------------------- MASTER ------------------------------------------------

void gartou::EVALUATE_GENERATION ( long start )
{
  long i, k ;

  long odkud ;
  long all, div, part ;
  long sum ;
  //long aid[max_processes] ; // slave id
  long position[max_processes] ;
  long amount[max_processes] ;
  
  all=ActualSize-(start*SelectedSize) ;
  part= all/M->size  ;
  div= all%M->size ;

  //
  // Index of beginning of aproppriate section of memmory
  sum=start*SelectedSize ;
  for ( i=( M->size )-1 ; i>-1 ; i-- )
    {
      amount[i]=part ;
      // adding one solution more for slaves if the total amount
      // is not dividable by no. processors
      if ( div>0 )
	{
	  amount[i]++ ;
	  div-- ;
	}
      position[i]=sum ;
      sum+=amount[i] ;
    }

  //
  // Sending

#if defined ( __PAR_DEBUG__ )
      printf( "Master: Starting sending data ... \n" );
      fflush( NULL ) ;
#endif

  // Starting from the most loaded one
  for( i=( M->size )-1 ; i>0 ; i-- )
    {
      //aid[i]=i ;
      MPI_Ssend( &M->mode, 1, MPI_LONG, i, tag_mode, MPI_COMM_WORLD ) ;
      MPI_Ssend( &amount[i], 1, MPI_LONG, i, tag_data, MPI_COMM_WORLD ) ;
      //MPI_Ssend( &position[i], 1, MPI_LONG, i, tag_data, MPI_COMM_WORLD ) ;
      for( k=0 ; k<amount[i] ; k++ )
	   MPI_Ssend( S[position[i]+k].CH , F->Dim, MPI_DOUBLE, i, tag_data, MPI_COMM_WORLD ) ;

#if defined ( __PAR_DEBUG__ )
      printf( "Master: Sent block from position %ld to %ld to node %ld \n", position[i], position[i]+k, i );
      fflush( NULL ) ;
#endif

    }

#if defined ( __PAR_DEBUG__ )
      printf( "Master: Sending done.\n" );
      fflush( NULL ) ;
#endif

  //
  // Computing
  EVALUATE_MASTER ( position[0] ) ;

#if defined ( __PAR_DEBUG__ )
      printf( "Master: Starting receiving data ... \n" );
      fflush( NULL ) ;
#endif
  //
  // Receiving
  // from any source
  for( i=1 ; i<M->size ; i++ )
    {
      MPI_Recv( &odkud, 1, MPI_LONG, MPI_ANY_SOURCE, tag_id, MPI_COMM_WORLD, &M->status ) ;

#if defined ( __PAR_DEBUG__ )
      printf( "Master: Receiving block from position %ld to %ld from node %ld \n", position[odkud], position[odkud]+amount[odkud], odkud );
      fflush( NULL ) ;
#endif


      for( k=0 ; k<amount[odkud] ; k++ )
	    MPI_Recv( S[ position[odkud]+k ].value,  F->Obj, MPI_DOUBLE, odkud, tag_value, MPI_COMM_WORLD, &M->status ) ;
    }  

#if defined ( __PAR_DEBUG__ )
      printf( "Master: Receiving done.\n" );
      fflush( NULL ) ;
#endif
}	

void gartou::EVALUATE_MASTER ( long start_id )
{
    long i;

    for ( i=start_id; i<ActualSize; i++ )
    {
      //S[i].id = i ;
      S[i].created = M->rank ;
      F->value ( S[i] );
      M->solved++;
    }
}	

void gartou::run_master ( my_mpi *oM )
{
  M=oM ;
  M->mode=tag_run ;

  run() ;

  //
  // MPI END
  M->mode=tag_end ;
  for( long i=1 ; i<M->size ; i++ ) 
    {
      MPI_Ssend( &M->mode, 1, MPI_LONG, i, tag_mode, MPI_COMM_WORLD ) ;
    }
}

// -------------------------------------- SLAVE ------------------------------------------------

void gartou::run_slave ( my_mpi *oM )
{
    M=oM ;
    CONFIGURATION_SLAVE () ;

    
    while ( to_continue_slave() )
      {
	EVALUATE_SLAVE() ;
      }
}

void gartou::CONFIGURATION_SLAVE ( void )
{
  if ( pool_rate > 0 ) {
    PoolSize = 2*pool_rate*F->Dim;
    SelectedSize = pool_rate*F->Dim;
  }
  ActualSize=0 ;

  S = new solution[PoolSize];

  for ( long i=0; i<PoolSize; i++ )
    {
      F->alloc( S[i] ) ;
    }
}

long gartou::to_continue_slave ( void )
{
  MPI_Recv( &M->mode, 1, MPI_LONG, 0, tag_mode, MPI_COMM_WORLD, &M->status ) ;
  return ( M->mode!=tag_end ) ;
}

void gartou::EVALUATE_SLAVE ( void )
{
  long i, j ;
  //long idx ;
  //long sid ;

#if defined ( __PAR_DEBUG__ )
      printf( "Slave %d: Starting receiving data ... \n",M->rank+1 );
      fflush( NULL ) ;
#endif

  //
  // Receiving
  MPI_Recv( &ActualSize, 1, MPI_LONG, 0, tag_data, MPI_COMM_WORLD, &M->status ) ;
  //MPI_Recv( &idx, 1, MPI_LONG, 0, tag_data, MPI_COMM_WORLD, &M->status ) ;
  for( i=0 ; i<ActualSize ; i++ )
    MPI_Recv( S[i].CH, F->Dim, MPI_DOUBLE, 0, tag_data, MPI_COMM_WORLD, &M->status ) ;

#if defined ( __PAR_DEBUG__ )
  printf( "Slave %d: Receiving done.\n",M->rank+1 );
  fflush( NULL ) ;

  for( i=0 ; i<ActualSize ; i++ ) {
    printf( "Slave %d: Received ",M->rank+1 );
    for ( j=0 ; j<F->Dim ; j++)
      printf( "%f ", S[i].CH[j] );
    printf( " \n" );
  }
  printf( " \n" );
      
  fflush( NULL ) ;
#endif

  //
  // Computing
  for ( i=0; i<ActualSize; i++ )
    {
      //S[i].id = idx+i ;
      S[i].created = M->rank; 
      F->value ( S[i] );
      M->solved++;
    }
  
#if defined ( __PAR_DEBUG__ )
      printf( "Slave %d: Starting sending data ... \n",M->rank+1 );
      fflush( NULL ) ;
#endif

  //
  // Sending
  MPI_Ssend( &M->rank, 1, MPI_LONG, 0, tag_id, MPI_COMM_WORLD ) ;

#if defined ( __PAR_DEBUG__ )
      printf( "Slave %d: My rank sent \n",M->rank+1 );
      fflush( NULL ) ;
      printf( "Slave %d: Sending ",M->rank+1 );
      for( i=0 ; i<ActualSize ; i++ )
	printf( "%f ", S[i].value[0] );

      printf( " \n" );

      fflush( NULL ) ;
#endif

  for( i=0 ; i<ActualSize ; i++ )
    MPI_Ssend( S[i].value, F->Obj, MPI_DOUBLE, 0, tag_value, MPI_COMM_WORLD ) ;

#if defined ( __PAR_DEBUG__ )
  printf( "Slave %d: Sending done.\n",M->rank+1 );
  fflush( NULL ) ;
#endif
}

