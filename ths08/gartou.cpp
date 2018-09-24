#include "vilog.cpp"

#if !defined ( __GARTOU__ )
#include "gartou.h"
#endif

#if defined ( __PARALLEL__ )
#include "parallel.cpp"
#endif

#include "cmlfile.h"

gartou::gartou ( void )
{
//------------------variables of SADE:

    //pool_rate = 1;
    //radioactivity = 0.1;
    //local_radioactivity = 0.1;
    //cross_rate = 0.3;
    //radioactivity = 0.3;
    //local_radioactivity = 0.00;
    //cross_rate = 0.2;

    //mutagen_rate=400;
    //mutation_rate=0.5;

    PoolSize=0;
    SelectedSize=0;

    cmlfile f( "parameters.cfg" ) ;
  
    f.set_labels( 9 ) ;
    f.set_label_string(  0,"SadePoolRate" ) ;
    f.set_label_string(  1,"SadeRadioactivity" ) ;
    f.set_label_string(  2,"SadeLocalRadioactivity" ) ;
    f.set_label_string(  3,"SadeGradient" ) ;
    f.set_label_string(  4,"SadeCrossRate" ) ;
    f.set_label_string(  5,"SadeMutationRate" ) ;
    f.set_label_string(  6,"SadeMutagenRate" ) ;
    f.set_label_string(  7,"SadeSelectedSize" ) ;
    f.set_label_string(  8,"SadePoolSize" ) ;
  
    f.get_value( 0,pool_rate ) ;
    f.get_value( 1,radioactivity ) ;
    f.get_value( 2,local_radioactivity ) ;
    f.get_value( 3,gradient ) ;
    f.get_value( 4,cross_rate ) ;
    f.get_value( 5,mutation_rate ) ;
    f.get_value( 6,mutagen_rate ) ;
    
    if ( pool_rate == 0 ) {
      f.get_value( 7,SelectedSize ) ;
      f.get_value( 8,PoolSize ) ;
      if ( PoolSize <= SelectedSize )
    {
      fprintf ( stderr, "\n PoolSize too small \n\n" ) ;
      exit(1) ;
    }
    }

    f.close() ;

    ranks = NULL ;


    ActualSize=0;
    generation=0;
    fitness_call=0;
    Force=NULL;
    mutagen=NULL;
    //CH=NULL;
    S=NULL ;
    //S.clear() ;
    //bsf=NULL;
    bsf_value=MinDouble;
    bsf_birth=0;
    //btg=NULL;
    btg_id=0;
    btg_value=MinDouble;

//------------------constants and variables of CERAF:

    deact_rate  = 0.995;
    //quiet = 1700000;

    quiet = 600;

    ExtremeSize = 1048576;
    radius_rate = 0.75;

    n_extreme=0;
    killed=0;
    monsters=0;
    Extremes=NULL;

    id=0;
    lastchange=0;
    lastbtg=MinDouble;
    quiet_generation=0.0;

    vilog=1;
    save=0;
    load=0;
}

gartou::~gartou ( void )
{
  //long i;
  //for ( i=0; i<PoolSize; i++ ) if ( CH[i] ) delete [] CH[i];
  //if ( CH ) delete [] CH;
  
  if ( S ) delete [] S ;
  //S.clear() ;
  if ( Force ) delete [] Force;
  //if ( bsf ) delete [] bsf;
  if ( mutagen ) delete [] mutagen;
  if ( ranks ) delete [] ranks ;
}

// -------------------------------------- SADE TECHNOLOGY --------------------------------------------------

//double* gartou::new_point ( void )
void gartou::new_point ( double* p )
{
    long i;

    for ( i=0; i<F->Dim; i++ )
    {
        p[i] = random_double ( F->Domain[i][0], F->Domain[i][1] );
    }
}

void gartou::configuration ( void )
{
    if ( pool_rate > 0 ) {
      PoolSize = 2*pool_rate*F->Dim;
      SelectedSize = pool_rate*F->Dim;
    }
    fprintf ( stderr, "\n SADE algorithm is used: \n\n" ) ;
    fprintf ( stderr, " pool_rate: %ld \n", pool_rate ) ;
    fprintf ( stderr, " SelectedSize: %ld \n", SelectedSize ) ;
    fprintf ( stderr, " PoolSize: %ld \n", PoolSize ) ;
    fprintf ( stderr, " radioactivity: %f \n", radioactivity ) ;
    fprintf ( stderr, " local_radioactivity: %f \n", local_radioactivity ) ;
    fprintf ( stderr, " gradient: %ld \n", gradient ) ;
    fprintf ( stderr, " cross_rate: %f \n", cross_rate ) ;
    fprintf ( stderr, " mutation_rate: %f \n", mutation_rate ) ;
    fprintf ( stderr, " mutagen_rate: %f \n\n", mutagen_rate ) ;

    long i;
    mutagen = new double [F->Dim];
    for ( i=0; i<F->Dim; i++ )
    {
        mutagen[i] = (F->Domain[i][1]-F->Domain[i][0])/mutagen_rate;
        if ( F->Domain[i][2]>mutagen[i] ) mutagen[i] = F->Domain[i][2];
    }

    
    if ( pool_rate > 0 ) {
      quiet_generation=quiet/pool_rate;
    }
    else
      quiet_generation=quiet ;
}

void gartou::FIRST_GENERATION ( void)
{
    long i,j;

    Force = new double [PoolSize];
    //CH = new p_double [PoolSize];
    //allocm ( PoolSize, F->Dim, CH ) ;
    S = new solution[PoolSize];
    //S.resize( PoolSize ) ;
    if ( F->Obj > 1 )
      ranks = new long [PoolSize] ;

    for ( i=0; i<PoolSize; i++ )
    {
      F->alloc( S[i] ) ;
      //CH[i] = new_point ();
      new_point ( S[i].CH ) ;
    }

    F->alloc( bsf ) ;    
    F->alloc( btg ) ;
    
    
    FILE *f;
    if( load )
    {
        ActualSize = SelectedSize;
        f=fopen( load_file, "rt" ); 
        fscanf( f, "generation %ld", &i );
    printf( "nactena %ld. generace\n", i );
    for ( i=0; i<SelectedSize; i++ )
    {
        for ( j=0; j<F->Dim; j++ ) fscanf( f, "%lf", &S[i].CH[j] );
        fscanf( f, "%lf", &Force[i] );
    }
    fclose( f );    
    }
    else
    {
    ActualSize = PoolSize;
    EVALUATE_GENERATION ( 0 );
    COMPUTE_FITNESS ( 0 );
    SELECT ();
    }
}   

void gartou::MUTATE ( void )
{
    double p, *x;
    long i, j, index;

    for ( i=0; i<SelectedSize; i++ )
    {
        if ( ActualSize == PoolSize ) break;
        p = random_double ( 0,1 );
        if ( p<=radioactivity )
        {
            index = random_long ( 0,SelectedSize-1 );
            mutation_rate=random_double ( 0,1 );
        // x = new_point() ;
            x = new double[F->Dim] ;
        new_point ( x ) ;
            for ( j=0;  j<F->Dim; j++ )
            {
                S[ActualSize].CH[j] = S[index].CH[j]+mutation_rate*( x[j]-S[index].CH[j] );
            }
        S[ActualSize].created = 0 ;
        S[ActualSize].computed = 0 ;
            delete [] x;
            ActualSize++;
        }
    }
}

void gartou::LOCAL_MUTATE ( void )
{
    double p,dCH ;
    long i,j,index ;
    for ( i=0; i<SelectedSize; i++ )
    {
        if ( ActualSize == PoolSize ) break;
        p = random_double( 0,1 );
        if ( p <= local_radioactivity )
        {
            index = random_long( 0, SelectedSize-1 );

            for ( j=0; j<F->Dim; j++ )
            {
                dCH = random_double( -mutagen[j],mutagen[j] );
                S[ActualSize].CH[j] = S[index].CH[j]+dCH;
            }
        S[ActualSize].created = 0 ;
        S[ActualSize].computed = 0 ;
            ActualSize++;
        }
    }
}

void gartou::CROSS ( void )
{
    long i1,i2,i3,j ;
    long direction = 1;
    double rate ;
    while ( ActualSize < PoolSize )
    {
        i1 = random_long( 0,SelectedSize-1 );
        i2 = random_long( 1,SelectedSize-1 );
        if ( i1==i2 ) i2--;
        i3 = random_long( 0,SelectedSize-1 );
        if ( gradient && (Force[i1]>Force[i2])) direction=-1;
    //rate =  random_double( 0.0, cross_rate ) ;
    rate =  cross_rate  ;

        for ( j=0 ; j<F->Dim ; j++ )
        {
            S[ActualSize].CH[j] = S[i3].CH[j]+rate*(double)(direction)*( S[i2].CH[j]-S[i1].CH[j] ) ;
        }
    S[ActualSize].created = 0 ;
    S[ActualSize].computed = 0 ;
        ActualSize++;
    }
}

#if !defined ( __PARALLEL__ )
void gartou::EVALUATE_GENERATION ( long start )
{
    solution *P ;
    
    if ( start ==0 )
    {

      for ( long i=0 ; i<PoolSize ; i++ )
    {
      S[i].id = i ;
    }

        P = &S[0] ; 
        F->evaluate ( P, SelectedSize );
        P = &S[SelectedSize] ;
        F->evaluate ( P, SelectedSize );
    }   
    else
    {
        P = &S[SelectedSize] ;

    for ( long i=0 ; i<SelectedSize ; i++ )
      {
        P[i].id = i ;
      }


        F->evaluate ( P, SelectedSize );
    }
        
    
//    for ( long i=start*SelectedSize; i<ActualSize; i++ )
//    {
//  F->value ( S[i] );
//    }
}
#endif

void gartou::COMPUTE_FITNESS ( long start )
{
    long i, k ; 
    // Btg and Bsf are useless for MO optimization but not fo SO! 
    btg_value = MinDouble ;                                                                                                                                                                              
    //bsf_value = MinDouble ;                                                                                                                                                                              

    for ( i=start*SelectedSize; i<ActualSize; i++ )
    {
    if ( F->Obj < 2 ) 
      {
        Force[i] = S[i].value[0] ;
        // For comparison with btg and bsf
        S[i].tmp = Force[i] ;
      }
    else
      {
        S[i].tmp = 0.0 ;
        for ( k=0; k<F->Obj; k++ ) {
          S[i].tmp += S[i].value[k] ;
        }
      }

        if ( vilog_gati ) printf ( "%ld\n", fitness_call );
        fitness_call++;
    }

    //
    // Average Ranking algorithm
    //
    if ( F->Obj > 1 ) {
      long max, tmp ;
      long idx ;
      //
      // Btg and Bsf are useless here! Maybe not! :o)
      //
      //btg_value = MinDouble ;
      //bsf_value = MinDouble ;

      //
      // Attempt for bubble-sort
      // 1 is the best
      //
      for ( k = 0 ; k<F->Obj ; k++ ) {
    idx = 0 ;
    //for ( i=start*SelectedSize; i<ActualSize; i++ ) {
    for ( i=0 ; i<ActualSize; i++ ) {
      ranks[i] = idx ;
      idx++ ;
    }
    for ( i=0; i<ActualSize-1; i++ ) {
      max = i ;
      for ( long j=i+1 ; j<ActualSize ; j++ ) {
        if ( S[ranks[max]].value[k] < S[ranks[j]].value[k] )
          max = j ;
      }
      tmp = ranks[i] ;
      ranks[i] = ranks[max] ;
      ranks[max] = tmp ;
    }
    idx = 1 ;
    for ( i=0; i<ActualSize; i++ ) {
      S[ranks[i]].ranks[k] = idx ;
      idx++ ;
    }
      }
      //
      // AR Fitness
      //
      for ( i=0; i<ActualSize; i++ ) {
    Force[i] = 0.0 ;
    for ( k = 0 ; k<F->Obj ; k++ ) {
      //
      // Count only objective functions above precision 
      //
      //if ( fabs(F->optimum[k]-S[i].value[k]) > F->precision ) {
        // - for maximization
        Force[i] -= (double) S[i].ranks[k] ;
        //}
    }
    Force[i] /= (double) F->Obj ;
      }
    }

    for ( i=start*SelectedSize; i<ActualSize; i++ )
    {
      // Not know how to measure "best-ness" - sum of all objectives used instead
      if ( S[i].tmp>btg_value )
    {
      btg_value = S[i].tmp;
      btg=S[i];            
      btg_id=i;
        }
    }

    if ( btg_value>bsf_value )
    {
      //memcpy ( bsf, S[btg].CH, F->Dim*sizeof ( double ) );
      bsf = btg ;
      bsf_value = btg_value;
      bsf_birth = fitness_call;
      //printf (" ** %12.8f %6ld \n", bsf_value, generation) ;
      printf (" ++ %8ld OBJs: ", bsf.id ) ;
      fprintf ( stderr, " ++ %8ld OBJs:", bsf.id ) ;
      for ( long m=0; m<F->Obj; m++ ){
    printf (" %8.5g", bsf.value[m] ) ;
    fprintf ( stderr, " %8.5g", bsf.value[m] ) ;
      }
      printf (" TMP: %8.5g", bsf.tmp ) ;
      fprintf ( stderr, " TMP: %8.5g \n", bsf.tmp ) ;

      printf( "X =  " ) ; 

      for ( long k=0; k<F->Dim; k++ ){
    printf (" %8.5g", bsf[k] ) ;
    fprintf ( stderr, " %8.5g", bsf[k] ) ;
      }
      printf( " \n" ) ; 
    
    }
}   

void gartou::SELECT ( void )
{
    long i1, i2, dead, last;

    while ( ActualSize > SelectedSize )
    {
        i1 = random_long ( 0,ActualSize-1 );
        i2 = random_long ( 1,ActualSize-1 );
        if ( i1==i2 ) i2--;
        if ( Force[i1] >= Force[i2] ) dead = i2;
        else dead = i1;
        last = ActualSize-1;
    
        //h = CH[last];
        //CH[last] = CH[dead];
        //CH[dead] = h;
    //for ( long i=0 ; i<F->Dim ; i++ )
    //  {
    //    CH[dead][i]=CH[last][i] ;
    //  }
    S[dead] = S[last] ;
    if ( btg_id==last ) 
      {
        btg_id=dead;
      }
        Force[dead] = Force[last];
        ActualSize--;
    }
}   

// -------------------------------------- CERAF TECHNOLOGY -------------------------------------------------

void gartou::init_ceraf ( void )
{
    Extremes = new extreme [ExtremeSize];
}

long gartou::new_ceraf ( void )
{
    long i;

    if ( btg_value-lastbtg > F->precision/2 )
    {
        lastchange=generation;
        lastbtg = btg_value;
    }

    if ( ( generation-lastchange ) < quiet_generation ) return 1;
    Extremes[id].location = new double [F->Dim];
    Extremes[id].radius = new double [F->Dim];
    memcpy ( Extremes[id].location, btg.CH, F->Dim*sizeof ( double ) );
    Extremes[id].value = btg_value;
    Extremes[id].birth = generation;

    for ( i=0; i<F->Dim; i++ )
    {
        Extremes[id].radius[i]=radius_rate*(F->Domain[i][1]-F->Domain[i][0]);
    }
    btg_value=MinDouble;
    lastbtg=MinDouble;
    lastchange=generation;
    id++;
    n_extreme++;

    for ( i=0; i<SelectedSize; i++ )
      {
    //CH[i]=new_point();
    new_point( S[i].CH ) ;
      }

    return 0;
}

void gartou::check_ceraf_distance ( long start )
{
    long i, j, k;
    double r;

    for ( i=start*SelectedSize; i<PoolSize; i++ )
    {
        for ( j=0; j<n_extreme; j++ )
        {
            r=0.0;
            for ( k=0; k<F->Dim; k++ )
            {
                r += pow ( 2*( S[i].CH[k]-Extremes[j].location[k] )/Extremes[j].radius[k], 2 );
            }
            if ( r < 1 )
            {
                if ( i-SelectedSize > monsters )
                {
                    for ( k=0; k<F->Dim; k++ ) Extremes[j].radius[k]*=deact_rate;
                }
                //delete [] CH[i];
                //CH[i] = new_point ();
        new_point( S[i].CH ) ;
                j=n_extreme;
                i--;
                killed++;
            }
        }
    }
}

void gartou::clear_ceraf ( void )
{
    long i;
    for ( i=0; i<n_extreme; i++ )
    {
        if ( Extremes[i].location ) delete [] Extremes[i].location;
        if ( Extremes[i].radius ) delete [] Extremes[i].radius;
    }
    if ( Extremes ) delete [] Extremes;
}


// -------------------------------------- STOPPING CRITERIA ------------------------------------------------

long gartou::to_continue ( void )
{
  //  long status ;
  long good = 1 ;

  if ( fitness_call > fitness_calls_limit ) return 0;
  if ( F->optimum ) {
    for( long i=0 ; i<F->Obj ; i++ ) {
      if ( fabs(F->optimum[i]-bsf.value[i])>F->precision ) {
    good = 0 ;
    break ;
      }
    }
    if ( good ) {
      /*    status = F->update() ;
    if ( status ) {
    
    #if defined ( __DEBUG__ )
    printf( "\n Re-evaluating population ...\n" );
    #endif
    for ( long i=0 ; i<PoolSize ; i++ )
    S[i].computed = 0 ;
    
    btg_value=MinDouble;
    bsf_value=MinDouble;
    
    EVALUATE_GENERATION ( 0 );
    return 1 ;
    }
    else */
      return 0 ;
    }
  }
  return 1;
}

// -------------------------------------- MAIN -------------------------------------------------------------

void gartou::run ( void )
{

    long stop=0, ceraf=1;
    static long save_turn=0;
    configuration ();
    init_vilog ();
    init_ceraf ();

    FIRST_GENERATION ();

    vilog_save ();

    while ( !stop )
    {
        MUTATE ();
        monsters=ActualSize-SelectedSize;
        LOCAL_MUTATE ();
        CROSS ();
        check_ceraf_distance ( ceraf );
        EVALUATE_GENERATION ( ceraf );
        COMPUTE_FITNESS ( ceraf );
    SELECT ();
        ceraf = new_ceraf ();
        if ( !vilog_news ( ceraf ) ) stop=1;
        if ( !to_continue () ) stop=1;
        generation++;
        //fprintf( stderr, "generation %ld\n", generation );
    if ( save_turn > 0 ) 
      {
      vilog_save ();
      save_turn = 0;
      }
    save_turn++;
    }
    vilog_results ();
    vilog_make_out ();

    //F->evaluate ( bsf );

    close_vilog ();
    clear_ceraf ();
}
