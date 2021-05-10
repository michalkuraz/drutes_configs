#include "general.h"

#include <sys/time.h>
#include <unistd.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//long rseed=-115889033;
//long rseed=59381828;
long rseed=-1;
//long rseed=-125780985;

void initialize_randoms ( void )
{
    if ( rseed==-1 )
    {
        struct timeval tv ;
        struct timezone tz ;
        gettimeofday( &tv,&tz ) ;
        rseed=tv.tv_usec+( tv.tv_sec%4194304 )*1000 ;

        FILE *f;
        f=fopen("rseed", "wt");
        fprintf ( f, "%ld", rseed );
        fclose ( f );
    }
    srand( rseed ) ;
}

long random_long ( long olow , long ohigh )
{
    return rand()%( ohigh-olow+1 )+olow ;
}

double random_double ( double olow , double ohigh )
{
    return (( double )rand()/RAND_MAX )*( ohigh-olow )+olow ;
}

void errmsg ( p_char omsg )
{
    printf( "\n %s \n" , omsg ) ;
    exit ( 1 )  ;
}

