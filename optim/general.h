#ifndef __GENERAL__
#define __GENERAL__

/*#define MinDouble -1.7976931348623157e+308
#define MaxDouble +1.7976931348623157e+308*/
#define MinDouble -1.0e+99
#define MaxDouble +1.0e+99

#define sqr( x ) x*x

typedef double* p_double ;
typedef p_double* pp_double ;
typedef long* p_long ;
typedef char* p_char;

void initialize_randoms ( void ) ;
long random_long ( long olow , long ohigh ) ;
double random_double ( double olow , double ohigh ) ;
void errmsg ( p_char omsg ) ;
    
#endif // __GENERAL__

