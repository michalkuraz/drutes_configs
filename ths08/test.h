#include "function.h"

#if defined ( __PARALLEL__ )
#include "my_mpi.h"
#endif

class test : public objective_function
{
    public:
#if defined ( __PARALLEL__ )
        test ( my_mpi *oM );
#else
        test ( void );
#endif

        ~test ( void );
	void value ( solution& S ) ;

    private:
	long idx ;

#if defined ( __PARALLEL__ )
	my_mpi *M ;
#endif

};

