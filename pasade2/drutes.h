#include "function.h"
#include "fradu.h"

#if defined ( __PARALLEL__ )
#include "my_mpi.h"
#endif

class drutes : public objective_function
{
    public:
#if defined ( __PARALLEL__ )
        drutes ( my_mpi *oM );
#else
        drutes ( void );
#endif

        ~drutes ( void );
	void value ( solution& S ) ;
    private:
	void run_oofem ( char *fname ) ;

    private:
	long idx ;
	fradu Grafy ;


#if defined ( __PARALLEL__ )
	my_mpi *M ;
#endif

};

