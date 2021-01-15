#if !defined ( __GARTOU__ )
#define        __GARTOU__

#if !defined ( __OBJECTIVE_FUNCTION__ )
#include "function.h"
#endif

#include <math.h>
#include <string.h>
#include <stdlib.h>

//#include "matrix.h"

#if defined ( __HAVE_KRESLITKO__ )

#include <kreslitko.h>

#define CloseSignal 101
#define NextSignal 102
#define StartStopSignal 103
#define PostScriptSignal 104

#endif

#if defined ( __PARALLEL__ )
#include "my_mpi.h"
#endif

struct extreme
{
    int birth;
    double *location, value, *radius;
};

class gartou
{
    public:
        gartou ( void );
        ~gartou ( void );
#if defined ( __PARALLEL__ )
        void run_master ( my_mpi *oM ) ;
	void run_slave ( my_mpi *oM );
        void CONFIGURATION_SLAVE ( void ) ;
	void EVALUATE_SLAVE ( void ) ;
	void EVALUATE_MASTER ( long start_id ) ;
        long to_continue_slave ( void ) ;
    private:
	long mpi_mode ;
	my_mpi *M ;
    public:
#endif
        void run ( void );
	long get_calls ( void ) { return fitness_call; };
        long vilog, save, load;
        long vilog_gati;
        long fitness_calls_limit;
	char load_file[32];

    private:
// -------------------------------------- SADE TECHNOLOGY --------------------------------------------------
        //double* new_point ( void );
	void new_point ( double* p );
        void configuration ( void );
        void FIRST_GENERATION ( void);
        void MUTATE ( void );
        void LOCAL_MUTATE ( void );
        void CROSS ( void );
        void EVALUATE_GENERATION ( long start );
        void COMPUTE_FITNESS ( long start );
        void SELECT ( void );
// -------------------------------------- CERAF TECHNOLOGY -------------------------------------------------
        void init_ceraf ( void );
        long new_ceraf ( void );
        void check_ceraf_distance ( long start );
        void clear_ceraf ( void );
// -------------------------------------- STOPPING CRITERIA ------------------------------------------------
        long to_continue ( void );

//--------------------------------------- VISUALISATION AND LOGGING ----------------------------------------
        void init_vilog ( void );
        long vilog_news ( long ceraf );
	void vilog_save ( void );
        void vilog_results ();
        void vilog_make_out ();
        void close_vilog ( void );
#if defined ( __HAVE_KRESLITKO__ )
        kreslitko *Ko ;
        long running ;
        void open_graphics ( void );
        void draw_chromos ( matrix &CH, long on );
        void draw_cerafs ( extreme *Extremes, long on );
        void draw_top ( void );
        void draw_peaks ( void );
        long until_next_step ( void );
        void close_graphics ( void );
#endif

//------------------variables of SADE:
    public:
        long pool_rate;
        double radioactivity, local_radioactivity;
        double mutation_rate, mutagen_rate;
        double cross_rate;
#if !defined ( __FSWE__ )
        objective_function *F;
#endif
#if defined ( __FSWE__ )
#include "fswe.h"
        fswe *F;
#endif

    private:
        long ActualSize, generation, fitness_call, PoolSize, SelectedSize;
        double *Force, *mutagen;
        //p_double *CH;
	//matrix CH ;
        long gradient;
	solution *S ;
	//CSolutionSet S ;
	long *ranks ;

	solution bsf ;
	solution btg ;
        double bsf_value, btg_value;
	//double *btg;
	long btg_id;
        long bsf_birth;
//------------------variables of CERAF:
    public:
        long ExtremeSize, quiet;
        double radius_rate, deact_rate;

    private:
        long n_extreme, killed, monsters;
        extreme *Extremes;
        long id, lastchange;
        double lastbtg, quiet_generation;
};
#endif
