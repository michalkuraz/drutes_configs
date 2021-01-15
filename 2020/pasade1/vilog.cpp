#if !defined ( __GARTOU__ )
#include "gartou.h"
#endif

#include <stdio.h>

#if defined ( __FSWE__ )
#include "fswe.h"
#endif

#define _X 0
#define _Y 1

void gartou::init_vilog ( void )
{
    if ( !vilog ) return;

#if defined ( __HAVE_KRESLITKO__ )
    running=0;
    open_graphics ();
#endif
}

void gartou::close_vilog ( void )
{
    if ( !vilog ) return;

#if defined ( __HAVE_KRESLITKO__ )
    close_graphics() ;
#endif
}

long gartou::vilog_news ( long ceraf )
{
    if ( !vilog ) return 1;

    static double sbtg=MinDouble;

    if ( !ceraf )
    {
        printf ( "fc=%ld   CERAF\n", fitness_call );
        sbtg=MinDouble;
    }
    long i;
    FILE *f;
    //f = fopen ( "mata.out", "at" );
    f = stdout ;
    if ( (btg_value-sbtg) > F->precision )
    {
        if ( F->optimum )
        {
            fprintf( f, "fc=%ld   BTG=%.10g   BSF=%.10g   opt=%.8f   extremes:%ld\n", fitness_call, btg_value, bsf_value, *F->optimum, n_extreme ) ;
	    fprintf( f, "btg= " );
	    for ( i=0; i<F->Dim; i++ ) fprintf( f, "%g ",btg[i] );
	    fprintf( f, "\n" );
        }
        else 
	{
		fprintf( f, "fc=%ld   BTG=%.10g   BSF=%.10g   extremes:%ld\n", fitness_call, btg_value, bsf_value, n_extreme ) ;
		fprintf( f, "btg= " );
        	for ( i=0; i<F->Dim; i++ ) fprintf( f, "%g ", btg[i] );
		fprintf( f, "\n" );
	}
	//printf( "killed:%ld\n", killed );
        sbtg=btg_value;
    }
    //fclose( f );
    f = NULL ;
#if defined ( __HAVE_KRESLITKO__ )
    Ko->delete_items( 1 ) ;
//    Ko->delete_items( 11 ) ;
    Ko->redraw_items( 3 ) ;
    draw_chromos( S,ActualSize ) ;
#if defined ( __FSWE__ )
    draw_peaks();
#endif
    draw_cerafs( Extremes, n_extreme );
    if ( !ceraf ) running=0;
    if ( !until_next_step () ) return 0;
#endif

//    specialitka pro prubeh fcni hodnoty behem jednoho vypoctu  ...........
/*
    static long k=0;
    char file[32];
    FILE *f;
    k++;
    if ( k==10 )
    {
        sprintf( file, "statistc/%ld.out", fitness_call);
        f = fopen( file, "at" );
        fprintf( f, "%.8g\n", bsf_value );
        k=0;
    }
*/
// .........................................


    return 1;
}

void gartou::vilog_save ()
{
    if ( !save ) return;
    
    long i, j;
    FILE *f;
    
    f = fopen( "gartou.sav", "wt" );
    
    fprintf( f, "generation %ld\n", generation );
    for ( i=0; i<SelectedSize; i++ )
    {
	
	for ( j=0; j<F->Dim; j++ ) fprintf( f, "%g ", S[i].CH[j] );
	fprintf( f, "%g\n", Force[i] );
    }
    
    fclose( f );
}

void gartou::vilog_results ()
{
    if ( !vilog ) return;

    long i,j;

    printf ( "\n" );
    for  ( i=0; i<n_extreme; i++ )
    {
        printf ( "extreme%ld=( %.5g", i, Extremes[i].location[0] );
        for ( j=1; j<F->Dim; j++ ) printf ( ", %.5g", Extremes[i].location[j] );
        printf ( " ) value=%.10g \n", Extremes[i].value );
    }
    printf ( "\n" );

    if ( F->optimum )
    {
        printf( "fc=%ld   BSF=%.10g   bsf_birth=%ld   opt=%.8f\n", fitness_call, bsf_value, bsf_birth, *F->optimum ) ;
    }
    else printf( "fc=%ld   BSF=%.10g bsf_birth=%ld\n", fitness_call, bsf_value, bsf_birth ) ;

    printf ( "x=( %.5g",bsf[0] );
    for ( i=1; i<F->Dim; i++ ) printf ( ", %.5g", bsf[i] );
    printf ( " )\n" );
    printf ( "f=( %.5g",bsf.value[0] );
    for ( i=1; i<F->Obj; i++ ) printf ( ", %.5g", bsf.value[i] );
    printf ( " )\n" );

#if defined ( __HAVE_KRESLITKO__ )
    draw_top();
    running=0;
    if ( !until_next_step () ) return;
#endif

}

void gartou::vilog_make_out ()
{
    if ( !vilog ) return;

    long i;
    char file[32];

    FILE *f;

    if ( F->outfile )
    {
        sprintf ( file, "%s.out", F->outfile );
        f = fopen ( file, "at" );
    }
    else f = fopen ( "results", "at" );

    fprintf ( f, "%ld\t%.10g\t%ld\t%ld\t", fitness_call, bsf_value, bsf_birth, n_extreme );
    if ( F->optimum ) fprintf ( f, "%.5g\t", fabs(*F->optimum-bsf_value) );
    for ( i=0; i<F->Dim; i++ ) fprintf ( f, "%.5g\t", bsf[i] );
    fprintf ( f, "\n" );
    fclose ( f );
}

#if defined ( __HAVE_KRESLITKO__ )
void gartou::open_graphics ( void )
{
  double delta=10. ;
    Ko=new kreslitko ;
    Ko->create_button( "Close",CloseSignal ) ;
    Ko->create_button( "Next",NextSignal ) ;
    Ko->create_button( "Start/Stop",StartStopSignal ) ;
    Ko->create_button( "PostScript",PostScriptSignal ) ;
    Ko->open_graphics( F->Domain[_X][0]-delta,F->Domain[_Y][0]-delta,F->Domain[_X][1]+delta,F->Domain[_Y][1]+delta,700,700 ) ;
    Ko->set_background_color( 200,200,200 ) ;
    Ko->assign_rgb_color( 1,0,0,150 ) ;
    Ko->assign_rgb_color( 2,0,100,0 ) ;
    Ko->assign_rgb_color( 3,200,0,0 ) ;
    Ko->assign_rgb_color( 4,0,0,0 ) ;
    Ko->set_color( 1 ) ;
    Ko->set_text_color( 3 ) ;
    Ko->assign_font( 0,12,"Times",0 ) ;
    //Ko->x_scale( F->Domain[_Y][0],F->Domain[_X][0],F->Domain[_X][1],0,8 ) ;
    //Ko->y_scale( F->Domain[_X][0],F->Domain[_Y][0],F->Domain[_Y][1],0,8 ) ;
}

void gartou::draw_chromos ( solution &S , long on )
{
    double xave=0.0,yave=0.0 ;
    long i ;

    
    Ko->begin_serial() ;
    Ko->set_color( 4 ) ;
    //Ko->mark_absolute( 0.0,0.0,KRESLITKO_MarkRectangle,350,1 );
    Ko->set_color( 1 ) ;
    for ( i=0 ; i<on ; i++ )
    {
        Ko->mark_percentual( S[i][_X],S[i][_Y],KRESLITKO_MarkCross,3,1 ) ;
        xave+=S[i][_X] ;
        yave+=S[i][_Y] ;
    }
    xave/=( double )on ;
    yave/=( double )on ;
//    Ko->mark_percentual( xave,yave,KRESLITKO_MarkRectangle,15,1 ) ;
    Ko->flush_serial() ;
    Ko->is_ready() ;
}

void gartou::draw_cerafs ( extreme *Extremes, long on )
{
    long i;

    Ko->set_color( 3 ) ;
    Ko->set_pen_width ( 2 );
    Ko->assign_font( 1,9,"Times",1 ) ;

    char *value=new char[32];

    double movea,moveb,movec;

    Ko->begin_serial() ;
    for ( i=0 ; i<on ; i++ )
    {
        movea=1.0;
        moveb=1.0;
        movec=1.0;

        if ( Extremes[i].location[_Y]+10.0 > F->Domain[_Y][1] ) movea=-2.0;
        if ( Extremes[i].location[_X]+20.0 > F->Domain[_X][1] ) moveb=-2.0;
        if ( Extremes[i].location[_Y]-10.0 < F->Domain[_Y][0] ) movec=-2.0;


        //        Ko->ellipse ( Extremes[i].location[_X], Extremes[i].location[_Y], Extremes[i].radius[_X]/2, Extremes[i].radius[_Y]/2,1 );
        Ko->mark_absolute( Extremes[i].location[_X], Extremes[i].location[_Y], KRESLITKO_MarkCircle,long (Extremes[i].radius[_X]/(F->Domain[_X][1]-F->Domain[_X][0])*350.0),1 );
        Ko->mark_percentual( Extremes[i].location[_X], Extremes[i].location[_Y],KRESLITKO_MarkLyingCross,2,1 ) ;
        sprintf ( value, "%.4g", Extremes[i].value );
        Ko->text( Extremes[i].location[_X]+5.0*movea/moveb, Extremes[i].location[_Y]+5.0*movea, value, 1 );
        sprintf ( value, "*%ld", Extremes[i].birth );
        Ko->text( Extremes[i].location[_X]+5.0/moveb*movec, Extremes[i].location[_Y]-10.0/movec, value, 1 );
    }
    Ko->flush_serial() ;
    Ko->is_ready() ;
    Ko->set_pen_width ( 1 );

    delete [] value;
}

void gartou::draw_top ( void )
{
    char *value=new char[32];
    double movea=1.0,moveb=1.0,movec=1.0;

    if ( bsf[_Y]+10.0 > F->Domain[_Y][1] ) movea=-2.0;
    if ( bsf[_X]+20.0 > F->Domain[_X][1] ) moveb=-2.0;
    if ( bsf[_Y]-10.0 < F->Domain[_Y][0] ) movec=-2.0;

    Ko->set_color( 1 ) ;
    Ko->assign_font( 1,9,"Times",1 ) ;
    Ko->set_pen_width ( 2 );

    sprintf ( value, "%.4g", bsf_value );
    Ko->text( bsf[_X]+5.0*movea*moveb, bsf[_Y]+5.0*movea, value, 1 );
    sprintf ( value, "*%ld", long(bsf_birth/pool_rate/F->Dim )-1);
    Ko->text( bsf[_X]-20.0*movea*moveb, bsf[_Y]+5.0*movec, value, 1 );
}
#if defined ( __FSWE__ )
void gartou::draw_peaks ( void )
{
    long i;

    Ko->set_color( 2 ) ;
    Ko->set_pen_width ( 2 );
    Ko->assign_font( 1,9,"Times",1 ) ;

    char *value=new char[32];

    double movea,moveb;

    Ko->begin_serial() ;
    for ( i=0 ; i<Members ; i++ )
    {
        movea=1.0;
        moveb=1.0;

        if ( F->Peaks[i].x_0[1]+10.0 > F->Domain[1][1] ) movea=-2.0;
        if ( F->Peaks[i].x_0[0]+20.0 > F->Domain[0][1] ) moveb=-2.0;

        Ko->mark_percentual( F->Peaks[i].x_0[0], F->Peaks[i].x_0[1],KRESLITKO_MarkLyingCross,2,1 ) ;
        sprintf ( value, "%.4g", F->Peaks[i].value );
        Ko->text( F->Peaks[i].x_0[0]+5.0*movea/moveb, F->Peaks[i].x_0[1]+5.0*movea, value, 1 );
    }
    Ko->mark_percentual( F->Peaks[F->top].x_0[0], F->Peaks[F->top].x_0[1],KRESLITKO_MarkCross,3,1 ) ;
    Ko->flush_serial() ;
    Ko->is_ready() ;
    Ko->set_pen_width ( 1 );

    delete [] value;
}
#endif

long gartou::until_next_step ( void )
{
    FILE *f;
    long i;

    long signal ;
    if ( running )
    {
        if ( !Ko->receive_signal( signal )) return 1 ;
        if ( signal==StartStopSignal ) running=0 ;
    }
    if ( !running )
    {
        Ko->wait_and_receive_signal( signal ) ;
        switch ( signal )
        {
            case CloseSignal:
                return 0 ;
                break ;
            case NextSignal:
                return 1 ;
                break ;
            case StartStopSignal:
                f=fopen( "CHs.dat", "wt" );
                for( i=0; i<SelectedSize; i++ ) fprintf( f, "%g %g\n", S[i].CH[0], S[i].CH[1] );
                fclose( f );
                running=1 ;
                break ;
            case PostScriptSignal:
                char *psfile=new char[64];
                sprintf ( psfile, "gartou.ps" );
                Ko->export_as_ps( psfile, 8.0, 8.0 );
                printf( "filename <%s> exported as PostScript\n", psfile ) ;
                delete [] psfile;
                break;
        }
    }
    return 1 ;
}

void gartou::close_graphics ( void )
{
    long signal ;
    do
    {
        Ko->wait_and_receive_signal( signal ) ;
    }
    while ( signal!=CloseSignal ) ;
    if ( Ko ) delete Ko ;
}
#endif

