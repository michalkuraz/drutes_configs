//
//  File:         logo.c
//
//  Description:  Class for logging output
//
//  Author:       Matej Leps
//

#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>

#include "logo.h"

const long KeywordLength=64 ;

struct logo_keyword
{
  logo_keyword ( void ) ;
 ~logo_keyword ( void ) ;
  char *word ;
  long binary ;
  long mode ;
  bool init ;
} ;

logo_keyword::logo_keyword ( void )
{
  binary=0 ;
  mode=0 ;
  word=new char[KeywordLength] ;
  init=false ;
}

logo_keyword::~logo_keyword ( void )
{
  if ( word )
    {
      delete [] word ;
    }
}

// -------------------------------------- C_LOGO --------------------------------------------------

logo::logo ( void ) 
{
  f=NULL ;
  printlevel=0 ;
  savelevel=0 ;
  files=0 ;
  levels=0 ;
  names=NULL ;
  filenames=NULL ;
}

logo::~logo ( void )
{
  for ( long i=0 ; i<files ; i++ )
    {
      if ( filenames[i].init ) fclose ( f[i] ) ;
    }
  if ( f ) delete [] f ;
  if ( names ) delete [] names ;
  if ( filenames ) delete [] filenames ;
}

void logo::set_levels ( long olevels )
{
  levels=olevels ;
  names=new logo_keyword[levels] ;
}

void logo::set_level_name ( long olevel, char *olevel_name )
{
  strcpy( names[olevel].word,olevel_name ) ;
  names[olevel].binary=give_binary( olevel+1 ) ;
}

void logo::set_files ( long ofiles )
{
  files=ofiles ;
  filenames=new logo_keyword[files] ;
  f=new PFILE[files] ;
}

void logo::set_file_name ( long ofile, char *ofile_name, long omode )
{
  strcpy( filenames[ofile].word,ofile_name ) ;
  filenames[ofile].mode=omode ;
}

void logo::load_files ( void )
{
  if ( savelevel > 0 )
    {
      for ( long i=0 ; i<files ; i++ )
	{
	  switch ( filenames[i].mode )
	    {
	    case mode_always_ask:
	      safe_load( i ) ;
	      break ;
	      
	    case mode_append:
	      append_load( i ) ;
	      break ;
	      
	    case mode_rewrite:
	      rewrite_load( i ) ;
	      break ;
	      
	    default:
	      fprintf ( stderr, "\n\n unknown open_file mode in function" ) ;
	      fprintf ( stderr, "\n logo::load_files (file %s, line %d).\n",__FILE__,__LINE__ ) ;
	      exit ( 1 ) ;
	    }
	}
    }
}

void logo::safe_load ( long ofile )
{
  int c ;
  if (( f[ofile]=fopen ( filenames[ofile].word,"rt" ) ) != NULL )
    {
      printf ( "File %s exists, rewrite, append, exit? [R/A/E]:", filenames[ofile].word ) ;
      c=getchar() ;
      while (getchar() != '\n')
	;
      if ( fclose( f[ofile] ) == EOF )
	{
	  f[ofile]=NULL ;
	  fprintf ( stderr,"\n\n Error during closing file %s", filenames[ofile].word ) ;
	  fprintf ( stderr,"\n in function logo::safe_load (%s, line %d).\n",__FILE__,__LINE__);
	  exit ( 1 ) ;
	}
      if ( !( c=='r' || c=='R' || c=='a' || c=='A'))
	{
	  f[ofile]=NULL ;
	  fprintf ( stderr,"\nUser exit.\n" ) ;
	  exit( 0 ) ;
	}
      else if ( c=='a' || c=='A' )
	{
	  f[ofile]=fopen ( filenames[ofile].word,"at" ) ;
	  if ( f[ofile]==NULL )
	    {
	      fprintf ( stderr,"\n\n Error during opening file %s", filenames[ofile].word ) ;
	      fprintf ( stderr,"\n in function logo::safe_load (%s, line %d).\n",__FILE__,__LINE__);
	      exit ( 1 ) ;
	    }
	  return ;
	}
    }

  f[ofile]=fopen ( filenames[ofile].word,"wt" ) ;
  if ( f[ofile]==NULL )
    {
      fprintf ( stderr,"\n\n Error during opening file %s", filenames[ofile].word ) ;
      fprintf ( stderr,"\n in function logo::safe_load (%s, line %d).\n",__FILE__,__LINE__);
      exit ( 1 ) ;
    }
  filenames[ofile].init=true ;
}

void logo::rewrite_load ( long ofile )
{
  f[ofile]=fopen ( filenames[ofile].word,"wt" ) ;
  if ( f[ofile]==NULL )
    {
      fprintf ( stderr,"\n\n Error during opening file %s", filenames[ofile].word ) ;
      fprintf ( stderr,"\n in function logo::rewrite_load (%s, line %d).\n",__FILE__,__LINE__);
      exit ( 1 ) ;
    }
  filenames[ofile].init=true ;
}

void logo::append_load ( long ofile )
{
  f[ofile]=fopen ( filenames[ofile].word,"at" ) ;
  if ( f[ofile]==NULL )
    {
      fprintf ( stderr,"\n\n Error during opening file %s", filenames[ofile].word ) ;
      fprintf ( stderr,"\n in function logo::append_load (%s, line %d).\n",__FILE__,__LINE__);
      exit ( 1 ) ;
    }
  filenames[ofile].init=true ;
}

void logo::start_time ( void )
{
  start=clock() ;
}

void logo::print_time ( long olevel, char *olabel, long otag, long osize )
{
  if ( names[olevel].binary & printlevel )
    {
      end=clock() ;
      the_time=(double)(end-start)/((double)osize*CLOCKS_PER_SEC) ;
      print ( olevel, olabel ) ;
      print ( olevel, " %f s \n", the_time );
      if ( otag==tag_stop_time )
	start=clock() ;
    }
}

void logo::save_time ( long olevel, long ofile, char *olabel, long otag, long osize )
{
  if ( filenames[ofile].init && ( names[olevel].binary & savelevel ))
    {
      end=clock() ;
      the_time=(double)(end-start)/((double)osize*CLOCKS_PER_SEC) ;
      save ( olevel, ofile, olabel ) ;
      save ( olevel, ofile, " %f s \n", the_time ) ;
      fflush( f[ofile] ) ;
      if ( otag==tag_stop_time )
	start=clock() ;
    }
}

void logo::print_sys_time ( long olevel )
{
  if ( names[olevel].binary & printlevel )
    {
      struct rusage usage ;

      getrusage( RUSAGE_SELF, &usage ) ;
      print( olevel, "CPU usage: User = %ld.%06ld, System = %ld.%06ld\n",
	     usage.ru_utime.tv_sec, usage.ru_utime.tv_usec,
	     usage.ru_stime.tv_sec, usage.ru_stime.tv_usec ) ;
    }
}

void logo::print ( long olevel, char *format, ... )
{
  if ( names[olevel].binary & printlevel )
    {
      va_list argumenty ;
      
      va_start ( argumenty, format ) ;
      vprintf ( format, argumenty ) ;
      va_end ( argumenty ) ;
    }
}

void logo::save ( long olevel, long ofile, char *format, ... )
{
  if ( filenames[ofile].init && ( names[olevel].binary & savelevel ))
    {
      va_list argumenty ;
      
      va_start ( argumenty, format ) ;
      vfprintf ( f[ofile], format, argumenty ) ;
      va_end ( argumenty ) ;
    }
}

void logo::print_help ( void )
{
  long i ;
  printf( "\nLevels:\n" ) ;
  for ( i=0 ; i<levels ; i++ )
    {
      printf ( "      %32s %10ld \n", names[i].word, names[i].binary ) ;
    }
  printf( "\nFiles:\n" ) ;
  for ( i=0 ; i<files ; i++ )
    {
      printf ( "      %32s \n", filenames[i].word ) ;
    }
}

//
// commands for command line
// by O. Hrstka, sometimes in 2001
//

long scan_cmd_line_direct ( long oargc , char *oargv[] , char omark , long &ovalue )
{
    long i ;
    for ( i=1 ; i<oargc ; i++ )
    {
        if (( oargv[i][0]=='-' ) && ( oargv[i][1]==omark ))
        {
            sscanf( oargv[i]+2,"%ld",&ovalue ) ;
            return 1 ;
        }
    }
    return 0 ;
}

long scan_cmd_line_direct ( long oargc , char *oargv[] , char omark , double &ovalue )
{
    long i ;
    for ( i=1 ; i<oargc ; i++ )
    {
        if (( oargv[i][0]=='-' ) && ( oargv[i][1]==omark ))
        {
            sscanf( oargv[i]+2,"%lf",&ovalue ) ;
            return 1 ;
        }
    }
    return 0 ;
}

long  scan_cmd_line_direct ( long oargc , char *oargv[] , char omark , char *ovalue )
{
    long i ;
    for ( i=1 ; i<oargc ; i++ )
    {
        if (( oargv[i][0]=='-' ) && ( oargv[i][1]==omark ))
        {
            sscanf( oargv[i]+2,"%s",ovalue ) ;
            return 1 ;
        }
    }
    return 0 ;
}

long scan_cmd_line ( long oargc , char *oargv[] , char omark , long &ovalue )
{
    long i ;
    for ( i=1 ; i<oargc ; i++ )
    {
        if (( oargv[i][0]=='-' ) && ( oargv[i][1]==omark ))
        {
            sscanf( oargv[i+1],"%ld",&ovalue ) ;
            return 1 ;
        }
    }
    return 0 ;
}

long scan_cmd_line ( long oargc , char *oargv[] , char omark , double &ovalue )
{
    long i ;
    for ( i=1 ; i<oargc ; i++ )
    {
        if (( oargv[i][0]=='-' ) && ( oargv[i][1]==omark ))
        {
            sscanf( oargv[i+1],"%lf",&ovalue ) ;
            return 1 ;
        }
    }
    return 0 ;
}

long scan_cmd_line ( long oargc , char *oargv[] , char omark , char *ovalue )
{
    long i ;
    for ( i=1 ; i<oargc ; i++ )
    {
        if (( oargv[i][0]=='-' ) && ( oargv[i][1]==omark ))
        {
            sscanf( oargv[i+1],"%s",ovalue ) ;
            return 1 ;
        }
    }
    return 0 ;
}
