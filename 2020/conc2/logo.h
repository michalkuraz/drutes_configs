//
//  File:         logo.h
//
//  Description:  Class for logging output
//
//  Author:       Matej Leps
//

#ifndef __logo_h__
#define __logo_h__

#include <stdio.h>
#include <stdarg.h>
#include <time.h>

#define tag_stop_time         0
#define tag_cont_time         1

#define mode_always_ask      0
#define mode_append           1
#define mode_rewrite          2

typedef FILE* PFILE ;
struct logo_keyword ;

class logo
{
  public:
	  logo ( void ) ;
	 ~logo ( void ) ;
	 void set_levels ( long olevels ) ;
	 void set_level_name ( long olevel, char *olevel_name ) ;
	 void set_files ( long ofiles ) ;
	 void set_file_name ( long ofile, char *ofile_name, long omode=mode_always_ask ) ;
	 void load_files ( void ) ;
	 void safe_load ( long ofile ) ;
	 void append_load ( long ofile ) ;
	 void rewrite_load ( long ofile ) ;
	 void start_time ( void ) ;
	 void print_time ( long olevel, char *olabel, long otag=tag_stop_time, long osize=1 ) ;
	 void save_time ( long olevel, long ofile, char *olabel, long otag=tag_stop_time, long osize=1 ) ;
	 void print_sys_time ( long olevel ) ;
	 void print ( long olevel, char *format, ... ) ;
	 void save ( long olevel, long ofile, char *format, ... ) ;
	 void print_help ( void ) ;
  public:
	 FILE **f ;
	 clock_t start, end ;
	 double the_time ;
	 long printlevel ;
	 long savelevel ;
	 long files ;
	 long levels ;
	 logo_keyword *names ;
	 logo_keyword *filenames ;
  public:
	 //
         //        1 ->      1
         //        2 ->      2
         //        4 ->      3
         //        8 ->      4
         //       16 ->      5
	 //
	 inline long give_long ( long obinary )
	   {
	     long num=obinary ;
	     long poc=1 ;
	     while ( !(num & 1) )
	       {
		 num >>= 1;
		 poc++ ;
	       }
	     return poc ;
	   } ;
	 //
         //        1 ->      1
         //        2 ->      2
         //        3 ->      4
         //        4 ->      8
         //        5 ->     16
	 //
	 inline long give_binary ( long olong )
	   {
	     long binary=1 ;
	     binary <<= (olong-1) ;
	     return binary ;
	   } ;
} ;



// for command line
long scan_cmd_line_direct ( long oargc , char *oargv[] , char omark , long &ovalue ) ;
long scan_cmd_line_direct ( long oargc , char *oargv[] , char omark , double &ovalue ) ;
long scan_cmd_line_direct ( long oargc , char *oargv[] , char omark , char *ovalue ) ;
long scan_cmd_line ( long oargc , char *oargv[] , char omark , long &ovalue ) ;
long scan_cmd_line ( long oargc , char *oargv[] , char omark , double &ovalue ) ;
long scan_cmd_line ( long oargc , char *oargv[] , char omark , char *ovalue ) ;

#endif // __logo_h__
