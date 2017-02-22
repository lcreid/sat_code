/*
    sat_id.cpp     8 March 2003,  with updates as listed below

   An example 'main' function illustrating how to find which satellite(s)
are within a given radius of a given RA/dec,  as seen from a given
point.  The code reads in a file of observations in MPC format (name
provided as the first command-line argument).  For example:

sat_id mpc_obs.txt

   would hunt through the file 'mpc_obs.txt' for MPC-formatted
observations.  It would then read the file 'alldat.tle',  looking
for corresponding satellites within .2 degrees of said observations.
It then spits out the original file,  with satellite IDs added
(when found) after each observation line.  For each IDed satellite,
the international and NORAD designations are given,  along with
its angular distance from the search point,  position angle of
motion,  and apparent angular rate of motion in arcminutes/second
(or,  equivalently,  degrees/minute). */

/* 2 July 2003:  fixed the day/month/year to JD part of 'get_mpc_data()'
so it will work for all years >= 0 (previously,  it worked for years
2000 to 2099... plenty for the practical purpose of ID'ing recently-found
satellites,  but this is also an 'example' program.) */

/* 3 July 2005:  revised the check on the return value for parse_elements().
Now elements with bad checksums won't be rejected. */

/* 23 June 2006:  after comment from Eric Christensen,  revised to use
names 'ObsCodes.html' or 'ObsCodes.htm',  with 'stations.txt' being a
third choice.  Also added the '-a' command line switch to cause the program
to show all lines from input (default is now that only MPC astrometric
input gets echoed.)   */

/* 30 June 2006:  further comment from Eric Christensen:  when computing
object motion from two consecutive observations,  if the second one has
a date/time preceding the first,  you get a negative rate of motion that's
off by 180 degrees.  Fixed this. */

/* 17 Nov 2006:  artificial satellite data is now being provided in a
file named 'ALL_TLE.TXT'.  I've modified the default TLE to match. */

/* 22 Oct 2012:  minor cosmetic changes,  such as making constant variables
of type 'const',  updating URL for the MPC station code file,  adding a
comment or two. */

/* 7 Jan 2013:  revised output to show satellite name if available,  plus
the eccentricity,  orbital period,  and inclination. */

/* 2013 Dec 8:  revised to pay attention to "# MJD" and "#Ephem start"
lines,  for files that contain many TLEs covering different time spans
for the same object.  I sometimes create such files;  when that happens,
for each observation,  only the TLE(s) covering that observation's time
should be used,  and the others are suppressed.       */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include "norad.h"
#include "observe.h"

#define MAX_N_MATCHES 50
#define OBSERVATION struct observation

OBSERVATION
   {
   char text[81];
   double jd, ra, dec;
   double lon, rho_cos_phi, rho_sin_phi;
   char matches[MAX_N_MATCHES][170];
   };

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define TIME_EPSILON (1./86400.)

static char *fgets_trimmed( char *buff, const int buffsize, FILE *ifile)
{
   char *rval = fgets( buff, buffsize, ifile);

   if( rval)
      {
      size_t i = 0;

      while( rval[i] != 10 && rval[i] != 13 && rval[i])
         i++;
      rval[i] = '\0';
      }
   return( rval);
}

static int get_mpc_data( OBSERVATION *obs, const char *buff)
{
   int i1, i2, i, year, month;
   double tval, day;
   static const char month_len[12] = { 31, 28, 31, 30, 31, 30,
                                       31, 31, 30, 31, 30, 31 };

   if( strlen( buff) != 80)
      return( -1);
   if( sscanf( buff + 32, "%d %d %lf", &i1, &i2, &tval) != 3)
      return( -2);
   obs->ra = ((double)i1 + (double)i2 / 60. + tval / 3600.) * (PI / 12.);

   if( sscanf( buff + 45, "%d %d %lf", &i1, &i2, &tval) != 3)
      return( -3);
   obs->dec = ((double)i1 + (double)i2 / 60. + tval / 3600.) * (PI / 180.);
   if( buff[44] == '-')
      obs->dec = -obs->dec;
   else if( buff[44] != '+')
      return( -4);

               /* Read in the day/month/year from the record... */
   if( sscanf( buff + 15, "%d %d %lf", &year, &month, &day) != 3)
      return( -5);
   if( month < 1 || month > 12 || day < 1.
                     || day > 1. + (double)month_len[month - 1])
      return( -6);
               /* ...and convert to a JD (NOTE:  non-negative years only) */
   obs->jd = 1721059.5
         + (double)( year * 365 + year / 4 - year / 100 + year / 400) + day;
   for( i = 0; i < month - 1; i++)
      obs->jd += (double)month_len[i];
   if( month < 3 && !(year % 4))    /* leap years,  January and February */
      if( !(year % 400) || (year % 100))
         (obs->jd)--;
   strcpy( obs->text, buff);
   return( 0);
}

/* This loads up the file 'ObsCodes.html' into memory on its first call.
Then,  given an MPC code,  it finds the corresponding line and copies
it into 'station_code_data'. */

int verbose = 0;

static int get_station_code_data( char *station_code_data,
                  const char *mpc_code)
{
   static char *cached_data, *cached_ptr;

   if( !mpc_code)       /* freeing memory */
      {
      if( cached_data)
         free( cached_data);
      cached_data = cached_ptr = NULL;
      return( 0);
      }
   *station_code_data = '\0';
   if( !cached_data)
      {
      FILE *ifile = fopen( "ObsCodes.html", "rb");
      size_t size;

      if( !ifile)              /* perhaps stored with truncated extension? */
         ifile = fopen( "ObsCodes.htm", "rb");
      if( !ifile)
         {
         printf( "Failed to find MPC station list 'ObsCodes.html'\n");
         printf( "This can be downloaded at:\n\n");
         printf( "http://www.minorplanetcenter.org/iau/lists/ObsCodes.html\n");
         exit( -3);
         }
      fseek( ifile, 0L, SEEK_END);
      size = (size_t)ftell( ifile);
      fseek( ifile, 0L, SEEK_SET);
      cached_data = (char *)malloc( size + 1);
      if( fread( cached_data, 1, size, ifile) != size)
         {
         printf( "Failed to read station file\n");
         exit( -4);
         }
      fclose( ifile);
      cached_data[size] = '\0';
      if( verbose)
         printf( "Station codes: %u bytes read\n", (unsigned)size);
      }
   if( !cached_ptr || memcmp( cached_ptr, mpc_code, 3))
      {
      char search_buff[5];

      sprintf( search_buff, "\n%.3s", mpc_code);
      cached_ptr = strstr( cached_data, search_buff);
      if( cached_ptr)
         cached_ptr++;
      }
   if( cached_ptr)
      {
      size_t i;

      for( i = 0; cached_ptr[i] >= ' '; i++)
         station_code_data[i] = cached_ptr[i];
      station_code_data[i] = '\0';
      }
   else
      {
      printf( "Station code '%s' not found.\n", mpc_code);
#ifdef ON_LINE_VERSION
      printf( "If this is a new MPC code,  it could be that this service needs to be\n");
      printf( "updated to know about it.  Please contact pluto at projectpluto.com so\n");
      printf( "I can fix this.\n");
#else
      printf( "If this is a new MPC code,  you may need to get this file:\n");
      printf( "http://www.minorplanetcenter.org/iau/lists/ObsCodes.html\n");
      printf( "and replace the existing ObsCodes.html.\n");
#endif
      }
   return( cached_ptr ? 0 : -1);
}

/* Loads up MPC-formatted 80-column observations from a file.  Makes
a pass to find out how many observations there are,  allocates space
for them,  then reads them again to actually load the observations. */

static OBSERVATION *get_observations_from_file( FILE *ifile, size_t *n_found)
{
   int pass;
   OBSERVATION *rval = NULL, obs;

   memset( &obs, 0, sizeof( OBSERVATION));
   for( pass = 0; pass < 2; pass++)
      {
      char buff[100];
      size_t count = 0;

      fseek( ifile, 0L, SEEK_SET);
      while( fgets_trimmed( buff, sizeof( buff), ifile))
         if( !get_mpc_data( &obs, buff))
            {
            if( rval)
               {
               char station_data[100];

               if( get_station_code_data( station_data, obs.text + 77))
                  printf( "FAILED to find MPC code %s\n", obs.text + 77);
               sscanf( station_data + 3, "%lf %lf %lf", &obs.lon,
                                        &obs.rho_cos_phi, &obs.rho_sin_phi);
               obs.lon *= PI / 180.;
//             if( fabs( obs.rho_cos_phi) > 1.1 || fabs( obs.rho_sin_phi) > 1.1)
//                {                                            /* lat and alt actually given; */
//                const double lat = obs.rho_cos_phi * PI / 180.;  /* cvt them to parallax data   */
//                const double alt_in_meters = obs.rho_sin_phi;
//
//                lat_alt_to_parallax( lat, alt_in_meters,
//                            &obs.rho_cos_phi, &obs.rho_sin_phi);
//                printf( "Parallax constants for %s: %.6f %.6f\n", obs.text + 77,
//                            obs.rho_cos_phi, obs.rho_sin_phi);
//                }
               rval[count] = obs;
               }
            count++;
            }
      if( !pass)
         rval = (OBSERVATION *)calloc( count, sizeof( OBSERVATION));
      *n_found = count;
      }
   return( rval);
}

static int id_compare( const OBSERVATION *a, const OBSERVATION *b)
{
   return( memcmp( a->text, b->text, 12));
}

static int compare_obs( const void *a, const void *b, void *context)
{
   const OBSERVATION *aptr = (const OBSERVATION *)a;
   const OBSERVATION *bptr = (const OBSERVATION *)b;
   int rval = id_compare( aptr, bptr);

   if( !rval)        /* same IDs?  Then sort by JD of observation */
      rval = (aptr->jd > bptr->jd ? 1 : -1);
   return( rval);
}

/* Copied straight from 'mpc_obs.cpp' in Find_Orb.  See comments there. */

void shellsort_r( void *base, const size_t n_elements, const size_t esize,
         int (*compare)(const void *, const void *, void *), void *context)
{
#if (defined _GNU_SOURCE || defined __GNU__ || defined __linux)
   qsort_r( base, n_elements, esize, compare, context);
#else
   size_t gap = 1;
   char *data = (char *)base;
   char *pivot = (char *)alloca( esize);

   while( gap < n_elements)
      gap = gap * 3 + 1;
   while( gap /= 3)
      {
      size_t j;
      const size_t spacing = esize * gap;

      for( j = gap; j < n_elements; j++)
         {
         char *tptr = data + j * esize;
         char *tptr2 = tptr - spacing;

         if( (compare)( tptr2, tptr, context) > 0)
            {
            memcpy( pivot, tptr, esize);
            memcpy( tptr, tptr2, esize);
            tptr = tptr2;
            tptr2 -= spacing;
            while( tptr2 >= base && (compare)( tptr2, pivot, context) > 0)
               {
               memcpy( tptr, tptr2, esize);
               tptr = tptr2;
               tptr2 -= spacing;
               }
            memcpy( tptr, pivot, esize);
            }
         }
      }
#endif
}

/* Quick and dirty computation of apparent motion,  good enough
   for our humble purposes.  If you want absolute perfection,  see
   the 'dist_pa.cpp' code at http://www.projectpluto.com/source.htm . */

static int compute_motion( const double delta_t,
                const double d_ra,  const double d_dec,
                double *arcmin_per_sec, double *posn_ang)
{
   int rval = 0;

   if( delta_t && (d_ra || d_dec))
      {
      double total_motion = sqrt( d_ra * d_ra + d_dec * d_dec);

      *posn_ang = atan2( d_ra, d_dec) * 180. / PI;
      if( *posn_ang < 0.)
         *posn_ang += 360.;
      if( delta_t < 0.)
         {
         *posn_ang = fmod( *posn_ang + 180., 360.);
         total_motion *= -1.;
         }
      *arcmin_per_sec = (total_motion * 180. / PI) / (delta_t * 1440.);
      }
   else     /* undefined or no motion */
      {
      rval = -1;
      *arcmin_per_sec = *posn_ang = 0.;
      }
   return( rval);
}

static int compute_motion_between_obs( const OBSERVATION *obs1,
             const OBSERVATION *obs2,
             double *arcmin_per_sec, double *posn_ang)
{
   *arcmin_per_sec = *posn_ang = 0.;
   int rval = 0;

   if( !id_compare( obs1, obs2) && fabs( obs1->jd - obs2->jd) < .3)
      {
      const double dt = obs1->jd - obs2->jd;
      const double ra1 = obs1->ra,  dec1 = obs1->dec;
      const double ra2 = obs2->ra,  dec2 = obs2->dec;

      rval = compute_motion( dt,
                  (ra1 - ra2) * cos( (dec1 + dec2) / 2.),
                  (dec1 - dec2), arcmin_per_sec, posn_ang);
      }
   return( rval);
}

/* Out of all observations for a given object,  this function will keep
the last one,  and the preceding observation within 0.3 days.  The idea is
to get a preceding observation that will give a good idea of the object's
motion.  Of course,  in some cases,  we'll only have one observation and
the match will have to be based on position only,  which is usually not
all that reliable. */

static int drop_extra_obs( OBSERVATION *obs, const int n_obs)
{
   int i = 0, j, rval = 0;

   while( i < n_obs)
      {
      int k;

      j = i + 1;
      while( j < n_obs && !id_compare( obs + i, obs + j))
         j++;
      j--;           /* j = idx of last observation of this object */
      k = j;
      while( k > i && obs[j].jd - obs[k - 1].jd < 0.3
                  && !memcmp( &obs[j].text[77], &obs[k].text[77], 3))
         k--;
      if( k != j)
         {
         obs[rval++] = obs[k];
         }
      obs[rval++] = obs[j];
      i = j + 1;
      }
   return( rval);
}

static void error_exit( const int exit_code)
{
   printf(
"sat_id takes the name of an input file of MPC-formatted (80-column)\n\
astrometry as a command-line argument.  It searches for matches between\n\
the observation data and satellites in 'ALL_TLE.TXT'.  By default,  matches\n\
within .2 degrees are shown.\n\n\
Additional command-line arguments are:\n\
   -r(radius)    Reset search distance from the default of .2 degrees.\n\
   -t(filename)  Reset the filename of the .tle file.\n\
   -a            Show all lines from input,  not just those with astrometry.\n");
   exit( exit_code);
}

static bool is_in_range( const double jd, const double tle_start,
                                             const double tle_range)
{
   return( !tle_range || !tle_start ||
            (jd >= tle_start && jd <= tle_start + tle_range));
}

/* Given a set of MPC observations and a TLE file,  this function looks at
each TLE in the file and checks to see if that satellite came close to any
of the observations.  We keep track of the MAX_N_MATCHES closest satellites
to any given observation.  The function is called for each TLE file.
*/

int norad_id = 0;

static int add_tle_to_obs( OBSERVATION *obs, const size_t n_obs,
             const char *tle_file_name, const double search_radius,
             const double max_revs_per_day)
{
   char line0[100], line1[100], line2[100];
   double tle_start = 0., tle_range = 0.;
   FILE *tle_file = fopen( tle_file_name, "rb");
   int rval = 0, n_tles_found = 0;
   bool force_sgp4_only = false;
   bool check_updates = true;

   if( !tle_file)
      {
      printf( "Couldn't open TLE file %s\n", tle_file_name);
      return( -1);
      }
   if( verbose)
      printf( "Looking through TLE file '%s', %u obs, radius %f, max %f revs/day\n",
                 tle_file_name, (unsigned)n_obs, search_radius, max_revs_per_day);
   *line0 = *line1 = '\0';
   while( fgets_trimmed( line2, sizeof( line2), tle_file))
      {
      tle_t tle;  /* Structure for two-line elements set for satellite */
      const double mins_per_day = 24. * 60.;

      if( parse_elements( line1, line2, &tle) >= 0
                 && (tle.ephemeris_type == 'H'
                 || tle.xno < 2. * PI * max_revs_per_day / mins_per_day)
                 && (!norad_id || norad_id == tle.norad_number))
         {                           /* hey! we got a TLE! */
         const int is_deep = (force_sgp4_only ? 0 : select_ephemeris( &tle));
//       const int is_deep = select_ephemeris( &tle);
         double sat_params[N_SAT_PARAMS];

         if( verbose > 1)
            printf( "TLE found:\n%s\n%s\n", line1, line2);
         n_tles_found++;
         if( is_deep)
            SDP4_init( sat_params, &tle);
         else
            SGP4_init( sat_params, &tle);
         for( size_t idx = 0; idx < n_obs; idx++)
            if( is_in_range( obs[idx].jd, tle_start, tle_range))
               {
               OBSERVATION *optr = obs + idx;
               const double target_ra  = optr->ra;
               const double target_dec = optr->dec;
               const double jd         = optr->jd;
               double observer_loc[3], observer_loc2[3];
               double radius, d_ra, d_dec;
               double ra, dec, dist_to_satellite, t_since;
               double pos[3]; /* Satellite position vector */
               double unused_delta2;
               int sxpx_rval;

               observer_cartesian_coords( jd, optr->lon,
                       optr->rho_cos_phi, optr->rho_sin_phi, observer_loc);
               observer_cartesian_coords( jd + TIME_EPSILON,
                      optr->lon, optr->rho_cos_phi, optr->rho_sin_phi, observer_loc2);

               t_since = (jd - tle.epoch) * 1440.;
               if( is_deep)
                  sxpx_rval = SDP4( t_since, &tle, sat_params, pos, NULL);
               else
                  sxpx_rval = SGP4( t_since, &tle, sat_params, pos, NULL);
               if( verbose > 2 && sxpx_rval)
                  printf( "TLE failed: %d\n", sxpx_rval);
               get_satellite_ra_dec_delta( observer_loc, pos,
                                    &ra, &dec, &dist_to_satellite);
               epoch_of_date_to_j2000( jd, &ra, &dec);
               d_ra = (ra - target_ra + PI * 4.);
               while( d_ra > PI)
                  d_ra -= PI + PI;
               d_dec = dec - target_dec;
               radius = sqrt( d_ra * d_ra + d_dec * d_dec) * 180. / PI;
               if( radius < search_radius)      /* good enough for us! */
                  {
                  double arcmin_per_sec, posn_ang;

                                 /* Compute position one second later,  so we */
                                 /* can show speed/PA of motion: */
                  t_since += TIME_EPSILON * 1440.;
                  if( is_deep)
                     SDP4( t_since, &tle, sat_params, pos, NULL);
                  else
                     SGP4( t_since, &tle, sat_params, pos, NULL);
                  get_satellite_ra_dec_delta( observer_loc2, pos,
                                       &d_ra, &d_dec, &unused_delta2);
                  epoch_of_date_to_j2000( jd, &d_ra, &d_dec);
                  if( verbose > 1)
                     printf( "Made a match to obs %d\n", (int)idx + 1);
                  d_ra -= ra;
                  d_dec -= dec;
                  while( d_ra > PI)
                     d_ra -= PI + PI;
                  while( d_ra < -PI)
                     d_ra += PI + PI;
                          /* Put RA into 0 to 2pi range: */
                  if( !compute_motion( TIME_EPSILON, d_ra * cos( dec), d_dec,
                              &arcmin_per_sec, &posn_ang))
                     {
                     char obuff[200];
                     char full_intl_desig[20];
                     size_t match_loc = 0;

                     line1[8] = line1[16] = '\0';
                     memcpy( line1 + 30, line1 + 11, 6);
                     line1[11] = '\0';
                     sprintf( full_intl_desig, "%s%.2s-%s",
                              (tle.intl_desig[0] < '5' ? "20" : "19"),
                              tle.intl_desig, tle.intl_desig + 2);
                     sprintf( obuff, "      %5dU = %-9s",
                           tle.norad_number, full_intl_desig);
                     sprintf( obuff + strlen( obuff),
                               "e=%.2f; P=%.1f min; i=%.1f",
                               tle.eo, 2. * PI / tle.xno,
                               tle.xincl * 180. / PI);
                     if( strlen( line0) < 30)         /* object name given... */
                        sprintf( obuff + strlen( obuff), ": %s", line0);
                     obuff[79] = '\0';    /* avoid buffer overrun */
                     strcat( obuff, "\n");
                     sprintf( obuff + strlen( obuff),
                        "             motion %6.3f'/sec at PA %.1f; dist=%8.1f km; offset=%5.2f deg\n",
                            arcmin_per_sec, posn_ang,
                            dist_to_satellite, radius);
                              /* "Speed" is displayed in arcminutes/second,
                                  or in degrees/minute */
                     while( match_loc < MAX_N_MATCHES &&
                                    obs[idx].matches[match_loc][0])
                        {
                        double r = atof( strstr( obs[idx].matches[match_loc], "offset=") + 7);

                        if( radius < r)
                           break;
                        match_loc++;
                        }
                     if( match_loc < MAX_N_MATCHES)
                        {
                        for( size_t i = MAX_N_MATCHES - 1; i > match_loc; i--)
                           strcpy( obs[idx].matches[i], obs[idx].matches[i - 1]);
                        assert( strlen( obuff) < 169);
                        strcpy( obs[idx].matches[match_loc], obuff);
                        }
                     }
                  else
                     if( verbose)
                        printf( "Motion check failed on obs %d\n", (int)idx + 1);
                  }
               }
         }
      else if( !memcmp( line2, "# SGP4 only", 11))
         force_sgp4_only = true;
      else if( !memcmp( line2, "# No updates", 12))
         check_updates = false;
      else if( !memcmp( line2, "# Ephem range:", 14))
         {
         const double mjd_1970 = 40587.;     /* MJD for 1970 Jan 1 */
         double mjd_start, mjd_end;
         double curr_mjd = mjd_1970 + (double)time( NULL) / 86400.;

         sscanf( line2 + 14, "%lf %lf %lf\n", &mjd_start, &mjd_end, &tle_range);
         if( check_updates && mjd_end < curr_mjd + 7.)
            printf( "WARNING:  Update TLEs in '%s'\n", tle_file_name);
         }
      else if( !memcmp( line2, "# MJD ", 6))
         tle_start = atof( line2 + 6) + 2400000.5;
      else if( !memcmp( line2, "# Include ", 10))
         rval = add_tle_to_obs( obs, n_obs, line2 + 10, search_radius,
                                    max_revs_per_day);
      strcpy( line0, line1);
      strcpy( line1, line2);
      }
   if( verbose)
      printf( "%d TLEs read from '%s'\n", n_tles_found, tle_file_name);
   fclose( tle_file);
   return( rval);
}


static void extract_motion( const char *match_text, double *xvel, double *yvel)
{
   const char *tptr = strstr( match_text, "'/sec at PA");
   double total_vel, posn_angle;

   assert( tptr);
   total_vel = atof( tptr - 7);
   posn_angle = atof( tptr + 11) * PI / 180.;
   *xvel = total_vel * sin( posn_angle);
   *yvel = total_vel * cos( posn_angle);
}

/* The "on-line version",  sat_id2,  gathers data from a CGI multipart form,
   puts it into a file,  possibly adds in some options,  puts together the
 command-line arguments,  and then calls sat_id_main.            */

#ifdef ON_LINE_VERSION
int sat_id_main( const int argc, const char **argv)
#else
int main( const int argc, const char **argv)
#endif
{
   const char *tle_file_name = "tle_list.txt";
   FILE *ifile = fopen( argv[1], "rb");
   OBSERVATION *obs;
   size_t n_obs;
   double search_radius = 4;     /* default to 4-degree search */
   double max_revs_per_day = 6.;
   double speed_cutoff = 0.;
   double motion_mismatch_limit = .015;  /* up to .015'/s motion discrepancy OK */
   double motion_mismatch_limit_squared;
// int debug_level = 0;
   int rval;

   if( argc == 1)
      error_exit( -2);

   if( !ifile)
      {
      printf( "Couldn't open input file %s\n", argv[1]);
      return( -1);
      }
   obs = get_observations_from_file( ifile, &n_obs);
   fclose( ifile);
   printf( "%d observations found\n", (int)n_obs);
   shellsort_r( obs, n_obs, sizeof( obs[0]), compare_obs, NULL);
   n_obs = drop_extra_obs( obs, n_obs);
   printf( "%d observations left after dropping extras\n", (int)n_obs);
   if( verbose)
      for( int i = 0; i < argc; i++)
         printf( "Arg %d: '%s'\n", i, argv[i]);
   if( !obs || !n_obs)
      return( -2);
   for( int i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'r':
               search_radius = atof( argv[i] + 2);
               break;
            case 'y':
               motion_mismatch_limit = atof( argv[i] + 2);
               break;
            case 'm':
               max_revs_per_day = atof( argv[i] + 2);
               break;
            case 'n':
               norad_id = atoi( argv[i] + 2);
               break;
            case 't':
               tle_file_name = argv[i] + 2;
               break;
               break;
            case 'v':
               verbose = atoi( argv[i] + 2) + 1;
               break;
//          case 'd':
//             debug_level = atoi( argv[i] + 2);
//             break;
            case 'z':
               speed_cutoff = atof( argv[i] + 2);
               break;
            default:
               printf( "Unrecognized command-line option '%s'\n", argv[i]);
               exit( -2);
               break;
            }

   rval = add_tle_to_obs( obs, n_obs, tle_file_name, search_radius,
                                    max_revs_per_day);
   motion_mismatch_limit_squared = motion_mismatch_limit * motion_mismatch_limit;
   for( size_t idx = 0; idx < n_obs; idx++)
      if( idx == n_obs - 1 || id_compare( obs + idx, obs + idx + 1))
         {
         double motion = 0., posn_ang = 0.;

         if( idx)
            compute_motion_between_obs( obs + idx, obs + idx - 1,
                        &motion, &posn_ang);
         if( motion >= speed_cutoff)
            {
            const double xvel0 = motion * sin( posn_ang * PI / 180.);
            const double yvel0 = motion * cos( posn_ang * PI / 180.);
            int n_matches_shown = 0;

            for( size_t i = 0; i < MAX_N_MATCHES && obs[idx].matches[i][0]; i++)
               {
               double xvel, yvel;

               extract_motion( obs[idx].matches[i], &xvel, &yvel);
               xvel -= xvel0;
               yvel -= yvel0;
               if( xvel * xvel + yvel * yvel < motion_mismatch_limit_squared)
                  {
                  if( !n_matches_shown)
                     {
                     printf( "\n%s\n", obs[idx].text);
                     printf( "    Object motion is %.3f'/sec at PA %.1f\n",
                        motion, posn_ang);
                     }
                  printf( "%s", obs[idx].matches[i]);
                  n_matches_shown++;
                  }
               }
            }
         }
   free( obs);
   get_station_code_data( NULL, NULL);
   return( rval);
}     /* End of main() */

