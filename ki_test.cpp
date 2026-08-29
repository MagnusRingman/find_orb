/* ki_test.cpp: check keplerian_link.cpp against the worked examples in
   G. F. Gronchi,  'Orbit determination with the Keplerian integrals',
   arXiv:2111.02406,  sections 'Numerical test with Link2' and
   'Numerical test with Link3'.

Copyright (C) 2026, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.

Compile with

g++ -Wall -O3 -Wextra -pedantic -I$(PREFIX)/include -o ki_test ki_test.cpp \
              keplerian_link.cpp -L$(PREFIX)/lib -llunar

   The Link2 tracklets of (4542) Mossotti are 921 days apart and the Link3
   tracklets of (4628) Laplace span 564 days;  no method based on a Taylor
   expansion of the equations of motion could link either.  Everything here
   is done in equatorial J2000,  the frame the paper's attributables are
   given in;  the algebra is frame-agnostic,  so within find_orb the same
   code runs in ecliptic coordinates instead.

   The mean epochs tabulated in the paper are TT,  not UTC.  That is not
   stated there,  but it is what the numbers require:  reading them as UTC
   shifts the observer by 66 seconds of Earth motion and moves the Link3
   ranges by about 0.004 AU,  roughly twenty times the scatter the
   published rounding of the attributables can account for.  Link2 is far
   less sensitive,  moving only about 0.0003 AU.        */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "watdefs.h"
#include "afuncs.h"
#include "mpc_obs.h"
#include "keplerian_link.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)

/* Heliocentric equatorial J2000 states of Pan-STARRS 1 (MPC code F51) at
   the mean epochs,  from the JPL ephemeris via astropy,  with the epochs
   read as TT.  Baked in so that this test needs no ephemeris of its own. */

static const double mossotti_observer[2][6] = {
   { -7.961988931154208e-01, -5.653585969674157e-01, -2.450639310609458e-01,
     +1.048470199264797e-02, -1.263413320302988e-02, -5.440647950162496e-03 },
   { +7.370728045806806e-01, +6.088109966947408e-01, +2.639293975103021e-01,
     -1.199161296069380e-02, +1.183584414111527e-02, +5.061583527440962e-03 } };

static const double laplace_observer[3][6] = {
   { +8.562470590914019e-01, -4.945194825243194e-01, -2.143587085981919e-01,
     +9.094022194166683e-03, +1.344206274853188e-02, +5.764238203370036e-03 },
   { +8.303688315167107e-01, +5.011952363622151e-01, +2.172709745626089e-01,
     -9.975142139072338e-03, +1.319957590037597e-02, +5.689355845663648e-03 },
   { -9.655655776516187e-01, +2.104402198248298e-01, +9.122577474585647e-02,
     -4.511261694088860e-03, -1.545137910687363e-02, -6.684156427062916e-03 } };

/* alpha, delta (radians),  alphadot, deltadot (radians/day),  as printed
   in the paper,  and the mean epochs in MJD (TT).      */

static const double mossotti_attr[2][4] = {
   { 4.127242, -0.094234, -0.00316982,  0.00064761 },
   { 0.896144,  0.078622, -0.00364403, -0.00065882 } };
static const double mossotti_epoch[2] = { 55679.52985, 56600.45442 };

static const double laplace_attr[3][4] = {
   { 5.497266, -0.067965, -0.00379969, -0.00072536 },
   { 0.715891,  0.542071, -0.00422693, -0.00136864 },
   { 0.831367,  0.390747,  0.00622482,  0.00054073 } };
static const double laplace_epoch[3] = { 55794.36935, 56226.53746, 56358.24760 };

static void set_attributable( ATTRIBUTABLE *attr, const double *elems,
                            const double epoch, const double *observer)
{
   const double alpha = elems[0], delta = elems[1];
   const double adot = elems[2], ddot = elems[3];
   const double ca = cos( alpha), sa = sin( alpha);
   const double cd = cos( delta), sd = sin( delta);
   const double e_alpha[3] = { -sa, ca, 0. };
   const double e_delta[3] = { -sd * ca, -sd * sa, cd };
   int i;

   attr->t = epoch + 2400000.5;
   attr->u[0] = cd * ca;
   attr->u[1] = cd * sa;
   attr->u[2] = sd;
   for( i = 0; i < 3; i++)
      {
      attr->eta[i] = adot * cd * e_alpha[i] + ddot * e_delta[i];
      attr->q[i] = observer[i];
      attr->qdot[i] = observer[i + 3];
      }
}

/* Semimajor axis,  eccentricity,  and the inclination and node measured
   in the ecliptic (which is how the paper reports them).           */

static void elements_from_state( const double *state, double *a, double *ecc,
                                       double *incl, double *node)
{
   const double r = vector3_length( state);
   const double v2 = dot_product( state + 3, state + 3);
   double h[3], evec[3], hlen;
   int i;

   *a = 1. / (2. / r - v2 / SOLAR_GM);
   vector_cross_product( h, state, state + 3);
   vector_cross_product( evec, state + 3, h);
   for( i = 0; i < 3; i++)
      evec[i] = evec[i] / SOLAR_GM - state[i] / r;
   *ecc = vector3_length( evec);
   equatorial_to_ecliptic( h);
   hlen = vector3_length( h);
   *incl = acos( h[2] / hlen) * 180. / PI;
   *node = atan2( h[0], -h[1]) * 180. / PI;
   if( *node < 0.)
      *node += 360.;
}

static int check( const char *what, const double found, const double expected,
                                                    const double tolerance)
{
   const int failed = (fabs( found - expected) > tolerance);

   if( failed)
      printf( "  FAILED %s: got %.5f,  expected %.5f +/- %.5f\n",
                              what, found, expected, tolerance);
   return( failed);
}

static int test_link2( const bool verbose)
{
   ATTRIBUTABLE a1, a2;
   LINK2_ROOT roots[20];
   int i, n_roots, n_bound = 0, rval = 0;

   printf( "Link2,  (4542) Mossotti,  two tracklets 921 days apart:\n");
   set_attributable( &a1, mossotti_attr[0], mossotti_epoch[0],
                                            mossotti_observer[0]);
   set_attributable( &a2, mossotti_attr[1], mossotti_epoch[1],
                                            mossotti_observer[1]);
   n_roots = link2( roots, 20, &a1, &a2, .05, 60.);
   if( n_roots < 0)
      {
      printf( "  FAILED: link2() could not set up\n");
      return( 1);
      }
   for( i = 0; i < n_roots; i++)
      {
      double a, ecc, incl, node;

      elements_from_state( roots[i].state1, &a, &ecc, &incl, &node);
      if( verbose || a > 0.)
         printf( "  rho1=%9.4f rho2=%9.4f  |c1-c2|=%.2e  "
                 "a=%9.4f e=%7.4f I=%8.4f Om=%9.4f\n",
                 roots[i].rho1, roots[i].rho2, roots[i].max_c_err,
                 a, ecc, incl, node);
      if( a > 0.)
         {
         n_bound++;
         rval += check( "rho1", roots[i].rho1,   1.8802, .002);
         rval += check( "rho2", roots[i].rho2,   2.1774, .002);
         rval += check( "a",    a,               3.03055, .005);
         rval += check( "e",    ecc,             0.06436, .001);
         rval += check( "I",    incl,           11.22246, .01);
         rval += check( "Om",   node,          104.80204, .01);
                     /* Angular momentum must be conserved to rounding.
                        This is what catches the rhodot1 / rhodot2
                        assignment being interchanged.      */
         rval += check( "|c1-c2|", roots[i].max_c_err, 0., 1e-12);
         }
      }
   printf( "  paper:   1.8802     2.1774                  "
           "  3.0306   0.0644  11.2225  104.8020\n");
   if( n_bound != 1)
      {
      printf( "  FAILED: expected exactly one bound solution,  got %d\n",
                                                            n_bound);
      rval++;
      }
   return( rval);
}

/* The paper's two surviving Link3 triplets.  The first has badly
   inconsistent elements from epoch to epoch and would be thrown out by the
   compatibility test;  the second is the real orbit.  Both must appear
   among our roots.                                        */

static const double laplace_expected[2][3] = {
   { 2.1955, 1.9028, 2.9200 },
   { 1.9379, 1.8279, 2.8870 } };

static int test_link3( const bool verbose)
{
   ATTRIBUTABLE attr[3];
   LINK3_ROOT roots[20];
   int i, k, n_roots, n_accepted = 0, rval = 0;
   int found[2] = { 0, 0 };
   double sl[3];

   printf( "\nLink3,  (4628) Laplace,  three tracklets spanning 564 days:\n");
   for( k = 0; k < 3; k++)
      set_attributable( attr + k, laplace_attr[k], laplace_epoch[k],
                                                   laplace_observer[k]);
   for( k = 0; k < 3; k++)
      sl[k] = straight_line_rho( attr + k);
   if( verbose)
      printf( "  straight-line rho (closed form):  %.4f %.4f %.4f\n",
                                             sl[0], sl[1], sl[2]);
   n_roots = link3( roots, 20, attr, attr + 1, attr + 2, .05, 60.);
   if( n_roots < 0)
      {
      printf( "  FAILED: link3() could not set up\n");
      return( 1);
      }
   for( i = 0; i < n_roots; i++)
      {
      const LINK3_ROOT *root = roots + i;
      double a[3], ecc[3], incl, node;
      bool bound = true, straight_line = (root->c_len < 1e-8);
      bool positive = (root->rho[0] > 0. && root->rho[1] > 0.
                                         && root->rho[2] > 0.);

      for( k = 0; k < 3; k++)
         {
         elements_from_state( root->state[k], a + k, ecc + k, &incl, &node);
         if( a[k] <= 0.)
            bound = false;
         }
      printf( "  rho=(%8.4f,%8.4f,%8.4f) |c|=%.2e |dc|=%.1e  %s\n",
              root->rho[0], root->rho[1], root->rho[2],
              root->c_len, root->max_c_err,
              straight_line ? "straight line" :
                    (!positive ? "rho <= 0" : (!bound ? "unbound" : "ACCEPTED")));
      rval += check( "|c_i-c_j|", root->max_c_err, 0., 1e-12);
      if( straight_line)     /* must agree with the closed form */
         for( k = 0; k < 3; k++)
            rval += check( "straight-line rho", root->rho[k], sl[k], 1e-4);
      if( bound && positive && !straight_line)
         {
         n_accepted++;
         printf( "        a = %.5f / %.5f / %.5f   e = %.5f / %.5f / %.5f\n",
                          a[0], a[1], a[2], ecc[0], ecc[1], ecc[2]);
         for( k = 0; k < 2; k++)
            if( fabs( root->rho[0] - laplace_expected[k][0]) < .005
             && fabs( root->rho[1] - laplace_expected[k][1]) < .005
             && fabs( root->rho[2] - laplace_expected[k][2]) < .005)
               found[k] = 1;
         }
      }
   printf( "  paper:  (  2.1955,  1.9028,  2.9200)"
           "  a = 2.86808 / 2.64520 / 2.59619\n");
   printf( "          (  1.9379,  1.8279,  2.8870)"
           "  a = 2.64614 / 2.64562 / 2.64427\n");
   if( n_accepted != 2)
      {
      printf( "  FAILED: expected two accepted solutions,  got %d\n",
                                                       n_accepted);
      rval++;
      }
   for( k = 0; k < 2; k++)
      if( !found[k])
         {
         printf( "  FAILED: did not recover the paper's solution %d\n", k + 1);
         rval++;
         }
   return( rval);
}

int main( const int argc, const char **argv)
{
   const bool verbose = (argc > 1 && !strcmp( argv[1], "-v"));
   int n_failures = test_link2( verbose) + test_link3( verbose);

   printf( n_failures ? "\n%d CHECK(S) FAILED\n" : "\nAll checks passed\n",
                                                            n_failures);
   return( n_failures);
}
