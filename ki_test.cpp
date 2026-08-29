/* ki_test.cpp: check keplerian_link.cpp against the worked example in
   G. F. Gronchi,  'Orbit determination with the Keplerian integrals',
   arXiv:2111.02406,  section 'Numerical test with Link2'.

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

   The paper's two tracklets of (4542) Mossotti are 921 days apart,  which
   no method based on a Taylor expansion of the equations of motion could
   link.  Everything here is done in equatorial J2000,  the frame the
   paper's attributables are given in;  the algebra is frame-agnostic,  so
   within find_orb the same code runs in ecliptic coordinates instead.  */

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

/* Heliocentric equatorial J2000 state of Pan-STARRS 1 (MPC code F51) at
   the two mean epochs,  from the JPL ephemeris via astropy.  Baked in so
   that this test needs no ephemeris of its own.       */

static const double observer_state[2][6] = {
   { -7.961908614146815e-01, -5.653682744617511e-01, -2.450680986773588e-01,
     +1.048528815311778e-02, -1.263286487027348e-02, -5.440593985799470e-03 },
   { +7.370634796207131e-01, +6.088201997318874e-01, +2.639333333361885e-01,
     -1.199256527771150e-02, +1.183474375826706e-02, +5.061521869374504e-03 } };

/* The attributables as printed in the paper: alpha, delta (radians),
   alphadot, deltadot (radians/day),  and the mean epoch (MJD).     */

static const double paper_attributable[2][4] = {
   { 4.127242, -0.094234, -0.00316982,  0.00064761 },
   { 0.896144,  0.078622, -0.00364403, -0.00065882 } };

static const double paper_epoch[2] = { 55679.52985, 56600.45442 };

static void set_attributable( ATTRIBUTABLE *attr, const int idx)
{
   const double alpha = paper_attributable[idx][0];
   const double delta = paper_attributable[idx][1];
   const double adot  = paper_attributable[idx][2];
   const double ddot  = paper_attributable[idx][3];
   const double ca = cos( alpha), sa = sin( alpha);
   const double cd = cos( delta), sd = sin( delta);
   const double e_alpha[3] = { -sa, ca, 0. };
   const double e_delta[3] = { -sd * ca, -sd * sa, cd };
   int i;

   attr->t = paper_epoch[idx] + 2400000.5;
   attr->u[0] = cd * ca;
   attr->u[1] = cd * sa;
   attr->u[2] = sd;
   for( i = 0; i < 3; i++)
      {
      attr->eta[i] = adot * cd * e_alpha[i] + ddot * e_delta[i];
      attr->q[i] = observer_state[idx][i];
      attr->qdot[i] = observer_state[idx][i + 3];
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

int main( const int argc, const char **argv)
{
   ATTRIBUTABLE a1, a2;
   LINK2_DATA ld;
   LINK2_ROOT roots[20];
   int i, n_roots, n_bound = 0, n_failures = 0;
   const bool verbose = (argc > 1 && !strcmp( argv[1], "-v"));

   set_attributable( &a1, 0);
   set_attributable( &a2, 1);
   if( link2_setup( &ld, &a1, &a2))
      {
      fprintf( stderr, "link2_setup() failed: D1 parallel to D2\n");
      return( -1);
      }
   if( verbose)
      printf( "conic:  q20=%+.6e q10=%+.6e q02=%+.6e q01=%+.6e q00=%+.6e\n",
                     ld.q20, ld.q10, ld.q02, ld.q01, ld.q00);

   n_roots = link2( roots, 20, &a1, &a2, .05, 60.);
   printf( "%d root(s) of the Link2 system;  bound solutions:\n", n_roots);
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
                     /* Compare to Table 2 of the paper.  The tolerances
                        allow for the attributables having been printed
                        to only six or seven digits.        */
         if( fabs( roots[i].rho1 - 1.8802) > .002) n_failures++;
         if( fabs( roots[i].rho2 - 2.1774) > .002) n_failures++;
         if( fabs( a    -   3.03055) > .005)       n_failures++;
         if( fabs( ecc  -   0.06436) > .001)       n_failures++;
         if( fabs( incl -  11.22246) > .01)        n_failures++;
         if( fabs( node - 104.80204) > .01)        n_failures++;
                     /* Angular momentum must be conserved to rounding.
                        This is what catches the rhodot1 / rhodot2
                        assignment being interchanged.      */
         if( roots[i].max_c_err > 1e-12)           n_failures++;
         }
      }
   if( n_bound != 1)
      {
      printf( "FAILED: expected exactly one bound solution,  got %d\n", n_bound);
      n_failures++;
      }
   printf( "paper:    rho1=   1.8802 rho2=   2.1774  "
           "                a=   3.0306 e= 0.0644 I= 11.2225 Om= 104.8020\n");
   printf( n_failures ? "\n%d CHECK(S) FAILED\n" : "\nAll checks passed\n",
                                                            n_failures);
   return( n_failures);
}
