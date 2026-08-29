/* keplerian_link.cpp: linkage of too-short arcs via the Keplerian integrals

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
02110-1301, USA.    */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "watdefs.h"
#include "afuncs.h"
#include "comets.h"
#include "mpc_obs.h"
#include "keplerian_link.h"

/* Initial orbit determination for the 'linkage' problem:  given two
   tracklets ("too-short arcs",  TSAs) that may or may not belong to the
   same object,  find the orbits consistent with both.  This implements
   the 'Link2' method of

   G. F. Gronchi,  'Orbit determination with the Keplerian integrals',
   arXiv:2111.02406 (2021),  reviewing Gronchi,  Bau & Maro (2015).

   Unlike Gauss/Laplace/Herget/Vaisala,  which expand the _equations of
   motion_ about a central epoch and therefore need observations close
   together in time,  this uses the _conservation laws_ of the two-body
   problem.  There is accordingly no constraint on the time between the
   two tracklets;  they can be separated by years,  or by more than one
   orbital period.  That is the entire point of the method.

   An 'attributable' compresses a tracklet to the line of sight and its
   time derivative at the mean epoch.  Gronchi writes this as
   (alpha, delta, alphadot, deltadot) in spherical coordinates,  with a
   frame (e^rho, e^alpha, e^delta).  We store instead the unit vector
   u = e^rho and its derivative eta = du/dt,  which is what the algebra
   actually uses:

      eta = alphadot cos(delta) e^alpha + deltadot e^delta

   That is the same quantity,  but avoids the cos(delta) factors and the
   coordinate singularity at the poles,  and matches how find_orb already
   stores a line of sight (OBSERVE.vect,  a unit vector in ecliptic
   coordinates).  Everything below is frame-agnostic:  it only requires
   that u,  eta,  q and qdot are expressed in one consistent frame.

   With  r = q + rho u  and  rdot = qdot + rhodot u + rho eta,  the
   angular momentum r x rdot expands to

      c(rho, rhodot) = D rhodot + E rho^2 + F rho + G

   with the four vectors (Gronchi's eq. 'DEFG',  rewritten frame-free)

      D = q x u
      E = u x eta
      F = q x eta + u x qdot
      G = q x qdot

   Note that c is linear in rhodot and quadratic in rho,  which is what
   makes the elimination possible.                                    */

#define J2000 2451545.
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)
      /* AU_PER_DAY -- the speed of light in AU/day -- comes from afuncs.h */

/* Step used when differencing the map from attributables to Delta,  as a
   fraction of the standard deviation of the component being varied.  A
   full sigma of the rate is a large perturbation -- for a half-hour
   tracklet it is most of a percent of the proper motion -- and Link3 is
   stiff enough that the root then moves onto another branch:  on the
   paper's (4628) Laplace case,  steps of a fifth of a sigma or more do
   not reconverge at the first attempt,  and rely on the step-shrinking
   fallback below.

   Measured on that case,  the chi^2 of the _good_ solution is 0.0155,
   holding to within a few tenths of a percent for steps between a
   two-hundredth and a thirtieth of a sigma,  and to about three percent
   out to a tenth.  The chi^2 of the bad solution is not stable at all:  it
   wanders over 179 to 1676 across the same range of steps,  because with
   Delta some eighty degrees from zero the map is nowhere near linear over
   a perturbation of this size.  That does not matter for the purpose --
   the two stay at least four orders of magnitude apart however the step is
   chosen -- but it does mean the absolute chi^2 of a clearly rejected
   solution must not be read as a probability.  Only the chi^2 of a
   solution that nearly passes is quantitatively meaningful.        */

#ifndef FD_STEP_IN_SIGMAS
   #define FD_STEP_IN_SIGMAS .03
#endif

/* Least-squares fit of a straight line to each component of the observed
   unit vectors,  giving the line of sight and its rate at the mean epoch.
   For a real tracklet the arc is short enough that a linear fit is ample;
   Gronchi allows a quadratic fit,  which would matter only for a fast
   mover.  We renormalise u and force eta perpendicular to it,  since
   |u| = 1 implies u . du/dt = 0.                                     */

int compute_attributable( ATTRIBUTABLE *attr, const OBSERVE FAR *obs,
                                                        const int n_obs)
{
   const double arcsec = PI / (180. * 3600.);
   double t_mean = 0., sum_dt2 = 0., dot;
   double sigma_east = 0., sigma_north = 0.;
   int i, j;

   if( n_obs < 2)
      return( -1);
   for( i = 0; i < n_obs; i++)
      t_mean += obs[i].jd;
   t_mean /= (double)n_obs;
   attr->t = t_mean;
   for( j = 0; j < 3; j++)
      {
      attr->u[j] = 0.;
      attr->eta[j] = 0.;
      attr->q[j] = 0.;
      attr->qdot[j] = 0.;
      }
   for( i = 0; i < n_obs; i++)
      {
      const double dt = obs[i].jd - t_mean;

      sum_dt2 += dt * dt;
      for( j = 0; j < 3; j++)
         {
         attr->u[j] += obs[i].vect[j];
         attr->eta[j] += obs[i].vect[j] * dt;
         attr->q[j] += obs[i].obs_posn[j];
         attr->qdot[j] += obs[i].obs_vel[j];
         }
      }
   if( sum_dt2 == 0.)         /* all observations at the same instant */
      return( -2);
   for( j = 0; j < 3; j++)
      {
      attr->u[j] /= (double)n_obs;
      attr->eta[j] /= sum_dt2;
      attr->q[j] /= (double)n_obs;
      attr->qdot[j] /= (double)n_obs;
      }
   normalize_vect3( attr->u);
   dot = dot_product( attr->u, attr->eta);
   for( j = 0; j < 3; j++)      /* enforce eta perpendicular to u */
      attr->eta[j] -= dot * attr->u[j];
                  /* Uncertainty of the fit.  Because the epochs are taken
                     relative to their mean,  the intercept and the slope
                     are uncorrelated,  and the covariance is diagonal:
                     sigma^2 / n for the line of sight and sigma^2 / sum
                     dt^2 for its rate.  We take a representative sigma per
                     direction,  which is exact when the observations are
                     equally weighted -- the usual case within one
                     tracklet.  posn_sigma_1 is taken as the east (RA)
                     uncertainty and posn_sigma_2 as the north (dec) one,
                     as elsewhere in find_orb;  either being unset falls
                     back to one arcsecond.               */
   for( i = 0; i < n_obs; i++)
      {
      const double s1 = (obs[i].posn_sigma_1 > 0. ? obs[i].posn_sigma_1 : 1.);
      const double s2 = (obs[i].posn_sigma_2 > 0. ? obs[i].posn_sigma_2 : 1.);

      sigma_east += s1 * s1;
      sigma_north += s2 * s2;
      }
   sigma_east = sqrt( sigma_east / (double)n_obs) * arcsec;
   sigma_north = sqrt( sigma_north / (double)n_obs) * arcsec;
   for( i = 0; i < 4; i++)
      for( j = 0; j < 4; j++)
         attr->covar[i][j] = 0.;
   attr->covar[0][0] = sigma_east * sigma_east / (double)n_obs;
   attr->covar[1][1] = sigma_north * sigma_north / (double)n_obs;
   attr->covar[2][2] = sigma_east * sigma_east / sum_dt2;
   attr->covar[3][3] = sigma_north * sigma_north / sum_dt2;
   return( 0);
}

/* The vectors D, E, F, G of Gronchi's eq. (DEFG),  in the frame-free form
   described at the top of this file.                                 */

static void compute_defg( const ATTRIBUTABLE *attr, double *dvect,
                     double *evect, double *fvect, double *gvect)
{
   double tvect[3];
   int i;

   vector_cross_product( dvect, attr->q, attr->u);
   vector_cross_product( evect, attr->u, attr->eta);
   vector_cross_product( fvect, attr->q, attr->eta);
   vector_cross_product( tvect, attr->u, attr->qdot);
   for( i = 0; i < 3; i++)
      fvect[i] += tvect[i];
   vector_cross_product( gvect, attr->q, attr->qdot);
}

/* Set up everything that depends only on the two attributables,  so that
   evaluating the constraint polynomials at a trial (rho1, rho2) is cheap.

   The conservation of angular momentum c1 = c2 reads

      D1 rhodot1 - D2 rhodot2 = J(rho1, rho2)

   with J = E2 rho2^2 - E1 rho1^2 + F2 rho2 - F1 rho1 + G2 - G1.  Since
   D1 x D2 is perpendicular to both D1 and D2,  projecting onto it
   eliminates rhodot1 and rhodot2 entirely and leaves

      q(rho1, rho2) = J . (D1 x D2) = 0,

   a plain quadratic -- a conic in the (rho1, rho2) plane -- with no
   cross term in rho1 rho2.  This is the whole of 'Tier 0':  the
   condition that the two attributables share an angular momentum
   direction,  obtained with nothing but vector algebra.

   To recover the radial velocities we project onto a direction
   perpendicular to one of the D's.  D2 x (D1 x D2) is perpendicular to
   D2,  so it isolates rhodot1;  D1 x (D1 x D2) isolates rhodot2.  Note
   this is the opposite assignment to the one printed in eq. (13)-(14)
   of arXiv:2111.02406, which appears to have the two labels
   interchanged;  link2_verify_angular_momentum() below checks the
   assignment numerically rather than trusting either source.        */

int link2_setup( LINK2_DATA *ld, const ATTRIBUTABLE *a1,
                                 const ATTRIBUTABLE *a2)
{
   double dcross[3], len2;
   int i;

   ld->a1 = a1;
   ld->a2 = a2;
   compute_defg( a1, ld->d1, ld->e1, ld->f1, ld->g1);
   compute_defg( a2, ld->d2, ld->e2, ld->f2, ld->g2);
   vector_cross_product( dcross, ld->d1, ld->d2);
   len2 = dot_product( dcross, dcross);
   if( len2 == 0.)      /* D1 parallel to D2;  the projection basis */
      return( -1);      /* collapses and the method has no leverage */
   memcpy( ld->dcross, dcross, 3 * sizeof( double));
   ld->dcross_len2 = len2;
   vector_cross_product( ld->d1_x_dcross, ld->d1, dcross);
   vector_cross_product( ld->d2_x_dcross, ld->d2, dcross);
   for( i = 0; i < 3; i++)
      ld->dg[i] = ld->g2[i] - ld->g1[i];
                     /* coefficients of the conic q(rho1, rho2) = 0 : */
   ld->q20 = -dot_product( ld->e1, dcross);
   ld->q10 = -dot_product( ld->f1, dcross);
   ld->q02 =  dot_product( ld->e2, dcross);
   ld->q01 =  dot_product( ld->f2, dcross);
   ld->q00 =  dot_product( ld->dg, dcross);
   return( 0);
}

/* J(rho1, rho2),  the right-hand side of the angular momentum equation. */

static void compute_j( const LINK2_DATA *ld, const double rho1,
                                 const double rho2, double *jvect)
{
   int i;

   for( i = 0; i < 3; i++)
      jvect[i] = ld->e2[i] * rho2 * rho2 - ld->e1[i] * rho1 * rho1
               + ld->f2[i] * rho2        - ld->f1[i] * rho1
               + ld->dg[i];
}

double link2_conic( const LINK2_DATA *ld, const double rho1, const double rho2)
{
   return( ld->q20 * rho1 * rho1 + ld->q10 * rho1
         + ld->q02 * rho2 * rho2 + ld->q01 * rho2 + ld->q00);
}

/* Given a trial (rho1, rho2),  recover the radial velocities and hence
   the full heliocentric state at each epoch.                         */

void link2_states( const LINK2_DATA *ld, const double rho1, const double rho2,
                                            double *state1, double *state2)
{
   double jvect[3];
   const double *u1 = ld->a1->u, *u2 = ld->a2->u;
   double rhodot1, rhodot2;
   int i;

   compute_j( ld, rho1, rho2, jvect);
   rhodot1 = dot_product( jvect, ld->d2_x_dcross) / ld->dcross_len2;
   rhodot2 = dot_product( jvect, ld->d1_x_dcross) / ld->dcross_len2;
   for( i = 0; i < 3; i++)
      {
      state1[i    ] = ld->a1->q[i] + rho1 * u1[i];
      state1[i + 3] = ld->a1->qdot[i] + rhodot1 * u1[i] + rho1 * ld->a1->eta[i];
      state2[i    ] = ld->a2->q[i] + rho2 * u2[i];
      state2[i + 3] = ld->a2->qdot[i] + rhodot2 * u2[i] + rho2 * ld->a2->eta[i];
      }
}

/* Internal consistency check:  on the conic q = 0,  the two angular
   momentum vectors must agree in all three components.  Returns the
   largest absolute discrepancy.  This is what settles the rhodot1 /
   rhodot2 assignment discussed above -- with the labels swapped, this
   returns a large number instead of a value at rounding level.       */

double link2_verify_angular_momentum( const LINK2_DATA *ld,
                              const double rho1, const double rho2)
{
   double state1[6], state2[6], c1[3], c2[3], worst = 0.;
   int i;

   link2_states( ld, rho1, rho2, state1, state2);
   vector_cross_product( c1, state1, state1 + 3);
   vector_cross_product( c2, state2, state2 + 3);
   for( i = 0; i < 3; i++)
      if( worst < fabs( c1[i] - c2[i]))
         worst = fabs( c1[i] - c2[i]);
   return( worst);
}

/* Gronchi's vector xi,  eq. (xi).  It is built from the energy and
   Laplace-Lenz differences in such a way that the non-polynomial terms
   mu / |r_j| cancel identically,  so no square roots and no auxiliary
   variable z are needed:

      xi = (|rdot2|^2 - |rdot1|^2) / 2 * (r1 x r2)
           - (rdot1 . r1) rdot1 x (r1 - r2)
           + (rdot2 . r2) rdot2 x (r1 - r2)

   We only ever evaluate it numerically,  so none of the elimination
   theory used to reduce this to a univariate degree-9 polynomial is
   required here.                                                     */

static void compute_xi( const LINK2_DATA *ld, const double rho1,
                                 const double rho2, double *xi)
{
   double state1[6], state2[6], r1_x_r2[3], dr[3], t1[3], t2[3];
   const double *r1, *r2, *v1, *v2;
   double coeff;
   int i;

   link2_states( ld, rho1, rho2, state1, state2);
   r1 = state1;   v1 = state1 + 3;
   r2 = state2;   v2 = state2 + 3;
   vector_cross_product( r1_x_r2, r1, r2);
   for( i = 0; i < 3; i++)
      dr[i] = r1[i] - r2[i];
   vector_cross_product( t1, v1, dr);
   vector_cross_product( t2, v2, dr);
   coeff = .5 * (dot_product( v2, v2) - dot_product( v1, v1));
   for( i = 0; i < 3; i++)
      xi[i] = coeff * r1_x_r2[i]
            - dot_product( v1, r1) * t1[i]
            + dot_product( v2, r2) * t2[i];
}

/* p1 = xi . u1.  Gronchi notes that the highest-degree monomials of xi
   are all multiplied by u1 x u2,  so projecting onto u1 (or u2) drops
   the total degree from 6 to 5.                                      */

double link2_p( const LINK2_DATA *ld, const double rho1, const double rho2,
                                                        const int which)
{
   double xi[3];

   compute_xi( ld, rho1, rho2, xi);
   return( dot_product( xi, which ? ld->a2->u : ld->a1->u));
}

/* Solve the conic for rho1 given rho2.  Because q has no rho1 rho2 cross
   term this is an ordinary quadratic,  with at most two real branches. */

static int solve_conic_for_rho1( const LINK2_DATA *ld, const double rho2,
                                                       double *rho1)
{
   const double c = ld->q02 * rho2 * rho2 + ld->q01 * rho2 + ld->q00;
   int n_found = 0;

   if( fabs( ld->q20) < 1e-14)         /* degenerates to a linear eqn */
      {
      if( fabs( ld->q10) > 1e-14)
         rho1[n_found++] = -c / ld->q10;
      }
   else
      {
      const double discr = ld->q10 * ld->q10 - 4. * ld->q20 * c;

      if( discr >= 0.)
         {
         const double sqrt_discr = sqrt( discr);

         rho1[n_found++] = (-ld->q10 - sqrt_discr) / (2. * ld->q20);
         rho1[n_found++] = (-ld->q10 + sqrt_discr) / (2. * ld->q20);
         }
      }
   return( n_found);
}

/* Refine a bracketed root of p1 along one branch of the conic by
   bisection on rho2,  re-solving the conic at each step so that we stay
   on the curve.                                                       */

static int refine_root( const LINK2_DATA *ld, double rho2_lo, double rho2_hi,
               const int branch, const int which_p, double *rho1, double *rho2)
{
   double r1_lo[2], r1_hi[2], p_lo;
   int i;

   if( solve_conic_for_rho1( ld, rho2_lo, r1_lo) <= branch)
      return( -1);
   p_lo = link2_p( ld, r1_lo[branch], rho2_lo, which_p);
   for( i = 0; i < 80; i++)
      {
      const double mid = (rho2_lo + rho2_hi) * .5;
      double r1_mid[2], p_mid;

      if( solve_conic_for_rho1( ld, mid, r1_mid) <= branch)
         return( -1);
      p_mid = link2_p( ld, r1_mid[branch], mid, which_p);
      if( (p_lo < 0.) == (p_mid < 0.))
         {
         rho2_lo = mid;
         p_lo = p_mid;
         }
      else
         rho2_hi = mid;
      }
   *rho2 = (rho2_lo + rho2_hi) * .5;
   if( solve_conic_for_rho1( ld, *rho2, r1_hi) <= branch)
      return( -1);
   *rho1 = r1_hi[branch];
   return( 0);
}

/* 'Tier 1':  find the solutions of the Link2 system without ever forming
   the degree-9 resultant.  We march along each branch of the conic
   q = 0 -- parameterised by rho2 -- and bracket sign changes of
   p1 = xi . u1.  Each such crossing is an intersection of the curves
   q = 0 and p1 = 0,  i.e. a solution of the linkage problem.

   Solutions are rejected here if either range is non-positive or if the
   resulting orbit is unbound;  Gronchi additionally applies a
   covariance-based compatibility test on (delta a, delta l),  which is
   not implemented yet.                                               */

int link2( LINK2_ROOT *roots, const int max_roots, const ATTRIBUTABLE *a1,
        const ATTRIBUTABLE *a2, const double rho_min, const double rho_max)
{
   LINK2_DATA ld;
   const int n_steps = 4000;
   int i, branch, n_roots = 0;

   if( link2_setup( &ld, a1, a2))
      return( -1);
   for( branch = 0; branch < 2; branch++)
      {
      double prev_p = 0., prev_rho2 = 0.;
      int have_prev = 0;

      for( i = 0; i <= n_steps; i++)
         {
         const double frac = (double)i / (double)n_steps;
         const double rho2 = rho_min * exp( frac * log( rho_max / rho_min));
         double r1[2], p;

         if( solve_conic_for_rho1( &ld, rho2, r1) <= branch || r1[branch] <= 0.)
            {
            have_prev = 0;
            continue;
            }
         p = link2_p( &ld, r1[branch], rho2, 0);
         if( have_prev && (p < 0.) != (prev_p < 0.) && n_roots < max_roots)
            {
            double found_rho1, found_rho2;

            if( !refine_root( &ld, prev_rho2, rho2, branch, 0,
                                          &found_rho1, &found_rho2))
               {
               LINK2_ROOT *root = roots + n_roots;

               root->rho1 = found_rho1;
               root->rho2 = found_rho2;
               link2_states( &ld, found_rho1, found_rho2,
                                          root->state1, root->state2);
               root->max_c_err =
                     link2_verify_angular_momentum( &ld, found_rho1, found_rho2);
               n_roots++;
               }
            }
         prev_p = p;
         prev_rho2 = rho2;
         have_prev = 1;
         }
      }
   return( n_roots);
}

/* 'Link3':  joining three attributables.

   Gronchi's second method (section 5 of arXiv:2111.02406,  from Gronchi,
   Bau & Maro 2017) needs only conservation of angular momentum.  There is
   no energy,  no Laplace-Lenz vector,  no auxiliary variable z and no
   radical anywhere in it,  which makes it markedly simpler than Link2.

   Writing c_i = c_j for each of the three pairs and projecting onto
   D_i x D_j eliminates both radial velocities,  exactly as in Link2,  and
   leaves one conic per pair,  each in only two of the three ranges:

      q_01(rho0, rho1) = 0,  q_12(rho1, rho2) = 0,  q_20(rho2, rho0) = 0

   Three quadratics in three unknowns have at most 2 * 2 * 2 = 8 common
   solutions,  which is exactly the degree of the univariate polynomial
   the paper reaches by taking two resultants in succession.  We do not
   form that polynomial:  marching in rho1 gives rho0 from the first conic
   and rho2 from the second -- each an ordinary quadratic with two
   branches -- and we look for sign changes of the third conic over the
   four branch combinations.  That avoids building resultant coefficients,
   which involve heavy cancellation.

   The system has a spurious family of solutions with c_j = 0 at every
   epoch:  an object moving in a straight line through the sun.  Those
   satisfy c_1 = c_2 = c_3 trivially and must be discarded.  They can be
   written down in closed form,  one range per attributable,  which is
   what straight_line_rho() below provides.        */

static void compute_pair( LINK3_PAIR *pair, const double *di, const double *ei,
             const double *fi, const double *gi, const double *dj,
             const double *ej, const double *fj, const double *gj)
{
   double dg[3];
   int i;

   vector_cross_product( pair->w, di, dj);
   pair->wlen2 = dot_product( pair->w, pair->w);
   vector_cross_product( pair->di_x_w, di, pair->w);
   for( i = 0; i < 3; i++)
      dg[i] = gj[i] - gi[i];
   pair->ca =  dot_product( ej, pair->w);
   pair->cc =  dot_product( fj, pair->w);
   pair->cb = -dot_product( ei, pair->w);
   pair->cd = -dot_product( fi, pair->w);
   pair->ce =  dot_product( dg, pair->w);
}

int link3_setup( LINK3_DATA *ld, const ATTRIBUTABLE *a0,
                 const ATTRIBUTABLE *a1, const ATTRIBUTABLE *a2)
{
   int k;

   ld->a[0] = a0;
   ld->a[1] = a1;
   ld->a[2] = a2;
   for( k = 0; k < 3; k++)
      compute_defg( ld->a[k], ld->d[k], ld->e[k], ld->f[k], ld->g[k]);
   for( k = 0; k < 3; k++)
      {
      const int j = (k + 1) % 3;

      compute_pair( ld->pair + k, ld->d[k], ld->e[k], ld->f[k], ld->g[k],
                                  ld->d[j], ld->e[j], ld->f[j], ld->g[j]);
      if( ld->pair[k].wlen2 == 0.)     /* D_i parallel to D_j */
         return( -1);
      }
                  /* Gronchi's condition (D1 x D2) . D3 != 0 :  without it
                     the six projections no longer imply c1 = c2 = c3.   */
   ld->triple_product = dot_product( ld->pair[0].w, ld->d[2]);
   if( ld->triple_product == 0.)
      return( -2);
   return( 0);
}

/* The conic belonging to pair 'pair_idx',  which joins attributable
   pair_idx to attributable (pair_idx + 1) mod 3.  'ri' is the range at
   the first of those,  'rj' the range at the second.        */

double link3_conic( const LINK3_DATA *ld, const int pair_idx,
                                const double ri, const double rj)
{
   const LINK3_PAIR *p = ld->pair + pair_idx;

   return( p->ca * rj * rj + p->cc * rj + p->cb * ri * ri + p->cd * ri + p->ce);
}

static int solve_quadratic( const double a, const double b, const double c,
                                                            double *roots)
{
   if( fabs( a) < 1e-300)
      {
      if( fabs( b) < 1e-300)
         return( 0);
      roots[0] = -c / b;
      return( 1);
      }
   else
      {
      const double discr = b * b - 4. * a * c;

      if( discr < 0.)
         return( 0);
      else
         {
         const double sqrt_discr = sqrt( discr);

         roots[0] = (-b - sqrt_discr) / (2. * a);
         roots[1] = (-b + sqrt_discr) / (2. * a);
         return( 2);
         }
      }
}

/* Given rho at the second attributable of a pair,  solve that pair's conic
   for rho at the first;  and vice versa.                    */

static int solve_pair_for_i( const LINK3_DATA *ld, const int pair_idx,
                                     const double rj, double *ri)
{
   const LINK3_PAIR *p = ld->pair + pair_idx;

   return( solve_quadratic( p->cb, p->cd,
                            p->ca * rj * rj + p->cc * rj + p->ce, ri));
}

static int solve_pair_for_j( const LINK3_DATA *ld, const int pair_idx,
                                     const double ri, double *rj)
{
   const LINK3_PAIR *p = ld->pair + pair_idx;

   return( solve_quadratic( p->ca, p->cc,
                            p->cb * ri * ri + p->cd * ri + p->ce, rj));
}

/* rhodot at the _second_ attributable of a pair,  which is what the
   projection onto D_i x (D_i x D_j) gives.  Each of the three radial
   velocities therefore comes from a different pair:  rhodot_1 from pair
   (0,1),  rhodot_2 from (1,2),  rhodot_0 from (2,0).

   Note that the Link3 section of arXiv:2111.02406 gives this assignment,
   while its Link2 section gives the opposite one;  the two sections
   contradict each other,  and this is the version that conserves angular
   momentum.  See the comment on link2_setup().     */

static double pair_rhodot( const LINK3_DATA *ld, const int pair_idx,
                                     const double ri, const double rj)
{
   const LINK3_PAIR *p = ld->pair + pair_idx;
   const int i = pair_idx, j = (pair_idx + 1) % 3;
   double jvect[3];
   int k;

   for( k = 0; k < 3; k++)
      jvect[k] = ld->e[j][k] * rj * rj - ld->e[i][k] * ri * ri
               + ld->f[j][k] * rj      - ld->f[i][k] * ri
               + ld->g[j][k]           - ld->g[i][k];
   return( dot_product( jvect, p->di_x_w) / p->wlen2);
}

void link3_states( const LINK3_DATA *ld, const double *rho, double *states)
{
   double rhodot[3];
   int k;

   for( k = 0; k < 3; k++)
      rhodot[(k + 1) % 3] = pair_rhodot( ld, k, rho[k], rho[(k + 1) % 3]);
   for( k = 0; k < 3; k++)
      {
      const ATTRIBUTABLE *a = ld->a[k];
      double *state = states + k * 6;
      int i;

      for( i = 0; i < 3; i++)
         {
         state[i]     = a->q[i] + rho[k] * a->u[i];
         state[i + 3] = a->qdot[i] + rhodot[k] * a->u[i] + rho[k] * a->eta[i];
         }
      }
}

/* The zero-angular-momentum ('straight line') solution for one
   attributable,  in closed form.  With
   u = q - (q.e_rho) e_rho - (q.eta) eta / |eta|^2  orthogonal to both
   e_rho and eta,  we get lambda = (qdot.u)/|u|^2 and then rho directly.
   A root matching this at all three epochs is spurious.       */

double straight_line_rho( const ATTRIBUTABLE *attr)
{
   const double eta2 = dot_product( attr->eta, attr->eta);
   const double q_dot_u = dot_product( attr->q, attr->u);
   const double q_dot_eta = dot_product( attr->q, attr->eta);
   double uvect[3], tvect[3], lambda;
   int i;

   if( eta2 == 0.)
      return( 0.);
   for( i = 0; i < 3; i++)
      uvect[i] = attr->q[i] - q_dot_u * attr->u[i]
                           - q_dot_eta * attr->eta[i] / eta2;
   lambda = dot_product( attr->qdot, uvect) / dot_product( uvect, uvect);
   for( i = 0; i < 3; i++)
      tvect[i] = lambda * attr->q[i] - attr->qdot[i];
   return( dot_product( tvect, attr->eta) / eta2);
}

static void fill_in_root( const LINK3_DATA *ld, LINK3_ROOT *root,
                                              const double *rho)
{
   double c[3][3];
   int k;

   memcpy( root->rho, rho, 3 * sizeof( double));
   link3_states( ld, rho, root->state[0]);
   for( k = 0; k < 3; k++)
      vector_cross_product( c[k], root->state[k], root->state[k] + 3);
   root->c_len = vector3_length( c[0]);
   root->max_c_err = 0.;
   for( k = 0; k < 3; k++)
      {
      const int j = (k + 1) % 3;
      int i;

      for( i = 0; i < 3; i++)
         if( root->max_c_err < fabs( c[k][i] - c[j][i]))
            root->max_c_err = fabs( c[k][i] - c[j][i]);
      }
}

/* Residual of the third conic,  for a given rho1 and a choice of branch
   on each of the first two conics.  Returns zero if that branch pair does
   not exist at this rho1.                                    */

static int link3_residual( const LINK3_DATA *ld, const double rho1,
          const int branch0, const int branch2, double *resid, double *rho)
{
   double r0[2], r2[2];

   if( solve_pair_for_i( ld, 0, rho1, r0) <= branch0)
      return( 0);
   if( solve_pair_for_j( ld, 1, rho1, r2) <= branch2)
      return( 0);
   rho[0] = r0[branch0];
   rho[1] = rho1;
   rho[2] = r2[branch2];
   *resid = link3_conic( ld, 2, rho[2], rho[0]);
   return( 1);
}

int link3( LINK3_ROOT *roots, const int max_roots, const ATTRIBUTABLE *a0,
        const ATTRIBUTABLE *a1, const ATTRIBUTABLE *a2,
        const double rho_min, const double rho_max)
{
   LINK3_DATA ld;
   const int n_steps = 8000;
   int i, branch0, branch2, n_roots = 0;

   if( link3_setup( &ld, a0, a1, a2))
      return( -1);
   for( branch0 = 0; branch0 < 2; branch0++)
      for( branch2 = 0; branch2 < 2; branch2++)
         {
         double prev_resid = 0., prev_rho1 = 0.;
         int have_prev = 0;

         for( i = 0; i <= n_steps; i++)
            {
            const double frac = (double)i / (double)n_steps;
            const double rho1 = rho_min * exp( frac * log( rho_max / rho_min));
            double resid, rho[3];

            if( !link3_residual( &ld, rho1, branch0, branch2, &resid, rho))
               {
               have_prev = 0;
               continue;
               }
            if( have_prev && (resid < 0.) != (prev_resid < 0.)
                                          && n_roots < max_roots)
               {
               double lo = prev_rho1, hi = rho1;
               int iter;

               for( iter = 0; iter < 100; iter++)
                  {
                  const double mid = (lo + hi) * .5;
                  double mid_resid, mid_rho[3];

                  if( !link3_residual( &ld, mid, branch0, branch2,
                                                &mid_resid, mid_rho))
                     break;
                  if( (mid_resid < 0.) == (prev_resid < 0.))
                     lo = mid;
                  else
                     hi = mid;
                  }
               if( link3_residual( &ld, (lo + hi) * .5, branch0, branch2,
                                                        &resid, rho))
                  {
                  LINK3_ROOT troot;
                  int is_dup = 0, k;

                  for( k = 0; k < n_roots; k++)
                     if( fabs( roots[k].rho[0] - rho[0]) < 1e-4
                      && fabs( roots[k].rho[1] - rho[1]) < 1e-4
                      && fabs( roots[k].rho[2] - rho[2]) < 1e-4)
                        is_dup = 1;
                  fill_in_root( &ld, &troot, rho);
                     /* Where the two branches of a conic merge -- at a zero
                        of its discriminant -- the branch ordering swaps and
                        the residual jumps.  Bisection then converges on that
                        discontinuity rather than on a root.  A true root
                        conserves angular momentum to rounding,  so test that
                        rather than the residual,  which has no natural
                        scale.                                          */
                  if( !is_dup && troot.max_c_err < 1e-8 * troot.c_len + 1e-12)
                     roots[n_roots++] = troot;
                  }
               }
            prev_resid = resid;
            prev_rho1 = rho1;
            have_prev = 1;
            }
         }
   return( n_roots);
}

/* The compatibility test.

   Conservation laws cannot pin down every element.  Link2 imposes angular
   momentum,  energy and the Laplace-Lenz vector,  which leaves only the
   semimajor axis and the mean anomaly free;  Link3 imposes angular
   momentum alone,  which additionally leaves the argument of perihelion
   free.  If the tracklets really do belong to one object those leftover
   differences must vanish,  so they make a natural test of the
   identification.

   The catch is that they are not small in absolute terms.  The mean
   anomaly is exactly what the integrals fail to constrain,  and in the
   paper's own accepted Link2 solution the two epochs disagree by about
   ten degrees in it.  So a raw threshold is useless;  what matters is the
   size of the discrepancy relative to its own uncertainty.  We therefore
   propagate the attributable covariances through to Delta and use the
   resulting quadratic form.                                        */

void attributable_basis( const ATTRIBUTABLE *attr, double *east, double *north)
{
   const double pole[3] = { 0., 0., 1. };
   const double xaxis[3] = { 1., 0., 0. };
   double len;
   int i;

   vector_cross_product( east, pole, attr->u);
   len = vector3_length( east);
   if( len < 1e-8)         /* looking along the pole;  any basis will do */
      {
      vector_cross_product( east, xaxis, attr->u);
      len = vector3_length( east);
      }
   for( i = 0; i < 3; i++)
      east[i] /= len;
   vector_cross_product( north, attr->u, east);
}

static void perturb_attributable( ATTRIBUTABLE *out, const ATTRIBUTABLE *in,
                                                     const double *delta)
{
   double east[3], north[3], dot;
   int i;

   *out = *in;
   attributable_basis( in, east, north);
   for( i = 0; i < 3; i++)
      {
      out->u[i]   += delta[0] * east[i] + delta[1] * north[i];
      out->eta[i] += delta[2] * east[i] + delta[3] * north[i];
      }
   normalize_vect3( out->u);
   dot = dot_product( out->u, out->eta);
   for( i = 0; i < 3; i++)
      out->eta[i] -= dot * out->u[i];
}

/* Solve the n x n system a x = b in place by Gauss-Jordan elimination with
   partial pivoting.  n is at most six here.  Returns nonzero if singular. */

static int solve_linear( double *a, double *b, const int n)
{
   int i, j, k;

   for( i = 0; i < n; i++)
      {
      int best = i;
      double pivot;

      for( j = i + 1; j < n; j++)
         if( fabs( a[j * n + i]) > fabs( a[best * n + i]))
            best = j;
      if( best != i)
         {
         double temp;

         for( j = 0; j < n; j++)
            {
            temp = a[i * n + j];
            a[i * n + j] = a[best * n + j];
            a[best * n + j] = temp;
            }
         temp = b[i];
         b[i] = b[best];
         b[best] = temp;
         }
      pivot = a[i * n + i];
      if( pivot == 0.)
         return( -1);
      for( j = 0; j < n; j++)
         a[i * n + j] /= pivot;
      b[i] /= pivot;
      for( k = 0; k < n; k++)
         if( k != i)
            {
            const double factor = a[k * n + i];

            for( j = 0; j < n; j++)
               a[k * n + j] -= factor * a[i * n + j];
            b[k] -= factor * b[i];
            }
      }
   return( 0);
}

static double wrap_angle( double angle)
{
   angle = fmod( angle, PI + PI);
   if( angle > PI)
      angle -= PI + PI;
   if( angle < -PI)
      angle += PI + PI;
   return( angle);
}

static int state_to_elements( ELEMENTS *elem, const double *state,
                                              const double epoch)
{
   elem->gm = SOLAR_GM;
   elem->epoch = epoch;
   elem->central_obj = 0;
   calc_classical_elements( elem, state, epoch, 1);
   if( elem->ecc >= 1. || elem->major_axis <= 0.)
      return( -1);            /* unbound;  no mean anomaly to compare */
   return( 0);
}

/* Newton refinement of the ranges,  used to follow a root as the
   attributables are perturbed.  For Link2 the two constraints are the
   conic and xi . u1;  for Link3 they are the three conics.  The radial
   velocities are not unknowns here -- they were eliminated analytically --
   so the systems are 2 x 2 and 3 x 3 rather than the 4 x 4 and 6 x 6 of
   the paper.                                                       */

static int link_residuals( double *resid, const ATTRIBUTABLE *attr,
                              const int n_attr, const double *rho)
{
   if( n_attr == 2)
      {
      LINK2_DATA ld;

      if( link2_setup( &ld, attr, attr + 1))
         return( -1);
      resid[0] = link2_conic( &ld, rho[0], rho[1]);
      resid[1] = link2_p( &ld, rho[0], rho[1], 0);
      }
   else
      {
      LINK3_DATA ld;
      int k;

      if( link3_setup( &ld, attr, attr + 1, attr + 2))
         return( -1);
      for( k = 0; k < 3; k++)
         resid[k] = link3_conic( &ld, k, rho[k], rho[(k + 1) % 3]);
      }
   return( 0);
}

static int refine_ranges( const ATTRIBUTABLE *attr, const int n_attr,
                                                    double *rho)
{
   const int n = n_attr;      /* one constraint per unknown range */
   int iter;

   for( iter = 0; iter < 50; iter++)
      {
      double f[3], jac[9], step[3], worst = 0., scale;
      int i, j;

      if( link_residuals( f, attr, n_attr, rho))
         return( -1);
      for( j = 0; j < n; j++)
         {
         const double h = fabs( rho[j]) * 1e-6 + 1e-10;
         double saved = rho[j], fp[3], fm[3];

         rho[j] = saved + h;
         if( link_residuals( fp, attr, n_attr, rho))
            return( -1);
         rho[j] = saved - h;
         if( link_residuals( fm, attr, n_attr, rho))
            return( -1);
         rho[j] = saved;
         for( i = 0; i < n; i++)
            jac[i * n + j] = (fp[i] - fm[i]) / (2. * h);
         }
      for( i = 0; i < n; i++)
         step[i] = -f[i];
      if( solve_linear( jac, step, n))
         return( -2);
                  /* Damp the step.  Link3 in particular is stiff enough
                     that an undamped Newton step can overshoot into
                     negative ranges,  or onto a different branch of a
                     conic.  Limit each range to a half of its current
                     value per iteration,  and never let one go
                     non-positive.                        */
      scale = 1.;
      for( i = 0; i < n; i++)
         if( fabs( step[i]) > .5 * rho[i])
            {
            const double limit = .5 * rho[i] / fabs( step[i]);

            if( scale > limit)
               scale = limit;
            }
      for( i = 0; i < n; i++)
         {
         const double applied = scale * step[i];

         rho[i] += applied;
         if( worst < fabs( applied))
            worst = fabs( applied);
         if( rho[i] <= 0.)
            return( -3);
         }
      if( worst < 1e-13)
         return( 0);
      }
   return( -4);         /* failed to converge */
}

/* The leftover element differences.  Link2 compares the two epochs;  Link3
   compares epochs 1 and 3 against epoch 2,  which is Gronchi's choice of
   reference.                                            */

static int compute_delta( double *delta, const ATTRIBUTABLE *attr,
                              const int n_attr, const double *rho)
{
   ELEMENTS elem[3];
   double t[3], n_ref;
   int k;

   if( n_attr == 2)
      {
      LINK2_DATA ld;
      double state[2][6];

      if( link2_setup( &ld, attr, attr + 1))
         return( -1);
      link2_states( &ld, rho[0], rho[1], state[0], state[1]);
      for( k = 0; k < 2; k++)
         {
         t[k] = attr[k].t - rho[k] / AU_PER_DAY;   /* light-time corrected */
         if( state_to_elements( elem + k, state[k], t[k]))
            return( -1);
         }
      n_ref = GAUSS_K / (elem[1].major_axis * sqrt( elem[1].major_axis));
      delta[0] = elem[0].major_axis - elem[1].major_axis;
      delta[1] = wrap_angle( elem[0].mean_anomaly - elem[1].mean_anomaly
                                          - n_ref * (t[0] - t[1]));
      }
   else
      {
      LINK3_DATA ld;
      double states[18];
      int out = 0;

      if( link3_setup( &ld, attr, attr + 1, attr + 2))
         return( -1);
      link3_states( &ld, rho, states);
      for( k = 0; k < 3; k++)
         {
         t[k] = attr[k].t - rho[k] / AU_PER_DAY;
         if( state_to_elements( elem + k, states + k * 6, t[k]))
            return( -1);
         }
      n_ref = GAUSS_K / (elem[1].major_axis * sqrt( elem[1].major_axis));
      for( k = 0; k < 3; k += 2)       /* epochs 0 and 2,  against 1 */
         {
         delta[out++] = elem[k].major_axis - elem[1].major_axis;
         delta[out++] = wrap_angle( elem[k].arg_per - elem[1].arg_per);
         delta[out++] = wrap_angle( elem[k].mean_anomaly - elem[1].mean_anomaly
                                             - n_ref * (t[k] - t[1]));
         }
      }
   return( 0);
}

/* Propagate the attributable covariances onto Delta and form the quadratic
   form.  The Jacobian is taken by central differences of the whole map
   from attributables to Delta,  re-converging the ranges at each perturbed
   point;  that is the total derivative,  including the implicit dependence
   of the ranges on the data,  without having to assemble the partials of
   the constraint system by hand.  Steps are one standard deviation of the
   component being varied,  which keeps them well clear of rounding while
   staying in the linear regime.                                     */

static int compatibility( double *chi2, double *delta_out,
              const ATTRIBUTABLE *attr, const int n_attr, const double *rho)
{
   const int n_delta = (n_attr == 2 ? 2 : 6);
   const int n_params = 4 * n_attr;
   double delta0[6], jac[6 * 12], gamma[36], rhs[6];
   double rho0[3];
   int i, j, k, b;

   memcpy( rho0, rho, n_attr * sizeof( double));
   if( compute_delta( delta0, attr, n_attr, rho0))
      return( -1);
   if( delta_out)
      memcpy( delta_out, delta0, n_delta * sizeof( double));
   for( j = 0; j < n_params; j++)
      {
      const int which = j / 4, comp = j % 4;
      const double sigma = sqrt( attr[which].covar[comp][comp]);
      ATTRIBUTABLE perturbed[3];
      double offset[4], dplus[6], dminus[6], trial[3];
      double h = FD_STEP_IN_SIGMAS * sigma;
      int sign, attempt, converged = 0;

      for( i = 0; i < n_delta; i++)
         jac[i * n_params + j] = 0.;
      if( sigma <= 0.)           /* no uncertainty quoted;  no sensitivity */
         continue;
                  /* If a perturbed system will not reconverge -- the root
                     having moved onto another branch -- shrink the step
                     and try again rather than giving up on the whole
                     test.                                         */
      for( attempt = 0; attempt < 6 && !converged; attempt++, h *= .25)
         {
         converged = 1;
         for( sign = -1; sign <= 1 && converged; sign += 2)
            {
            for( k = 0; k < 4; k++)
               offset[k] = (k == comp ? sign * h : 0.);
            for( k = 0; k < n_attr; k++)
               if( k == which)
                  perturb_attributable( perturbed + k, attr + k, offset);
               else
                  perturbed[k] = attr[k];
            memcpy( trial, rho0, n_attr * sizeof( double));
            if( refine_ranges( perturbed, n_attr, trial)
                     || compute_delta( sign > 0 ? dplus : dminus, perturbed,
                                                          n_attr, trial))
               converged = 0;
            }
         if( converged)
            for( i = 0; i < n_delta; i++)
               jac[i * n_params + j] = (dplus[i] - dminus[i]) / (2. * h);
         }
      if( !converged)
         return( -2);
      }
                  /* Gamma_Delta = J Gamma_A J^T,  with Gamma_A block
                     diagonal:  one 4 x 4 block per attributable.    */
   for( i = 0; i < n_delta * n_delta; i++)
      gamma[i] = 0.;
   for( b = 0; b < n_attr; b++)
      for( i = 0; i < n_delta; i++)
         for( j = 0; j < n_delta; j++)
            {
            double sum = 0.;
            int m, n;

            for( m = 0; m < 4; m++)
               for( n = 0; n < 4; n++)
                  sum += jac[i * n_params + b * 4 + m]
                       * attr[b].covar[m][n]
                       * jac[j * n_params + b * 4 + n];
            gamma[i * n_delta + j] += sum;
            }
   memcpy( rhs, delta0, n_delta * sizeof( double));
   if( solve_linear( gamma, rhs, n_delta))
      return( -3);
   *chi2 = 0.;
   for( i = 0; i < n_delta; i++)
      *chi2 += delta0[i] * rhs[i];
   return( 0);
}

int link2_compatibility( double *chi2, double *delta, const ATTRIBUTABLE *a1,
                              const ATTRIBUTABLE *a2, const double *rho)
{
   ATTRIBUTABLE attr[2];

   attr[0] = *a1;
   attr[1] = *a2;
   return( compatibility( chi2, delta, attr, 2, rho));
}

int link3_compatibility( double *chi2, double *delta, const ATTRIBUTABLE *a0,
       const ATTRIBUTABLE *a1, const ATTRIBUTABLE *a2, const double *rho)
{
   ATTRIBUTABLE attr[3];

   attr[0] = *a0;
   attr[1] = *a1;
   attr[2] = *a2;
   return( compatibility( chi2, delta, attr, 3, rho));
}

/* Driving the above from a set of observations.

   Everything up to this point works on attributables.  This turns a plain
   list of observations into them,  runs whichever method the number of
   tracklets allows,  and returns the best-scoring orbit as a state vector.

   Note that nothing here calls into the rest of find_orb:  the whole file
   links against the 'lunar' library alone,  which is what lets ki_test
   check the mathematics without dragging in an ephemeris,  a config
   directory or an integrator.  Configuration and reporting are therefore
   the caller's business -- see initial_orbit() in orb_func.cpp.       */

/* Split the observations into tracklets on gaps larger than
   'max_tracklet_gap' days.  A tracklet is a handful of observations over
   minutes to hours;  the next one is typically a night or longer away,
   so almost any threshold between about an hour and a day gives the same
   answer.                                                            */

static int find_tracklets( const OBSERVE FAR *obs, const int n_obs,
       int *starts, int *lengths, const int max_tracklets,
       const double max_tracklet_gap)
{
   int i, n = 0;

   if( !n_obs)
      return( 0);
   starts[0] = 0;
   for( i = 1; i < n_obs && n < max_tracklets; i++)
      if( obs[i].jd - obs[i - 1].jd > max_tracklet_gap)
         {
         lengths[n] = i - starts[n];
         n++;
         if( n < max_tracklets)
            starts[n] = i;
         }
   if( n < max_tracklets)
      {
      lengths[n] = n_obs - starts[n];
      n++;
      }
   return( n);
}

static int is_bound( const double *state)
{
   ELEMENTS elem;

   return( !state_to_elements( &elem, state, J2000));
}

/* Score a candidate:  reject it outright if any range is non-positive,  if
   any epoch gives an unbound orbit,  or if the angular momentum vanishes
   (the straight-line solution),  and otherwise return the compatibility
   chi^2.  Lower is better;  a negative return means 'rejected'.      */

#define REJECTED (-1.)

static double score_link3_root( const LINK3_ROOT *root, const ATTRIBUTABLE *attr)
{
   double chi2;
   int k;

   if( root->c_len < 1e-8)                /* straight-line solution */
      return( REJECTED);
   for( k = 0; k < 3; k++)
      if( root->rho[k] <= 0. || !is_bound( root->state[k]))
         return( REJECTED);
   if( link3_compatibility( &chi2, NULL, attr, attr + 1, attr + 2, root->rho))
      return( REJECTED);
   return( chi2);
}

static double score_link2_root( const LINK2_ROOT *root, const ATTRIBUTABLE *attr)
{
   double chi2, rho[2];

   if( root->rho1 <= 0. || root->rho2 <= 0.)
      return( REJECTED);
   if( !is_bound( root->state1) || !is_bound( root->state2))
      return( REJECTED);
   rho[0] = root->rho1;
   rho[1] = root->rho2;
   if( link2_compatibility( &chi2, NULL, attr, attr + 1, rho))
      return( REJECTED);
   return( chi2);
}

/* Try one triple of tracklets,  and keep it if it beats what we have. */

#define MAX_WINDOWS_TRIED 16

static void try_link3_window( KEPLERIAN_LINK_RESULT *result,
       const ATTRIBUTABLE *attr, const int *starts, const int *lengths,
       const int a, const int b, const int c, const double rho_max)
{
   LINK3_ROOT roots[20];
   ATTRIBUTABLE chosen[3];
   int i, n_roots;

   chosen[0] = attr[a];
   chosen[1] = attr[b];
   chosen[2] = attr[c];
   n_roots = link3( roots, 20, chosen, chosen + 1, chosen + 2, .05, rho_max);
   for( i = 0; i < n_roots; i++)
      {
      const double chi2 = score_link3_root( roots + i, chosen);

      if( chi2 >= 0. && (!result->epoch || chi2 < result->chi2))
         {
         result->chi2 = chi2;
         result->n_tracklets = 3;
         result->first_obs = starts[a];
         result->n_obs_used = starts[c] + lengths[c] - starts[a];
         memcpy( result->orbit, roots[i].state[0], 6 * sizeof( double));
         result->epoch = chosen[0].t - roots[i].rho[0] / AU_PER_DAY;
         }
      }
}

static void try_link2_window( KEPLERIAN_LINK_RESULT *result,
       const ATTRIBUTABLE *attr, const int *starts, const int *lengths,
       const int a, const int c, const double rho_max)
{
   LINK2_ROOT roots[20];
   ATTRIBUTABLE chosen[2];
   int i, n_roots;

   chosen[0] = attr[a];
   chosen[1] = attr[c];
   n_roots = link2( roots, 20, chosen, chosen + 1, .05, rho_max);
   for( i = 0; i < n_roots; i++)
      {
      const double chi2 = score_link2_root( roots + i, chosen);

      if( chi2 >= 0. && (!result->epoch || chi2 < result->chi2))
         {
         result->chi2 = chi2;
         result->n_tracklets = 2;
         result->first_obs = starts[a];
         result->n_obs_used = starts[c] + lengths[c] - starts[a];
         memcpy( result->orbit, roots[i].state1, 6 * sizeof( double));
         result->epoch = chosen[0].t - roots[i].rho1 / AU_PER_DAY;
         }
      }
}

int keplerian_link_orbit( KEPLERIAN_LINK_RESULT *result,
             const OBSERVE FAR *obs, const int n_obs,
             const double max_tracklet_gap, const double max_span,
             const double rho_max)
{
   int starts[MAX_KI_TRACKLETS], lengths[MAX_KI_TRACKLETS];
   int t_start[MAX_KI_TRACKLETS], t_len[MAX_KI_TRACKLETS];
   ATTRIBUTABLE attr[MAX_KI_TRACKLETS];
   int n_tracklets, i, n_attr = 0, n_windows = 0;

   memset( result, 0, sizeof( KEPLERIAN_LINK_RESULT));
   n_tracklets = find_tracklets( obs, n_obs, starts, lengths,
                                 MAX_KI_TRACKLETS, max_tracklet_gap);
   for( i = 0; i < n_tracklets; i++)
      if( lengths[i] >= 2)          /* need two obs to get a rate of motion */
         if( !compute_attributable( attr + n_attr, obs + starts[i], lengths[i]))
            {
            t_start[n_attr] = starts[i];
            t_len[n_attr] = lengths[i];
            n_attr++;
            }
   if( n_attr < 2)
      return( -1);
                  /* Slide a window of at most 'max_span' days along the
                     arc.  Within each placement we take the two ends and
                     the middle:  adjacent tracklets would throw away the
                     separation these methods exist to exploit,  while the
                     whole arc would break the two-body assumption they
                     rest on.  Anchors are spread evenly so that a long
                     arc costs no more than a short one.       */
   for( i = 0; i < n_attr - 1 && n_windows < MAX_WINDOWS_TRIED; i++)
      {
      const int step = (n_attr <= MAX_WINDOWS_TRIED ? 1
                              : n_attr / MAX_WINDOWS_TRIED);
      const int a = (i * step < n_attr ? i * step : n_attr - 1);
      int c = a;

      while( c + 1 < n_attr && attr[c + 1].t - attr[a].t <= max_span)
         c++;
      if( c == a)             /* nothing else within the span */
         continue;
      n_windows++;
      if( c - a >= 2)
         try_link3_window( result, attr, t_start, t_len, a, (a + c) / 2, c,
                                                            rho_max);
      else
         try_link2_window( result, attr, t_start, t_len, a, c, rho_max);
      if( i * step >= n_attr - 1)
         break;
      }
                  /* If Link3 never produced anything acceptable,  fall
                     back on Link2 over the widest admissible baseline. */
   if( !result->epoch)
      {
      int c = 0;

      while( c + 1 < n_attr && attr[c + 1].t - attr[0].t <= max_span)
         c++;
      if( c)
         try_link2_window( result, attr, t_start, t_len, 0, c, rho_max);
      }
   return( result->epoch ? 0 : -2);
}
