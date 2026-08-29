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

/* Least-squares fit of a straight line to each component of the observed
   unit vectors,  giving the line of sight and its rate at the mean epoch.
   For a real tracklet the arc is short enough that a linear fit is ample;
   Gronchi allows a quadratic fit,  which would matter only for a fast
   mover.  We renormalise u and force eta perpendicular to it,  since
   |u| = 1 implies u . du/dt = 0.                                     */

int compute_attributable( ATTRIBUTABLE *attr, const OBSERVE FAR *obs,
                                                        const int n_obs)
{
   double t_mean = 0., sum_dt2 = 0., dot;
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
