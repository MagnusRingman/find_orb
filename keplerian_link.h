/* keplerian_link.h: linkage of too-short arcs via the Keplerian integrals

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

#ifndef KEPLERIAN_LINK_H_INCLUDED
#define KEPLERIAN_LINK_H_INCLUDED

/* An 'attributable':  a tracklet compressed to a line of sight and its
   time derivative at the mean epoch,  plus the observer's heliocentric
   state at that epoch.  See keplerian_link.cpp for why we store the
   unit vector and its derivative rather than Gronchi's
   (alpha, delta, alphadot, deltadot).   All vectors are in whatever
   frame the observations use (ecliptic,  for find_orb),  with distances
   in AU and times in days.                                          */

typedef struct
{
   double t;               /* mean epoch of the tracklet,  JD */
   double u[3];            /* unit vector toward the object */
   double eta[3];          /* du/dt,  radians/day;  perpendicular to u */
   double q[3];            /* heliocentric observer position,  AU */
   double qdot[3];         /* heliocentric observer velocity,  AU/day */
} ATTRIBUTABLE;

/* Everything that depends only on the pair of attributables,  so that
   evaluating the constraint polynomials at a trial (rho1, rho2) is
   cheap.                                                            */

typedef struct
{
   const ATTRIBUTABLE *a1, *a2;
   double d1[3], e1[3], f1[3], g1[3];
   double d2[3], e2[3], f2[3], g2[3];
   double dg[3];                    /* G2 - G1 */
   double dcross[3];                /* D1 x D2 */
   double dcross_len2;              /* |D1 x D2|^2 */
   double d1_x_dcross[3];           /* D1 x (D1 x D2);  isolates rhodot2 */
   double d2_x_dcross[3];           /* D2 x (D1 x D2);  isolates rhodot1 */
   double q20, q10, q02, q01, q00;  /* coefficients of the conic q = 0 */
} LINK2_DATA;

typedef struct
{
   double rho1, rho2;
   double state1[6], state2[6];     /* heliocentric,  AU and AU/day */
   double max_c_err;                /* |c1 - c2|,  should be at rounding level */
} LINK2_ROOT;

int compute_attributable( ATTRIBUTABLE *attr, const OBSERVE FAR *obs,
                                                    const int n_obs);
int link2_setup( LINK2_DATA *ld, const ATTRIBUTABLE *a1,
                                 const ATTRIBUTABLE *a2);
double link2_conic( const LINK2_DATA *ld, const double rho1, const double rho2);
double link2_p( const LINK2_DATA *ld, const double rho1, const double rho2,
                                                    const int which);
void link2_states( const LINK2_DATA *ld, const double rho1, const double rho2,
                                        double *state1, double *state2);
double link2_verify_angular_momentum( const LINK2_DATA *ld,
                              const double rho1, const double rho2);
int link2( LINK2_ROOT *roots, const int max_roots, const ATTRIBUTABLE *a1,
        const ATTRIBUTABLE *a2, const double rho_min, const double rho_max);

#endif        /* #ifndef KEPLERIAN_LINK_H_INCLUDED */
