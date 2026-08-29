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
   double covar[4][4];     /* uncertainty of (u, eta);  see below */
} ATTRIBUTABLE;

/* 'covar' is the covariance of the attributable in the local orthonormal
   basis (east, north) perpendicular to the line of sight,  ordered

      (du.east, du.north, deta.east, deta.north)

   in radians and radians/day.  Those are the same four quantities as
   Gronchi's (alpha cos delta, delta, alphadot cos delta, deltadot),  but
   defined without reference to a pole.  attributable_basis() returns the
   two basis vectors.                                              */

void attributable_basis( const ATTRIBUTABLE *attr, double *east, double *north);

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

/* 'Link3' joins three attributables using conservation of angular momentum
   alone.  Each of the three pairs (0,1), (1,2), (2,0) contributes one conic
   in two of the ranges;  the three conics in three unknowns admit at most
   2 * 2 * 2 = 8 solutions,  which is the degree of the univariate
   polynomial Gronchi obtains by a double resultant.               */

typedef struct
{
   double w[3];              /* D_i x D_j */
   double wlen2;             /* |D_i x D_j|^2 */
   double di_x_w[3];         /* D_i x (D_i x D_j);  yields rhodot_j */
   double ca, cb, cc, cd, ce;   /* conic: ca rj^2 + cc rj + cb ri^2 + cd ri + ce */
} LINK3_PAIR;

typedef struct
{
   const ATTRIBUTABLE *a[3];
   double d[3][3], e[3][3], f[3][3], g[3][3];
   LINK3_PAIR pair[3];       /* pair[k] joins attributable k to k + 1 mod 3 */
   double triple_product;    /* D1 x D2 . D3;  must be nonzero */
} LINK3_DATA;

typedef struct
{
   double rho[3];
   double state[3][6];       /* heliocentric,  AU and AU/day */
   double max_c_err;         /* largest |c_i - c_j| */
   double c_len;             /* |c|;  ~0 for the spurious straight-line root */
} LINK3_ROOT;

int compute_attributable( ATTRIBUTABLE *attr, const OBSERVE FAR *obs,
                                                    const int n_obs);
int link3_setup( LINK3_DATA *ld, const ATTRIBUTABLE *a0,
                 const ATTRIBUTABLE *a1, const ATTRIBUTABLE *a2);
double link3_conic( const LINK3_DATA *ld, const int pair_idx,
                                const double ri, const double rj);
void link3_states( const LINK3_DATA *ld, const double *rho, double *states);
double straight_line_rho( const ATTRIBUTABLE *attr);
int link3( LINK3_ROOT *roots, const int max_roots, const ATTRIBUTABLE *a0,
        const ATTRIBUTABLE *a1, const ATTRIBUTABLE *a2,
        const double rho_min, const double rho_max);

/* The compatibility test that selects among the solutions.  The Keplerian
   integrals leave some elements unconstrained:  Link2 matches everything
   except the semimajor axis and the mean anomaly,  Link3 (using angular
   momentum alone) leaves the semimajor axis,  argument of perihelion and
   mean anomaly free.  If the attributables really do belong to one object,
   those leftover differences must vanish.  We form them into a vector
   Delta -- two components for Link2,  six for Link3 -- propagate the
   attributable covariances onto it,  and return

      chi2 = Delta . Gamma_Delta^-1 . Delta

   Gronchi's control is chi2 <= chi_max^2,  with chi_max chosen from
   simulation.  Note that Delta is not small in absolute terms:  the mean
   anomaly is precisely what the integrals do not constrain,  so a
   perfectly good solution can be tens of degrees out and still pass.  The
   covariance is doing all the work.

   'delta' receives the raw Delta vector (radians and AU) if not NULL.
   Returns 0 on success,  nonzero if the elements or the covariance could
   not be formed (an unbound solution,  say,  or a singular Gamma).   */

int link2_compatibility( double *chi2, double *delta, const ATTRIBUTABLE *a1,
                              const ATTRIBUTABLE *a2, const double *rho);
int link3_compatibility( double *chi2, double *delta, const ATTRIBUTABLE *a0,
       const ATTRIBUTABLE *a1, const ATTRIBUTABLE *a2, const double *rho);

/* Driving all of the above from a list of observations:  split them into
   tracklets,  pick a well-separated triple,  run Link3 (or Link2 if only
   two tracklets are available),  and keep the surviving orbit with the
   smallest compatibility chi^2.

   The separation is deliberately bounded.  These methods can link across
   arbitrarily long gaps in principle,  but only because they assume pure
   two-body motion,  and that assumption is what a long baseline destroys:
   the node regresses,  the inclination oscillates,  and the conserved
   quantities we are matching stop being conserved.  Gronchi's own test
   cases span 921 and 564 days.  Handing this the first and last tracklets
   of a twenty-five year arc produces a confidently wrong answer,  so we
   cap the baseline and slide the window along the arc instead,  keeping
   whichever placement gives the best chi^2.

   'result->epoch' is the epoch of the returned state vector -- corrected
   for light time,  so it is when the light left the object -- or zero if
   no orbit was found.  'first_obs' and 'n_obs_used' report which
   observations went into it,  so that the caller can score the orbit
   against those rather than against an entire multi-apparition arc no
   two-body orbit could ever fit.                                     */

#define MAX_KI_TRACKLETS 1024

typedef struct
{
   double epoch;              /* zero if nothing was found */
   double orbit[6];           /* heliocentric,  AU and AU/day */
   double chi2;               /* compatibility chi^2 of the chosen root */
   double baseline;           /* days between the first and last tracklet */
   int n_tracklets;           /* three for Link3,  two for Link2 */
   int first_obs, n_obs_used; /* the span from the first to the last */
   int trk_start[3], trk_len[3];   /* the linked tracklets themselves */
} KEPLERIAN_LINK_RESULT;

/* 'chi2_max' is the threshold below which a linkage counts as acceptable.
   Among acceptable windows we take the one with the longest baseline,  not
   the one with the smallest chi^2.  Those are very different choices:  two
   tracklets a night apart are compatible with almost any orbit and score a
   superb chi^2 while telling us nothing,  so selecting on chi^2 alone picks
   out precisely the linkages with no leverage.  Ranking by baseline and
   using chi^2 only as an admission test asks the question we actually care
   about -- what is the longest arc these observations will support?     */

int keplerian_link_orbit( KEPLERIAN_LINK_RESULT *result,
             const OBSERVE FAR *obs, const int n_obs,
             const double max_tracklet_gap, const double max_span,
             const double rho_max, const double chi2_max);
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
