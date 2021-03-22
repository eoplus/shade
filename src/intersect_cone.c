
/* Line segment and cone intersection:
 *
 * In the folling, all quantities are vector quantities (x, y, z). "pos" is the 
 * vector cordinates of an arbitrary point, "ppos" is the starting position of 
 * the ray, "cdir" is the vector of cosine directions of the ray, "s" is the
 * path length of the ray. For the cone, "vertex" is the vertex of the cone, 
 * "theta" is the halfangle of the cone and "axis" is a unit vector giving the 
 * orientation of the cone. "phi" is the azimuthal angle.
 *
 * The line equation is defined as:
 * pos = ppos + cdir * s                                                     (1)
 *
 * The cone's surface equation is given by:
 * pos = vertex +- (dir · axis) * length,                                    (2)
 *
 * with "dir" is given by the polar coordinates {theta, phi} and phi being in 
 * the range [0, 2*pi). We will restrict our self to a single sided cone, taken 
 * here as the positive part with regard to its vertex. Any point on this cone's 
 * surface must satisfy:
 * (pos - vertex) · axis = ||pos - vertex|| * cos(theta)                     (3)
 *
 * Intersections with the surface of the cone are "pos" equal for the line and 
 * for the cone's surface. By squaring Eq. 3 and substituting pos for the line 
 * equation we have:
 * ((ppos + cdir * s - vertex) · axis)^2 = 
 *   ||ppos + cdir * s - vertex||^2 * cos(theta)^2 ->
 *
 * ((s * cdir + ppos - vertex) · axis)^2 = 
 *   ||s * cdir + ppos - vertex||^2 * cos(theta)^2 ->
 * 
 * ((s * cdir + ppos_v) · axis)^2 = ||s * cdir + ppos_v||^2 * mu^2,          (4)
 * 
 * Where we are using ppos_v to represent ppos centered by the vertex and mu
 * to represent cos(theta). Expanding the left hand side we have:
 * ((s * cdir + ppos_v) · axis)^2 ->
 *
 * (s * cdir · axis + ppos_v · axis)^2 ->
 *
 * s^2 * (cdir · axis)^2 + 2 * s * (cdir · axis) * (ppos_v · axis) + 
 *   (ppos_v · axis)^2 ->
 * 
 * s^2 * cdir[2]^2 + 2 * s * cdir[2] * ppos_v[2] + ppos_v[2]^2.              (5)
 *
 * In the last step to reach Eq. 5 we have assumed the less generic case in 
 * which the axis of the cone is allinged with the z axis of the reference frame 
 * (or equivalently, that cdir has been projected into the arbitrary axis 
 * direction of the cone). Additionally, "[2]" represents the z axis position 
 * (since indexes in C strat at 0, 2 is the third index. Now we can also expand 
 * the right hand side of Eq. 4:
 * ||s * cdir + ppos_v||^2 * mu^2 ->
 *
 * ((s * cdir + ppos_v) · (s * cdir + ppos_v)) * mu^2 ->
 *
 * (s^2 * (cdir · cdir) + 2 * s * (cdir · ppos_v) + (ppos_v · ppos_v)^2) * 
 *   mu^2 ->
 * 
 * (s^2 + 2 * s * (cdir · ppos_v) + (ppos_v · ppos_v)^2) * mu^2 -> 
 * 
 * s^2 * mu^2 + 2 * s * (cdir · ppos_v) * mu^2 + (ppos_v · ppos_v)^2 * mu^2, (6)
 *
 * where since cdir is a unit vector, cdir · cdir = 1. We can rearrange the 
 * expansions for the righ and left side of Eq. 4 by combining equation 5 and 6:
 * s^2 * cdir[2]^2 + 2 * s * cdir[2] * ppos_v[2] + ppos_v[2]^2 - 
 *   (s^2 * mu^2 + 2 * s * (cdir · ppos_v) * mu^2 + 
 *   (ppos_v · ppos_v)^2 * mu^2) = 0 ->
 * 
 * s^2 * cdir[2]^2 - s^2 * mu^2 + 2 * s * cdir[2] * ppos_v[2] - 
 *   2 * s * (cdir · ppos_v) * mu^2 + ppos_v[2]^2 - 
 *   (ppos_v · ppos_v)^2 * mu^2 ->
 *
 * s^2 * (cdir[2]^2 - mu^2) + 2 * s * (cdir[2] * ppos_v[2] - (cdir · ppos_v) * 
 *   mu^2) + ppos_v[2]^2 - (ppos_v · ppos_v)^2 * mu^2,                       (7)
 * 
 * which is a quadratic equation of the form A * x^2 + B * x + C = 0, with:
 * A = cdir[2]^2 - mu^2,
 * B = 2 * cdir[2] * ppos_v[2] - (cdir · ppos_v) * mu^2,
 * C = ppos_v[2]^2 - (ppos_v · ppos_v)^2 * mu^2
 * delta = B^2 - 4 * A * C
 * s = - B +- sqrt(delta) / (2 * A)
 *
 * if A == 0, the inclination of the ray is the same as of the cone's surface 
 * and the solution is linear with:
 * s = - C / (2 * B)
 *
 * and if A == 0 and B == 0, the line segment passes through the vertex.
 *
 * if A != 0, and if delta < 0, no intersection occours. If delta = 0, the line 
 * segment intersect the cone at one point and if delta > 0, the line segment 
 * intersect the cone in two points (or many points if it lies on the surface of 
 * the cone). For a finite and optionally truncated cone, it is necessary to 
 * test it the point of intersecton is in the range of the cone height relative 
 * to vertex.
 * 
 * If it is a hollow cone, those tests are sufficient. If it is a closed or
 * semi-closed cone, it is also necessary to check for intersection at the 
 * circular caps.
 *
 * It is important to remember that this function assumes that the cone's axis 
 * is alligned with the Z axis, therefore is necessary to project the ray cosine 
 * direction into the reference frame of the cone.
 * 
 */
/*
 int intrs_cone (const double *ppos, 
                 const double *cpos, 
                 const double *cdir, 
                 const double s, 
                 const struct str_cone *cone,
                 const int    segment)
 {
   double ppos_v[3];
   double cdotpv, pvdotpv, cdota, pvdota;
   double delta, sdelta;
   double s_intrs;
   double p_intrs[3];
   double A, B, C; 
   double p_cdir;

   ppos_v[0] = ppos[0] - (*cone).vertex[0];
   ppos_v[1] = ppos[1] - (*cone).vertex[1];
   ppos_v[2] = ppos[2] - (*cone).vertex[2];

   cdotpv  = cdir[0] * ppos_v[0] + 
             cdir[1] * ppos_v[1] + 
             cdir[2] * ppos_v[2];

   pvdotpv = ppos_v[0] * ppos_v[0] + 
             ppos_v[1] * ppos_v[1] + 
             ppos_v[2] * ppos_v[2];

/*   If ray cosine directions are not projected, then the equation would be:

   cdota   = cdir[0] * (*cone).caxis[0] + 
             cdir[1] * (*cone).caxis[1] + 
             cdir[2] * (*cone).caxis[2];

   pvdota  = ppos_v[0] * (*cone).caxis[0] + 
             ppos_v[1] * (*cone).caxis[1] + 
             ppos_v[2] * (*cone).caxis[2];

   A = cdota * cdota - (*cone).mu_sq;
   B = 2.0 * (cdota * pvdota - cdotpv * (*cone).mu_sq);
   C = pvdota * pvdota - pvdotpv * (*cone).mu_sq;
*/
/*
   A = cdir[2] * cdir[2] - (*cone).mu_sq;
   B = 2.0 * (cdir[2] * ppos_v[2] - cdotpv * (*cone).mu_sq);
   C = ppos_v[2] * ppos_v[2] - pvdotpv * (*cone).mu_sq;

   if ( NUM_EQU(A, 0.0, 1E-12) )
   {
     // If A == 0, quadratic term is zero and root is linear: B * x + C = 0
     // -> x = -C / B.
     if ( NUM_EQU(B, 0.0, 1E-12) )
     {
       // If A and B are zero, line would pass through the vertex. Line segment
       // migt not, but still intersection occours as path is through the 
       // surface of the cone. All is needed is that the z ranges intersect.
       if ( NUM_EQU((*cone).height[0], 0.0, 1E-12) )
       { 
         if (cdir[2] > 0)
         {
           if ( ( (ppos[2] > ((*cone).vertex[2] + (*cone).height[0])) && 
                  (ppos[2] < ((*cone).vertex[2] + (*cone).height[1])) ) || 
                ( (ppos[2] < ((*cone).vertex[2] + (*cone).height[0])) && 
                  (cpos[2] > ((*cone).vertex[2] + (*cone).height[0])) ) )
           {
             return 1;
           }
         } else {
           if ( ( (ppos[2] < ((*cone).vertex[2] + (*cone).height[1])) && 
                  (ppos[2] > ((*cone).vertex[2] + (*cone).height[0])) ) || 
                ( (ppos[2] > ((*cone).vertex[2] + (*cone).height[1])) && 
                  (cpos[2] < ((*cone).vertex[2] + (*cone).height[1])) ) )
           {
             return 1;
           }
         }
       }
     } else {
       // If B is not zero, the root (= -C / B) gives the single point of 
       // intersection. If this root is positive and lower than the panthlength 
       // s of the ray, intersection is possible. However, we whant only the 
       // positive part of the cone, so we can project the vertex centered 
       // vector of the intersection point into the cone axis and this 
       // projection has to be positive.  
       s_intrs = -C / B;
       if ( (s_intrs > 0.0) && (s_intrs <= s) )
       {
         p_intrs[0] = ppos_v[0] + cdir[0] * s_intrs; 
         p_intrs[1] = ppos_v[1] + cdir[1] * s_intrs;
         p_intrs[2] = ppos_v[2] + cdir[2] * s_intrs;
         p_cdir = p_intrs[0] * (*cone).caxis[0] + 
                  p_intrs[1] * (*cone).caxis[1] +
                  p_intrs[2] * (*cone).caxis[2];
         if ( p_cdir > 0.0 && (p_intrs[2] > (*cone).height[0]) && 
           (p_intrs[2] < (*cone).height[1]))
         {
           return 1;
         }
       }
     }    
   } else {
     delta = B * B - 4 * A * C;
     if (delta >= 0.0) 
     {
       sdelta = sqrt(delta);

       s_intrs = ((-B - sdelta) / (2 * A));
       if ( (s_intrs > 0.0) && (s_intrs <= s) )
       {
         p_intrs[0] = ppos_v[0] + cdir[0] * s_intrs; 
         p_intrs[1] = ppos_v[1] + cdir[1] * s_intrs;
         p_intrs[2] = ppos_v[2] + cdir[2] * s_intrs;
         p_cdir = p_intrs[0] * (*cone).caxis[0] + 
                  p_intrs[1] * (*cone).caxis[1] +
                  p_intrs[2] * (*cone).caxis[2];
         if ( (p_cdir > 0.0) && (p_intrs[2] > (*cone).height[0]) && 
           (p_intrs[2] < (*cone).height[1]) )
         {
           return 1;
         }
       }

       s_intrs = ((-B + sdelta) / (2 * A));
       if ( (s_intrs > 0.0) && (s_intrs <= s) )
       {
         p_intrs[0] = ppos_v[0] + cdir[0] * s_intrs; 
         p_intrs[1] = ppos_v[1] + cdir[1] * s_intrs;
         p_intrs[2] = ppos_v[2] + cdir[2] * s_intrs;
         p_cdir = p_intrs[0] * (*cone).caxis[0] + 
                  p_intrs[1] * (*cone).caxis[1] +
                  p_intrs[2] * (*cone).caxis[2];
         if ( (p_cdir > 0.0) && (p_intrs[2] > (*cone).height[0]) && 
           (p_intrs[2] < (*cone).height[1]) )
         {
           return 1;
         }
       }
     }
   }

   // If no intersections happened at the side surface, it might still be 
   // necessary to test against circular caps if semi-closed or closed:
   if ( (*cone).closed[0] )
   {
     s_intrs = ((*cone).height[0] - ppos_v[2]) / cdir[2];
     if (s_intrs > 0 && s_intrs < s)
     {
       p_intrs[0] = ppos_v[0] + cdir[0] * s_intrs; 
       p_intrs[1] = ppos_v[1] + cdir[1] * s_intrs;
       if ( sqrt(p_intrs[0] * p_intrs[0] + p_intrs[1] * p_intrs[1]) <= 
            (*cone).radius[0] )
       {
         return 1;
       }
     }
   }

   if ( (*cone).closed[1] )
   {
     s_intrs = ((*cone).height[1] - ppos_v[2]) / cdir[2];
     if (s_intrs > 0 && s_intrs < s)
     {
       p_intrs[0] = ppos_v[0] + cdir[0] * s_intrs; 
       p_intrs[1] = ppos_v[1] + cdir[1] * s_intrs;
       if ( sqrt(p_intrs[0] * p_intrs[0] + p_intrs[1] * p_intrs[1]) <= 
            (*cone).radius[1] )
       {
         return 1;
       }
     }
   }

   return 0;
 }

 void setup_cone (struct str_cone *cone,
                  const double *cn_vertex,
                  const double *cn_axis,
                  const double *cn_height,
                  const double cn_theta,
                  const int    *cn_closed)
 {
   (*cone).vertex[0] = cn_vertex[0];
   (*cone).vertex[1] = cn_vertex[1];
   (*cone).vertex[2] = cn_vertex[2];

   (*cone).paxis[0] = cn_axis[0];
   (*cone).paxis[1] = cn_axis[1];

   (*cone).caxis[0] = sin(cn_axis[0]) * cos(cn_axis[1]);
   (*cone).caxis[1] = sin(cn_axis[0]) * sin(cn_axis[1]);
   (*cone).caxis[2] = cos(cn_axis[0]);

   (*cone).psi   = cn_theta;
   (*cone).mu    = cos(cn_theta);
   (*cone).mu_sq = (*cone).mu * (*cone).mu;

   (*cone).height[0] = cn_height[0];
   (*cone).height[1] = cn_height[1];

   (*cone).radius[0] = (*cone).height[0] * tan(cn_theta);
   (*cone).radius[1] = (*cone).height[1] * tan(cn_theta);

   (*cone).closed[0] = cn_closed[0];
   (*cone).closed[1] = cn_closed[1];

   (*cone).rotate_f = ((*cone).caxis[2] < 1.0)? 1 : 0;
 }

*/

