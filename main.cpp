/**
     Shadow equation is given by:
                  [x^2 + y^2 + z^2] cos(epsilon)^2 == x^2 - 2 x R sin(epsilon) + R^2
     Transforms into quartic equation by parametric substitution:
                  coord_i = c_i + a_i cos(fi) + b_i sin(fi)  (fi -> 0 to 2pi)
     The trigonometric substitution
                  cos(fi) = (1-t^2)/(1+t^2), sin = 2t/(1+t^2)
     on both sides of the equation leads to:

                  a1 t^4 + b1 t^3 + c1 t^2 + d1 t + e1 = a2 t^4 + b2 t^3 + c2 t^2 + d2 t + e2
     Bringing the right and the left terms together we get quartic equation with t = tg(fi/2):
                  A t^4 + B t^3 + C t^2 + D t + E = 0
**/

#include <iostream>
#include "calc_shadow_intersections.h"


/**test example**/
int main() {

    //vec7 = {rx, ry, rz, vx, vy, vz, time}
    auto orbit = orbitFromVec7(
            {2.66306016e+07, 9.70516730e+06, 0.00000000e+00, -2.34551642e+03, 1.33025620e+03, 0.00000000e+00, 0});

    //epsilon is equal earthR/earth2sunR and depends on time
    calc_shadow_intersection_times(orbit, 4.2666e-6);

    return 0;
}