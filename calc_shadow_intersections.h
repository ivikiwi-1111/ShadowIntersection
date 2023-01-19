//
// Created by is on 18.12.22.
//

#ifndef SHADOW_INTERSECTION_MODIFIED_CALC_SHADOW_INTERSECTIONS_H
#define SHADOW_INTERSECTION_MODIFIED_CALC_SHADOW_INTERSECTIONS_H

#endif //SHADOW_INTERSECTION_MODIFIED_CALC_SHADOW_INTERSECTIONS_H

#include <vector>
#include "quartic.h"
#include "orbit_params.h"

using std::sin;
using std::cos;

scalar sq(const scalar &a) {
    return a * a;
}

void calc_shadow_intersection_times(Orbit const &orb, scalar const epsilon) {
    auto const sin_eps = sin(epsilon);
    auto const cos_eps = cos(epsilon);
    auto c_ = -orb.a * orb.e;


    auto const a1 = orb.a.x, a2 = orb.a.y, a3 = orb.a.z;
    auto const b1 = orb.b.x, b2 = orb.b.y, b3 = orb.b.z;
    auto const c1 = c_.x, c2 = c_.y, c3 = c_.z;

    complex const A = (
            (sq(a1) + sq(a2) + sq(a3))*sq(cos_eps) -
            sq(a1)  + sq(b1) -
            (sq(b1) + sq(b2) + sq(b3))*sq(cos_eps) -
            (a1*b1 + a2*b2 + a3*b3)*sq(cos_eps)*2i + a1*b1*2i
            )/4.;

    complex const B = (
            (a1*c1 + a2*c2 + a3*c3)*sq(cos_eps)    -
            (b1*c1 + b2*c2 + b3*c3)*sq(cos_eps)*1i +
            (earthR*a1*sin_eps - a1*c1) -
            (earthR*b1*sin_eps - b1*c1)*1i
            );

    complex const C = (
            (sq(c1) + sq(c2) + sq(c3))*sq(cos_eps)*2 +
            (sq(a1) + sq(a2) + sq(a3))*sq(cos_eps)   +
            (sq(b1) + sq(b2) + sq(b3))*sq(cos_eps)   -
            2*sq(earthR) + 4*earthR*c1*sin_eps -
            (sq(a1) + sq(b1) + 2*sq(c1))
            )/2.;

    complex const D = (
            (b1*c1 + b2*c2 + b3*c3)*sq(cos_eps)*1i +
            (a1*c1 + a2*c2 + a3*c3)*sq(cos_eps)    +
            (earthR*b1*sin_eps - b1*c1)*1i            +
            (earthR*a1*sin_eps - a1*c1)
            );


    complex const E =(
            (sq(a1) + sq(a2) + sq(a3))*sq(cos_eps)   -
            sq(a1)  + sq(b1) -
            (sq(b1) + sq(b2) + sq(b3))*sq(cos_eps)   +
            (a1*b1 + a2*b2 + a3*b3)*sq(cos_eps)*2i - a1*b1*2i
            )/4.;;


    auto v = quartic_solver(B / A, C / A, D / A, E / A);
    for(auto it:v){
        std::cout << abs(it) << " ";
    }

}