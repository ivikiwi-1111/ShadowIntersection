//
// Created by is on 18.12.22.
//

#ifndef SHADOW_INTERSECTION_MODIFIED_CALC_SHADOW_INTERSECTIONS_H
#define SHADOW_INTERSECTION_MODIFIED_CALC_SHADOW_INTERSECTIONS_H

#endif //SHADOW_INTERSECTION_MODIFIED_CALC_SHADOW_INTERSECTIONS_H

#include <vector>
#include <algorithm>
#include "quartic.h"
#include "orbit_params.h"

using std::sin;
using std::cos;

scalar sq(const scalar &a) {
    return a * a;
}


std::pair<std::array<scalar, 4>, unsigned short> calc_shadow_intersection_anomaly_sorted(Orbit const &orb, scalar const epsilon) {
    scalar const sin_eps = sin(epsilon);
    scalar const cos_eps = cos(epsilon);
    vec3 const c_ = -orb.a * orb.e;

    //Substitution params coord_i = c_i + a_i cos(fi) + b_i sin(fi)  (fi -> 0 to 2pi)
    scalar const a1 = orb.a.x, a2 = orb.a.y, a3 = orb.a.z;
    scalar const b1 = orb.b.x, b2 = orb.b.y, b3 = orb.b.z;
    scalar const c1 = c_.x, c2 = c_.y, c3 = c_.z;

    //Quartic eq. A t^4 + B t^3 + C t^2 + D t + E = 0
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
            )/4.;


    //Quartic solver call -> returns complex solutions of quartic eq:
    auto const solution = quartic_solver(B / A, C / A, D / A, E / A);


    unsigned short counter = 0;
    std::array<scalar, 4> anomaly{};
    //Get the intersection anomalies
    for(auto it:solution){
        //abs(exp) = 1 <=> fi is real
        if( fabs(fabs(it) - 1) < std::numeric_limits<scalar>::epsilon()){
            scalar E0 = atan2(it.imag(), it.real());
            //x>0 <=> point is in shadow
            if(a1*cos(E0)+b1*sin(E0)+c1 > 0){
                anomaly[counter] = E0>orb.E0 ? E0 : E0+2*PI;
                counter ++;
            }
        }
    }
    std::sort(anomaly.begin(), anomaly.begin()+counter);
    return std::make_pair(anomaly, counter);
}

//already sorted
std::pair<std::array<scalar, 4>, unsigned short> get_intersection_times(Orbit const &orb, scalar const epsilon){
    auto const v = calc_shadow_intersection_anomaly_sorted(orb, epsilon);
    std::array<scalar, 4> times{};

    for(auto i = 0; i<v.second; ++i){
        scalar const ndt = (v.first[i] - orb.E0) - orb.e * (sin(v.first[i]) - sin(orb.E0));
        scalar const n = sqrt(pow(dot(orb.a, orb.a), -1.5f) * earthMu);
        times[i] = orb.t0 + ndt / n;
    }
    return make_pair(times, v.second);
}
