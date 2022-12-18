//
// Created by is on 18.12.22.
//

#ifndef SHADOW_INTERSECTION_MODIFIED_CALC_SHADOW_INTERSECTIONS_H
#define SHADOW_INTERSECTION_MODIFIED_CALC_SHADOW_INTERSECTIONS_H

#endif //SHADOW_INTERSECTION_MODIFIED_CALC_SHADOW_INTERSECTIONS_H

#include <vector>
#include "quartic.h"
#include "orbit_params.h"

using complex = std::complex<scalar>;
using std::sin;
using std::cos;

scalar sq(const scalar &a) {
    return a * a;
}

std::vector<scalar> calc_shadow_intersection_times(Orbit const &orb, scalar const epsilon) {
    auto const sin_eps = sin(epsilon);
    auto const cos_eps = cos(epsilon);
    auto c_ = -orb.a * orb.e;


    auto const a0 = orb.a.x, a1 = orb.a.y, a2 = orb.a.z;
    auto const b0 = orb.b.x, b1 = orb.b.y, b2 = orb.b.z;
    auto const c0 = c_.x, c1 = c_.y, c2 = c_.z;

    auto a_1 = sq(c0 - a0) + sq(c1 - a1) + sq(c2 - a2);
    a_1 = a_1 * sq(cos_eps);

    scalar const a_2 = sq(c0 - a0) - 2. * earthR * sin_eps * (c0 - a0) + sq(earthR);


    auto b_1 = 4. * (b0 * (c0 - a0) + b1 * (c1 - a1) + b2 * (c2 - a2));
    b_1 = b_1 * sq(cos_eps);

    scalar const b_2 = 4. * b0 * (c0 - a0) - 4. * b0 * earthR * sin_eps;


    scalar c_1 =
            (2. * sq(c0) + 4. * sq(b0) - 2. * sq(a0)) +
            (2. * sq(c1) + 4. * sq(b1) - 2. * sq(a1)) +
            (2. * sq(c2) + 4. * sq(b2) - 2. * sq(a2));

    c_1 = c_1 * sq(cos_eps);
    scalar const c_2 = (2. * sq(c0) + 4. * sq(b0) - 2. * sq(a0)) - 4. * c0 * earthR * sin_eps + 2. * sq(earthR);


    scalar d_1 = 4. * b0 * (a0 + c0) + 4. * b1 * (a1 + c1) + 4. * b2 * (a2 + c2);
    d_1 = d_1 * sq(cos_eps);
    scalar const d_2 = 4. * b0 * (a0 + c0) - 4. * b0 * earthR * sin_eps;


    auto e_1 = sq(c0 + a0) + sq(c1 + a1) + sq(c2 + a2);
    e_1 = e_1 * sq(cos_eps);
    scalar const e_2 = sq(c0 + a0) - 2. * (c0 + a0) * earthR * sin_eps + sq(earthR);


    scalar const A = a_1 - a_2, B = b_1 - b_2, C = c_1 - c_2, D = d_1 - d_2, E = e_1 - e_2;

    auto v = quartic_solver(B / A, C / A, D / A, E / A);

    std::vector<scalar> res;
    //inf case not solved yet
    for (auto &it: v) {
        //roots are real
        if (std::abs(it.imag()) < 1e-12) {
            scalar const t = it.real();
            scalar const sin_ = 2 * t / (1 + t * t);
            scalar const cos_ = (1 - t * t) / (1 + t * t);

            //point is in shadow
            if (c0 + a0 * cos_ + b0 * sin_ >= 0) {
                scalar const E1 = atan2(sin_, cos_);
                scalar const ndt = (E1 - orb.E0) - orb.e * (std::sin(E1) - std::sin(orb.E0));
                scalar const n = sqrt(pow(dot(orb.a, orb.a), -1.5f) * earthMu);
                res.insert(res.end(), ndt > 0 ? orb.t0 + ndt / n : (orb.t0 + ndt + 2 * PI) / n);
            }
        }
    }
    return res;
}