//
// Created by is on 15.12.22.
//

#ifndef SHADOW_INTERSECTION_MODIFIED_STRUCTS_H
#define SHADOW_INTERSECTION_MODIFIED_STRUCTS_H

#endif //SHADOW_INTERSECTION_MODIFIED_STRUCTS_H

#include "vec3.h"

using std::sqrt;
using std::atan2;


using Time = scalar;
struct vec7 {
    vec3 r;
    vec3 v;
    Time t;
};

struct Orbit {
    vec3 a, b;
    scalar e;
    scalar E0;
    Time t0;
};


constexpr scalar earthMu = 3.986004418e14;
constexpr scalar PI = 3.1415926535;
constexpr scalar earthR = 6371e3;

Orbit orbitFromVec7(vec7 const point) {
    vec3 c = cross(point.r, point.v);

    auto const r = normalize(point.r);
    auto const e = cross(point.v, c) / earthMu - r;
    auto const e2 = dot(e, e);
    vec3 const ex = e2 > 1e-7f ? normalize(e) : r;
    vec3 const ey = cross(normalize(c), ex);

    auto const p = dot(c, c) / earthMu;
    auto const a = p / (1.f - e2);
    scalar const b = a * sqrt(1.f - e2);

    scalar const sinE0 = dot(point.r, ey) / b;
    scalar const cosE0 = dot(point.r, ex) / a + sqrt(e2);

    return Orbit
            {
                    ex * a,
                    ey * b,
                    sqrt(e2),
                    atan2(sinE0, cosE0),
                    point.t
            };
}
