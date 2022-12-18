//
// Created by is on 15.12.22.
//

#ifndef SHADOW_INTERSECTION_MODIFIED_STRUCTS_H
#define SHADOW_INTERSECTION_MODIFIED_STRUCTS_H

#endif //SHADOW_INTERSECTION_MODIFIED_STRUCTS_H

using scalar = double;

struct vec3 {
    scalar x, y, z;
};

inline vec3 operator+(vec3 const &v1, vec3 const &v2) {
    return
            {
                    v1.x + v2.x,
                    v1.y + v2.y,
                    v1.z + v2.z,
            };
}

inline vec3 operator-(vec3 const &v1, vec3 const &v2) {
    return
            {
                    v1.x - v2.x,
                    v1.y - v2.y,
                    v1.z - v2.z,
            };
}

inline vec3 operator-(vec3 const &v) {
    return {-v.x, -v.y, -v.z};
}

inline vec3 operator*(vec3 const &v, scalar const f) {
    return
            {
                    v.x * f,
                    v.y * f,
                    v.z * f,
            };
}


inline vec3 operator/(vec3 const &v, scalar const f) {
    return
            {
                    v.x / f,
                    v.y / f,
                    v.z / f,
            };
}

inline vec3 operator/(vec3 const &v1, vec3 const &v2) {
    return
            {
                    v1.x / v2.x,
                    v1.y / v2.y,
                    v1.z / v2.z,
            };
}

inline scalar dot(vec3 const &v1, vec3 const &v2) {
    return v1.x * v2.x
           + v1.y * v2.y
           + v1.z * v2.z;
}

inline vec3 cross(vec3 const &v1, vec3 const &v2) {
    return
            {
                    v1.y * v2.z - v1.z * v2.y,
                    v1.z * v2.x - v1.x * v2.z,
                    v1.x * v2.y - v1.y * v2.x,
            };
}

inline scalar length(vec3 const &v) {
    return std::sqrt(dot(v, v));
}

inline vec3 normalize(vec3 const &v) {
    return v / length(v);
}