//
// Created by is on 15.12.22.
//

#ifndef SHADOW_INTERSECTION_MODIFIED_QUARTIC_H
#define SHADOW_INTERSECTION_MODIFIED_QUARTIC_H

#endif //SHADOW_INTERSECTION_MODIFIED_QUARTIC_H

#include <cmath>
#include <complex>
#include <array>
#include <vector>

using scalar = double;
using complex = std::complex<scalar>;


complex complex_cbrt(complex const &z) {
    return std::pow(z, 1. / 3.);
}


std::array<complex, 4> quartic_solver(scalar const &b, scalar const &c, scalar const &d, scalar const &e) {

    complex const Q1 = c * c - 3. * b * d + 12. * e;
    complex const Q2 = 2. * c * c * c - 9. * b * c * d + 27. * d * d + 27. * b * b * e - 72. * c * e;
    complex const Q3 = 8. * b * c - 16. * d - 2. * b * b * b;
    complex const Q4 = 3. * b * b - 8. * c;

    complex const Q5 = complex_cbrt(Q2 / 2. + std::sqrt(Q2 * Q2 / 4. - Q1 * Q1 * Q1));
    complex const Q6 = (Q1 / Q5 + Q5) / 3.;
    complex const Q7 = 2. * std::sqrt(Q4 / 12. + Q6);


    return std::array<complex, 4>{
            (-b - Q7 - std::sqrt(4. * Q4 / 6. - 4. * Q6 - Q3 / Q7)) / 4.,
            (-b - Q7 + std::sqrt(4. * Q4 / 6. - 4. * Q6 - Q3 / Q7)) / 4.,
            (-b + Q7 - std::sqrt(4. * Q4 / 6. - 4. * Q6 + Q3 / Q7)) / 4.,
            (-b + Q7 + std::sqrt(4. * Q4 / 6. - 4. * Q6 + Q3 / Q7)) / 4.
    };

}