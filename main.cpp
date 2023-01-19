#include <iostream>
#include "calc_shadow_intersections.h"


/**test example**/
int main() {

    //vec7 = {rx, ry, rz, vx, vy, vz, time}
    auto orbit = orbitFromVec7(
            {2.66306016e+07, 9.70516730e+06, 0.00000000e+00, -2.34551642e+03, 1.33025620e+03, 0.00000000e+00, 0});

    //epsilon is equal earthR/earth2sunR and depends on time
    auto v = calc_shadow_intersection_times(orbit, 4.2666e-6);

    std::cout << "Times till intersection:" << std::endl;
    for( auto it: v){
        if(it>0)
            std::cout << it << " ";
    }
    return 0;
}