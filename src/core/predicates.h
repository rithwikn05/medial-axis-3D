#pragma once

#include "vec3.h"

#ifdef __cplusplus
extern "C" {
#endif

void exactinit();

double orient3d(const double* pa,
                const double* pb,
                const double* pc,
                const double* pd);

double insphere(const double* pa,
                const double* pb,
                const double* pc,
                const double* pd,
                const double* pe);

#ifdef __cplusplus
}
#endif

namespace medial_axis_3d {

inline void ensure_predicates_initialized() {
    static bool initialized = false;
    if (!initialized) {
        exactinit();
        initialized = true;
    }
}

inline double orient3d(const Vec3& a,
                       const Vec3& b,
                       const Vec3& c,
                       const Vec3& d) {
    ensure_predicates_initialized();

    const double pa[3] = {a.x, a.y, a.z};
    const double pb[3] = {b.x, b.y, b.z};
    const double pc[3] = {c.x, c.y, c.z};
    const double pd[3] = {d.x, d.y, d.z};

    return ::orient3d(pa, pb, pc, pd);
}

inline double insphere(const Vec3& a,
                       const Vec3& b,
                       const Vec3& c,
                       const Vec3& d,
                       const Vec3& e) {
    ensure_predicates_initialized();

    const double pa[3] = {a.x, a.y, a.z};
    const double pb[3] = {b.x, b.y, b.z};
    const double pc[3] = {c.x, c.y, c.z};
    const double pd[3] = {d.x, d.y, d.z};
    const double pe[3] = {e.x, e.y, e.z};

    return ::insphere(pa, pb, pc, pd, pe);
}

inline bool insphere_contains(const Vec3& a,
                              const Vec3& b,
                              const Vec3& c,
                              const Vec3& d,
                              const Vec3& point) {
    const double orientation = orient3d(a, b, c, d);
    if (orientation > 0.0) {
        return insphere(a, b, c, d, point) > 0.0;
    } else if (orientation < 0.0) {
        return insphere(a, b, c, d, point) < 0.0;
    }
    return false;
}

}  // namespace medial_axis_3d
