#pragma once

#include <cmath>
#include <ostream>

namespace medial_axis_3d {

struct Vec3 {
    double x{};
    double y{};
    double z{};

    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

inline Vec3 operator+(const Vec3& a, const Vec3& b) {
    return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Vec3 operator-(const Vec3& a, const Vec3& b) {
    return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline Vec3 operator*(const Vec3& a, double scalar) {
    return Vec3(a.x * scalar, a.y * scalar, a.z * scalar);
}

inline Vec3 operator*(double scalar, const Vec3& a) {
    return a * scalar;
}

inline Vec3 operator/(const Vec3& a, double scalar) {
    return Vec3(a.x / scalar, a.y / scalar, a.z / scalar);
}

inline bool operator==(const Vec3& a, const Vec3& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

inline double squared_norm(const Vec3& a) {
    return dot(a, a);
}

inline double norm(const Vec3& a) {
    return std::sqrt(squared_norm(a));
}

inline Vec3 normalized(const Vec3& a) {
    const double length = norm(a);
    return length > 0.0 ? a / length : Vec3{};
}

inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    os << '(' << v.x << ", " << v.y << ", " << v.z << ')';
    return os;
}

}  // namespace medial_axis_3d
