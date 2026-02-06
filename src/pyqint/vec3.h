#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

typedef double mat33[3][3];

/**
 * @brief Custom 3-vector class.
 */
class Vec3 {
public:
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    /**
     * @brief Default constructor.
     */
    Vec3() {}

    /**
     * @brief Initialization constructor.
     *
     * @param _x  x component
     * @param _y  y component
     * @param _z  z component
     */
    Vec3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

    /**
     * @brief Access vector components by index.
     *
     * @param n  component index (0=x, 1=y, 2=z)
     *
     * @return reference to the selected component
     */
    inline double& operator[](int n) {
        switch(n) {
            case 0:
                return this->x;
            case 1:
                return this->y;
            case 2:
                return this->z;
            default:    // should never reach this
                return this->z;
        }
    }

    /**
     * @brief Access vector components by index (const).
     *
     * @param n  component index (0=x, 1=y, 2=z)
     *
     * @return const reference to the selected component
     */
    inline const double& operator[](int n) const {
        switch(n) {
            case 0:
                return this->x;
            case 1:
                return this->y;
            case 2:
                return this->z;
            default:    // should never reach this
                return this->z;
        }
    }

    /**
     * @brief Multiply a matrix with a vector.
     *
     * @param lhs  3x3 matrix
     * @param rhs  vector
     *
     * @return result of matrix-vector multiplication
     */
    friend Vec3 operator*(const mat33& lhs, const Vec3& rhs) {
        Vec3 r;
        r.x = lhs[0][0] * rhs.x + lhs[1][0] * rhs.y + lhs[2][0] * rhs.z;
        r.y = lhs[0][1] * rhs.x + lhs[1][1] * rhs.y + lhs[2][1] * rhs.z;
        r.z = lhs[0][2] * rhs.x + lhs[1][2] * rhs.y + lhs[2][2] * rhs.z;

        return r;
    }

    /**
     * @brief Multiply a scalar with a vector.
     *
     * @param v    scalar
     * @param rhs  vector
     *
     * @return scaled vector
     */
    friend Vec3 operator*(double v, const Vec3& rhs) {
        return Vec3(v * rhs.x, v * rhs.y, v * rhs.z);
    }

    /**
     * @brief Multiply a vector with a scalar.
     *
     * @param rhs  vector
     * @param v    scalar
     *
     * @return scaled vector
     */
    friend Vec3 operator*(const Vec3& rhs, double v) {
        return Vec3(v * rhs.x, v * rhs.y, v * rhs.z);
    }

    /**
     * @brief Return normalized vector.
     *
     * @return normalized vector
     */
    Vec3 normalized() const {
        double l = std::sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
        return Vec3(this->x / l, this->y / l, this->z / l);
    }

    /**
     * @brief Return dot product between two vectors.
     *
     * @param rhs  right-hand-side vector
     *
     * @return dot product
     */
    double dot(const Vec3& rhs) const {
        return this->x * rhs.x + this->y * rhs.y + this->z * rhs.z;
    }

    /**
     * @brief Addition operation between two vectors.
     *
     * @param lhs  left-hand-side vector
     * @param rhs  right-hand-side vector
     *
     * @return summed vector
     */
    friend Vec3 operator+(const Vec3& lhs, const Vec3& rhs) {
        return Vec3(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
    }

    /**
     * @brief Subtraction operation between two vectors.
     *
     * @param lhs  left-hand-side vector
     * @param rhs  right-hand-side vector
     *
     * @return difference vector
     */
    friend Vec3 operator-(const Vec3& lhs, const Vec3& rhs) {
        return Vec3(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
    }

    /**
     * @brief Addition assignment operation.
     *
     * @param rhs  right-hand-side vector
     */
    void operator+=(const Vec3& rhs) {
        this->x += rhs.x;
        this->y += rhs.y;
        this->z += rhs.z;
    }

    /**
     * @brief Subtraction assignment operation.
     *
     * @param rhs  right-hand-side vector
     */
    void operator-=(const Vec3& rhs) {
        this->x -= rhs.x;
        this->y -= rhs.y;
        this->z -= rhs.z;
    }

    /**
     * @brief Divide vector by a scalar operation.
     *
     * @param rhs  vector
     * @param v    scalar
     *
     * @return scaled vector
     */
    friend Vec3 operator/(const Vec3& rhs, double v) {
        return Vec3(rhs.x / v, rhs.y / v, rhs.z / v);
    }

    /**
     * @brief Calculate cross product between two vectors.
     *
     * @param rhs  right-hand-side vector
     *
     * @return cross product
     */
    Vec3 cross(const Vec3& rhs) const {
        return Vec3(
            this->y * rhs.z - this->z * rhs.y,
            this->z * rhs.x - this->x * rhs.z,
            this->x * rhs.y - this->y * rhs.x
        );
    }

    /**
     * @brief Calculate squared sum of coefficients.
     *
     * @return squared norm
     */
    double norm2() const {
        return (this->x * this->x) + (this->y * this->y) + (this->z * this->z);
    }

    /**
     * @brief Calculate product of coefficients.
     *
     * @return product of x, y, and z components
     */
    double prod() const {
        return this->x * this->y * this->z;
    }
};
