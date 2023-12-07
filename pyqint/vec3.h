#ifndef _VEC3_H
#define _VEC3_H

typedef double mat33[3][3];

/**
 * Custom 3-vector class
*/
class Vec3 {
public:
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    /**
     * Default constructor 
    */
    Vec3() {}

    /**
     * Initialization constructor
    */
    Vec3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

    double& operator[](int n) {
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

    const double& operator[](int n) const {
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
     * Multiplication operation between matrix and vector
    */
    friend Vec3 operator*(const mat33& lhs, const Vec3& rhs) {
        Vec3 r;
        r.x = lhs[0][0] * rhs.x + lhs[1][0] * rhs.y + lhs[2][0] * rhs.z;
        r.y = lhs[0][1] * rhs.x + lhs[1][1] * rhs.y + lhs[2][1] * rhs.z;
        r.z = lhs[0][2] * rhs.x + lhs[1][2] * rhs.y + lhs[2][2] * rhs.z;

        return r;
    }

    /**
     * Multiplication operation between scalar and vector
    */
    friend Vec3 operator*(double v, const Vec3& rhs) {
        return Vec3(v * rhs.x, v * rhs.y, v * rhs.z);
    }

    /**
     * Multiplication operation between vector and scalar
    */
    friend Vec3 operator*(const Vec3& rhs, double v) {
        return Vec3(v * rhs.x, v * rhs.y, v * rhs.z);
    }

    /**
     * Return normalized vector
    */
    Vec3 normalized() const {
        double l = std::sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
        return Vec3(this->x / l, this->y / l, this->z / l);
    }

    /**
     * Return dot product between two vectors
    */
    double dot(const Vec3& rhs) const {
        return this->x * rhs.x + this->y * rhs.y + this->z * rhs.z;
    }

    /**
     * Addition operation between two vectors
    */
    friend Vec3 operator+(const Vec3& lhs, const Vec3& rhs) {
        return Vec3(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
    }

    /**
     * Subtraction operation between two vectors
    */
    friend Vec3 operator-(const Vec3& lhs, const Vec3& rhs) {
        return Vec3(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
    }

    /**
     * Addition assignment operation
    */
    void operator+=(const Vec3& rhs) {
        this->x += rhs.x;
        this->y += rhs.y;
        this->z += rhs.z;
    }

    /**
     * Subtraction assignment operation
    */
    void operator-=(const Vec3& rhs) {
        this->x -= rhs.x;
        this->y -= rhs.y;
        this->z -= rhs.z;
    }

    /**
     * Divide vector by a scalar operation
    */
    friend Vec3 operator/(const Vec3& rhs, double v) {
        return Vec3(rhs.x / v, rhs.y / v, rhs.z / v);
    }

    /**
     * Calculate cross product between two vectors
    */
    Vec3 cross(const Vec3& rhs) const {
        return Vec3(
            this->y * rhs.z - this->z * rhs.y,
            this->z * rhs.x - this->x * rhs.z,
            this->x * rhs.y - this->y * rhs.x
        );
    }

    /**
     * Calculate squared sum of coefficients
    */
    double norm2() const {
        return (this->x * this->x) + (this->y * this->y) + (this->z * this->z);
    }

    /**
     * Calculate product of coefficients
    */
    double prod() const {
        return this->x * this->y * this->z;
    }
};

#endif // _VEC3_H
