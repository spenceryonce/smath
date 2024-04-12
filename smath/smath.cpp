 #include "smath.hpp"
 
 namespace smath 
 {
    // Vector2 + operator overload
    Vector2 operator+(const Vector2& v1, const Vector2& v2) {
        return {v1.x + v2.x, v1.y + v2.y};
    }

    // Vector2 - operator overload
    Vector2 operator-(const Vector2& v1, const Vector2& v2) {
        return {v1.x - v2.x, v1.y - v2.y};
    }

    // Vector2 * operator overload
    Vector2 operator*(const Vector2& v, double scalar) {
        return {v.x * scalar, v.y * scalar};
    }

    // Vector2 * operator overload (2 vectors element-wise multiplication)
    Vector2 operator*(const Vector2& v1, const Vector2& v2) {
        return {v1.x * v2.x, v1.y * v2.y};
    }

    // Vector2 / operator overload
    Vector2 operator/(const Vector2& v, double scalar) {
        return {v.x / scalar, v.y / scalar};
    }

    // Vector2 / operator overload (2 vectors element-wise division)
    Vector2 operator/(const Vector2& v1, const Vector2& v2) {
        return {v1.x / v2.x, v1.y / v2.y};
    }


    // Vector3 + operator overload
    Vector3 operator+(const Vector3& v1, const Vector3& v2) {
        return {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
    }

    // Vector3 - operator overload
    Vector3 operator-(const Vector3& v1, const Vector3& v2) {
        return {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
    }

    // Vector3 * operator overload
    Vector3 operator*(const Vector3& v, double scalar) {
        return {v.x * scalar, v.y * scalar, v.z * scalar};
    }

    // Vector3 * operator overload (2 vectors element-wise multiplication)
    Vector3 operator*(const Vector3& v1, const Vector3& v2) {
        return {v1.x * v2.x, v1.y * v2.y, v1.z * v2.z};
    }

    // Vector3 / operator overload
    Vector3 operator/(const Vector3& v, double scalar) {
        return {v.x / scalar, v.y / scalar, v.z / scalar};
    }

    // Vector3 / operator overload (2 vectors element-wise division)
    Vector3 operator/(const Vector3& v1, const Vector3& v2) {
        return {v1.x / v2.x, v1.y / v2.y, v1.z / v2.z};
    }


    // Vector4 + operator overload
    Vector4 operator+(const Vector4& v1, const Vector4& v2) {
        return {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, v1.w + v2.w};
    }

    // Vector4 - operator overload
    Vector4 operator-(const Vector4& v1, const Vector4& v2) {
        return {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.w - v2.w};
    }

    // Vector4 * operator overload
    Vector4 operator*(const Vector4& v, double scalar) {
        return {v.x * scalar, v.y * scalar, v.z * scalar, v.w * scalar};
    }

    // Vector4 * operator overload (2 vectors element-wise multiplication)
    Vector4 operator*(const Vector4& v1, const Vector4& v2) {
        return {v1.x * v2.x, v1.y * v2.y, v1.z * v2.z, v1.w * v2.w};
    }

    // Vector4 / operator overload
    Vector4 operator/(const Vector4& v, double scalar) {
        return {v.x / scalar, v.y / scalar, v.z / scalar, v.w / scalar};
    }

    // Vector4 / operator overload (2 vectors element-wise division)
    Vector4 operator/(const Vector4& v1, const Vector4& v2) {
        return {v1.x / v2.x, v1.y / v2.y, v1.z / v2.z, v1.w / v2.w};
    }


    // Vector2 Cout operator overload
    std::ostream& operator<<(std::ostream& os, const Vector2& v) {
        os << "(" << v.x << ", " << v.y << ")" << std::endl;
        return os;
    }

    // Vector3 Cout operator overload
    std::ostream& operator<<(std::ostream& os, const Vector3& v) {
        os << "(" << v.x << ", " << v.y << ", " << v.z << ")" << std::endl;
        return os;
    }

    // Vector4 Cout operator overload
    std::ostream& operator<<(std::ostream& os, const Vector4& v) {
        os << "(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")" << std::endl;
        return os;
    }

    std::istream& operator>>(std::istream& is, Vector2& v) {
        is >> v.x >> v.y;
        return is;
    }

    std::istream& operator>>(std::istream& is, Vector3& v) {
        is >> v.x >> v.y >> v.z;
        return is;
    }

    std::istream& operator>>(std::istream& is, Vector4& v) {
        is >> v.x >> v.y >> v.z >> v.w;
        return is;
    }


    // Matrix2 + operator overload
    Matrix2 operator+(const Matrix2& m1, const Matrix2& m2) {
        return {m1.m00 + m2.m00, m1.m01 + m2.m01, m1.m10 + m2.m10, m1.m11 + m2.m11};
    }

    // Matrix2 - operator overload
    Matrix2 operator-(const Matrix2& m1, const Matrix2& m2) {
        return {m1.m00 - m2.m00, m1.m01 - m2.m01, m1.m10 - m2.m10, m1.m11 - m2.m11};
    }

    // Matrix2 * operator overload
    Matrix2 operator*(const Matrix2& m, double scalar) {
        return {m.m00 * scalar, m.m01 * scalar, m.m10 * scalar, m.m11 * scalar};
    }

    // Matrix2 * operator overload (2 matrices element-wise multiplication)
    Matrix2 operator*(const Matrix2& m1, const Matrix2& m2) {
        return {m1.m00 * m2.m00, m1.m01 * m2.m01, m1.m10 * m2.m10, m1.m11 * m2.m11};
    }

    // Matrix2 / operator overload
    Matrix2 operator/(const Matrix2& m, double scalar) {
        return {m.m00 / scalar, m.m01 / scalar, m.m10 / scalar, m.m11 / scalar};
    }

    // Matrix2 / operator overload (2 matrices element-wise division)
    Matrix2 operator/(const Matrix2& m1, const Matrix2& m2) {
        return {m1.m00 / m2.m00, m1.m01 / m2.m01, m1.m10 / m2.m10, m1.m11 / m2.m11};
    }


    // Matrix3 + operator overload
    Matrix3 operator+(const Matrix3& m1, const Matrix3& m2) {
        return {m1.m00 + m2.m00, m1.m01 + m2.m01, m1.m02 + m2.m02, m1.m10 + m2.m10, m1.m11 + m2.m11, m1.m12 + m2.m12, m1.m20 + m2.m20, m1.m21 + m2.m21, m1.m22 + m2.m22};
    }

    // Matrix3 - operator overload
    Matrix3 operator-(const Matrix3& m1, const Matrix3& m2) {
        return {m1.m00 - m2.m00, m1.m01 - m2.m01, m1.m02 - m2.m02, m1.m10 - m2.m10, m1.m11 - m2.m11, m1.m12 - m2.m12, m1.m20 - m2.m20, m1.m21 - m2.m21, m1.m22 - m2.m22};
    }

    // Matrix3 * operator overload
    Matrix3 operator*(const Matrix3& m, double scalar) {
        return {m.m00 * scalar, m.m01 * scalar, m.m02 * scalar, m.m10 * scalar, m.m11 * scalar, m.m12 * scalar, m.m20 * scalar, m.m21 * scalar, m.m22 * scalar};
    }

    // Matrix3 * operator overload (2 matrices element-wise multiplication)
    Matrix3 operator*(const Matrix3& m1, const Matrix3& m2) {
        return {m1.m00 * m2.m00, m1.m01 * m2.m01, m1.m02 * m2.m02, m1.m10 * m2.m10, m1.m11 * m2.m11, m1.m12 * m2.m12, m1.m20 * m2.m20, m1.m21 * m2.m21, m1.m22 * m2.m22};
    }

    // Matrix3 / operator overload
    Matrix3 operator/(const Matrix3& m, double scalar) {
        return {m.m00 / scalar, m.m01 / scalar, m.m02 / scalar, m.m10 / scalar, m.m11 / scalar, m.m12 / scalar, m.m20 / scalar, m.m21 / scalar, m.m22 / scalar};
    }

    // Matrix3 / operator overload (2 matrices element-wise division)
    Matrix3 operator/(const Matrix3& m1, const Matrix3& m2) {
        return {m1.m00 / m2.m00, m1.m01 / m2.m01, m1.m02 / m2.m02, m1.m10 / m2.m10, m1.m11 / m2.m11, m1.m12 / m2.m12, m1.m20 / m2.m20, m1.m21 / m2.m21, m1.m22 / m2.m22};
    }


    // Matrix4 + operator overload
    Matrix4 operator+(const Matrix4& m1, const Matrix4& m2) {
        return {m1.m00 + m2.m00, m1.m01 + m2.m01, m1.m02 + m2.m02, m1.m03 + m2.m03, m1.m10 + m2.m10, m1.m11 + m2.m11, m1.m12 + m2.m12, m1.m13 + m2.m13, m1.m20 + m2.m20, m1.m21 + m2.m21, m1.m22 + m2.m22, m1.m23 + m2.m23, m1.m30 + m2.m30, m1.m31 + m2.m31, m1.m32 + m2.m32, m1.m33 + m2.m33};
    }

    // Matrix4 - operator overload
    Matrix4 operator-(const Matrix4& m1, const Matrix4& m2) {
        return {m1.m00 - m2.m00, m1.m01 - m2.m01, m1.m02 - m2.m02, m1.m03 - m2.m03, m1.m10 - m2.m10, m1.m11 - m2.m11, m1.m12 - m2.m12, m1.m13 - m2.m13, m1.m20 - m2.m20, m1.m21 - m2.m21, m1.m22 - m2.m22, m1.m23 - m2.m23, m1.m30 - m2.m30, m1.m31 - m2.m31, m1.m32 - m2.m32, m1.m33 - m2.m33};
    }

    // Matrix4 * operator overload
    Matrix4 operator*(const Matrix4& m, double scalar) {
        return {m.m00 * scalar, m.m01 * scalar, m.m02 * scalar, m.m03 * scalar, m.m10 * scalar, m.m11 * scalar, m.m12 * scalar, m.m13 * scalar, m.m20 * scalar, m.m21 * scalar, m.m22 * scalar, m.m23 * scalar, m.m30 * scalar, m.m31 * scalar, m.m32 * scalar, m.m33 * scalar};
    }

    // Matrix4 * operator overload (2 matrices element-wise multiplication)
    Matrix4 operator*(const Matrix4& m1, const Matrix4& m2) {
        return {m1.m00 * m2.m00, m1.m01 * m2.m01, m1.m02 * m2.m02, m1.m03 * m2.m03, m1.m10 * m2.m10, m1.m11 * m2.m11, m1.m12 * m2.m12, m1.m13 * m2.m13, m1.m20 * m2.m20, m1.m21 * m2.m21, m1.m22 * m2.m22, m1.m23 * m2.m23, m1.m30 * m2.m30, m1.m31 * m2.m31, m1.m32 * m2.m32, m1.m33 * m2.m33};
    }

    // Matrix4 / operator overload
    Matrix4 operator/(const Matrix4& m, double scalar) {
        return {m.m00 / scalar, m.m01 / scalar, m.m02 / scalar, m.m03 / scalar, m.m10 / scalar, m.m11 / scalar, m.m12 / scalar, m.m13 / scalar, m.m20 / scalar, m.m21 / scalar, m.m22 / scalar, m.m23 / scalar, m.m30 / scalar, m.m31 / scalar, m.m32 / scalar, m.m33 / scalar};
    }

    // Matrix4 / operator overload (2 matrices element-wise division)
    Matrix4 operator/(const Matrix4& m1, const Matrix4& m2) {
        return {m1.m00 / m2.m00, m1.m01 / m2.m01, m1.m02 / m2.m02, m1.m03 / m2.m03, m1.m10 / m2.m10, m1.m11 / m2.m11, m1.m12 / m2.m12, m1.m13 / m2.m13, m1.m20 / m2.m20, m1.m21 / m2.m21, m1.m22 / m2.m22, m1.m23 / m2.m23, m1.m30 / m2.m30, m1.m31 / m2.m31, m1.m32 / m2.m32, m1.m33 / m2.m33};
    }

    // Matrix2 Cout operator overload
    std::ostream& operator<<(std::ostream& os, const Matrix2& m) {
        os << m.m00 << " " << m.m01 << std::endl;
        os << m.m10 << " " << m.m11 << std::endl;
        return os;
    }

    // Matrix3 Cout operator overload
    std::ostream& operator<<(std::ostream& os, const Matrix3& m) {
        os << m.m00 << " " << m.m01 << " " << m.m02 << std::endl;
        os << m.m10 << " " << m.m11 << " " << m.m12 << std::endl;
        os << m.m20 << " " << m.m21 << " " << m.m22 << std::endl;
        return os;
    }

    // Matrix4 Cout operator overload
    std::ostream& operator<<(std::ostream& os, const Matrix4& m) {
        os << m.m00 << " " << m.m01 << " " << m.m02 << " " << m.m03 << std::endl;
        os << m.m10 << " " << m.m11 << " " << m.m12 << " " << m.m13 << std::endl;
        os << m.m20 << " " << m.m21 << " " << m.m22 << " " << m.m23 << std::endl;
        os << m.m30 << " " << m.m31 << " " << m.m32 << " " << m.m33 << std::endl;
        return os;
    }
















    double factorial(int n) {
        double result = 1;
        for (int i = 2; i <= n; ++i) {
            result *= i;
        }
        return result;
    }

    double abs(double x) {
        return x < 0 ? -x : x;
    }

    double fmod(double x, double y) {
        return x - y * floor(x / y);
    }

    double fabs(double x) {
        return x < 0 ? -x : x;
    }

    double deg2rad(double x) {
        return x * PI / 180;
    }

    double rad2deg(double x) {
        return x * 180 / PI;
    }


    double pow(double x, double y) {
        double result = 1;
        for (int i = 0; i < y; ++i) {
            result *= x;
        }
        return result;
    }

    double sqrt(double x) {
        double guess = x / 2;
        double error = 0.0001;
        while (fabs(guess * guess - x) > error) {
            guess = (guess + x / guess) / 2;
        }
        return guess;
    }

    double cbrt(double x) {
        double guess = x / 3;
        double error = 0.0001;
        while (fabs(guess * guess * guess - x) > error) {
            guess = (2 * guess + x / (guess * guess)) / 3;
        }
        return guess;
    }

    double hypot(double x, double y) {
        return sqrt(x * x + y * y);
    }


    double exp(double x) {
        double sum = 1;
        double term = 1;
        int terms = 10; // Number of terms in the Taylor series
        for (int i = 1; i < terms; ++i) {
            term *= x / i;
            sum += term;
        }
        return sum;
    }

    double log(double x) {
        double sum = 0;
        double term = (x - 1) / (x + 1);
        int terms = 10; // Number of terms in the Taylor series
        for (int i = 0; i < terms; ++i) {
            sum += term / (2 * i + 1);
            term *= term * term;
        }
        return 2 * sum;
    }

    double log10(double x) {
        return log(x) / log(10);
    }

    double log2(double x) {
        return log(x) / log(2);
    }

    double expm1(double x) {
        double sum = 1;
        double term = 1;
        int terms = 10; // Number of terms in the Taylor series
        for (int i = 1; i < terms; ++i) {
            term *= x / i;
            sum += term;
        }
        return sum - 1;
    }

    double log1p(double x) {
        double sum = 0;
        double term = x / (1 - x);
        int terms = 10; // Number of terms in the Taylor series
        for (int i = 0; i < terms; ++i) {
            sum += term / (i + 1);
            term *= x / (1 - x);
        }
        return sum;
    }


    double round(double x) {
        return static_cast<int>(x + 0.5);
    }

    double floor(double x) {
        return static_cast<int>(x);
    }

    double ceil(double x) {
        return static_cast<int>(x) + 1;
    }

    double trunc(double x) {
        return x < 0 ? ceil(x) : floor(x);
    }

    double frac(double x) {
        return x - floor(x);
    }


    double min(double x, double y) {
        return x < y ? x : y;
    }

    double max(double x, double y) {
        return x > y ? x : y;
    }


    double sign(double x) {
        return x > 0 ? 1 : x < 0 ? -1 : 0;
    }

    double delta(double x, double y) {
        return fabs(x - y);
    }

    double average(double x, double y) {
        return (x + y) / 2;
    }

    double square(double x) {
        return x * x;
    }

    double cube(double x) {
        return x * x * x;
    }

    double mod(double x, double y) {
        return x - y * floor(x / y);
    }


    // Approximate sin function using Taylor series
    double sin(double x){
        // Normalize x to the range [-PI, PI]
        x = fmod(x, 2 * PI);
        if (x > PI) {
            x -= 2 * PI;
        } else if (x < -PI) {
            x += 2 * PI;
        }

        double sum = 0;
        int terms = 10; // Number of terms in the Taylor series
        for (int i = 0; i < terms; ++i) {
            double term = pow(-1, i) * pow(x, 2 * i + 1) / factorial(2 * i + 1);
            sum += term;
        }
        return sum;
    }

    // Approximate cos function using Taylor series
    double cos(double x){
        // Normalize x to the range [-PI, PI]
        x = fmod(x, 2 * PI);
        if (x > PI) {
            x -= 2 * PI;
        } else if (x < -PI) {
            x += 2 * PI;
        }

        double sum = 0;
        int terms = 10; // Number of terms in the Taylor series
        for (int i = 0; i < terms; ++i) {
            double term = pow(-1, i) * pow(x, 2 * i) / factorial(2 * i);
            sum += term;
        }
        return sum;
    }

    // Approximate tan function using Taylor series
    double tan(double x){
        return sin(x) / cos(x);
    }

    // Approximate asin function using Taylor series
    double asin(double x){
        double sum = 0;
        int terms = 10; // Number of terms in the Taylor series
        for (int i = 0; i < terms; ++i) {
            double term = factorial(2 * i) * pow(x, 2 * i + 1) / (pow(4, i) * pow(factorial(i), 2) * (2 * i + 1));
            sum += term;
        }
        return sum;
    }

    // Approximate acos function using Taylor series
    double acos(double x){
        return PI / 2 - asin(x);
    }

    // Approximate atan function using Taylor series
    double atan(double x){
        double sum = 0;
        int terms = 10; // Number of terms in the Taylor series
        for (int i = 0; i < terms; ++i) {
            double term = pow(-1, i) * pow(x, 2 * i + 1) / (2 * i + 1);
            sum += term;
        }
        return sum;
    }

    // Approximate atan2 function using atan
    double atan2(double y, double x){
        if (x > 0) {
            return atan(y / x);
        } else if (x < 0 && y >= 0) {
            return atan(y / x) + PI;
        } else if (x < 0 && y < 0) {
            return atan(y / x) - PI;
        } else if (x == 0 && y > 0) {
            return PI / 2;
        } else if (x == 0 && y < 0) {
            return -PI / 2;
        } else {
            return 0; // Undefined
        }
    }
 }
 
