#include <iostream>

#pragma once
    // Simple math library to replace cmath.
    // This library is not optimized for performance, but for simplicity and readability.
    // Made with the intention of being used in educational purposes.
    // Author: Spencer Yonce
namespace smath 
{
    
    // Constants

    // PI approximated to 20 decimal places
    const double PI = 3.14159265358979323846;
    // Euler's number
    const double E = 2.71828182845904523536;
    // Golden ratio
    const double PHI = 1.61803398874989484820;

    struct Vector2
    {
        double x;
        double y;
    };

    struct Vector3
    {
        double x;
        double y;
        double z;
    };

    struct Vector4
    {
        double x;
        double y;
        double z;
        double w;
    };

    struct Matrix2
    {
        double m00;
        double m01;
        double m10;
        double m11;
    };

    struct Matrix3
    {
        double m00;
        double m01;
        double m02;
        double m10;
        double m11;
        double m12;
        double m20;
        double m21;
        double m22;
    };

    struct Matrix4
    {
        double m00;
        double m01;
        double m02;
        double m03;
        double m10;
        double m11;
        double m12;
        double m13;
        double m20;
        double m21;
        double m22;
        double m23;
        double m30;
        double m31;
        double m32;
        double m33;
    };

    // Vector Opeartor Overloads

    Vector2 operator+(const Vector2& a, const Vector2& b);
    Vector2 operator-(const Vector2& a, const Vector2& b);
    Vector2 operator*(const Vector2& a, double b);
    Vector2 operator*(const Vector2& a, const Vector2& b);
    Vector2 operator/(const Vector2& a, double b);
    Vector2 operator/(const Vector2& a, const Vector2& b);

    Vector3 operator+(const Vector3& a, const Vector3& b);
    Vector3 operator-(const Vector3& a, const Vector3& b);
    Vector3 operator*(const Vector3& a, double b);
    Vector3 operator*(const Vector3& a, const Vector3& b);
    Vector3 operator/(const Vector3& a, double b);
    Vector3 operator/(const Vector3& a, const Vector3& b);

    Vector4 operator+(const Vector4& a, const Vector4& b);
    Vector4 operator-(const Vector4& a, const Vector4& b);
    Vector4 operator*(const Vector4& a, double b);
    Vector4 operator*(const Vector4& a, const Vector4& b);
    Vector4 operator/(const Vector4& a, double b);
    Vector4 operator/(const Vector4& a, const Vector4& b);

    // Vector Cout and Cin Overloads

    std::ostream& operator<<(std::ostream& os, const Vector2& v);
    std::ostream& operator<<(std::ostream& os, const Vector3& v);
    std::ostream& operator<<(std::ostream& os, const Vector4& v);

    // Vector2 Cin operator overload
    // Example input: 1 2 (requires a space between the two numbers)
    std::istream& operator>>(std::istream& is, Vector2& v);
    // Vector3 Cin operator overload
    // Example input: 1 2 3 (requires a space between the three numbers)
    std::istream& operator>>(std::istream& is, Vector3& v);
    // Vector4 Cin operator overload
    // Example input: 1 2 3 4 (requires a space between the four numbers)
    std::istream& operator>>(std::istream& is, Vector4& v);

    // Matrix Operator Overloads

    Matrix2 operator+(const Matrix2& a, const Matrix2& b);
    Matrix2 operator-(const Matrix2& a, const Matrix2& b);
    Matrix2 operator*(const Matrix2& a, double b);
    Matrix2 operator*(const Matrix2& a, const Matrix2& b);
    Matrix2 operator/(const Matrix2& a, double b);
    Matrix2 operator/(const Matrix2& a, const Matrix2& b);

    Matrix3 operator+(const Matrix3& a, const Matrix3& b);
    Matrix3 operator-(const Matrix3& a, const Matrix3& b);
    Matrix3 operator*(const Matrix3& a, double b);
    Matrix3 operator*(const Matrix3& a, const Matrix3& b);
    Matrix3 operator/(const Matrix3& a, double b);
    Matrix3 operator/(const Matrix3& a, const Matrix3& b);

    Matrix4 operator+(const Matrix4& a, const Matrix4& b);
    Matrix4 operator-(const Matrix4& a, const Matrix4& b);
    Matrix4 operator*(const Matrix4& a, double b);
    Matrix4 operator*(const Matrix4& a, const Matrix4& b);
    Matrix4 operator/(const Matrix4& a, double b);
    Matrix4 operator/(const Matrix4& a, const Matrix4& b);

    // Matrix Cout and Cin Overloads

    std::ostream& operator<<(std::ostream& os, const Matrix2& m);
    std::ostream& operator<<(std::ostream& os, const Matrix3& m);
    std::ostream& operator<<(std::ostream& os, const Matrix4& m);





    // Helpers

    // Returns the factorial of a number.
     // A factorial is the product of an integer and all the integers below it.
     // For example, the factorial of 5 is 5 * 4 * 3 * 2 * 1 = 120
    double factorial(int n);
    // Returns the absolute value of a number.
    double abs(double x);
    // Returns the remainder of a division operation.
    double fmod(double x, double y);
    // Returns the absolute value of a number.
    double fabs(double x);
    double deg2rad(double x);
    double rad2deg(double x);


    // Power Functions

    // Returns the result of raising a number to a power.
    double pow(double x, double y);
    // Returns the square root of a number.
    double sqrt(double x);
    // Returns the cube root of a number.
    double cbrt(double x);
    // Returns the hypotenuse of a right-angled triangle given the lengths of the other two sides.
    double hypot(double x, double y);

    // Exponential and Logarithmic Functions

    // Returns e raised to the power of a number.
    double exp(double x);
    // Returns the natural logarithm of a number.
    double log(double x);
    // Returns the base 10 logarithm of a number.
    double log10(double x);
    // Returns the base 2 logarithm of a number.
    double log2(double x);
    // Returns the result of raising a number to the power of e.
    // This is the inverse of the natural logarithm function.
    double expm1(double x);
    // Returns the natural logarithm of a number plus 1.
    // This is the inverse of the exp function.
    double log1p(double x);

    // Rounding Functions

    // Returns the nearest integer to a number.
    // Rounds halfway cases away from zero.
    double round(double x);
    // Returns the largest integer less than or equal to a given number.
    double floor(double x);
    // Returns the smallest integer greater than or equal to a given number.
    double ceil(double x);
    // Returns the integer part of a number by removing any fractional digits.
    double trunc(double x);
    // Returns the fractional part of a number by removing any integer digits.
    double frac(double x);

    // Min and Max Functions

    // Returns the smaller of two numbers.
    double min(double x, double y);
    // Returns the larger of two numbers.
    double max(double x, double y);

    // Miscellaneous Functions

    // Returns the sign of a number.
    // Returns 1 if the number is positive, -1 if the number is negative, and 0 if the number is zero.
    double sign(double x);
    // Returns the absolute value of the difference between two numbers.
    double delta(double x, double y);
    // Returns the average of two numbers.
    double average(double x, double y);
    // Returns the square of a number.
    double square(double x);
    // Returns the cube of a number.
    double cube(double x);
    // Returns the remainder of a division operation.
    double mod(double x, double y);



    // Trigonometry

    // Returns the sine of an angle.
    // The angle must be in radians.
    // To convert degrees to radians, multiply by PI / 180 or use the deg2rad function.
    double sin(double x);
    // Returns the cosine of an angle.
    // The angle must be in radians.
    // To convert degrees to radians, multiply by PI / 180 or use the deg2rad function.
    double cos(double x);
    // Returns the tangent of an angle.
    // The angle must be in radians.
    // To convert degrees to radians, multiply by PI / 180 or use the deg2rad function.
    double tan(double x);
    // Returns the arcsine of a number.
    double asin(double x);
    // Returns the arccosine of a number.
    double acos(double x);
    // Returns the arctangent of a number.
    double atan(double x);
    // Returns the angle whose tangent is the quotient of two specified numbers.
    // For example, atan2(1, 1) returns PI / 4 or 45 degrees.
    double atan2(double y, double x);



}

