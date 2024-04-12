# SMATH Library Documentation

## Overview

The `smath` library is a simple mathematical library designed to replace `cmath`. It prioritizes simplicity and readability over performance, making it ideal for educational purposes.

### Author
Spencer Yonce

---

## Constants

### Mathematical Constants

- `PI`: π approximated to 20 decimal places (`3.14159265358979323846`).
- `E`: Euler's number (`2.71828182845904523536`).
- `PHI`: Golden ratio (`1.61803398874989484820`).

---

## Data Structures

### Vectors

1. **Vector2**
   - Components: `double x, double y`
   
2. **Vector3**
   - Components: `double x, double y, double z`
   
3. **Vector4**
   - Components: `double x, double y, double z, double w`

### Matrices

1. **Matrix2**
   - Components: `m00, m01, m10, m11`
   
2. **Matrix3**
   - Components: `m00 to m22`
   
3. **Matrix4**
   - Components: `m00 to m33`

---

## Functions

### Vector Operations

#### Overloads for Vector2, Vector3, and Vector4

- Addition (`+`)
- Subtraction (`-`)
- Multiplication (`*`) - Supports scalar multiplication and vector-wise multiplication.
- Division (`/`) - Supports scalar division and vector-wise division.

### Matrix Operations

#### Overloads for Matrix2, Matrix3, and Matrix4

- Addition (`+`)
- Subtraction (`-`)
- Multiplication (`*`) - Supports scalar multiplication and matrix multiplication.
- Division (`/`) - Supports scalar division.

### I/O Operations for Vectors and Matrices

- Output stream overloads for `Vector2`, `Vector3`, `Vector4`, `Matrix2`, `Matrix3`, and `Matrix4`.
- Input stream overloads for `Vector2`, `Vector3`, and `Vector4`.

### Mathematical Functions

#### Power Functions

- `pow(double x, double y)`
- `sqrt(double x)`
- `cbrt(double x)`
- `hypot(double x, double y)`

#### Exponential and Logarithmic Functions

- `exp(double x)`
- `log(double x)`
- `log10(double x)`
- `log2(double x)`
- `expm1(double x)`
- `log1p(double x)`

#### Rounding Functions

- `round(double x)`
- `floor(double x)`
- `ceil(double x)`
- `trunc(double x)`
- `frac(double x)`

#### Min and Max Functions

- `min(double x, double y)`
- `max(double x, double y)`

#### Miscellaneous Functions

- `sign(double x)`
- `delta(double x, double y)`
- `average(double x, double y)`
- `square(double x)`
- `cube(double x)`
- `mod(double x, double y)`

### Trigonometry Functions

- `sin(double x)`
- `cos(double x)`
- `tan(double x)`
- `asin(double x)`
- `acos(double x)`
- `atan(double x)`
- `atan2(double y, double x)`

### Helper Functions

- `factorial(int n)`
- `abs(double x)`
- `fmod(double x, double y)`
- `fabs(double x)`
- `deg2rad(double x)`
- `rad2deg(double x)`

---

## Future Directions

### Potential Math Types and Functions

1. **Complex Numbers**
   - Addition, subtraction, multiplication, division, and other complex number operations.

2. **Enhanced Vector and Matrix Operations**
   - Include functions for dot product, cross product, determinant, and inverse.

3. **Geometric Functions**
   - Functions for computing distances, angles, and intersections in both 2D and 3D spaces.

4. **Statistics Functions**
   - Mean, median, mode, variance, standard deviation.

5. **Discrete Mathematics Functions**
   - Functions related to graph theory, combinatorics, number theory, etc.

This documentation aims to be comprehensive and modular, allowing for easy updates and expansion as new features are developed in the `smath` library.