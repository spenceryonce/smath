#include <iostream>
#include "smath.hpp"

int main() {
    double angle = 3.14159265358979323846 / 2; // 90 degrees in radians
    
    // Vector2
    smath::Vector2 v1 = {1, 2};
    smath::Vector2 v2 = {3, 4};

    // Add two vectors
    smath::Vector2 v3 = v1 + v2;
    std::cout << "v1 + v2 = " << v3 << std::endl;

    // Matrix2
    smath::Matrix2 m1 = {1, 2, 3, 4};
    smath::Matrix2 m2 = {5, 6, 7, 8};

    // Multiply two matrices
    smath::Matrix2 m3 = m1 * m2;

    std::cout << "m1 * m2 = " << std::endl;
    std::cout << m3 << std::endl;

    int exit;
    std::cout << "Enter 0 to exit: ";
    std::cin >> exit;
    if(exit == 0)
    {
        return 0;
    } 
    return 0;
}
