#ifndef GBC_MATH_HPP
#define GBC_MATH_HPP

#include "VertexR2.hpp"
#include <vector>
namespace gbc
{
    const double PI = 3.14159265358979323846; // 自定义的 PI 常量

    double dot(VertexR2 a, VertexR2 b);

    double dot(double a, double b, double c, double d);

    double cross(VertexR2 a, VertexR2 b);

    double cross(double a, double b, double c, double d);

    VertexR2 diff(VertexR2 a, VertexR2 b);

    double sum(std::vector<double>& nums);

    void div(std::vector<double>& nums, double s);

    void add(std::vector<double>& nums1, std::vector<double>& nums2);

    void rotateSomeAngel(VertexR2 & a);
}

#endif // GBC_MATH_HPP