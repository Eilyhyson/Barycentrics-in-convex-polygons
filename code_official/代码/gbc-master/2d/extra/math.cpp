#include "math.h"
using namespace gbc;

double gbc::dot(VertexR2 a, VertexR2 b)
{
    return a[0] * b[0] + a[1] * b[1];
}

double gbc::dot(double a, double b, double c, double d)
{
    return a * c + b * d;
}

double gbc::cross(VertexR2 a, VertexR2 b)
{
    return a[0] * b[1] - a[1] * b[0];
}

double gbc::cross(double a, double b, double c, double d)
{
    return a * d - b * c;
}

VertexR2 gbc::diff(VertexR2 a, VertexR2 b)
{
    return VertexR2(a[0] - b[0], a[1] - b[1]);
}

double gbc::sum(std::vector<double>& nums)
{
    double ret = 0.0;
    for (int i = 0; i < nums.size(); ++i)
    {
        ret += nums[i];
    }
    return ret;
}

void gbc::div(std::vector<double>& nums, double s)
{
    for (int i = 0; i < nums.size(); ++i)
    {
        nums[i] /= s;
    }
}

void gbc::add(std::vector<double>& nums1, std::vector<double>& nums2)
{
    for (int i = 0; i < nums1.size(); ++i)
    {
        nums1[i] += nums2[i];
    }
}

//旋转90度
void gbc::rotateSomeAngel(VertexR2 & a)
{
    double tmp = a[0];
    a[0] = -a[1];
    a[1] = tmp;
}