#include <vector>
#include "PointwiseIterativeR2.h"
#include "../extra/math.h"

using namespace gbc;

void gbc::PointwiseIterativeR2::compute(const VertexR2 &p, std::vector<double> &b)
{
    b.clear();
    b.resize(n_v, 0.0);
    // Boundary.
    if (computeBoundaryCoordinates(p, b)) return;
    // Interior.
    bcInterior(p, b);
    return;
}


void gbc::PointwiseIterativeR2::bcInterior(const VertexR2 &p, std::vector<double> &b)
{
    std::vector<VertexR2> p1;
    for (int i = 0; i < n_v; ++i)
    {
        p1.push_back(_v[i]);
    }

    // 平移
    for (auto &p1_i : p1)
    {
        p1_i -= p;
    }

    // 计算p1
    std::vector<double> p1_coe(n_v, 0.0);
    double tmp;
    for (int i = 0; i < n_v; ++i)
    {
        double tmp = p1[i].length();
        p1[i] *= 1 / tmp;
        p1_coe[i] += 1 / tmp;
    }

    // 计算p2
    std::vector<VertexR2> p2(n_v);
    std::vector<double> p2_coe(n_v, 0.0);
    int ip;
    for (int i = 0; i < n_v; ++i)
    {   
        ip = (i + 1) % n_v;
        p2[i] = VertexR2(p1[i][0] + p1[ip][0], p1[i][1] + p1[ip][1]);
        tmp = p2[i].length();
        p2[i] *= 1 / tmp;
        p2_coe[i] += 1 / tmp;
    }

    // 迭代
    for (int i = 0; i < n_v; ++i)
    {
        VertexR2 delta = p1[i];
        std::vector<double> delta_coe(n_v, 0.0);
        delta_coe[i] = p1_coe[i];
        update_delta(p1, p1_coe, p2, p2_coe, delta, delta_coe, i);
        double sum_delta = gbc::sum(delta_coe);
        gbc::div(delta_coe, sum_delta);
        gbc::add(b, delta_coe);    
    }
    double sum_b = gbc::sum(b);
    gbc::div(b, sum_b);
}

void gbc::PointwiseIterativeR2::update_delta(std::vector<VertexR2>& p1, std::vector<double>& p1_coe, 
            std::vector<VertexR2>& p2,std::vector<double>& p2_coe,
            VertexR2& delta, std::vector<double>& delta_coe, int whichp)
{
    double delta_norm = delta.length();
    if (delta_norm > 1.0)
    {
        delta_norm = 1.0;
    }

    while (delta_norm > _tol)
    {
        double alpha = std::asin(delta_norm);
        // alpha = 0;
        std::vector<double> rho1_coe(n_v, 0.0);
        std::vector<double> rho2_coe(n_v, 0.0);

        get_rhoi_information(p1, alpha, delta, delta_coe, whichp, rho1_coe);
        get_rhoi_information(p2, alpha, delta, delta_coe, -1, rho2_coe);
        
        double sum = 0.0;
        for (int i = 0; i < n_v; ++i)
        {
            sum += rho1_coe[i] + rho2_coe[i];
        }

        std::vector<double> p_sum_coe(n_v, 0.0);
        VertexR2 p_sum;
        double ip;
        for (int i = 0; i < n_v; ++i)
        {   
            ip = (i + 1) % n_v;
            p_sum_coe[i] += (rho1_coe[i] * p1_coe[i] + rho2_coe[i] * p2_coe[i] * p1_coe[i]) / sum;
            p_sum_coe[ip] += (rho2_coe[i] * p2_coe[i] * p1_coe[ip]) / sum;
            p_sum[0] += (rho1_coe[i] * p1[i][0] + rho2_coe[i] * p2[i][0]) / sum;
            p_sum[1] += (rho1_coe[i] * p1[i][1] + rho2_coe[i] * p2[i][1]) / sum;
        }
        double p_sum_norm = p_sum.length();
        for (int i = 0; i < n_v; ++i)
        {
            p_sum_coe[i] /= p_sum_norm;
        }
        p_sum[0] /= p_sum_norm;
        p_sum[1] /= p_sum_norm;

        
        VertexR2 del_p = gbc::diff(delta, p_sum);
        double d_coe = 0.0 - gbc::dot(p_sum, del_p) / gbc::dot(del_p, del_p);

        delta =  d_coe * (delta - p_sum) + p_sum;
        for (int i = 0; i < n_v; ++i)
        {
            delta_coe[i] = d_coe * delta_coe[i] + (1 - d_coe) * p_sum_coe[i];
        }
        delta_norm = delta.length();
    }
}

void gbc::PointwiseIterativeR2::get_rhoi_information(std::vector<VertexR2>& pi, double alpha,
        VertexR2& delta, std::vector<double>& delta_coe, int whichp,
        std::vector<double>& rhoi_coe)
{
    double tmp, theta;
    double delta_norm = delta.length();
    for (int i = 0; i < n_v; ++i)
    {
        if (i != whichp)
        {
            tmp = (pi[i][0] * delta[0] + pi[i][1] * delta[1]) / (pi[i].length() * delta_norm);
            if (tmp > 1.0) 
            {
                tmp = 1.0;
            }
            if (tmp < -1.0)
            {
                tmp = -1.0;
            }
            theta = std::acos(tmp);
            if (theta > PI / 2 - alpha)
            {
                rhoi_coe[i] = f1((theta - PI/2 + alpha) / (PI/2 + alpha));
            }
        }
    }
}

