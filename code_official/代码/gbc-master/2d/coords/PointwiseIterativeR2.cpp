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
    p0.clear();
    std::vector<VertexR2> p1;
    for (int i = 0; i < n_v; ++i)
    {
        p1.push_back(_v[i]);
    }

    // 平移
    for (auto &p1_i : p1)
    {
        p1_i -= p;
        p0.push_back(p1_i);
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
    std::vector<std::vector<double>> p2_coe(n_v, std::vector<double>(2, 0.0));
    int ip;
    for (int i = 0; i < n_v; ++i)
    {   
        ip = (i + 1) % n_v;
        p2[i] = VertexR2(p1[i][0] + p1[ip][0], p1[i][1] + p1[ip][1]);
        tmp = p2[i].length();
        p2[i] *= 1 / tmp;
        p2_coe[i][0] += p1_coe[i] / tmp;
        p2_coe[i][1] += p1_coe[ip] / tmp;
    }

    // for (int i = 0; i < n_v; i++)
    // {
    //     ip = (i + 1) % n_v;
    //     VertexR2 p2_test = p2_coe[i][0] * p0[i] + p2_coe[i][1] * p0[ip];

    //     if ((p2_test - p2[i]).length() > 1e-10)
    //     {
    //         std::cout << "p2计算错误" << std::endl;
    //     }
    // }
    //  std::cout << "p2计算正确" << std::endl;

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
            std::vector<VertexR2>& p2,std::vector<std::vector<double>> & p2_coe,
            VertexR2& delta, std::vector<double>& delta_coe, int whichp)
{
    double delta_norm = delta.length();
    if (delta_norm > 1.0)
    {
        delta_norm = 1.0;
    }
    int times = 0;
    while (delta_norm > _tol)
    {
        times += 1;
        double alpha = std::asin(delta_norm);
        // alpha = 0;
        std::vector<double> rho1_coe(n_v, 0.0);
        std::vector<double> rho2_coe(n_v, 0.0);
        std::vector<double> d1_coe(n_v, 0.0);
        std::vector<double> d2_coe(n_v, 0.0);
        std::vector<double> new_delta_coe(n_v, 0.0);

        get_rhoi_information(p1, alpha, delta, delta_coe, whichp, rho1_coe, d1_coe);
        get_rhoi_information(p2, alpha, delta, delta_coe, -1, rho2_coe, d2_coe);
        
        double sum = 0.0;
        double sum_delta_rho = 0.0; // delta前的系数
        for (int i = 0; i < n_v; ++i)
        {
            sum += rho1_coe[i] + rho2_coe[i];
            sum_delta_rho += rho1_coe[i] * d1_coe[i] + rho2_coe[i] * d2_coe[i];
        }

        double ip;
        std::vector<VertexR2> d1(n_v);
        std::vector<VertexR2> d2(n_v);
        for (int i = 0; i < n_v; ++i)
        {
            ip = (i - 1 + n_v) % n_v;
            new_delta_coe[i] = sum_delta_rho * delta_coe[i] + rho1_coe[i] * (1 - d1_coe[i]) * p1_coe[i] 
                            + rho2_coe[i] * (1 - d2_coe[i]) * p2_coe[i][0]
                            + rho2_coe[ip] * (1 - d2_coe[ip]) * p2_coe[ip][1];
            d1[i] = rho1_coe[i] * (d1_coe[i] * delta + (1 - d1_coe[i]) * p1[i]);
            d2[i] = rho2_coe[i] * (d2_coe[i] * delta + (1 - d2_coe[i]) * p2[i]);
        }

        delta[0] = 0.0;
        delta[1] = 0.0;
        for (int i = 0; i < n_v; ++i)
        {
            delta_coe[i] = new_delta_coe[i] / sum;
            delta += d1[i] + d2[i];
        }
        delta = delta / sum;
        delta_norm = delta.length();

        // VertexR2 delta_tmp;
        // for (int i = 0; i < n_v; ++i)
        // {
        //     delta_tmp += delta_coe[i] * p0[i];
        // }
        // if ((delta_tmp - delta).length() < 1e-10)
        // {
        //     std::cout << "delta计算正确" << std::endl;
        // }
        // else{
        //     std::cout << "delta计算错误" << std::endl;
        // }
        if (times > 100)
        {
            int a =1;
        }
   }
}

void gbc::PointwiseIterativeR2::get_rhoi_information(std::vector<VertexR2>& pi, double alpha,
        VertexR2& delta, std::vector<double>& delta_coe, int whichp,
        std::vector<double>& rhoi_coe, std::vector<double>& di_coe)
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
                VertexR2 del_p = delta - pi[i];
                di_coe[i] = 0.0 - gbc::dot(pi[i], del_p) / gbc::dot(del_p, del_p);
            }
        }
    }
}

