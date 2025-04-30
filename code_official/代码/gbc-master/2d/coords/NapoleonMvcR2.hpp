
#ifndef GBC_NAPOLEONMVC_HPP
#define GBC_NAPOLEONMVC_HPP
// STL includes.
#include <vector>
#include <cassert>
#include <cmath>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/math.h"
#include "../extra/BarycentricCoordinatesR2.hpp"

namespace gbc {

    // PointwiseIterativeCoordinates in R2.
    class NapoleonMvcR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        NapoleonMvcR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol) 
        {
            n_v = _v.size();
            _theta = 2 * gbc::PI / n_v; 
            _tan_theta = std::tan(_theta);
            _tan_theta = 1;
        }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "NapoleonMvcR2";
        }
        
        void compute(const VertexR2 &p, std::vector<double> &b)
        {
            b.clear();
            b.resize(n_v, 0.0);
            // Boundary.
            if (computeBoundaryCoordinates(p, b)) return;
            // Interior.
            bcInterior(p, b);
            return;
        }

        // Compute the coordinates at p using the internal storage from the VertexR2 class.
        inline void compute(VertexR2 &p){
            compute(p, p.b());
        }

        // Compute coordinates bb at all points p in the vector.
        void compute(const std::vector<VertexR2> &p, std::vector<std::vector<double> > &bb)
        {
            const size_t numP = p.size();
            
            bb.resize(numP);
            for (size_t i = 0; i < numP; ++i) 
            {   
                compute(p[i], bb[i]);
            }
        }

        // Compute coordinates at all points p in the vector using the internal storage from the VertexR2 class.
        void compute(std::vector<VertexR2> &p){

            const size_t numP = p.size();
            for (size_t i = 0; i < numP; ++i) 
            {
                compute(p[i], p[i].b());
            }
        }

        // Implementation of the virtual function to compute all coordinates.
        inline void bc(std::vector<VertexR2> &p) {
            compute(p);
        }

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;
        double _theta = 0.0;
        double _tan_theta = 0.0;

        void bcInterior(const VertexR2 &p, std::vector<double> &b)
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

            // 计算r_i
            std::vector<double> p1_coe(n_v, 0.0);
            double tmp;
            for (int i = 0; i < n_v; ++i)
            {
                double tmp = p1[i].length();
                p1[i] *= 1 / tmp;
                p1_coe[i] += 1 / tmp;
            }

            // 计算t_i
            std::vector<double> t(n_v, 0.0);
            for (int i = 0; i < n_v; ++i)
            {
                int ip = (i + 1) % n_v;
                // int sign = -1;
                // double side = gbc::cross(p1[i], p1[ip]);
                // if (side < 0)
                // {
                //     sign *= -1;
                // }
                // double ti_2 = (p1[i] + p1[ip]).length();
                // double cos_a = gbc::dot(p1[i], p1[ip]);
                // double ti_1 = std::sqrt((1 + cos_a) / 2) + sign * std::sqrt((1 - cos_a) / 2) * _tan_theta;
                // t[i] = ti_1 / ti_2;

                VertexR2 v1 = (p1[i] - p1[ip]) / 2;
                VertexR2 v2 = (p1[i] + p1[ip]) / 2;
                gbc::rotateSomeAngel(v1);
                v2 += v1;
                t[i] = v2.length() / (p1[i] + p1[ip]).length();
                
            }

            // 计算重心坐标
            double sum = 0.0;
            for (int i = 0; i < n_v; ++i)
            {
                int ip = (i - 1 + n_v) % n_v;
                sum += (t[i] + t[ip] - 1) * p1_coe[i];
                b[i] = (t[i] + t[ip] - 1) * p1_coe[i];
            }
            div(b, sum);
        }
    

    };

} // namespace gbc

#endif // GBC_GBC_NAPOLEONMVC_HPP
