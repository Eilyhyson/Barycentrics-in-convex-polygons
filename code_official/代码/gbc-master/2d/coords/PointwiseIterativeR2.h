
// README:
/*

    Mean value coordinates.

    This class depends on:
    1. BarycentricCoordinatesR2.hpp
    2. SegmentCoordinatesR2.hpp
    3. VertexExpressionsR2.hpp
    4. VertexR2.hpp

*/

#ifndef POINTWISEITERATIVER2_HPP
#define POINTWISEITERATIVER2_HPP

// STL includes.
#include <vector>
#include <cassert>
#include <cmath>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"

namespace gbc {

    // PointwiseIterativeCoordinates in R2.
    class PointwiseIterativeR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        PointwiseIterativeR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol) 
        {
            n_v = _v.size();
        }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "PointwiseIterativeR2";
        }
        
        void compute(const VertexR2 &p, std::vector<double> &b);

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
                // std::cout << i << std::endl;
                compute(p[i], p[i].b());
            }
            // compute(p[1044], p[1044].b());
        }

        // Implementation of the virtual function to compute all coordinates.
        inline void bc(std::vector<VertexR2> &p) {
            compute(p);
        }

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;

        void bcInterior(const VertexR2 &p, std::vector<double> &b);
    
        void update_delta(std::vector<VertexR2>& p1, std::vector<double>& p1_coe, std::vector<VertexR2>& p2, std::vector<std::vector<double>>& p2_coe,
            VertexR2& delta, std::vector<double>& delta_coe, int whichp);
    
        void get_rhoi_information(std::vector<VertexR2>& pi, double alpha,
            VertexR2& delta, std::vector<double>& delta_coe, int whichp, 
            std::vector<double>& rhoi_coe, std::vector<double>& d_coei);

        inline double f1(double x)
        {
            return x * x * x + 3 * x * x * (1 - x);
        };

        std::vector<VertexR2> p0;
    };

} // namespace gbc

#endif // GBC_POINTWISEITERATIVER2_HPP
