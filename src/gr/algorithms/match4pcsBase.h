//
// Created by Sandra Alfaro on 24/05/18.
//

#ifndef OPENGR_MATCH4PCSBASE_H
#define OPENGR_MATCH4PCSBASE_H

#include <vector>

#ifdef SUPER4PCS_USE_OPENMP
#include <omp.h>
#endif

#include "../shared.h"
#include "../sampling.h"
#include "../accelerators/kdtree.h"
#include "../utils/logger.h"
#include "matchBase.h"

#ifdef TEST_GLOBAL_TIMINGS
#   include "super4pcs/utils/timer.h"
#endif

namespace gr {


    // ----- 4PCS Options -----
    struct Match4PCSOptions : public MatchOptions{
        using Scalar = typename Point3D::Scalar;
        Match4PCSOptions() {}

        /// Maximum normal difference.
        Scalar max_normal_difference = -1;
        /// Maximum translation distance. Set negative to ignore
        Scalar max_translation_distance = -1;
        /// Maximum color RGB distance between corresponding vertices. Set negative to ignore
        Scalar max_color_distance = -1;
    };

    struct Traits4pcs {
        static constexpr int size() { return 4; }
        using Base = std::array<int,4>;
        using Set = std::vector<Base>;
        using Coordinates = std::array<Point3D, 4>;
    };

    /// Class for the computation of the 4PCS algorithm.
    /// \param Functor use to determinate the use of Super4pcs or 4pcs algorithm.
    template <typename Functor>
    class Match4pcsBase : public MatchBase<Traits4pcs> {
    protected:
        Functor fun_;

    public:

        Match4pcsBase (const Match4PCSOptions& options
                , const Utils::Logger& logger);

        virtual ~Match4pcsBase();

        inline const Functor& getFunctor() const { return fun_; }

        /// Takes quadrilateral as a base, computes robust intersection point
        /// (approximate as the lines might not intersect) and returns the invariants
        /// corresponding to the two selected lines. The method also updates the order
        /// of the base base_3D_.
        bool TryQuadrilateral(Scalar &invariant1, Scalar &invariant2,
                              int &id1, int &id2, int &id3, int &id4);

        /// Selects a random triangle in the set P (then we add another point to keep the
        /// base as planar as possible). We apply a simple heuristic that works in most
        /// practical cases. The idea is to accept maximum distance, computed by the
        /// estimated overlap, multiplied by the diameter of P, and try to have
        /// a triangle with all three edges close to this distance. Wide triangles helps
        /// to make the transformation robust while too large triangles makes the
        /// probability of having all points in the inliers small so we try to trade-off.
        bool SelectQuadrilateral(Scalar &invariant1, Scalar &invariant2,
                                 int& base1, int& base2, int& base3, int& base4);

        /// Initializes the data structures and needed values before the match
        /// computation.
        /// @param [in] point_P First input set.
        /// @param [in] point_Q Second input set.
        /// expected to be in the inliers.
        /// This method is called once the internal state of the Base class as been
        /// set.
        void Initialize(const std::vector<Point3D>& /*P*/,
                        const std::vector<Point3D>& /*Q*/) override;

        /// Find all the congruent set similar to the base in the second 3D model (Q).
        /// It could be with a 3 point base or a 4 point base.
        /// \param base use to find the similar points congruent in Q.
        /// \param congruent_set a set of all point congruent found in Q.
        bool generateCongruents (Base& base,Set& congruent_quads) override;

    };
}

#include "match4pcsBase.hpp"

#endif //OPENGR_MATCH4PCSBASE_H
