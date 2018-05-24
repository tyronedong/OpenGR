//
// Created by Sandra alfaro on 24/05/18.
//

#ifndef OPENGR_MATCH4PCSBASE_H
#define OPENGR_MATCH4PCSBASE_H

#include <vector>

#ifdef SUPER4PCS_USE_OPENMP
#include <omp.h>
#endif

#include "../shared4pcs.h"
#include "../sampling.h"
#include "../accelerators/kdtree.h"
#include "../utils/logger.h"
#include "matchBase.h"

#ifdef TEST_GLOBAL_TIMINGS
#   include "super4pcs/utils/timer.h"
#endif

namespace gr {

    template <typename Functor>
    class Match4pcsBase : public MatchBase {
        using Base = MatchBase;

    protected:
        /// The transformation matrix by wich we transform Q to P
        Eigen::Matrix<Scalar, 4, 4> transform_;
        /// Quad centroids in first and second clouds
        /// They are used temporarily and makes the transformations more robust to
        /// noise. At the end, the direct transformation applied as a 4x4 matrix on
        /// every points in Q is computed and returned.
        Eigen::Matrix<Scalar, 3, 1> qcentroid1_, qcentroid2_;
        /// The points in the base (indices to P). It is being updated in every
        /// RANSAC iteration.
        int base_[4];
        /// The current congruent 4 points from Q. Every RANSAC iteration the
        /// algorithm examines a set of such congruent 4-points from Q and retains
        /// the best from them (the one that realizes the best LCP).
        int current_congruent_[4];
        Functor fun_;

    public:

        inline Match4pcsBase<Functor>::Match4pcsBase (const Match4PCSOptions& options
                , const Utils::Logger& logger);

        virtual ~Match4pcsBase();

        inline const Functor& Match4pcsBase<Functor>::getFunctor() const { return fun_; }

        /// Takes quadrilateral as a base, computes robust intersection point
        /// (approximate as the lines might not intersect) and returns the invariants
        /// corresponding to the two selected lines. The method also updates the order
        /// of the base base_3D_.
        inline bool Match4pcsBase<Functor>::TryQuadrilateral(Scalar &invariant1, Scalar &invariant2,
                                                      int &base1, int &base2, int &base3, int &base4);

        /// Tries one base and finds the best transformation for this base.
        /// Returns true if the achieved LCP is greater than terminate_threshold_,
        /// else otherwise.
        template <typename Visitor>
        inline bool Match4pcsBase<Functor>::TryOneBase(const Visitor &v);

        /// Selects a quadrilateral from P and returns the corresponding invariants
        /// and point indices. Returns true if a quadrilateral has been found, false
        /// otherwise.
        inline bool Match4pcsBase<Functor>::SelectQuadrilateral(Scalar &invariant1, Scalar &invariant2,
                                                                int& base1, int& base2, int& base3, int& base4);

        /// Loop over the set of congruent 4-points and test the compatibility with the
        /// input base.
        /// \param [out] Nb Number of quads corresponding to valid configurations
        template <typename Visitor>
        inline bool Match4pcsBase<Functor>::TryCongruentSet(int base_id1,
                                                            int base_id2,
                                                            int base_id3,
                                                            int base_id4,
                                                            const std::vector<Quadrilateral> &congruent_quads,
                                                            const Visitor &v,
                                                            size_t &nbCongruent);
    };
}
#endif //OPENGR_MATCH4PCSBASE_H
