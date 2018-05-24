//
// Created by Sandra Alfaro on 24/05/18.
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
                                                      int &base1, int &base2, int &base3, int &base4);


        /// Computes an approximation of the best LCP (directional) from Q to P
        /// and the rigid transformation that realizes it. The input sets may or may
        /// not contain normal information for any point.
        /// @param [in] P The first input set.
        /// @param [in] Q The second input set.
        /// as a fraction of the size of P ([0..1]).
        /// @param [out] transformation Rigid transformation matrix (4x4) that brings
        /// Q to the (approximate) optimal LCP. Initial value is considered as a guess
        /// @return the computed LCP measure.
        /// The method updates the coordinates of the second set, Q, applying
        /// the found transformation.
        template <typename Sampler = DefaultSampler,
                typename Visitor = DummyTransformVisitor>
        Scalar ComputeTransformation(const std::vector<Point3D>& P,
                                     std::vector<Point3D>* Q,
                                     Eigen::Ref<MatrixType> transformation,
                                     const Sampler& sampler = Sampler(),
                                     const Visitor& v = Visitor());

        /// Performs n RANSAC iterations, each one of them containing base selection,
        /// finding congruent sets and verification. Returns true if the process can be
        /// terminated (the target LCP was obtained or the maximum number of trials has
        /// been reached), false otherwise.
        template <typename Visitor>
        bool Perform_N_steps(int n,
                             Eigen::Ref<MatrixType> transformation,
                             std::vector<Point3D>* Q,
                             const Visitor& v);

        /// Tries one base and finds the best transformation for this base.
        /// Returns true if the achieved LCP is greater than terminate_threshold_,
        /// else otherwise.
        template <typename Visitor>
        bool TryOneBase(const Visitor &v);

        /// Selects a quadrilateral from P and returns the corresponding invariants
        /// and point indices. Returns true if a quadrilateral has been found, false
        /// otherwise.
        bool SelectQuadrilateral(Scalar &invariant1, Scalar &invariant2,
                                 int &id1, int &id2, int &id3, int &id4);

        /// Loop over the set of congruent 4-points and test the compatibility with the
        /// input base.
        /// \param [out] Nb Number of quads corresponding to valid configurations
        template <typename Visitor>
        bool TryCongruentSet(int base_id1,
                                                            int base_id2,
                                                            int base_id3,
                                                            int base_id4,
                                                            const std::vector<Quadrilateral> &congruent_quads,
                                                            const Visitor &v,
                                                            size_t &nbCongruent);
        // Initialize all internal data structures and data members.
        void Initialize(const std::vector<Point3D>& /*P*/,
                               const std::vector<Point3D>& /*Q*/) override;

    };
}

#include "match4pcsBase.hpp"

#endif //OPENGR_MATCH4PCSBASE_H
