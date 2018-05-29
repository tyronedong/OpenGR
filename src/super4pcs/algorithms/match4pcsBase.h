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

    struct Traits4pcs {
        static constexpr int size() { return 4; }
        using Base = std::array<int,4>;
        using Set = std::vector<Base>;
        using Coordinates = std::array<Point3D, 4>;
    /* .... */
    };

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

        /// Selects a quadrilateral from P and returns the corresponding invariants
        /// and point indices. Returns true if a quadrilateral has been found, false
        /// otherwise.
        bool SelectQuadrilateral(Scalar &invariant1, Scalar &invariant2,
                                 int& base1, int& base2, int& base3, int& base4);

        /// Initialize all internal data structures and data members.
        void Initialize(const std::vector<Point3D>& /*P*/,
                        const std::vector<Point3D>& /*Q*/) override;


        //TODO New Fonctions

        //Redefinition TryOneBase
        bool generateCongruents (Base& base,Set& congruent_quads) override;

    };
}

#include "match4pcsBase.hpp"

#endif //OPENGR_MATCH4PCSBASE_H
