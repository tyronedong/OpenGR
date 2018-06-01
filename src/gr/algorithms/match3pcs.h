//
// Created by Sandra Alfaro on 30/05/18.
//

#ifndef OPENGR_MATCH3PCSBASE_H
#define OPENGR_MATCH3PCSBASE_H

#include <vector>
#include "../shared.h"
#include "../sampling.h"
#include "../utils/logger.h"
#include "matchBase.h"

#ifdef TEST_GLOBAL_TIMINGS
#   include "super4pcs/utils/timer.h"
#endif

namespace gr {
    struct Traits3pcs {
        static constexpr int size() { return 3; }
        using Base = std::array<int,3>;
        using Set = std::vector<Base>;
        using Coordinates = std::array<Point3D, 3>;
    };

    // ----- 3PCS Options -----
    struct Match3PCSOptions : public MatchOptions{
        using Scalar = typename Point3D::Scalar;
        Match3PCSOptions() {}

        /// Maximum normal difference.
        Scalar max_normal_difference = -1;
        /// Maximum translation distance. Set negative to ignore
        Scalar max_translation_distance = -1;
        /// Maximum color RGB distance between corresponding vertices. Set negative to ignore
        Scalar max_color_distance = -1;
    };

    class Match3pcs : public MatchBase<Traits3pcs> {
    public:
        Match3pcs (const Match3PCSOptions& options
                , const Utils::Logger& logger);

        virtual ~Match3pcs();

        bool generateCongruents (Base& base, Set& congruent_quads) override;

        void Initialize(const std::vector<Point3D>& /*P*/,
                        const std::vector<Point3D>& /*Q*/) override {};


    };
}

#include "match3pcs.hpp"
#endif //OPENGR_MATCH3PCSBASE_H
