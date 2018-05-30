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

    class Match3pcs : public MatchBase<Traits3pcs> {
        Match3pcs (const MatchOptions& options
                , const Utils::Logger& logger);

        virtual ~Match3pcs();

        bool generateCongruents (Base& base, Set& congruent_quads) override;


    };
}

//#include "match3pcsBase.hpp"
#endif //OPENGR_MATCH3PCSBASE_H
