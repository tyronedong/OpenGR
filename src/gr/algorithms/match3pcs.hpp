//
// Created by Sandra Alfaro on 30/05/18.
//

#include <vector>
#include <atomic>
#include <chrono>
#include "../shared.h"
#include "../sampling.h"
#include "../utils/logger.h"
#include "match3pcs.h"

namespace gr {

    Match3pcs::Match3pcs(const gr::Match3PCSOptions &options,
                                 const gr::Utils::Logger &logger)
                                 : MatchBase(options,logger)
    {
        base_3D_.resize(3);
    }

    Match3pcs::~Match3pcs() {};

    inline bool Match3pcs::generateCongruents (Base& base, Set& congruent_set) {

        //Find base in P (random triangle)
        if (!SelectRandomTriangle(base[0], base[1], base[2]))
            return false;
        base_3D_ [0] = sampled_P_3D_[base[0]];
        base_3D_ [1] = sampled_P_3D_[base[1]];
        base_3D_ [2] = sampled_P_3D_[base[2]];


        // Computes distance between points.
        const Scalar d1 = (base_3D_[0].pos()- base_3D_[1].pos()).norm();
        const Scalar d2 = (base_3D_[0].pos()- base_3D_[2].pos()).norm();
        const Scalar d3 = (base_3D_[1].pos()- base_3D_[2].pos()).norm();

       /*
        // Compute normal angles.
        const Scalar normal_angle_AB;
        const Scalar normal_angle_AC;
        const Scalar normal_angle_BC;
       */

        // Find all 3pcs in Q
        for (int i=0; i<sampled_Q_3D_.size(); ++i) {
            const Point3D& a = sampled_Q_3D_[i];
            for (int j=i+1; j<sampled_Q_3D_.size(); ++j) {
                const Point3D& b = sampled_Q_3D_[j];
                const Scalar dAB = (b.pos() - a.pos()).norm();
                if (std::abs(dAB - d1) > distance_factor * options_.delta) continue;
                for (int k=j+1; k<sampled_Q_3D_.size(); ++k) {
                    const Point3D& c = sampled_Q_3D_[k];
                    const Scalar dAC = (c.pos() - a.pos()).norm();
                    const Scalar dBC = (c.pos() - b.pos()).norm();
                    if (std::abs(dAC - d2) > distance_factor * options_.delta) continue;
                    if (std::abs(dBC - d3) > distance_factor * options_.delta) continue;
                    congruent_set.push_back({i,j,k});
                }
            }
        }

        //TODO add filter points
        /* //Change the normal_angle
         // Compute normal angles.
        const Scalar normal_angle1 = (base_3D_[0].normal() - base_3D_[1].normal()).norm();
        const Scalar normal_angle2 = (base_3D_[2].normal() - base_3D_[3].normal()).norm();

          PointFilterFunctor fun(myOptions_, myBase_3D_);
                    std::pair<bool,bool> res = fun(p,q, pair_normals_angle, base_point1,base_point2);
                    if (res.first)
                        pairs->emplace_back(i, j);
                    if (res.second)
                        pairs->emplace_back(j, i);
        */
        return congruent_set.size()!=0;
    }

}