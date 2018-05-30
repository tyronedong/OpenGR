//
// Created by Sandra Alfaro on 24/05/18.
//

#include <vector>
#include <chrono>
#include <atomic>
#include <Eigen/Geometry>                 // MatrixBase.homogeneous()
#include <Eigen/SVD>
#include <Eigen/Core>                     // Transform.computeRotationScaling()


#ifdef SUPER4PCS_USE_OPENMP
#include <omp.h>
#endif

#include "super4pcs/shared.h"
#include "../sampling.h"
#include "../accelerators/kdtree.h"
#include "../utils/logger.h"
#include "match4pcsBase.h"

#ifdef TEST_GLOBAL_TIMINGS
#   include "super4pcs/utils/timer.h"
#endif

template <typename VectorType, typename Scalar>
static Scalar distSegmentToSegment(const VectorType& p1, const VectorType& p2,
                                   const VectorType& q1, const VectorType& q2,
                                   Scalar& invariant1, Scalar& invariant2) {

    static const Scalar kSmallNumber = 0.0001;
    VectorType u = p2 - p1;
    VectorType v = q2 - q1;
    VectorType w = p1 - q1;
    Scalar a = u.dot(u);
    Scalar b = u.dot(v);
    Scalar c = v.dot(v);
    Scalar d = u.dot(w);
    Scalar e = v.dot(w);
    Scalar f = a * c - b * b;
    // s1,s2 and t1,t2 are the parametric representation of the intersection.
    // they will be the invariants at the end of this simple computation.
    Scalar s1 = 0.0;
    Scalar s2 = f;
    Scalar t1 = 0.0;
    Scalar t2 = f;

    if (f < kSmallNumber) {
        s1 = 0.0;
        s2 = 1.0;
        t1 = e;
        t2 = c;
    } else {
        s1 = (b * e - c * d);
        t1 = (a * e - b * d);
        if (s1 < 0.0) {
            s1 = 0.0;
            t1 = e;
            t2 = c;
        } else if (s1 > s2) {
            s1 = s2;
            t1 = e + b;
            t2 = c;
        }
    }

    if (t1 < 0.0) {
        t1 = 0.0;
        if (-d < 0.0)
            s1 = 0.0;
        else if (-d > a)
            s1 = s2;
        else {
            s1 = -d;
            s2 = a;
        }
    } else if (t1 > t2) {
        t1 = t2;
        if ((-d + b) < 0.0)
            s1 = 0;
        else if ((-d + b) > a)
            s1 = s2;
        else {
            s1 = (-d + b);
            s2 = a;
        }
    }
    invariant1 = (std::abs(s1) < kSmallNumber ? 0.0 : s1 / s2);
    invariant2 = (std::abs(t1) < kSmallNumber ? 0.0 : t1 / t2);

    return ( w + (invariant1 * u) - (invariant2 * v)).norm();
}

namespace gr {
    template <typename Functor>
    Match4pcsBase<Functor>::Match4pcsBase (const Match4PCSOptions& options
            , const Utils::Logger& logger)
            : MatchBase(options,logger)
            , fun_(sampled_Q_3D_,base_3D_,options) {}

    template <typename Functor>
    Match4pcsBase<Functor>::~Match4pcsBase() {}

    template <typename Functor>
    bool Match4pcsBase<Functor>::TryQuadrilateral(Scalar &invariant1, Scalar &invariant2,
                                                  int &id1, int &id2, int &id3, int &id4) {

        Scalar min_distance = std::numeric_limits<Scalar>::max();
        int best1, best2, best3, best4;
        best1 = best2 = best3 = best4 = -1;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (i == j) continue;
                int k = 0;
                while (k == i || k == j) k++;
                int l = 0;
                while (l == i || l == j || l == k) l++;
                double local_invariant1;
                double local_invariant2;
                // Compute the closest points on both segments, the corresponding
                // invariants and the distance between the closest points.
                Scalar segment_distance = distSegmentToSegment(
                        base_3D_[i].pos(), base_3D_[j].pos(),
                        base_3D_[k].pos(), base_3D_[l].pos(),
                        local_invariant1, local_invariant2);
                // Retail the smallest distance and the best order so far.
                if (segment_distance < min_distance) {
                    min_distance = segment_distance;
                    best1 = i;
                    best2 = j;
                    best3 = k;
                    best4 = l;
                    invariant1 = local_invariant1;
                    invariant2 = local_invariant2;
                }
            }
        }

        if(best1 < 0 || best2 < 0 || best3 < 0 || best4 < 0 ) return false;

        std::vector<Point3D> tmp = base_3D_;
        base_3D_[0] = tmp[best1];
        base_3D_[1] = tmp[best2];
        base_3D_[2] = tmp[best3];
        base_3D_[3] = tmp[best4];

        std::array<int, 4> tmpId = {id1, id2, id3, id4};
        id1 = tmpId[best1];
        id2 = tmpId[best2];
        id3 = tmpId[best3];
        id4 = tmpId[best4];

        return true;
    }

    template <typename Functor>
    bool Match4pcsBase<Functor>::SelectQuadrilateral(Scalar &invariant1, Scalar &invariant2,
                                                     int& base1, int& base2, int& base3, int& base4)  {

        const Scalar kBaseTooSmall (0.2);
        int current_trial = 0;

        // Try fix number of times.
        while (current_trial < kNumberOfDiameterTrials) {
            // Select a triangle if possible. otherwise fail.
            if (!SelectRandomTriangle(base1, base2, base3)){
                return false;
            }

            base_3D_[0] = sampled_P_3D_[base1];
            base_3D_[1] = sampled_P_3D_[base2];
            base_3D_[2] = sampled_P_3D_[base3];

            // The 4th point will be a one that is close to be planar to the other 3
            // while still not too close to them.
            const double x1 = base_3D_[0].x();
            const double y1 = base_3D_[0].y();
            const double z1 = base_3D_[0].z();
            const double x2 = base_3D_[1].x();
            const double y2 = base_3D_[1].y();
            const double z2 = base_3D_[1].z();
            const double x3 = base_3D_[2].x();
            const double y3 = base_3D_[2].y();
            const double z3 = base_3D_[2].z();

            // Fit a plan.
            Scalar denom = (-x3 * y2 * z1 + x2 * y3 * z1 + x3 * y1 * z2 - x1 * y3 * z2 -
                            x2 * y1 * z3 + x1 * y2 * z3);

            if (denom != 0) {
                Scalar A =
                        (-y2 * z1 + y3 * z1 + y1 * z2 - y3 * z2 - y1 * z3 + y2 * z3) / denom;
                Scalar B =
                        (x2 * z1 - x3 * z1 - x1 * z2 + x3 * z2 + x1 * z3 - x2 * z3) / denom;
                Scalar C =
                        (-x2 * y1 + x3 * y1 + x1 * y2 - x3 * y2 - x1 * y3 + x2 * y3) / denom;
                base4 = -1;
                Scalar best_distance = std::numeric_limits<Scalar>::max();
                // Go over all points in P.
                const Scalar too_small = std::pow(max_base_diameter_ * kBaseTooSmall, 2);
                for (unsigned int i = 0; i < sampled_P_3D_.size(); ++i) {
                    if ((sampled_P_3D_[i].pos()- sampled_P_3D_[base1].pos()).squaredNorm() >= too_small &&
                        (sampled_P_3D_[i].pos()- sampled_P_3D_[base2].pos()).squaredNorm() >= too_small &&
                        (sampled_P_3D_[i].pos()- sampled_P_3D_[base3].pos()).squaredNorm() >= too_small) {
                        // Not too close to any of the first 3.
                        const Scalar distance =
                                std::abs(A * sampled_P_3D_[i].x() + B * sampled_P_3D_[i].y() +
                                         C * sampled_P_3D_[i].z() - 1.0);
                        // Search for the most planar.
                        if (distance < best_distance) {
                            best_distance = distance;
                            base4 = int(i);
                        }
                    }
                }
                // If we have a good one we can quit.
                if (base4 != -1) {
                    base_3D_[3] = sampled_P_3D_[base4];
                    if(TryQuadrilateral(invariant1, invariant2, base1, base2, base3, base4))
                        return true;
                }
            }
            current_trial++;
        }

        // We failed to find good enough base..
        return false;
    }

    template <typename Functor>
    // Initialize all internal data structures and data members.
    inline void Match4pcsBase<Functor>::Initialize(const std::vector<Point3D>& P,
                           const std::vector<Point3D>& Q) {
        fun_.Initialize(P,Q);
    }


    template <typename Functor>
    inline bool Match4pcsBase<Functor>::generateCongruents (Base& base, Set& congruent_quads) {

        Scalar invariant1, invariant2;
//#define STATIC_BASE

#ifdef STATIC_BASE
        static bool first_time = true;

  if (first_time){
      base[0] = 0;
      base[1] = 3;
      base[2] = 1;
      base[3] = 4;

      base_3D_[0] = sampled_P_3D_ [base[0]];
      base_3D_[1] = sampled_P_3D_ [base[1]];
      base_3D_[2] = sampled_P_3D_ [base[2]];
      base_3D_[3] = sampled_P_3D_ [base[3]];
      TryQuadrilateral(&invariant1, &invariant2, base[0], base[1], base[2], base[3]);

      first_time = false;
  }
  else
      return false;

#else
        if (!SelectQuadrilateral(invariant1, invariant2, base[0], base[1],
                                 base[2], base[3])) {
            return false;
        }
#endif

        // Computes distance between pairs.
        const Scalar distance1 = (base_3D_[0].pos()- base_3D_[1].pos()).norm();
        const Scalar distance2 = (base_3D_[2].pos()- base_3D_[3].pos()).norm();

        std::vector<std::pair<int, int>> pairs1, pairs2;

        // Compute normal angles.
        const Scalar normal_angle1 = (base_3D_[0].normal() - base_3D_[1].normal()).norm();
        const Scalar normal_angle2 = (base_3D_[2].normal() - base_3D_[3].normal()).norm();

        fun_.ExtractPairs(distance1, normal_angle1, distance_factor * options_.delta, 0, 1, &pairs1);
        fun_.ExtractPairs(distance2, normal_angle2, distance_factor * options_.delta, 2, 3, &pairs2);

//  Log<LogLevel::Verbose>( "Pair creation ouput: ", pairs1.size(), " - ", pairs2.size());

        if (pairs1.size() == 0 || pairs2.size() == 0) {
            return false;
        }

        if (!fun_.FindCongruentQuadrilaterals(invariant1, invariant2,
                                         distance_factor * options_.delta,
                                         distance_factor * options_.delta,
                                         pairs1,
                                         pairs2,
                                         &congruent_quads)) {
            return false;
        }

        return true;
    }
}

