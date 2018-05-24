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

#include "../shared4pcs.h"
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
            : Base(options,logger)
            , fun_(sampled_Q_3D_,base_3D_,options) {}

    template <typename Functor>
    Match4pcsBase<Functor>::~Match4pcsBase() {}

    template <typename Functor>
    // Initialize all internal data structures and data members.
    void Match4pcsBase<Functor>::Initialize(const std::vector<Point3D>& P,
                           const std::vector<Point3D>& Q) {
        fun_.Initialize(P,Q);
    }

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


    // The main 4PCS function. Computes the best rigid transformation and transfoms
    // Q toward P by this transformation
    template <typename Functor>
    template <typename Sampler, typename Visitor>
    typename Match4pcsBase<Functor>::Scalar
    Match4pcsBase<Functor>::ComputeTransformation(const std::vector<Point3D>& P,
                                     std::vector<Point3D>* Q,
                                     Eigen::Ref<MatrixType> transformation,
                                     const Sampler& sampler,
                                     const Visitor& v) {

        if (Q == nullptr) return kLargeNumber;
        if (P.empty() || Q->empty()) return kLargeNumber;

        init(P, *Q, sampler);

        if (best_LCP_ != Scalar(1.))
            Perform_N_steps(number_of_trials_, transformation, Q, v);

#ifdef TEST_GLOBAL_TIMINGS
        Log<LogLevel::Verbose>( "----------- Timings (msec) -------------" );
  Log<LogLevel::Verbose>( " Total computation time  : ", totalTime   );
  Log<LogLevel::Verbose>( " Total verify time       : ", verifyTime  );
  Log<LogLevel::Verbose>( "    Kdtree query         : ", kdTreeTime  );
  Log<LogLevel::Verbose>( "----------------------------------------" );
#endif

        return best_LCP_;
    }


    // Performs N RANSAC iterations and compute the best transformation. Also,
    // transforms the set Q by this optimal transformation.
    template <typename Functor>
    template <typename Visitor>
    bool Match4pcsBase<Functor>::Perform_N_steps(int n,
                                    Eigen::Ref<MatrixType> transformation,
                                    std::vector<Point3D>* Q,
                                    const Visitor &v) {
        using std::chrono::system_clock;
        if (Q == nullptr) return false;

#ifdef TEST_GLOBAL_TIMINGS
        Timer t (true);
#endif


        // The transformation has been computed between the two point clouds centered
        // at the origin, we need to recompute the translation to apply it to the original clouds
        auto getGlobalTransform = [this](Eigen::Ref<MatrixType> transformation){
            Eigen::Matrix<Scalar, 3, 3> rot, scale;
            Eigen::Transform<Scalar, 3, Eigen::Affine> (transform_).computeRotationScaling(&rot, &scale);
            transformation = transform_;
            transformation.col(3) = (qcentroid1_ + centroid_P_ - ( rot * scale * (qcentroid2_ + centroid_Q_))).homogeneous();
        };

        Scalar last_best_LCP = best_LCP_;
        v(0, best_LCP_, transformation);

        bool ok = false;
        std::chrono::time_point<system_clock> t0 = system_clock::now(), end;
        for (int i = current_trial_; i < current_trial_ + n; ++i) {
            ok = TryOneBase(v);

            Scalar fraction_try  = Scalar(i) / Scalar(number_of_trials_);
            Scalar fraction_time =
                    std::chrono::duration_cast<std::chrono::seconds>
                            (system_clock::now() - t0).count() /
                    options_.max_time_seconds;
            Scalar fraction = std::max(fraction_time, fraction_try);

            if (v.needsGlobalTransformation()) {
                getGlobalTransform(transformation);
            } else {
                transformation = transform_;
            }

            v(fraction, best_LCP_, transformation);

            // ok means that we already have the desired LCP.
            if (ok || i > number_of_trials_ || fraction >= 0.99 || best_LCP_ == 1.0) break;
        }

        current_trial_ += n;
        if (best_LCP_ > last_best_LCP) {
            *Q = Q_copy_;

            getGlobalTransform(transformation);

            // Transforms Q by the new transformation.
            for (size_t i = 0; i < Q->size(); ++i) {
                (*Q)[i].pos() = (transformation * (*Q)[i].pos().homogeneous()).head<3>();
            }
        }
#ifdef TEST_GLOBAL_TIMINGS
        totalTime += Scalar(t.elapsed().count()) / Scalar(CLOCKS_PER_SEC);
#endif

        return ok || current_trial_ >= number_of_trials_;
    }


    template <typename Functor>
    template <typename Visitor>
    bool Match4pcsBase<Functor>::TryOneBase(const Visitor &v) {
        Scalar invariant1, invariant2;
        int base_id1, base_id2, base_id3, base_id4;

//#define STATIC_BASE

#ifdef STATIC_BASE
        static bool first_time = true;

  if (first_time){
      base_id1 = 0;
      base_id2 = 3;
      base_id3 = 1;
      base_id4 = 4;

      base_3D_[0] = sampled_P_3D_ [base_id1];
      base_3D_[1] = sampled_P_3D_ [base_id2];
      base_3D_[2] = sampled_P_3D_ [base_id3];
      base_3D_[3] = sampled_P_3D_ [base_id4];
        // fun_. TryQuadrilateral(&invariant1, &invariant2, base_id1, base_id2, base_id3, base_id4);
      TryQuadrilateral(&invariant1, &invariant2, base_id1, base_id2, base_id3, base_id4);

      first_time = false;
  }
  else
      return false;

#else
        if (!SelectQuadrilateral(invariant1, invariant2, base_id1, base_id2,
                                 base_id3, base_id4)) {
            return false;
        }
#endif

        // Computes distance between pairs.
        const Scalar distance1 = (base_3D_[0].pos()- base_3D_[1].pos()).norm();
        const Scalar distance2 = (base_3D_[2].pos()- base_3D_[3].pos()).norm();

        std::vector<std::pair<int, int>> pairs1, pairs2;
        std::vector<Quadrilateral> congruent_quads;

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

        size_t nb = 0;

        bool match = TryCongruentSet(base_id1, base_id2, base_id3, base_id4,
                                     congruent_quads,
                                     v,
                                     nb);

        //if (nb != 0)
        //  Log<LogLevel::Verbose>( "Congruent quads: (", nb, ")    " );

        return match;
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
    template <typename Visitor>
    bool Match4pcsBase<Functor>::TryCongruentSet(int base_id1,
                                                 int base_id2,
                                                 int base_id3,
                                                 int base_id4,
                                                 const std::vector<Quadrilateral> &congruent_quads,
                                                 const Visitor &v,
                                                 size_t &nbCongruent) {
        static const double pi = std::acos(-1);

        // get references to the basis coordinates
        const Point3D& b1 = sampled_P_3D_[base_id1];
        const Point3D& b2 = sampled_P_3D_[base_id2];
        const Point3D& b3 = sampled_P_3D_[base_id3];
        const Point3D& b4 = sampled_P_3D_[base_id4];

        // set the basis coordinates in the congruent quad array
        const std::array<Point3D, 4> congruent_base {{b1, b2, b3, b4}};


        // Centroid of the basis, computed once and using only the three first points
        Eigen::Matrix<Scalar, 3, 1> centroid1 = (b1.pos() + b2.pos() + b3.pos()) / Scalar(3);


        std::atomic<size_t> nbCongruentAto(0);

#ifdef SUPER4PCS_USE_OPENMP
#pragma omp parallel for num_threads(omp_nthread_congruent_)
#endif
        for (int i = 0; i < int(congruent_quads.size()); ++i) {
            std::array<Point3D, 4> congruent_candidate;

            Eigen::Matrix<Scalar, 4, 4> transform;

            // Centroid of the sets, computed in the loop using only the three first points
            Eigen::Matrix<Scalar, 3, 1> centroid2;

            const int a = congruent_quads[i].vertices[0];
            const int b = congruent_quads[i].vertices[1];
            const int c = congruent_quads[i].vertices[2];
            const int d = congruent_quads[i].vertices[3];
            congruent_candidate[0] = sampled_Q_3D_[a];
            congruent_candidate[1] = sampled_Q_3D_[b];
            congruent_candidate[2] = sampled_Q_3D_[c];
            congruent_candidate[3] = sampled_Q_3D_[d];

#ifdef STATIC_BASE
            Log<LogLevel::Verbose>( "Ids: ", base_id1, "\t", base_id2, "\t", base_id3, "\t", base_id4);
      Log<LogLevel::Verbose>( "     ", a, "\t", b, "\t", c, "\t", d);
#endif

            centroid2 = (congruent_candidate[0].pos() +
                         congruent_candidate[1].pos() +
                         congruent_candidate[2].pos()) / Scalar(3.);

            Scalar rms = -1;

            const bool ok =
                    ComputeRigidTransformation(congruent_base,     // input congruent quad
                                               congruent_candidate,// tested congruent quad
                                               centroid1,          // input: basis centroid
                                               centroid2,          // input: candidate quad centroid
                                               options_.max_angle * pi / 180.0, // maximum per-dimension angle, check return value to detect invalid cases
                                               transform,          // output: transformation
                                               rms,                // output: rms error of the transformation between the basis and the congruent quad
#ifdef MULTISCALE
                            true
#else
                                               false
#endif
                    );             // state: compute scale ratio ?

            if (ok && rms >= Scalar(0.)) {

                // We give more tolerant in computing the best rigid transformation.
                if (rms < distance_factor * options_.delta) {

                    nbCongruentAto++;
                    // The transformation is computed from the point-clouds centered inn [0,0,0]

                    // Verify the rest of the points in Q against P.
                    Scalar lcp = Verify(transform);

                    // transformation has been computed between the two point clouds centered
                    // at the origin, we need to recompute the translation to apply it to the original clouds
                    auto getGlobalTransform =
                            [this, transform, centroid1, centroid2]
                                    (Eigen::Ref<MatrixType> transformation){
                                Eigen::Matrix<Scalar, 3, 3> rot, scale;
                                Eigen::Transform<Scalar, 3, Eigen::Affine> (transform).computeRotationScaling(&rot, &scale);
                                transformation = transform;
                                transformation.col(3) = (centroid1 + centroid_P_ - ( rot * scale * (centroid2 + centroid_Q_))).homogeneous();
                            };

                    if (v.needsGlobalTransformation())
                    {
                        Eigen::Matrix<Scalar, 4, 4> transformation = transform;
                        getGlobalTransform(transformation);
                        v(-1, lcp, transformation);
                    }
                    else
                        v(-1, lcp, transform);

#pragma omp critical
                    if (lcp > best_LCP_) {
                        // Retain the best LCP and transformation.
                        base_[0] = base_id1;
                        base_[1] = base_id2;
                        base_[2] = base_id3;
                        base_[3] = base_id4;

                        current_congruent_[0] = a;
                        current_congruent_[1] = b;
                        current_congruent_[2] = c;
                        current_congruent_[3] = d;

                        best_LCP_    = lcp;
                        transform_   = transform;
                        qcentroid1_  = centroid1;
                        qcentroid2_  = centroid2;
                    }
                    // Terminate if we have the desired LCP already.
                    if (best_LCP_ > options_.getTerminateThreshold()){
                        continue;
                    }
                }
            }
        }

        nbCongruent = nbCongruentAto;

        // If we reached here we do not have yet the desired LCP.
        return best_LCP_ > options_.getTerminateThreshold() /*false*/;
    }
}

