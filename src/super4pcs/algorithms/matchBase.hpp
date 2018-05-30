//
// Created by Sandra Alfaro on 24/05/18.
//

#include <vector>
#include <atomic>
#include <chrono>

#ifdef SUPER4PCS_USE_OPENMP
#include <omp.h>
#endif

#include "super4pcs/shared.h"
#include "super4pcs/sampling.h"
#include "super4pcs/accelerators/kdtree.h"
#include "super4pcs/utils/logger.h"

#ifdef TEST_GLOBAL_TIMINGS
#   include "super4pcs/utils/timer.h"
#endif


namespace gr {
    template <typename Traits>
    MatchBase<Traits>::MatchBase(  const MatchOptions& options
            , const Utils::Logger& logger
    )
            :number_of_trials_(0)
            , max_base_diameter_(-1)
            , P_mean_distance_(1.0)
            , best_LCP_(0.0)
            , options_(options)
            , randomGenerator_(options.randomSeed)
            , logger_(logger)
#ifdef SUPER4PCS_USE_OPENMP
    , omp_nthread_congruent_(1)
#endif
    {}

    template <typename Traits>
    MatchBase<Traits>::~MatchBase(){}


// The main 4PCS function. Computes the best rigid transformation and transfoms
// Q toward P by this transformation
 template <typename Traits>
template <typename Sampler, typename Visitor>
typename MatchBase<Traits>::Scalar
MatchBase<Traits>::ComputeTransformation(const std::vector<Point3D>& P,
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
 template <typename Traits>
template <typename Visitor>
bool
MatchBase<Traits>::Perform_N_steps(int n,
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


    template <typename Traits>
    typename MatchBase<Traits>::Scalar
    MatchBase<Traits>::MeanDistance() {
        const Scalar kDiameterFraction = 0.2;
        using RangeQuery = gr::KdTree<Scalar>::RangeQuery<>;

        int number_of_samples = 0;
        Scalar distance = 0.0;

        for (size_t i = 0; i < sampled_P_3D_.size(); ++i) {

            RangeQuery query;
            query.sqdist = P_diameter_ * kDiameterFraction;
            query.queryPoint = sampled_P_3D_[i].pos().template cast<Scalar>();

            auto resId = kd_tree_.doQueryRestrictedClosestIndex(query , i).first;

            if (resId != gr::KdTree<Scalar>::invalidIndex()) {
                distance += (sampled_P_3D_[i].pos() - sampled_P_3D_[resId].pos()).norm();
                number_of_samples++;
            }
        }

        return distance / number_of_samples;
    }

    template <typename Traits>
    bool MatchBase<Traits>::SelectRandomTriangle(int &base1, int &base2, int &base3) {
        int number_of_points = sampled_P_3D_.size();
        base1 = base2 = base3 = -1;

        // Pick the first point at random.
        int first_point = randomGenerator_() % number_of_points;

        const Scalar sq_max_base_diameter_ = max_base_diameter_*max_base_diameter_;

        // Try fixed number of times retaining the best other two.
        Scalar best_wide = 0.0;
        for (int i = 0; i < kNumberOfDiameterTrials; ++i) {
            // Pick and compute
            const int second_point = randomGenerator_() % number_of_points;
            const int third_point = randomGenerator_() % number_of_points;
            const VectorType u =
                    sampled_P_3D_[second_point].pos() -
                    sampled_P_3D_[first_point].pos();
            const VectorType w =
                    sampled_P_3D_[third_point].pos() -
                    sampled_P_3D_[first_point].pos();
            // We try to have wide triangles but still not too large.
            Scalar how_wide = (u.cross(w)).norm();
            if (how_wide > best_wide &&
                u.squaredNorm() < sq_max_base_diameter_ &&
                w.squaredNorm() < sq_max_base_diameter_) {
                best_wide = how_wide;
                base1 = first_point;
                base2 = second_point;
                base3 = third_point;
            }
        }
        return base1 != -1 && base2 != -1 && base3 != -1;
    }

    template <typename Traits>
    void MatchBase<Traits>::initKdTree(){
        size_t number_of_points = sampled_P_3D_.size();

        // Build the kdtree.
        kd_tree_ = gr::KdTree<Scalar>(number_of_points);

        for (size_t i = 0; i < number_of_points; ++i) {
            kd_tree_.add(sampled_P_3D_[i].pos());
        }
        kd_tree_.finalize();
    }

    template <typename Traits>
    bool MatchBase<Traits>::ComputeRigidTransformation(
            const Coordinates& ref,
            const Coordinates& candidate,
            const Eigen::Matrix<Scalar, 3, 1>& centroid1,
            Eigen::Matrix<Scalar, 3, 1> centroid2,
            Scalar max_angle,
            Eigen::Ref<MatrixType> transform,
            Scalar& rms_,
            bool computeScale ) const {

        rms_ = kLargeNumber;

        Scalar kSmallNumber = 1e-6;

        // We only use the first 3 pairs. This simplifies the process considerably
        // because it is the planar case.

        const VectorType& p0 = ref[0].pos();
        const VectorType& p1 = ref[1].pos();
        const VectorType& p2 = ref[2].pos();
        VectorType  q0 = candidate[0].pos();
        VectorType  q1 = candidate[1].pos();
        VectorType  q2 = candidate[2].pos();

        Scalar scaleEst (1.);

        // Compute scale factor if needed
        if (computeScale){
            const VectorType& p3 = ref[3].pos();
            const VectorType& q3 = candidate[3].pos();

            const Scalar ratio1 = (p1 - p0).norm() / (q1 - q0).norm();
            const Scalar ratio2 = (p3 - p2).norm() / (q3 - q2).norm();

            const Scalar ratioDev  = std::abs(ratio1/ratio2 - Scalar(1.));  // deviation between the two
            const Scalar ratioMean = (ratio1+ratio2)/Scalar(2.);            // mean of the two

            if ( ratioDev > Scalar(0.1) )
                return kLargeNumber;


            //Log<LogLevel::Verbose>( ratio1, " ", ratio2, " ", ratioDev, " ", ratioMean);
            scaleEst = ratioMean;

            // apply scale factor to q
            q0 = q0*scaleEst;
            q1 = q1*scaleEst;
            q2 = q2*scaleEst;
            centroid2 *= scaleEst;
        }

        VectorType vector_p1 = p1 - p0;
        if (vector_p1.squaredNorm() == 0) return kLargeNumber;
        vector_p1.normalize();
        VectorType vector_p2 = (p2 - p0) - ((p2 - p0).dot(vector_p1)) * vector_p1;
        if (vector_p2.squaredNorm() == 0) return kLargeNumber;
        vector_p2.normalize();
        VectorType vector_p3 = vector_p1.cross(vector_p2);

        VectorType vector_q1 = q1 - q0;
        if (vector_q1.squaredNorm() == 0) return kLargeNumber;
        vector_q1.normalize();
        VectorType vector_q2 = (q2 - q0) - ((q2 - q0).dot(vector_q1)) * vector_q1;
        if (vector_q2.squaredNorm() == 0) return kLargeNumber;
        vector_q2.normalize();
        VectorType vector_q3 = vector_q1.cross(vector_q2);

        //cv::Mat rotation = cv::Mat::eye(3, 3, CV_64F);
        Eigen::Matrix<Scalar, 3, 3> rotation = Eigen::Matrix<Scalar, 3, 3>::Identity();

        Eigen::Matrix<Scalar, 3, 3> rotate_p;
        rotate_p.row(0) = vector_p1;
        rotate_p.row(1) = vector_p2;
        rotate_p.row(2) = vector_p3;

        Eigen::Matrix<Scalar, 3, 3> rotate_q;
        rotate_q.row(0) = vector_q1;
        rotate_q.row(1) = vector_q2;
        rotate_q.row(2) = vector_q3;

        rotation = rotate_p.transpose() * rotate_q;


        // Discard singular solutions. The rotation should be orthogonal.
        if (((rotation * rotation).diagonal().array() - Scalar(1) > kSmallNumber).any())
            return false;

        //FIXME
        if (max_angle >= 0) {
            // Discard too large solutions (todo: lazy evaluation during boolean computation
            if (! (
                    std::abs(std::atan2(rotation(2, 1), rotation(2, 2)))
                    <= max_angle &&

                    std::abs(std::atan2(-rotation(2, 0),
                                        std::sqrt(std::pow(rotation(2, 1),2) +
                                                  std::pow(rotation(2, 2),2))))
                    <= max_angle &&

                    std::abs(atan2(rotation(1, 0), rotation(0, 0)))
                    <= max_angle
            ))
                return false;
        }


        //FIXME
        // Compute rms and return it.
        rms_ = Scalar(0.0);
        {
            VectorType first, transformed;

            //cv::Mat first(3, 1, CV_64F), transformed;
            for (int i = 0; i < 3; ++i) {
                first = scaleEst*candidate[i].pos() - centroid2;
                transformed = rotation * first;
                rms_ += (transformed - ref[i].pos() + centroid1).norm();
            }
        }

        rms_ /= Scalar(ref.size());

        Eigen::Transform<Scalar, 3, Eigen::Affine> etrans (Eigen::Transform<Scalar, 3, Eigen::Affine>::Identity());
        transform = etrans
                .scale(scaleEst)
                .translate(centroid1)
                .rotate(rotation)
                .translate(-centroid2)
                .matrix();

        return true;
    }


// Verify a given transformation by computing the number of points in P at
// distance at most (normalized) delta from some point in Q. In the paper
// we describe randomized verification. We apply deterministic one here with
// early termination. It was found to be fast in practice.
    template <typename Traits>
    typename  MatchBase<Traits>::Scalar
    MatchBase<Traits>::Verify(const Eigen::Ref<const MatrixType> &mat) const {
        using RangeQuery = gr::KdTree<Scalar>::RangeQuery<>;

#ifdef TEST_GLOBAL_TIMINGS
        Timer t_verify (true);
#endif

        // We allow factor 2 scaling in the normalization.
        const Scalar epsilon = options_.delta;
#ifdef OPENGR_USE_WEIGHTED_LCP
        std::atomic<float> good_points(0);

        auto kernel = [](Scalar x) {
            return std::pow(std::pow(x,4) - Scalar(1), 2);
        };

        auto computeWeight = [kernel](Scalar sqx, Scalar th) {
            return kernel( std::sqrt(sqx) / th );
        };
#else
        std::atomic_uint good_points(0);
#endif
        const size_t number_of_points = sampled_Q_3D_.size();
        const size_t terminate_value = best_LCP_ * number_of_points;

        const Scalar sq_eps = epsilon*epsilon;
#ifdef OPENGR_USE_WEIGHTED_LCP
        const Scalar    eps = std::sqrt(sq_eps);
#endif

        for (size_t i = 0; i < number_of_points; ++i) {

            // Use the kdtree to get the nearest neighbor
#ifdef TEST_GLOBAL_TIMINGS
            Timer t (true);
#endif

            RangeQuery query;
            query.queryPoint = (mat * sampled_Q_3D_[i].pos().homogeneous()).template head<3>();
            query.sqdist     = sq_eps;

            auto result = kd_tree_.doQueryRestrictedClosestIndex( query );

#ifdef TEST_GLOBAL_TIMINGS
            kdTreeTime += Scalar(t.elapsed().count()) / Scalar(CLOCKS_PER_SEC);
#endif

            if ( result.first != gr::KdTree<Scalar>::invalidIndex() ) {
//      Point3D& q = sampled_P_3D_[near_neighbor_index[0]];
//      bool rgb_good =
//          (p.rgb()[0] >= 0 && q.rgb()[0] >= 0)
//              ? cv::norm(p.rgb() - q.rgb()) < options_.max_color_distance
//              : true;
//      bool norm_good = norm(p.normal()) > 0 && norm(q.normal()) > 0
//                           ? fabs(p.normal().ddot(q.normal())) >= cos_dist
//                           : true;
//      if (rgb_good && norm_good) {
#ifdef OPENGR_USE_WEIGHTED_LCP
                assert (result.second <= query.sqdist);
                good_points = good_points + computeWeight(result.second, eps);
#else
                good_points++;
#endif
//      }
            }

            // We can terminate if there is no longer chance to get better than the
            // current best LCP.
            if (number_of_points - i + good_points < terminate_value) {
                break;
            }
        }

#ifdef TEST_GLOBAL_TIMINGS
        verifyTime += Scalar(t_verify.elapsed().count()) / Scalar(CLOCKS_PER_SEC);
#endif
        return Scalar(good_points) / Scalar(number_of_points);
    }

    template <typename Traits>
    template <typename Sampler>
    void MatchBase<Traits>::init(const std::vector<Point3D>& P,
                         const std::vector<Point3D>& Q,
                         const Sampler& sampler){

#ifdef TEST_GLOBAL_TIMINGS
        kdTreeTime = 0;
    totalTime  = 0;
    verifyTime = 0;
#endif

        const Scalar kSmallError = 0.00001;
        const int kMinNumberOfTrials = 4;
        const Scalar kDiameterFraction = 0.3;

        centroid_P_ = VectorType::Zero();
        centroid_Q_ = VectorType::Zero();

        sampled_P_3D_.clear();
        sampled_Q_3D_.clear();

        // prepare P
        if (P.size() > options_.sample_size){
            sampler(P, options_, sampled_P_3D_);
        }
        else
        {
            Log<LogLevel::ErrorReport>( "(P) More samples requested than available: use whole cloud" );
            sampled_P_3D_ = P;
        }



        // prepare Q
        if (Q.size() > options_.sample_size){
            std::vector<Point3D> uniform_Q;
            sampler(Q, options_, uniform_Q);


            std::shuffle(uniform_Q.begin(), uniform_Q.end(), randomGenerator_);
            size_t nbSamples = std::min(uniform_Q.size(), options_.sample_size);
            auto endit = uniform_Q.begin(); std::advance(endit, nbSamples );
            std::copy(uniform_Q.begin(), endit, std::back_inserter(sampled_Q_3D_));
        }
        else
        {
            Log<LogLevel::ErrorReport>( "(Q) More samples requested than available: use whole cloud" );
            sampled_Q_3D_ = Q;
        }


        // center points around centroids
        auto centerPoints = [](std::vector<Point3D>&container,
                               VectorType& centroid){
            for(const auto& p : container) centroid += p.pos();
            centroid /= Scalar(container.size());
            for(auto& p : container) p.pos() -= centroid;
        };
        centerPoints(sampled_P_3D_, centroid_P_);
        centerPoints(sampled_Q_3D_, centroid_Q_);


        initKdTree();
        // Compute the diameter of P approximately (randomly). This is far from being
        // Guaranteed close to the diameter but gives good results for most common
        // objects if they are densely sampled.
        P_diameter_ = 0.0;
        for (int i = 0; i < kNumberOfDiameterTrials; ++i) {
            int at = randomGenerator_() % sampled_Q_3D_.size();
            int bt = randomGenerator_() % sampled_Q_3D_.size();

            Scalar l = (sampled_Q_3D_[bt].pos() - sampled_Q_3D_[at].pos()).norm();
            if (l > P_diameter_) {
                P_diameter_ = l;
            }
        }

        // Mean distance and a bit more... We increase the estimation to allow for
        // noise, wrong estimation and non-uniform sampling.
        P_mean_distance_ = MeanDistance();

        // Normalize the delta (See the paper) and the maximum base distance.
        // delta = P_mean_distance_ * delta;
        max_base_diameter_ = P_diameter_;  // * estimated_overlap_;

        // RANSAC probability and number of needed trials.
        Scalar first_estimation =
                std::log(kSmallError) / std::log(1.0 - pow(options_.getOverlapEstimation(),
                                                           static_cast<Scalar>(kMinNumberOfTrials)));
        // We use a simple heuristic to elevate the probability to a reasonable value
        // given that we don't simply sample from P, but instead, we bound the
        // distance between the points in the base as a fraction of the diameter.
        number_of_trials_ =
                static_cast<int>(first_estimation * (P_diameter_ / kDiameterFraction) /
                                 max_base_diameter_);
        if (number_of_trials_ < kMinNumberOfTrials)
            number_of_trials_ = kMinNumberOfTrials;

        Log<LogLevel::Verbose>( "norm_max_dist: ", options_.delta );
        current_trial_ = 0;
        best_LCP_ = 0.0;

        Q_copy_ = Q;
        for (int i = 0; i < 4; ++i) {
            base_[i] = 0;
            current_congruent_[i] = 0;
        }
        transform_ = Eigen::Matrix<Scalar, 4, 4>::Identity();

        // call Virtual handler

        Initialize(P,Q);

        best_LCP_ = Verify(transform_);
        Log<LogLevel::Verbose>( "Initial LCP: ", best_LCP_ );
    }

    template <typename Traits>
    template <typename Visitor>
    bool MatchBase<Traits>::TryOneBase(const Visitor &v) {
            Base base;
            Set congruent_quads;
            if (!generateCongruents(base,congruent_quads))
                return false;

            size_t nb = 0;

            bool match = TryCongruentSet(base,congruent_quads,v,nb);

            //if (nb != 0)
            //  Log<LogLevel::Verbose>( "Congruent quads: (", nb, ")    " );

            return match;
    }

    template <typename Traits>
    template <typename Visitor>
    bool MatchBase<Traits>::TryCongruentSet(Base& base, Set& set, Visitor &v,size_t &nbCongruent) {
        static const double pi = std::acos(-1);

        // get references to the basis coordinate
         Coordinates references;
        for (int i = 0; i!= Traits::size(); ++i)
            references[i] = sampled_P_3D_[base[i]];
        const Coordinates& ref = references;

        // Centroid of the basis, computed once and using only the three first points
        Eigen::Matrix<Scalar, 3, 1> centroid1 = (ref[0].pos() + ref[1].pos() + ref[2].pos()) / Scalar(3);


        std::atomic<size_t> nbCongruentAto(0);

#ifdef SUPER4PCS_USE_OPENMP
#pragma omp parallel for num_threads(omp_nthread_congruent_)
#endif
        Coordinates congruent_candidate;
        for (int i = 0; i < int(set.size()); ++i) {
            const Base& congruent_ids = set[i];
            for (int j = 0; j!= Traits::size(); ++j)
                congruent_candidate[j] = sampled_Q_3D_[congruent_ids[j]];

            Eigen::Matrix<Scalar, 4, 4> transform;

            // Centroid of the sets, computed in the loop using only the three first points
            Eigen::Matrix<Scalar, 3, 1> centroid2;


#ifdef STATIC_BASE
            Log<LogLevel::Verbose>( "Ids: ");
            for (int j = 0; j!= Traits::size(); ++j)
               Log<LogLevel::Verbose>( base[j], "\t");
            Log<LogLevel::Verbose>( "     ");
            for (int j = 0; j!= Traits::size(); ++j)
                Log<LogLevel::Verbose>( congruent_ids[j], "\t");
#endif

            centroid2 = (congruent_candidate[0].pos() +
                         congruent_candidate[1].pos() +
                         congruent_candidate[2].pos()) / Scalar(3.);

            Scalar rms = -1;

            const bool ok =
                    ComputeRigidTransformation(ref,     // input congruent quad
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
                        for (int j = 0; j!= Traits::size(); ++j)
                            base_[j] = base[j];


                        for (int j = 0; j!= Traits::size(); ++j)
                            current_congruent_[j] = congruent_ids[j];

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