//
// Created by Sandra Alfaro on 26/04/18.
//

#ifndef OPENGR_FUNCTORFEATUREPOINTTEST_H
#define OPENGR_FUNCTORFEATUREPOINTTEST_H

#include <vector>

namespace GlobalRegistration {
    struct FilterTests {
    public :
        using TypeBase = std::vector<Point3D>;
        using Scalar      = typename Point3D::Scalar;
        using PairsVector = std::vector< std::pair<int, int> >;
        using VectorType  = typename Point3D::VectorType;
        using OptionType  = Match4PCSOptions;

    private :
        OptionType myOptions_;
        TypeBase myBase_3D_;

    public :
        inline FilterTests (OptionType options, TypeBase base) : myOptions_(options), myBase_3D_(base) {}

        inline std::pair<bool,bool> operator() (const Point3D& p, const Point3D& q, Scalar pair_normals_angle, int base_point1, int base_point2) {
            std::pair<bool,bool> res;
            res.first = false;
            res.second = false;

            VectorType segment1 = (myBase_3D_[base_point2].pos() -
                    myBase_3D_[base_point1].pos()).normalized();

            if ( myOptions_.max_normal_difference > 0 &&
                 q.normal().squaredNorm() > 0 &&
                 p.normal().squaredNorm() > 0) {
                const Scalar norm_threshold =
                        0.5 * myOptions_.max_normal_difference * M_PI / 180.0;
                const double first_normal_angle = (q.normal() - p.normal()).norm();
                const double second_normal_angle = (q.normal() + p.normal()).norm();
                // Take the smaller normal distance.
                const Scalar first_norm_distance =
                        std::min(std::abs(first_normal_angle - pair_normals_angle),
                                 std::abs(second_normal_angle - pair_normals_angle));
                // Verify appropriate angle between normals and distance.

                if (first_norm_distance > norm_threshold) return res;
            }
            // Verify restriction on the rotation angle, translation and colors.
            if (myOptions_.max_color_distance > 0) {
                const bool use_rgb = (p.rgb()[0] >= 0 && q.rgb()[0] >= 0 &&
                                      myBase_3D_[base_point1].rgb()[0] >= 0 &&
                                      myBase_3D_[base_point2].rgb()[0] >= 0);
                bool color_good = (p.rgb() - myBase_3D_[base_point1].rgb()).norm() <
                                  myOptions_.max_color_distance &&
                                  (q.rgb() - myBase_3D_[base_point2].rgb()).norm() <
                                  myOptions_.max_color_distance;

                if (use_rgb && ! color_good) return res;
            }

            if (myOptions_.max_translation_distance > 0) {
                const bool dist_good = (p.pos() - myBase_3D_[base_point1].pos()).norm() <
                                       myOptions_.max_translation_distance &&
                                       (q.pos() - myBase_3D_[base_point2].pos()).norm() <
                                       myOptions_.max_translation_distance;
                if (! dist_good) return res;
            }

            // need cleaning here
            if (myOptions_.max_angle > 0){
                VectorType segment2 = (q.pos() - p.pos()).normalized();
                if (std::acos(segment1.dot(segment2)) <= myOptions_.max_angle * M_PI / 180.0) {
                    res.second = true;
                }

                if (std::acos(segment1.dot(- segment2)) <= myOptions_.max_angle * M_PI / 180.0) {
                    // Add ordered pair.
                    res.first = true;
                }
            }else {
                res.first = true;
                res.second = true;
            }
            return res;
        }
    };
}

#endif //OPENGR_FUNCTORFEATUREPOINTTEST_H
