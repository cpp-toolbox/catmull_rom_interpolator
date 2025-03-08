#include "catmull_rom_interpolator.hpp"

#include <stdexcept>
#include <iostream>
#include <algorithm>

CatmullRomInterpolator::CatmullRomInterpolator(double duration, float tau)
    : points{}, polynomials{}, duration{duration}, accumulatedTime{0.0}, tau{tau} {}

glm::vec3 CatmullRomInterpolator::interpolate(double deltaTime) {
    if (points.size() < 4) {
        // TODO: or if you have overlapping adjacent points
        throw std::runtime_error("CatmullRomInterpolator needs at least 4 control points!");
    }

    accumulatedTime = std::clamp(accumulatedTime + deltaTime, 0.0, duration);

    float total_arc_length = 0.0f;
    for (const auto &polynomial : polynomials) {
        total_arc_length += polynomial.arc_length;
    }

    float curr_arc_length = total_arc_length * accumulatedTime / duration;

    float arc_length = 0.0f;
    int i = 0;
    while (arc_length + polynomials[i].arc_length < curr_arc_length) {
        arc_length += polynomials[i].arc_length;
        i++;
    }

    float t = 0.0f;
    if (deltaTime >= 0.0 && i == cache.i) {
        // use cached variables to numerically integrate for arc length, if possible
        t = cache.t;
        arc_length = cache.arc_length;
    }
    while (arc_length < curr_arc_length) {
        arc_length += glm::length(polynomials[i].coef_linear + 2.0f * polynomials[i].coef_quadratic * t +
                                  3.0f * polynomials[i].coef_cubic * t * t) *
                      dt;
        t += dt;
    }

    cache.i = i;
    cache.t = t;
    cache.arc_length = arc_length;

    return polynomials[i].coef_constant + polynomials[i].coef_linear * t + polynomials[i].coef_quadratic * t * t +
           polynomials[i].coef_cubic * t * t * t;
}

int CatmullRomInterpolator::get_num_points() const { return points.size(); }

glm::vec3 CatmullRomInterpolator::get_point(int i) const { return points[i].coords; }

void CatmullRomInterpolator::set_point(int i, glm::vec3 coords) {
    Point point;
    point.coords = coords;

    points[i] = point;

    if (points.size() >= 4) {
        init_polynomials(std::max(i - 3, 0), std::min(i + 1, (int)polynomials.size()));
    }

    cache.i = -1;
}

void CatmullRomInterpolator::insert_point(int i, glm::vec3 coords) {
    points.emplace(points.begin() + i, coords);

    if (points.size() >= 4) {
        int j = std::clamp(i - 1, 0, (int)polynomials.size());
        polynomials.emplace(polynomials.begin() + j);
        init_polynomials(std::max(j - 3, 0), std::min(j + 1, (int)polynomials.size()));
    }

    cache.i = -1;
}

void CatmullRomInterpolator::append_point(glm::vec3 coords) { insert_point(points.size(), coords); }

void CatmullRomInterpolator::delete_point(int i) {
    points.erase(points.begin() + i);

    if (points.size() >= 4) {
        int j = std::clamp(i - 1, 0, (int)polynomials.size());
        polynomials.emplace(polynomials.begin() + j);
        init_polynomials(std::max(j - 3, 0), std::min(j, (int)polynomials.size()));
    }

    cache.i = -1;
}

double CatmullRomInterpolator::get_duration() { return duration; }

void CatmullRomInterpolator::set_duration(double duration_) { duration = duration_; }

float CatmullRomInterpolator::get_tension() { return tau; }

void CatmullRomInterpolator::set_tension(float tau_) { tau = tau_; }

void CatmullRomInterpolator::init_polynomials(int start, int end) {
    for (int i = start; i < end; i++) {
        // for derivation, see https://www.cs.cmu.edu/~fp/courses/graphics/asst5/catmullRom.pdf
        polynomials[i].coef_constant = points[i + 1].coords;
        polynomials[i].coef_linear = -tau * points[i].coords + tau * points[i + 2].coords;
        polynomials[i].coef_quadratic = 2.0f * tau * points[i].coords + (tau - 3.0f) * points[i + 1].coords +
                                        (3.0f - 2.0f * tau) * points[i + 2].coords - tau * points[i + 3].coords;
        polynomials[i].coef_cubic = -tau * points[i].coords + (2.0f - tau) * points[i + 1].coords +
                                    (tau - 2.0f) * points[i + 2].coords + tau * points[i + 3].coords;

        polynomials[i].arc_length = 0.0f;

        float t = 0.0f;
        while (t < 1.0f) {
            polynomials[i].arc_length +=
                glm::length(polynomials[i].coef_linear + 2.0f * polynomials[i].coef_quadratic * t +
                            3.0f * polynomials[i].coef_cubic * t * t) *
                dt;

            t += dt;
        }
    }
}

// #include "catmull_rom_interpolator.hpp"

// #include <stdexcept>
// #include <iostream>
// #include <algorithm>

// CatmullRomInterpolator::CatmullRomInterpolator(std::vector<glm::vec3> points, double ms_start_time, double
// ms_end_time,
//                                                float tau)
//     : points{points}, coefs_constant{}, coefs_linear{}, coefs_quadratic{}, coefs_cubic{},
//     ms_start_time{ms_start_time},
//       ms_end_time{ms_end_time}, tau{tau} {
//     compute_coefs();
// }

// int CatmullRomInterpolator::get_num_points() { return points.size(); }

// glm::vec3 CatmullRomInterpolator::interpolate(double ms_curr_time) {
//     if (ms_curr_time <= ms_start_time) {
//         return points[1];
//     }

//     if (ms_curr_time >= ms_end_time) {
//         return points[points.size() - 2];
//     }

//     float total_arc_length = cummulative_arc_lengths[cummulative_arc_lengths.size() - 1];
//     float curr_arc_length = total_arc_length * (ms_curr_time - ms_start_time) / (ms_end_time - ms_start_time);

//     int i = 0;
//     while (cummulative_arc_lengths[i + 1] < curr_arc_length)
//         i++;

//     float t = 0.0f;
//     float cummulative_arc_length = cummulative_arc_lengths[i];
//     if (ms_curr_time >= cache.ms_curr_time && i == cache.i) {
//         // use cached variables to numerically integrate for arc length, if possible
//         t = cache.t;
//         cummulative_arc_length = cache.cummulative_arc_length;
//     }
//     while (cummulative_arc_length < curr_arc_length) {
//         cummulative_arc_length +=
//             glm::length(coefs_linear[i] + 2.0f * coefs_quadratic[i] * t + 3.0f * coefs_cubic[i] * t * t) * dt;
//         t += dt;
//     }

//     cache.ms_curr_time = ms_curr_time;
//     cache.i = i;
//     cache.t = t;
//     cache.cummulative_arc_length = cummulative_arc_length;

//     return coefs_constant[i] + coefs_linear[i] * t + coefs_quadratic[i] * t * t + coefs_cubic[i] * t * t * t;
// }

// void CatmullRomInterpolator::insert_point(int i, glm::vec3 point) {

//     points.insert(points.begin() + i, point);
//     int coefs_size = coefs_constant.size();
//     coefs_constant.insert(coefs_constant.begin() + std::clamp(i - 1, 0, coefs_size), glm::vec3(0.0));
//     coefs_linear.insert(coefs_linear.begin() + std::clamp(i - 1, 0, coefs_size), glm::vec3(0.0));
//     coefs_quadratic.insert(coefs_quadratic.begin() + std::clamp(i - 1, 0, coefs_size), glm::vec3(0.0));
//     coefs_cubic.insert(coefs_cubic.begin() + std::clamp(i - 1, 0, coefs_size), glm::vec3(0.0));
//     cummulative_arc_lengths.insert(cummulative_arc_lengths.begin() + std::clamp(i - 1, 0, coefs_size + 1), 0.0);

//     recompute_coefs(i);

//     cache.i = -1;
// }

// void CatmullRomInterpolator::delete_point(int i) {
//     points.erase(points.begin() + i);
//     int coefs_size = coefs_constant.size();
//     coefs_constant.erase(coefs_constant.begin() + std::clamp(i - 1, 0, coefs_size));
//     coefs_linear.erase(coefs_linear.begin() + std::clamp(i - 1, 0, coefs_size));
//     coefs_quadratic.erase(coefs_quadratic.begin() + std::clamp(i - 1, 0, coefs_size));
//     coefs_cubic.erase(coefs_cubic.begin() + std::clamp(i - 1, 0, coefs_size));
//     cummulative_arc_lengths.erase(cummulative_arc_lengths.begin() + std::clamp(i - 1, 0, coefs_size + 1));

//     recompute_coefs(i);

//     cache.i = -1;
// }

// void CatmullRomInterpolator::update_point(int i, glm::vec3 point) {
//     points[i] = point;

//     recompute_coefs(i);

//     cache.i = -1;
// }

// void CatmullRomInterpolator::compute_coefs() {
//     cummulative_arc_lengths.clear();
//     coefs_constant.clear();
//     coefs_linear.clear();
//     coefs_quadratic.clear();
//     coefs_cubic.clear();

//     float cummulative_arc_length = 0.0f;
//     cummulative_arc_lengths.push_back(cummulative_arc_length);

//     /**
//      *      points: 0 1 2 3 4 5
//      *       coefs: X 0 1 2 X X
//      * cummulative: X 0 1 2 3 X
//      */

//     for (int i = 1; i < points.size() - 2; i++) {
//         // for derivation, see https://www.cs.cmu.edu/~fp/courses/graphics/asst5/catmullRom.pdf
//         coefs_constant.push_back(points[i]);
//         coefs_linear.push_back(-tau * points[i - 1] + tau * points[i + 1]);
//         coefs_quadratic.push_back(2.0f * tau * points[i - 1] + (tau - 3.0f) * points[i] +
//                                   (3.0f - 2.0f * tau) * points[i + 1] - tau * points[i + 2]);
//         coefs_cubic.push_back(-tau * points[i - 1] + (2.0f - tau) * points[i] + (tau - 2.0f) * points[i + 1] +
//                               tau * points[i + 2]);

//         float t = 0.0f;
//         while (t < 1.0f) {
//             cummulative_arc_length += glm::length(coefs_linear[i - 1] + 2.0f * coefs_quadratic[i - 1] * t +
//                                                   3.0f * coefs_cubic[i - 1] * t * t) *
//                                       dt;

//             t += dt;
//         }

//         cummulative_arc_lengths.push_back(cummulative_arc_length);
//     }
// }

// void CatmullRomInterpolator::recompute_coefs(int i) {
//     if (points.size() < 4) {
//         // or if you have overlapping points
//         throw std::runtime_error("CatmullRomInterpolator needs at least 4 control points!");
//     }

//     // should be from point i - 2 to i + 2
//     /**
//      *      points: 0 1 2 3 4 5
//      *       coefs: X 0 1 2 X X
//      * cummulative: X 0 1 2 3 X
//      */
//     for (int j = i - 3; j < i + 1; j++) { // TODO: boundary

//         // for derivation, see https://www.cs.cmu.edu/~fp/courses/graphics/asst5/catmullRom.pdf
//         coefs_constant[j] = points[j + 1];
//         coefs_linear[j] = -tau * points[j] + tau * points[j + 2];
//         coefs_quadratic[j] = 2.0f * tau * points[j] + (tau - 3.0f) * points[j + 1] +
//                              (3.0f - 2.0f * tau) * points[j + 2] - tau * points[j + 3];
//         coefs_cubic[j] =
//             -tau * points[j] + (2.0f - tau) * points[j + 1] + (tau - 2.0f) * points[j + 2] + tau * points[j + 3];
//     }

//     // update cummulative arc lengths
//     float copy = cummulative_arc_lengths[i + 1];
//     float cummulative_arc_length = cummulative_arc_lengths[i - 3];
//     for (int j = i - 2; j < i + 2; j++) {
//         float t = 0.0f;
//         while (t < 1.0f) {
//             cummulative_arc_length += glm::length(coefs_linear[j - 1] + 2.0f * coefs_quadratic[j - 1] * t +
//                                                   3.0f * coefs_cubic[j - 1] * t * t) *
//                                       dt;

//             t += dt;
//         }

//         cummulative_arc_lengths[j] = cummulative_arc_length;
//     }

//     float change_in_arc_length = cummulative_arc_length - copy;
//     for (int j = i + 2; j < points.size() - 2; j++) {
//         cummulative_arc_lengths[j] += change_in_arc_length;
//     }
// }
