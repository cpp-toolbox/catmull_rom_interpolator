#include "catmull_rom_interpolator.hpp"

#include <stdexcept>
#include <iostream>

CatmullRomInterpolator::CatmullRomInterpolator(std::vector<glm::vec3> points, double ms_start_time, double ms_end_time,
                                               float tau)
    : points{points}, coefs_constant{}, coefs_linear{}, coefs_quadratic{}, coefs_cubic{}, ms_start_time{ms_start_time},
      ms_end_time{ms_end_time}, tau{tau} {
    compute_coefs();
}

int CatmullRomInterpolator::get_num_points() { return points.size(); }

glm::vec3 CatmullRomInterpolator::interpolate(double ms_curr_time) {
    if (ms_curr_time <= ms_start_time) {
        return points[1];
    }

    if (ms_curr_time >= ms_end_time) {
        return points[points.size() - 2];
    }

    float total_arc_length = cummulative_arc_lengths[cummulative_arc_lengths.size() - 1];
    float curr_arc_length = total_arc_length * (ms_curr_time - ms_start_time) / (ms_end_time - ms_start_time);

    int i = 0;
    while (cummulative_arc_lengths[i + 1] < curr_arc_length)
        i++;

    float t = 0.0f;
    float cummulative_arc_length = cummulative_arc_lengths[i];
    if (ms_curr_time >= cache.ms_curr_time && i == cache.i) {
        // use cached variables to numerically integrate for arc length, if possible
        t = cache.t;
        cummulative_arc_length = cache.cummulative_arc_length;
    }
    while (cummulative_arc_length < curr_arc_length) {
        cummulative_arc_length +=
            glm::length(coefs_linear[i] + 2.0f * coefs_quadratic[i] * t + 3.0f * coefs_cubic[i] * t * t) * dt;
        t += dt;
    }

    cache.ms_curr_time = ms_curr_time;
    cache.i = i;
    cache.t = t;
    cache.cummulative_arc_length = cummulative_arc_length;

    return coefs_constant[i] + coefs_linear[i] * t + coefs_quadratic[i] * t * t + coefs_cubic[i] * t * t * t;
}

void CatmullRomInterpolator::insert_point(int i, glm::vec3 point) {
    points.insert(points.begin() + i, point);
    coefs_constant.insert(coefs_constant.begin() + i - 1, glm::vec3(0.0));
    coefs_linear.insert(coefs_linear.begin() + i - 1, glm::vec3(0.0));
    coefs_quadratic.insert(coefs_quadratic.begin() + i - 1, glm::vec3(0.0));
    coefs_cubic.insert(coefs_cubic.begin() + i - 1, glm::vec3(0.0));
    cummulative_arc_lengths.insert(cummulative_arc_lengths.begin() + i - 1, 0.0);

    recompute_coefs(i);

    cache.i = -1;
}

void CatmullRomInterpolator::delete_point(int i) {
    points.erase(points.begin() + i);
    coefs_constant.erase(coefs_constant.begin() + i - 1);
    coefs_linear.erase(coefs_linear.begin() + i - 1);
    coefs_quadratic.erase(coefs_quadratic.begin() + i - 1);
    coefs_cubic.erase(coefs_cubic.begin() + i - 1);
    cummulative_arc_lengths.erase(cummulative_arc_lengths.begin() + i - 1);

    recompute_coefs(i);

    cache.i = -1;
}

void CatmullRomInterpolator::update_point(int i, glm::vec3 point) {
    points[i] = point;

    recompute_coefs(i);

    cache.i = -1;
}

void CatmullRomInterpolator::compute_coefs() {
    if (points.size() < 4) {
        // or if you have overlapping points
        throw std::runtime_error("CatmullRomInterpolator needs at least 4 control points!");
    }

    cummulative_arc_lengths.clear();
    coefs_constant.clear();
    coefs_linear.clear();
    coefs_quadratic.clear();
    coefs_cubic.clear();

    float cummulative_arc_length = 0.0f;
    cummulative_arc_lengths.push_back(cummulative_arc_length);

    for (int i = 1; i < points.size() - 2; i++) {
        // for derivation, see https://www.cs.cmu.edu/~fp/courses/graphics/asst5/catmullRom.pdf
        coefs_constant.push_back(points[i]);
        coefs_linear.push_back(-tau * points[i - 1] + tau * points[i + 1]);
        coefs_quadratic.push_back(2.0f * tau * points[i - 1] + (tau - 3.0f) * points[i] +
                                  (3.0f - 2.0f * tau) * points[i + 1] - tau * points[i + 2]);
        coefs_cubic.push_back(-tau * points[i - 1] + (2.0f - tau) * points[i] + (tau - 2.0f) * points[i + 1] +
                              tau * points[i + 2]);

        float t = 0.0f;
        while (t < 1.0f) {
            cummulative_arc_length += glm::length(coefs_linear[i - 1] + 2.0f * coefs_quadratic[i - 1] * t +
                                                  3.0f * coefs_cubic[i - 1] * t * t) *
                                      dt;

            t += dt;
        }

        cummulative_arc_lengths.push_back(cummulative_arc_length);
    }
}

void CatmullRomInterpolator::recompute_coefs(int i) {
    if (points.size() < 4) {
        // or if you have overlapping points
        throw std::runtime_error("CatmullRomInterpolator needs at least 4 control points!");
    }

    for (int j = i - 3; j < i + 1; j++) { // TODO: boundary
        // for derivation, see https://www.cs.cmu.edu/~fp/courses/graphics/asst5/catmullRom.pdf
        coefs_constant[j] = points[j + 1];
        coefs_linear[j] = -tau * points[j] + tau * points[j + 2];
        coefs_quadratic[j] = 2.0f * tau * points[j] + (tau - 3.0f) * points[j + 1] +
                             (3.0f - 2.0f * tau) * points[j + 2] - tau * points[j + 3];
        coefs_cubic[j] =
            -tau * points[j] + (2.0f - tau) * points[j + 1] + (tau - 2.0f) * points[j + 2] + tau * points[j + 3];
    }

    // update cummulative arc lengths
    float copy = cummulative_arc_lengths[i + 1];
    float cummulative_arc_length = cummulative_arc_lengths[i - 3];
    for (int j = i - 2; j < i + 2; j++) {
        float t = 0.0f;
        while (t < 1.0f) {
            cummulative_arc_length += glm::length(coefs_linear[j - 1] + 2.0f * coefs_quadratic[j - 1] * t +
                                                  3.0f * coefs_cubic[j - 1] * t * t) *
                                      dt;

            t += dt;
        }

        cummulative_arc_lengths[j] = cummulative_arc_length;
    }

    float change_in_arc_length = cummulative_arc_length - copy;
    for (int j = i + 2; j < points.size() - 2; j++) {
        cummulative_arc_lengths[j] += change_in_arc_length;
    }
}
