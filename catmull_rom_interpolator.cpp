#include "catmull_rom_interpolator.hpp"

#include <stdexcept>
#include <iostream>
#include <algorithm>

CatmullRomInterpolator::CatmullRomInterpolator(double duration, float tau)
    : points{}, polynomials{}, duration{duration}, accumulated_time{0.0}, tau{tau} {}

std::vector<float> CatmullRomInterpolator::get_accumulated_arc_lengths() {
    std::vector<float> arc_lengths;
    float arc_length = 0.0f;
    for (const auto &polynomial : polynomials) {
        arc_lengths.push_back(arc_length);
        arc_length += polynomial.arc_length;
    }
    arc_lengths.push_back(arc_length);
    return arc_lengths;
}

glm::vec3 CatmullRomInterpolator::interpolate(double delta_time) {
    if (points.size() < 4) {
        // TODO: or if you have overlapping adjacent points
        throw std::runtime_error("CatmullRomInterpolator needs at least 4 control points!");
    }

    accumulated_time = std::clamp(accumulated_time + delta_time, 0.0, duration);

    float total_arc_length = 0.0f;
    for (const auto &polynomial : polynomials) {
        total_arc_length += polynomial.arc_length;
    }

    float curr_arc_length = total_arc_length * accumulated_time / duration;

    float arc_length = 0.0f;
    int i = 0;
    while (arc_length + polynomials[i].arc_length < curr_arc_length) {
        arc_length += polynomials[i].arc_length;
        i++;
    }

    float t = 0.0f;
    if (delta_time >= 0.0 && i == cache.i) {
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

void CatmullRomInterpolator::reset() {
    accumulated_time = 0.0;
    cache.i = -1;
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
        polynomials.erase(polynomials.begin() + j);
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
