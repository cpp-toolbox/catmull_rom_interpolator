#ifndef CATMULL_ROM_INTERPOLATOR
#define CATMULL_ROM_INTERPOLATOR

#include <glm/glm.hpp> // TODO: make generic
#include <vector>

class CatmullRomInterpolator {
  public:
    // TODO: needs to take in duration (seconds) instead
    CatmullRomInterpolator(std::vector<glm::vec3> points, double ms_start_time, double ms_end_time, float tau);

    // TODO: needs to take in a delta time
    glm::vec3 interpolate(double ms_curr_time);

    // TODO: we need a function to update the duration and tension value

    glm::vec3 get_point(int i) const;
    void insert_point(int i, glm::vec3 point);
    void delete_point(int i);
    void update_point(int i, glm::vec3 point);

    int get_num_points() const;

  private:
    std::vector<glm::vec3> points;
    std::vector<glm::vec3> coefs_constant;
    std::vector<glm::vec3> coefs_linear;
    std::vector<glm::vec3> coefs_quadratic;
    std::vector<glm::vec3> coefs_cubic;
    std::vector<float> cummulative_arc_lengths; // arc length traversed until i-th control point
    double ms_start_time;
    double ms_end_time;
    float tau;

    struct Cache {
        double ms_curr_time = 0.0;
        int i = 0;
        float t = 0.0f;
        float cummulative_arc_length = 0.0f;
    } cache; // use cached variables to numerically integrate for arc length

    float dt = 0.001f;

    void compute_coefs();
    void recompute_coefs(int i); // i is the index of the inserted/deleted/updated point
};

#endif // CATMULL_ROM_INTERPOLATOR
