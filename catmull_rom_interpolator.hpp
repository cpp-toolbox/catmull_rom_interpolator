#ifndef CATMULL_ROM_INTERPOLATOR
#define CATMULL_ROM_INTERPOLATOR

#include <glm/glm.hpp> // TODO: make generic
#include <vector>

class CatmullRomInterpolator {
  public:
    CatmullRomInterpolator(double duration, float tau);

    glm::vec3 interpolate(double delta_time);
    void reset();

    int get_num_points() const;
    glm::vec3 get_point(int i) const;
    void set_point(int i, glm::vec3 coords);
    void insert_point(int i, glm::vec3 coords);
    void append_point(glm::vec3 coords);
    void delete_point(int i);

    double get_duration();
    void set_duration(double duration);

    float get_tension();
    void set_tension(float tau);

    std::vector<float> get_accumulated_arc_lengths();

    // bad - make private
    double accumulated_time;

  private:
    struct Polynomial {
        glm::vec3 coef_constant;
        glm::vec3 coef_linear;
        glm::vec3 coef_quadratic;
        glm::vec3 coef_cubic;
        float arc_length;
    };

    struct Point {
        glm::vec3 coords;
    };

    std::vector<Point> points;
    std::vector<Polynomial> polynomials;
    double duration;

    float tau;

    struct Cache {
        int i = 0;
        float t = 0.0f;
        float arc_length = 0.0f;
    } cache; // use cached variables to numerically integrate for arc length

    float dt = 0.001f;

    void init_polynomials(int start, int end);
};

#endif // CATMULL_ROM_INTERPOLATOR
