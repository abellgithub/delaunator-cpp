#include "delaunator_test_main.hpp"

#include <cmath>

#include "delaunator-header-only.hpp"

using namespace std;

struct Point2D {
    double x, y;

    Point2D()
        : x(std::numeric_limits<double>::quiet_NaN())
        , y(std::numeric_limits<double>::quiet_NaN())
    {}

    Point2D(double x, double y)
        : x(x)
        , y(y)
    {}
};

Point2D rotate(Point2D const& p_in, Point2D const& center,
    double const& angle_deg)
{
    const double pi        = 3.14159265358979323846;
    const double angle_rad = angle_deg * pi / 180.0;
    const double s         = sin(angle_rad);
    const double c         = cos(angle_rad);
    Point2D p;

    // translate point back to origin:
    p.x = p_in.x - center.x;
    p.y = p_in.y - center.y;

    // rotate point
    double tmp = p.x * c - p.y * s;
    p.y = p.x * s + p.y * c;
    p.x = tmp;

    // translate point back:
    p.x += center.x;
    p.y += center.y;

    return p;
}

double distance(Point2D const& p1, Point2D const& p2)
{
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) +
        (p1.y - p2.y) * (p1.y - p2.y));
}

bool approx(double v1, double v2)
{
    const double eps = 0.01;
    return fabs(v1 - v2) < eps;
}

TEST(External_Lib_Delaunator_Test, rotate_grid_6_8_8_207_0)
{
    const size_t  grid_size = 6;
    const Point2D center(8, 8);
    const double  angle = 207.0;

    std::vector<double> points;

    // Generate rotated point grid of 1-by-1 squares.
    for (size_t x = 0; x < grid_size; ++x) {
        for (size_t y = 0; y < grid_size; ++y) {
            Point2D p = rotate(Point2D(x, y), center, angle);
            points.push_back(p.x);
            points.push_back(p.y);
        }
    }

    delaunator::Delaunator delaunator(points);

    // Validate.
    // Assume:
    // (1) since the grid is 1-by-1 square, that each triangle has to
    //     have two edges with length 1 and one edge with length sqrt(2) -
    //     give or take epsilon.
    // (2) due to rounding errors during rotation the hull points of the
    //     grid are not in a numerically perfect line, thus degenerate
    //     triangles will appear on the outer edges. Their area might be
    //     rather small (< 0.01) compared to a regular grid triangle (== 0.5).

    for (size_t i = 0; i < delaunator.triangles.size(); i += 3)
    {
        Point2D p1(points[2 * delaunator.triangles[i]],
                   points[1 + 2 * delaunator.triangles[i]]);
        Point2D p2(points[2 * delaunator.triangles[i + 1]],
                   points[1 + 2 * delaunator.triangles[i + 1]]);
        Point2D p3(points[2 * delaunator.triangles[i + 2]],
                   points[1 + 2 * delaunator.triangles[i + 2]]);

        // Determine edge lengths and triangle area.
        double len_1 = distance(p1, p2);
        double len_2 = distance(p2, p3);
        double len_3 = distance(p3, p1);
        double area  = fabs((p1.x * (p2.y - p3.y) +
                             p2.x * (p3.y - p1.y) +
                             p3.x * (p1.y - p2.y)) / 2.0);

        // Check edge lengths according to (1).
        bool cond_1 = false;
        if (approx(len_1, sqrt(2.0)))
        {
            cond_1 = approx(len_2, 1.0) && approx(len_3, 1.0);
        }
        else if (approx(len_2, sqrt(2.0)))
        {
            cond_1 = approx(len_1, 1.0) && approx(len_3, 1.0);
        }
        else if (approx(len_3, sqrt(2.0)))
        {
            cond_1 = approx(len_1, 1.0) && approx(len_2, 1.0);
        }

        // Check triangle area according to (2).
        bool cond_2 = false;
        if (cond_1)
        {
            // Triangle seems regular half 1-b-1 square. Just to be sure,
            // check the area.
            cond_2 = approx(area, 0.5);
        }
        else {
            // Irregular triangle. Accept degenerate ones. (Assume they
            // are on the outer egde, but do not check.)
            cond_2 = area < 0.01;
        }

        EXPECT_TRUE(cond_2) << "triangle = "
            "[" << p1.x << ", " << p1.y << "]; "
            "[" << p2.x << ", " << p2.y << "]; "
            "[" << p3.x << ", " << p3.y << "]";
    }
}
