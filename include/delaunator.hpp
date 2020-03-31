#pragma once

#ifdef DELAUNATOR_HEADER_ONLY
#define INLINE inline
#else
#define INLINE
#endif

#include <limits>
#include <vector>

namespace delaunator {

constexpr std::size_t INVALID_INDEX =
    (std::numeric_limits<std::size_t>::max)();

class Point
{
public:
    Point(double x, double y) : m_x(x), m_y(y)
    {}
    Point() : m_x(0), m_y(0)
    {}


    double x() const
    { return m_x; }

    double y() const
    { return m_y; }

private:
    double m_x;
    double m_y;
};

class Delaunator {

public:
    std::vector<double> const& coords;
    // 'triangles' stores the indices to the 'X's of the input
    // 'coords'.
    std::vector<std::size_t> triangles;

    // 'halfedges' store indices into 'triangles'.  If halfedges[X] = Y,
    // It says that there's an edge from X to Y where a) X and Y are
    // both indices into triangles and b) X and Y are indices into different
    // triangles in the array.  This allows you to get from a triangle to
    // its adjacent triangle.  If the a triangle edge has no adjacent triangle,
    //  it half edge will be INVALID_INDEX.
    std::vector<std::size_t> halfedges;

    std::vector<std::size_t> hull_prev;
    std::vector<std::size_t> hull_next;

    // This contains indexes into the triangles array.
    std::vector<std::size_t> hull_tri;
    std::size_t hull_start;

    INLINE Delaunator(std::vector<double> const& in_coords);
    INLINE double get_hull_area();

private:
    std::vector<std::size_t> m_hash;
    Point m_center;
    std::size_t m_hash_size;
    std::vector<std::size_t> m_edge_stack;

    INLINE std::size_t legalize(std::size_t a);
    INLINE std::size_t hash_key(double x, double y) const;
    INLINE std::size_t add_triangle(
        std::size_t i0,
        std::size_t i1,
        std::size_t i2,
        std::size_t a,
        std::size_t b,
        std::size_t c);
    INLINE void link(std::size_t a, std::size_t b);
};

} //namespace delaunator

#undef INLINE
