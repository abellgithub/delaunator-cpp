#include <delaunator.hpp>
#include <iostream>
#include <string>

int main() {
    using namespace delaunator;

    /* x0, y0, x1, y1, ... */
    std::vector<double> coords = {-3, 1, -1, 1, 1, 1, 3, 1,
                                  3, -1, 1, -1, -1, -1, -3, -1};

    //triangulation happens here
    Delaunator d(coords);

    auto xcoord = [&d](size_t tri, int coord) -> double
    {
        size_t index = d.triangles[tri + coord];
        return d.coords[2 * index];
    };

    auto ycoord = [&d](size_t tri, int coord) -> double
    {
        size_t index = d.triangles[tri + coord];
        return d.coords[2 * index + 1];
    };

    for (std::size_t i = 0; i < d.triangles.size(); i += 3) {
        std::cout << "Triangle " << (i / 3) << " points: [[" <<
            xcoord(i, 0) << ", " << ycoord(i, 0) << "], [" <<
            xcoord(i, 1) << ", " << ycoord(i, 1) << "], [" <<
            xcoord(i, 2) << ", " << ycoord(i, 2) << "]]\n";
    }

    auto adjTriangle = [&d](size_t edge) -> std::string
    {
        size_t adj = d.halfedges[edge];
        return adj == INVALID_INDEX ? "None" : std::to_string(adj / 3);
    };

    for (std::size_t i = 0; i < d.triangles.size(); i += 3)
    {
        std::cerr << "Adjacent triangles for triangle number " << (i / 3) << ": " <<
            adjTriangle(i) << ", " << adjTriangle(i + 1) << " and " << adjTriangle(i + 2) << "\n";
    }
}
