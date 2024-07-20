#pragma once

#include <string>
#include <vector>

/** d-dimensional point in R^d */
class Point {

private:

    std::vector<double>             coordinates;    /** array of scalars representing the point */

public:

    /** construct a point from the given coordinates */
    Point(std::initializer_list<double> const & newCoordinates);

    /** create point owning the contents of provided vector */
    Point(std::vector<double> && newCoordinates);

    /** j-th coordinate */
    double operator[](int j) const;

};

/** extracts list of points from a text file
 *  returns the dimension of the points together with the list of points
 */
std::pair<int, std::vector<Point>> readPointsFromFile(std::string const & file);

/** rc(X, Y) problem instance */
struct RCInstance {

    int                             d;              /** dimension of the problem */

    std::vector<Point>              X;              /** set for which each point needs to be separated by all hyperplanes */
    std::vector<Point>              Y;              /** set for which each point needs to be separated by at least one hyperplane */

    RCInstance(
        int                         d,              /** dimension of the problem */
        std::vector<Point> const &  newX,           /** points defining X */
        std::vector<Point> const &  newY            /** points defining Y */
    );

    RCInstance(
        std::string const &         fileX,          /** file containing X */
        std::string const &         fileY           /** file containing Y */
    );

    /** prints problem size to std::cout and optionally enumerates all points */
    void print(bool enumeratePoints) const;

};