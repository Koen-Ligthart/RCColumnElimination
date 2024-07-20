#include <fstream>
#include <iostream>
#include "rc_instance.h"

Point::Point(std::initializer_list<double> const & newCoordinates) : coordinates { newCoordinates } {}

Point::Point(std::vector<double> && newCoordinates) : coordinates { newCoordinates } {}

double Point::operator[](int const j) const { return coordinates[j]; }

std::pair<int, std::vector<Point>> readPointsFromFile(std::string const & file) {
    std::ifstream in(file);
    if (!in)
        throw std::invalid_argument("file could not be found");
    int n, d;
    in >> n >> d;
    std::vector<Point> points;
    points.reserve(n);
    for (int i = 0; i < n; ++i) {
        std::vector<double> coordinates;
        coordinates.reserve(d);
        for (int j = 0; j < d; ++j) {
            double coordinate;
            in >> coordinate;
            coordinates.push_back(coordinate);
        }
        points.emplace_back(std::move(coordinates));
    }
    return {d, points};
}

RCInstance::RCInstance(int const newD, std::vector<Point> const & newX, std::vector<Point> const & newY) : d{ newD }, X{ newX }, Y{ newY } {}

RCInstance::RCInstance(std::string const & fileX, std::string const & fileY) {
    auto const [dX, newX] = readPointsFromFile(fileX);
    auto const [dY, newY] = readPointsFromFile(fileY);
    if (dX != dY)
        throw std::runtime_error("X and Y do not share same dimension");
    d = dX;
    X = newX;
    Y = newY;
}

void printPointSet(int const d, std::vector<Point> const & points) {
    for (Point const & point : points) {
        std::cout << "   ";
        for (int i = 0; i < d; ++i)
            std::cout << " " << point[i];
        std::cout << "\n";
    }
}

void RCInstance::print(bool const enumeratePoints) const {
    std::cout << "problem instance of size (|X| + |Y|) * d = (m + n) * d = (" << X.size() << " + " << Y.size() << ") * " << d << "\n";
    if (enumeratePoints) {
        std::cout << "point set X\n";
        printPointSet(d, X);
        std::cout << "point set Y\n";
        printPointSet(d, Y);
    }
    std::cout.flush();
}