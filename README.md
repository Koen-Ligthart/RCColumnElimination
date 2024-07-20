# RCColumnElimination
This is the source code of the column elimination algorithm that is used for my master thesis "Computing The Relaxation Complexity for Finite Point Sets Using Column Elimination" based on the work "Graph coloring with decision diagrams" by W.-J. van Hoeve.

It computes the least number of facets needed for a polyhedron to contain a given finite point set X and exclude all the points of a finite point set Y.

This depends on:
- C++20 standard
- CLI11 (https://github.com/CLIUtils/CLI11)
- cudd-3.0.0 (Colorado University Decision Diagram library by F. Somenzi)
- relaxation-complexity-1.0 (https://github.com/christopherhojny/relaxation_complexity, used as "Efficient MIP techniques for computing the relaxation complexity" by G. Averkov, C. Hojny and M. Schymura) together with its dependencies
- scipoptsuite-9.0.0 (https://scipopt.org/)

The algorithm can be called through the command line and the input file format is the format used in "relaxation_complexity".