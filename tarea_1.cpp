#include <iostream>
#include "numericalanalysis.h"

int main (int argc, char** argv) {

    if(argc < 6)
    {
        std::cerr  << "Usage: " << argv[0] << " function" << " point_a" << " point_b" << " tolerance" << " iterations" << std::endl;
        return EXIT_FAILURE;
    }

    std::string function_str = argv[1];
    double point_a = std::stod(argv[2]);
    double point_b = std::stod(argv[3]);
    int tolerance = std::stoi(argv[4]);
    int iterations = std::stoi(argv[5]);
    NumericalAnalysis::Polynomial function;
    function.extract_expression(function_str);

    if (point_a > 0 && point_b > 0) {
    std::cerr << "For the bisectional method to work one point must be positive and one negative." << std::endl;
    return EXIT_FAILURE;
    }
    else if (point_a < 0 && point_b < 0) {
    std::cerr << "For the bisectional method to work one point must be positive and one negative." << std::endl;
    return EXIT_FAILURE;
    }
    if (point_a == 0 || point_b == 0) {
    std::cerr << "For the bisectional method to work one point must be positive and one negative." << std::endl;
    return EXIT_FAILURE;
    }

    double result = NumericalAnalysis::bisection(function, point_a, point_b, tolerance, iterations);

    std::cout << "Resultado del metodo de biseccion: " << result << std::endl;

    return EXIT_SUCCESS;
}