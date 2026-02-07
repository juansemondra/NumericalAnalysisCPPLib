#ifndef NUMERICALANALYSIS_H
#define NUMERICALANALYSIS_H

#include <map>
#include <iostream>

namespace NumericalAnalysis {

    class Polynomial{
        private:
        std::map<int, double> coeff;

        public:
        Polynomial();
        double evaluate(double value);
        double get(int degree);
        void update(int degree, double val);
        void add(int degree, double val);
        void extract_expression(std::string expression);
    };

    double bisection(Polynomial pol, double point_a, double point_b, int tolerance, int iterations);
    bool evaluate_tolerance(double value, int tolerance);

    

}

#endif