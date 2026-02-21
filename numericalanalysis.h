#ifndef NUMERICALANALYSIS_H
#define NUMERICALANALYSIS_H

#include <map>
#include <string>
#include <iostream>

namespace NumericalAnalysis {

    class Function {
    private:
        std::map<std::string, float> coeff;
        bool valid_key(const std::string& key) const;
        double evaluate_key(const std::string& key, double x) const;

    public:
        Function();
        double evaluate(double x) const;
        float get(const std::string& key) const;
        void update(const std::string& key, float val);
        void add(const std::string& key, float val);
        void extract_expression(const std::string& expression);
        void print() const;
    };

    double bisection(Function func, double point_a, double point_b,
                     int tolerance, int iterations);

    bool evaluate_tolerance(double value, int tolerance);

}

#endif 