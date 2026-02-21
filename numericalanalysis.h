#ifndef NUMERICALANALYSIS_H
#define NUMERICALANALYSIS_H

#include <map>
#include <string>
#include <iostream>

namespace NumericalAnalysis {

    class Function {
    private:
        std::map<std::string, double> coeff;
        bool valid_key(const std::string& key) const;
        double evaluate_key(const std::string& key, double x) const;

    public:
        Function();
        double evaluate(double x) const;
        double derivate_evaluate(double x) const;
        float get(const std::string& key) const;
        void update(const std::string& key, float val);
        void add(const std::string& key, float val);
        void extract_expression(const std::string& expression);
        void print() const;
    };

    static double eval_arg(const std::string &arg, double x);
    static double eval_arg_deriv(const std::string &arg, double x);

    bool evaluate_tolerance(double xn, double xnp1, double tolerance);

    double bisection(Function func, double point_a, double point_b, double tolerance, int iterations);
    double fixed_point(Function func, double initial_point, double tolerance, int iterations);
    double fake_position(Function func, double point_a, double point_b, double tolerance, int iterations);
    double newton_raphson(Function func, double initial_point, double tolerance, int iterations);
    double secant_method(Function func, double point_a, double point_b, double tolerance, int iterations);
}

#endif 