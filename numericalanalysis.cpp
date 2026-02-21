#include "numericalanalysis.h"
#include <cmath>
#include <regex>
#include <string>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <stdexcept>

namespace NumericalAnalysis
{

    bool Function::valid_key(const std::string &key) const
    {
        if (key == "sin" || key == "cos" || key == "tan")
            return true;
        std::regex poly_key_re(R"(x\^\d+)");
        return std::regex_match(key, poly_key_re);
    }

    /*
     *  "x^N"  →  x raised to the power N
     *  "sin"  →  sin(x)
     *  "cos"  →  cos(x)
     *  "tan"  →  tan(x)
     */

    double Function::evaluate_key(const std::string &key, double x) const
    {
        if (key == "sin")
            return std::sin(x);
        if (key == "cos")
            return std::cos(x);
        if (key == "tan")
            return std::tan(x);

        if (key.size() > 2 && key[0] == 'x' && key[1] == '^')
        {
            int degree = std::stoi(key.substr(2));
            return std::pow(x, degree);
        }

        std::cerr << "[Function::evaluate_key] Unknown key: " << key << "\n";
        return 0.0;
    }

    Function::Function() {}

    double Function::evaluate(double x) const
    {
        double result = 0.0;
        for (const auto &[key, coef] : coeff)
            result += static_cast<double>(coef) * evaluate_key(key, x);
        return result;
    }

    float Function::get(const std::string &key) const
    {
        auto it = coeff.find(key);
        return (it != coeff.end()) ? it->second : 0.0f;
    }

    void Function::update(const std::string &key, float val)
    {
        if (!valid_key(key))
        {
            std::cerr << "[Function::update] Invalid key: " << key << "\n";
            return;
        }
        auto it = coeff.find(key);
        if (it != coeff.end())
            it->second = val;
        else
            std::cout << "[Function::update] Term \"" << key
                      << "\" does not exist; use add() instead.\n";
    }

    void Function::add(const std::string &key, float val)
    {
        if (!valid_key(key))
        {
            std::cerr << "[Function::add] Invalid key: " << key << "\n";
            return;
        }
        if (coeff.count(key))
            std::cout << "[Function::add] Term \"" << key
                      << "\" already exists; use update() to change it.\n";
        else
            coeff.emplace(key, val);
    }

    void Function::extract_expression(const std::string &expression)
    {

        std::string expr = expression;
        expr.erase(std::remove_if(expr.begin(), expr.end(), ::isspace), expr.end());
        std::regex token_re(R"([+-]?[^+-]+)");
        std::regex trig_re(R"(^([+-]?)(\d+\.?\d*|\.?\d+)?(sin|cos|tan)$)");
        std::regex poly_re(R"(^([+-]?)(\d+\.?\d*|\.?\d+)?(x?)(?:\^(\d+))?$)");

        for (std::sregex_iterator it(expr.begin(), expr.end(), token_re), end;
             it != end; ++it)
        {
            std::string term = it->str();
            std::smatch m;

            if (std::regex_match(term, m, trig_re))
            {
                std::string sign = m[1].str();
                std::string digits = m[2].str();
                std::string func = m[3].str();

                float value = digits.empty() ? 1.0f : std::stof(digits);
                if (sign == "-")
                    value = -value;

                if (coeff.count(func))
                    coeff[func] += value;
                else
                    coeff.emplace(func, value);

                continue;
            }

            if (std::regex_match(term, m, poly_re))
            {
                std::string sign = m[1].str();
                std::string digits = m[2].str();
                bool has_x = !m[3].str().empty();
                std::string exp_str = m[4].str();

                if (!has_x && digits.empty())
                    continue;

                float value = digits.empty() ? 1.0f : std::stof(digits);
                if (sign == "-")
                    value = -value;

                int degree = has_x ? (exp_str.empty() ? 1 : std::stoi(exp_str)) : 0;
                std::string key = "x^" + std::to_string(degree);

                if (coeff.count(key))
                    coeff[key] += value;
                else
                    coeff.emplace(key, value);

                continue;
            }

            std::cerr << "[Function::extract_expression] Unrecognised term: \""
                      << term << "\"\n";
        }
    }

    void Function::print() const
    {
        std::cout << "f(x) =";
        bool first = true;
        for (const auto &[key, coef] : coeff)
        {
            if (first)
            {
                std::cout << " ";
                first = false;
            }
            else
            {
                std::cout << (coef >= 0 ? " + " : " - ");
            }

            float display_coef = (coef < 0 && !first) ? -coef : coef;
            if (first && coef < 0)
            {
                std::cout << "-";
                display_coef = -coef;
            }

            // Print coefficient only when it is not exactly 1 (or -1 for sign)
            if (std::fabs(display_coef - 1.0f) > 1e-6f || key == "x^0")
                std::cout << display_coef;

            if (key == "x^0")
            { /* constant, already printed */
            }
            else
            {
                std::cout << key;
            }
        }
        std::cout << "\n";
    }

    bool evaluate_tolerance(double xn, double xnp1, double tolerance)
    {
        double result = (xnp1 - xn) / xnp1;
        result = std::abs(result);
        if (result < tolerance) return true;
        return false;
    }

    double bisection(Function func, double point_a, double point_b, double tolerance, int iterations)
    {
        double fa = func.evaluate(point_a);
        double fb = func.evaluate(point_b);
        double p;
        double fp;
        for (int i = 0; i < iterations; i++){
            p = point_a + ((point_b - point_a) / 2);
            fp = func.evaluate(p);
            if (fp == 0 || ((point_b - point_a) / 2) < tolerance) return p;
            if ( (fa * fp) > 0 ) point_a = p;
            else point_b = p;
        }
        return -1;
    }

    double fixed_point(Function func, double initial_point, double tolerance, int iterations)
    {

        return -1;
    }

    double fake_position(Function func, double point_a, double point_b, double tolerance, int iterations)
    {

        return -1;
    }

    double newton_raphson(Function func, double initial_point, double tolerance, int iterations){

        return -1;
    }

    double secant_method(Function func, double point_a, double point_b, double tolerance, int iterations)
    {

        return -1;
    }