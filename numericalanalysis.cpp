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

    static double eval_arg(const std::string &arg, double x)
    {
        if (arg.size() > 2 && arg[0] == 'x' && arg[1] == '^')
            return std::pow(x, std::stoi(arg.substr(2)));

        std::cerr << "[eval_arg] Unsupported argument: \"" << arg << "\"\n";
        return 0.0;
    }

    static double eval_arg_deriv(const std::string &arg, double x)
    {
        if (arg.size() > 2 && arg[0] == 'x' && arg[1] == '^')
        {
            int n = std::stoi(arg.substr(2));
            if (n == 0) return 0.0;
            return n * std::pow(x, n - 1);
        }
        return 0.0;
    }

    bool Function::valid_key(const std::string &key) const
    {
        static const std::regex poly_key_re(R"(^x\^\d+$)");
        if (std::regex_match(key, poly_key_re)) return true;

        // Trig key must have normalised arg, e.g. "sin(x^2)" — never bare "sin(x)"
        static const std::regex trig_key_re(R"(^(sin|cos|tan)\(x\^\d+\)$)");
        return std::regex_match(key, trig_key_re);
    }

    double Function::evaluate_key(const std::string &key, double x) const
    {
        if (key.size() > 2 && key[0] == 'x' && key[1] == '^')
            return std::pow(x, std::stoi(key.substr(2)));

        static const std::regex trig_key_re(R"(^(sin|cos|tan)\((x\^\d+)\)$)");
        std::smatch m;
        if (std::regex_match(key, m, trig_key_re))
        {
            double av = eval_arg(m[2].str(), x);
            const std::string &func = m[1].str();
            if (func == "sin") return std::sin(av);
            if (func == "cos") return std::cos(av);
            if (func == "tan") return std::tan(av);
        }

        std::cerr << "[Function::evaluate_key] Unknown key: \"" << key << "\"\n";
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

    // -------------------------------------------------------------------------
    //  derivate_evaluate
    //
    //    d/dx [ c * sin(g(x)) ] =  c * cos(g(x)) * g'(x)
    //    d/dx [ c * cos(g(x)) ] = -c * sin(g(x)) * g'(x)
    //    d/dx [ c * tan(g(x)) ] =  c * sec²(g(x)) * g'(x)
    //    d/dx [ c * x^N ]       =  c * N * x^(N-1)
    // -------------------------------------------------------------------------

    double Function::derivate_evaluate(double x) const
    {
        double result = 0.0;

        static const std::regex trig_key_re(R"(^(sin|cos|tan)\((x\^\d+)\)$)");

        for (const auto &[key, coef] : coeff)
        {
            double c = static_cast<double>(coef);

            std::smatch m;
            if (std::regex_match(key, m, trig_key_re))
            {
                const std::string &func = m[1].str();
                const std::string &arg  = m[2].str();

                double gx  = eval_arg(arg, x);      
                double gpx = eval_arg_deriv(arg, x); 

                if (func == "sin")
                    result += c * std::cos(gx) * gpx;
                else if (func == "cos")
                    result += c * -std::sin(gx) * gpx;
                else if (func == "tan")
                {
                    double sec = 1.0 / std::cos(gx);
                    result += c * sec * sec * gpx;
                }
                continue;
            }

            if (key.size() > 2 && key[0] == 'x' && key[1] == '^')
            {
                int degree = std::stoi(key.substr(2));
                if (degree == 0) continue; 
                result += c * degree * std::pow(x, degree - 1);
                continue;
            }

            std::cerr << "[Function::derivate_evaluate] Unknown key: \"" << key << "\"\n";
        }

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
            std::cerr << "[Function::update] Invalid key: \"" << key << "\"\n";
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
            std::cerr << "[Function::add] Invalid key: \"" << key << "\"\n";
            return;
        }
        if (coeff.count(key))
            std::cout << "[Function::add] Term \"" << key
                      << "\" already exists; use update() to change it.\n";
        else
            coeff.emplace(key, val);
    }

    // -------------------------------------------------------------------------
    //  extract_expression
    //
    //  Parses compact expressions such as:
    //    "3x^2 + 2sin(x) - 5cos(x^2) + tan(x^3) - 7"
    //
    //  RULES
    //  -----
    //  • Trig functions MUST be followed by a parenthesised monomial argument.
    //    Writing bare "sin", "cos", or "tan" without parens throws std::invalid_argument.
    //  • The argument inside the parens must be  x  or  x^N.
    //  • "x" inside parens is normalised to "x^1" so map keys are always "func(x^N)".
    //  • Arguments with inner +/- (e.g. sin(x+1)) are not supported.
    // -------------------------------------------------------------------------

    void Function::extract_expression(const std::string &expression)
    {
        std::string expr = expression;
        expr.erase(std::remove_if(expr.begin(), expr.end(), ::isspace), expr.end());

        static const std::regex bare_trig_re(R"((sin|cos|tan)(?!\())");
        if (std::regex_search(expr, bare_trig_re))
            throw std::invalid_argument(
                "[Function::extract_expression] Trigonometric functions require "
                "a parenthesised argument, e.g. sin(x) or cos(x^2).");

        std::vector<std::string> tokens;
        {
            std::string current;
            int depth = 0;
            for (std::size_t i = 0; i < expr.size(); ++i)
            {
                char c = expr[i];
                if (c == '(') ++depth;
                else if (c == ')') --depth;

                if ((c == '+' || c == '-') && depth == 0 && i != 0)
                {
                    tokens.push_back(current);
                    current.clear();
                }
                current += c;
            }
            if (!current.empty()) tokens.push_back(current);
        }

        static const std::regex trig_re(
            R"(^([+-]?)(\d+\.?\d*|\.?\d+)?(sin|cos|tan)\((x(?:\^(\d+))?)\)$)");

        static const std::regex poly_re(
            R"(^([+-]?)(\d+\.?\d*|\.?\d+)?(x?)(?:\^(\d+))?$)");

        for (const std::string &term : tokens)
        {
            std::smatch m;

            if (std::regex_match(term, m, trig_re))
            {
                std::string sign   = m[1].str();
                std::string digits = m[2].str();
                std::string func   = m[3].str(); 
                std::string arg    = m[4].str(); 
                std::string exp    = m[5].str(); 

                std::string norm_arg = exp.empty() ? "x^1" : ("x^" + exp);

                float value = digits.empty() ? 1.0f : std::stof(digits);
                if (sign == "-") value = -value;

                std::string key = func + "(" + norm_arg + ")";
                if (coeff.count(key)) coeff[key] += value;
                else                  coeff.emplace(key, value);
                continue;
            }

            if (std::regex_match(term, m, poly_re))
            {
                std::string sign    = m[1].str();
                std::string digits  = m[2].str();
                bool        has_x   = !m[3].str().empty();
                std::string exp_str = m[4].str();

                if (!has_x && digits.empty()) continue; 

                float value = digits.empty() ? 1.0f : std::stof(digits);
                if (sign == "-") value = -value;

                int degree = has_x ? (exp_str.empty() ? 1 : std::stoi(exp_str)) : 0;
                std::string key = "x^" + std::to_string(degree);

                if (coeff.count(key)) coeff[key] += value;
                else                  coeff.emplace(key, value);
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
            std::string prefix;
            if (first) { prefix = " "; first = false; }
            else        prefix = (coef >= 0) ? " + " : " - ";
            std::cout << prefix;

            float display = (coef < 0 && prefix != " ") ? -coef : coef;
            if (prefix == " " && coef < 0) { std::cout << "-"; display = -coef; }

            if (std::fabs(display - 1.0f) > 1e-6f || key == "x^0")
                std::cout << display;

            if (key != "x^0") std::cout << key;
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

    double bisection(Function func, double point_a, double point_b, double tolerance, int iterations) {
        double fa;
        double fb;
        double p;
        double fp;
        for (int i = 0; i < iterations; i++){
            fa = func.evaluate(point_a);
            fb = func.evaluate(point_b);
            p = point_a + ((point_b - point_a) / 2);
            fp = func.evaluate(p);
            if (fp == 0 || ((point_b - point_a) / 2) < tolerance) return p;
            if ( (fa * fp) > 0 ) point_a = p;
            else point_b = p;
        }
        return -1;
    }

    double fixed_point(Function func, double initial_point, double tolerance, int iterations) {
        double point = initial_point;
        double next_point;
        double f_next;
        for (int i = 0; i < iterations; i++){
            next_point = point - func.evaluate(point);
            f_next = func.evaluate(next_point);
            if (f_next == 0 || std::abs(next_point - point) < tolerance) return next_point;
            point = next_point;
        }
        return -1;
    }

    double fake_position(Function func, double point_a, double point_b, double tolerance, int iterations) {
        double fa;
        double fb;
        double p;
        double fp;
        for (int i = 0; i < iterations; i++) {
            fa = func.evaluate(point_a);
            fb = func.evaluate(point_b);
            p = ((point_a * fb) - (point_b * fa) ) / (fb - fa);
            fp = func.evaluate(p);
            if ( (fp * fa) < 0) {
                if (NumericalAnalysis::evaluate_tolerance(point_b, p, tolerance)) return p;
                point_b = p;
            }
            else if ( (fp * fb) < 0) {
                if (NumericalAnalysis::evaluate_tolerance(point_a, p, tolerance)) return p;
                point_a = p;
            }
        }
        return -1;
    }

    double newton_raphson(Function func, double initial_point, double tolerance, int iterations){
        double point = initial_point;
        double next_point;
        double fp;
        for (int i = 0; i < iterations; i++){
            next_point = point - (func.evaluate(point) / func.derivate_evaluate(point));
            if (std::abs(next_point - point) < tolerance) return next_point;
            point = next_point;
        }
        return -1;
    }

    double secant_method(Function func, double point_a, double point_b, double tolerance, int iterations) {
        // Entendemos como point_a = x(n-1), point_b = x(n) y p = x(n+1)
        double fa;
        double fb;
        double p;
        double fp;
        for (int i = 0; i < iterations; i++){
            fa = func.evaluate(point_a);
            fb = func.evaluate(point_b);
            p = ( (point_a * fb) - (point_b * fa) ) / (fb - fa);
            if (NumericalAnalysis::evaluate_tolerance(point_b, p, tolerance)) return p;
            point_a = point_b;
            point_b = p;
        }
        return -1;
    }
}