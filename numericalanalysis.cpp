#include "numericalanalysis.h"
#include <iostream>
#include <map>
#include <regex>
#include <iomanip>
#include <sstream>

namespace NumericalAnalysis {

    double bisection(Polynomial pol, double point_a, double point_b, int tolerance, int iterations){

        double current_point_a = point_a;
        double current_point_b = point_b;
        double next_point;
        double result;

        for (int i = 0; i < iterations; i++){

            next_point = (current_point_a + current_point_b) / 2;
            double evaluate_a = pol.evaluate(current_point_a);
            double evaluate_b = pol.evaluate(current_point_b);
            double evaluate_p = pol.evaluate(next_point);
            if (evaluate_a > 0 && evaluate_p > 0)
                current_point_a = next_point;
            else if (evaluate_b > 0 && evaluate_p > 0)
                current_point_b = next_point;
            else if (evaluate_a < 0 && evaluate_p < 0)
                current_point_a = next_point;
            else if (evaluate_b < 0 && evaluate_b < 0)
                current_point_b = next_point;
         
            if (evaluate_tolerance(evaluate_p, tolerance)) return result;
            
        }

        return result;
    }

    bool evaluate_tolerance(double value, int tolerance){

        std::ostringstream oss;
        oss << std::fixed << std::setprecision(15) << value;
        std::string s = oss.str();
        auto position = s.find('.');
        if (position == std::string::npos) return false;
        std::string fraction = s.substr(position + 1);
        if ((int)fraction.size() < tolerance) return false;

        int counter = 0;
        for(char x : fraction) if (x == '0') counter++;

        if (counter == tolerance) return true;
        return false;
    }

    Polynomial::Polynomial(){

    }

    double Polynomial::get(int degree){

        if (coeff.find(degree) == coeff.end()) return 0;
        else return coeff[degree];

    }

    void Polynomial::update(int degree, double val){

        if (auto search = coeff.find(degree); search != coeff.end()) coeff[degree] = val;
        else 
            std::cout << "Polynomial does not include degree x^" << degree << " in this function." << std::endl;
        
    }

    double Polynomial::evaluate(double value){

        double result = 0;

        for (auto it = coeff.begin(); it != coeff.end(); it++) result = result + (it->first * it->second);
        
        return result;
    }

    void Polynomial::add(int degree, double val){

        if (auto search = coeff.find(degree); search != coeff.end())
            std::cout << "Polynomial already includes degree x^" << degree << " in this function." << std::endl;
        else
            coeff.insert({degree, val});

    }

    void Polynomial::extract_expression(std::string expression){

        std::regex term_re(R"([+-]?[^+-]+)");
        std::regex parse_re(R"(^([+-]?)(\d*)(x?)(?:\^(\d+))?$)");
        
        for (std::sregex_iterator it(expression.begin(), expression.end(), term_re), end; it != end; ++it){
            std::string term = it->str();
            std::smatch m;

            if (!std::regex_match(term, m, parse_re)) {
                std::cerr << "Bad term: " << term << "\n";
                continue;
            }

            std::string sign = m[1].str();
            std::string digits = m[2].str();
            bool has_x = !m[3].str().empty();
            std::string exp_str = m[4].str();

            double value = 0;
            int degree = 0;

            if (has_x) {
                value = digits.empty() ? 1 : std::stoll(digits);
                degree = exp_str.empty() ? 1 : std::stoi(exp_str);
            } else {
                value = digits.empty() ? 1 : std::stoll(digits);
                degree = 0;
            }
            if (sign == "-") value = -value;
            if (value == 0) return;
            coeff.insert({degree, value});
        }   

    }

}