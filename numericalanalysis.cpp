#include "numericalanalysis.h"
#include <cmath>
#include <regex>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <fstream>
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

 
    // =====================================================================
    //  Matrix class implementation
    // =====================================================================

    Matrix::Matrix() : rows(0), columns(0) {}

    Matrix::Matrix(const std::string& filename) : rows(0), columns(0)
    {
        read_from_file(filename);
    }

    Matrix::Matrix(int rows, int columns)
        : data(rows, std::vector<double>(columns, 0.0)),
          rows(rows), columns(columns) {}

    Matrix::Matrix(const std::vector<std::vector<double>>& data)
        : data(data),
          rows(static_cast<int>(data.size())),
          columns(data.empty() ? 0 : static_cast<int>(data[0].size())) {}

    void Matrix::set(int row, int column, double value)
    {
        if (row < 0 || row >= rows || column < 0 || column >= columns)
        {
            std::cerr << "[Matrix::set] Index out of bounds ("
                      << row << ", " << column << ")\n";
            return;
        }
        data[row][column] = value;
    }

    double Matrix::get(int row, int column) const
    {
        if (row < 0 || row >= rows || column < 0 || column >= columns)
        {
            std::cerr << "[Matrix::get] Index out of bounds ("
                      << row << ", " << column << ")\n";
            return 0.0;
        }
        return data[row][column];
    }

    int Matrix::getRows() const { return rows; }
    int Matrix::getCols() const { return columns; }

    void Matrix::print() const
    {
        for (int i = 0; i < rows; i++)
        {
            std::cout << "| ";
            for (int j = 0; j < columns; j++)
                std::cout << std::setw(10) << std::fixed
                          << std::setprecision(4) << data[i][j] << " ";
            std::cout << "|\n";
        }
    }

    void Matrix::add(const Matrix& other)
    {
        if (rows != other.rows || columns != other.columns)
        {
            std::cerr << "[Matrix::add] Dimension mismatch\n";
            return;
        }
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++)
                data[i][j] += other.data[i][j];
    }

    void Matrix::subtract(const Matrix& other)
    {
        if (rows != other.rows || columns != other.columns)
        {
            std::cerr << "[Matrix::subtract] Dimension mismatch\n";
            return;
        }
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++)
                data[i][j] -= other.data[i][j];
    }

    void Matrix::multiply(const Matrix& other)
    {
        if (columns != other.rows)
        {
            std::cerr << "[Matrix::multiply] Incompatible dimensions ("
                      << rows << "x" << columns << ") * ("
                      << other.rows << "x" << other.columns << ")\n";
            return;
        }
        std::vector<std::vector<double>> result(
            rows, std::vector<double>(other.columns, 0.0));
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < other.columns; j++)
                for (int k = 0; k < columns; k++)
                    result[i][j] += data[i][k] * other.data[k][j];
        data = result;
        columns = other.columns;
    }

    void Matrix::divide(const Matrix& other)
    {
        if (rows != other.rows || columns != other.columns)
        {
            std::cerr << "[Matrix::divide] Dimension mismatch\n";
            return;
        }
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++)
            {
                if (std::abs(other.data[i][j]) < 1e-12)
                {
                    std::cerr << "[Matrix::divide] Division by zero at ("
                              << i << ", " << j << ")\n";
                    return;
                }
                data[i][j] /= other.data[i][j];
            }
    }

    void Matrix::transpose()
    {
        std::vector<std::vector<double>> result(
            columns, std::vector<double>(rows));
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++)
                result[j][i] = data[i][j];
        data = result;
        std::swap(rows, columns);
    }

    void Matrix::inverse()
    {
        if (rows != columns)
        {
            std::cerr << "[Matrix::inverse] Matrix must be square\n";
            return;
        }
        int n = rows;
        std::vector<std::vector<double>> aug(n, std::vector<double>(2 * n, 0.0));
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                aug[i][j] = data[i][j];
            aug[i][n + i] = 1.0;
        }

        for (int i = 0; i < n; i++)
        {
            int maxRow = i;
            for (int k = i + 1; k < n; k++)
                if (std::abs(aug[k][i]) > std::abs(aug[maxRow][i]))
                    maxRow = k;
            std::swap(aug[i], aug[maxRow]);

            if (std::abs(aug[i][i]) < 1e-12)
            {
                std::cerr << "[Matrix::inverse] Singular matrix, cannot invert\n";
                return;
            }

            double pivot = aug[i][i];
            for (int j = 0; j < 2 * n; j++)
                aug[i][j] /= pivot;

            for (int k = 0; k < n; k++)
            {
                if (k != i)
                {
                    double factor = aug[k][i];
                    for (int j = 0; j < 2 * n; j++)
                        aug[k][j] -= factor * aug[i][j];
                }
            }
        }

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                data[i][j] = aug[i][n + j];
    }

    double Matrix::determinant()
    {
        if (rows != columns)
        {
            std::cerr << "[Matrix::determinant] Matrix must be square\n";
            return 0.0;
        }
        int n = rows;
        auto temp = data;
        double det = 1.0;
        int sign = 1;

        for (int i = 0; i < n; i++)
        {
            int maxRow = i;
            for (int k = i + 1; k < n; k++)
                if (std::abs(temp[k][i]) > std::abs(temp[maxRow][i]))
                    maxRow = k;
            if (maxRow != i)
            {
                std::swap(temp[i], temp[maxRow]);
                sign = -sign;
            }
            if (std::abs(temp[i][i]) < 1e-12)
                return 0.0;
            det *= temp[i][i];
            for (int k = i + 1; k < n; k++)
            {
                double factor = temp[k][i] / temp[i][i];
                for (int j = i + 1; j < n; j++)
                    temp[k][j] -= factor * temp[i][j];
            }
        }
        return sign * det;
    }

    int Matrix::rank()
    {
        auto temp = data;
        int r = 0;
        for (int col = 0; col < columns && r < rows; col++)
        {
            int pivot = -1;
            for (int row = r; row < rows; row++)
            {
                if (std::abs(temp[row][col]) > 1e-12)
                {
                    pivot = row;
                    break;
                }
            }
            if (pivot == -1) continue;
            std::swap(temp[r], temp[pivot]);
            for (int row = r + 1; row < rows; row++)
            {
                if (std::abs(temp[r][col]) < 1e-12) continue;
                double factor = temp[row][col] / temp[r][col];
                for (int j = col; j < columns; j++)
                    temp[row][j] -= factor * temp[r][j];
            }
            r++;
        }
        return r;
    }

    void Matrix::read_from_file(const std::string& filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "[Matrix::read_from_file] No se pudo abrir: "
                      << filename << "\n";
            return;
        }

        data.clear();
        std::string line;
        while (std::getline(file, line))
        {
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::vector<double> row;
            std::string token;
            while (iss >> token)
            {
                size_t slash = token.find('/');
                if (slash != std::string::npos)
                {
                    double num = std::stod(token.substr(0, slash));
                    double den = std::stod(token.substr(slash + 1));
                    row.push_back(num / den);
                }
                else
                {
                    row.push_back(std::stod(token));
                }
            }
            if (!row.empty())
                data.push_back(row);
        }

        rows = static_cast<int>(data.size());
        columns = rows > 0 ? static_cast<int>(data[0].size()) : 0;
    }

    void Matrix::write_to_file(const std::string& filename) const
    {
        std::ofstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "[Matrix::write_to_file] No se pudo abrir: "
                      << filename << "\n";
            return;
        }

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                if (j > 0) file << " ";
                file << data[i][j];
            }
            file << "\n";
        }
    }

    // =====================================================================

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
        for (int i = 0; i < iterations; i++){
            if (std::abs(func.derivate_evaluate(point)) < 1e-12) return -1;
            next_point = point - (func.evaluate(point) / func.derivate_evaluate(point));
            if (std::abs(next_point - point) < tolerance) return next_point;
            point = next_point;
        }
        return -1;
    }

    double secant_method(Function func, double point_a, double point_b, double tolerance, int iterations) {
        // Entendemos como point_a = x(n-1), point_b = x(n) y p = x(n+1)
        double fa, fb, p;
        for (int i = 0; i < iterations; i++){
            fa = func.evaluate(point_a);
            fb = func.evaluate(point_b);
            if (std::abs(fb - fa) < 1e-12) return -1;
            p = ( (point_a * fb) - (point_b * fa) ) / (fb - fa);
            if (NumericalAnalysis::evaluate_tolerance(point_b, p, tolerance)) return p;
            point_a = point_b;
            point_b = p;
        }
        return -1;
    }

    // =====================================================================
    //  Segundo Corte — Sistemas de Ecuaciones Lineales
    // =====================================================================

    // -----------------------------------------------------------------
    //  Sustitución Regresiva (Tarea 1)
    //
    //  Entrada: Matriz aumentada [R|c] de n×(n+1), donde R es
    //           triangular superior con r_ii ≠ 0.
    //  Salida:  Vector solución x como Matrix de n×1.
    //
    //  Algoritmo:
    //    x_n = c_n / r_nn
    //    Para i = (n-1) hasta 1:
    //      suma = Σ_{j=i+1}^{n} r_ij * x_j
    //      x_i = (c_i - suma) / r_ii
    // -----------------------------------------------------------------

    Matrix regressive_substitution(Matrix matrix)
    {
        int n = matrix.getRows();
        Matrix x(n, 1);

        if (matrix.getCols() != n + 1)
        {
            std::cerr << "[regressive_substitution] La matriz debe ser aumentada "
                      << n << "x" << (n + 1) << ", se recibió "
                      << n << "x" << matrix.getCols() << "\n";
            return x;
        }

        double rnn = matrix.get(n - 1, n - 1);
        if (std::abs(rnn) < 1e-12)
        {
            std::cerr << "[regressive_substitution] r_nn = 0, no se puede resolver\n";
            return x;
        }
        x.set(n - 1, 0, matrix.get(n - 1, n) / rnn);

        for (int i = n - 2; i >= 0; i--)
        {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++)
                sum += matrix.get(i, j) * x.get(j, 0);

            double rii = matrix.get(i, i);
            if (std::abs(rii) < 1e-12)
            {
                std::cerr << "[regressive_substitution] r_" << i+1 << i+1
                          << " = 0, no se puede resolver\n";
                return x;
            }
            x.set(i, 0, (matrix.get(i, n) - sum) / rii);
        }

        return x;
    }

    // -----------------------------------------------------------------
    //  Eliminación Gaussiana — solo eliminación hacia adelante
    //
    //  Acepta cualquier matriz (cuadrada o aumentada).
    //  Retorna la matriz en forma escalonada (triangular superior).
    // -----------------------------------------------------------------

    Matrix gaussian_elimination_step(Matrix matrix)
    {
        int n = matrix.getRows();
        int cols = matrix.getCols();

        for (int i = 0; i < n - 1; i++)
        {
            int p = -1;
            for (int k = i; k < n; k++)
            {
                if (std::abs(matrix.get(k, i)) > 1e-12)
                {
                    p = k;
                    break;
                }
            }

            if (p == -1) continue;

            if (p != i)
            {
                for (int j = 0; j < cols; j++)
                {
                    double temp = matrix.get(i, j);
                    matrix.set(i, j, matrix.get(p, j));
                    matrix.set(p, j, temp);
                }
            }

            for (int j = i + 1; j < n; j++)
            {
                double mji = matrix.get(j, i) / matrix.get(i, i);
                for (int k = i; k < cols; k++)
                    matrix.set(j, k, matrix.get(j, k) - mji * matrix.get(i, k));
            }
        }

        return matrix;
    }

    // -----------------------------------------------------------------
    //  Eliminación Gaussiana con sustitución hacia atrás (Tarea 2)
    //
    //  Entrada: Matriz aumentada [A|b] de n×(n+1).
    //  Salida:  Vector solución x como Matrix de n×1.
    // -----------------------------------------------------------------

    Matrix gaussian_elimination_with_regressive_substitution(Matrix matrix)
    {
        int n = matrix.getRows();

        if (matrix.getCols() != n + 1)
        {
            std::cerr << "[gaussian_elimination] La matriz debe ser aumentada "
                      << n << "x" << (n + 1) << ", se recibió "
                      << n << "x" << matrix.getCols() << "\n";
            return Matrix(n, 1);
        }

        for (int i = 0; i < n - 1; i++)
        {
            // Buscar el menor p ≥ i tal que a[p][i] ≠ 0
            int p = -1;
            for (int k = i; k < n; k++)
            {
                if (std::abs(matrix.get(k, i)) > 1e-12)
                {
                    p = k;
                    break;
                }
            }

            if (p == -1)
            {
                std::cerr << "[gaussian_elimination] No existe solución única\n";
                return Matrix(n, 1);
            }

            // Intercambio de filas si p ≠ i
            if (p != i)
            {
                for (int j = 0; j <= n; j++)
                {
                    double temp = matrix.get(i, j);
                    matrix.set(i, j, matrix.get(p, j));
                    matrix.set(p, j, temp);
                }
            }

            // Eliminación: E_j ← E_j - m_ji * E_i
            for (int j = i + 1; j < n; j++)
            {
                double mji = matrix.get(j, i) / matrix.get(i, i);
                for (int k = i; k <= n; k++)
                    matrix.set(j, k, matrix.get(j, k) - mji * matrix.get(i, k));
            }
        }

        if (std::abs(matrix.get(n - 1, n - 1)) < 1e-12)
        {
            std::cerr << "[gaussian_elimination] No existe solución única (a_nn = 0)\n";
            return Matrix(n, 1);
        }

        return regressive_substitution(matrix);
    }

    // -----------------------------------------------------------------
    //  Factorización LU con pivoteo parcial — PA = LU
    //
    //  Entrada: Matriz cuadrada A de n×n.
    //  Salida:  Matrices L y U tales que PA = LU.
    // -----------------------------------------------------------------

    void lu_factorization(Matrix A, Matrix& L, Matrix& U)
    {
        int n = A.getRows();

        if (A.getCols() != n)
        {
            std::cerr << "[lu_factorization] Se esperaba una matriz cuadrada, se recibió "
                      << n << "x" << A.getCols() << "\n";
            L = Matrix(); U = Matrix();
            return;
        }

        Matrix work(n, n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                work.set(i, j, A.get(i, j));

        for (int j = 0; j < n; j++)
        {
            int r = j;
            double maxVal = std::abs(work.get(j, j));
            for (int i = j + 1; i < n; i++)
            {
                double val = std::abs(work.get(i, j));
                if (val > maxVal) { maxVal = val; r = i; }
            }

            if (maxVal < 1e-12)
            {
                std::cerr << "[lu_factorization] Matriz singular, no se puede factorizar\n";
                L = Matrix(); U = Matrix();
                return;
            }

            if (r != j)
            {
                for (int k = 0; k < n; k++)
                {
                    double temp = work.get(j, k);
                    work.set(j, k, work.get(r, k));
                    work.set(r, k, temp);
                }
            }

            for (int i = j + 1; i < n; i++)
            {
                double mij = work.get(i, j) / work.get(j, j);
                work.set(i, j, mij);
                for (int k = j + 1; k < n; k++)
                    work.set(i, k, work.get(i, k) - mij * work.get(j, k));
            }
        }

        L = Matrix(n, n);
        U = Matrix(n, n);
        for (int i = 0; i < n; i++)
        {
            L.set(i, i, 1.0);
            for (int j = 0; j < i; j++)
                L.set(i, j, work.get(i, j));
            for (int j = i; j < n; j++)
                U.set(i, j, work.get(i, j));
        }
    }

    // -----------------------------------------------------------------
    //  Sustitución LU — Resuelve Ax = b usando PA = LU (Tarea 3)
    //
    //  Entrada: Matriz aumentada [A|b] de n×(n+1).
    //  Salida:  Vector solución x como Matrix de n×1.
    // -----------------------------------------------------------------

    Matrix lu_substitution(Matrix matrix)
    {
        int n = matrix.getRows();

        if (matrix.getCols() != n + 1)
        {
            std::cerr << "[lu_substitution] La matriz debe ser aumentada "
                      << n << "x" << (n + 1) << ", se recibió "
                      << n << "x" << matrix.getCols() << "\n";
            return Matrix(n, 1);
        }

        std::vector<int> perm(n);
        for (int i = 0; i < n; i++) perm[i] = i;

        for (int j = 0; j < n; j++)
        {
            // 2.1 Pivoteo parcial: buscar r con |a_rj| máximo
            int r = j;
            double maxVal = std::abs(matrix.get(j, j));
            for (int i = j + 1; i < n; i++)
            {
                double val = std::abs(matrix.get(i, j));
                if (val > maxVal)
                {
                    maxVal = val;
                    r = i;
                }
            }

            // 2.2 Chequeo de singularidad
            if (maxVal < 1e-12)
            {
                std::cerr << "[lu_substitution] Matriz singular, no se puede factorizar\n";
                return Matrix(n, 1);
            }

            // 2.3 Intercambio de filas
            if (r != j)
            {
                std::swap(perm[j], perm[r]);
                for (int k = 0; k <= n; k++)
                {
                    double temp = matrix.get(j, k);
                    matrix.set(j, k, matrix.get(r, k));
                    matrix.set(r, k, temp);
                }
            }

            // 2.4 Eliminación
            for (int i = j + 1; i < n; i++)
            {
                double mij = matrix.get(i, j) / matrix.get(j, j);
                matrix.set(i, j, mij);
                matrix.set(i, n, matrix.get(i, n) - mij * matrix.get(j, n));
                for (int k = j + 1; k < n; k++)
                    matrix.set(i, k, matrix.get(i, k) - mij * matrix.get(j, k));
            }
        }

        // 3. Chequeo final
        if (std::abs(matrix.get(n - 1, n - 1)) < 1e-12)
        {
            std::cerr << "[lu_substitution] Matriz singular (a_nn = 0)\n";
            return Matrix(n, 1);
        }

        // 5. Resolver Ux = b_modificado con sustitución regresiva
        Matrix upper(n, n + 1);
        for (int i = 0; i < n; i++)
        {
            for (int jj = i; jj < n; jj++)
                upper.set(i, jj, matrix.get(i, jj));
            upper.set(i, n, matrix.get(i, n));
        }

        return regressive_substitution(upper);
    }

    // -----------------------------------------------------------------
    //  Método Iterativo de Gauss-Seidel (Tarea 4)
    //
    //  Entrada: Matriz aumentada [A|b] de n×(n+1),
    //           vector inicial x0 (Matrix n×1),
    //           tolerancia E, máximo de iteraciones N.
    //  Salida:  Vector solución x como Matrix de n×1.
    //
    //  Algoritmo (pg. 44 de los apuntes):
    //    k = 1
    //    Mientras k ≤ N:
    //      Para i = 1..n:
    //        x_i = (1/a_ii)(b_i - Σ_{j<i} a_ij*x_j - Σ_{j>i} a_ij*x0_j)
    //        (usa valores nuevos para j<i, valores viejos para j>i)
    //      Si ||x - x0||∞ < E → Éxito
    //      x0 ← x
    //    Fracaso
    //
    //  Convergencia garantizada si A es diagonal estrictamente
    //  dominante sobre filas: |a_ii| > Σ_{j≠i} |a_ij|
    // -----------------------------------------------------------------

    Matrix gauss_seidel(Matrix matrix, Matrix initial, double tolerance, int iterations)
    {
        int n = matrix.getRows();

        if (matrix.getCols() != n + 1)
        {
            std::cerr << "[gauss_seidel] La matriz debe ser aumentada "
                      << n << "x" << (n + 1) << ", se recibió "
                      << n << "x" << matrix.getCols() << "\n";
            return Matrix(n, 1);
        }

        Matrix x0(n, 1);
        Matrix x(n, 1);

        for (int i = 0; i < n; i++)
            x0.set(i, 0, initial.get(i, 0));

        std::cout << "\n--- Iteraciones de Gauss-Seidel ---\n";
        std::cout << std::fixed << std::setprecision(8);

        for (int k = 0; k < iterations; k++)
        {
            for (int i = 0; i < n; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < n; j++)
                {
                    if (j < i)
                        sum += matrix.get(i, j) * x.get(j, 0);
                    else if (j > i)
                        sum += matrix.get(i, j) * x0.get(j, 0);
                }

                double aii = matrix.get(i, i);
                if (std::abs(aii) < 1e-12)
                {
                    std::cerr << "[gauss_seidel] a_" << i+1 << i+1
                              << " = 0, no se puede resolver\n";
                    return x;
                }
                x.set(i, 0, (matrix.get(i, n) - sum) / aii);
            }

            double norm = 0.0;
            for (int i = 0; i < n; i++)
            {
                double diff = std::abs(x.get(i, 0) - x0.get(i, 0));
                if (diff > norm) norm = diff;
            }

            std::cout << "  Iter " << std::setw(3) << (k + 1) << ":  x = [";
            for (int i = 0; i < n; i++)
            {
                if (i > 0) std::cout << ", ";
                std::cout << std::setw(14) << x.get(i, 0);
            }
            std::cout << " ]  ||e|| = " << norm << "\n";

            if (norm < tolerance)
            {
                std::cout << "  Convergencia alcanzada en " << (k + 1) << " iteraciones.\n";
                return x;
            }

            for (int i = 0; i < n; i++)
                x0.set(i, 0, x.get(i, 0));
        }

        std::cerr << "[gauss_seidel] Se excedió el máximo de iteraciones ("
                  << iterations << ")\n";
        return x;
    }

    double inferior_sums(Function func, double a, double b, int n)
    {
        double dx = (b - a) / n;
        double sum = 0.0;

        for (int i = 1; i <= n; i++)
        {
            double f_prev = func.evaluate(a + (i - 1) * dx);
            double f_curr = func.evaluate(a + i * dx);
            sum += std::min(f_prev, f_curr) * dx;
        }

        return sum;
    }

    double superior_sums(Function func, double a, double b, int n)
    {
        double dx = (b - a) / n;
        double sum = 0.0;

        for (int i = 1; i <= n; i++)
        {
            double f_prev = func.evaluate(a + (i - 1) * dx);
            double f_curr = func.evaluate(a + i * dx);
            sum += std::max(f_prev, f_curr) * dx;
        }

        return sum;
    }

    double trapezoidal_rule(Function func, double a, double b, int n)
    {
        double h = (b - a) / n;
        double s0 = func.evaluate(a) + func.evaluate(b);
        double s1 = 0.0;

        for (int i = 1; i < n; i++)
            s1 += func.evaluate(a + i * h);

        return (h / 2.0) * (s0 + 2.0 * s1);
    }

    double simpson_rule(Function func, double a, double b, int n)
    {
        if (n <= 0 || n % 2 != 0)
        {
            std::cerr << "[simpson_rule] n debe ser un entero positivo par\n";
            return -1;
        }

        double h = (b - a) / n;
        double s0 = func.evaluate(a) + func.evaluate(b);
        double s1 = 0.0; // Terminos con indice impar
        double s2 = 0.0; // Terminos con indice par

        for (int i = 1; i < n; i++)
        {
            double fx = func.evaluate(a + i * h);
            if (i % 2 == 0) s2 += fx;
            else            s1 += fx;
        }

        return (h / 3.0) * (s0 + 4.0 * s1 + 2.0 * s2);
    }
}