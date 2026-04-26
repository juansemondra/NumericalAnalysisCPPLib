#include <iostream>
#include "numericalanalysis.h"
#include "menu.h"

int main(int argc, char **argv)
{

    for (int i = 0; i < argc; i++)
    {
        std::string temporal = argv[i];
        if (temporal == "--help")
        {
            helper_function();
            return EXIT_SUCCESS;
        }
    }
    bool menu_continue = true;
    int menu_option;

    while (menu_continue){
        print_menu();
        std::cin >> menu_option;
        std::cout << std::endl;
        switch (menu_option){
            case 1:
                check_error(call_bisection());
            break;

            case 2:
                check_error(call_fixed_point());
            break;

            case 3:
                check_error(call_fake_position());
            break;

            case 4:
                check_error(call_newton_raphson());
            break;

            case 5:
                check_error(call_secant_method());
            break;

            case 6:
                call_regressive_substitution();
            break;

            case 7:
                call_gaussian_elimination();
            break;

            case 8:
                call_lu_substitution();
            break;

            case 9:
                call_gauss_seidel();
            break;

            case 10:
                call_inferior_sums();
            break;

            case 11:
                call_superior_sums();
            break;

            case 12:
                call_trapezoidal_rule();
            break;

            case 13:
                call_simpson_rule();
            break;

            case 0:
                std::cout << "Gracias por usar el sistema!" << std::endl;
                menu_continue = false;
            break;

            default:
                std::cout << "Opción inválida, por favor intente nuevamente" << std::endl;
            break;
        }
    }

    return EXIT_SUCCESS;
}