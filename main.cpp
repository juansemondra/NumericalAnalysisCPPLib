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

            case 0:
                std::cout << "Gracias por usar el sistema!" << std::endl;
                menu_continue = false;
            break;

            default:
            break;
        }
    }

    return EXIT_SUCCESS;
}