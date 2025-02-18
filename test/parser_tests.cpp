#include "high_exponent_real.hpp"
#include "view_parameters.hpp"
#include "real_number.hpp"
#include <iomanip>
#include <sstream>

auto test1 = "Location:\n"
             "Re: "
             "-1."
             "47994622332507888020258065344256383359082887482853327232891946750"
             "45014280415514581021231577152136510355459435420781673489538857873"
             "41902612509987523844493790997569942182815961274\n"
             "Im: "
             "0."
             "00090139732902035398019779186619717356625156681804510241106763038"
             "64886922871890491456215844360269342187635275772906318094547966181"
             "107674245803000517934891341563129581287143605336i\n"
             "Magnification: 2.2e141 (where 1.0 fits entire Mandelbrot set "
             "into the window)\n"
             "Scale: 6.6e143 (pixels/unit, for renders seen here)\n"
             "Iterations: 80000\n";

int main() {
  using H512 = fractals::real_number<512, 0, 0>;

  // Numerical exponents on high_precision_real
  {
    H512 n1;
    std::stringstream test1("1.234e-50");
    test1 >> n1;
    std::cout << std::setprecision(60) << n1 << std::endl;
  }

  // Test numerical parsing

  {
    H512 n1 = 3.1415e-50;
    std::stringstream ss;
    ss << n1;
    std::cout << std::setprecision(60) << n1;
  }

  {
    std::stringstream ss(test1);
    fractals::view_parameters p1;
    try_import(ss, p1);
    ss >> p1;
    std::cout << p1 << std::endl;
  }
}
