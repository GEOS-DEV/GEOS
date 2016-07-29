#include "FractalSurface.h"
#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{
  //Instantiate aperturator class
  FractalSurface fs;

  //Read parameters from the first file
  {
    std::ifstream input;
    input.open(argv[1]);
    if(!input){
      std::cout << "Error 1" << std::endl;
      return 1;
    }

    R1TensorT<2> lower(0.0);
    R1TensorT<2> upper(1.0);
    localIndex nlevels, n0, n1;
    realT hfct;
    input >> lower(0) >> lower(1) >> upper(0) >> upper(1)
        >> n0 >> n1 >> hfct >> nlevels;

    if(nlevels == 0)
    {
      realT hurst = 1.3, mean = 0.0, stdev = 1.0;
      input >> nlevels >> mean >> stdev >> hurst;
      fs.InitializeFractal(mean, stdev, lower, upper, hurst, nlevels, n0, n1, hfct);
    }
    else
    {
      Array2dT<realT> parameters;
      parameters.resize2(nlevels, 2);
      for(localIndex i = 0; i < nlevels; i++)
        input >> parameters(i,0) >> parameters(i,1);
      fs.Initialize(parameters, lower, upper, n0, n1, hfct);
    }

    input.close();
  }

  //Read points in from the second file and output apertures
  {
    std::ifstream input;
    input.open(argv[2]);
    if(!input){
      std::cout << "Error 2" << std::endl;
      return 1;
    }

    int num = 0;
    input >> num;

    R1TensorT<2> x(0.0);
    for(int i = 0; i < num; i++) {
      input >> x(0) >> x(1);
      realT aa = fs.Value(x);
      std::cout << x(0) << " " << x(1) << " 0 " << aa << "\n";
    }

    input.close();
  }

  return 0;
}
