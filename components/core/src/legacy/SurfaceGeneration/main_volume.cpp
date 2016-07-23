#include "FractalVolume.h"
#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{
  //Instantiate
  FractalVolume fv;

  //Read parameters from the first file
  {
    std::ifstream input;
    input.open(argv[1]);
    if(!input){
      std::cout << "Error 1" << std::endl;
      return 1;
    }

    R1Tensor lower(0.0);
    R1Tensor upper(1.0);
    localIndex nlevels, n0, n1, n2;
    realT hfct;
    input >> lower(0) >> lower(1) >> lower(2)
        >> upper(0) >> upper(1) >> upper(2)
        >> n0 >> n1 >> n2 >> hfct >> nlevels;

    if(nlevels == 0)
    {
      realT hurst = 1.3, mean = 0.0, stdev = 1.0;
      input >> nlevels >> mean >> stdev >> hurst;
      fv.InitializeFractal(mean, stdev, lower, upper, hurst, nlevels, n0, n1, n2, hfct);
    }
    else
    {
      Array2dT<realT> parameters;
      parameters.resize2(nlevels, 2);
      for(localIndex i = 0; i < nlevels; i++)
        input >> parameters(i,0) >> parameters(i,1);
      fv.Initialize(parameters, lower, upper, n0, n1, n2, hfct);
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

    R1Tensor x(0.0);
    for(int i = 0; i < num; i++) {
      input >> x(0) >> x(1) >> x(2);
      const realT val = fv.Value(x);
      std::cout << x(0) << " " << x(1) << " " << x(2) << " " << val << "\n";
    }

    input.close();
  }

  return 0;
}
