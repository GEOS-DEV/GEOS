class LinearSolverParameters
{
public:

  int          verbosity = 0;               //!< Output level [0=none, 1=basic, 2=everything]
  std::string  solverType = "cg";           //!< Solver type [direct, cg, gmres, bicgstab]
  std::string  preconditionerType = "ilut"; //!< Preconditioner type [none, ilu, ilut, icc, amg]
  int          dofsPerNode = 1;             //!< Can be used to enable dense-block algorithms if available

  struct
  {
    double  tolerance = 1e-6;
    int     maxIterations = 200;
    int     maxRestart = 200;
  }
  krylov;

  struct
  {
    bool useRowScaling = false;
    bool useRowColScaling = false;
  }
  scaling; //TODO: not implemented

  struct
  {
    int          maxLevels = 20;
    std::string  cycleType = "V";
    std::string  smootherType = "gaussSeidel";
    std::string  coarseType = "direct";
    int          numSweeps = 2;
    bool         isSymmetric = true;
    std::string  nullSpaceType = "constantModes";
  }
  amg;
  
  struct
  {
    int     fill = 0;
    double  threshold = 0.0;
  }
  ilu;

  // constructor
  LinearSolverParameters() = default;

  // destructor
  ~LinearSolverParameters() = default;
};
