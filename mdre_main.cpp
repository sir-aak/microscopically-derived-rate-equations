// version: 18.09.2020

#include "mdre_parameters.hpp"
#include "mdre_command_line_options.hpp"
#include "mdre_node_equations.hpp"
#include "mdre_network.hpp"
#include "mdre_solver.hpp"
#include "mdre_time_series.hpp"
#include "mdre_auxiliary.hpp"
#include "mdre_analyzer.hpp"
#include "mdre_TISI_finder.hpp"
#include "mdre_features.hpp"
#include "mdre_function_switches.hpp"


using namespace std;


// this program provides features to analyze networks of the 
// microscopically derived rate equations (MDRE) based on LIN14
int main (int argc, char* argv[])
{
  
  setCommandLineOptions(argc, argv);
  enableFunctionSwitches(argc, argv);
  
  return (0);
  
}

