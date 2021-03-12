#include "chi_runtime.h"
#include "ChiConsole/chi_console.h"
 
int main(int argc, char* argv[])
{
  ChiTech::Initialize(argc,argv);
  auto& console = ChiConsole::GetInstance();

  auto L  = console.consoleState;
  #include "LinearBoltzmannSolver/AffineDecomposition/lua/lua_register.h"

  ChiTech::RunBatch(argc, argv);

  ChiTech::Finalize();
  return 0;
}
