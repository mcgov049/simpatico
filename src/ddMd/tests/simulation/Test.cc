#include "SimulationTest.h"
//#include "SimulationTest2.h"
int main()
{
   #ifdef UTIL_MPI 
   MPI::Init();
   IntVector::commitMpiType();
   Vector::commitMpiType();
   #endif 

   TEST_RUNNER(SimulationTest) runner;
   runner.run();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

} 
