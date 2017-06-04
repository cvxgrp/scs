/*! \mainpage SuperSCS Documentation
 *
 * \section intro_sec Installation
 *
 * What SuperSCS is about
 * 
 * Installation runs as simple as 
 * 
 *     make
 *  
 * If you want to run the tests, do 
 *
 *     make run-test
 *
 * 
 * If, additionally, you want to run the tests and perform a memory check
 * using valgrind, do 
 * 
 *     make run-test-mem
 * 
 * 
 * For more advanced options, type in your terminal 
 * 
 *     make help
 * 
 * 
 * 
 * \section examples-of-use Examples of use 
 * 
 * ~~~~~~{.c}
 * sol = initSol();
 * info = initInfo();
 * data->stgs->do_super_scs = 1;
 * scs(data, cone, sol, info);
 * freeData(data, cone);
 * freeSol(sol);
 * scs_free(info);
 * ~~~~~~
 */

