/*! \mainpage SuperSCS Documentation
 *
 * \section intro_sec Installation
 * 
 * \subsection installation_in_c For use in C
 *
 * 
 * First, you need to [download SuperSCS](https://github.com/kul-forbes/scs/archive/master.zip)
 * from the [github repo](https://github.com/kul-forbes/scs) of SuperSCS, or use the 
 * command:
 * 
 *     git clone https://github.com/kul-forbes/scs.git
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
 * using [valgrind](http://valgrind.org), do 
 * 
 *     make run-test-mem
 * 
 * 
 * For more advanced options, type in your terminal 
 * 
 *     make help
 * 
 * 
 * \subsection installation_in_matlab For use in MATLAB
 * 
 * In order to use SuperSCS from MATLAB, you need to compile a MEX interface.
 * 
 *     cd matlab
 *     make_scs; 
 *  
 * \section sec-first-steps First steps
 * 
 * \subsection main-examples Documentation
 * 
 * Examples of use:
 * 
 * - \ref page_socp
 * - \ref page_benchmarks
 * - \ref page_directions
 * - \ref page_dev
 * 
 * \subsection sec-github-page Source code
 * 
 * The source code of [SuperSCS](https://github.com/kul-forbes/scs) is available on [github](https://github.com/kul-forbes/scs).
 * 
 * This project was originally forked from [cvxgrp/scs](https://github.com/cvxgrp/scs).
 * 
 * \subsection sec-verification Verification
 * 
 * Quality assurance:
 * - [unit tests](https://github.com/kul-forbes/scs/tree/master/tests)
 * - [memory management tests](http://valgrind.org) using valgrind
 * - [continuous integration](https://travis-ci.org/kul-forbes/scs) on Travis CI
 * - [coverage reports](https://codecov.io/gh/kul-forbes/scs) on Codecov
 * - [code quality reports](https://www.codacy.com/app/alphaville/scs/dashboard) on codacy
 * - [lcov report] using <code>make cov</code>
 */

