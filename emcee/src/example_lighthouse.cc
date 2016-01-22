/**
 *  Application to the Light House problem.
 *  from Sivia and  Skilling 2006, p 192
 *
 *  "The aim of this application module is to solve the 'lighthouse' problem of
 *  Section 2.4, using the locations of flashes observed along the coastline to
 *  locate the lighthouse that emitted them in random directions. The lighthouse is
 *  here assumed to be somewhere in the rectangle −2 < x < 2 , 0 < y < 2 , with
 *  uniform prior.
 *  	
 *  	              u=0                                 u=1
 *  	               -------------------------------------
 *  	          y=2 |:::::::::::::::::::::::::::::::::::::| v=1
 *  	              |::::::::::::::::::::::LIGHT::::::::::|
 *  	         north|::::::::::::::::::::::HOUSE::::::::::|
 *  	              |:::::::::::::::::::::::::::::::::::::|
 *  	              |:::::::::::::::::::::::::::::::::::::|
 *  	          y=0 |:::::::::::::::::::::::::::::::::::::| v=0
 *  	 --*--------------*----*--------*-**--**--*-*-------------*--------
 *  	             x=-2          coastline -->east      x=2
 *
 *  Problem:
 *
 *  	  Lighthouse at (x,y) emitted n flashes observed at D[.] on coast.
 *
 *  	  Position    is 2-dimensional -2 < x < 2, 0 < y < 2 with flat prior
 *
 *  	  Likelihood  is L(x,y) = PRODUCT[k] (y/pi) / ((D[k] - x)^2 + y^2)
 *
 *  	  Evidence    is Z = INTEGRAL L(x,y) Prior(x,y) dxdy
 *
 *  	  Posterior   is P(x,y) = L(x,y) / Z estimating lighthouse position
 *
 *  	  Information is H = INTEGRAL P(x,y) log(P(x,y)/Prior(x,y)) dxdy
 *
 *  MCMC will only need the log-likelihood function (without specific prior)
 *
 *  Meanwhile, the position of the lighthouse as estimated from the given data
 *  was x = 1.24 +/- 0.18, y = 1.00 +/- 0.19. 
 *  The current code runs the sampling and display the x, y, lnp values on stdout
 */

#include <stdlib.h>
#include <random>
#include <iostream>
#include <iomanip>
#include <vector>
#include "emcee.hh"


/** 
 * likelihood of the lighthouse problem.
 *
 * @param pos  predicted position x,y
 *
 * @return lnp  ln-likelihood as defined above
 */
double lighthouse(std::vector<double>& pos){

    std::vector<double>
    D = { 4.73,  0.45, -1.73,  1.09,  2.19,  0.12,
          1.31,  1.00,  1.32,  1.07,  0.86, -0.49, -2.59,  1.73,  2.11,
          1.61,  4.98,  1.71,  2.23,-57.20,  0.96,  1.25, -1.56,  2.45,
          1.19,  2.17,-10.66,  1.91, -4.16,  1.92,  0.10,  1.98, -2.51,
          5.55, -0.47,  1.91,  0.95, -0.78, -0.84,  1.72, -0.01,  1.48,
          2.70,  1.21,  4.41, -4.79,  1.33,  0.81,  0.20,  1.58,  1.29,
          16.19,  2.75, -2.38, -1.79,  6.50,-18.53,  0.72,  0.94,  3.64,
          1.94, -0.11,  1.57,  0.57};

    double pi = 3.14159265;

    double lnp = 0.;

    for (size_t k=0; k < D.size(); ++k){
        lnp += std::log( (pos[1] / pi) / ( (D[k] - pos[0]) * (D[k] - pos[0]) + pos[1] * pos[1]));
    }

    return lnp;
}


int main(int argc, char *argv[]){

    size_t nwalkers = 40; /** number of Goodman & Weare walkers" */
    double a = 2.0;       /** proposal scale parameter */
    size_t ndim = 2;      /** number of dimensions in parameter space */
    size_t burn = 2000;    /** burning number of steps */
    size_t nsteps = 2000;  /** number of steps */

    // define function to optimize
    std::function<double (std::vector<double>&)> lnpfunc;
    // lnpfunc = test_lnprob;
    // lnpfunc = fn2;
    lnpfunc = lighthouse;
    
    // guess center of population
    std::vector<double> pos0(ndim);
    for (size_t dim = 0; dim < ndim; ++dim) {
        pos0[dim] = 1.; 
    }

    // init. walkers
    std::vector<std::vector<double> > walkers;
    std::vector<double> lnp(nwalkers);
    
    emcee::init_walkers_from_ball(nwalkers, walkers, pos0, 0.1);
    for (size_t k = 0; k < nwalkers; ++k) {
        lnp[k] = lnpfunc(walkers[k]);
    }

    // init. accept vector
    std::vector<bool> accept(nwalkers);
    for (size_t k = 0; k < nwalkers; ++k) {
        accept[k] = false;
    }

    // start by running the burn-in
    emcee::step(walkers, lnp, accept, lnpfunc, burn, a);
    //
    // run a production chain
    for(size_t step=0; step < nsteps; ++step){
        emcee::step(walkers, lnp, accept, lnpfunc);
        emcee::cout_walkers_lnp(walkers, lnp);
    }

    std::cout << "# acceptance ratio: " << emcee::compute_acceptance_ratio(accept) << std::endl;


    std::cout << "# Mean position: " ;
    for (size_t dim = 0;  dim < ndim ; ++dim) {
        double sum = 0.;
        for (size_t i = 0; i < nwalkers; ++i) {
            sum += walkers[i][dim];
        }
        double mean = sum / (double) nwalkers;
        double stddev = 0.;
        sum = 0.;
        for (size_t i = 0; i < nwalkers; ++i) {
            sum += (walkers[i][dim] - mean) * (walkers[i][dim] - mean);
        }
        stddev = std::sqrt(sum / ((double) nwalkers - 1));

        std::cout << mean << " +/- " << stddev << "   ";
    }
    std::cout << std::endl;
    return 0;
}
		
// vim: expandtab:ts=4:softtabstop=4:shiftwidth=4
