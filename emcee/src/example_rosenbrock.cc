/**
 *  Example of fitting a cardinal sinus sin(t) / t
 */

#include <stdlib.h>
#include <random>
#include <iostream>
#include <iomanip>
#include <vector>
#include "emcee.hh"


/** 
 * Rosenbrock function with global minimum.
 *
 * non-convex function used as a performance test problem for
 * optimization algorithms introduced by Howard H. Rosenbrock
 * in 1960
 *
 * @param pos  predicted position (ndim)
 *
 * @return lnp  ln-likelihood as defined above
 */
double rosenbrock(std::vector<double>& pos){
    double b = 100.;
    double a = 20;

    if ( (pos[0] <= -5) || (pos[0] >= 5) || (pos[1] <= -5) || (pos[1] >= 5)){
        return -1e30;
    }

    double val =  (b * (pos[1] - pos[0] * pos[0]) * (pos[1] - pos[0] * pos[0]) +
        (1 - pos[0]) * (1 - pos[0]) ) / a;

    return -val;
}


int main(int argc, char *argv[]){

    size_t nwalkers = 40; /** number of Goodman & Weare walkers" */
    double a = 2.0;       /** proposal scale parameter */
    size_t ndim = 2;      /** number of dimensions in parameter space */
    size_t burn = 1000;    /** burning number of steps */
    size_t nsteps = 2000;  /** number of steps */

    // define function to optimize
    std::function<double (std::vector<double>&)> lnpfunc;
    lnpfunc = rosenbrock;
    
    // guess center of population
    std::vector<double> pos0(ndim);
    for (size_t dim = 0; dim < ndim; ++dim) {
        pos0[dim] = 0.; 
    }

    // init. walkers
    std::vector<std::vector<double> > walkers;
    std::vector<double> lnp(nwalkers);
    
    emcee::init_walkers_from_ball(nwalkers, walkers, pos0, 1.);
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
