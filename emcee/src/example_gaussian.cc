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
 * likelihood
 *
 * @param pos  predicted position (ndim)
 *
 * @return lnp  ln-likelihood as defined above
 */
double gaussian(std::vector<double>& pos){
    double lp = 0;
    for(size_t dim=0; dim < pos.size(); ++dim){
        lp += ((pos[dim]) * (pos[dim]));
    }
    return -0.5 * lp;
}

/** 
 * likelihood with arguments
 *
 * @param pos  predicted position (ndim)
 *
 * @param mean mean of the gaussian
 *
 * @param std  inverse variance
 *
 * @return lnp  ln-likelihood as defined above
 */
double gaussian_fn(std::vector<double>& pos, double mean, double ivar){
    double lp = 0;
    for(size_t dim=0; dim < pos.size(); ++dim){
        lp += ivar * ((pos[dim] - mean) * (pos[dim] - mean));
    }
    return -0.5 * lp;
}


int main(int argc, char *argv[]){

    size_t nwalkers = 40; /** number of Goodman & Weare walkers" */
    double a = 2.0;       /** proposal scale parameter */
    size_t ndim = 2;      /** number of dimensions in parameter space */
    size_t burn = 2000;    /** burning number of steps */
    size_t nsteps = 2000;  /** number of steps */

    // define function to optimize
    std::function<double (std::vector<double>&)> lnpfunc;
    auto partial_gaussian = std::bind(gaussian_fn,
            std::placeholders::_1, 2, 0.3);
    lnpfunc = partial_gaussian;
    
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
