#pragma once
#include <stdlib.h>
#include <random>
#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>

namespace emcee { 

    /**
     *  This subroutine advances an ensemble of walkers using the
     *  Goodman & Weare stretch move.
     * 
     *  @param ndim 	The dimension of the parameter space.
     * 
     *  @param nwalkers	The number of walkers.
     * 
     *  @param a 		The proposal scale (a tuning parameter). 
     *                  Using `a=2` is almost always the right move.
     * 
     *  @param pin      starting positions (ndim, nwalkers)]
     * 
     *  @param lpin     log-probability at positions `pin` (nwalkers,)
     * 
     *  @param pout    	The new positions of the walkers in parameter space.
     * 
     *  @param lpout 	log-proba at `pout`
     * 
     *  @param accept 	A binary list indicating whether or not each proposal
     *                  was accepted.
     *
     *  @param lnpfunc  function to sample
     */
    void advance(size_t ndim, size_t nwalkers, 
        double a, 
        std::vector<std::vector<double> >& pin, 
        std::vector<double>& lpin,
        std::vector<std::vector<double> >& pout, 
        std::vector<double>& lpout,
        std::vector<bool>& accept,
        std::function<double (std::vector<double>& )>& lnpfunc
        );

    /**
     *  This subroutine advances an ensemble of walkers using the
     *  Goodman & Weare stretch move.
     * 
     *  @param walkers  starting positions [nwalkers][ndim]
     * 
     *  @param lnp      log-probability at positions `walkers` [nwalkers]
     * 
     *  @param accept 	A binary list indicating whether or not each proposal
     *                  was accepted.
     *
     *  @param lnpfunc  function to sample
     *
     *  @param nsteps   number of steps to advance (default = 1)
     *
     *  @param a 		The proposal scale (a tuning parameter). 
     *                  Using `a=2` is almost always the right move.
     */
    void step(
        std::vector<std::vector<double> >& walkers, 
        std::vector<double>& lnp,
        std::vector<bool>& accept,
        std::function<double (std::vector<double>& )>& lnpfunc,
        size_t nstep = 1,
        double a = 2.);

    /**
     * convenient function that initializes a set of walkers from a given
     * position randomly within a Gaussian ball.
     *
     * @param nwalkers  number of walkers to generate
     *
     * @param walkers   vector to populate
     *
     * @param pos0      central position [ndim]
     *
     * @param radius   size of the ball
     */
    void init_walkers_from_ball(size_t nwalkers,
        std::vector<std::vector<double> >& walkers,
        std::vector<double>& pos0,
        double radius);

    /** convenient function that computes the acceptance ratio given the accept
     * vector
     *
     *  @param accept  accept vector [nwalkers]
     *
     *  @return val   acceptance ratio
     *
     */

    float compute_acceptance_ratio(std::vector<bool>& accept);

    /** convenient function that cout walker positions and lnp
     *
     *  @param walkers  starting positions [nwalkers][ndim]
     * 
     *  @param lnp      log-probability at positions `walkers` [nwalkers]
     *
     */
    void cout_walkers_lnp(
        std::vector<std::vector<double> >& walkers,
        std::vector<double>& lnp);
}

// vim: expandtab:ts=4:softtabstop=4:shiftwidth=4
