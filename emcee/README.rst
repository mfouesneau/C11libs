emcee -- The MCMC hammer c-version
==================================


`emcee`_ is a pure-Python implementation of `Goodman & Weare’s Affine Invariant
Markov chain Monte Carlo (MCMC) Ensemble sampler <http://msp.berkeley.edu/camcos/2010/5-1/p04.xhtml>`_.

This repository is a quick C++11 implementation.



Basic Usage
-----------

If you wanted to draw samples from a 10 dimensional Gaussian, you would do something like:

.. code:: cpp

        include "emcee.hh"

        /** 
         * likelihood with arguments
         *
         * @param pos  predicted position (ndim)
         * @param mean mean of the gaussian
         * @param std  inverse variance
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
              
            size_t nwalkers = 40; /** number of Goodman & Weare walkers"      */
            double a = 2.0;       /** proposal scale parameter                */
            size_t ndim = 10;     /** number of dimensions in parameter space */
            size_t nsteps = 2000; /** number of steps                         */

            // define function to optimize
            std::function<double (std::vector<double>&)> lnpfunc;
            auto partial_gaussian = std::bind(gaussian_fn, std::placeholders::_1, 2, 0.3);
            lnpfunc = partial_gaussian;

            // guess center of population
            std::vector<double> pos0(ndim);
            for (size_t dim = 0; dim < ndim; ++dim) { pos0[dim] = 1.; }
      
            // init. walkers
            std::vector<std::vector<double> > walkers;
            std::vector<double> lnp(nwalkers);

            emcee::init_walkers_from_ball(nwalkers, walkers, pos0, 0.1);
            for (size_t k = 0; k < nwalkers; ++k) { lnp[k] = lnpfunc(walkers[k]); }

            // init. accept vector
            std::vector<bool> accept(nwalkers);
            for (size_t k = 0; k < nwalkers; ++k) { accept[k] = false; }

            // run a production chain
            for(size_t step=0; step < nsteps; ++step){
                emcee::step(walkers, lnp, accept, lnpfunc);
                // convenient cout shortcut.
                emcee::cout_walkers_lnp(walkers, lnp);
            }
            return 0;
        }


.. code:: bash

        > g++ -std=c++11 -O3 emcee.cc example.cc -o example
        > ./example


Note the use of `std::bind <http://en.cppreference.com/w/cpp/utility/functional/bind>`_ 
to allow us the use of partial functions.


Check the examples for more. (the Makefile will compile them for you)

.. _emcee: http://dan.iel.fm/emcee/
