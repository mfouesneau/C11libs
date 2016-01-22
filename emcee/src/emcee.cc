#include "emcee.hh"

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
        ){

        // random generator
        std::random_device rd;
        std::mt19937 gen(rd());
        // 
        std::uniform_real_distribution<double> randu(0, 1);
        std::uniform_int_distribution<size_t> rand(0, nwalkers - 1);

        std::vector<double> q(ndim);

        for (size_t k = 0; k < nwalkers; ++k) {
            // compute the random stretch factor 
            double z = (a - 1) * randu(gen) + 1;
            z = z * z / a;

            // select the helper walker
            size_t ri = rand(gen);
            if (ri >= k){
                for (size_t dim=0; dim < ndim; ++dim)
                    q[dim] = pin[ri][dim];
            } else {
                q = pout[ri];
                for (size_t dim=0; dim < ndim; ++dim)
                    q[dim] = pout[ri][dim];
            }
            
            // compute the proposal
            for (size_t dim=0; dim < ndim; ++dim){
                q[dim] = (1 - z) * q[dim] + z * pin[k][dim];
            }
            double lp = lnpfunc(q);
            double diff = (ndim - 1) * std::log(z) + lp - lpin[k];

            // Accept or reject
            if (diff >= 0){
                accept[k] = true;
            } else {
                if( diff >= std::log(randu(gen))){
                    accept[k] = true;
                } else {
                    accept[k] = false;
                }
            }

            // update
            if (accept[k] == true){
                for (size_t dim=0; dim < ndim; ++dim)
                    pout[k][dim] = q[dim];
                lpout[k] = lp;
            } else {
                for (size_t dim=0; dim < ndim; ++dim)
                    pout[k][dim] = pin[k][dim];
                lpout[k] = lpin[k];
            }
            
        }

    }

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
        size_t nstep,
        double a){

        size_t ndim = walkers[0].size();
        size_t nwalkers = walkers.size();

        for( size_t step=0; step < nstep ; ++ step){
            advance(ndim, nwalkers, a, 
                    walkers, lnp, 
                    walkers, lnp, 
                    accept, lnpfunc);
        }
    }

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
        double radius=1){
        
        size_t ndim = pos0.size();

        // define random generator
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> randn(0, radius);

        // init. walkers
        for (size_t k=0; k < nwalkers; ++k){
            std::vector<double> pos(ndim);
            for (size_t dim = 0; dim < ndim; ++dim) {
                pos[dim] = randn(gen) + pos0[dim];
            }
            walkers.push_back(pos);
        }
    }

    /** convenient function that computes the acceptance ratio given the accept
     * vector
     *
     *  @param accept  accept vector [nwalkers]
     *
     *  @return val   acceptance ratio
     *
     */

    float compute_acceptance_ratio(std::vector<bool>& accept){
        float a_ratio = 0;
        size_t nwalkers = accept.size();

        for(size_t k=0; k < nwalkers; ++k){
            a_ratio += static_cast<float>(accept[k]) / static_cast<float>(nwalkers);
        }
        return a_ratio;
    }

    /** convenient function that cout walker positions and lnp
     *
     *  @param walkers  starting positions [nwalkers][ndim]
     * 
     *  @param lnp      log-probability at positions `walkers` [nwalkers]
     *
     */
    void cout_walkers_lnp(
        std::vector<std::vector<double> >& walkers,
        std::vector<double>& lnp){

        size_t nwalkers = walkers.size();
        size_t ndim = walkers[0].size();

        for(size_t k=0; k < nwalkers; ++k){
            for (size_t dim = 0; dim < ndim; ++dim)
                std::cout << std::setw(8) << walkers[k][dim] << " ";
            std::cout << lnp[k] << std::endl;
        }
    }
}

// vim: expandtab:ts=4:softtabstop=4:shiftwidth=4
		
