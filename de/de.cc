/**----------------------------------------------------------------------------
 * @file de.cc
 *
 * Differential evolution optimizer 
 * ================================
 *
 * Differential evolution is a stochastic population based method that is useful
 * for global optimization problems. At each pass through the population the
 * algorithm mutates each candidate solution by mixing with other candidate
 * solutions to create a trial candidate. There are several strategies for
 * creating trial candidates, which suit some problems more than others. The
 * 'best1bin' strategy is a good starting point for many systems. In this
 * strategy two members of the population are randomly chosen. Their difference
 * is used to mutate the best member (the `best` in `best1bin`), :math:`b_0`, so
 * far:
 *
 * .. math::
 *
 *     b' = b_0 + mutation * (population[rand0] - population[rand1])
 *
 * A trial vector is then constructed. Starting with a randomly chosen 'i'th
 * parameter the trial is sequentially filled (in modulo) with parameters from
 * `b'` or the original candidate. The choice of whether to use `b'` or the
 * original candidate is made with a binomial distribution (the 'bin' in
 * 'best1bin') - a random number in [0, 1) is generated.  If this number is less
 * than the `recombination` constant then the parameter is loaded from `b'`,
 * otherwise it is loaded from the original candidate.  The final parameter is
 * always loaded from `b'`.  Once the trial candidate is built its fitness is
 * assessed. If the trial is better than the original candidate then it takes
 * its place. If it is also better than the best overall candidate it also
 * replaces that.
 *
 * To improve your chances of finding a global minimum use higher `popsize`
 * values, with higher `mutation` and (dithering), but lower `recombination`
 * values. This has the effect of widening the search radius, but slowing
 * convergence.
 *
 * References
 *   ----------
 *  .. [1] Storn, R and Price, K, Differential Evolution - a Simple and
 *         Efficient Heuristic for Global Optimization over Continuous Spaces,
 *         Journal of Global Optimization, 1997, 11, 341 - 359.
 *  .. [2] http://www1.icsi.berkeley.edu/~storn/code.html
 *  .. [3] http://en.wikipedia.org/wiki/Differential_evolution
 *
 *-----------------------------------------------------------------------------*/
 
#include "de.hh"

/**
 * ===  FUNCTION  ======================================================================
 *         Name:  initialize_population
 *  Description:  Random initialization of a DE population uniformly within the
 *                boundaries
 *
 * @param gen 	        Random generator
 * @param ndim          number of dimensions to explore
 * @param npop          size of the population
 * @param population    population vectors
 * @param bounds        vector of boundaries (a, b) for each dimension
 * =====================================================================================
 */
void initialize_population(std::mt19937& gen, 
                size_t ndim,
                size_t npop, 
                std::vector<std::vector<double> >& population,
                std::vector<std::vector<double> >& bounds)
{
    std::uniform_real_distribution<double> randu(0, 1);

    for (size_t pop=0; pop<npop; ++pop){
            std::vector<double> pos(ndim);
            for (size_t dim=0; dim<ndim; ++dim) {
                double bl = bounds[dim][0];
                double bw = bounds[dim][1];
                pos[dim] = bl + randu(gen) * bw;
            }
            population.push_back(pos);
    }
}


/**
 * ===  FUNCTION  ======================================================================
 *         Name:  initialize_population_lhs
 *  Description:  Initializes the population with Latin Hypercube Sampling.
 *                Latin Hypercube Sampling ensures that each parameter is
 *                uniformly sampled over its range.
 * =====================================================================================
 */
void initialize_population_lhs(std::mt19937& gen, 
                size_t ndim,
                size_t npop, 
                std::vector<std::vector<double> >& population,
                std::vector<std::vector<double> >& bounds){

    // make population vector
    for (size_t pop=0; pop<npop; ++pop){
            std::vector<double> pos(ndim);
            population.push_back(pos);
    }
    // default population initialization is a latin hypercube design, but
    // there are other population initializations possible.
    double num_population_members = npop * ndim;
    // Each parameter range needs to be sampled uniformly. The scaled  parameter
    // range ([0, 1)) needs to be split into  `num_population_members`
    // segments, each of which has the following size:
    double segsize = 1.0 / num_population_members;
    // Within each segment we sample from a uniform random distribution.
    // We need to do this sampling for each parameter.
    std::uniform_real_distribution<double> randu(0, 1);
    for (size_t i = 0; i < npop; ++i) {
        for (size_t dim = 0;  dim < ndim; ++dim) {
           population[i][dim] = (bounds[dim][1] - bounds[dim][0]) * (segsize * randu(gen) + i); 
        }
    }
    //Initialize population of candidate solutions by permutation of the random
    //samples.
    std::vector<double> v(npop);
    for (size_t dim = 0;  dim < ndim; ++dim) {
        v.clear();
        for (size_t i = 0; i < npop; ++i) {
            v.push_back(population[i][dim]);
        }
        std::random_shuffle(v.begin(), v.end());
        for (size_t i = 0; i < npop; ++i) {
            population[i][dim] = v[i];
        }
    }
}

/**
 * ===  FUNCTION  ======================================================================
 *         Name:  print_population
 *  Description:  Display on STDOUT the population positions
 *
 * @param population    population vectors
 * @param id            First column content (useful to indicate DE stage)
 * =====================================================================================
 */
void print_population(std::vector< std::vector<double> >& pos, std::string id){

	using namespace std;

    size_t npop = pos.size();
    size_t ndim = pos[0].size();
	
	for (size_t k=0; k<npop; ++k){
		cout << setw(12) << id;
		for (size_t i=0; i<ndim; ++i){
			cout << setw(18) << pos[k][i] << "  ";
		}
        cout << endl;
	}
}


/**
 * ===  FUNCTION  ======================================================================
 *         Name:  print_population
 *  Description:  Display on STDOUT the population positions
 *
 * @param population    population vectors
 * @param fitness       fitness vector
 * @param id            First column content (useful to indicate DE stage)
 * =====================================================================================
 */
void print_population(std::vector< std::vector<double> >& pos,
        std::vector<double>& fitness,
		std::string id){

	using namespace std;

    size_t npop = pos.size();
    size_t ndim = pos[0].size();
	
	for (size_t k=0; k<npop; ++k){
		cout << setw(12) << id;
		for (size_t i=0; i<ndim; ++i){
			cout << setw(18) << pos[k][i] << "  ";
		}
        cout << " | " << setw(18) << fitness[k] << endl;
	}
}


/**
 * ===  FUNCTION  ======================================================================
 *         Name:  de_step
 *  Description:  Make a step in DE optimization algorithm
 *
 * @param fun           function to optimize
 * @param ndim          number of dimensions to explore
 * @param npop          size of the population
 * @param population    population vectors
 * @param fitness       population fitness array (function values)
 * @param bounds        boundaries of the parameter space to explore
 * @param generator     random generator instance
 * @param F             the difference amplification factor.
 *                      Values of 0.5-0.8 are good in most cases.
 * @param C             The cross-over probability. 
 *                      use 0.9 to test for fast convergence, and smaller (~0.1)
 *                      for elaborated search.
 * @param maximize      set to maximize fun instead of minimize.
 * =====================================================================================
 */
void de_step(std::function<double (std::vector<double>&)>& fun,
             std::vector<std::vector<double> >& population,
             std::vector<double>& fitness,
             std::vector<std::vector<double> >& bounds,
             std::mt19937& generator,
             double F, 
             double C,
             bool maximize){

     size_t npop = population.size();
     size_t ndim = population[0].size();

    std::uniform_int_distribution<size_t> randint(0, npop - 1);
    std::uniform_int_distribution<size_t> randdim(0, ndim - 1);
    std::uniform_real_distribution<double> randu(0, 1);

    for (size_t i=0; i<npop; ++i){
        // find population vectors randomly
        size_t t0=i, t1=i, t2=i;
         while (t0 == i)
            t0 = randint(generator);
        while (t1 == i)
            t1 = randint(generator);
        while (t2 == i)
            t2 = randint(generator);

        std::vector<double> pt0 = population[t0];
        std::vector<double> pt1 = population[t1];
        std::vector<double> pt2 = population[t2];

        std::vector<double> v(ndim);
        for (size_t dim=0; dim < ndim; ++dim){
            if (randu(generator) > C){
                // crossover 
                v[dim] = population[i][dim];
            } else {
                v[dim] = pt0[dim] + F * (pt1[dim] - pt2[dim]);
            }
        }
        // shuffle dimensions
        std::vector<double> u(ndim);
        for (size_t dim=0; dim < ndim; ++dim){
            size_t ri = randdim(generator);
            u[ri] = v[ri];
        }

        //compute new fitness
        double ufit = fun(u);
        if (maximize){
                ufit *= -1;
        }

        // acceptance criterion
        if (ufit < fitness[i]){
            for (size_t dim=0; dim < ndim; ++dim){
                population[i][dim] = u[dim];
            }
            fitness[i] = ufit;
        }
    }
}


/**
 * ===  FUNCTION  ======================================================================
 *         Name:  get_result
 *  Description:  Get result from optimization
 *
 *  @param population   population vectors
 *  @param fitness      fitness vector
 *  @return result      optimization_result structure
 * =====================================================================================
 */
struct optimization_result get_result(
        std::vector<std::vector<double> >& population,
        std::vector<double>& fitness){

    // best fit 
    size_t npop = population.size();
    size_t ndim = population[0].size();
    struct optimization_result result;
    result.x = population[0];
    result.fitness = fitness[0];
    result.ndim = ndim;
    for (size_t i = 1; i < npop; ++i) {
        if (fitness[i] < result.fitness) {
            result.x = population[i];
            result.fitness = fitness[i];
        }
    }
    result.fitnesses = fitness;
    result.population = population;
    return result;
}


/**
 * ===  FUNCTION  ======================================================================
 *         Name:  optimization_convergence
 *  Description:  calculate the convergence status of the population
 *
 *  @param fitness  energies of the population
 *  @return val     deviation over mean ratio
 * =====================================================================================
 */
double optimization_convergence(std::vector<double>& fitness){
    size_t count = 0;
    double mean = 0.;
    for (size_t i=0; i< fitness.size(); ++i){
        if (!std::isnan(fitness[i])){
            mean += fitness[i];
            count++;
        }
    }
    mean = mean / (double) count;

    double var = 0.;
    count = 0;
    for (size_t i=0; i<fitness.size(); ++i){
        if (!std::isnan(fitness[i])){
            var += (fitness[i] - mean) * (fitness[i] - mean);
            count++;
        }
    }
    var = var / (double) (count - 1);

    double convergence = std::sqrt(var) / (1e-15 + (std::abs(mean)));
    return convergence;
}


/**
 * ===  FUNCTION  ======================================================================
 *         Name:  optimize
 *  Description:  Optimize function
 *
 *  @param model        function to optimize
 *  @param bounds       boundary of the parameter space ([a,b], ...[a,b])
 *  @param ndim         number of dimensions of the parameter space
 *  @param npop         size of the population
 *  @param nsteps       number of steps during optimization
 *  @param generator    random generator
 *  @param F            the difference amplification factor.
 *                      Values of 0.5-0.8 are good in most cases.
 *  @param C            The cross-over probability. 
 *                      use 0.9 to test for fast convergence, and smaller (~0.1)
 *                      for elaborated search.
 *  @param maximize     set to maximize fun instead of minimize.
 *  @param tol          convergence is considered if the fitness relative std < tol
 *  @return result      result of the optimization
 * =====================================================================================
 */
struct optimization_result optimize(
        std::function<double (std::vector<double>&)>& model,
        std::vector< std::vector<double> >& bounds, 
        size_t ndim, size_t npop, size_t nsteps, 
        std::mt19937& generator,
        double F, 
        double C,
        bool maximize,
        double tol){

    // initialize population
    std::vector<std::vector<double> > population;
    // initialize_population(generator, ndim, npop, population, bounds);
    initialize_population_lhs(generator, ndim, npop, population, bounds);


    // initialize fitness
    std::vector<double> fitness(npop);
    for (size_t k=0; k<npop; ++k){
        fitness[k] = model(population[k]);
    }

    // print_population(population, fitness, "0");

    // optimize
    size_t number_of_steps = nsteps;
    double convergence = 0;
    for (size_t step=0; step < nsteps; ++step){
        de_step(model, population, fitness, bounds, generator);
        convergence = optimization_convergence(fitness);
        // stop when the fractional s.d. of the population is less than tol of
        // the mean energy
        if (convergence < tol){
            number_of_steps = step;
            break;
        }
    }

    // best fit 
    struct optimization_result result = get_result(population, fitness);

    std::cout << "Optimization ended" << std::endl;
    std::cout << "  + population of size " << npop << std::endl;
    std::cout << "  + after nsteps: " << number_of_steps << std::endl;
    std::cout << "  + s.d / mean energy pop: " << convergence << std::endl;
    std::cout << "  + Best fit value:" << std::endl;
    for (size_t i = 0; i < result.ndim; ++i) {
        std::cout << std::setw(12) << result.x[i] << " ";
    } 
    std::cout << std::setw(12) << result.fitness << std::endl;

    return result;
}

/**
 * ===  FUNCTION  ======================================================================
 *         Name:  optimize
 *  Description:  Optimize function
 *
 *  @param model        function to optimize
 *  @param bounds       boundary of the parameter space ([a,b], ...[a,b])
 *  @param ndim         number of dimensions of the parameter space
 *  @param npop         size of the population
 *  @param nsteps       number of steps during optimization
 *  @param F            the difference amplification factor.
 *                      Values of 0.5-0.8 are good in most cases.
 *  @param C            The cross-over probability. 
 *                      use 0.9 to test for fast convergence, and smaller (~0.1)
 *                      for elaborated search.
 *  @param maximize     set to maximize fun instead of minimize.
 *  @return result      result of the optimization
 * =====================================================================================
 */
struct optimization_result optimize(
        std::function<double (std::vector<double>&)>& model,
        std::vector< std::vector<double> >& bounds, 
        size_t ndim, size_t npop, size_t nsteps,
        double F, 
        double C,
        bool maximize){

    // define random generator
    std::random_device rd;
    std::mt19937 generator(rd());

    return optimize(model, bounds, ndim, npop, nsteps,
            generator);
}


// vim: expandtab:ts=4:softtabstop=4:shiftwidth=4
