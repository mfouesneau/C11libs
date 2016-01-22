/** ---------------------------------------------------------------------------
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
#ifndef __DE_H__
#define __DE_H__
#include <stdlib.h>
#include <random>
#include <iostream>
#include <iomanip>
#include <vector>

/** 
 *  best position x and value fitness.
 */
struct optimization_result {
    size_t ndim;                                  /**< number of dimensions in parameter space*/
    std::vector<double> x;                        /**< best fit value */
    double fitness;                               /**< energy of the best fit */
    std::vector<std::vector<double> > population; /**< final population positions*/
    std::vector<double> fitnesses;                /**< energy of each positions  */
};

/**
 * Make a step in DE optimization algorithm
 *
 *  @param fun           function to optimize
 *  @param population    population vectors
 *  @param fitness       population fitness array (function values)
 *  @param bounds        boundaries of the parameter space to explore
 *  @param generator     random generator instance
 *  @param F             the difference amplification factor.
 *                       Values of 0.5-0.8 are good in most cases.
 *  @param C             The cross-over probability. 
 *                       use 0.9 to test for fast convergence, and smaller (~0.1)
 *                       for elaborated search.
 *  @param maximize      set to maximize fun instead of minimize.
 */
void de_step(std::function<double (std::vector<double>&)>& fun,
             std::vector<std::vector<double> >& population,
             std::vector<double>& fitness,
             std::vector<std::vector<double> >& bounds,
             std::mt19937& generator,
             double F=0.5, 
             double C=0.5,
             bool maximize=false);

/**
 * Make a step in DE optimization algorithm
 *
 * @param fun           function to optimize
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
 */
void de_step(std::function<double (std::vector<double>&)>& fun,
             std::vector<std::vector<double> >& population,
             std::vector<double>& fitness,
             std::vector<std::vector<double> >& bounds,
             double F=0.5, 
             double C=0.5,
             bool maximize=false);

/**
 *  Get result from optimization
 *
 *  @param population   population vectors
 *  @param fitness      fitness vector
 *  @return result      optimization_result structure
 */
struct optimization_result get_result(
        std::vector<std::vector<double> >& population,
        std::vector<double>& fitness);
/**
 * calculate the convergence status of the population
 *
 *  @param fitness  energies of the population
 *  @return val     deviation over mean ratio
 */
double optimization_convergence(std::vector<double>& fitness);
/**
 * Random initialization of a DE population uniformly within the boundaries
 *
 * @param gen 	        Random generator
 * @param ndim          number of dimensions to explore
 * @param npop          size of the population
 * @param population    population vectors
 * @param bounds        vector of boundaries (a, b) for each dimension
 */
void initialize_population(std::mt19937& gen, 
                size_t ndim,
                size_t npop, 
                std::vector<std::vector<double> >& population,
                std::vector<std::vector<double> >& bounds);
/**
 * Initializes the population with Latin Hypercube Sampling. 
 *
 * Latin Hypercube Sampling ensures that each parameter is uniformly sampled
 * over its range.
 */
void initialize_population_lhs(std::mt19937& gen, 
                size_t ndim,
                size_t npop, 
                std::vector<std::vector<double> >& population,
                std::vector<std::vector<double> >& bounds);
/**
 *  Optimize function.
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
 */
struct optimization_result optimize(
        std::function<double (std::vector<double>&)>& model,
        std::vector< std::vector<double> >& bounds, 
        size_t ndim, size_t npop, size_t nsteps, 
        std::mt19937& generator,
        double F=0.5, 
        double C=0.5,
        bool maximize=false,
        double tol=1e-3);
/**
 *  Optimize function
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
 */
struct optimization_result optimize(
        std::function<double (std::vector<double>&)>& model,
        std::vector< std::vector<double> >& bounds, 
        size_t ndim, size_t npop, size_t nsteps,
        double F=0.5, 
        double C=0.5,
        bool maximize=false);
/**
 * Display on STDOUT the population positions
 *
 * @param population    population vectors
 * @param id            First column content (useful to indicate DE stage)
 */
void print_population(std::vector< std::vector<double> >& pos, std::string id="0");
/**
 *  Display on STDOUT the population positions
 *
 * @param population    population vectors
 * @param fitness       fitness vector
 * @param id            First column content (useful to indicate DE stage)
 */
void print_population(std::vector< std::vector<double> >& pos,
        std::vector<double>& fitness, std::string id="0");

#endif
