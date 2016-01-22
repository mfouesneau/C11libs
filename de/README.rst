Differential evolution optimizer 
================================

Differential evolution is a stochastic population based method that is useful
for global optimization problems. At each pass through the population the
algorithm mutates each candidate solution by mixing with other candidate
solutions to create a trial candidate. There are several strategies for
creating trial candidates, which suit some problems more than others. The
'best1bin' strategy is a good starting point for many systems. In this
strategy two members of the population are randomly chosen. Their difference
is used to mutate the best member (the `best` in `best1bin`), :math:`b_0`, so
far:

.. math::

    b^\prime = b_0 + mutation \times (population[rand_0] - population[rand_1])


A trial vector is then constructed. Starting with a randomly chosen 'i'th
parameter the trial is sequentially filled (in modulo) with parameters from
:math:`b^\prime` or the original candidate. The choice of whether to use
:math:`b^\prime`  or the original candidate is made with a binomial
distribution (the 'bin' in 'best1bin') - a random number in [0, 1) is
generated.  If this number is less than the `recombination` constant then the
parameter is loaded from `b'`, otherwise it is loaded from the original
candidate.  The final parameter is always loaded from `b'`.  Once the trial
candidate is built its fitness is assessed. If the trial is better than the
original candidate then it takes its place. If it is also better than the
best overall candidate it also replaces that.

To improve your chances of finding a global minimum use higher `popsize`
values, with higher `mutation` and (dithering), but lower `recombination`
values. This has the effect of widening the search radius, but slowing
convergence.

References
----------

 .. [1] Storn, R and Price, K, Differential Evolution - a Simple and
        Efficient Heuristic for Global Optimization over Continuous Spaces,
        Journal of Global Optimization, 1997, 11, 341 - 359.

 .. [2] http://www1.icsi.berkeley.edu/~storn/code.html

 .. [3] http://en.wikipedia.org/wiki/Differential_evolution


