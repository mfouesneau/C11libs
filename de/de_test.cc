/*-----------------------------------------------------------------------------
 *  Test and examples for DE
 *-----------------------------------------------------------------------------*/

#include "de.hh"

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  rosenbrock_fn
 *  Description:  Rosenbrock function with global minimum 
 *                   x = [a, ..., a], and f(x) = 0
 *                non-convex function used as a performance test problem for
 *                optimization algorithms introduced by Howard H. Rosenbrock
 *                in 1960
 *
 * @param x     vector position (n, )
 * @return val  f(x)
 * =====================================================================================
 */
double rosenbrock_fn(std::vector<double>& x){
    size_t ndim = x.size();
    double a = 1.;
    double b = 100.;
    double val = 0;
    for (size_t dim=0; dim < ndim - 1; ++dim){
        val += b * (x[dim + 1] - x[dim] * x[dim]) * (x[dim + 1] - x[dim] * x[dim]) +
            (a - x[dim]) * (a - x[dim]);
    }
    return val;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  ackley_fn
 *  Description:  Ackley function with global minimum at x = [0, ..., 0], f(x) = 0
 * =====================================================================================
 */
double ackley_fn(std::vector<double>& x){
    size_t ndim = x.size();
    double invndim = 1. / (double) ndim;

    double sumx2 =0;
    double sumcosx = 0;
    double pi = std::atan(1) * 4.; 
    for (size_t i = 0; i < ndim; ++i) {
        sumx2 += x[i] * x[i];
        sumcosx += std::cos(2. * pi * x[i]);
    }
    return ( -20. * std::exp(-0.2  * std::sqrt(invndim * sumx2)) 
            - std::exp(invndim * sumcosx) + 20. + std::exp(1));
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  styblinsky_fn
 *  Description:  Styblinski & Tang function 
                  global minimum at x = [-2.903534, ..., -2.903534], and
                  f(x) = -39.16599 * len(n)
 * =====================================================================================
 */
double styblinsky_fn(std::vector<double>& x){
    size_t ndim = x.size();
    double val = 0.;
    for (size_t i = 0; i < ndim; ++i) {
        val += 0.5 * (std::pow(x[i], 4) - 16. * x[i] * x[i] + 5 * x[i]);
    }
    return val;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Main
 *  Description:  
 * =====================================================================================
 */
int main(int argc, char *argv[])
{
    size_t ndim = 3;
    size_t npop = 200;
    size_t nsteps = 2000;
    
    // define function to optimize
    std::function<double (std::vector<double>&)> model;
    model = rosenbrock_fn;

    // limits of parameter space
    std::vector<std::vector<double> > bounds;
    for (size_t dim=0; dim<ndim; ++dim){
            std::vector<double> b = {-5, 5};
            bounds.push_back(b);
    }
    
    std::cout << "Rosenbrock: expecting [1, ..., 1], f(x) = 0" << std::endl;
    optimize(model, bounds, ndim, 10 * npop, nsteps);
    std::cout << "Ackley: expecting [0, ..., 0], f(x) = 0" << std::endl;
    model = ackley_fn;
    optimize(model, bounds, ndim, npop, nsteps);
    std::cout << "Styblinski: expecting [-2.903534, ..., -2.903534], f(x) = -39.16599 * len(n)" 
              << std::endl;
    model = styblinsky_fn;
    optimize(model, bounds, ndim, npop, nsteps);

    return 0;
}

// vim: expandtab:ts=4:softtabstop=4:shiftwidth=4
