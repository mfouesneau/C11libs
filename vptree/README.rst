A VP-Tree implementation
========================

VP-Tree is similar to that a k-d tree, but uses circular (or hyperspherical)
rather than rectilinear partitioning scheme. 

This implementation is based on Peter Yianilos (1993) `[1]`_ "Data Structures and
Algorithms for Nearest Neighbor Search" and Steve Hanovs's blog post 
`VP trees - A data structure for finding stuff fast`_

Also inspired by blog post
http://stevehanov.ca/blog/index.php?id=130

A `vantage-point tree` (or `VP-tree`) is a tree that segregates data in a
metric space by choosing a position in the space (the *vantage point*) and
partitioning the data points into two parts: those points that are nearer to the
vantage point than a threshold, and those points that are further away. 
By recursively applying this procedure to partition the data into smaller and
smaller sets, a tree data structure is created where neighbors in the tree are
likely to be neighbors in the space.




With high dimensional data, the benefits of the K-D tree are soon lost. As
the number of dimensions increase, the points tend to scatter and it becomes
difficult to pick a good splitting dimension. 

The authors of "Data Mining: Practical machine Learning Tools and Techniques"
suggests using Ball Trees. Each node of a Ball tree describes a bounding
sphere, using a centre and a radius. To make the search efficient, the nodes
should use the minimal sphere that completely contains all of its children,
and overlaps the least with other sibling spheres in the tree.

References
----------

.. _[1]: Peter Yianilos (1993) http://pnylab.com/pny/papers/vptree/vptree.pdf

.. _VP trees - A data structure for finding stuff fast: http://stevehanov.ca/blog/index.php?id=130



Basic Usage
-----------


.. code:: cpp


        // This example is from `vptree.cc`
        #include "vptree.hh"

        /** Example of item that keeps also an index number of reference.
         * Shows how to keep meta data
         */
        struct IndexedVector {
                std::vector<double> x;
                size_t index;
                IndexedVector(size_t index, std::vector<double> x) :  x(x), index(index) {};
        };

        /* Corresponding L1, L2 norms, original functions are in vptree */
        double L1(const IndexedVector& p1, const IndexedVector& p2 ) { return L1(p1.x, p2.x); }
        double L2(const IndexedVector& p1, const IndexedVector& p2 ) { return L2(p1.x, p2.x); }


        int main( int argc, char* argv[] ) {

                // generate some random data
                std::vector<IndexedVector> X;
                size_t N = 1e7;   // how many data points
                size_t ndim = 5;  // data dimensionality

                for (size_t i = 0; i < N; ++i) {
                        std::vector<double> x(ndim);
                        for (size_t dim = 0;  dim < ndim; ++dim) {
                                 x[dim] = (double)rand() / RAND_MAX;
                        }
                        IndexedVector v (i, x);
                        X.push_back(v);
                }

                // create tree: takes about 20 seconds 
                // for N=1e7, ndim=5
                VPTree<IndexedVector, L2> tree;
                tree.create(X);

                // search takes a few micro-seconds
                std::vector<IndexedVector> results;
                std::vector<double> distances;
                tree.search(X[0], 5, results, distances);

                // show results
                std::cout << "query: " << std::endl << " x = ";
                for (auto v: X[0].x){ std::cout << v << " "; }
                std::cout << std::endl;
                for (size_t i = 0; i < results.size(); ++i) {
                        auto x = results[i];
                        std::cout << x.index << " | ";
                        for (auto v: x.x){ std::cout << v << " "; }
                        std::cout << " | " << distances[i] << std::endl;
                }
        }
