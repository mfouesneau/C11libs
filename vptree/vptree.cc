/** A VP-Tree implementation.
 * Based on "Data Structures and Algorithms for Nearest Neighbor Search" by
 * Peter N. Yianilos [http://pnylab.com/pny/papers/vptree/vptree.pdf]
 *
 * Also inspired by blog post
 * http://stevehanov.ca/blog/index.php?id=130
 *
 * With high dimensional data, the benefits of the K-D tree are soon lost. As
 * the number of dimensions increase, the points tend to scatter and it becomes
 * difficult to pick a good splitting dimension. 
 *
 * The authors of "Data Mining: Practical machine Learning Tools and Techniques"
 * suggests using Ball Trees. Each node of a Ball tree describes a bounding
 * sphere, using a centre and a radius. To make the search efficient, the nodes
 * should use the minimal sphere that completely contains all of its children,
 * and overlaps the least with other sibling spheres in the tree.
 */
#include "vptree.hh"

/** Example of item that keeps also an index number of reference.
 *
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


/** Example usage */
int main( int argc, char* argv[] ) {

	std::vector<IndexedVector> X;
	size_t N = 1e7;
	size_t ndim = 5;

	std::cout << "VPTree example " << std::endl
		  << "   Using a dataset of size " << N << std::endl
		  << "   of dimensionality " << ndim << std::endl;

	// timers
	std::chrono::system_clock::time_point start, end;

	// create data
	start = std::chrono::system_clock::now();
	for (size_t i = 0; i < N; ++i) {
		std::vector<double> x(ndim);
		for (size_t dim = 0;  dim < ndim; ++dim) {
	    		 x[dim] = (double)rand() / RAND_MAX;
		}
		IndexedVector v (i, x);
		X.push_back(v);
	}
	end = std::chrono::system_clock::now();
	std::cout << "Create data took " 
		  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
	          << " seconds"
		  << std::endl;

	// create tree
	start = std::chrono::system_clock::now();
	VPTree<IndexedVector, L2> tree;
	tree.create(X);
	end = std::chrono::system_clock::now();
	std::cout << "Create tree took " 
		  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
	          << " seconds"
		  << std::endl;

	// search
	std::vector<IndexedVector> results;
	std::vector<double> distances;
	// 
	start = std::chrono::system_clock::now();
	tree.search(X[0], 5, results, distances);
	end = std::chrono::system_clock::now();
	std::cout << "Search took "
		  << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
	          << " us"
		  << std::endl;

	// show results
	std::cout << "query: " << std::endl << " x = ";
	for (auto v: X[0].x){
		std::cout << v << " ";
	}
	std::cout << std::endl;
	for (size_t i = 0; i < results.size(); ++i) {
		auto x = results[i];
		std::cout << x.index << " | ";
		for (auto v: x.x){
			std::cout << v << " ";
		}
		std::cout << " | " << distances[i] << std::endl;
	}

}
