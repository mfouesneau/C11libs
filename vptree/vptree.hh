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
#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>
#include <limits>
#include <cmath>
#include <chrono>


/** L2 norm distance.
 *
 * L2 = sqrt (sum_k  (p1[k] - p2[k])^2 )
 *
 * @param p1 first vector
 * @param p2 second vector
 * @return distance between p1 and p2
 */
template<typename T>
T L2(const std::vector<T>& p1, const std::vector<T>& p2 ) {
	size_t ndim = p1.size();
	T d2 = 0;
	for (size_t i = 0; i < ndim; ++i) {
		d2 += (p1[i] - p2[i]) * (p1[i] - p2[i]);
	}
    return sqrt(d2);
}


/** L1 norm distance.
 *
 * L1 = sum_k  | p1[k] - p2[k] |
 *
 * @param p1 first vector
 * @param p2 second vector
 * @return distance between p1 and p2
 */
template<typename T>
T L1(const std::vector<T>& p1, const std::vector<T>& p2 ) {
	size_t ndim = p1.size();
	T d1 = 0;
	for (size_t i = 0; i < ndim; ++i) {
		d1 += std::abs(p1[i] - p2[i]);
	}
    return d1;
}


/** Node structure */
struct Node {
	size_t index;          //< data index number
	double threshold;
	Node* left;
	Node* right;
	// constructor
	Node() : index(0), threshold(0.), left(0), right(0) {}
	// destructor
	~Node() { delete left; delete right; }
};

/** queue Item */
struct HeapItem {
	int index;
	double dist;
	// constructor
	HeapItem( int index, double dist) : index(index), dist(dist) {}
	// comparator
	bool operator<( const HeapItem& o ) const { return dist < o.dist;}
};

/** Distance comparator on items.
 *
 * Allow to compare generic items with a given distance function
 */
template<typename T, double (*distance)( const T&, const T& )>
struct DistanceComparator {
	const T& item;
	// constructor
	DistanceComparator( const T& item ) : item(item) {}
	// Comparator
	bool operator()(const T& a, const T& b) {
	    return distance( item, a ) < distance( item, b );
	}
};


/**
 * Main tree structure
 */
template <typename T, double (*distance)( const T&, const T&) >
class VPTree {

public:
	// Constructor
	VPTree() {this->_root = new Node;};
	// destructor
	~VPTree() {delete this->_root;};
	// initialize the tree from items of type T
	void create(const std::vector<T>& items);
	// search k nearest items
	void search(const T& target, size_t k, 
			std::vector<T>& results, 
			std::vector<double>& distances);
	void search_ball(const T& target, 
		double tau, 
		std::vector<T>& results, 
		std::vector<double>& distances);

private:
	Node* _root;             //< primary node
	double _tau;             //< Ball size (lower): distance to neighbors
	std::vector<T> _items;   //< items in the tree

	// Recursively building the tree nodes
	Node* buildFromPoints(size_t lower, size_t upper );
	// Recursively search for a target
	void  search( Node* node, const T& target, int k, 
			std::priority_queue<HeapItem>& heap );
	
		
};

/** Create the tree from a vector of items.
 *
 * @param items values to index
 */
template<typename T, double (*distance)( const T&, const T& )>
void  VPTree<T, distance>::create(const std::vector<T>& items ) {
        delete _root;
        _items = items;
        _root = buildFromPoints(0, items.size());
    }


/** Build recursively the tree.
 *
 * @param lower  lower index value
 * @param upper  upper index value 
 *
 * @return new Node for the tree
 */
template<typename T, double (*distance)( const T&, const T& )>
Node* VPTree<T, distance>::buildFromPoints(size_t lower, size_t upper ){
	if ( upper == lower ) {
	    return NULL;
	}

	Node* node = new Node();
	node->index = lower;

	if ( upper - lower > 1 ) {

	    // Select the vantage point
	    // choose an arbitrary point and move it to the start
	    size_t i = (size_t)((double)rand() / RAND_MAX * (upper - lower - 1) ) + lower;
	    // swap item places
	    std::swap( this->_items[lower], this->_items[i] );

	    size_t median = ( upper + lower ) / 2;

	    // partitian around the median distance
	    std::nth_element( 
		this->_items.begin() + lower + 1, 
		this->_items.begin() + median,
		this->_items.begin() + upper,
		DistanceComparator<T, distance>( this->_items[lower] ));

	    // what was the median?
	    node->threshold = distance( this->_items[lower], this->_items[median] );

	    node->index = lower;
	    node->left  = this->buildFromPoints( lower + 1, median );
	    node->right = this->buildFromPoints( median, upper );
	}
	return node;
}

/** Search for k nearest items.
 *
 * @param target     object of reference
 * @param k          number of neighbors to extract
 * @param results    vector of items that were found
 * @param distances  distance between results and target
 */
template<typename T, double (*distance)( const T&, const T& )>
void  VPTree<T, distance>::search(const T& target, 
		size_t k, 
		std::vector<T>& results, 
		std::vector<double>& distances) 
{
	std::priority_queue<HeapItem> heap;

	this->_tau = std::numeric_limits<double>::max();
	// search from the root node
	search(_root, target, k, heap );

	// Clean variables
	results.clear(); distances.clear();

	// set results
	while( !heap.empty() ) {
	    results.push_back( _items[heap.top().index] );
	    distances.push_back( heap.top().dist );
	    heap.pop();
	}

	std::reverse( results.begin(), results.end() );
	std::reverse( distances.begin(), distances.end() );
}


/** Search for nearest items within a distance.
 *
 * @param target     object of reference
 * @param tau        maximum distance from target
 * @param results    vector of items that were found
 * @param distances  distance between results and target
 */
template<typename T, double (*distance)( const T&, const T& )>
void  VPTree<T, distance>::search_ball(const T& target, 
		double tau, 
		std::vector<T>& results, 
		std::vector<double>& distances) 
{
	// TODO
}

/**
 * Recursive search between nodes
 *
 * @param node       node to search in
 * @param target     object of reference
 * @param k          number of neighbors to extract
 * @param heap       queue to put new items into
 */
template<typename T, double (*distance)( const T&, const T& )>
void  VPTree<T, distance>::search( Node* node, 
		const T& target, int k, 
		std::priority_queue<HeapItem>& heap )
{
	if ( node == NULL ) return;

	double dist = distance( _items[node->index], target );

	if ( dist < _tau ) {
	    if ( heap.size() == k ) heap.pop();
	    heap.push( HeapItem(node->index, dist) );
	    if ( heap.size() == k ) _tau = heap.top().dist;
	}

	if ( node->left == NULL && node->right == NULL ) {
	    return;
	}

	if ( dist < node->threshold ) {
	    if ( dist - _tau <= node->threshold ) {
		search( node->left, target, k, heap );
	    }

	    if ( dist + _tau >= node->threshold ) {
		search( node->right, target, k, heap );
	    }

	} else {
	    if ( dist + _tau >= node->threshold ) {
		search( node->right, target, k, heap );
	    }

	    if ( dist - _tau <= node->threshold ) {
		search( node->left, target, k, heap );
	    }
	}
}
// vim: expandtab:ts=4:softtabstop=4:shiftwidth=4
