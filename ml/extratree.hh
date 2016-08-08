/**
 * An extra-trees regressor. 
 * 
 * This class implements a meta-estimator that fits a number of randomized
 * decision trees (a.k.a. extra-trees) on various sub-samples of the dataset and
 * use averaging to improve the predictive accuracy and control over-fitting.
 * 
 * Extra-tree algorithm is a perturb-and-combine technique [B1998] specifically
 * designed for trees. This means a diverse set of classifiers is created by
 * introducing randomness in the classifier construction. The prediction of the
 * ensemble is given as the averaged prediction of the individual classifiers.
 * 
 * In extremely randomized trees thresholds are drawn at random for each
 * candidate feature and the best of these randomly-generated thresholds is
 * picked as the splitting rule. This usually allows to reduce the variance of
 * the model a bit more, at the expense of a slightly greater increase in bias.
 * 
 * The main parameters to adjust when using these methods is `n_estimators` and
 * `max_features`. The former is the number of trees in the forest. The larger
 * the better, but also the longer it will take to compute. In addition, note
 * that results will stop getting significantly better beyond a critical number
 * of trees. The latter is the size of the random subsets of features to
 * consider when splitting a node. The lower the greater the reduction of
 * variance, but also the greater the increase in bias. Empirical good default
 * values are max_features=n_features for regression problems, and
 * max_features=sqrt(n_features) for classification tasks (where n_features is
 * the number of features in the data). Good results are often achieved when
 * setting max_depth=None in combination with min_samples_split=1 (i.e., when
 * fully developing the trees). Bear in mind though that these values are
 * usually not optimal, and might result in models that consume a lot of ram.
 * The best parameter values should always be cross-validated. 
 *
 *
 * The Class uses templating to allow for arbitrary precision (float or double)
 * And can also be used for classification if the output is of interger type
 * (e.g. int,long, size_t)
 *
 * 	extratree::ExtraTreesRegressor<double, size_t> reg(...);
 * or
 * 	extratree::ExtraTreesRegressor<float, short> reg(...);
 *
 * 
 * FIXME: add bootstrapping option during learning
 * FIXME: fix names in the header doc.
 * 
 * References
 * ----------
 * 
 * [B1998]	Breiman, “Arcing Classifiers”, Annals of Statistics 1998.
 * [GEW2006]	P. Geurts, D. Ernst., and L. Wehenkel, “Extremely randomized trees”, 
 *                  Machine Learning, 63(1), 3-42, 2006.
 *
 * @author fouesneau
 */
#pragma once
#include <stdlib.h>
#include <random>
#include <iostream>
#include <vector>
#include <cmath>
#include <memory>

namespace extratrees { 

// std::vector copies the content on storage. 
// Not best in our case, we need to store references/pointers so 
// modifications are persistent.
// 
// vector of shared pointers 
template <typename T>
using sptr_vector = std::vector<std::shared_ptr<T>>;
// Note that we must create node pointers directly when using make_shared as it
// corresponds to calling new.


////////////////////////////////////////////////////////////////////////////////////////
/** Defines a Node for a Regressor Tree.
 *
 * It is mainly a container with nested structure. 
 *
 * @author fouesneau
 *
 */
template <typename T, typename U>
class ExtraTreeNode
{
public:
	ExtraTreeNode(const std::vector<std::vector<T>>& X, const std::vector<U>& y);
	bool hasChildren();
	void clean();
	std::vector<std::vector<T>> getX();
	std::vector<U> gety();
	void setAsRoot();
	bool isRoot();
	void setCutValue(const U cut);
	U getCutValue();
	void setChildrenNodes(const sptr_vector<ExtraTreeNode>& children);
	sptr_vector<ExtraTreeNode> getChildrenNodes();
	bool stopSplit(size_t nmin);
	bool checkAreLabelsConstant();
	bool checkAreAttributesConstant();
	U predict();
	U getScore();
	std::vector<T> getCutAttributes(size_t d);
	void setCutProjection(size_t dimensionOfCut);
	size_t getCutProjection();

private:
	std::vector<std::vector<T>> X;
	std::vector<U> y;
	size_t ndimX = 0;
	size_t n_members = 0;
	U y_mean = 0.;
	T cutvalue = 0;
	bool is_root = false;
	sptr_vector<ExtraTreeNode> children;
	size_t dimensionOfCut = -1;
};


/**
 * Constructor
 * @param X   input features to store
 * @param y   input labels to store
 * @throws RuntimeError when X and y have incompatible shapes
 */
template <typename T, typename U>
ExtraTreeNode<T,U>::ExtraTreeNode(const std::vector<std::vector<T>>& X, 
				  const std::vector<U>& y){
	// Check that number of labels and input vectors are identical.
	if (y.size() != X.size()) {
		throw std::runtime_error("Numbers of input vectors and labels is not identical.");
	}
	this->X = X;
	this->y = y;
	this->n_members = y.size();
	this->ndimX = X[0].size();
	// Calculate mean value of given labels in this->node.
	this->y_mean = 0;
	for (U yk : y){ this->y_mean += yk; }
	this->y_mean /= this->n_members;
}
	
/**
 * Whether the Node has sub-nodes
 * @return  true if has children
 */
template <typename T, typename U>
bool ExtraTreeNode<T,U>::hasChildren(){
	if (this->children.size() > 0)
		return true;
	else
		return false;
}
	
/**
 * Drop references to X, y data.
 */
template <typename T, typename U>
void ExtraTreeNode<T,U>::clean() {
	this->X.clear();
	this->y.clear();
}
	
/**
 * Return input features
 * @return X
 */
template <typename T, typename U>
std::vector<std::vector<T>> ExtraTreeNode<T,U>::getX() {
	return this->X;
}

/**
 * Return input labels
 * @return y
 */
template <typename T, typename U>
std::vector<U> ExtraTreeNode<T,U>::gety() {
	return this->y;
}
	
/**
 * Set root flag
 */
template <typename T, typename U>
void ExtraTreeNode<T,U>::setAsRoot() {
	this->is_root = true;
}
	
	
/**
 * Check if root
 * @return True if root is set.
 */
template <typename T, typename U>
bool ExtraTreeNode<T,U>::isRoot() {
	return this->is_root;
}
	
/**
 * Set the cutting value
 * @param cut  threshold for split
 */
template <typename T, typename U>
void ExtraTreeNode<T,U>::setCutValue(const U cut) {
	this->cutvalue = cut;
}
	
	
/**
 * Get the split value
 * @return cutvalue
 */
template <typename T, typename U>
U ExtraTreeNode<T,U>::getCutValue() {
	return this->cutvalue;
}
	
/**
 * Set children to a pair of nodes
 * @param children  vector of nodes
 */
template <typename T, typename U>
void ExtraTreeNode<T,U>::setChildrenNodes(const sptr_vector<ExtraTreeNode>& children){
	// this->children = children;
	for (auto child : children)
		this->children.push_back(std::shared_ptr<ExtraTreeNode>(child));
}

/**
 * Returns the children nodes (could be null)
 * @return children
 */
template <typename T, typename U>
sptr_vector<ExtraTreeNode<T,U>> ExtraTreeNode<T,U>::getChildrenNodes() {
	return this->children;
}
	
	
/** 
 * Check if the node can split 
 * 
 * stop splitting if at least one of the following conditions is true:
 *  - number of node elements is less than nmin
 *  - all attributes are constant in this node
 *  - all labels are constant in this node
 *  
 * @param nmin  minimum number of object to allow a split
 * @return true if split is not allowed
 */
template <typename T, typename U>
bool ExtraTreeNode<T,U>::stopSplit(size_t nmin) {
	if      (this->n_members < nmin)             return true;
	else if (this->checkAreLabelsConstant())     return true;
	else if (this->checkAreAttributesConstant()) return true;
	else                                         return false;
}
	
/** Check if all labels are constant. 
 * @return true if all y are identical
 */
template <typename T, typename U>
bool ExtraTreeNode<T,U>::checkAreLabelsConstant() {
	bool areConstant = true;
	U yref = this->y[0];
	for (U yk : this->y) 
		if (yk != yref) {
			areConstant = false;
			break;
		}
	return areConstant;
}
	
/** Check if all attributes X are constant. 
 * @return true if all X are identical
 */
template <typename T, typename U>
bool ExtraTreeNode<T,U>::checkAreAttributesConstant() {
	bool areConstant = true;
	std::vector<T> Xref = this->X[0];
	for (std::vector<T> Xk : this->X) 
		for (size_t d=0; d < ndimX; ++d) 
			if (Xk[d] != Xref[d]) {
				areConstant = false;
				break;
			}
	return areConstant;
}

/** Predict according to arithmetic mean of label in this node. 
 * @return y_mean
 */
template <typename T, typename U>
U ExtraTreeNode<T,U>::predict() {
	return this->y_mean;
}
	
/** Get RMS error score of predictions from this node. 
 * @return rms of the labels to its mean.
 */
template <typename T, typename U>
U ExtraTreeNode<T,U>::getScore() {
	U chi2 = 0.0;
	U dev = 0.0;
	for (T yk: this->y){
		dev = yk - this->y_mean;
		chi2 += dev * dev;
	}
	return chi2;
}
	
/** Get range of attribute values in some dimension d. */
template <typename T, typename U>
std::vector<T> ExtraTreeNode<T,U>::getCutAttributes(size_t d) {
	std::vector<T> cutattributes(this->n_members);
	for (size_t n=0; n<n_members; ++n) cutattributes[n] = this->X[n][d];
	return cutattributes;
}

/** Set the dimension of the cut.
 *
 * @param dimensionOfCut dimension
 * */
template <typename T, typename U>
void ExtraTreeNode<T,U>::setCutProjection(size_t dimensionOfCut){
	this->dimensionOfCut = dimensionOfCut;
}

/** Get the dimension of the cut.
 *
 * @return dimensionOfCut dimension
 * */
template <typename T, typename U>
size_t ExtraTreeNode<T,U>::getCutProjection() {
	return this->dimensionOfCut;
}

////////////////////////////////////////////////////////////////////////////////////////

/**
 * Extremely randomized Tree
 * 
 * This class implements one single tree.
 * 
 * @author fouesneau
 */
template <typename T, typename U>
class ExtraTree
{
public:
	ExtraTree(size_t n_min, size_t n_splits);
	ExtraTree(size_t n_min, size_t n_splits, long seed);
	ExtraTree(size_t n_min, size_t n_splits, const std::mt19937& gen);
	void fit(const std::vector<std::vector<T>>& X, const std::vector<U>& y);
	U predict(const std::vector<T>& X);
	sptr_vector<ExtraTreeNode<T,U>> splitNode(const std::shared_ptr<ExtraTreeNode<T,U>>& node, 
					          const std::vector<T>& cutattributes, T cut);
	void clean();


private:
	/* data */
	std::mt19937 gen;                /**< Random generator for that tree */
	size_t DimAttributeSpace = -1;   /**< Dimension of the X/feature space */
	size_t n_min = 1;                /**< Number of minimal data to split a node */
	size_t n_splits = 1;             /**< Number of random split during learning */
	long seed = 1234;                /**< Random generator seed */

	/** Data */
	sptr_vector<ExtraTreeNode<T,U>> allnodes;
	sptr_vector<ExtraTreeNode<T,U>> leafnodes;
};


/**
 * Constructor 
 * @param n_min     Number of minimal data to split a node
 * @param n_splits  Number of random split during learning
 */
template <typename T, typename U>
ExtraTree<T,U>::ExtraTree(size_t n_min, size_t n_splits) {
	this->n_min = n_min;
	this->n_splits = n_splits;
	// Set random seed.
	std::random_device rd;
	std::mt19937 gen(rd());
	this->gen = gen;
}

/**
 * Constructor 
 * @param n_min     Number of minimal data to split a node
 * @param n_splits  Number of random split during learning
 * @param seed      Random generator seed
 */
template <typename T, typename U>
ExtraTree<T,U>::ExtraTree(size_t n_min, size_t n_splits, long seed) {
	this->n_min = n_min;
	this->n_splits = n_splits;
	this->seed = seed;
	// Set random seed.
	std::mt19937 gen(seed);
	this->gen = gen;
}

/**
 * Constructor 
 * @param n_min     Number of minimal data to split a node
 * @param n_splits  Number of random split during learning
 * @param rng      Random generator
 */
template <typename T, typename U>
ExtraTree<T,U>::ExtraTree(size_t n_min, size_t n_splits, const std::mt19937& gen) {
	this->n_min = n_min;
	this->n_splits = n_splits;
	this->seed = -1;
	// Set random seed.
	this->gen = gen;
}


/**
 * Make the tree from a set of Features/Label pairs
 * @param X    Input training samples
 * @param y    Output training samples
 * @throws RuntimeError when X and y have incompatible shapes
 */
template <typename T, typename U>
void ExtraTree<T,U>::fit(const std::vector<std::vector<T>>& X, const std::vector<U>& y) {
	// Validate dimensions
	if (y.size() != X.size()) {
		throw std::runtime_error("Numbers of input vectors and labels is not identical.");
	}
	// Get dimension of attribute space
	this->DimAttributeSpace = X[0].size();

	// Random dimension generator
	std::uniform_int_distribution<size_t> rand_dim(0, this->DimAttributeSpace - 1);

	// Initialize all data as a single node.
	std::shared_ptr<ExtraTreeNode<T,U>> node = std::make_shared<ExtraTreeNode<T,U>>(X, y);
	node->setAsRoot();
	this->allnodes.push_back(node);
	this->leafnodes.push_back(node);

	// Continue to grow until all nodes refuse further splitting.
	bool growing = true;
	while (growing) {

		// initialize empty node vector
		sptr_vector<ExtraTreeNode<T,U>> newleafnodes;

		// Stop growing if no node wants to further split.
		growing = false;

		// Iterate over all nodes and try to split them.
		for (std::shared_ptr<ExtraTreeNode<T,U>> leaf : this->leafnodes){
			// Ask if the current node wants to get split.
			if (leaf->stopSplit(this->n_min)) {
				// If it wants to stop splitting, append the new set of nodes.
				newleafnodes.push_back(leaf);
			} else {
				// This node wants to split.
				growing = true;
				// Select K attributes a_1,...,a_K and their cut values from node.
				size_t bestDimension      = -1;
				U bestScore          = std::numeric_limits<T>::infinity();
				T bestCut            = std::numeric_limits<T>::infinity();
				sptr_vector<ExtraTreeNode<T,U>> bestSplit;
				for (size_t k=0; k < n_splits; ++k) {
					// Randomly select a dimension in attribute space where you will apply the cut.
					size_t d = rand_dim(this->gen);
					// Get cut attributes.
					std::vector<T> cutattributes = leaf->getCutAttributes(d);
					// Randomly draw a cut.
					std::uniform_int_distribution<size_t> rand_length(0, cutattributes.size() - 1);
					size_t i1 = rand_length(gen);
					size_t i2 = rand_length(gen);
					T cut = 0.5 * (cutattributes[i1] + cutattributes[i2]);
					// Try to split the current node according to this cut.
					try {
						sptr_vector<ExtraTreeNode<T,U>> split = this->splitNode(leaf, cutattributes, cut);
						// Get scores of both new nodes.
						U chi2  = split[0]->getScore() + split[1]->getScore();
						// Is this the best split?
						if (chi2 < bestScore) {
							bestScore     = chi2;
							bestDimension = d;
							bestCut       = cut;
							bestSplit     = split;
						}
					} catch (const std::exception& e) {
						continue;
					}
				}

				if (bestSplit.empty() == 1) {  // no split was found, just add it.
					newleafnodes.push_back(leaf);
				} else { // add both nodes to the pile
					// Set cut value and cut dimension of parent node.
					leaf->setCutProjection(bestDimension);
					leaf->setCutValue(bestCut);
					// Set new leaf nodes as children of parent node.
					leaf->setChildrenNodes(bestSplit);
					// Add new leaf nodes.
					newleafnodes.push_back(bestSplit[0]);
					newleafnodes.push_back(bestSplit[1]);
					this->allnodes.push_back(bestSplit[0]);
					this->allnodes.push_back(bestSplit[1]);
				}
			}
		}
		// Check if set of new nodes is larger than set of old nodes.
		if (newleafnodes.size() > this->leafnodes.size()) {
			growing = true;
		} else {
			growing = false;
		}
		// Replace old set of nodes by new.
		this->leafnodes = newleafnodes;
	}
	for (auto _node : this->allnodes){ _node->clean(); }
}


/** 
 * Predict label from Extra Tree given a test vector of attributes.
 * 
 * @param X   Input test data
 * @return    prediction
 * @throws RuntimeError when X has different shape from expected
 */
template <typename T, typename U>
U ExtraTree<T,U>::predict(const std::vector<T>& X) {
	if (X.size() != this->DimAttributeSpace) {
		throw std::runtime_error("Shape of input vector incorrect");
	}
	T ypred = std::numeric_limits<T>::quiet_NaN();
	bool hasChildren = true;
	// Start at root node.
	std::shared_ptr<ExtraTreeNode<T,U>> node = this->allnodes[0];
	while (hasChildren) {
		// Check if the current node still has children.
		hasChildren = node->hasChildren();
		if (!hasChildren) {
			// If it has no children, return its prediction.
			ypred = node->predict();
		} else {
			sptr_vector<ExtraTreeNode<T,U>> children = node->getChildrenNodes();
			// If it does have children, follow its decision to the correct child.
			size_t d     = node->getCutProjection();
			T cut = node->getCutValue();
			// Project attributes.
			if (X[d] < cut) 
				node = children[0];
			else 
				node = children[1];
		}
	}
	return ypred;
}

/**
 * Split a node given the values and threshold for the cut
 * 
 * @param node             Tree node to split
 * @param cutattributes    values from the input features for the cut X[:, d]
 * @param cut              threshold value
 * @return Vector of 2 nodes
 * @throws Exception  Raise an exception if it does not pass cut tests.
 */
template <typename T, typename U>
sptr_vector<ExtraTreeNode<T,U>> ExtraTree<T,U>::splitNode(const std::shared_ptr<ExtraTreeNode<T,U>>& node, 
   		                                        const std::vector<T>& cutattributes, 
						        T cut) {
	// Get attributes and labels of current node.
	std::vector<std::vector<T>> attributes = node->getX();
	std::vector<U> labels = node->gety();
	size_t N = labels.size();
	
	// Count how many examples will go into first node.
	size_t NA = 0;
	for (double value : cutattributes) if (value < cut) NA++;

	// Split attributes and labels.
	size_t NB = N - NA;
	if (NA<1 || NB<1) throw std::runtime_error("Each new node must at least contain one example.");
	std::vector<std::vector<T>> attA;
	std::vector<std::vector<T>> attB;
	std::vector<U> labelA;
	std::vector<U> labelB;
	for (size_t n=0; n<N; ++n) {
		if (cutattributes[n] < cut) {
			attA.push_back(attributes[n]);
			labelA.push_back(labels[n]);
		} else {
			attB.push_back(attributes[n]);
			labelB.push_back(labels[n]);
		}
	}
	// Create new nodes.
	sptr_vector<ExtraTreeNode<T,U>> newnodes = {std::make_shared<ExtraTreeNode<T,U>>(attA, labelA),
						    std::make_shared<ExtraTreeNode<T,U>>(attB, labelB)};
	return newnodes;
}

/**
 * Clean memory from left over data. 
 */
template <typename T, typename U>
void ExtraTree<T,U>::clean() {
	for (auto _node : this->allnodes){
		_node->clean();
	}
	for (auto _node : this->leafnodes){
		_node->clean();
	}
}

////////////////////////////////////////////////////////////////////////////////////////

/**
 * A forest of ExtraTree objects for regression
 * @author fouesneau
 */
template <typename T, typename U>
class ExtraTreesRegressor 
{
public:
	ExtraTreesRegressor(size_t n_estimators, size_t n_min, size_t n_splits);
	ExtraTreesRegressor(size_t n_estimators, size_t n_min, size_t n_splits, long seed);
	void fit(const std::vector<std::vector<T>>& X, const std::vector<U>& y);
	void fit_bootstrap(const std::vector<std::vector<T>>& X, const std::vector<U>& y);
	std::vector<U> predict(const std::vector<T>& X);

private:

	sptr_vector<ExtraTree<T,U>> forest;  /**< set of trees */
	std::mt19937 gen;                  /**< Random generator for that tree */
	size_t n_min = 1;                  /**< Number of minimal data to split a node */
	size_t n_splits = 1;               /**< Number of random split during learning */
	long seed = 1234;                  /**< Random generator seed */
};

/**
 * Constructor 
 * @param n_estimators   number of trees in the forest
 * @param n_min          Minimum number of objects to split a leaf
 * @param n_splits       Number of random splits during learning
 */
template <typename T, typename U>
ExtraTreesRegressor<T,U>::ExtraTreesRegressor(size_t n_estimators, size_t n_min, 
				              size_t n_splits) {

	// keep pars
	this->n_min = n_min;
	this->n_splits = n_splits;
	
	// Set random seed.
	std::random_device rd;
	std::mt19937 gen(rd());
	this->gen = gen;

	// generate the forest
	for (size_t i = 0; i < n_estimators; ++i) {
		forest.push_back(std::make_shared<ExtraTree<T,U>>(n_min, n_splits, this->gen));
	}
}

/**
 * Constructor 
 * @param n_min     Number of minimal data to split a node
 * @param n_splits  Number of random split during learning
 * @param seed      Random generator seed
 */
template <typename T, typename U>
ExtraTreesRegressor<T,U>::ExtraTreesRegressor(size_t n_estimators, size_t n_min, 
		  		              size_t n_splits, long seed) {
	// keep pars
	this->n_min = n_min;
	this->n_splits = n_splits;
	this->seed = seed;
	
	// Set random seed.
	std::random_device rd;
	std::mt19937 gen(seed);
	this->gen = gen;

	// generate the forest
	for (size_t i = 0; i < n_estimators; ++i) {
		forest.push_back(std::make_shared<ExtraTree<T,U>>(n_min, n_splits, this->gen));
	}
}
	
/**
* Make the trees from a set of Features/Label pairs
*
* @param X    Input training samples
* @param y    Output training samples
*/
template <typename T, typename U>
void ExtraTreesRegressor<T,U>::fit(const std::vector<std::vector<T>>& X, 
				   const std::vector<U>& y){
	for (std::shared_ptr<ExtraTree<T,U>> tree: this->forest){
		tree->fit(X, y);
	}
}

/**
* Make the trees from a set of Features/Label pairs with bootstrapping the
* training set between each tree.
*
* @param X    Input training samples
* @param y    Output training samples
*/
template <typename T, typename U>
void ExtraTreesRegressor<T,U>::fit_bootstrap(const std::vector<std::vector<T>>& X,
			   	             const std::vector<U>& y){
	
    // Randomly draw a cut.
	size_t N = y.size();
        std::uniform_int_distribution<size_t> rand_length(0, N - 1);
	auto gen = this->gen;

	for (std::shared_ptr<ExtraTree<T,U>> tree: this->forest){
		std::vector<std::vector<T>> Xs(N);
		std::vector<U> ys(N);
		for (size_t i = 0; i < N; ++i) {
				size_t idx = rand_length(gen);
				Xs[i] = X[idx];
				ys[i] = y[idx];
		}
		tree->fit(Xs, ys);
	}
}
	
/** 
 * Predict label from Extra Tree given a test vector of attributes.
 * 
 * @param X   Input test data
 * @return    predictions (one per tree)
 */
template <typename T, typename U>
std::vector<U> ExtraTreesRegressor<T,U>::predict(const std::vector<T>& X){
	std::vector<U> ypred(this->forest.size());
	for (size_t i = 0; i < ypred.size(); ++i) {
		ypred[i] = this->forest[i]->predict(X);
	}
	return ypred;
}

}
