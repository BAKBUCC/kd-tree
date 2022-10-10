// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Point.hpp"

template <typename value_type>
struct KDNode
{
	value_type val;
	KDNode<value_type>* nodes[2];
	KDNode(const value_type& value) {
	  nodes[0] = nodes[1] = nullptr;
    val = value;
	}
	KDNode(value_type value, KDNode<value_type>* izq, KDNode<value_type>* der)
	{
    nodes[0] = izq;
		nodes[1] = der;
		val = value;
	}
};

template <size_t N, typename ElemType>
class KDTree {
public:
	typedef std::pair<Point<N>, ElemType> value_type;

	KDTree();

	~KDTree();

	KDTree(const KDTree& rhs);
	KDTree& operator=(const KDTree& rhs);

	size_t dimension() const;

	size_t size() const;
	bool empty() const;

	bool find(const Point<N>& pt, KDNode<value_type>**& p) const;

	bool contains(const Point<N>& pt) const;

	void insert(const Point<N>& pt, const ElemType& value);

	ElemType& operator[](const Point<N>& pt);

	ElemType& at(const Point<N>& pt);
	const ElemType& at(const Point<N>& pt) const;




	void knn_helper(Point<N> key, KDNode<value_type>* currentNode, int nivel, std::vector<std::pair<ElemType, double> >& tipo) const;





	ElemType knn_value(const Point<N>& key, size_t k) const;

	std::vector<ElemType> knn_query(const Point<N>& key, size_t k) const;

private:
	size_t dimension_;
	size_t size_;

	KDNode<value_type>* root = nullptr;
};

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
	dimension_ = N;
	size_ = 0;
}

template <typename value_type>
void deleter(KDNode<value_type>* node)
{
	if (node != nullptr) {
		deleter((node->nodes)[1]);
		deleter((node->nodes)[0]);
		delete node;
	}
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
	deleter(root);
}

template <typename value_type>
KDNode<value_type>* CopyAndPutNodes(const KDNode<value_type>* node)
{
	if (node != nullptr) {
		KDNode<value_type>* nodeCopy = new KDNode<value_type>(node->val, CopyAndPutNodes((node->nodes)[0]), CopyAndPutNodes((node->nodes)[1]));
		return nodeCopy;
	}
	return nullptr;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& rhs) {
	root = CopyAndPutNodes(rhs.root);
	dimension_ = rhs.dimension_;
	size_ = rhs.size_;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>& KDTree<N, ElemType>::operator=(const KDTree& rhs) {
	root = CopyAndPutNodes(rhs.root);
	dimension_ = rhs.dimension_;
	size_ = rhs.size_;
	return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
	return N;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
	return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
	return !root;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::find(const Point<N>& pt, KDNode<value_type>**& p) const {
	size_t index = 0;

	//referencia al find del profesor alex cuadros
	for (p = const_cast<KDNode<value_type>**> (&root); *p && ((*p)->val).first != pt; p = &((*p)->nodes[pt[index % N] > (((*p)->val).first)[index % N]])) 
	{
		index++;
	}
	return *p != 0;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N>& pt) const {
	KDNode<value_type>** p;
	if (!find(pt, p)) {
		return false;
	}
	return true;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value) {
	KDNode<value_type>** p;
	if (!find(pt, p)) {
		value_type temp;
		temp.first = pt;
		temp.second = value;
		*p = new KDNode<value_type>(temp);
		size_ += 1;
	}
	((*p)->val).second = value;
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::operator[](const Point<N>& pt) {
	KDNode<value_type>** p;
	if (!find(pt, p)) {
		value_type temp;
		temp.first = pt;
		temp.second = size_;
		*p = new KDNode<value_type>(temp);
		size_ += 1;
		return ((*p)->val).second;
	}
	return ((*p)->val).second;
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) {
	KDNode<value_type>** p;
	if (find(pt, p)) {
		return ((*p)->val).second;
	}
	throw std::out_of_range("");
}

template <size_t N, typename ElemType>
const ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) const {
	KDNode<value_type>** p;
	if (find(pt, p)) {
		return ((*p)->val).second;
	}
	throw std::out_of_range("");
}


template <size_t N, typename ElemType>
void KDTree<N, ElemType>::knn_helper(Point<N> key, KDNode<value_type>* cNode, int nivel, std::vector<std::pair<ElemType, double> >& tipo) const
{
	if (cNode == nullptr)
		return;
	double d = distance((cNode->val).first, key);
	tipo.push_back(std::make_pair((cNode->val).second, d));

	knn_helper(key, (cNode->nodes)[0], ++nivel, tipo);
	knn_helper(key, (cNode->nodes)[1], ++nivel, tipo);
}


template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N>& key, size_t k) const {
	if (k > size_)
		k = size_;

	std::vector<ElemType> vecElem;
	std::vector<ElemType> vals(k);
	vals = knn_query(key, k);
	vecElem.push_back(vals[0]);


	for (size_t i = 1; i < k; i++) {
		bool aux = false;
		for (size_t j = 0; j < vecElem.size(); j++) {
			if (vals[i] == vecElem[j])
				aux = true;
		}
		if (!aux)
		{
			vecElem.push_back(vals[i]);
		}
	}

	std::vector<std::pair<ElemType, size_t> > cont;
	for (size_t i = 0; i < vecElem.size(); i++) 
	{
		cont.push_back(std::make_pair(vecElem[i], 0));
	}


	for (size_t i = 0; i < k; i++) {
		for (size_t j = 0; j < cont.size(); j++) {
			if (vals[i] == cont[j].first)
				cont[j].second += 1;
		}
	}


	std::pair<ElemType, size_t> maximo = cont[0];
	for (size_t i = 1; i < cont.size(); i++) {
		if (maximo.second < cont[i].second)
			maximo = cont[i];
	}
	return maximo.first;
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N>& key,
	size_t k) const {
	std::vector<std::pair<ElemType, double> > tipo;
	KDNode<value_type>* cNode = root;
	size_t nivel = 0;

	knn_helper(key, cNode, nivel, tipo);
	std::sort(tipo.begin(), tipo.end(), [](const auto& x, const auto& y) { return x.second < y.second; });

	std::vector<ElemType> val;
	for (size_t i = 0; i < k; i++) {
		val.push_back(tipo[i].first);
	}
	return val;
}

// TODO(me): finish the implementation of the rest of the KDTree class

//gracias por la ayuda a muchos de mis compañeros que me ayudaron a corregir y completar el trabajo
#endif  // SRC_KDTREE_HPP_




