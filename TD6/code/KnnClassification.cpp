
#include "KnnClassification.hpp"
#include <iostream>
#include <ANN/ANN.h>

 
KnnClassification::KnnClassification(int k, Dataset* dataset, int col_class)
: Classification(dataset, col_class) {
    // TODO Exercise 1.1
    m_k = k;
    m_dataPts = annAllocPts(m_dataset->getNbrSamples(), m_dataset->getDim()-1);

    for (int i = 0; i < m_dataset->getNbrSamples(); i++) {
        int z = 0;
        for (int j = 0; j < m_dataset->getDim(); j++) {
            if (j == m_col_class)continue;
            m_dataPts[i][z] = m_dataset->getInstance(i)[j];
            z++;
        }
    }

    m_kdTree = new ANNkd_tree(m_dataPts, m_dataset->getNbrSamples(), m_dataset->getDim()-1);

}

KnnClassification::~KnnClassification() {
    // TODO Exercise 1.1
    delete m_kdTree;
    annDeallocPts(m_dataPts);
    annClose();
}

int KnnClassification::Estimate(const ANNpoint & x, double threshold) {
    ANNidxArray nnIdx; // near neighbor indices
    ANNdistArray dists; // near neighbor distances
    nnIdx = new ANNidx[m_k]; // allocate near neigh indices
    dists = new ANNdist[m_k]; // allocate near neighbor dists
    m_kdTree->annkSearch(x,m_k, nnIdx, dists);
    double p = 0.0;
    for (int i = 0; i < m_k; i++) {
        p += m_dataset->getInstance(nnIdx[i])[m_col_class];
    }
    p /= m_k;

    delete[] nnIdx;
    delete[] dists;

    return (p>threshold);
}

int KnnClassification::getK() {
    return m_k;
}

ANNkd_tree* KnnClassification::getKdTree() {
    return m_kdTree;
}
