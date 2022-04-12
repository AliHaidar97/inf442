
#include "KnnRegression.hpp"
#include<iostream>
#include <ANN/ANN.h>


KnnRegression::KnnRegression( int k, Dataset* dataset, int col_regr)
: Regression(dataset, col_regr) {
	m_k = k;
    m_k = k;
    m_dataPts = annAllocPts(m_dataset->GetNbrSamples(), m_dataset->GetDim() - 1);

    for (int i = 0; i < m_dataset->GetNbrSamples(); i++) {
        int z = 0;
        for (int j = 0; j < m_dataset->GetDim(); j++) {
            if (j ==m_col_regr)continue;
            m_dataPts[i][z] = m_dataset->GetInstance(i)[j];
            z++;
        }
    }

    m_kdTree = new ANNkd_tree(m_dataPts, m_dataset->GetNbrSamples(), m_dataset->GetDim() - 1);
}

KnnRegression::~KnnRegression() {
    delete m_kdTree;
    annDeallocPts(m_dataPts);
    annClose();
}

double KnnRegression::Estimate(const Eigen::VectorXd & x) const {
	assert(x.size()==m_dataset->GetDim()-1);
    ANNidxArray nnIdx; // near neighbor indices
    ANNdistArray dists; // near neighbor distances
    nnIdx = new ANNidx[m_k]; // allocate near neigh indices
    dists = new ANNdist[m_k]; // allocate near neighbor dists
    ANNpoint xtemp = annAllocPt(x.size());
    for (int i = 0; i < x.size(); i++) {
        xtemp[i] = x[i];
    }
    m_kdTree->annkSearch(xtemp, m_k, nnIdx, dists);
    double p = 0.0;
    for (int i = 0; i < m_k; i++) {
        p += m_dataset->GetInstance(nnIdx[i])[m_col_regr];
    }
    p /= m_k;

    delete[] nnIdx;
    delete[] dists;
    annDeallocPt(xtemp);
    return p;
}

int KnnRegression::GetK() const {
	return m_k;
}

ANNkd_tree* KnnRegression::GetKdTree() const {
	return m_kdTree;
}
