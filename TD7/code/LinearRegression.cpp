#include "LinearRegression.hpp"
#include "Dataset.hpp"
#include "Regression.hpp"
#include<iostream>
#include<cassert>


LinearRegression::LinearRegression( Dataset* dataset, int col_regr ) 
: Regression(dataset, col_regr) {
	SetCoefficients();
}

LinearRegression::~LinearRegression() {
	m_beta->resize(0);
	delete m_beta;
}

void LinearRegression::SetCoefficients() {
	m_beta = new Eigen::VectorXd;
	int n = m_dataset->GetNbrSamples();
	int d = m_dataset->GetDim();
	Eigen::MatrixXd X(n, d);
	for (int i = 0; i < n; i++) {
		X(i, 0) = 1;
	}
	for (int i = 0; i < n; i++) {
		int z = 1;
		for (int j = 0; j < d; j++) {
			if (j == m_col_regr)continue;
			X(i, z) = m_dataset->GetInstance(i)[j];
			z++;
		}
	}
	Eigen::VectorXd y(n);
	for (int i = 0; i < n; i++) {
		y[i] = m_dataset->GetInstance(i)[m_col_regr];
	}
	Eigen::MatrixXd A = X.transpose() * X;
	Eigen::VectorXd rightPart = X.transpose() * y;
	*(m_beta) = A.colPivHouseholderQr().solve(rightPart);

}

const Eigen::VectorXd* LinearRegression::GetCoefficients() const {
	if (!m_beta) {
		std::cout<<"Coefficients have not been allocated."<<std::endl;
		return NULL;
	}
	return m_beta;
}

void LinearRegression::ShowCoefficients() const {
	if (!m_beta) {
		std::cout<<"Coefficients have not been allocated."<<std::endl;
		return;
	}
	
	if (m_beta->size() != m_dataset->GetDim()) {  // ( beta_0 beta_1 ... beta_{d} )
		std::cout<< "Warning, unexpected size of coefficients vector: " << m_beta->size() << std::endl;
	}
	
	std::cout<< "beta = (";
	for (int i=0; i<m_beta->size(); i++) {
		std::cout<< " " << (*m_beta)[i];
	}
	std::cout << " )"<<std::endl;
}

void LinearRegression::SumsOfSquares( Dataset* dataset, double& ess, double& rss, double& tss ) const {
	assert(dataset->GetDim()==m_dataset->GetDim());
	int n = dataset->GetNbrSamples();
	int d = dataset->GetDim() ;
	double ybar = 0;

	Eigen::MatrixXd X(n, d-1);
	for (int i = 0; i < n; i++) {
		int z = 0;
		for (int j = 0; j < d ; j++) {
			if (j == m_col_regr)continue;
			X(i, z) = dataset->GetInstance(i)[z];
			z++;
		}
	}
	Eigen::VectorXd y(n);
	for (int i = 0; i < n; i++) {
		y[i] = dataset->GetInstance(i)[m_col_regr];
		ybar += y[i];
	}
	ybar /= n;
	Eigen::VectorXd y_estimation(n);
	for (int i = 0; i < n; i++) {
		y_estimation[i] = Estimate(X.row(i));
	}
	tss = 0;
	rss = 0;
	ess = 0;
	for (int i = 0; i < n; i++) {
		tss += (y[i] - ybar)*(y[i]-ybar);
		ess += (y_estimation[i] - ybar)*(y_estimation[i]-ybar);
		rss += (y_estimation[i] - y[i]) * (y_estimation[i] - y[i]);
	}
}

double LinearRegression::Estimate( const Eigen::VectorXd & x ) const {
	double is = (*m_beta)[0];
	for (int i = 0; i < x.size(); i++) {
		is += ((*m_beta)[i + 1]) * x[i];
	}
	return is;
}
