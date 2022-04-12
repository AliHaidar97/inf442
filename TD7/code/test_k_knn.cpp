
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <ANN/ANN.h>
#include "Dataset.hpp"
#include "KnnRegression.hpp"

/** @file
 * Test suite for the impact of K in KNN regression.
 *  This executable will put the provided CSV files in objects of class Dataset 
 *  and perform KNN regression for k=1 to 40
 *  for the samples in the test file using instances
 *  of the train file. It then prints the resulting MSE for each value of k.
*/
int main(int argc, char* argv[]) {
	if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <train_file> <test_file> [ <column_for_regression> ]" << std::endl;
        return 1;
    }

	Dataset train_dataset(argv[1]);
	Dataset regr_dataset(argv[2]);
    
    int col_regr;
    if (argc == 4) {
    	col_regr = atoi(argv[3]);
    } else {
		col_regr = train_dataset.GetDim()-1;
		std::cout<< "No column specified for regression, assuming last column of dataset ("<< col_regr <<")."<<std::endl;    
    }

	train_dataset.Show(false);

	assert(train_dataset.GetDim() == regr_dataset.GetDim());
    // otherwise doesn't make sense

	int dim = train_dataset.GetDim();
	
	std::ofstream out("output.txt");
	for (int k = 1; k <= 100; k++) {
		KnnRegression knn_reg(k, &train_dataset, col_regr);

        double mse = 0;
        for (int i = 0; i < regr_dataset.GetNbrSamples(); i++) {
            std::vector<double> sample = regr_dataset.GetInstance(i);
            // extract column for regression
            Eigen::VectorXd query(regr_dataset.GetDim() - 1);
            for (int j = 0, j2 = 0; j < train_dataset.GetDim() - 1 && j2 < train_dataset.GetDim(); j++, j2++) {
                if (j == col_regr && j2 == col_regr) {
                    j--;
                    continue;
                }
                query[j] = sample[j2];
            }
            double estim = knn_reg.Estimate(query);
            mse += (estim - sample[col_regr]) * (estim - sample[col_regr]) / regr_dataset.GetNbrSamples();
        }
        out << k << ' ' << mse << '\n';
	}
    out.close();
	annClose();
	return 0;
}
