#include "Regression.hpp"
#include "Dataset.hpp"

Regression::Regression(Dataset* dataset, int col_regr) {
    m_dataset = dataset;
    m_col_regr = col_regr;
}

int Regression::GetColRegr() const {
    return m_col_regr;
}

Dataset* Regression::GetDataset() const {
    return m_dataset;
}
