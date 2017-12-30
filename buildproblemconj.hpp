//
//  buildproblemconj.hpp
//  FEM
//
//  Created by Aidan Hamilton on 11/14/17.
//  Copyright © 2017 Aidan Hamilton. All rights reserved.
//

#ifndef buildproblemconj_hpp
#define buildproblemconj_hpp

#include <stdio.h>
#include <tuple>
#include <iostream>
using namespace std;
#include <vector>
using namespace std;
#include "Eigen/Dense"
using namespace Eigen;
#include "Eigen/Sparse"
using namespace Eigen;
#include <Eigen/LU>
#include <Eigen/SparseCore>
using namespace Eigen;
#include <fstream>

std::tuple<SparseMatrix<double,RowMajor>,VectorXd> buildproblemconj(vector<vector<int>>& nodes, vector<vector<double>>& nodepoints, vector<int>& boundarynodes, double (*fun)(double,double), double (*bc)(double,double));

#endif /* buildproblemconj_hpp */
