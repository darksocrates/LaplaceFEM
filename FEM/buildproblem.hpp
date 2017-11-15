//
//  buildproblem.hpp
//  FEM
//
//  Created by Aidan Hamilton on 10/30/17.
//  Copyright © 2017 Aidan Hamilton. All rights reserved.
//

#ifndef buildproblem_hpp
#define buildproblem_hpp

#include <stdio.h>
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



std::tuple<SparseMatrix<double,RowMajor>,VectorXd> buildproblem(vector<vector<int>>& nodes, vector<vector<double>>& nodepoints, vector<int>& boundarynodes, double (*fun)(double,double), double (*bc)(double,double)) ;

#endif /* buildproblem_hpp */