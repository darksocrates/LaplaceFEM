//
//  StiffnessMatrix.hpp
//  FEM
//
//  Created by Aidan Hamilton on 10/18/17.
//  Copyright © 2017 Aidan Hamilton. All rights reserved.
//


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
using namespace Eigen;
#include <fstream>


#ifndef StiffnessMatrix_hpp
#define StiffnessMatrix_hpp
MatrixXd StiffnessMatrix(vector<vector<int>>& nodes, vector<vector<double>>& nodepoints);

#endif /* StiffnessMatrix_hpp */
