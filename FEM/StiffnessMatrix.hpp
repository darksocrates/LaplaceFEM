//
//  StiffnessMatrix.hpp
//  FEM
//
//  Created by Aidan Hamilton on 10/18/17.
//  Copyright © 2017 Aidan Hamilton. All rights reserved.
//

#ifndef StiffnessMatrix_hpp
#define StiffnessMatrix_hpp

#include <stdio.h>
#include <iostream>
using namespace std;
#include "Eigen/Dense"
using namespace Eigen;
#include "Eigen/Sparse"
using namespace Eigen;
#include <Eigen/LU>
using namespace Eigen;
#include <fstream>
#include <vector> 


MatrixXd StiffnessMatrix(vector<vector<int>>& nodes, vector<vector<double>>& nodepoints)

#endif /* StiffnessMatrix_hpp */
