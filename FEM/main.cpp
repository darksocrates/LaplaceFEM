//
//  main.cpp
//  FEM
//
//  Created by Aidan Hamilton on 10/4/17.
//  Copyright © 2017 Aidan Hamilton. All rights reserved.
//

#include <iostream>
using namespace std;
#include <vector>
using namespace std;
#include "StiffnessMatrix.hpp"
#include <Eigen/Dense>
#include <vector>



//function to read mesh data. Returns pointer to array.
template<class T>
vector<vector<T>> readMatrix(string path){
    
    // The result of the read is placed in here
    // In C++, we use a vector like an array but vectors can dynamically grow
    // as required when we get more data.
    std::vector<std::vector<T>>  Data;
    
    // Replace 'Plop' with your file name.
    std::ifstream  file(path);
    
    if(!file.is_open()){
        cout << "file did not open";
    }
    else{
    
    std::string line;
    // Read one line at a time into the variable line:
    while(std::getline(file, line))
    {
        std::vector<T>   lineData;
        std::stringstream  lineStream(line);
        
        T value;
        // Read an integer at a time from the line
        while(lineStream >> value)
        {
            // Add the integers from a line to a 1D array (vector)
            lineData.push_back(value);
        }
        // When all the integers have been read, add the 1D array
        // into a 2D array (as one line in the 2D array)
        Data.push_back(lineData);
    }
    }
    return Data;
        
    
        
}


int main(int argc, const char * argv[]) {
    // defining paths to load Mesh data
    
    string basepath  = "/Users/AidanHamilton/Desktop/Math_Research/AdvNumericalAnalysis/MeshData/";
    
    string lvl = "2d_laplace_lvl1/";
    string nodepath  = basepath + lvl + "element_node";
    
    string cordpath = basepath + lvl + "cordinati";
    
    
    // load node data as vector vector 'array'
    vector<vector<int>> nodes = readMatrix<int>(nodepath);
    // load node cordinates as vector vector 'array'
    vector<vector<double>> cord = readMatrix<double>(cordpath);
    
    // test the stiffness matrix function
    
    MatrixXd blah  = StiffnessMatrix(nodes,cord);
    
    
    return 0;
}