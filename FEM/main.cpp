//
//  main.cpp
//  FEM
//
//  Created by Aidan Hamilton on 10/4/17.
//  Copyright © 2017 Aidan Hamilton. All rights reserved.
//

#include <iostream>
using namespace std;
#include <Eigen/Dense>
#include <vector>
using namespace std;
#include <fstream>
#include "buildproblem.hpp"
#include "buildproblem2.hpp"
#include <Eigen/SparseCore>
#include "StiffnessMatrix.hpp"


const double pi = 3.141592653589793;


//function to read mesh data. Returns pointer to array.
template<class T>
vector<vector<T>> readMatrix(string path){
    
    // The result of the read is placed in here
    // In C++, we use a vector like an array but vectors can dynamically grow
    // as required when we get more data.
    std::vector<std::vector<T> >  Data;
    
    std::ifstream file(path);
    
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


//creates the sparse relation tables for the given mesh, in particular gives (indirectly) the boundary
//nodes

std::tuple<Eigen::SparseMatrix<int,RowMajor>,Eigen::SparseMatrix<int,RowMajor>,Eigen::SparseMatrix<int,RowMajor>,std::vector<int>> createrelations(vector<vector<int>>& nodes,int size){
    
    //create element to nodes matrix in sparse format.
    
    int len = nodes.size();
    
    typedef Eigen::Triplet<int> T;
    std::vector<T> enholder;
    enholder.reserve(size*3);
    
    for (int i = 0; i<len; ++i){
        for (int j = 0; j<3; ++j){
            enholder.push_back(T(i,nodes[i][j],1));
        }
    }
    
    Eigen::SparseMatrix<int,RowMajor> elementnode(len,size);
    
    /*
    for(int i =0; i<enholder.size();++i){
        cout<<enholder[i].row()<<" "<<enholder[i].col() << " "<<enholder[i].value()<<"\n";
    }
    */
    elementnode.setFromTriplets(enholder.begin(),enholder.end());
    elementnode.prune(1);
  
    
    //we now create node to element matrix
    
    SparseMatrix<int,RowMajor> nodeelement = elementnode.transpose();
    
    // node to node matrix
    
    SparseMatrix<int,RowMajor> nodenode = (nodeelement*elementnode);
    
    //now iterate through and find all node-node pairs with value 1. That is the edge node1 to node2 belongs to a single element. So that edge is on the boundary. which means that node1 is a boundary node. This method is much improved on the method tried in the below comment.
    
    std::vector<int> boundarynodes;
    
    int test = 0;
    
    for (int k=0; k<nodenode.outerSize(); ++k){
        for (SparseMatrix<int,RowMajor>::InnerIterator it(nodenode,k); it; ++it)
        {
            if (it.value() == 1 ){
                test = 1;
            }
            
            //cout<<it.value()<<"\n";
        }
        if (test == 1){
            boundarynodes.push_back(k);
            test = 0;
        }
    }
  
 
    return std::make_tuple(elementnode,nodeelement,nodenode,boundarynodes);
}

// function describing boundary data.
double bc(double x, double y){
    double f = 1+x+y;
    return f;
}

double f(double x, double y){
    double c = 0;
    return c;
}

double thetruth(double x, double y){
    double f = 1+x+y;
    return f;
}

int main(int argc, const char * argv[]) {
    
    //user input, give me the mesh files!!!!
    string basepath  = "/Users/AidanHamilton/Desktop/Math_Research/AdvNumericalAnalysis/MeshData/";
    
    string lvl = "2d_laplace_lvl1/";
    
    
    
    string nodepath  = basepath + lvl + "element_node";
    
    string pointspath = basepath + lvl + "coordinati";
    
    //load data
    
    vector<vector<int>> nodes = readMatrix<int>(nodepath);
    
    vector<vector<double>> cord = readMatrix<double>(pointspath);


    // create mesh relations
    SparseMatrix<int,RowMajor> a,b,c;
    std::vector<int> boundarynodes;
    int size = cord.size();

    
    std::tie(a,b,c,boundarynodes) = createrelations(nodes,size);
    
    
    FILE * pfile;
    
    pfile = fopen("/Users/AidanHamilton/Desktop/testingFEM/boundarynodes.txt","w");

    for(int i =0; i<boundarynodes.size();++i){
        fprintf(pfile,"%5d \n",boundarynodes[i]);
    }

    //declare function pointers
    
    SparseMatrix<double,RowMajor> A;
    VectorXd RHS;

    
    std::tie(A,RHS) = buildproblem2(nodes, cord, boundarynodes, &f, &bc);

    ConjugateGradient<SparseMatrix<double,RowMajor>,Lower|Upper> solver;
    //BiCGSTAB<SparseMatrix<double,RowMajor>> solver;
    solver.compute(A);
    double tol = 1.0e-10;
    solver.setTolerance(tol);
    VectorXd u  = solver.solve(RHS);
    cout<<solver.error()<<" "<<solver.iterations(); 
    VectorXd truesol(cord.size());
    for(int i =0; i<cord.size(); ++i){
        truesol(i) = thetruth(cord[i][0],cord[i][1]);
    }
    
    //cout<<RHS;
    auto error = (u-truesol).cwiseAbs();
    
    FILE * afile;
    
    afile = fopen("/Users/AidanHamilton/Desktop/testingFEM/debug.txt","w");
    
    for (int i =0; i<u.size(); ++i){
        
        fprintf(afile,"%.7g %.7g %.7g \n", u(i),truesol(i),error(i));
        
    }
    return 0;
}