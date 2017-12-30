//
//  main.cpp
//  FEM
//
//  Created by Aidan Hamilton on 10/4/17.
//  Copyright © 2017 Aidan Hamilton. All rights reserved.
//

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include "buildproblemconj.hpp"
#include <Eigen/SparseCore>

using namespace std;

const double pi = 3.141592653589793;


//function to read mesh data. Returns array.
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

//function to read mesh data. Returns array.
template<class T>
vector<vector<T>> readMatrix2(string path){
    
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
                lineData.push_back(value-1);
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
    
    //creates the element to nodes matrix in sparse format.
    
    // gets the number of elements in the mesh
    int len = nodes.size();
    
    //Eigen::Triplets is a template of the std::vector function from the standard library that stores three pieces of information at every index of the vector, i,j and value. Where it is assumed that i,j are the i and j positions of the value in some given matrix. This is used for filling am Eigen Sparse matrix.
    
    typedef Eigen::Triplet<int> T;
    std::vector<T> enholder;
    //reserve memory for the templated vector enholder
    enholder.reserve(size*3);
    
    //converts the given element to node matrix in vector<vecotr<int>> format to this Eigen::Triplet format.
    for (int i = 0; i<len; ++i){
        for (int j = 0; j<3; ++j){
            enholder.push_back(T(i,nodes[i][j],1));
        }
    }
    
    //declare the elementnode Sparse matrix.
    
    Eigen::SparseMatrix<int,RowMajor> elementnode(len,size);
    
    //fill that sparse matrix with the Triplet elementnode data.
    
    elementnode.setFromTriplets(enholder.begin(),enholder.end());
    //remove excess memory.
    elementnode.prune(1);
  
    
    //we now create the node to element matrix. Transpose is an optmized function for generating sparse matrix transpose of a spare matrix
    
    SparseMatrix<int,RowMajor> nodeelement = elementnode.transpose();
    
    // node to node matrix. The * operator is overloaded to function a linear algebra operator in Eigen and will output exactly what you think it should output. Eigen has dense matrixs and vectors so if you use * between two compatible dense matrices it will output a dense matrix and will do the dense * operation. In this case as both objects are sparse matrices the * operator outputs a sparse matrix and uses the sparse matrix muliplication alogrithm, which I believe is optimized.
    
    SparseMatrix<int,RowMajor> nodenode = (nodeelement*elementnode);
    
    //now iterate through and find all node-node pairs with value 1. This means that the edge node1 to node2 belongs to a single element. So that edge is on the boundary. which means that node1 is a boundary node. Now this also means that node2 is on the  boundary as well but because I'm taking a lazy iteration procedure I don't add that to the list yet, and wait for that node2 to be the node1 in the iteration.
    
    std::vector<int> boundarynodes;
    
    int test = 0;
    
    //outersize is the number of rows in the sparse matrix.
    for (int k=0; k<nodenode.outerSize(); ++k){
        //the InnerIterator is simply Eigen's method of accessing the elements of a Sparse matrices given row.
        for (SparseMatrix<int,RowMajor>::InnerIterator it(nodenode,k); it; ++it)
        {
            if (it.value() == 1 ){
                test = 1;
            }
        }
        if (test == 1){
            boundarynodes.push_back(k);
            test = 0;
        }
    }
  // this std::make_tuple thing is just a creating a tuple object that contains references to the data you actually want to pass through. It's jsut the way C++ handles returning multiple outputs.
 
    return std::make_tuple(elementnode,nodeelement,nodenode,boundarynodes);
}

//function describing boundary data of the Poisson Equation
double bc(double x, double y){
    //change f to the function you want
    double f = sin(x)*sin(y);
    return f;
}

//Function describing the RHS of the Poisson Equation
double f(double x, double y){
    //change c to the function you want
    double c = 2*sin(x)*sin(y);
    return c;
}

//function describing the true solution if known. Used for test problems.
double thetruth(double x, double y){
    //change f to the function you want.
    double f = sin(x)*sin(y);
    return f;
}


int main(int argc, const char * argv[]) {
    
    //these four lines are for generating the path to the mesh files for the node cordinate data and the element to node matrix. It is split up to take advantage of the specific file strucutre I was working with. Change as needed.
    string basepath  = "/Users/AidanHamilton/desktop/Math_Research/AdvNumericalAnalysis/FEM_Laplace/MeshData/";

    string lvl = "2d_anis_lvl2/";
    
    string nodepath  = basepath + lvl + "element_node";
    
    string pointspath = basepath + lvl + "coordinati";
    
    //load data, using the ReadMatrix function.
    
    vector<vector<int>> nodes = readMatrix2<int>(nodepath);
    
    vector<vector<double>> cord = readMatrix<double>(pointspath);
    
    // create mesh relations using the create relations function. I only actually use the boundarynodes data, the rest was used for debugging purposes. The other stuff is still kept because why not.
    SparseMatrix<int,RowMajor> a,b,c;
    std::vector<int> boundarynodes;
    int size = cord.size();

    
    std::tie(a,b,c,boundarynodes) = createrelations(nodes,size);
   
    
    //the following commmented out code was used for debugging purposes and wrote the boundarynode data to a text file for later viewing. I keep it in but commented out because it can be useful.
    
    /*
    FILE * pfile;
    
    pfile = fopen("/Users/AidanHamilton/Desktop/testingFEM/boundarynodes.txt","w");

    for(int i =0; i<boundarynodes.size();++i){
        fprintf(pfile,"%5d \n",boundarynodes[i]);
    }
    */
    
    //declare the stiffness matrix A and the RHS vecotr RHS.
    SparseMatrix<double,RowMajor> A;
    VectorXd RHS;

    // generate the stiffness matrix A and the RHS vector RHS for the Poisson problem. Note that &f and &bc are function pointers.
    std::tie(A,RHS) = buildproblemconj(nodes, cord, boundarynodes, &f, &bc);

    //Eigen has its own implementation of a conjugate gradient solver for sparse matrices. It is a class so below I create an object of it.
    ConjugateGradient<SparseMatrix<double,RowMajor>,Lower|Upper> solver;
    
    //.compute(A) is neccessary to call before calling solve, it does some preconditioning an generates a sparsity pattern used internally in the ConjugateGradient object.
    solver.compute(A);
    //set the tolerance for the conjugategradient solution, computed in the euclidean norm.
    double tol = 1.0e-10;
    solver.setTolerance(tol);
    // solve the Au = RHS using conjugate gradient.
    VectorXd u  = solver.solve(RHS);
    //output the error and number of iterations in the conjugate gradient.
    cout<<solver.error()<<" "<<solver.iterations();
    
    //calculate the true solution using the known truesolution
    VectorXd truesol(cord.size());
    for(int i =0; i<cord.size(); ++i){
        truesol(i) = thetruth(cord[i][0],cord[i][1]);
    }
    //calculate the absolute difference between this true solution and the computed solution.
    auto error = (u-truesol).cwiseAbs();
    
    
    // the following outputs the true solution, computed solution and computed error to a text file with path in the fopen function.
    FILE * afile;
    
    afile = fopen("/Users/AidanHamilton/Desktop/testingFEM/debug.txt","w");
    
    for (int i =0; i<u.size(); ++i){
        
        // this is the c method for writing to a file, what this means is that I output everything
        fprintf(afile,"%.7g %.7g %.7g \n", u(i),truesol(i),error(i));
        
    }
    
    return 0;
}