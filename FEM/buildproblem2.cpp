//
//  buildproblem.cpp
//  FEM
//
//  Created by Aidan Hamilton on 10/30/17.
//  Copyright © 2017 Aidan Hamilton. All rights reserved.
//

#include "buildproblem2.hpp"


/*
 NOTES
 
 This function utilizes the EIGEN C++ library to handle the sparse matrices.
 
 
 INPUTS
 
 nodes - nx3 vector<vector<int>> of element to node map -- passed by reference
 nodepoints - m x 2 vector<vector<double>> of node to x,y cordinate map -- passed by reference
 boundary nodes -  vecotr<int> of the nodes on the boundary -- passed by reference
 fun - a function corresponding to the dirichlet boundary data -- passed by function pointer
 
 
 OUTPUT
 
 A - The assembled stiffness matrix for the FEM for poissons with linear interpolation in CSR format. Contains only the the interior nodes.
 
 RHS - The assembled RHS for the FEM with enforced boundary conditions, contains only the interior nodes. The boundary nodes have been removed and had their boundary conditions enforced. 
 
 BoundaryValue - VectorXd containing the weights forthe boundary nodes.
 */




std::tuple<SparseMatrix<double,RowMajor>,VectorXd> buildproblem(vector<vector<int>>& nodes, vector<vector<double>>& nodepoints, vector<int>& boundarynodes, double (*fun)(double,double), double (*bc)(double,double))
{
    //create necessary constants in more readable form
    const unsigned int numberelements = nodes.size();
    const unsigned int numnodes = nodepoints.size();
    
    //create some other constant vectors
    
    //local stiffness matrix in triplet format, essentially assembles itself this way
    typedef Eigen::Triplet<double> T;
    std::vector<T> localstiff;
    localstiff.reserve(numberelements*9);
    
    
    
    const int enumtriangle[3][2]= {{1,2},{2,0},{0,1}};
    
    //preallocation
    
    VectorXd RHS(numnodes);
    
    double trianglepoints[3][2];
    //double a[3];
    double b[3];
    double c[3];
    
    //creation of local stiffness matrices
    for (int globalindex = 0; globalindex < numberelements; globalindex++){
        // calculation of constant vectors neccessary below. Using simplified method particular to Poisson's problem. Reference is An
        // introduction to the finite element method with applications to non-linear problems by R.E.
        // White, John Wiley & Sons.
        
        
        //create array contining x,y cordinates of the nodes of the current element.
        for(int i = 0; i<=2;i++){
            trianglepoints[i][0] = nodepoints[nodes[globalindex][i]][0];
            trianglepoints[i][1] = nodepoints[nodes[globalindex][i]][1];
        }
        
        //create array containing the xy cordinates of the midpoints of the nodes of the current element
        
        double midp[3][2];
        
        midp[2][0] = (trianglepoints[0][0] + trianglepoints[1][0])/2; midp[2][1] = (trianglepoints[0][1] + trianglepoints[1][1])/2;
        midp[0][0] = (trianglepoints[1][0] + trianglepoints[2][0])/2; midp[0][1] = (trianglepoints[1][1] + trianglepoints[2][1])/2;
        midp[1][0] = (trianglepoints[2][0] + trianglepoints[0][0])/2; midp[1][1] = (trianglepoints[2][1] + trianglepoints[0][1])/2;
        
        //populate a,b,c vectors
        //NOTE: for now not forming a as not needed in stiffness matrix formation
        for(int i = 0; i<=2;i++){
            int jj = enumtriangle[i][0];
            int mm = enumtriangle[i][1];
            
            //a[i] = trianglepoints[jj][0]*trianglepoints[mm][1] - trianglepoints[mm][0]*trianglepoints[jj][1];
            b[i] = trianglepoints[jj][1] - trianglepoints[mm][1];
            c[i] =trianglepoints[mm][0] - trianglepoints[jj][0];
            
            
        }
        
        //calculate area of elements triangle.
        Matrix3d triholder;
        triholder << 1,trianglepoints[0][0],trianglepoints[0][1],1,trianglepoints[1][0],trianglepoints[1][1],1,trianglepoints[2][0],trianglepoints[2][1];
        double trianglearea =0.5*std::abs(triholder.determinant());
        
        // local stiffness matrix calculation and RHS calculation
        
        for(int i = 0; i<=2;++i){
            int ii = nodes[globalindex][i];
            
            //stiffness
            for(int j = 0; j<=2;++j){
                
                int jj = nodes[globalindex][j];
                
                localstiff.push_back(T(ii,jj, (b[i]*b[j]+c[i]*c[j])/(4*trianglearea)));
                
            }
            //RHS
            
            //other quadrature, both work.
            
            //RHS(ii) = RHS(ii)+ trianglearea*(fun(trianglepoints[i][0],trianglepoints[i][1])/6+fun(trianglepoints[enumtriangle[i][1]][0],trianglepoints[enumtriangle[i][1]][1])/12+fun(trianglepoints[enumtriangle[i][0]][0],trianglepoints[enumtriangle[i][0]][1])/12);
            
            //midpoint quadrature
            
            double f1 = fun(midp[enumtriangle[i][0]][0],midp[enumtriangle[i][0]][1]);
            double f2 = fun(midp[enumtriangle[i][1]][0],midp[enumtriangle[i][1]][1]);
            
            RHS(ii) = RHS(ii) + trianglearea/6*(f1+f2);
            
        }
        
    }
    
    //assemble matrix in EIGEN Sparse format
    
    Eigen::SparseMatrix<double,RowMajor> A(numnodes,numnodes);
    A.setFromTriplets(localstiff.begin(), localstiff.end());
    
    //now iterate through sparse matrix and enforce direchlet boundary conditions. I remove the row containing the boundary node in A and RHS and then add the boundary value to a double variable which acucumalates all the boundary values. I then at the en of iterating through the matrix subtract this value from all the remaining RHS entries.
    
    int bnxt= 0;
    
    for (int k=0; k<A.outerSize(); ++k){
        if (k == boundarynodes[bnxt]) {
            bnxt++;
            // for all nodes on boundary set only the entry corresponding to that nodes basis to 1 in A, all others 0.
            for (SparseMatrix<double,RowMajor>::InnerIterator it(A,k); it; ++it)
            {
                if (it.col() == k){
                    it.valueRef() = 1;
                }
                else{
                    it.valueRef() = 0;
                }
                
            }
            // now change the RHS
            
            RHS(k) = bc(nodepoints[k][0],nodepoints[k][1]);
        }
    }
    A.prune(.1);
    
    //now I iterate through the non-zero entries again and zero out the columns of the boundary nodes and change the RHS's appropriately

    bnxt = 0;
    int iter = 0;
    for (int k=0; k<A.outerSize(); ++k){
        
        //avoid working on boundary rows.
        if (k == boundarynodes[bnxt]) {
            bnxt++;
        }
        //non-boundary rows
        else {
        iter = 0;
        
        for (SparseMatrix<double,RowMajor>::InnerIterator it(A,k); it; ++it){
            //check if any column in the row is the same index as a boundary rows, if so zero out. Takes advantaged of how everything is ordered from smallest to largest (both it.col() values and boundary nodes indexs).
            
            while(boundarynodes[iter]<it.col()){
                iter++;
            }
            if (boundarynodes[iter] == it.col()){
                
                RHS(k)  = RHS(k) - it.value()*RHS(boundarynodes[iter]);
                it.valueRef() = 0;
                
            }
            
        }
            
        }
    }
    
    
    return std::make_tuple(A,RHS);
    
}