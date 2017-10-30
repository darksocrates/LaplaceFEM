//
//  StiffnessMatrix.cpp
//  FEM
//
//  Created by Aidan Hamilton on 10/18/17.
//  Copyright © 2017 Aidan Hamilton. All rights reserved.
//

#include "StiffnessMatrix.hpp"

// following build then assemble method for generating stiffness matrix. After developing this matrix its converted into sparse format to then be
// used in sparse linear solver operations

/*
 INPUTS 
 
 nodes - nx3 matrix of element to node map -- passed by reference
 nodepoints - m x 2 matrix of node to x,y cordinate map -- passed by reference
 
 
 FUTURE [need to create] [going to create using sparse matrix algebra using the eigen library (THEY HAVE TRIPLETS TO SPARSE CONVERSION!!!!)
 -------------------------------------------
 edges - cx3 matrix of element to edges map
 boundary - matrix containing boundary edges
 -------------------------------------------
 
 OUTPUT
 
 A - The assembled stiffness matrix for the FEM for poissons with linear interpolation.
 
 */


MatrixXd StiffnessMatrix(vector<vector<int>>& nodes, vector<vector<double>>& nodepoints)
{
    //create necessary constants in more readable form
    unsigned int numberelements = nodes.size();
    
    //create some other constant vectors
    
    const int8_t enumtriangle[3][2]= {{2,3},{3,1},{1,2}};
    
    //preallocation
    double localstiffness[numberelements-1][3][3];
    double trianglepoints[3][2];
    //double a[3];
    double b[3];
    double c[3];
    
    //creation of local stiffness matrices
    for (int globalindex = 0; globalindex < numberelements; globalindex++){
        // calculation of constant vectors neccessary below. Using simplified method particular to Poisson's problem. Reference is An
        // introduction to the finite element method with applications to non-linear problems by R.E.
        // White, John Wiley & Sons.
        
        
        //create arrary contining x,y cordinates of the nodes of the current element.
            for(int i = 0; i<=2;i++){
                trianglepoints[i][0] = nodepoints[nodes[globalindex][i]][0];
                trianglepoints[i][1] = nodepoints[nodes[globalindex][i]][1];
                }
        
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
        
        
        // local stiffness matrix calculation
        // abuses fact matrix is symmetric so only compute upper diagonal portion, then copy for bottom portion
            for(int i = 0; i<=2;i++){
                //populate upper triangle portion
                for(int j = i; j<=2;j++){
                    
                    localstiffness[globalindex][i][j] = (b[i]*b[j]+c[i]*c[j])/(4*trianglearea);
                
                }
                //now copy to fill bottom 3 empty entries
            
                localstiffness[globalindex][1][0] = localstiffness[globalindex][0][1];
                
                localstiffness[globalindex][2][0] = localstiffness[globalindex][0][2];
                
                localstiffness[globalindex][2][1] = localstiffness[globalindex][1][2];
                
            }
        
        
    }
    
    
    
    //write local stiffness to file for debuggging
    
    std::ofstream output;
    output.open("localstiffnesstset.txt");
    
    for(int i = 0; i<numberelements;i++){
        for(int index1 = 0; index1 <=2; index1++){
            for(int index2 = 0; index2 <=2;index2++){
                
                output << localstiffness[i][index1][index2] << " ";
                
                
            }
            output << "\n";
            
        }
        output << "\n";
        output << "\n";
    }
    output.close();
    // crap because I hate warning errors while working...
    MatrixXd A;
    return A;
    
    
    
}