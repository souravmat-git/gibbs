//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FourthRankTensorVoigtNotation36.h"
#include "DenseMatrix.h"
#include <fstream>
using namespace std;

registerMooseObject("gibbsApp", FourthRankTensorVoigtNotation36);

template<>
InputParameters
validParams<FourthRankTensorVoigtNotation36>()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Converts a fourth rank tensor into a 6 by 6 array.");
  params.addRequiredParam<MaterialPropertyName>("C_tensor_name", "Fourth rank tensor");
  params.addRequiredParam<MaterialPropertyName>("C_matrix_name", "Stiffness_matrix name");
  return params;
}

FourthRankTensorVoigtNotation36::FourthRankTensorVoigtNotation36(const InputParameters & parameters)
 : Material(parameters),
  _C(getMaterialProperty<RankFourTensor>
               (getParam<MaterialPropertyName>("C_tensor_name"))),
  _C_matrix_name(getParam<MaterialPropertyName>("C_matrix_name")),
  _V(declareProperty<DenseMatrix<Real>>(_C_matrix_name)) 
{
}

void
FourthRankTensorVoigtNotation36::computeQpProperties()
{
  
   _V[_qp].resize(6,6);
   _V[_qp].zero();
    
  //Diagonal terms: 11, 22, 33, 44, 55, 66
  _V[_qp](0,0)  = _C[_qp](0,0,0,0);
  _V[_qp](1,1)  = _C[_qp](1,1,1,1);
  _V[_qp](2,2)  = _C[_qp](2,2,2,2);
  _V[_qp](3,3)  = _C[_qp](1,2,1,2);
  _V[_qp](4,4)  = _C[_qp](0,2,0,2);
  _V[_qp](5,5)  = _C[_qp](0,1,0,1);

  // First row begins with index 0 5 elements
  _V[_qp](0,1) = _C[_qp](0,0,1,1);
  _V[_qp](0,2) = _C[_qp](0,0,2,2);
  _V[_qp](0,3) = _C[_qp](0,0,1,2);
  _V[_qp](0,4) = _C[_qp](0,0,0,2);
  _V[_qp](0,5) = _C[_qp](0,0,0,1);
  
  // No major symmetry C_ijkl != C_klij 
  _V[_qp](1,0) = _C[_qp](1,1,0,0);
  _V[_qp](2,0) = _C[_qp](2,2,0,0);
  _V[_qp](3,0) = _C[_qp](1,2,0,0);
  _V[_qp](4,0) = _C[_qp](0,2,0,0);
  _V[_qp](5,0) = _C[_qp](0,1,0,0);
  
  // Second row begins with index 1 and has 4 elements
  _V[_qp](1,2) = _C[_qp](1,1,2,2);
  _V[_qp](1,3) = _C[_qp](1,1,1,2);
  _V[_qp](1,4) = _C[_qp](1,1,0,2);
  _V[_qp](1,5) = _C[_qp](1,1,0,1); 
  
  // No major symmetry _C_ijkl != C_klij 
  _V[_qp](2,1) = _C[_qp](2,2,1,1);
  _V[_qp](3,1) = _C[_qp](1,2,1,1);
  _V[_qp](4,1) = _C[_qp](0,2,1,1);
  _V[_qp](5,1) = _C[_qp](0,1,1,1);
 
  // Third row 3 elements
  _V[_qp](2,3) = _C[_qp](2,2,1,2);
  _V[_qp](2,4) = _C[_qp](2,2,0,2);
  _V[_qp](2,5) = _C[_qp](2,2,0,1);
  
  //No major symmetry _C_ijkl != C_klij 
  _V[_qp](3,2) = _C[_qp](1,2,2,2);
  _V[_qp](4,2) = _C[_qp](0,2,2,2);
  _V[_qp](5,2) = _C[_qp](0,1,2,2);
 
  // Fourth row 2 elements
  _V[_qp](3,4) = _C[_qp](1,2,0,2);
  _V[_qp](3,5) = _C[_qp](1,2,0,1); 
  
  //No major symmetry _C_ijkl = C_klij 
  _V[_qp](4,3) = _C[_qp](0,2,1,2);
  _V[_qp](5,3) = _C[_qp](0,1,1,2);
  
   // Fifth row 1 elements
   _V[_qp](4,5) = _C[_qp](0,2,0,1);
  
  // Due to major symmetry _C = C_klij 
  _V[_qp](5,4) = _C[_qp](0,1,0,2);
  
  
   //Generate a CSV file to store the array
   if (_t_step)
   {
     ofstream C_matrix;
     C_matrix.open(_C_matrix_name + ".csv");
  
     //store the array V is a m by n matrix
     for (unsigned int i = 0; i< 6 ; i++){
        for (unsigned int j = 0; j< 6; j++){
          C_matrix << setprecision(5)  << scientific << _V[_qp](i,j) << ", ";
        }
     C_matrix << "\n"; 
    }
   } 
}
