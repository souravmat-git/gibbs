//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
class FourthRankTensorVoigtNotation;

//MOOSE includes
#include "Material.h"

template <>
InputParameters validParams<FourthRankTensorVoigtNotation>();

class FourthRankTensorVoigtNotation : public Material
{
public:
  FourthRankTensorVoigtNotation(const InputParameters & parameters);

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:

    const MaterialProperty<RankFourTensor> & _C;
    //Given stress components
    std::string _C_matrix_name;
    MaterialProperty<DenseMatrix<Real>> & _V;
   
};
