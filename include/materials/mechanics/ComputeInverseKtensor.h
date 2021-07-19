//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
// Forward Declarations
class ComputeInverseKtensor;

//MOOSe includes
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

template <>
InputParameters validParams<ComputeInverseKtensor>();

class ComputeInverseKtensor : public Material
{
public:
  ComputeInverseKtensor(const InputParameters & parameters);
  
  //This class calculates the inverse of K tensor

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:

    const VariableGradient & _grad_eta ;

    //unit normal
    //const MaterialProperty<RealVectorValue> & _n;
 
    //Stiffness tensor of alpha and beta phase
    const MaterialProperty<RankFourTensor> & _stiffness_alpha;
    const MaterialProperty<RankFourTensor> & _stiffness_beta;
    const MaterialProperty<Real> & _h;
    
    MaterialProperty<RankTwoTensor> & _inv_K_val;
    MaterialProperty<RankTwoTensor> & _dK_dh_val;
};
