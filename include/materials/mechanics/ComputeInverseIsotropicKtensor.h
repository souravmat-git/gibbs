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
class ComputeInverseIsotropicKtensor;

//MOOSe includes
#include "Material.h"
#include "RankTwoTensor.h"

template <>
InputParameters validParams<ComputeInverseIsotropicKtensor>();

class ComputeInverseIsotropicKtensor : public Material
{
public:
  ComputeInverseIsotropicKtensor(const InputParameters & parameters);
  
  //This class calculates the jump in strain or the strain difference
  //and its derivative wrt to overall strain

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:

    //unit normal
    const MaterialProperty<RealVectorValue> & _n;
 
     //Lames parameters for alpha and beta phases
    const MaterialProperty<Real> & _alpha_lambda;
    const MaterialProperty<Real> & _beta_lambda;
  
    const MaterialProperty<Real> & _alpha_mu;
    const MaterialProperty<Real> & _beta_mu;
    
    //interpolation function
    const MaterialProperty<Real> & _h;
    
    MaterialProperty<RankTwoTensor> & _inv_K_val;
    MaterialProperty<RankTwoTensor> & _dK_dh_val;
};
