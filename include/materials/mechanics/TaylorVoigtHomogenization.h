//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html#pragma once
// Forward Declarations

class TaylorVoigtHomogenization;

//MOOSe includes
#include "Material.h"
#include "RankTwoTensor.h"

template <>
InputParameters validParams<TaylorVoigtHomogenization>();

class TaylorVoigtHomogenization : public Material
{
public:
 TaylorVoigtHomogenization(const InputParameters & parameters);
  
  //This class calculates the difference in compatible elastic strain
  //and supplies this value to other material class
  //Its current value depends on the local value of the jump in displacement
  //and the unit normal at that point

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:
   
    //Compute the second-rank tensor
    MaterialProperty<RankTwoTensor> & _strain_jump;  
    MaterialProperty<RankTwoTensor> & _dstrainjump_dphi;        
};
