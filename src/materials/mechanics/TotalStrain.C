//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "TotalStrain.h"
registerMooseObject("gibbsApp", TotalStrain);

template <>
InputParameters
validParams<TotalStrain>()
{
  InputParameters params = validParams<ComputeStrainBase>();
  return params;
}

TotalStrain::TotalStrain(const InputParameters & parameters)
  : ComputeStrainBase(parameters)
{
} 

void
TotalStrain::computeQpProperties(){  
   
   //Define the displacement gradient tensor which is in general not symmetric
    RankTwoTensor _disp_grad_tensor((*_grad_disp[0])[_qp], 
                                    (*_grad_disp[1])[_qp],  
                                    (*_grad_disp[2])[_qp]);
   
  //Assuming infinitesimal deformation, the total strain is given by 
  _total_strain[_qp] = (_disp_grad_tensor + _disp_grad_tensor.transpose()) / 2.0;
  
  if (_global_strain)
      _total_strain[_qp] += (*_global_strain)[_qp];

    _mechanical_strain[_qp] = _total_strain[_qp];

}
