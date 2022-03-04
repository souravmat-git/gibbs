//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "ComputeSmallStrainAlpha.h"
registerMooseObject("gibbsApp", ComputeSmallStrainAlpha);

//template <>
InputParameters
ComputeSmallStrainAlpha::validParams()
{
  InputParameters params = ComputeStrainBase::validParams();
  params.addRequiredParam<std::string>("base_name","The phase name");
  //params.addRequiredParam<MaterialPropertyName>("strain_jump", "Jump in strain");
  //params.addRequiredParam<MaterialPropertyName>("h", "Interpolation Function");
  return params;
}

ComputeSmallStrainAlpha::ComputeSmallStrainAlpha(const InputParameters & parameters)
  : ComputeStrainBase(parameters),
   _strain_jump(getMaterialProperty<RankTwoTensor>("strain_jump")),
    //interpolation functions
   _h(getMaterialProperty<Real>("h"))
{
}

void
ComputeSmallStrainAlpha::computeQpProperties()
{
  //Define the displacement gradient tensor which is in general not symmetric
  const auto _disp_grad_tensor = RankTwoTensor::initializeFromRows((*_grad_disp[0])[_qp],
                                    (*_grad_disp[1])[_qp],
                                    (*_grad_disp[2])[_qp]);

  //Assuming infinitesimal deformation
  // total strain_in alpha phase = (total) strain +  strain jump * interpolation
  _total_strain[_qp] = (_disp_grad_tensor + _disp_grad_tensor.transpose()) / 2.0;


  if (_global_strain)
      _total_strain[_qp] += (*_global_strain)[_qp];

    _mechanical_strain[_qp] = _total_strain[_qp];

    // Remove the Eigen strain
    for (auto es : _eigenstrains)
      _mechanical_strain[_qp] -= (*es)[_qp];

  // + _h[_qp] * _strain_jump[_qp];
  _mechanical_strain[_qp] += _h[_qp] * _strain_jump[_qp];
}
