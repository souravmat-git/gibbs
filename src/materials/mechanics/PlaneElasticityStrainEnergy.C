//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "PlaneElasticityStrainEnergy.h"
registerMooseObject("gibbsApp", PlaneElasticityStrainEnergy);

//template <>
InputParameters
PlaneElasticityStrainEnergy::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>("fel", "Elastic strain energy density");
  params.addRequiredParam<std::string>("phase_name", "phase for which the elastic strain energy");
  return params;
}

PlaneElasticityStrainEnergy::PlaneElasticityStrainEnergy(const InputParameters & parameters)
 : Material(parameters),
  _phase_name(getParam<std::string>("phase_name")),
  _fel_name(getParam<MaterialPropertyName>("fel")),
  _sxx(getMaterialProperty<Real>("sx_" + _phase_name)),
  _syy(getMaterialProperty<Real>("sy_" + _phase_name)),
  _sxy(getMaterialProperty<Real>("sxy_" + _phase_name)),
  _exx(getMaterialProperty<Real>("ex_" + _phase_name)),
  _eyy(getMaterialProperty<Real>("ey_" + _phase_name)),
  _exy(getMaterialProperty<Real>("exy_" + _phase_name)),
  _fel_val(declareProperty<Real>(_fel_name))
{
}

void
PlaneElasticityStrainEnergy::computeQpProperties(){

  _fel_val[_qp] = 0.5*(_sxx[_qp] * _exx[_qp]
                     + _syy[_qp] * _eyy[_qp]
               + 2.0 * _sxy[_qp] * _exy[_qp]);

}
