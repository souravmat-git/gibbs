//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//* This Material class supplies the interpolation function,
//* and its first and second derivatives

#include "ParabolicPhaseMaterial.h"
registerMooseObject("gibbsApp", ParabolicPhaseMaterial);

//template <>
InputParameters
ParabolicPhaseMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>("free_energy", "Free energy of the given phase");
  params.addRequiredParam<MaterialPropertyName>("B_diff_pot", "Diffusion potential of comp B");
  params.addRequiredParam<MaterialPropertyName>("B_therm_factor","Thermodyanmic factor");
  params.addRequiredCoupledVar("mol_fraction_B", "Mole fraction of the solute component");
  params.addRequiredParam<Real>("curvature", "curvature of the parabola");
  params.addRequiredParam<Real>("horiz_shift", "horizontal shift of the parabola");
  params.addParam<Real>("vert_shift", 0.0,"vertical shift of the parabola");
  params.addClassDescription("This class returns the above declared prop..."
                              "assuming a vertex form of the parabola");
  return params;
}

ParabolicPhaseMaterial::ParabolicPhaseMaterial(const InputParameters & parameters)
  : Material(parameters),
  _f_phase_name(getParam<MaterialPropertyName>("free_energy")),
  _B_diff_pot_phase_name(getParam<MaterialPropertyName>("B_diff_pot")),
  _B_therm_factor_phase_name(getParam<MaterialPropertyName>("B_therm_factor")),
  _f_phase_val(declareProperty<Real>(_f_phase_name)),
  _B_diff_pot_phase_val(declareProperty<Real>(_B_diff_pot_phase_name)),
  _B_therm_factor_phase_val(declareProperty<Real>(_B_therm_factor_phase_name)),
  _mol_fraction_B(coupledValue("mol_fraction_B")),
  _A(getParam<Real>("curvature")),
  _h(getParam<Real>("horiz_shift")),
  _k(getParam<Real>("vert_shift"))
{
}

void
ParabolicPhaseMaterial::computeQpProperties()
{
  //Parbolic free energy
  _f_phase_val[_qp] = 0.5*_A*std::pow((_mol_fraction_B[_qp] - _h),2.0) + _k;

  //Linear diffusion potential
  _B_diff_pot_phase_val[_qp] = _A*(_mol_fraction_B[_qp] - _h);

  //Constant therm factor
  _B_therm_factor_phase_val[_qp] = _A;

}
