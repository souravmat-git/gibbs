//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EigenStrainPrefactor.h"
registerMooseObject("gibbsApp", EigenStrainPrefactor);

template<>
InputParameters
validParams<EigenStrainPrefactor>()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Calculates the prefactor required by" 
                             "EigenStrainPhaseMaterial");
  params.addRequiredParam<std::string>("base_name", "Phase name");
  params.addRequiredParam<Real>("R", "Gas constant");
  params.addRequiredParam<Real>("T", "Simulation temperature");
  params.addRequiredParam<Real>("overall_xB", "Overall mole fraction");
  params.addRequiredParam<Real>("xB_eqm", "Equilibrium mole fraction in any given phase");
  params.addRequiredParam<Real>("B_tf_eqm", "Eqm therm. factor of comp B in any given phase");
  params.addRequiredParam<Real>("B_diff_pot_eqm", "Equilibrium diffusion potential of comp B");
  params.addRequiredCoupledVar("B_diff_pot", "diffusion potential");
  return params;
}

EigenStrainPrefactor::EigenStrainPrefactor(const InputParameters & parameters)
 :Material(parameters),
 _base_name(getParam<std::string>("base_name")),
 //constant material parameters
 _R(getParam<Real>("R")),
 _T(getParam<Real>("T")),
 _xB_eqm(getParam<Real>("xB_eqm")),
 _B_tf_eqm(getParam<Real>("B_tf_eqm")),
 _B_diff_pot_eqm(getParam<Real>("B_diff_pot_eqm")),
 _xB_o(getParam<Real>("overall_xB")),
 _prefactor(declareProperty<Real>(_base_name + "_prefactor")),
 _dprefactor(declareProperty<Real>(_base_name + "_dprefactor")),
 _d2prefactor(declareProperty<Real>(_base_name + "_d2prefactor")),
  // coupled variables
 _B_diff_pot(coupledValue("B_diff_pot"))
{
}

void
EigenStrainPrefactor::computeQpProperties(){ 
  
  Real _Ec = _R*_T;
  
  Real _xB  = (_B_diff_pot[_qp] - (_B_diff_pot_eqm/_Ec))/(_B_tf_eqm/_Ec) + _xB_eqm;  
  Real _chi_B = (1.0/(_B_tf_eqm/_Ec));
  
  //Note that the phase-compositions are obtained as functions of diffusion potentaial
  _prefactor[_qp] = (_xB - _xB_o);
  _dprefactor[_qp] = _chi_B;
  _d2prefactor[_qp] = 0;

}
