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

#include "GPTaylorApproximation.h"
registerMooseObject("gibbsApp", GPTaylorApproximation);

//template <>
InputParameters
GPTaylorApproximation::validParams()
{
  InputParameters params = TabulatedPhaseMaterial::validParams();
  params.addRequiredParam<MaterialPropertyName>("xB", "Mole fraction of comp B");
  params.addRequiredParam<MaterialPropertyName>("inv_B_tf", "Thermodynamic factor of comp B");
  params.addRequiredParam<MaterialPropertyName>("A_chem_pot", "Chemical potential of dependent comp");
  params.addRequiredParam<Real>("xB_eqm", "Equilibrium mole fraction in alpha phase");
  params.addRequiredParam<Real>("B_tf_eqm", "Eqm therm. factor of comp B");
  params.addRequiredParam<Real>("B_diff_pot_eqm", "Equilibrium diffusion potential of comp B");
  params.addRequiredParam<Real>("A_chem_pot_eqm", "Equilibrium chemical potential of comp A");
  params.addRequiredCoupledVar("B_diff_pot", "diffusion potential");
  params.addClassDescription("This class approximates the free energy of a phase..."
                              "by a Taylor approximation");
  return params;
}

GPTaylorApproximation::GPTaylorApproximation(const InputParameters & parameters)
  : TabulatedPhaseMaterial(parameters),
    _xB_name(getParam<MaterialPropertyName>("xB")),
    _inv_B_tf_name(getParam<MaterialPropertyName>("inv_B_tf")),
    _A_chem_pot_name(getParam<MaterialPropertyName>("A_chem_pot")),
    //constant material parameters
    _xB_eqm(getParam<Real>("xB_eqm")),
    _B_tf_eqm(getParam<Real>("B_tf_eqm")),
    _B_diff_pot_eqm(getParam<Real>("B_diff_pot_eqm")),
    _A_chem_pot_eqm(getParam<Real>("A_chem_pot_eqm")),
    //material properties
    _xB_val(declareProperty<Real>(_xB_name)),
    _inv_B_tf_val(declareProperty<Real>(_inv_B_tf_name)),
    _A_chem_pot_val(declareProperty<Real>(_A_chem_pot_name)),
    // coupled variables
    _B_diff_pot(coupledValue("B_diff_pot"))
{
}

void
GPTaylorApproximation::computeQpProperties()
{

    //phase comp. alpha as functon of diffusion potential
    _xB_val[_qp] = (_B_diff_pot[_qp] - (_B_diff_pot_eqm/_Ec))/(_B_tf_eqm/_Ec) + _xB_eqm;

    //First derivative of phase composition alpha with respect to mu
    _inv_B_tf_val[_qp] = (1.0/(_B_tf_eqm/_Ec));

    //Grand-potential of alpha phase as function of diff_pot
    _A_chem_pot_val[_qp] =  ((_A_chem_pot_eqm/_Ec) - _xB_eqm*(_B_diff_pot[_qp] - (_B_diff_pot_eqm/_Ec))
                            -(0.5/(_B_tf_eqm/_Ec))*std::pow((_B_diff_pot[_qp] - (_B_diff_pot_eqm/_Ec)),2.0)) ;

}
