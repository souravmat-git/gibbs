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

#include "GPTernaryTaylorApproximation.h"
registerMooseObject("gibbsApp",GPTernaryTaylorApproximation);

//template <>
InputParameters
GPTernaryTaylorApproximation::validParams()
{
  InputParameters params = GPTaylorApproximation::validParams();
  params.addRequiredParam<MaterialPropertyName>("xC", "Mole fraction of comp C");
  params.addRequiredParam<MaterialPropertyName>("inv_C_tf", "Thermodynamic factor of comp C");
  params.addRequiredParam<Real>("xC_eqm", "Equilibrium mole fraction of comp C");
  params.addRequiredParam<Real>("C_tf_eqm", "Eqm therm. factor of comp C");
  params.addRequiredParam<Real>("C_diff_pot_eqm", "Equilibrium diffusion potential of comp B");
  params.addRequiredCoupledVar("C_diff_pot", "diffusion potential");
  params.addClassDescription("This class approximates the free energy of a phase..."
                              "by a Taylor approximation");
  return params;
}

GPTernaryTaylorApproximation::GPTernaryTaylorApproximation(const InputParameters & parameters)
  : GPTaylorApproximation(parameters),
    _xC_name(getParam<MaterialPropertyName>("xC")),
    _inv_C_tf_name(getParam<MaterialPropertyName>("inv_C_tf")),
    //constant material parameters
    _xC_eqm(getParam<Real>("xC_eqm")),
    _C_tf_eqm(getParam<Real>("C_tf_eqm")),
    _C_diff_pot_eqm(getParam<Real>("C_diff_pot_eqm")),
    //material properties
    _xC_val(declareProperty<Real>(_xC_name)),
    _inv_C_tf_val(declareProperty<Real>(_inv_C_tf_name)),
    // coupled variables
    _C_diff_pot(coupledValue("C_diff_pot"))
{
}

void
GPTernaryTaylorApproximation::computeQpProperties()
{

    //phase comp. alpha as functon of diffusion potential
    _xB_val[_qp] = (_B_diff_pot[_qp] - (_B_diff_pot_eqm/_Ec))/(_B_tf_eqm/_Ec) + _xB_eqm;
    _xC_val[_qp] = (_C_diff_pot[_qp] - (_C_diff_pot_eqm/_Ec))/(_C_tf_eqm/_Ec) + _xC_eqm;

    //First derivative of phase composition alpha with respect to mu
    _inv_B_tf_val[_qp] = (1.0/(_B_tf_eqm/_Ec));
    _inv_C_tf_val[_qp] = (1.0/(_C_tf_eqm/_Ec));

    //Grand-potential of alpha phase as function of diff_pot
    _A_chem_pot_val[_qp] =  ((_A_chem_pot_eqm/_Ec) - _xB_eqm*(_B_diff_pot[_qp] - (_B_diff_pot_eqm/_Ec))
                                                   - _xC_eqm*(_C_diff_pot[_qp] - (_C_diff_pot_eqm/_Ec))
                            -(0.5/(_B_tf_eqm/_Ec))*std::pow((_B_diff_pot[_qp] - (_B_diff_pot_eqm/_Ec)),2.0)
                            -(0.5/(_C_tf_eqm/_Ec))*std::pow((_C_diff_pot[_qp] - (_C_diff_pot_eqm/_Ec)),2.0)) ;
    //std::cout << _A_chem_pot_val[_qp] << "\n";
}
