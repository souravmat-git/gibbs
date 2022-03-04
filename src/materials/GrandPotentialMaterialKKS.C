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

#include "GrandPotentialMaterialKKS.h"
registerMooseObject("gibbsApp",GrandPotentialMaterialKKS);

//template <>
InputParameters
GrandPotentialMaterialKKS::validParams()
{
  InputParameters params = Material::validParams();
  params.addCoupledVar("B_diff_pot", "diffusion potential");
  return params;
}

GrandPotentialMaterialKKS::GrandPotentialMaterialKKS(const InputParameters & parameters)
  : Material(parameters),
    // Declare that this material is going to provide a Real
    // valued properties named "_f,_df,_d2f" that Kernels can use.
    _xB_alpha(declareProperty<Real>("xB_alpha")),
    _xB_beta(declareProperty<Real>("xB_beta")),
    _inv_B_tf_alpha(declareProperty<Real>("inv_B_tf_alpha")),
    _inv_B_tf_beta(declareProperty<Real>("inv_B_tf_beta")),
    _A_chem_pot_alpha(declareProperty<Real>("A_chem_pot_alpha")),
    _A_chem_pot_beta(declareProperty<Real>("A_chem_pot_beta")),
    // coupled variables
    _mu(coupledValue("B_diff_pot"))
{
}

//Function

void
GrandPotentialMaterialKKS::computeQpProperties()
{
    Real A_alpha = 2.0;
    Real eqm_comp_alpha = 0.1;

    Real A_beta = 2.0;
    Real eqm_comp_beta = 0.9;

    //phase comp. alpha as functon of diffusion potential
    _xB_alpha[_qp] = (_mu[_qp]/A_alpha + eqm_comp_alpha);

    //phase comp.beta as function of diffusion potential
    _xB_beta[_qp] = (_mu[_qp]/A_beta + eqm_comp_beta);

    //First derivative of phase composition alpha with respect to mu
    _inv_B_tf_alpha[_qp] = (1.0/A_alpha);

    //First derivative of phase composition beta with respect to mu
    _inv_B_tf_beta[_qp] = (1.0/A_beta);

    //Grand-potential of alpha phase as function of diff_pot
    _A_chem_pot_alpha[_qp] = -(0.5*(_mu[_qp]*_mu[_qp]/A_alpha) + _mu[_qp]*eqm_comp_alpha);

    //Grand-potential of beta phase as function of diff_pot
    _A_chem_pot_beta[_qp] = -(0.5*(_mu[_qp]*_mu[_qp]/A_beta) + _mu[_qp]*eqm_comp_beta);


}
