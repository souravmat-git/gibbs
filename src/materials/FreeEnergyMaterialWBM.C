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

#include "FreeEnergyMaterialWBM.h"
registerMooseObject("gibbsApp",FreeEnergyMaterialWBM);

template <>
InputParameters
validParams<FreeEnergyMaterialWBM>()
{
  InputParameters params = validParams<Material>();
  params.addCoupledVar("c", "Composition of the solute");
  return params;
}

FreeEnergyMaterialWBM::FreeEnergyMaterialWBM(const InputParameters & parameters)
  : Material(parameters),
    // Declare that this material is going to provide a Real
    // valued properties named "_f,_df,_d2f" that Kernels can use.
    _free_energy_alpha(declareProperty<Real>("f_alpha")),
    _dfalpha_dc(declareProperty<Real>("dfalpha_dc")),
    _d2falpha_dc2(declareProperty<Real>("d2falpha_dc2")),
    _free_energy_beta(declareProperty<Real>("f_beta")),
    _dfbeta_dc(declareProperty<Real>("dfbeta_dc")),
    _d2fbeta_dc2(declareProperty<Real>("d2fbeta_dc2")),
    // coupled variables
    _comp(coupledValue("c"))
{
}

//Function

void 
FreeEnergyMaterialWBM::computeQpProperties()
{
    Real A_alpha = 2.0;
    Real eqm_comp_alpha = 0.1;
    
    Real A_beta = 2.0;
    Real eqm_comp_beta = 0.9;


  _free_energy_alpha[_qp] = (A_alpha/2.0)*(_comp[_qp] - eqm_comp_alpha)*(_comp[_qp] - eqm_comp_alpha);
  _dfalpha_dc[_qp] = A_alpha*(_comp[_qp] - eqm_comp_alpha);
  _d2falpha_dc2[_qp] = A_alpha;
  
   _free_energy_beta[_qp] = (A_beta/2.0)*(_comp[_qp] - eqm_comp_beta)*(_comp[_qp] - eqm_comp_beta);
  _dfbeta_dc[_qp] = A_beta*(_comp[_qp] - eqm_comp_beta);
  _d2fbeta_dc2[_qp] = A_beta;
  
}
