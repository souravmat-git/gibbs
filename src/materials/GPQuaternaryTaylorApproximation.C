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

#include "GPQuaternaryTaylorApproximation.h"
registerMooseObject("gibbsApp",GPQuaternaryTaylorApproximation);

template <>
InputParameters
validParams<GPQuaternaryTaylorApproximation>()
{
  InputParameters params = validParams<GPTernaryTaylorApproximation>();
  params.addRequiredParam<MaterialPropertyName>("xD", "Mole fraction of comp C");
  params.addRequiredParam<MaterialPropertyName>("inv_D_tf", "Thermodynamic factor of comp C");
  params.addRequiredParam<Real>("xD_eqm", "Equilibrium mole fraction of comp C");
  params.addRequiredParam<Real>("D_tf_eqm", "Eqm therm. factor of comp C");
  params.addRequiredParam<Real>("D_diff_pot_eqm", "Equilibrium diffusion potential of comp B");
  params.addRequiredCoupledVar("D_diff_pot", "diffusion potential");
  params.addClassDescription("This class approximates the free energy of a phase..."
                              "by a Taylor approximation");
  return params;
}

GPQuaternaryTaylorApproximation::GPQuaternaryTaylorApproximation(const InputParameters & parameters)
  : GPTernaryTaylorApproximation(parameters),
    _xD_name(getParam<MaterialPropertyName>("xD")),
    _inv_D_tf_name(getParam<MaterialPropertyName>("inv_D_tf")),
    //constant material parameters
    _xD_eqm(getParam<Real>("xD_eqm")),
    _D_tf_eqm(getParam<Real>("D_tf_eqm")),
    _D_diff_pot_eqm(getParam<Real>("D_diff_pot_eqm")),
    //material properties
    _xD_val(declareProperty<Real>(_xD_name)),
    _inv_D_tf_val(declareProperty<Real>(_inv_D_tf_name)),
    // coupled variables
    _D_diff_pot(coupledValue("D_diff_pot"))
{
}

void 
GPQuaternaryTaylorApproximation::computeQpProperties()
{
    
    //phase comp. alpha as functon of diffusion potential
    _xB_val[_qp] = (_B_diff_pot[_qp] - (_B_diff_pot_eqm/_Ec))/(_B_tf_eqm/_Ec) + _xB_eqm;
    _xC_val[_qp] = (_C_diff_pot[_qp] - (_C_diff_pot_eqm/_Ec))/(_C_tf_eqm/_Ec) + _xC_eqm;
    _xD_val[_qp] = (_D_diff_pot[_qp] - (_D_diff_pot_eqm/_Ec))/(_D_tf_eqm/_Ec) + _xD_eqm;
    
    //First derivative of phase composition alpha with respect to mu   
    _inv_B_tf_val[_qp] = (1.0/(_B_tf_eqm/_Ec));
    _inv_C_tf_val[_qp] = (1.0/(_C_tf_eqm/_Ec));
    _inv_D_tf_val[_qp] = (1.0/(_D_tf_eqm/_Ec));
    
    //Grand-potential of alpha phase as function of diff_pot
    _A_chem_pot_val[_qp] =  ((_A_chem_pot_eqm/_Ec) - _xB_eqm*(_B_diff_pot[_qp] - (_B_diff_pot_eqm/_Ec))
                                                   - _xC_eqm*(_C_diff_pot[_qp] - (_C_diff_pot_eqm/_Ec))
                                                   - _xD_eqm*(_D_diff_pot[_qp] - (_D_diff_pot_eqm/_Ec))                          
                            -(0.5/(_B_tf_eqm/_Ec))*std::pow((_B_diff_pot[_qp] - (_B_diff_pot_eqm/_Ec)),2.0)
                            -(0.5/(_C_tf_eqm/_Ec))*std::pow((_C_diff_pot[_qp] - (_C_diff_pot_eqm/_Ec)),2.0)
                            -(0.5/(_D_tf_eqm/_Ec))*std::pow((_D_diff_pot[_qp] - (_D_diff_pot_eqm/_Ec)),2.0));

}
