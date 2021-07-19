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

#include "GPFittedMaterial.h"
registerMooseObject("gibbsApp",GPFittedMaterial);

template <>
InputParameters
validParams<GPFittedMaterial>()
{
  InputParameters params = validParams<TabulatedPhaseMaterial>();
  params.addRequiredParam<MaterialPropertyName>("xB", "Mole fraction of comp B");
  params.addRequiredParam<MaterialPropertyName>("inv_B_tf", "Thermodynamic factor of comp B");
  params.addRequiredParam<MaterialPropertyName>("A_chem_pot", "Chemical potential of dependent comp");
  params.addParam<MaterialPropertyName>("inv_B_td", 0.0, "Third derivative of comp B");
  params.addRequiredParam<Real>("a", "Parameter 1");
  params.addRequiredParam<Real>("b", "Parameter 2");
  params.addRequiredParam<Real>("c2", "Parameter 3");
  params.addRequiredParam<Real>("c3", "Parameter 4");
  params.addRequiredParam<Real>("c4", "Parameter 5");
  params.addRequiredParam<Real>("c5", "Parameter 6");
  params.addRequiredParam<Real>("xB_eqm", "Equilibrium mole fraction of comp B");
  params.addRequiredParam<Real>("B_diff_pot_eqm", "Eqm. diffusion potential");
  params.addRequiredCoupledVar("B_diff_pot", "diffusion potential");
  params.addClassDescription("This class approximates the free energy of a phase..."
                              "by a Taylor approximation");
  return params;
}

GPFittedMaterial::GPFittedMaterial(const InputParameters & parameters)
  : TabulatedPhaseMaterial(parameters),
    _xB_name(getParam<MaterialPropertyName>("xB")),
    _inv_B_tf_name(getParam<MaterialPropertyName>("inv_B_tf")),
    _inv_B_td_name(getParam<MaterialPropertyName>("inv_B_td")),
    _A_chem_pot_name(getParam<MaterialPropertyName>("A_chem_pot")),
    //constant material parameters
    _a(getParam<Real>("a")),
    _b(getParam<Real>("b")),
    _c2(getParam<Real>("c2")),
    _c3(getParam<Real>("c3")),
    _c4(getParam<Real>("c4")),
    _c5(getParam<Real>("c5")),
    _xB_eqm(getParam<Real>("xB_eqm")),
    _B_diff_pot_eqm(getParam<Real>("B_diff_pot_eqm")),
    //material properties
    _xB_val(declareProperty<Real>(_xB_name)),
    _inv_B_tf_val(declareProperty<Real>(_inv_B_tf_name)),
    _inv_B_td_val(declareProperty<Real>(_inv_B_td_name)),
    _A_chem_pot_val(declareProperty<Real>(_A_chem_pot_name)),
    // coupled variables
    _B_diff_pot(coupledValue("B_diff_pot"))
{
}

void 
GPFittedMaterial::computeQpProperties()
{
    //Note that _B_diff_pot is assumed to be non-dimensional
    
    Real x = (_B_diff_pot[_qp] - (_B_diff_pot_eqm/_Ec));
    
    //phase composition as functon of diffusion potential
    _xB_val[_qp] = -(-_xB_eqm - _b*(x/(std::sqrt(1 + x*x))) 
                   - (_c2/_Ec)*x - (_c3/_Ec)*x*x - (_c4/_Ec)*std::pow(x,3.0) 
                   - (_c5/_Ec)*std::pow(x,4.0));
    
    //First derivative of phase composition alpha with respect to mu   
    _inv_B_tf_val[_qp] = -( -_b*(1.0/((1.0+x*x)*sqrt(1.0 + x*x))) -(_c2/_Ec)
                          -2.0*(_c3/_Ec)*x - 3.0*(_c4/_Ec)*x*x - 4.0*(_c5/_Ec)*std::pow(x,3.0));
    
    _inv_B_td_val[_qp] = 0.0;
    
    //Grand-potential of any phase as a function of diffusion potential
    //Note that _Ec is for non-dimensionalization
    _A_chem_pot_val[_qp] =  (_a/_Ec - _xB_eqm*(x) -_b*std::sqrt(1.0 + x*x) 
                          - 0.5*(_c2/_Ec)*std::pow(x,2.0) -(1/3)*(_c3/_Ec)*std::pow(x,3.0)
                          - (1.0/4.0)*(_c4/_Ec)*std::pow(x,4.0) - (1/5.0)*(_c5/_Ec)*std::pow(x,5.0)) ;
        
}
