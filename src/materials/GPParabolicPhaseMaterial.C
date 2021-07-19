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

#include "GPParabolicPhaseMaterial.h"
registerMooseObject("gibbsApp",GPParabolicPhaseMaterial);

template <>
InputParameters
validParams<GPParabolicPhaseMaterial>()
{
  InputParameters params = validParams<TabulatedPhaseMaterial>();
  params.addRequiredParam<MaterialPropertyName>("xB", "Mole fraction of comp B");
  params.addRequiredParam<MaterialPropertyName>("inv_B_tf", "Thermodynamic factor of comp B");
  params.addRequiredParam<MaterialPropertyName>("A_chem_pot", "Chemical potential of dependent comp");
  params.addRequiredParam<Real>("curvature", "Curvature of the parabola");
  params.addRequiredParam<Real>("horiz_shift", "Horizontal shift of the parabola");
  params.addParam<Real>("vert_shift", 0.0,"Vertical shift of the parabola");
  params.addRequiredCoupledVar("B_diff_pot", "diffusion potential");
  params.addClassDescription("This class approximates the grand-potential of a phase..."
                              "by a parabolic experession");
  return params;
}

GPParabolicPhaseMaterial::GPParabolicPhaseMaterial(const InputParameters & parameters)
  : TabulatedPhaseMaterial(parameters),
    _xB_name(getParam<MaterialPropertyName>("xB")),
    _inv_B_tf_name(getParam<MaterialPropertyName>("inv_B_tf")),
    _A_chem_pot_name(getParam<MaterialPropertyName>("A_chem_pot")),
    //constant material parameters
    _A(getParam<Real>("curvature")),
    _h(getParam<Real>("horiz_shift")),
    _k(getParam<Real>("vert_shift")),
    //material properties
    _xB_val(declareProperty<Real>(_xB_name)),
    _inv_B_tf_val(declareProperty<Real>(_inv_B_tf_name)),
    _A_chem_pot_val(declareProperty<Real>(_A_chem_pot_name)),
    // coupled variables
    _B_diff_pot(coupledValue("B_diff_pot"))
{
}

void 
GPParabolicPhaseMaterial::computeQpProperties()
{
    Real _mu = _B_diff_pot[_qp];
    
    //phase comp. alpha as functon of diffusion potential
    _xB_val[_qp] = (_mu/_A) + _h;
       
    //First derivative of phase composition alpha with respect to mu   
    _inv_B_tf_val[_qp] = (1.0/_A);
    
    //Grand-potential of alpha phase as function of diff_pot
    _A_chem_pot_val[_qp] = -0.5*((_mu*_mu)/_A) -_mu*_h + _k ; 
}
