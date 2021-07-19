//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html 

//* This kernel implements the driving force 
//* which is the difference in elastic strain energy between the two phases

#include "NFDrivingForce.h"
registerMooseObject("gibbsApp", NFDrivingForce);

template<>
InputParameters
validParams<NFDrivingForce>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: dh*[f_beta_el -f_alpha_el] = 0");
  params.addRequiredCoupledVar("ex_alpha", "Strain in the alpha phase");
  params.addRequiredCoupledVar("ex_beta" ,  "Strain in the beta phase");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
  params.addRequiredParam<MaterialPropertyName>("nd_factor", "Non-dimensional factor");
  return params;
}

NFDrivingForce::NFDrivingForce(const InputParameters & parameters)
  :Kernel(parameters),
  //Coupled variables
  _ex_alpha(coupledValue("ex_alpha")),
  _ex_alpha_var(coupled("ex_alpha")),
  _ex_beta(coupledValue("ex_beta")),
  _ex_beta_var(coupled("ex_beta")),
  //Interpolation function
  _dh(getMaterialProperty<Real>("dh")),
  _d2h(getMaterialProperty<Real>("d2h")),
  //Requisite material properties
  _fel_alpha(getMaterialProperty<Real>("fel_alpha")),
  _fel_beta(getMaterialProperty<Real>("fel_beta")),
  _sx_alpha(getMaterialProperty<Real>("sx_alpha")),
  _sx_beta(getMaterialProperty<Real>("sx_beta")),
  _mat_const_alpha(getMaterialProperty<Real>("mat_const_alpha")),
  _mat_const_beta(getMaterialProperty<Real>("mat_const_beta")),
  _L(getMaterialProperty<Real>("mob_name")),
  _nd_factor(getMaterialProperty<Real>("nd_factor"))
{
}

Real
NFDrivingForce::QuantMechDrivingForce() const{
  
  const Real lc_beta  = _fel_beta[_qp]  - _sx_beta[_qp]  * _ex_beta[_qp];
  const Real lc_alpha = _fel_alpha[_qp] - _sx_alpha[_qp] * _ex_alpha[_qp];
  
  return (lc_beta -lc_alpha); 
}


Real
NFDrivingForce::computeQpResidual(){
  //Difference in elastic complementary strain density is the driving force of PT
  //Residual: N_{i} * h_{\phi} (f_el_beta - f_el_alpha) = 0 ;
  return _test[_i][_qp] *_L[_qp] * _dh[_qp] * 
          NFDrivingForce::QuantMechDrivingForce() * _nd_factor[_qp];

}

Real
NFDrivingForce::computeQpJacobian(){
  // This will be non-zero since the 
  // free energies of each phase is dependent strain
  return _test[_i][_qp] * _L[_qp] * _d2h[_qp] * 
          NFDrivingForce::QuantMechDrivingForce() * _phi[_j][_qp] * _nd_factor[_qp];
}

Real
NFDrivingForce::computeQpOffDiagJacobian(unsigned int jvar){

  if (jvar == _ex_alpha_var){
    return _test[_i][_qp] * _L[_qp] * _dh[_qp] *
           ( _ex_alpha[_qp] * _mat_const_alpha[_qp])* _phi[_j][_qp] * _nd_factor[_qp];
  }
  else if (jvar == _ex_beta_var){
    return -_test[_i][_qp] * _L[_qp] * _dh[_qp] *
            ( _ex_beta[_qp] * _mat_const_beta[_qp])* _phi[_j][_qp] * _nd_factor[_qp];
  }
  else
    return 0;
}
