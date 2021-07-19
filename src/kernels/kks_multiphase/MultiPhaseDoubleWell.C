//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//* This was written by S.Chatterjee
//f(phi) = phi^(2) * (1-phi)^(2)

#include "MultiPhaseDoubleWell.h"

registerMooseObject("gibbsApp", MultiPhaseDoubleWell);

template <>
InputParameters
validParams<MultiPhaseDoubleWell>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the barrier term"
                              "Eqn: H*dg/dphi");
  params.addRequiredCoupledVar("eta2", "Coupled phase field variable 2");
  params.addCoupledVar("eta3", 0.0, "Coupled phase field variable 3");
  params.addCoupledVar("eta4", 0.0, "Coupled phase field variable 4");
  params.addCoupledVar("eta5", 0.0, "Coupled phase field variable 5");
  params.addRequiredParam<MaterialPropertyName>("gamma","Interface energy parameter");
  params.addRequiredParam<MaterialPropertyName>("barrier_height", "Height of the double well potential");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "The mobility used with the kernel");

  return params;
}

MultiPhaseDoubleWell::MultiPhaseDoubleWell(const InputParameters & parameters)
  : Kernel(parameters),
   _eta2(coupledValue("eta2")),
   _eta2_var(coupled("eta2")),
   _eta3(coupledValue("eta3")),
   _eta3_var(coupled("eta3")),
   _eta4(coupledValue("eta4")),
   _eta4_var(coupled("eta4")),
   _eta5(coupledValue("eta5")),
   _eta5_var(coupled("eta5")),
   _gamma(getMaterialProperty<Real>("gamma")),
   _BH(getMaterialProperty<Real>("barrier_height")),
   _L(getMaterialProperty<Real>("mob_name"))
{
}

Real
MultiPhaseDoubleWell::computeQpResidual()
{
  return _L[_qp] * _test[_i][_qp] * _BH[_qp]*( (std::pow(_u[_qp],3.0) -_u[_qp] 
            + 2.0*_u[_qp] * _gamma[_qp] * (_eta2[_qp]  * _eta2[_qp] 
                                          + _eta3[_qp] * _eta3[_qp] 
                                          + _eta4[_qp] * _eta4[_qp] 
                                          + _eta5[_qp] * _eta5[_qp]) ) );
}

Real
MultiPhaseDoubleWell::computeQpJacobian()
{
  return  _L[_qp] * _test[_i][_qp] *_BH[_qp] *_phi[_j][_qp]*((3.0*std::pow(_u[_qp],2.0) - 1.0
            + 2.0* _gamma[_qp] * (_eta2[_qp] * _eta2[_qp] + _eta3[_qp] * _eta3[_qp] 
                                  + _eta4[_qp]* _eta4[_qp]+ _eta5[_qp] * _eta5[_qp])));
}

Real
MultiPhaseDoubleWell::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _eta2_var)
  {
    return _test[_i][_qp]* _BH[_qp]* _L[_qp]*2.0* _u[_qp]* _gamma[_qp]* (2.0* _eta2[_qp]) * _phi[_j][_qp];
  }
  
  else if (jvar == _eta3_var)
  {
    return _test[_i][_qp]* _BH[_qp]* _L[_qp]*2.0* _u[_qp]* _gamma[_qp]*(2.0* _eta3[_qp]) * _phi[_j][_qp]; 
  
  }
  else if (jvar == _eta4_var)
  {
    return _test[_i][_qp]* _BH[_qp]* _L[_qp]*2.0* _u[_qp]* _gamma[_qp]*(2.0* _eta4[_qp]) * _phi[_j][_qp]; 
  
  }
  else if (jvar == _eta5_var)
  {
    return _test[_i][_qp]* _BH[_qp]* _L[_qp]*2.0* _u[_qp]* _gamma[_qp]*(2.0* _eta5[_qp]) * _phi[_j][_qp]; 
  
  }
  else 
  return 0.0; 
}
