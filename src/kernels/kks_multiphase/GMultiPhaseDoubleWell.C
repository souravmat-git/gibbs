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

#include "GMultiPhaseDoubleWell.h"

registerMooseObject("gibbsApp", GMultiPhaseDoubleWell);

template <>
InputParameters
validParams<GMultiPhaseDoubleWell>()
{
  InputParameters params = ACBulk<Real>::validParams();
  params.addClassDescription("Multi-double well barrier term");
  params.addRequiredCoupledVar("eta2", "Coupled phase field variable 2");
  params.addCoupledVar("eta3", 0.0, "Coupled phase field variable 3");
  params.addCoupledVar("eta4", 0.0, "Coupled phase field variable 4");
  params.addCoupledVar("eta5", 0.0, "Coupled phase field variable 5");
  params.addRequiredParam<MaterialPropertyName>("gamma","Interface energy parameter");
  params.addRequiredParam<MaterialPropertyName>("barrier_height", "Height of the double well potential");
  return params;
}

GMultiPhaseDoubleWell::GMultiPhaseDoubleWell(const InputParameters & parameters)
  : ACBulk<Real>(parameters),
   _eta2(coupledValue("eta2")),
   _eta2_var(coupled("eta2")),
   _eta3(coupledValue("eta3")),
   _eta3_var(coupled("eta3")),
   _eta4(coupledValue("eta4")),
   _eta4_var(coupled("eta4")),
   _eta5(coupledValue("eta5")),
   _eta5_var(coupled("eta5")),
   _gamma(getMaterialProperty<Real>("gamma")),
   _BH(getMaterialProperty<Real>("barrier_height"))
{
}

Real
GMultiPhaseDoubleWell::computeDFDOP(PFFunctionType type)
{

  switch (type)
  {
    case Residual:
    {
      return (_BH[_qp]*( (std::pow(_u[_qp],3.0) -_u[_qp] 
              + 2.0*_u[_qp] * _gamma[_qp] * (_eta2[_qp]  * _eta2[_qp] 
                                          + _eta3[_qp] * _eta3[_qp] 
                                          + _eta4[_qp] * _eta4[_qp] 
                                          + _eta5[_qp] * _eta5[_qp]) ) ));
                                          
    }   
    case Jacobian:
    {
      return (_phi[_j][_qp]* (_BH[_qp]*((3.0*std::pow(_u[_qp],2.0) - 1.0
                 + 2.0* _gamma[_qp] * (_eta2[_qp] * _eta2[_qp] + _eta3[_qp] * _eta3[_qp] 
                                  + _eta4[_qp]* _eta4[_qp]+ _eta5[_qp] * _eta5[_qp]) )) ));                                   
    }
  }
  
  mooseError("Internal error");
}

Real
GMultiPhaseDoubleWell::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _eta2_var)
  {
    return ACBulk<Real>::computeQpOffDiagJacobian(jvar) + 
           _test[_i][_qp]* _BH[_qp]* _L[_qp]*2.0* _u[_qp]* _gamma[_qp]* (2.0* _eta2[_qp]) * _phi[_j][_qp];
  }
  
  else if (jvar == _eta3_var)
  {
    return ACBulk<Real>::computeQpOffDiagJacobian(jvar) + 
          _test[_i][_qp]* _BH[_qp]* _L[_qp]*2.0* _u[_qp]* _gamma[_qp]*(2.0* _eta3[_qp]) * _phi[_j][_qp]; 
  
  }
  else if (jvar == _eta4_var)
  {
    return ACBulk<Real>::computeQpOffDiagJacobian(jvar) + 
    _test[_i][_qp]* _BH[_qp]* _L[_qp]*2.0* _u[_qp]* _gamma[_qp]*(2.0* _eta4[_qp]) * _phi[_j][_qp]; 
  
  }
  else if (jvar == _eta5_var)
  {
    return ACBulk<Real>::computeQpOffDiagJacobian(jvar) + 
      _test[_i][_qp]* _BH[_qp]* _L[_qp]*2.0* _u[_qp]* _gamma[_qp]*(2.0* _eta5[_qp]) * _phi[_j][_qp]; 
  
  }
  else 
  return 0.0; 
}
