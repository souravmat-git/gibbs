//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//This kernel is based on StressDivergenceTensors.C

#include "QMomentumBalance3DFun.h"
#include "Function.h"
registerMooseObject("gibbsApp", QMomentumBalance3DFun);

template <>
InputParameters
validParams<QMomentumBalance3DFun>(){
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the momentum balance");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction "
                                        "the variable this kernel acts in. (0 for x, "
                                        "1 for y, 2 for z)");
  params.addRequiredCoupledVar("displacements",
                   "The string of displacements suitable for the problem statement");
  params.addRequiredParam<FunctionName>("eta", "Phase field variable");
  return params;
}

QMomentumBalance3DFun::QMomentumBalance3DFun(const InputParameters & parameters)
  :Kernel(parameters),
  _eta(getFunction("eta")),
  //Global stress
  _stress(getMaterialProperty<RankTwoTensor>("stress")),
  //Jacobian_mult
  _Jacobian_mult(getMaterialProperty<RankFourTensor>("Jacobian_mult")),
  _component(getParam<unsigned int>("component")),
  _ndisp(coupledComponents("displacements")),
  _disp_var(_ndisp)
{
  for (unsigned int i=0; i < _ndisp; i++)
    _disp_var[i] = coupled("displacements", i);

}

Real
QMomentumBalance3DFun::computeQpResidual()
{
   return _grad_test[_i][_qp] * _stress[_qp].row(_component);     
}

Real
QMomentumBalance3DFun::computeQpJacobian()
{ 
  Real jacobian = 0.0;
  // B^T_i * C * B_j
  jacobian += ElasticityTensorTools::elasticJacobian(
      _Jacobian_mult[_qp], _component, _component, _grad_test[_i][_qp], _grad_phi[_j][_qp]);
      
   return jacobian;
}     

Real
QMomentumBalance3DFun::computeQpOffDiagJacobian(unsigned int jvar)
{   //Based on StressDivergeceTensors.C in MOOSE
   for (unsigned int disp_comp = 0; disp_comp < _ndisp; disp_comp++)
    if (jvar == _disp_var[disp_comp])
    {
    
        Real jacobian = 0.0;
        // B^T_i * C * B_j
        jacobian += ElasticityTensorTools::elasticJacobian(_Jacobian_mult[_qp],
                                                         _component,
                                                          disp_comp,
                                                         _grad_test[_i][_qp],
                                                         _grad_phi[_j][_qp]);
      
        return jacobian;
      } 
  return 0; 
}
