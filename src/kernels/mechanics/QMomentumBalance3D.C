//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* This was written by S.Chatterjee

#include "QMomentumBalance3D.h"
#include "ElasticityTensorTools.h"

registerMooseObject("gibbsApp", QMomentumBalance3D);

template <>
InputParameters
validParams<QMomentumBalance3D>(){
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Based on stress diveregence kernel");
  params.addRequiredCoupledVar("eta", "Phase field variable");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction "
                                        "the variable this kernel acts in. (0 for x, "
                                        "1 for y, 2 for z)");
  params.addRequiredCoupledVar("displacements",
                   "The string of displacements suitable for the problem statement");
  return params;
}

QMomentumBalance3D::QMomentumBalance3D(const InputParameters & parameters)
  :Kernel(parameters),
  //_eta(coupledValue("eta")),
  _eta_var(coupled("eta")),
   //Global stress
  _stress(getMaterialProperty<RankTwoTensor>("stress")),
  _dsigma_dphi(getMaterialProperty<RankTwoTensor>("dsigma_dphi")),
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
QMomentumBalance3D::computeQpResidual()
{
    //_stress[_qp].row(_component=0) returns s01, s02, s03
    //The * operator takes a dot product between the gradient of 
    // the test function and the stress tensor vector
    return _grad_test[_i][_qp] * _stress[_qp].row(_component);              
}

Real
QMomentumBalance3D::computeQpJacobian()
{  
  Real jacobian = 0.0;
  // B^T_i * C * B_j
  jacobian += ElasticityTensorTools::elasticJacobian(
      _Jacobian_mult[_qp], _component, _component, _grad_test[_i][_qp], _grad_phi[_j][_qp]);
      
   return jacobian;
}     

Real
QMomentumBalance3D::computeQpOffDiagJacobian(unsigned int jvar)
{
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
   
  if (jvar == _eta_var)
    return _grad_test[_i][_qp] * _dsigma_dphi[_qp].row(_component) * _phi[_j][_qp];  
 
  return 0;
}
