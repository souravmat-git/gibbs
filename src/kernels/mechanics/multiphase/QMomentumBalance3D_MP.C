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

#include "QMomentumBalance3D_MP.h"
#include "ElasticityTensorTools.h"

registerMooseObject("gibbsApp", QMomentumBalance3D_MP);

template <>
InputParameters
validParams<QMomentumBalance3D_MP>(){
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Based on stress diveregence kernel");
  params.addRequiredCoupledVar("phase_alpha", "Phase field variable for alpha phase");
  params.addRequiredCoupledVar("phase_beta",  "Phase field variable for beta phase");
  params.addRequiredCoupledVar("phase_gamma", "Phase field variable for gamma phase");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction "
                                        "the variable this kernel acts in. (0 for x, "
                                        "1 for y, 2 for z)");
  params.addRequiredCoupledVar("displacements",
                   "The string of displacements suitable for the problem statement");
  return params;
}

QMomentumBalance3D_MP::QMomentumBalance3D_MP(const InputParameters & parameters)
  :Kernel(parameters),
  _phase_alpha(coupledValue("phase_alpha")),
  _phase_beta(coupledValue("phase_beta")),
  _phase_gamma(coupledValue("phase_gamma")),
  _phase_alpha_var(coupled("phase_alpha")),
  _phase_beta_var(coupled("phase_beta")),
  _phase_gamma_var(coupled("phase_beta")),
   //Global stress
  _stress(getMaterialProperty<RankTwoTensor>("stress")),
  //Its derivatives
  _dstress_dphialpha(getMaterialProperty<RankTwoTensor>("dstress_dphialpha")),
  _dstress_dphibeta(getMaterialProperty<RankTwoTensor>("dstress_dphibeta")),
  _dstress_dphigamma(getMaterialProperty<RankTwoTensor>("dstress_dphigamma")),
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
QMomentumBalance3D_MP::computeQpResidual()
{
    //_stress[_qp].row(_component=0) returns s01, s02, s03
    //The * operator takes a dot product between the gradient of
    // the test function and the stress tensor vector
    return _grad_test[_i][_qp] * _stress[_qp].row(_component);
}

Real
QMomentumBalance3D_MP::computeQpJacobian()
{
  Real jacobian = 0.0;
  // B^T_i * C * B_j
  jacobian += ElasticityTensorTools::elasticJacobian(
      _Jacobian_mult[_qp], _component, _component, _grad_test[_i][_qp], _grad_phi[_j][_qp]);

   return jacobian;
}

Real
QMomentumBalance3D_MP::computeQpOffDiagJacobian(unsigned int jvar)
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

  //if (jvar == _phase_alpha_var){
  //  return _grad_test[_i][_qp] * _dstress_dphialpha[_qp].row(_component) * _phi[_j][_qp];
  //}
  //else if (jvar == _phase_beta_var){
  //  return _grad_test[_i][_qp] * _dstress_dphibeta[_qp].row(_component) * _phi[_j][_qp];
  //}
  //else if (jvar == _phase_gamma_var){
  //  return _grad_test[_i][_qp] * _dstress_dphigamma[_qp].row(_component) * _phi[_j][_qp];
  //}
  //else
  //   return 0;

  return 0;
}
