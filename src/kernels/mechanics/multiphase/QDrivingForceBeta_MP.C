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

#include "QDrivingForceBeta_MP.h"
registerMooseObject("gibbsApp", QDrivingForceBeta_MP);

template<>
InputParameters
validParams<QDrivingForceBeta_MP>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Mechanical driving traction of phase transformation");
  //params.addRequiredCoupledVar("ux", "Displacecment in the x-direction");
  params.addRequiredCoupledVar("phase_alpha", "phase-field variable representing alpha phase");
  params.addRequiredCoupledVar("phase_gamma", "phase-field variable representing gamma phase");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
  params.addRequiredParam<MaterialPropertyName>("nd_factor", "Non-dimensional factor");
  return params;
}

QDrivingForceBeta_MP::QDrivingForceBeta_MP(const InputParameters & parameters)
  :Kernel(parameters),
  _phase_alpha(coupledValue("phase_alpha")),
  _phase_alpha_var(coupled("phase_alpha")),
  _phase_gamma(coupledValue("phase_gamma")),
  _phase_gamma_var(coupled("phase_gamma")),
  //First derivatives of interpolation function
  _dhalpha_dphibeta(getMaterialProperty<Real>("dhalpha_dphibeta")),
  _dhgamma_dphibeta(getMaterialProperty<Real>("dhgamma_dphibeta")),
  //Their second derivatives with resect to beta
  _d2halpha_dphibeta2(getMaterialProperty<Real>("d2halpha_dphibeta2")),
  _d2hgamma_dphibeta2(getMaterialProperty<Real>("d2hgamma_dphibeta2")),
  //Their second derivatives with resect to alpha
  _d2halpha_dphibeta_dphialpha(getMaterialProperty<Real>("d2halpha_dphibeta_dphialpha")),
  _d2hgamma_dphibeta_dphialpha(getMaterialProperty<Real>("d2hgamma_dphibeta_dphialpha")),
  //Their second derivatives with resect to gamma
  _d2halpha_dphibeta_dphigamma(getMaterialProperty<Real>("d2halpha_dphibeta_dphigamma")),
  _d2hgamma_dphibeta_dphigamma(getMaterialProperty<Real>("d2hgamma_dphibeta_dphigamma")),
  //Strain energies
  _alpha_strain_energy(getMaterialProperty<Real>("alpha_strain_energy_density")),
  _beta_strain_energy(getMaterialProperty<Real>("beta_strain_energy_density")),
  _gamma_strain_energy(getMaterialProperty<Real>("gamma_strain_energy_density")),
  //overall stress
  _stress(getMaterialProperty<RankTwoTensor>("stress")),
  //Strain jumps
  _strain_jump_alpha_beta(getMaterialProperty<RankTwoTensor>("strain_jump_alpha_beta")),
  _strain_jump_beta_gamma(getMaterialProperty<RankTwoTensor>("strain_jump_beta_gamma")),
  //mobility and nd_factor
  _L(getMaterialProperty<Real>(getParam<MaterialPropertyName>("mob_name"))),
  _nd_factor(getMaterialProperty<Real>
             (getParam<MaterialPropertyName>("nd_factor")))
{
}

Real
QDrivingForceBeta_MP::driving_force() const{

  return  _dhalpha_dphibeta[_qp] * (_alpha_strain_energy[_qp]  - _beta_strain_energy[_qp])
        + _dhgamma_dphibeta[_qp] * (_gamma_strain_energy[_qp]  - _beta_strain_energy[_qp])
        - _dhalpha_dphibeta[_qp] *  _stress[_qp].doubleContraction(_strain_jump_alpha_beta[_qp])
        + _dhgamma_dphibeta[_qp] *  _stress[_qp].doubleContraction(_strain_jump_beta_gamma[_qp]);
}

Real
QDrivingForceBeta_MP::derivative_df() const{

  return  _d2halpha_dphibeta2[_qp] * (_alpha_strain_energy[_qp]  - _beta_strain_energy[_qp])
        + _d2hgamma_dphibeta2[_qp] * (_gamma_strain_energy[_qp]  - _beta_strain_energy[_qp])
        - _d2halpha_dphibeta2[_qp] *  _stress[_qp].doubleContraction(_strain_jump_alpha_beta[_qp])
        - _d2hgamma_dphibeta2[_qp] *  _stress[_qp].doubleContraction(_strain_jump_beta_gamma[_qp]);
}

Real
QDrivingForceBeta_MP::computeQpResidual(){
  return _test[_i][_qp] * _L[_qp] *  driving_force() * _nd_factor[_qp];
}

Real
QDrivingForceBeta_MP::computeQpJacobian(){
  return _test[_i][_qp] * _L[_qp] * derivative_df() * _nd_factor[_qp] * _phi[_j][_qp];
}

Real
QDrivingForceBeta_MP::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _phase_alpha_var)
  {
     Real df_dphialpha = _d2halpha_dphibeta_dphialpha[_qp] * (_alpha_strain_energy[_qp]  - _beta_strain_energy[_qp])
                       + _d2hgamma_dphibeta_dphialpha[_qp] * (_gamma_strain_energy[_qp]  - _beta_strain_energy[_qp])
                       - _d2halpha_dphibeta_dphialpha[_qp] *  _stress[_qp].doubleContraction(_strain_jump_alpha_beta[_qp])
                       + _d2hgamma_dphibeta_dphialpha[_qp] *  _stress[_qp].doubleContraction(_strain_jump_beta_gamma[_qp]);

      return _test[_i][_qp] * _L[_qp] * df_dphialpha * _nd_factor[_qp] * _phi[_j][_qp];
  }
  else if (jvar == _phase_gamma_var)
  {
    Real df_dphigamma = _d2halpha_dphibeta_dphigamma[_qp] * (_alpha_strain_energy[_qp]  - _beta_strain_energy[_qp])
                      + _d2hgamma_dphibeta_dphigamma[_qp] * (_gamma_strain_energy[_qp]  - _beta_strain_energy[_qp])
                      - _d2halpha_dphibeta_dphigamma[_qp] *  _stress[_qp].doubleContraction(_strain_jump_alpha_beta[_qp])
                      + _d2hgamma_dphibeta_dphigamma[_qp] *  _stress[_qp].doubleContraction(_strain_jump_beta_gamma[_qp]);

      return _test[_i][_qp] * _L[_qp] * df_dphigamma * _nd_factor[_qp] * _phi[_j][_qp];

  }

    return 0;
}
