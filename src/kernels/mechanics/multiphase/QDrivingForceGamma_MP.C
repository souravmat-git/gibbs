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

#include "QDrivingForceGamma_MP.h"
registerMooseObject("gibbsApp", QDrivingForceGamma_MP);

template<>
InputParameters
validParams<QDrivingForceGamma_MP>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Mechanical driving traction of phase transformation");
  params.addRequiredCoupledVar("phase_beta", "phase-field variable for beta phase");
  params.addRequiredCoupledVar("phase_alpha", "phase-field variable for alpha phase");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
  params.addRequiredParam<MaterialPropertyName>("nd_factor", "Non-dimensional factor");
  return params;
}

QDrivingForceGamma_MP::QDrivingForceGamma_MP(const InputParameters & parameters)
  :Kernel(parameters),
  _phase_beta(coupledValue("phase_beta")),
  _phase_beta_var(coupled("phase_beta")),
  _phase_alpha(coupledValue("phase_alpha")),
  _phase_alpha_var(coupled("phase_alpha")),
  //First derivatives of interpolation function
  _dhalpha_dphigamma(getMaterialProperty<Real>("dhalpha_dphigamma")),
  _dhbeta_dphigamma(getMaterialProperty<Real>("dhbeta_dphigamma")),
  _dhgamma_dphigamma(getMaterialProperty<Real>("dhgamma_dphigamma")),
  //Their second derivatives with resect to gamma
  _d2halpha_dphigamma2(getMaterialProperty<Real>("d2halpha_dphigamma2")),
  _d2hbeta_dphigamma2(getMaterialProperty<Real>("d2halpha_dphigamma2")),
  //Their second derivative with respect to beta
  _d2halpha_dphigamma_dphibeta(getMaterialProperty<Real>("d2halpha_dphigamma_dphibeta")),
  _d2hbeta_dphigamma_dphibeta(getMaterialProperty<Real>("d2hbeta_dphigamma_dphibeta")),
  _d2hgamma_dphigamma_dphibeta(getMaterialProperty<Real>("d2hgamma_dphibeta_dphigamma")),
  //Second derivative with resect to alpha
  _d2halpha_dphigamma_dphialpha(getMaterialProperty<Real>("d2halpha_dphigamma_dphialpha")),
  _d2hbeta_dphigamma_dphialpha(getMaterialProperty<Real>("d2hbeta_dphigamma_dphialpha")),
  _d2hgamma_dphigamma_dphialpha(getMaterialProperty<Real>("d2hgamma_dphialpha_dphigamma")),
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
QDrivingForceGamma_MP::driving_force() const{

  return  _dhalpha_dphigamma[_qp] * (_alpha_strain_energy[_qp]  - _gamma_strain_energy[_qp])
        + _dhbeta_dphigamma[_qp]  * (_beta_strain_energy[_qp]   - _gamma_strain_energy[_qp])
        - _dhalpha_dphigamma[_qp] *  _stress[_qp].doubleContraction(_strain_jump_alpha_beta[_qp])
        + _dhgamma_dphigamma[_qp] *  _stress[_qp].doubleContraction(_strain_jump_beta_gamma[_qp]);
}

Real
QDrivingForceGamma_MP::derivative_df() const{

  Real  _d2hgamma_dphigamma2 = -(_d2hbeta_dphigamma2[_qp] + _d2halpha_dphigamma2[_qp]);

  return  _d2halpha_dphigamma2[_qp] * (_alpha_strain_energy[_qp]  - _gamma_strain_energy[_qp])
        + _d2hbeta_dphigamma2[_qp]  * (_beta_strain_energy[_qp]   - _gamma_strain_energy[_qp])
        - _d2halpha_dphigamma2[_qp] *  _stress[_qp].doubleContraction(_strain_jump_alpha_beta[_qp])
        + _d2hgamma_dphigamma2      *  _stress[_qp].doubleContraction(_strain_jump_beta_gamma[_qp]);
}

Real
QDrivingForceGamma_MP::computeQpResidual(){
  return _test[_i][_qp] *_L[_qp] *  driving_force() * _nd_factor[_qp];
}

Real
QDrivingForceGamma_MP::computeQpJacobian(){
  return _test[_i][_qp]* _L[_qp] * derivative_df()* _nd_factor[_qp]* _phi[_j][_qp];
}

Real
QDrivingForceGamma_MP::computeQpOffDiagJacobian(unsigned int jvar){

  if (jvar == _phase_beta_var)
  {
     Real df_dphibeta =  _d2halpha_dphigamma_dphibeta[_qp] * (_alpha_strain_energy[_qp]  - _gamma_strain_energy[_qp])
                       + _d2hbeta_dphigamma_dphibeta[_qp]  * (_beta_strain_energy[_qp]   - _gamma_strain_energy[_qp])
                       - _d2halpha_dphigamma_dphibeta[_qp] *  _stress[_qp].doubleContraction(_strain_jump_alpha_beta[_qp])
                       + _d2hgamma_dphigamma_dphibeta[_qp] *  _stress[_qp].doubleContraction(_strain_jump_beta_gamma[_qp]);

      return (_test[_i][_qp]* _L[_qp] * df_dphibeta * _nd_factor[_qp]* _phi[_j][_qp]);
   }
   else if  (jvar == _phase_alpha_var)
   {
     Real  df_dphialpha = _d2halpha_dphigamma_dphialpha[_qp] * (_alpha_strain_energy[_qp]  - _gamma_strain_energy[_qp])
                        + _d2hbeta_dphigamma_dphialpha[_qp]  * (_beta_strain_energy[_qp]   - _gamma_strain_energy[_qp])
                        - _d2halpha_dphigamma_dphialpha[_qp] *  _stress[_qp].doubleContraction(_strain_jump_alpha_beta[_qp])
                        + _d2hgamma_dphigamma_dphialpha[_qp] *  _stress[_qp].doubleContraction(_strain_jump_beta_gamma[_qp]);

     return (_test[_i][_qp]* _L[_qp] * df_dphialpha * _nd_factor[_qp]* _phi[_j][_qp]);

   }
   else
    return 0;
 }
