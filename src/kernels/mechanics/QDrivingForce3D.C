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

#include "QDrivingForce3D.h"
registerMooseObject("gibbsApp", QDrivingForce3D);

template<>
InputParameters
validParams<QDrivingForce3D>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Mechanical driving traction of phase transformation");
  //params.addRequiredCoupledVar("ux", "Displacecment in the x-direction");
  params.addRequiredParam<MaterialPropertyName>("strain_jump_name",
                                            "Name of strain jump material");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
  params.addRequiredParam<MaterialPropertyName>("nd_factor", "Non-dimensional factor");
  return params;
}

QDrivingForce3D::QDrivingForce3D(const InputParameters & parameters)
  :Kernel(parameters),
  //Interpolation function
  _h(getMaterialProperty<Real>("h")),
  _dh(getMaterialProperty<Real>("dh")),
  _d2h(getMaterialProperty<Real>("d2h")),
  //Requisite material properties
  _alpha_strain_energy(getMaterialProperty<Real>("alpha_strain_energy_density")),
  _beta_strain_energy(getMaterialProperty<Real>("beta_strain_energy_density")),
  _alpha_stress(getMaterialProperty<RankTwoTensor>("alpha_stress")),
  _beta_stress(getMaterialProperty<RankTwoTensor>("beta_stress")),
  //Strain jump is required for this kernel
  _strain_jump(getMaterialProperty<RankTwoTensor>
                      (getParam<MaterialPropertyName>("strain_jump_name"))),
  //modulus values
  _alpha_stiffness(getMaterialProperty<RankFourTensor>("alpha_elasticity_tensor")),
  _beta_stiffness(getMaterialProperty<RankFourTensor>("beta_elasticity_tensor")),
  //mobility and nd_factor
  _L(getMaterialProperty<Real>
             (getParam<MaterialPropertyName>("mob_name"))),
  _nd_factor(getMaterialProperty<Real>
             (getParam<MaterialPropertyName>("nd_factor")))
{
}

Real
QDrivingForce3D::driving_force() const{

  const RankTwoTensor _avg_stress =  (_beta_stress[_qp]* _h[_qp]
                                    + _alpha_stress[_qp]* (1.0- _h[_qp]) );


  return  (_beta_strain_energy[_qp] - _alpha_strain_energy[_qp])
         + (_avg_stress.doubleContraction(_strain_jump[_qp]));
}

Real
QDrivingForce3D::computeQpResidual(){
  return _test[_i][_qp] *_L[_qp] * _dh[_qp] * driving_force() * _nd_factor[_qp];
}

Real
QDrivingForce3D::computeQpJacobian(){
  return _test[_i][_qp]* _L[_qp] * _d2h[_qp] * driving_force()
                                  * _nd_factor[_qp]* _phi[_j][_qp];
}

Real
QDrivingForce3D::computeQpOffDiagJacobian(unsigned int /*jvar*/){
    return 0;
}
