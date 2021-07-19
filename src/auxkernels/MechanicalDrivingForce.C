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

#include "MechanicalDrivingForce.h"

registerMooseObject("gibbsApp", MechanicalDrivingForce);

template<>
InputParameters
validParams<MechanicalDrivingForce>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<MaterialPropertyName>("strain_jump_name",
                                            "Name of strain jump material");
  params.addRequiredParam<MaterialPropertyName>("nd_factor", "Non-dimensional factor"); 
  return params;
}

MechanicalDrivingForce::MechanicalDrivingForce(const InputParameters & parameters)
  :AuxKernel(parameters),
  _h(getMaterialProperty<Real>("h")),
  _dh(getMaterialProperty<Real>("dh")),
  _alpha_strain_energy(getMaterialProperty<Real>("alpha_strain_energy_density")),
  _beta_strain_energy(getMaterialProperty<Real>("beta_strain_energy_density")), 
  _alpha_stress(getMaterialProperty<RankTwoTensor>("alpha_stress")),
  _beta_stress(getMaterialProperty<RankTwoTensor>("beta_stress")),
  //Strain jump is required for this kernel
  _strain_jump(getMaterialProperty<RankTwoTensor>
                      (getParam<MaterialPropertyName>("strain_jump_name"))),
   _nd_factor(getMaterialProperty<Real>
             (getParam<MaterialPropertyName>("nd_factor")))
{
}

Real
MechanicalDrivingForce::driving_force_aux() const{
      
   RankTwoTensor _avg_stress =  (_beta_stress[_qp]* _h[_qp] 
                                    + _alpha_stress[_qp]* (1.0- _h[_qp]) );
        
      
  return  (_beta_strain_energy[_qp] - _alpha_strain_energy[_qp]) 
         + (_avg_stress.doubleContraction(_strain_jump[_qp]));
}
 
Real
MechanicalDrivingForce::computeValue(){
  return _dh[_qp]*driving_force_aux()*_nd_factor[_qp];
}
