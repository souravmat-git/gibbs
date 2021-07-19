//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
// Forward Declarations
class DrivingForceGradPhiTheta;

//MOOSe includes
#include "VariableUnitNormalBase.h"

template <>
InputParameters validParams<DrivingForceGradPhiTheta>();

class DrivingForceGradPhiTheta : public VariableUnitNormalBase
{
public:
  DrivingForceGradPhiTheta(const InputParameters & parameters);

  //This class calculates the bulk driving force due to gradient in phase-field

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

   const MaterialProperty<RankTwoTensor> & _alpha_stress;
   const MaterialProperty<RankTwoTensor> & _beta_stress;
   const MaterialProperty<RankTwoTensor> & _gamma_stress;

   //Interpolation functions
   const MaterialProperty<Real> & _h_alpha;
   const MaterialProperty<Real> & _h_beta;
   const MaterialProperty<Real> & _h_gamma;

   //Jump vectors
   const MaterialProperty<RealVectorValue> & _avec_alpha_beta;
   const MaterialProperty<RealVectorValue> & _avec_beta_gamma;

   //bulk driving force for alpha phase
   MaterialProperty<RealVectorValue>  & _alpha_dfbulk_dgradphi;
   MaterialProperty<RealVectorValue>  & _alpha_d2fbulk_dgradphi_dphi;

   //bulk driving force for beta phase
   MaterialProperty<RealVectorValue>  & _beta_dfbulk_dgradphi;
   MaterialProperty<RealVectorValue>  & _beta_d2fbulk_dgradphi_dphi;

   //bulk driving force for gamma phase
   MaterialProperty<RealVectorValue>  & _gamma_dfbulk_dgradphi;
   MaterialProperty<RealVectorValue>  & _gamma_d2fbulk_dgradphi_dphi;

};
