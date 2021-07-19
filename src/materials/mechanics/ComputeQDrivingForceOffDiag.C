//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "ComputeQDrivingForceOffDiag.h"
registerMooseObject("gibbsApp", ComputeQDrivingForceOffDiag);

template <>
InputParameters
validParams<ComputeQDrivingForceOffDiag>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<MaterialPropertyName>("strain_jump_name",
                                            "Name of strain jump material");
  params.addRequiredParam<MaterialPropertyName>("nd_factor", "Non-dimensional factor");
  return params;
}

ComputeQDrivingForceOffDiag::ComputeQDrivingForceOffDiag(const InputParameters & parameters)
  : Material(parameters),
     //Stresses of alpha and beta phase
   _alpha_stress(getMaterialProperty<RankTwoTensor>("alpha_stress")),
   _beta_stress(getMaterialProperty<RankTwoTensor>("beta_stress")),
   //Stiffness of alpha and beta phase
   _alpha_stiffness(getMaterialProperty<RankFourTensor>("alpha_elasticity_tensor")),
   _beta_stiffness(getMaterialProperty<RankFourTensor>("beta_elasticity_tensor")),
   //Strain jump
   _strain_jump(getMaterialProperty<RankTwoTensor>
                      (getParam<MaterialPropertyName>("strain_jump_name"))),
   //and its derivatives wrt to strain and phi
   _dstrainjump_dphi(getMaterialProperty<RankTwoTensor>("dstrainjump_dphi")),
   _ds_de(getMaterialProperty<RankFourTensor>("ds_de")),
   //Interpolation function
   _h(getMaterialProperty<Real>("h")),  
   _dh(getMaterialProperty<Real>("dh")),
   //A non-dimensional factor which is equal to RT/Vm*(barrier_height)
   _nd_factor(getMaterialProperty<Real>
             (getParam<MaterialPropertyName>("nd_factor"))),
   //Compute the following properties
   _d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain"))
{
} 

void
ComputeQDrivingForceOffDiag::computeQpProperties()
{
     RankTwoTensor _store = (1.0 - 2.0 * _h[_qp])*
                                    (_alpha_stress[_qp] - _beta_stress[_qp]);
                                    
     RankFourTensor _M = (_alpha_stiffness[_qp] - _beta_stiffness[_qp])*_ds_de[_qp];                                  
            
    //Driving force with respect to strain
    _d2Fdcdstrain[_qp] = _nd_factor[_qp] * _dh[_qp] * ((_beta_stress[_qp] - _alpha_stress[_qp])
                                                      + _ds_de[_qp].innerProductTranspose(_store)
                        + (_h[_qp] * _beta_stiffness[_qp] + (1-_h[_qp]) * _alpha_stiffness[_qp]) * _strain_jump[_qp]
                        + _M.innerProductTranspose(_strain_jump[_qp]) * _h[_qp] * (1.0 - _h[_qp]));        
           
}
