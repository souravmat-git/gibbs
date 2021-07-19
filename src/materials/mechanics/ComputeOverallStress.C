//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "ComputeOverallStress.h"
#include <fstream>
using namespace std;
registerMooseObject("gibbsApp", ComputeOverallStress);

template <>
InputParameters
validParams<ComputeOverallStress>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<MaterialPropertyName>("strain_jump_name",
                                            "Name of strain jump material");
  return params;
}

ComputeOverallStress::ComputeOverallStress(const InputParameters & parameters)
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
   //Compute the following properties
   _sigma_val(declareProperty<RankTwoTensor>("stress")),
   _dsigma_dphi_val(declareProperty<RankTwoTensor>("dsigma_dphi")),
   _dsigma_de_val(declareProperty<RankFourTensor>("Jacobian_mult"))
{
}

void
ComputeOverallStress::computeQpProperties()
{

    const RankFourTensor avg_C = (_h[_qp]* _beta_stiffness[_qp] +
                              (1.0 - _h[_qp])* _alpha_stiffness[_qp]);

    const RankFourTensor inv_avgC  = (_h[_qp]  * _alpha_stiffness[_qp] +
                                (1.0 - _h[_qp]) * _beta_stiffness[_qp]);

    //overall stress
    _sigma_val[_qp] = _beta_stress[_qp] * _h[_qp]
                    + _alpha_stress[_qp]* (1.0 - _h[_qp]);

    //Declare two rank two tensors to store the term1 and term2
    RankTwoTensor Z1 = (_alpha_stiffness[_qp] - _beta_stiffness[_qp]) * _strain_jump[_qp];
    RankTwoTensor Z2 = inv_avgC * _dstrainjump_dphi[_qp];

    //Derivative wrt phase-field variable (Second rank tensor)
   _dsigma_dphi_val[_qp] =
          _dh[_qp]*(_beta_stress[_qp] - _alpha_stress[_qp])
        + _dh[_qp]*avg_C *_strain_jump[_qp]
        + _h[_qp] * (1.0 - _h[_qp]) * (_alpha_stiffness[_qp] - _beta_stiffness[_qp]) * _dstrainjump_dphi[_qp]
        + _h[_qp] * (1.0 - _h[_qp]) * _ds_de[_qp].transposeMajor() * Z1 
        + _h[_qp] * (1.0 - _h[_qp]) * _ds_de[_qp].transposeMajor() * Z2; 

   //print the value
   // if (_t_step)
   // {
   //  ofstream  file;
   // file.open("dsigma_dphi.csv");
    
   //  for (unsigned int i = 0; i <3; ++i){
   //	for (unsigned int j = 0; j <3; ++j){
   //          file << setprecision(6) << scientific << _dsigma_dphi_val[_qp](i,j) << ",";
   //     }
   //   file << "\n";
   //  }
   //}

   //Derivative wrt strain A_{lmpq} = (avg_C)_lmpq + (diff_C)_lmki*(ds_de)_kipq + term2 + term3
   //where term2_lmpq = h*(1-h)*ds_de(j,i,p,q)*(diff_C)_{j,i,l,m}
   //where term3_lmpq = h*(1-h)*ds_de(j,i,p,q)*(inv_avgC)_{j,i,r,s}*_ds_de(r,s,l,m)

   _dsigma_de_val[_qp] = avg_C + _h[_qp]* (1.0-_h[_qp])*(_alpha_stiffness[_qp] - _beta_stiffness[_qp]) * _ds_de[_qp]
                        + _h[_qp]* (1.0-_h[_qp])*_ds_de[_qp].transposeMajor()*(_alpha_stiffness[_qp] - _beta_stiffness[_qp])
                        + _h[_qp]* (1.0-_h[_qp])*_ds_de[_qp].transposeMajor()* inv_avgC * _ds_de[_qp];

}
