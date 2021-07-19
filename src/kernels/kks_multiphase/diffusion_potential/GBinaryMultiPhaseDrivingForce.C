//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*

//*This Kernel implements the bulk part of the
//*Allen Cahn equation which tells us that at
//*equilibrium the grandpotentials must be equal
//*In this code both the free energy is obtained 
//* from the material class AnalyticalFreeEnergyMaterial

#include "GBinaryMultiPhaseDrivingForce.h"

registerMooseObject("gibbsApp", GBinaryMultiPhaseDrivingForce);

template <>
InputParameters
validParams<GBinaryMultiPhaseDrivingForce>()
{
  InputParameters params = ACBulk<Real>::validParams();
  params.addClassDescription("Eqn:Eqn: dh*(mu_A^{\beta} - mu_A^{\alpha})"
                             "This kernel operates on eta.");
  params.addRequiredCoupledVar("phase_2", "Phase_2, Note: (Phase2 - Phase1)");
  params.addCoupledVar("phase_3", 0.0, "Phase_3, Note: (Phase3 - Phase1)");
  params.addCoupledVar("phase_4", 0.0, "phase_4");
  params.addCoupledVar("phase_5", 0.0, "Phase_5");
  params.addRequiredCoupledVar("B_diff_pot", "Diffusion potential of component B");
  params.addRequiredParam<MaterialPropertyName>("A_chem_pot_1", "Grand-potential phase 1");
  params.addRequiredParam<MaterialPropertyName>("A_chem_pot_2", "Grand-potential phase 2");
  params.addRequiredParam<MaterialPropertyName>("xB_1", "Mole fraction of B in phase 1");
  params.addRequiredParam<MaterialPropertyName>("xB_2", "Mole fraction of B in phase 2");
  params.addRequiredParam<MaterialPropertyName>("dh", "Derivative of phase2 wr.t 1");
  params.addRequiredParam<MaterialPropertyName>("d2h","Derivative of dh w.r.t 1");
  params.addRequiredParam<MaterialPropertyName>("d2h_2", "Derivative of dh w.r.t 2");
  params.addParam<MaterialPropertyName>("d2h_3", 0.0,"Derivative of dh w.r.t 3");
  params.addParam<MaterialPropertyName>("d2h_4", 0.0,"Derivative of dh w.r.t 4");
  params.addParam<MaterialPropertyName>("d2h_5", 0.0, "Deriavtive of dh w.r.t 5");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
  params.addParam<MaterialPropertyName>("nd_factor", 1.0, "RT/Vm*barrier_height");
  return params;
}

GBinaryMultiPhaseDrivingForce::GBinaryMultiPhaseDrivingForce(const InputParameters & parameters)
  : ACBulk<Real>(parameters),
   _phase_2(coupledValue("phase_2")),
   _phase_2_var(coupled("phase_2")),
   //Extend to multiple phases
   _phase_3(coupledValue("phase_3")),
   _phase_3_var(coupled("phase_3")),
   //Multiple phases
   _phase_4(coupledValue("phase_4")),
   _phase_4_var(coupled("phase_4")),
   _phase_5(coupledValue("phase_5")),
   _phase_5_var(coupled("phase_5")),
   //for an A-B alloy
   _B_diff_pot(coupledValue("B_diff_pot")),
   _B_diff_pot_var(coupled("B_diff_pot")),
   //Get all the declared names for material property
   _A_chem_pot_1_name(getParam<MaterialPropertyName>("A_chem_pot_1")),
   _A_chem_pot_2_name(getParam<MaterialPropertyName>("A_chem_pot_2")),
   _xB_1_name(getParam<MaterialPropertyName>("xB_1")),
   _xB_2_name(getParam<MaterialPropertyName>("xB_2")),
   //Get all the declared names for interpolation function
   _dh_name(getParam<MaterialPropertyName>("dh")),
   _d2h_name(getParam<MaterialPropertyName>("d2h")),
   _d2h_2_name(getParam<MaterialPropertyName>("d2h_2")),
   _d2h_3_name(getParam<MaterialPropertyName>("d2h_3")),
   _d2h_4_name(getParam<MaterialPropertyName>("d2h_4")),
   _d2h_5_name(getParam<MaterialPropertyName>("d2h_5")),
   //Material specific properties
   _A_chem_pot_1(getMaterialProperty<Real>(_A_chem_pot_1_name)),
   _A_chem_pot_2(getMaterialProperty<Real>(_A_chem_pot_2_name)), 
   _xB_1(getMaterialProperty<Real>(_xB_1_name)),
   _xB_2(getMaterialProperty<Real>(_xB_2_name)),
   //interpolation function
   _dh(getMaterialProperty<Real>(_dh_name)),
   _d2h(getMaterialProperty<Real>(_d2h_name)),
   _d2h_2(getMaterialProperty<Real>(_d2h_2_name)),
   _d2h_3(getMaterialProperty<Real>(_d2h_3_name)),
   _d2h_4(getMaterialProperty<Real>(_d2h_4_name)),
   _d2h_5(getMaterialProperty<Real>(_d2h_5_name)),
   //Phase-field mobility is now obtained from ACBulk
   _nd_factor(getMaterialProperty<Real>("nd_factor"))
{
} 
        
Real
GBinaryMultiPhaseDrivingForce::computeDFDOP(PFFunctionType type)
{
  switch (type)
  {
    case Residual:
    {
      //_dF_dop = (RT/Vm*m)*h'*(Difference of GPs)
      return (_nd_factor[_qp]* _dh[_qp] * (_A_chem_pot_2[_qp] -_A_chem_pot_1[_qp]));
    }
    case Jacobian:
    {
      //_JdF_dop = (RT/Vm*m)*h'*(Difference of GPs)
      return (_nd_factor[_qp] * _d2h[_qp] * (_A_chem_pot_2[_qp] - _A_chem_pot_1[_qp])* _phi[_j][_qp]); 
    }
  }  
  mooseError("Internal error"); 
}   

Real
GBinaryMultiPhaseDrivingForce::computeQpOffDiagJacobian(unsigned int jvar)
{
    
  if (jvar == _B_diff_pot_var)
  {
    return ACBulk<Real>::computeQpOffDiagJacobian(jvar) +  
      (_L[_qp]*(_test[_i][_qp] * _nd_factor[_qp] * _dh[_qp] * (_xB_1[_qp] - _xB_2[_qp])* _phi[_j][_qp])); 
  }
  else if (jvar == _phase_2_var)
  {
    return  ACBulk<Real>::computeQpOffDiagJacobian(jvar) + 
      (_L[_qp]*(_test[_i][_qp] * _nd_factor[_qp] * _d2h_2[_qp]* (_A_chem_pot_2[_qp] 
                                                   - _A_chem_pot_1[_qp])* _phi[_j][_qp]));
  }
  else if (jvar == _phase_3_var)
  {
    return  ACBulk<Real>::computeQpOffDiagJacobian(jvar) + 
      (_L[_qp]*(_test[_i][_qp] * _nd_factor[_qp] * _d2h_3[_qp]* (_A_chem_pot_2[_qp] 
                                                   - _A_chem_pot_1[_qp])* _phi[_j][_qp]));
  }
  else if (jvar == _phase_4_var)
  {
    return  ACBulk<Real>::computeQpOffDiagJacobian(jvar) + 
      (_L[_qp]*(_test[_i][_qp] * _nd_factor[_qp] * _d2h_4[_qp]* (_A_chem_pot_2[_qp] 
                                                   - _A_chem_pot_1[_qp])* _phi[_j][_qp]));
  }
  else if (jvar == _phase_5_var)
  {
    return  ACBulk<Real>::computeQpOffDiagJacobian(jvar) + 
      (_L[_qp]*(_test[_i][_qp] * _nd_factor[_qp] * _d2h_5[_qp]* (_A_chem_pot_2[_qp] 
                                                   - _A_chem_pot_1[_qp])* _phi[_j][_qp]));
  }
  else    
      return 0.0; 
}
