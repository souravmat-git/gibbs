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

#include "TernaryMultiPhaseBase.h"

registerMooseObject("gibbsApp", TernaryMultiPhaseBase);

template <>
InputParameters
validParams<TernaryMultiPhaseBase>()
{
  InputParameters params = validParams<MultiPhaseBase>();
  params.addClassDescription("Base class for continuity eqn for a A-B-C alloy");
  //B thermodynamic factor
  params.addParam<MaterialPropertyName>("inv_B_tf_gamma",  0.0, "Thermodynamic factor B in gamma phase");
  params.addParam<MaterialPropertyName>("inv_B_tf_delta",  0.0, "Thermodynamic factor B in delta phase");
  params.addParam<MaterialPropertyName>("inv_B_tf_epsilon",0.0, "Thermodynamic factor B in epsilon phase");
  //BC thermodynamic factor
  params.addParam<MaterialPropertyName>("inv_BC_tf_gamma", 0.0, "Thermodynamic factor BC in gamma phase");
  params.addParam<MaterialPropertyName>("inv_BC_tf_delta", 0.0, "Thermodynamic factor BC in delta phase");
  params.addParam<MaterialPropertyName>("inv_BC_tf_epsilon", 0.0, "Thermodynamic factor BC in epsilon phase");
  //C thermodynamic factor
  params.addParam<MaterialPropertyName>("inv_C_tf_gamma", 0.0, "Thermodynamic factor C in gamma phase");
  params.addParam<MaterialPropertyName>("inv_C_tf_delta", 0.0, "Thermodynamic factor C in delta phase");
  params.addParam<MaterialPropertyName>("inv_C_tf_epsilon", 0.0, "Thermodynamic factor C in epsilon phase");                              
  //***************************************Onsager mobility assumed only for two phases***************//
  //L_BB
  params.addParam<MaterialPropertyName>("L_BB_gamma", 0.0, "Onsager mobility of comp B in gamma phase");
  params.addParam<MaterialPropertyName>("L_BB_delta", 0.0, "Onsager mobility of comp B in delta phase");
  params.addParam<MaterialPropertyName>("L_BB_epsilon", 0.0, "TOnsager mobility of comp B in epsilon phase");  
  //L_BC
  params.addParam<MaterialPropertyName>("L_BC_gamma", 0.0, "Onsager mobility of comp BC in gamma phase");
  params.addParam<MaterialPropertyName>("L_BC_delta", 0.0, "Onsager mobility of comp BC in delta phase");
  params.addParam<MaterialPropertyName>("L_BC_epsilon", 0.0, "Onsager mobility of comp BC in epsilon phase");  
  //L_CC
  params.addParam<MaterialPropertyName>("L_CC_gamma", 0.0, "Onsager mobility of comp C in gamma phase");
  params.addParam<MaterialPropertyName>("L_CC_delta", 0.0, "Onsager mobility of comp C in delta phase");
  params.addParam<MaterialPropertyName>("L_CC_epsilon", 0.0, "Onsager mobility of comp C in epsilon phase"); 
  //*************************Derivatives of L wr.t B********************************************************// 
  //dL_BB_muB
  params.addParam<MaterialPropertyName>("dL_BB_muB_gamma", 0.0, "Onsager mobility of comp B in gamma phase");
  params.addParam<MaterialPropertyName>("dL_BB_muB_delta", 0.0, "Onsager mobility of comp B in delta phase");
  params.addParam<MaterialPropertyName>("dL_BB_muB_epsilon", 0.0, "TOnsager mobility of comp B in epsilon phase");
  //dL_BC_muB
  params.addParam<MaterialPropertyName>("dL_BC_muB_gamma", 0.0, "Onsager mobility of comp BC in gamma phase");
  params.addParam<MaterialPropertyName>("dL_BC_muB_delta", 0.0, "Onsager mobility of comp BC in delta phase");
  params.addParam<MaterialPropertyName>("dL_BC_muB_epsilon", 0.0, "Onsager mobility of comp BC in epsilon phase"); 
  //dL_CC_muB
  params.addParam<MaterialPropertyName>("dL_CC_muB_gamma", 0.0, "Onsager mobility of comp C in gamma phase");
  params.addParam<MaterialPropertyName>("dL_CC_muB_delta", 0.0, "Onsager mobility of comp C in delta phase");
  params.addParam<MaterialPropertyName>("dL_CC_muB_epsilon", 0.0, "Onsager mobility of comp C in epsilon phase");  
  //*************************Derivatives of L w.r.t muC********************************************************//  
  //dL_BB_muC
  params.addParam<MaterialPropertyName>("dL_BB_muC_gamma", 0.0, "Onsager mobility of comp B in gamma phase");
  params.addParam<MaterialPropertyName>("dL_BB_muC_delta", 0.0, "Onsager mobility of comp B in delta phase");
  params.addParam<MaterialPropertyName>("dL_BB_muC_epsilon", 0.0, "TOnsager mobility of comp B in epsilon phase");  
  //dL_BC_muC
  params.addParam<MaterialPropertyName>("dL_BC_muC_gamma", 0.0, "Onsager mobility of comp BC in gamma phase");
  params.addParam<MaterialPropertyName>("dL_BC_muC_delta", 0.0, "Onsager mobility of comp BC in delta phase");
  params.addParam<MaterialPropertyName>("dL_BC_muC_epsilon", 0.0, "Onsager mobility of comp BC in epsilon phase");
  //dL_CC_muC
  params.addParam<MaterialPropertyName>("dL_CC_muC_gamma", 0.0, "Onsager mobility of comp C in gamma phase");
  params.addParam<MaterialPropertyName>("dL_CC_muC_delta", 0.0, "Onsager mobility of comp C in delta phase");
  params.addParam<MaterialPropertyName>("dL_CC_muC_epsilon", 0.0, "Onsager mobility of comp C in epsilon phase");
  return params; 
}

TernaryMultiPhaseBase::TernaryMultiPhaseBase(const InputParameters & parameters)
  : MultiPhaseBase(parameters),
  //coeff B of the susceptibility matrix
  _inv_B_tf_alpha(getMaterialProperty<Real>("inv_B_tf_alpha")),
  _inv_B_tf_beta(getMaterialProperty<Real>("inv_B_tf_beta")),
  _inv_B_tf_gamma(getMaterialProperty<Real>("inv_B_tf_gamma")),
  _inv_B_tf_delta(getMaterialProperty<Real>("inv_B_tf_delta")),
  _inv_B_tf_epsilon(getMaterialProperty<Real>("inv_B_tf_epsilon")),
  //coeff BC of the susceptibility matrix
  _inv_BC_tf_alpha(getMaterialProperty<Real>("inv_BC_tf_alpha")),
  _inv_BC_tf_beta(getMaterialProperty<Real>("inv_BC_tf_beta")),
  _inv_BC_tf_gamma(getMaterialProperty<Real>("inv_BC_tf_gamma")),
  _inv_BC_tf_delta(getMaterialProperty<Real>("inv_BC_tf_delta")),
  _inv_BC_tf_epsilon(getMaterialProperty<Real>("inv_BC_tf_epsilon")),
  //coeff C of the susceptibility matrix
  _inv_C_tf_alpha(getMaterialProperty<Real>("inv_C_tf_alpha")),
  _inv_C_tf_beta(getMaterialProperty<Real>("inv_C_tf_beta")),
  _inv_C_tf_gamma(getMaterialProperty<Real>("inv_C_tf_gamma")),
  _inv_C_tf_delta(getMaterialProperty<Real>("inv_C_tf_delta")),
  _inv_C_tf_epsilon(getMaterialProperty<Real>("inv_C_tf_epsilon")),
  //Onsager matrix
  _L_BB_alpha(getMaterialProperty<Real>("L_BB_alpha")),
  _L_BB_beta(getMaterialProperty<Real>("L_BB_beta")),
  //assumed to be zero..see above
  _L_BB_gamma(getMaterialProperty<Real>("L_BB_gamma")),
  _L_BB_delta(getMaterialProperty<Real>("L_BB_delta")),
  _L_BB_epsilon(getMaterialProperty<Real>("L_BB_epsilon")),
    //L_BC within each phase
  _L_BC_alpha(getMaterialProperty<Real>("L_BC_alpha")),
  _L_BC_beta(getMaterialProperty<Real>("L_BC_beta")),
  _L_BC_gamma(getMaterialProperty<Real>("L_BC_gamma")),
  _L_BC_delta(getMaterialProperty<Real>("L_BC_delta")),
  _L_BC_epsilon(getMaterialProperty<Real>("L_BC_epsilon")),
  //L_CC within each phase
  _L_CC_alpha(getMaterialProperty<Real>("L_CC_alpha")),
  _L_CC_beta(getMaterialProperty<Real>("L_CC_beta")),
  _L_CC_gamma(getMaterialProperty<Real>("L_CC_gamma")),
  _L_CC_delta(getMaterialProperty<Real>("L_CC_delta")),
  _L_CC_epsilon(getMaterialProperty<Real>("L_CC_epsilon")),
  //dL_BB_muB
  _dL_BB_muB_alpha(getMaterialProperty<Real>("dL_BB_muB_alpha")),
  _dL_BB_muB_beta(getMaterialProperty<Real>("dL_BB_muB_beta")),
  _dL_BB_muB_gamma(getMaterialProperty<Real>("dL_BB_muB_gamma")),
  _dL_BB_muB_delta(getMaterialProperty<Real>("dL_BB_muB_delta")),
  _dL_BB_muB_epsilon(getMaterialProperty<Real>("dL_BB_muB_epsilon")),
  //dL_BC_muB
  _dL_BC_muB_alpha(getMaterialProperty<Real>("dL_BC_muB_alpha")),
  _dL_BC_muB_beta(getMaterialProperty<Real>("dL_BC_muB_beta")),
  _dL_BC_muB_gamma(getMaterialProperty<Real>("dL_BC_muB_gamma")),
  _dL_BC_muB_delta(getMaterialProperty<Real>("dL_BC_muB_delta")),
  _dL_BC_muB_epsilon(getMaterialProperty<Real>("dL_BC_muB_epsilon")),
  //dL_CC_muB
  _dL_CC_muB_alpha(getMaterialProperty<Real>("dL_CC_muB_alpha")),
  _dL_CC_muB_beta(getMaterialProperty<Real>("dL_CC_muB_beta")),
  _dL_CC_muB_gamma(getMaterialProperty<Real>("dL_CC_muB_gamma")),
  _dL_CC_muB_delta(getMaterialProperty<Real>("dL_CC_muB_delta")),
  _dL_CC_muB_epsilon(getMaterialProperty<Real>("dL_CC_muB_epsilon")),
   //dL_BB_muC
  _dL_BB_muC_alpha(getMaterialProperty<Real>("dL_BB_muC_alpha")),
  _dL_BB_muC_beta(getMaterialProperty<Real>("dL_BB_muC_beta")),
  _dL_BB_muC_gamma(getMaterialProperty<Real>("dL_BB_muC_gamma")),
  _dL_BB_muC_delta(getMaterialProperty<Real>("dL_BB_muC_delta")),
  _dL_BB_muC_epsilon(getMaterialProperty<Real>("dL_BB_muC_epsilon")),
  //dL_BC_muC
  _dL_BC_muC_alpha(getMaterialProperty<Real>("dL_BC_muC_alpha")),
  _dL_BC_muC_beta(getMaterialProperty<Real>("dL_BC_muC_beta")),
  _dL_BC_muC_gamma(getMaterialProperty<Real>("dL_BC_muC_gamma")),
  _dL_BC_muC_delta(getMaterialProperty<Real>("dL_BC_muC_delta")),
  _dL_BC_muC_epsilon(getMaterialProperty<Real>("dL_BC_muC_epsilon")),
  //dL_CC_muC
  _dL_CC_muC_alpha(getMaterialProperty<Real>("dL_CC_muC_alpha")),
  _dL_CC_muC_beta(getMaterialProperty<Real>("dL_CC_muC_beta")),
  _dL_CC_muC_gamma(getMaterialProperty<Real>("dL_CC_muC_gamma")),
  _dL_CC_muC_delta(getMaterialProperty<Real>("dL_CC_muC_delta")),
  _dL_CC_muC_epsilon(getMaterialProperty<Real>("dL_CC_muC_epsilon"))
{
}

//****************************************************************************//
//The following steps is followed here:
//1) Create the coefficents of the overall susceptibility matrix
//2) Generate the determinant of the susceptibility matrix
//3) Create the coefficients of the overall thermodynamic factor matrix
//4) Interpolate the coeff. of the Onsager mobilities
//5) Perform a matrix product to obtain the Overall diffusioon coeff. 
//****************************************************************************//

Real
TernaryMultiPhaseBase::chi_BB() const
{
   return (_h_alpha[_qp]  * _inv_B_tf_alpha[_qp] 
         + _h_beta[_qp]   * _inv_B_tf_beta[_qp]
         + _h_gamma[_qp]  * _inv_B_tf_gamma[_qp]
         + _h_delta[_qp]  * _inv_B_tf_delta[_qp]
         + _h_epsilon[_qp]* _inv_B_tf_epsilon[_qp]);
 
}

Real
TernaryMultiPhaseBase::chi_BC() const
{
   return (_h_alpha[_qp]  * _inv_BC_tf_alpha[_qp] 
         + _h_beta[_qp]   * _inv_BC_tf_beta[_qp]
         + _h_gamma[_qp]  * _inv_BC_tf_gamma[_qp]
         + _h_delta[_qp]  * _inv_BC_tf_delta[_qp]
         + _h_epsilon[_qp]* _inv_BC_tf_epsilon[_qp]);
 
}

Real
TernaryMultiPhaseBase::chi_CC() const
{
   return (_h_alpha[_qp]  * _inv_C_tf_alpha[_qp] 
         + _h_beta[_qp]   * _inv_C_tf_beta[_qp]
         + _h_gamma[_qp]  * _inv_C_tf_gamma[_qp]
         + _h_delta[_qp]  * _inv_C_tf_delta[_qp]
         + _h_epsilon[_qp]* _inv_C_tf_epsilon[_qp]);
 
}

Real
TernaryMultiPhaseBase::det_chi() const
{
  return (TernaryMultiPhaseBase::chi_BB()* TernaryMultiPhaseBase::chi_CC()
         -TernaryMultiPhaseBase::chi_BC()* TernaryMultiPhaseBase::chi_BC());
}


//All coeffecients of the overall thermodynamic factor matrix

Real
TernaryMultiPhaseBase::thermodynamic_factorBB() const
{
   //The coeff BB of the overall TF matrix is:
   //1/det(overall_chi)*chi_CC;
   return ( (std::pow(TernaryMultiPhaseBase::det_chi(),-1)) 
                     *TernaryMultiPhaseBase::chi_CC() ); 
}

Real
TernaryMultiPhaseBase::thermodynamic_factorBC() const
{
   //The coeff BB of the overall TF matrix is:
   //1/det(overall_chi)*chi_BC;
   return -( (1.0/TernaryMultiPhaseBase::det_chi()) 
                 *TernaryMultiPhaseBase::chi_BC() );
}

Real
TernaryMultiPhaseBase::thermodynamic_factorCC() const
{                     
   //The coeff CC of the overall TF matrix is:
   //1/det(overall_chi)*chi_CC;
   return ( (1.0/TernaryMultiPhaseBase::det_chi()) 
                *TernaryMultiPhaseBase::chi_BB() ); 
}

//Interpolated Onsager Monilities
 
Real
TernaryMultiPhaseBase::L_BB_interp() const
{
  return (_L_BB_alpha[_qp] * _h_alpha[_qp] 
        + _L_BB_beta[_qp]  * _h_beta[_qp] 
        + _L_BB_gamma[_qp] * _h_gamma[_qp] 
        + _L_BB_delta[_qp] * _h_delta[_qp]
        + _L_BB_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
TernaryMultiPhaseBase::L_BC_interp() const
{
  return (_L_BC_alpha[_qp] * _h_alpha[_qp] 
        + _L_BC_beta[_qp]  * _h_beta[_qp] 
        + _L_BC_gamma[_qp] * _h_gamma[_qp] 
        + _L_BC_delta[_qp] * _h_delta[_qp]
        + _L_BC_epsilon[_qp]* _h_epsilon[_qp]);
}


Real
TernaryMultiPhaseBase::L_CC_interp() const
{
  return (_L_CC_alpha[_qp]  * _h_alpha[_qp] 
        + _L_CC_beta[_qp]   * _h_beta[_qp] 
        + _L_CC_gamma[_qp]  * _h_gamma[_qp] 
        + _L_CC_delta[_qp]  * _h_delta[_qp]
        + _L_CC_epsilon[_qp]* _h_epsilon[_qp]);
}


//From the Onsager mobility and the thermodynamic factor, we can derive the 
//overall interdiffusion coefficents, note that these need not be symmetric
//First is for the flux of component of B
Real
TernaryMultiPhaseBase::DC_BB_interp() const
{ 
  return (TernaryMultiPhaseBase::L_BB_interp()*TernaryMultiPhaseBase::thermodynamic_factorBB()
         +TernaryMultiPhaseBase::L_BC_interp()*TernaryMultiPhaseBase::thermodynamic_factorBC());
}

Real
TernaryMultiPhaseBase::DC_BC_interp() const
{
  //Note thermodynamic factor CD = thermodynamic factor DC
  return (TernaryMultiPhaseBase::L_BB_interp()*TernaryMultiPhaseBase::thermodynamic_factorBC()
         +TernaryMultiPhaseBase::L_BC_interp()*TernaryMultiPhaseBase::thermodynamic_factorCC());
}

//First is for the flux of component of C
Real
TernaryMultiPhaseBase::DC_CB_interp() const
{
  //Note Onsager Mobility BC = Onsager Mobility BC
  return (TernaryMultiPhaseBase::L_BC_interp()*TernaryMultiPhaseBase::thermodynamic_factorBB()
         +TernaryMultiPhaseBase::L_CC_interp()*TernaryMultiPhaseBase::thermodynamic_factorBC());
}

Real
TernaryMultiPhaseBase::DC_CC_interp() const
{
  //Note Onsager Mobility BC = Onsager Mobility BC
  return (TernaryMultiPhaseBase::L_BC_interp()*TernaryMultiPhaseBase::thermodynamic_factorBC()
         +TernaryMultiPhaseBase::L_CC_interp()*TernaryMultiPhaseBase::thermodynamic_factorCC());
}

//Interpolated first derivatives
Real
TernaryMultiPhaseBase::dL_BB_muB_interp() const
{
         
  return ( _dL_BB_muB_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_BB_muB_beta[_qp]   * _h_beta[_qp] 
         + _dL_BB_muB_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_BB_muB_delta[_qp]  * _h_delta[_qp]
         + _dL_BB_muB_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
TernaryMultiPhaseBase::dL_BC_muB_interp() const
{
  return ( _dL_BC_muB_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_BC_muB_beta[_qp]   * _h_beta[_qp] 
         + _dL_BC_muB_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_BC_muB_delta[_qp]  * _h_delta[_qp]
         + _dL_BC_muB_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
TernaryMultiPhaseBase::dL_CC_muB_interp() const
{
  return ( _dL_CC_muB_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_CC_muB_beta[_qp]   * _h_beta[_qp] 
         + _dL_CC_muB_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_CC_muB_delta[_qp]  * _h_delta[_qp]
         + _dL_CC_muB_epsilon[_qp]* _h_epsilon[_qp]);
}

/// w.r.t component C
Real
TernaryMultiPhaseBase::dL_BB_muC_interp() const
{
  return ( _dL_BB_muC_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_BB_muC_beta[_qp]   * _h_beta[_qp] 
         + _dL_BB_muC_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_BB_muC_delta[_qp]  * _h_delta[_qp]
         + _dL_BB_muC_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
TernaryMultiPhaseBase::dL_BC_muC_interp() const
{
  return ( _dL_BC_muC_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_BC_muC_beta[_qp]   * _h_beta[_qp] 
         + _dL_BC_muC_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_BC_muC_delta[_qp]  * _h_delta[_qp]
         + _dL_BC_muC_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
TernaryMultiPhaseBase::dL_CC_muC_interp() const
{
  return ( _dL_CC_muC_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_CC_muC_beta[_qp]   * _h_beta[_qp] 
         + _dL_CC_muC_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_CC_muC_delta[_qp]  * _h_delta[_qp]
         + _dL_CC_muC_epsilon[_qp]* _h_epsilon[_qp]);
}

//First derivatives  of each material property sum over all interpolation function
//This will be (# of property*# of phases)
//Onsger mobility L_BB

Real
TernaryMultiPhaseBase::sum_dh_L_BB_alpha() const
{
    return (_dhbeta_dphialpha[_qp]      * (_L_BB_beta[_qp]  - _L_BB_alpha[_qp])
             + _dhgamma_dphialpha[_qp]  * (_L_BB_gamma[_qp] - _L_BB_alpha[_qp])
             + _dhdelta_dphialpha[_qp]  * (_L_BB_delta[_qp] - _L_BB_alpha[_qp])
             + _dhepsilon_dphialpha[_qp]* (_L_BB_epsilon[_qp]- _L_BB_alpha[_qp]));
}

Real
TernaryMultiPhaseBase::sum_dh_L_BB_beta() const
{
    return (_dhalpha_dphibeta[_qp]     * (_L_BB_alpha[_qp] - _L_BB_beta[_qp])
             + _dhgamma_dphibeta[_qp]  * (_L_BB_gamma[_qp] - _L_BB_beta[_qp])
             + _dhdelta_dphibeta[_qp]  * (_L_BB_delta[_qp] - _L_BB_beta[_qp])
             + _dhepsilon_dphibeta[_qp]*(_L_BB_epsilon[_qp]- _L_BB_beta[_qp]));
}

Real
TernaryMultiPhaseBase::sum_dh_L_BB_gamma() const
{
    return (_dhalpha_dphigamma[_qp]     * (_L_BB_alpha[_qp] - _L_BB_gamma[_qp])
             + _dhbeta_dphigamma[_qp]   * (_L_BB_beta[_qp]  - _L_BB_gamma[_qp])
             + _dhdelta_dphigamma[_qp]  * (_L_BB_delta[_qp] - _L_BB_gamma[_qp])
             + _dhepsilon_dphigamma[_qp]* (_L_BB_epsilon[_qp]-_L_BB_gamma[_qp]));
}

Real
TernaryMultiPhaseBase::sum_dh_L_BB_delta() const
{
    return (_dhalpha_dphidelta[_qp]     * (_L_BB_alpha[_qp] - _L_BB_delta[_qp])
             + _dhbeta_dphidelta[_qp]   * (_L_BB_beta[_qp]  - _L_BB_delta[_qp])
             + _dhgamma_dphidelta[_qp]  * (_L_BB_gamma[_qp] - _L_BB_delta[_qp])
             + _dhepsilon_dphidelta[_qp]* (_L_BB_epsilon[_qp]-_L_BB_delta[_qp]));
}

Real
TernaryMultiPhaseBase::sum_dh_L_BB_epsilon() const
{
    return (_dhalpha_dphiepsilon[_qp]     * (_L_BB_alpha[_qp] - _L_BB_epsilon[_qp])
             + _dhbeta_dphiepsilon[_qp]   * (_L_BB_beta[_qp]  - _L_BB_epsilon[_qp])
             + _dhgamma_dphiepsilon[_qp]  * (_L_BB_gamma[_qp] - _L_BB_epsilon[_qp])
             + _dhdelta_dphiepsilon[_qp]*   (_L_BB_delta[_qp] - _L_BB_epsilon[_qp]));
}

//Onsger mobility L_BC

Real
TernaryMultiPhaseBase::sum_dh_L_BC_alpha() const
{
    return (_dhbeta_dphialpha[_qp]      * (_L_BC_beta[_qp]  - _L_BC_alpha[_qp])
             + _dhgamma_dphialpha[_qp]  * (_L_BC_gamma[_qp] - _L_BC_alpha[_qp])
             + _dhdelta_dphialpha[_qp]  * (_L_BC_delta[_qp] - _L_BC_alpha[_qp])
             + _dhepsilon_dphialpha[_qp]* (_L_BC_epsilon[_qp]- _L_BC_alpha[_qp]));
}

Real
TernaryMultiPhaseBase::sum_dh_L_BC_beta() const
{
    return (_dhalpha_dphibeta[_qp]     * (_L_BC_alpha[_qp] - _L_BC_beta[_qp])
             + _dhgamma_dphibeta[_qp]  * (_L_BC_gamma[_qp] - _L_BC_beta[_qp])
             + _dhdelta_dphibeta[_qp]  * (_L_BC_delta[_qp] - _L_BC_beta[_qp])
             + _dhepsilon_dphibeta[_qp]*(_L_BC_epsilon[_qp]- _L_BC_beta[_qp]));
}

Real
TernaryMultiPhaseBase::sum_dh_L_BC_gamma() const
{
    return  (_dhalpha_dphigamma[_qp]    * (_L_BC_alpha[_qp] - _L_BC_gamma[_qp])
             + _dhbeta_dphigamma[_qp]   * (_L_BC_beta[_qp]  - _L_BC_gamma[_qp])
             + _dhdelta_dphigamma[_qp]  * (_L_BC_delta[_qp] - _L_BC_gamma[_qp])
             + _dhepsilon_dphigamma[_qp]* (_L_BC_epsilon[_qp]-_L_BC_gamma[_qp]));
}

Real
TernaryMultiPhaseBase::sum_dh_L_BC_delta() const
{
    return (_dhalpha_dphidelta[_qp]     * (_L_BC_alpha[_qp] - _L_BC_delta[_qp])
             + _dhbeta_dphidelta[_qp]   * (_L_BC_beta[_qp]  - _L_BC_delta[_qp])
             + _dhgamma_dphidelta[_qp]  * (_L_BC_gamma[_qp] - _L_BC_delta[_qp])
             + _dhepsilon_dphidelta[_qp]* (_L_BC_epsilon[_qp]-_L_BC_delta[_qp]));
}

Real
TernaryMultiPhaseBase::sum_dh_L_BC_epsilon() const
{
    return (_dhalpha_dphiepsilon[_qp]     * (_L_BC_alpha[_qp] - _L_BC_epsilon[_qp])
             + _dhbeta_dphiepsilon[_qp]   * (_L_BC_beta[_qp]  - _L_BC_epsilon[_qp])
             + _dhgamma_dphiepsilon[_qp]  * (_L_BC_gamma[_qp] - _L_BC_epsilon[_qp])
             + _dhdelta_dphiepsilon[_qp]  * (_L_BC_delta[_qp] - _L_BC_epsilon[_qp]));
}

//Property::sum_dh_L_CC_phases()

Real
TernaryMultiPhaseBase::sum_dh_L_CC_alpha() const
{
    return (_dhbeta_dphialpha[_qp]      * (_L_CC_beta[_qp]  - _L_CC_alpha[_qp])
             + _dhgamma_dphialpha[_qp]  * (_L_CC_gamma[_qp] - _L_CC_alpha[_qp])
             + _dhdelta_dphialpha[_qp]  * (_L_CC_delta[_qp] - _L_CC_alpha[_qp])
             + _dhepsilon_dphialpha[_qp]* (_L_CC_epsilon[_qp]- _L_CC_alpha[_qp]));
}

Real
TernaryMultiPhaseBase::sum_dh_L_CC_beta() const
{
    return (_dhalpha_dphibeta[_qp]     * (_L_CC_alpha[_qp] - _L_CC_beta[_qp])
             + _dhgamma_dphibeta[_qp]  * (_L_CC_gamma[_qp] - _L_CC_beta[_qp])
             + _dhdelta_dphibeta[_qp]  * (_L_CC_delta[_qp] - _L_CC_beta[_qp])
             + _dhepsilon_dphibeta[_qp]*(_L_CC_epsilon[_qp]- _L_CC_beta[_qp]));
}

Real
TernaryMultiPhaseBase::sum_dh_L_CC_gamma() const
{
    return  (_dhalpha_dphigamma[_qp]    * (_L_CC_alpha[_qp] - _L_CC_gamma[_qp])
             + _dhbeta_dphigamma[_qp]   * (_L_CC_beta[_qp]  - _L_CC_gamma[_qp])
             + _dhdelta_dphigamma[_qp]  * (_L_CC_delta[_qp] - _L_CC_gamma[_qp])
             + _dhepsilon_dphigamma[_qp]* (_L_CC_epsilon[_qp]-_L_CC_gamma[_qp]));
}

Real
TernaryMultiPhaseBase::sum_dh_L_CC_delta() const
{
    return (_dhalpha_dphidelta[_qp]     * (_L_CC_alpha[_qp] - _L_CC_delta[_qp])
             + _dhbeta_dphidelta[_qp]   * (_L_CC_beta[_qp]  - _L_CC_delta[_qp])
             + _dhgamma_dphidelta[_qp]  * (_L_CC_gamma[_qp] - _L_CC_delta[_qp])
             + _dhepsilon_dphidelta[_qp]* (_L_CC_epsilon[_qp]-_L_CC_delta[_qp]));
}

Real
TernaryMultiPhaseBase::sum_dh_L_CC_epsilon() const
{
    return (_dhalpha_dphiepsilon[_qp]     * (_L_CC_alpha[_qp] - _L_CC_epsilon[_qp])
             + _dhbeta_dphiepsilon[_qp]   * (_L_CC_beta[_qp]  - _L_CC_epsilon[_qp])
             + _dhgamma_dphiepsilon[_qp]  * (_L_CC_gamma[_qp] - _L_CC_epsilon[_qp])
             + _dhdelta_dphiepsilon[_qp]*   (_L_CC_delta[_qp] - _L_CC_epsilon[_qp]));
}

Real
TernaryMultiPhaseBase::computeQpResidual(){
  return 0;
}

Real
TernaryMultiPhaseBase::computeQpJacobian(){
  return 0;
}

Real
TernaryMultiPhaseBase::computeQpOffDiagJacobian(unsigned int /*jvar*/){
  return 0;
}
