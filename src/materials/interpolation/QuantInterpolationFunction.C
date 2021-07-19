//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//* This Material class supplies the interpolation function,
//* for a three phase system

#include "QuantInterpolationFunction.h"
registerMooseObject("gibbsApp", QuantInterpolationFunction);

template <>
InputParameters
validParams<QuantInterpolationFunction>()
{
  InputParameters params = validParams<Material>();
  params.addCoupledVar("phase_alpha",0.0, "Phase field for alpha phase");
  params.addCoupledVar("phase_beta", 0.0, "Phase field for beta phase");
  params.addCoupledVar("phase_gamma",0.0, "Phase field for gamma phase");
  params.addCoupledVar("phase_delta", 0.0, "Phase field for delta phase");
  params.addCoupledVar("phase_epsilon", 0.0, "Phase field for epsilon phase");
  return params;
}

QuantInterpolationFunction::QuantInterpolationFunction(const InputParameters & parameters)
  : Material(parameters),
    // Declare that this material is going to provide a Real
    // valued properties named "_h_alpha,_h_beta,_h_gamma".
    _h_alpha(declareProperty<Real>("h_alpha")),
    _h_beta(declareProperty<Real>("h_beta")),
    _h_gamma(declareProperty<Real>("h_gamma")),
    _h_delta(declareProperty<Real>("h_delta")),
    _h_epsilon(declareProperty<Real>("h_epsilon")),
    
    //first derivatives of the interpolation function
    _dhalpha_dphialpha(declareProperty<Real>("dhalpha_dphialpha")),
    _dhbeta_dphibeta(declareProperty<Real>("dhbeta_dphibeta")),
    _dhgamma_dphigamma(declareProperty<Real>("dhgamma_dphigamma")),
    _dhdelta_dphidelta(declareProperty<Real>("dhdelta_dphidelta")),
    _dhepsilon_dphiepsilon(declareProperty<Real>("dhepsilon_dphiepsilon")),
    
    // phase_alpha
    _dhbeta_dphialpha(declareProperty<Real>("dhbeta_dphialpha")),
    _dhgamma_dphialpha(declareProperty<Real>("dhgamma_dphialpha")),
    _dhdelta_dphialpha(declareProperty<Real>("dhdelta_dphialpha")),
    _dhepsilon_dphialpha(declareProperty<Real>("dhepsilon_dphialpha")),
    
    // phase_beta
    _dhalpha_dphibeta(declareProperty<Real>("dhalpha_dphibeta")),
    _dhgamma_dphibeta(declareProperty<Real>("dhgamma_dphibeta")),
    _dhdelta_dphibeta(declareProperty<Real>("dhdelta_dphibeta")),
    _dhepsilon_dphibeta(declareProperty<Real>("dhepsilon_dphibeta")),
    
    // phase_gamma
    _dhalpha_dphigamma(declareProperty<Real>("dhalpha_dphigamma")),
    _dhbeta_dphigamma(declareProperty<Real>("dhbeta_dphigamma")),
    _dhdelta_dphigamma(declareProperty<Real>("dhdelta_dphigamma")),
    _dhepsilon_dphigamma(declareProperty<Real>("dhepsilon_dphigamma")),
    
    // phase_delta
    _dhalpha_dphidelta(declareProperty<Real>("dhalpha_dphidelta")),
    _dhbeta_dphidelta(declareProperty<Real>("dhbeta_dphidelta")),
    _dhgamma_dphidelta(declareProperty<Real>("dhgamma_dphidelta")), 
    _dhepsilon_dphidelta(declareProperty<Real>("dhepsilon_dphidelta")),
    
    // phase_epsilon
    _dhalpha_dphiepsilon(declareProperty<Real>("dhalpha_dphiepsilon")),
    _dhbeta_dphiepsilon(declareProperty<Real>("dhbeta_dphiepsilon")),
    _dhgamma_dphiepsilon(declareProperty<Real>("dhgamma_dphiepsilon")), 
    _dhdelta_dphiepsilon(declareProperty<Real>("dhdelta_dphiepsilon")),
              
    // phase_alpha (Diagonal terms)
    _d2hbeta_dphialpha2(declareProperty<Real>("d2hbeta_dphialpha2")),
    _d2hgamma_dphialpha2(declareProperty<Real>("d2hgamma_dphialpha2")),
    _d2hdelta_dphialpha2(declareProperty<Real>("d2hdelta_dphialpha2")),
    _d2hepsilon_dphialpha2(declareProperty<Real>("d2hepsilon_dphialpha2")),
    
    //phase_beta
    _d2halpha_dphibeta2(declareProperty<Real>("d2halpha_dphibeta2")),
    _d2hgamma_dphibeta2(declareProperty<Real>("d2hgamma_dphibeta2")),
    _d2hdelta_dphibeta2(declareProperty<Real>("d2hdelta_dphibeta2")),
    _d2hepsilon_dphibeta2(declareProperty<Real>("d2hepsilon_dphibeta2")),
    
    // phase_gamma
    _d2halpha_dphigamma2(declareProperty<Real>("d2halpha_dphigamma2")),
    _d2hbeta_dphigamma2(declareProperty<Real>("d2hbeta_dphigamma2")),
    _d2hdelta_dphigamma2(declareProperty<Real>("d2hdelta_dphigamma2")),
    _d2hepsilon_dphigamma2(declareProperty<Real>("d2hepsilon_dphigamma2")),
    
    // phase_delta
    _d2halpha_dphidelta2(declareProperty<Real>("d2halpha_dphidelta2")),
    _d2hbeta_dphidelta2(declareProperty<Real>("d2hbeta_dphidelta2")),
    _d2hgamma_dphidelta2(declareProperty<Real>("d2hgamma_dphidelta2")),
    _d2hepsilon_dphidelta2(declareProperty<Real>("d2hepsilon_dphidelta2")),
    
    // phase_epsilon
    _d2halpha_dphiepsilon2(declareProperty<Real>("d2halpha_dphiepsilon2")),
    _d2hbeta_dphiepsilon2(declareProperty<Real>("d2hbeta_dphiepsilon2")),
    _d2hgamma_dphiepsilon2(declareProperty<Real>("d2hgamma_dphiepsilon2")),
    _d2hdelta_dphiepsilon2(declareProperty<Real>("d2hdelta_dphiepsilon2")),
       
    // phase_alpha (off-diagonal) 
    //h_{beta,alpha}
    _d2hbeta_dphialpha_dphibeta(declareProperty<Real>("d2hbeta_dphialpha_dphibeta")),
    _d2hbeta_dphialpha_dphigamma(declareProperty<Real>("d2hbeta_dphialpha_dphigamma")),
    _d2hbeta_dphialpha_dphidelta(declareProperty<Real>("d2hbeta_dphialpha_dphidelta")),
    _d2hbeta_dphialpha_dphiepsilon(declareProperty<Real>("d2hbeta_dphialpha_dphiepsilon")),
    
    //h_{gamma,alpha}
    _d2hgamma_dphialpha_dphibeta(declareProperty<Real>("d2hgamma_dphialpha_dphibeta")),
    _d2hgamma_dphialpha_dphigamma(declareProperty<Real>("d2hgamma_dphialpha_dphigamma")),
    _d2hgamma_dphialpha_dphidelta(declareProperty<Real>("d2hgamma_dphialpha_dphidelta")),
    _d2hgamma_dphialpha_dphiepsilon(declareProperty<Real>("d2hgamma_dphialpha_dphiepsilon")),
    
    //h_{delta,alpha}
    _d2hdelta_dphialpha_dphibeta(declareProperty<Real>("d2hdelta_dphialpha_dphibeta")),
    _d2hdelta_dphialpha_dphigamma(declareProperty<Real>("d2hdelta_dphialpha_dphigamma")),
    _d2hdelta_dphialpha_dphidelta(declareProperty<Real>("d2hdelta_dphialpha_dphidelta")),
    _d2hdelta_dphialpha_dphiepsilon(declareProperty<Real>("d2hdelta_dphialpha_dphiepsilon")),
    
     //h_{epsilon,alpha}
    _d2hepsilon_dphialpha_dphibeta(declareProperty<Real>("d2hepsilon_dphialpha_dphibeta")),
    _d2hepsilon_dphialpha_dphigamma(declareProperty<Real>("d2hepsilon_dphialpha_dphigamma")),
    _d2hepsilon_dphialpha_dphidelta(declareProperty<Real>("d2hepsilon_dphialpha_dphidelta")),
    _d2hepsilon_dphialpha_dphiepsilon(declareProperty<Real>("d2hepsilon_dphialpha_dphiepsilon")),
        
    //phase_beta
    //h_{alpha,beta}
    _d2halpha_dphibeta_dphialpha(declareProperty<Real>("d2halpha_dphibeta_dphialpha")),
    _d2halpha_dphibeta_dphigamma(declareProperty<Real>("d2halpha_dphibeta_dphigamma")),
    _d2halpha_dphibeta_dphidelta(declareProperty<Real>("d2halpha_dphibeta_dphidelta")),
    _d2halpha_dphibeta_dphiepsilon(declareProperty<Real>("d2halpha_dphibeta_dphiepsilon")),
    
     //h_{gamma,beta}
    _d2hgamma_dphibeta_dphialpha(declareProperty<Real>("d2hgamma_dphibeta_dphialpha")),
    _d2hgamma_dphibeta_dphigamma(declareProperty<Real>("d2hgamma_dphibeta_dphigamma")),
    _d2hgamma_dphibeta_dphidelta(declareProperty<Real>("d2hgamma_dphibeta_dphidelta")),
    _d2hgamma_dphibeta_dphiepsilon(declareProperty<Real>("d2hgamma_dphibeta_dphiepsilon")),
    
     //h_{delta,beta}
    _d2hdelta_dphibeta_dphialpha(declareProperty<Real>("d2hdelta_dphibeta_dphialpha")),
    _d2hdelta_dphibeta_dphigamma(declareProperty<Real>("d2hdelta_dphibeta_dphigamma")),
    _d2hdelta_dphibeta_dphidelta(declareProperty<Real>("d2hdelta_dphibeta_dphidelta")),
    _d2hdelta_dphibeta_dphiepsilon(declareProperty<Real>("d2hdelta_dphibeta_dphiepsilon")),
    
    //h_{epsilon,beta}
    _d2hepsilon_dphibeta_dphialpha(declareProperty<Real>("d2hepsilon_dphibeta_dphialpha")),
    _d2hepsilon_dphibeta_dphigamma(declareProperty<Real>("d2hepsilon_dphibeta_dphigamma")),
    _d2hepsilon_dphibeta_dphidelta(declareProperty<Real>("d2hepsilon_dphibeta_dphidelta")),
    _d2hepsilon_dphibeta_dphiepsilon(declareProperty<Real>("d2hepsilon_dphibeta_dphiepsilon")),
    
    //phase_gamma
    //h_{alpha,gamma}
    _d2halpha_dphigamma_dphibeta(declareProperty<Real>("d2halpha_dphigamma_dphibeta")),
    _d2halpha_dphigamma_dphialpha(declareProperty<Real>("d2halpha_dphigamma_dphialpha")),
    _d2halpha_dphigamma_dphidelta(declareProperty<Real>("d2halpha_dphigamma_dphidelta")),
    _d2halpha_dphigamma_dphiepsilon(declareProperty<Real>("d2halpha_dphigamma_dphiepsilon")),
    
    //h_{beta,gamma}
    _d2hbeta_dphigamma_dphialpha(declareProperty<Real>("d2hbeta_dphigamma_dphialpha")),
    _d2hbeta_dphigamma_dphibeta(declareProperty<Real>("d2hbeta_dphigamma_dphibeta")),
    _d2hbeta_dphigamma_dphidelta(declareProperty<Real>("d2hbeta_dphigamma_dphidelta")),
    _d2hbeta_dphigamma_dphiepsilon(declareProperty<Real>("d2hbeta_dphigamma_dphiepsilon")),
    
    //h_{delta,gamma}
    _d2hdelta_dphigamma_dphialpha(declareProperty<Real>("d2hdelta_dphigamma_dphialpha")),
    _d2hdelta_dphigamma_dphibeta(declareProperty<Real>("d2hdelta_dphigamma_dphibeta")),
    _d2hdelta_dphigamma_dphidelta(declareProperty<Real>("d2hdelta_dphigamma_dphidelta")), 
    _d2hdelta_dphigamma_dphiepsilon(declareProperty<Real>("d2hdelta_dphigamma_dphiepsilon")),
    
    //h_{epsilon,gamma}
    _d2hepsilon_dphigamma_dphialpha(declareProperty<Real>("d2hepsilon_dphigamma_dphialpha")),
    _d2hepsilon_dphigamma_dphibeta(declareProperty<Real>("d2hepsilon_dphigamma_dphibeta")),
    _d2hepsilon_dphigamma_dphidelta(declareProperty<Real>("d2hepsilon_dphigamma_dphidelta")), 
    _d2hepsilon_dphigamma_dphiepsilon(declareProperty<Real>("d2hepsilon_dphigamma_dphiepsilon")),   
      
    //phase_delta
    //h_{alpha,delta}
    _d2halpha_dphidelta_dphialpha(declareProperty<Real>("d2halpha_dphidelta_dphialpha")),
    _d2halpha_dphidelta_dphibeta(declareProperty<Real>("d2halpha_dphidelta_dphibeta")),
    _d2halpha_dphidelta_dphigamma(declareProperty<Real>("d2halpha_dphidelta_dphigamma")),
    _d2halpha_dphidelta_dphiepsilon(declareProperty<Real>("d2halpha_dphidelta_dphiepsilon")),
    
    //h_{beta,delta}
    _d2hbeta_dphidelta_dphialpha(declareProperty<Real>("d2hbeta_dphidelta_dphialpha")),
    _d2hbeta_dphidelta_dphibeta(declareProperty<Real>("d2hbeta_dphidelta_dphibeta")),
    _d2hbeta_dphidelta_dphigamma(declareProperty<Real>("d2hbeta_dphidelta_dphigamma")),
    _d2hbeta_dphidelta_dphiepsilon(declareProperty<Real>("d2hbeta_dphidelta_dphiepsilon")),
        
     //h_{gamma,delta}
    _d2hgamma_dphidelta_dphialpha(declareProperty<Real>("d2hgamma_dphidelta_dphialpha")),
    _d2hgamma_dphidelta_dphibeta(declareProperty<Real>("d2hgamma_dphidelta_dphibeta")),
    _d2hgamma_dphidelta_dphigamma(declareProperty<Real>("d2hgamma_dphidelta_dphigamma")),
    _d2hgamma_dphidelta_dphiepsilon(declareProperty<Real>("d2hgamma_dphidelta_dphiepsilon")),
    
     //h_{epsilon,delta}
    _d2hepsilon_dphidelta_dphialpha(declareProperty<Real>("d2hepsilon_dphidelta_dphialpha")),
    _d2hepsilon_dphidelta_dphibeta(declareProperty<Real>("d2hepsilon_dphidelta_dphibeta")),
    _d2hepsilon_dphidelta_dphigamma(declareProperty<Real>("d2hepsilon_dphidelta_dphigamma")),
    _d2hepsilon_dphidelta_dphiepsilon(declareProperty<Real>("d2hepsilon_dphidelta_dphiepsilon")),
    
    //phase_epsilon
    //h_{alpha,epsilon}
    _d2halpha_dphiepsilon_dphialpha(declareProperty<Real>("d2halpha_dphiepsilon_dphialpha")),
    _d2halpha_dphiepsilon_dphibeta(declareProperty<Real>("d2halpha_dphiepsilon_dphibeta")),
    _d2halpha_dphiepsilon_dphigamma(declareProperty<Real>("d2halpha_dphiepsilon_dphigamma")),
    _d2halpha_dphiepsilon_dphidelta(declareProperty<Real>("d2halpha_dphiepsilon_dphidelta")),
    
    //h_{beta,epsilon}
    _d2hbeta_dphiepsilon_dphialpha(declareProperty<Real>("d2hbeta_dphiepsilon_dphialpha")),
    _d2hbeta_dphiepsilon_dphibeta(declareProperty<Real>("d2hbeta_dphiepsilon_dphibeta")),
    _d2hbeta_dphiepsilon_dphigamma(declareProperty<Real>("d2hbeta_dphiepsilon_dphigamma")),
    _d2hbeta_dphiepsilon_dphidelta(declareProperty<Real>("d2hbeta_dphiepsilon_dphidelta")),
        
     //h_{gamma,epsilon}
    _d2hgamma_dphiepsilon_dphialpha(declareProperty<Real>("d2hgamma_dphiepsilon_dphialpha")),
    _d2hgamma_dphiepsilon_dphibeta(declareProperty<Real>("d2hgamma_dphiepsilon_dphibeta")),
    _d2hgamma_dphiepsilon_dphigamma(declareProperty<Real>("d2hgamma_dphiepsilon_dphigamma")),
    _d2hgamma_dphiepsilon_dphidelta(declareProperty<Real>("d2hgamma_dphiepsilon_dphidelta")),
    
     //h_{delta, epsilon}
    _d2hdelta_dphiepsilon_dphialpha(declareProperty<Real>("d2hdelta_dphiepsilon_dphialpha")),
    _d2hdelta_dphiepsilon_dphibeta(declareProperty<Real>("d2hdelta_dphiepsilon_dphibeta")),
    _d2hdelta_dphiepsilon_dphigamma(declareProperty<Real>("d2hdelta_dphiepsilon_dphigamma")),
    _d2hdelta_dphiepsilon_dphidelta(declareProperty<Real>("d2hdelta_dphiepsilon_dphidelta")),
    
    _phase_alpha(coupledValue("phase_alpha")),
    _phase_beta(coupledValue("phase_beta")),
    _phase_gamma(coupledValue("phase_gamma")),
    _phase_delta(coupledValue("phase_delta")),
    _phase_epsilon(coupledValue("phase_epsilon"))
{
}

void
QuantInterpolationFunction::computeQpProperties()
{
    //Sum of the square of the interpolation function
    Real _getsumSquare = 1.0/(std::pow(_phase_alpha[_qp],2.0) 
    + std::pow(_phase_beta[_qp],2.0) + std::pow(_phase_gamma[_qp],2.0) + std::pow(_phase_delta[_qp], 2.0)
    + std::pow(_phase_epsilon[_qp],2.0));
           
    //The interpolation function for each phase
   _h_alpha[_qp] = std::pow(_phase_alpha[_qp], 2.0) * _getsumSquare;
   _h_beta[_qp]  = std::pow(_phase_beta[_qp], 2.0)  * _getsumSquare;
   _h_gamma[_qp] = std::pow(_phase_gamma[_qp], 2.0) * _getsumSquare;
   _h_delta[_qp] = std::pow(_phase_delta[_qp], 2.0) * _getsumSquare;
   _h_epsilon[_qp] = std::pow(_phase_epsilon[_qp], 2.0) *_getsumSquare;
   
   //The first derivative of the interpolation function
   //Note that h_{alpha}(\phi_{alpha}, \phi_{beta}, \phi_{gamma})
   //and h_{\alpha} + h_{\beta} + h_{\gamma} = 1
   
   //Diagonal terms of the derivative of the matrix
   _dhalpha_dphialpha[_qp] = 2.0 * _phase_alpha[_qp] * (1.0 - _h_alpha[_qp]) * _getsumSquare;
   _dhbeta_dphibeta[_qp] =   2.0 * _phase_beta[_qp]  * (1.0 - _h_beta[_qp])  * _getsumSquare;   
   _dhgamma_dphigamma[_qp] = 2.0 * _phase_gamma[_qp] * (1.0 - _h_gamma[_qp]) * _getsumSquare;
   _dhdelta_dphidelta[_qp] = 2.0 * _phase_delta[_qp] * (1.0 - _h_delta[_qp]) * _getsumSquare;
   _dhepsilon_dphiepsilon[_qp] = 2.0 * _phase_epsilon[_qp] * (1.0 - _h_epsilon[_qp]) * _getsumSquare;
      
   // OffDiagonals alpha column
   _dhbeta_dphialpha[_qp]  = -2.0 * _h_beta[_qp]  * _phase_alpha[_qp] * _getsumSquare;
   _dhgamma_dphialpha[_qp] = -2.0 * _h_gamma[_qp] * _phase_alpha[_qp] * _getsumSquare;
   _dhdelta_dphialpha[_qp] = -2.0 * _h_delta[_qp] * _phase_alpha[_qp] * _getsumSquare;
   _dhepsilon_dphialpha[_qp]= -2.0 * _h_epsilon[_qp] * _phase_alpha[_qp] * _getsumSquare;
    
   // OffDiagonals beta column
   _dhalpha_dphibeta[_qp] = -2.0 * _h_alpha[_qp] * _phase_beta[_qp] * _getsumSquare;
   _dhgamma_dphibeta[_qp] = -2.0 * _h_gamma[_qp] * _phase_beta[_qp] * _getsumSquare;
   _dhdelta_dphibeta[_qp] = -2.0 * _h_delta[_qp] * _phase_beta[_qp] * _getsumSquare;
   _dhepsilon_dphibeta[_qp] = -2.0* _h_epsilon[_qp] * _phase_beta[_qp] * _getsumSquare;
   
   // OffDiagonals gamma column
   _dhalpha_dphigamma[_qp] = -2.0 * _h_alpha[_qp] * _phase_gamma[_qp] * _getsumSquare;
   _dhbeta_dphigamma[_qp]  = -2.0 * _h_beta[_qp]  * _phase_gamma[_qp] * _getsumSquare;
   _dhdelta_dphigamma[_qp] = -2.0 * _h_delta[_qp] * _phase_gamma[_qp] * _getsumSquare;
   _dhepsilon_dphigamma[_qp] = -2.0 * _h_epsilon[_qp] * _phase_gamma[_qp] * _getsumSquare;
   
    // OffDiagonals delta column
   _dhalpha_dphidelta[_qp] = -2.0 * _h_alpha[_qp] * _phase_delta[_qp] * _getsumSquare;
   _dhbeta_dphidelta[_qp]  = -2.0 * _h_beta[_qp]  * _phase_delta[_qp] * _getsumSquare;
   _dhgamma_dphidelta[_qp] = -2.0 * _h_gamma[_qp] * _phase_delta[_qp] * _getsumSquare;
   _dhepsilon_dphidelta[_qp] = -2.0 * _h_epsilon[_qp] * _phase_delta[_qp] * _getsumSquare;
   
    // OffDiagonals epsilon column
   _dhalpha_dphiepsilon[_qp] = -2.0 * _h_alpha[_qp] * _phase_epsilon[_qp] * _getsumSquare;
   _dhbeta_dphiepsilon[_qp]  = -2.0 * _h_beta[_qp]  * _phase_epsilon[_qp] * _getsumSquare;
   _dhgamma_dphiepsilon[_qp] = -2.0 * _h_gamma[_qp] * _phase_epsilon[_qp] * _getsumSquare;
   _dhdelta_dphiepsilon[_qp] = -2.0 * _h_delta[_qp] * _phase_epsilon[_qp] * _getsumSquare;
   
   //*******************************Second derivatives************************************************************//
   
   // derivate of off-diagonal terms with respect to same phase alpha
   _d2hbeta_dphialpha2[_qp] = -2.0*(2.0* _dhbeta_dphialpha[_qp] * _phase_alpha[_qp] + _h_beta[_qp])* _getsumSquare;
   _d2hgamma_dphialpha2[_qp]= -2.0*(2.0* _dhgamma_dphialpha[_qp]* _phase_alpha[_qp] + _h_gamma[_qp])* _getsumSquare;
   _d2hdelta_dphialpha2[_qp]= -2.0*(2.0* _dhdelta_dphialpha[_qp]* _phase_alpha[_qp] + _h_delta[_qp])* _getsumSquare;
   _d2hepsilon_dphialpha2[_qp]= -2.0*(2.0* _dhepsilon_dphialpha[_qp]* _phase_alpha[_qp] + _h_epsilon[_qp])* _getsumSquare;
 
   // derivate of off-diagonal terms with respect to same phase beta
   _d2halpha_dphibeta2[_qp] = -2.0*(2.0* _dhalpha_dphibeta[_qp] * _phase_beta[_qp] + _h_alpha[_qp])* _getsumSquare;
   _d2hgamma_dphibeta2[_qp] = -2.0*(2.0* _dhgamma_dphibeta[_qp] * _phase_beta[_qp] + _h_gamma[_qp])* _getsumSquare;
   _d2hdelta_dphibeta2[_qp] = -2.0*(2.0* _dhdelta_dphibeta[_qp] * _phase_beta[_qp] + _h_delta[_qp])* _getsumSquare;
   _d2hepsilon_dphibeta2[_qp] = -2.0*(2.0* _dhepsilon_dphibeta[_qp] * _phase_beta[_qp] + _h_epsilon[_qp])* _getsumSquare;
     
   //derivate of off-diagonal terms with respect to same phase gamma 
   _d2halpha_dphigamma2[_qp] = -2.0*(2.0* _dhalpha_dphigamma[_qp]* _phase_gamma[_qp] + _h_alpha[_qp])*_getsumSquare;
   _d2hbeta_dphigamma2[_qp]  = -2.0*(2.0* _dhbeta_dphigamma[_qp] * _phase_gamma[_qp] + _h_beta[_qp])* _getsumSquare; 
   _d2hdelta_dphigamma2[_qp] = -2.0*(2.0* _dhdelta_dphigamma[_qp]* _phase_gamma[_qp] + _h_delta[_qp])*_getsumSquare;
   _d2hepsilon_dphigamma2[_qp] = -2.0*(2.0* _dhepsilon_dphigamma[_qp] * _phase_gamma[_qp] + _h_epsilon[_qp])* _getsumSquare;
   
   //derivate of off-diagonal terms with respect to same phase delta 
   _d2halpha_dphidelta2[_qp] = -2.0*(2.0* _dhalpha_dphidelta[_qp]* _phase_delta[_qp] + _h_alpha[_qp])*_getsumSquare;
   _d2hbeta_dphidelta2[_qp]  = -2.0*(2.0* _dhbeta_dphidelta[_qp] * _phase_delta[_qp] + _h_beta[_qp])* _getsumSquare; 
   _d2hgamma_dphidelta2[_qp] = -2.0*(2.0* _dhdelta_dphidelta[_qp]* _phase_delta[_qp] + _h_gamma[_qp])*_getsumSquare;
   _d2hepsilon_dphidelta2[_qp] = -2.0*(2.0* _dhepsilon_dphidelta[_qp] * _phase_delta[_qp] + _h_epsilon[_qp])* _getsumSquare;
   
   //derivate of off-diagonal terms with respect to same phase epsilon 
   _d2halpha_dphiepsilon2[_qp] = -2.0*(2.0* _dhalpha_dphiepsilon[_qp]* _phase_epsilon[_qp] + _h_alpha[_qp])*_getsumSquare;
   _d2hbeta_dphiepsilon2[_qp]  = -2.0*(2.0* _dhbeta_dphiepsilon[_qp] * _phase_epsilon[_qp] + _h_beta[_qp]) *_getsumSquare; 
   _d2hgamma_dphiepsilon2[_qp] = -2.0*(2.0* _dhgamma_dphiepsilon[_qp]* _phase_epsilon[_qp] + _h_gamma[_qp])*_getsumSquare;
   _d2hdelta_dphiepsilon2[_qp] = -2.0*(2.0* _dhdelta_dphiepsilon[_qp]* _phase_epsilon[_qp] + _h_delta[_qp])*_getsumSquare;     
   
   //****************************************For phase alpha********************************************************// 
   
   //  For h_{beta,alpha}
   _d2hbeta_dphialpha_dphibeta[_qp] = 
    -2.0*(_dhbeta_dphialpha[_qp] * _phase_beta[_qp]  + _dhbeta_dphibeta[_qp] * _phase_alpha[_qp]) * _getsumSquare;
   
   _d2hbeta_dphialpha_dphigamma[_qp] = 
    -2.0*(_dhbeta_dphigamma[_qp] * _phase_alpha[_qp] + _dhbeta_dphialpha[_qp] * _phase_gamma[_qp])* _getsumSquare;
   
   _d2hbeta_dphialpha_dphidelta[_qp] =
    -2.0*(_dhbeta_dphidelta[_qp] * _phase_alpha[_qp]  + _dhbeta_dphialpha[_qp] * _phase_delta[_qp])* _getsumSquare;
    
    _d2hbeta_dphialpha_dphiepsilon[_qp] =
    -2.0*(_dhbeta_dphiepsilon[_qp] * _phase_alpha[_qp]  + _dhbeta_dphialpha[_qp] * _phase_epsilon[_qp])* _getsumSquare;
    
    // For h_{gamma,alpha}
    
   _d2hgamma_dphialpha_dphibeta[_qp] = 
    -2.0*(_dhgamma_dphibeta[_qp] * _phase_alpha[_qp] + _dhgamma_dphialpha[_qp] *  _phase_beta[_qp]) * _getsumSquare;
    
   _d2hgamma_dphialpha_dphigamma[_qp] =
    -2.0*(_dhgamma_dphialpha[_qp] * _phase_gamma[_qp] + _dhgamma_dphigamma[_qp] * _phase_alpha[_qp]) *_getsumSquare;
        
    _d2hgamma_dphialpha_dphidelta[_qp] = 
    -2.0*(_dhgamma_dphidelta[_qp] * _phase_alpha[_qp] + _dhgamma_dphialpha[_qp] * _phase_delta[_qp]) * _getsumSquare;
    
    _d2hgamma_dphialpha_dphiepsilon[_qp] = 
    -2.0*(_dhgamma_dphiepsilon[_qp] * _phase_alpha[_qp] + _dhgamma_dphialpha[_qp] * _phase_epsilon[_qp]) * _getsumSquare;
        
    // For h_{delta,alpha}
    
   _d2hdelta_dphialpha_dphibeta[_qp] = 
    -2.0*(_dhdelta_dphibeta[_qp] * _phase_alpha[_qp] + _dhdelta_dphialpha[_qp]  * _phase_beta[_qp]) *_getsumSquare;
    
   _d2hdelta_dphialpha_dphigamma[_qp] =
    -2.0*(_dhdelta_dphigamma[_qp] * _phase_alpha[_qp] + _dhdelta_dphialpha[_qp] * _phase_gamma[_qp])*_getsumSquare;
        
    _d2hdelta_dphialpha_dphidelta[_qp] = 
    -2.0*(_dhdelta_dphialpha[_qp] * _phase_delta[_qp] + _dhdelta_dphidelta[_qp] * _phase_alpha[_qp])* _getsumSquare;
    
     _d2hdelta_dphialpha_dphiepsilon[_qp] = 
    -2.0*(_dhdelta_dphiepsilon[_qp] * _phase_alpha[_qp] + _dhdelta_dphialpha[_qp] * _phase_epsilon[_qp])* _getsumSquare;
    
    // For h_{epsilon,alpha}
    
     _d2hepsilon_dphialpha_dphibeta[_qp] = 
    -2.0*(_dhepsilon_dphibeta[_qp] * _phase_alpha[_qp] + _dhepsilon_dphialpha[_qp]  * _phase_beta[_qp]) *_getsumSquare;
    
     _d2hepsilon_dphialpha_dphigamma[_qp] = 
    -2.0*(_dhepsilon_dphigamma[_qp] * _phase_alpha[_qp] + _dhepsilon_dphialpha[_qp] * _phase_gamma[_qp]) *_getsumSquare;
    
     _d2hepsilon_dphialpha_dphidelta[_qp] = 
    -2.0*(_dhepsilon_dphidelta[_qp] * _phase_alpha[_qp] + _dhepsilon_dphialpha[_qp] * _phase_delta[_qp]) *_getsumSquare;
    
     _d2hepsilon_dphialpha_dphiepsilon[_qp] = 
    -2.0*(_dhepsilon_dphialpha[_qp] * _phase_epsilon[_qp] + _dhepsilon_dphiepsilon[_qp] * _phase_alpha[_qp]) *_getsumSquare;
    
   //****************************************For phase beta********************************************************//
   
   // For h_{alpha,beta}
   _d2halpha_dphibeta_dphialpha[_qp] =
    -2.0*(_dhalpha_dphibeta[_qp] * _phase_alpha[_qp] + _dhalpha_dphialpha[_qp] * _phase_beta[_qp]) * _getsumSquare;
    
   _d2halpha_dphibeta_dphigamma[_qp] = 
   -2.0*(_dhalpha_dphigamma[_qp] * _phase_beta[_qp] + _dhalpha_dphibeta[_qp] * _phase_gamma[_qp]) * _getsumSquare;
   
   _d2halpha_dphibeta_dphidelta[_qp] = 
   -2.0*(_dhalpha_dphidelta[_qp] * _phase_beta[_qp] + _dhalpha_dphibeta[_qp] * _phase_delta[_qp]) * _getsumSquare;
   
   _d2halpha_dphibeta_dphiepsilon[_qp] = 
   -2.0*(_dhalpha_dphiepsilon[_qp] * _phase_beta[_qp] + _dhalpha_dphibeta[_qp] * _phase_epsilon[_qp]) * _getsumSquare;
   
   //For h_{gamma,beta}
    
   _d2hgamma_dphibeta_dphialpha[_qp] = 
    -2.0*(_dhgamma_dphialpha[_qp] * _phase_beta[_qp] + _dhgamma_dphibeta[_qp] * _phase_alpha[_qp]) *_getsumSquare; 
    
   _d2hgamma_dphibeta_dphigamma[_qp] = 
    -2.0*(_dhgamma_dphibeta[_qp] * _phase_gamma[_qp] + _dhgamma_dphigamma[_qp] * _phase_beta[_qp]) * _getsumSquare;
    
   _d2hgamma_dphibeta_dphidelta[_qp] = 
    -2.0*(_dhgamma_dphidelta[_qp] * _phase_beta[_qp] + _dhgamma_dphibeta[_qp] * _phase_delta[_qp]) * _getsumSquare;
    
   _d2hgamma_dphibeta_dphiepsilon[_qp] = 
    -2.0*(_dhgamma_dphiepsilon[_qp] * _phase_beta[_qp] + _dhgamma_dphibeta[_qp] * _phase_epsilon[_qp]) * _getsumSquare;
    
  // For h_{delta,beta}
  
   _d2hdelta_dphibeta_dphialpha[_qp] = 
    -2.0*(_dhdelta_dphialpha[_qp] * _phase_beta[_qp] + _dhdelta_dphibeta[_qp] * _phase_alpha[_qp]) *_getsumSquare; 
    
   _d2hdelta_dphibeta_dphigamma[_qp] = 
    -2.0*(_dhdelta_dphigamma[_qp] * _phase_beta[_qp] + _dhdelta_dphibeta[_qp] * _phase_gamma[_qp]) * _getsumSquare;
    
   _d2hdelta_dphibeta_dphidelta[_qp] = 
    -2.0*(_dhdelta_dphibeta[_qp] * _phase_delta[_qp] + _dhdelta_dphidelta[_qp] * _phase_beta[_qp]) * _getsumSquare;
    
   _d2hdelta_dphibeta_dphiepsilon[_qp] = 
    -2.0*(_dhdelta_dphiepsilon[_qp] * _phase_beta[_qp] + _dhdelta_dphibeta[_qp] * _phase_epsilon[_qp]) * _getsumSquare;
    
    // For h_{epsilon,beta}
  
   _d2hepsilon_dphibeta_dphialpha[_qp] = 
    -2.0*(_dhepsilon_dphialpha[_qp] * _phase_beta[_qp] + _dhepsilon_dphibeta[_qp] * _phase_alpha[_qp]) *_getsumSquare; 
    
   _d2hepsilon_dphibeta_dphigamma[_qp] = 
    -2.0*(_dhepsilon_dphigamma[_qp] * _phase_beta[_qp] + _dhepsilon_dphibeta[_qp] * _phase_gamma[_qp]) * _getsumSquare;
    
   _d2hepsilon_dphibeta_dphidelta[_qp] = 
    -2.0*(_dhepsilon_dphidelta[_qp] * _phase_beta[_qp] + _dhepsilon_dphibeta[_qp] * _phase_delta[_qp]) * _getsumSquare;
    
   _d2hepsilon_dphibeta_dphiepsilon[_qp] = 
    -2.0*(_dhepsilon_dphibeta[_qp] * _phase_epsilon[_qp] + _dhepsilon_dphiepsilon[_qp] * _phase_beta[_qp])* _getsumSquare;
    
   
    //****************************************For phase gamma********************************************************//
    
   // For h_{alpha,gamma}
   
   _d2halpha_dphigamma_dphialpha[_qp] = 
     -2.0*(_dhalpha_dphigamma[_qp] * _phase_alpha[_qp] + _dhalpha_dphialpha[_qp] * _phase_gamma[_qp]) * _getsumSquare;
     
   _d2halpha_dphigamma_dphibeta[_qp] = 
    -2.0*(_dhalpha_dphibeta[_qp]* _phase_gamma[_qp] + _dhalpha_dphigamma[_qp] * _phase_beta[_qp]) *_getsumSquare;
   
   _d2halpha_dphigamma_dphidelta[_qp] = 
    -2.0*(_dhalpha_dphidelta[_qp]* _phase_gamma[_qp] + _dhalpha_dphigamma[_qp] * _phase_delta[_qp]) *_getsumSquare; 
    
   _d2halpha_dphigamma_dphiepsilon[_qp] = 
    -2.0*(_dhalpha_dphiepsilon[_qp]* _phase_gamma[_qp] + _dhalpha_dphigamma[_qp]*_phase_epsilon[_qp]) *_getsumSquare;  
        
   // For h_{beta,gamma}
   _d2hbeta_dphigamma_dphialpha[_qp] = 
     -2.0*(_dhbeta_dphialpha[_qp] * _phase_gamma[_qp] + _dhbeta_dphigamma[_qp] * _phase_alpha[_qp]) * _getsumSquare;
       
   _d2hbeta_dphigamma_dphibeta[_qp] = 
    -2.0*(_dhbeta_dphigamma[_qp] * _phase_beta[_qp] + _dhbeta_dphibeta[_qp] * _phase_gamma[_qp]) * _getsumSquare;
    
   _d2hbeta_dphigamma_dphidelta[_qp] = 
    -2.0*(_dhbeta_dphidelta[_qp] * _phase_gamma[_qp] + _dhbeta_dphigamma[_qp] * _phase_delta[_qp]) * _getsumSquare;
    
   _d2hbeta_dphigamma_dphiepsilon[_qp] = 
    -2.0*(_dhbeta_dphiepsilon[_qp] * _phase_gamma[_qp] + _dhbeta_dphigamma[_qp] * _phase_epsilon[_qp]) * _getsumSquare; 
    
     // For h_{delta,gamma}
   _d2hdelta_dphigamma_dphialpha[_qp] = 
     -2.0*(_dhdelta_dphialpha[_qp] * _phase_gamma[_qp] + _dhdelta_dphigamma[_qp] * _phase_alpha[_qp]) * _getsumSquare;
       
   _d2hdelta_dphigamma_dphibeta[_qp] = 
    -2.0*(_dhdelta_dphibeta[_qp] * _phase_gamma[_qp] + _dhdelta_dphigamma[_qp] * _phase_beta[_qp]) * _getsumSquare;
    
   _d2hdelta_dphigamma_dphidelta[_qp] = 
    -2.0*(_dhdelta_dphigamma[_qp] * _phase_delta[_qp] + _dhdelta_dphidelta[_qp] * _phase_gamma[_qp]) * _getsumSquare;
    
   _d2hdelta_dphigamma_dphiepsilon[_qp] = 
    -2.0*(_dhdelta_dphiepsilon[_qp] * _phase_gamma[_qp] + _dhdelta_dphigamma[_qp] * _phase_epsilon[_qp]) * _getsumSquare;
    
    // For h_{epsilon,gamma} 
    _d2hepsilon_dphigamma_dphialpha[_qp] = 
    -2.0*(_dhepsilon_dphialpha[_qp] * _phase_gamma[_qp] + _dhepsilon_dphigamma[_qp] * _phase_alpha[_qp]) * _getsumSquare;
   
    _d2hepsilon_dphigamma_dphibeta[_qp] = 
    -2.0*(_dhepsilon_dphibeta[_qp] * _phase_gamma[_qp] + _dhepsilon_dphigamma[_qp] * _phase_beta[_qp]) * _getsumSquare;
    
    _d2hepsilon_dphigamma_dphidelta[_qp] = 
    -2.0*(_dhepsilon_dphidelta[_qp] * _phase_gamma[_qp] + _dhepsilon_dphigamma[_qp] * _phase_delta[_qp]) * _getsumSquare;
    
    _d2hepsilon_dphigamma_dphiepsilon[_qp] = 
    -2.0*(_dhepsilon_dphigamma[_qp] * _phase_epsilon[_qp] + _dhepsilon_dphiepsilon[_qp] * _phase_gamma[_qp])* _getsumSquare;
    
   //****************************************For phase delta********************************************************//
    
   // For h_{alpha,delta}
   
    _d2halpha_dphidelta_dphialpha[_qp] = 
     -2.0*(_dhalpha_dphidelta[_qp] * _phase_alpha[_qp] + _dhalpha_dphialpha[_qp] * _phase_delta[_qp]) * _getsumSquare;
       
   _d2halpha_dphidelta_dphibeta[_qp] = 
    -2.0*(_dhalpha_dphibeta[_qp] * _phase_delta[_qp] + _dhalpha_dphidelta[_qp] * _phase_beta[_qp]) * _getsumSquare;
    
   _d2halpha_dphidelta_dphigamma[_qp] = 
    -2.0*(_dhalpha_dphigamma[_qp] * _phase_delta[_qp] + _dhalpha_dphidelta[_qp] * _phase_gamma[_qp]) * _getsumSquare; 
    
   _d2halpha_dphidelta_dphiepsilon[_qp] = 
    -2.0*(_dhalpha_dphiepsilon[_qp] * _phase_delta[_qp] + _dhalpha_dphidelta[_qp] * _phase_epsilon[_qp]) * _getsumSquare;  
     
    
    // For h_{beta,delta}
    
   _d2hbeta_dphidelta_dphialpha[_qp] = 
     -2.0*(_dhbeta_dphialpha[_qp] * _phase_delta[_qp] + _dhbeta_dphidelta[_qp] * _phase_alpha[_qp]) * _getsumSquare;
       
   _d2hbeta_dphidelta_dphibeta[_qp] = 
    -2.0*(_dhbeta_dphidelta[_qp] * _phase_beta[_qp] + _dhbeta_dphibeta[_qp] * _phase_delta[_qp]) * _getsumSquare;
    
   _d2hbeta_dphidelta_dphigamma[_qp] = 
    -2.0*(_dhbeta_dphigamma[_qp] * _phase_delta[_qp] + _dhbeta_dphidelta[_qp] * _phase_gamma[_qp]) * _getsumSquare;
    
    _d2hbeta_dphidelta_dphiepsilon[_qp] = 
    -2.0*(_dhbeta_dphiepsilon[_qp] * _phase_delta[_qp] + _dhbeta_dphidelta[_qp] * _phase_epsilon[_qp]) * _getsumSquare;  
    
    // For h_{gamma,delta}
    
   _d2hgamma_dphidelta_dphialpha[_qp] = 
     -2.0*(_dhgamma_dphialpha[_qp] * _phase_delta[_qp] + _dhgamma_dphidelta[_qp] * _phase_alpha[_qp]) * _getsumSquare;
       
   _d2hgamma_dphidelta_dphibeta[_qp] = 
    -2.0*(_dhgamma_dphibeta[_qp] * _phase_delta[_qp] + _dhgamma_dphidelta[_qp] * _phase_beta[_qp]) * _getsumSquare;
    
   _d2hgamma_dphidelta_dphigamma[_qp] = 
    -2.0*(_dhgamma_dphidelta[_qp] * _phase_gamma[_qp] + _dhgamma_dphigamma[_qp] * _phase_delta[_qp]) * _getsumSquare;
    
     _d2hgamma_dphidelta_dphiepsilon[_qp] = 
    -2.0*(_dhgamma_dphiepsilon[_qp] * _phase_delta[_qp] + _dhgamma_dphidelta[_qp] * _phase_epsilon[_qp]) * _getsumSquare;
    
    // For h_{epsilon,delta}
    
   _d2hepsilon_dphidelta_dphialpha[_qp] = 
     -2.0*(_dhepsilon_dphialpha[_qp] * _phase_delta[_qp] + _dhepsilon_dphidelta[_qp] * _phase_alpha[_qp]) * _getsumSquare;
       
   _d2hepsilon_dphidelta_dphibeta[_qp] = 
    -2.0*(_dhepsilon_dphibeta[_qp] * _phase_delta[_qp] + _dhepsilon_dphidelta[_qp] * _phase_beta[_qp]) * _getsumSquare;
    
   _d2hepsilon_dphidelta_dphigamma[_qp] = 
    -2.0*(_dhepsilon_dphigamma[_qp] * _phase_delta[_qp] + _dhepsilon_dphidelta[_qp] * _phase_gamma[_qp]) * _getsumSquare;
    
     _d2hepsilon_dphidelta_dphiepsilon[_qp] = 
    -2.0*(_dhepsilon_dphidelta[_qp] * _phase_epsilon[_qp] + _dhepsilon_dphiepsilon[_qp] * _phase_delta[_qp]) * _getsumSquare;
    
    //****************************************For phase epsilon********************************************************//
    
    // For h_{alpha,epsilon}
   
    _d2halpha_dphiepsilon_dphialpha[_qp] = 
     -2.0*(_dhalpha_dphiepsilon[_qp] * _phase_alpha[_qp] + _dhalpha_dphialpha[_qp] * _phase_epsilon[_qp]) * _getsumSquare;
       
   _d2halpha_dphiepsilon_dphibeta[_qp] = 
    -2.0*(_dhalpha_dphibeta[_qp] * _phase_epsilon[_qp] + _dhalpha_dphiepsilon[_qp] * _phase_beta[_qp]) *    _getsumSquare;
    
   _d2halpha_dphiepsilon_dphigamma[_qp] = 
    -2.0*(_dhalpha_dphigamma[_qp] * _phase_epsilon[_qp] + _dhalpha_dphiepsilon[_qp] * _phase_gamma[_qp]) * _getsumSquare; 
    
   _d2halpha_dphiepsilon_dphidelta[_qp] = 
    -2.0*(_dhalpha_dphidelta[_qp] * _phase_epsilon[_qp] + _dhalpha_dphiepsilon[_qp] * _phase_delta[_qp]) * _getsumSquare;
    
    // For h_{beta,epsilon}
   
    _d2hbeta_dphiepsilon_dphialpha[_qp] = 
     -2.0*(_dhbeta_dphialpha[_qp] * _phase_epsilon[_qp] + _dhbeta_dphiepsilon[_qp] * _phase_alpha[_qp]) * _getsumSquare;
       
   _d2hbeta_dphiepsilon_dphibeta[_qp] = 
    -2.0*(_dhbeta_dphiepsilon[_qp] * _phase_beta[_qp] + _dhbeta_dphibeta[_qp] * _phase_epsilon[_qp]) *    _getsumSquare;
    
   _d2hbeta_dphiepsilon_dphigamma[_qp] = 
    -2.0*(_dhbeta_dphigamma[_qp] * _phase_epsilon[_qp] + _dhbeta_dphiepsilon[_qp] * _phase_gamma[_qp]) * _getsumSquare; 
    
   _d2hbeta_dphiepsilon_dphidelta[_qp] = 
    -2.0*(_dhbeta_dphidelta[_qp] * _phase_epsilon[_qp] + _dhbeta_dphiepsilon[_qp] * _phase_delta[_qp]) * _getsumSquare;
    
    // For h_{gamma,epsilon}
    
   _d2hgamma_dphiepsilon_dphialpha[_qp] = 
     -2.0*(_dhgamma_dphialpha[_qp] * _phase_epsilon[_qp] + _dhgamma_dphiepsilon[_qp] * _phase_alpha[_qp]) * _getsumSquare;
       
   _d2hgamma_dphiepsilon_dphibeta[_qp] = 
    -2.0*(_dhgamma_dphibeta[_qp] * _phase_epsilon[_qp] + _dhgamma_dphiepsilon[_qp] * _phase_beta[_qp]) * _getsumSquare;
    
   _d2hgamma_dphiepsilon_dphigamma[_qp] = 
    -2.0*(_dhgamma_dphiepsilon[_qp] * _phase_gamma[_qp] + _dhgamma_dphigamma[_qp] * _phase_epsilon[_qp]) * _getsumSquare;
    
     _d2hgamma_dphiepsilon_dphidelta[_qp] = 
    -2.0*(_dhgamma_dphidelta[_qp] * _phase_epsilon[_qp] + _dhgamma_dphiepsilon[_qp] * _phase_delta[_qp]) * _getsumSquare; 
    
     // For h_{delta,epsilon}
    
   _d2hdelta_dphiepsilon_dphialpha[_qp] = 
     -2.0*(_dhdelta_dphialpha[_qp] * _phase_epsilon[_qp] + _dhdelta_dphiepsilon[_qp] * _phase_alpha[_qp]) * _getsumSquare;
       
   _d2hdelta_dphiepsilon_dphibeta[_qp] = 
    -2.0*(_dhdelta_dphibeta[_qp] * _phase_epsilon[_qp] + _dhdelta_dphiepsilon[_qp] * _phase_beta[_qp]) * _getsumSquare;
    
   _d2hdelta_dphiepsilon_dphigamma[_qp] = 
    -2.0*(_dhdelta_dphigamma[_qp] * _phase_epsilon[_qp] + _dhdelta_dphiepsilon[_qp] * _phase_gamma[_qp]) * _getsumSquare;
    
     _d2hdelta_dphiepsilon_dphidelta[_qp] = 
    -2.0*(_dhdelta_dphiepsilon[_qp] * _phase_delta[_qp] + _dhdelta_dphidelta[_qp] * _phase_epsilon[_qp]) * _getsumSquare;   
      
}
