//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//Include guard
#ifndef QUANTINTERPOLATIONFUNCTION_H
#define QUANTINTERPOLATIONFUNCTION_H

//Forward declaration
class QuantInterpolationFunction;

//Included dependencies
#include "Material.h"

template <>
InputParameters validParams<QuantInterpolationFunction>();

//Derived class: QuantInterpolationFunction
//Base class: Material.C

class QuantInterpolationFunction : public Material
{
  public:
    QuantInterpolationFunction(const InputParameters & parameters);
  
  protected:
    virtual void computeQpProperties() override;
  
  private:  
  
   /**
   * Define the member references that will hold the computed values
   * for the Real value properties in this class.
   * _h_alpha : holds the interpolation function for alpha phase
   *_h_beta : holds the interpolation function for the beta phase
   *_h_gamma : holds the inetrpolation function for the gamma phase
   */
    MaterialProperty<Real> & _h_alpha;
    MaterialProperty<Real> & _h_beta;
    MaterialProperty<Real> & _h_gamma;
    MaterialProperty<Real> & _h_delta;
    MaterialProperty<Real> & _h_epsilon;
    
    //Diagonal terms of the first derivatives.
    // Note that we do not require them in the residual
    // only the jacobian calculation we require them
    
    MaterialProperty<Real> & _dhalpha_dphialpha;
    MaterialProperty<Real> & _dhbeta_dphibeta;
    MaterialProperty<Real> & _dhgamma_dphigamma;
    MaterialProperty<Real> & _dhdelta_dphidelta;
    MaterialProperty<Real> & _dhepsilon_dphiepsilon;
    
    //Required in residual of kernels operating on phase_{alpha}
    MaterialProperty<Real> & _dhbeta_dphialpha;
    MaterialProperty<Real> & _dhgamma_dphialpha;
    MaterialProperty<Real> & _dhdelta_dphialpha;
    MaterialProperty<Real> & _dhepsilon_dphialpha;
    
    //Required in residual of kernels operating on phase_{beta}
    MaterialProperty<Real> & _dhalpha_dphibeta;
    MaterialProperty<Real> & _dhgamma_dphibeta;
    MaterialProperty<Real> & _dhdelta_dphibeta;
    MaterialProperty<Real> & _dhepsilon_dphibeta;
    
   //Required in residual of kernels operating on phase_{gamma}
    MaterialProperty<Real> & _dhalpha_dphigamma;
    MaterialProperty<Real> & _dhbeta_dphigamma;
    MaterialProperty<Real> & _dhdelta_dphigamma;
    MaterialProperty<Real> & _dhepsilon_dphigamma;
    
   //Required in residual of kernels operating on phase_{delta} 
    MaterialProperty<Real> & _dhalpha_dphidelta;
    MaterialProperty<Real> & _dhbeta_dphidelta;
    MaterialProperty<Real> & _dhgamma_dphidelta;
    MaterialProperty<Real> & _dhepsilon_dphidelta;
    
    //Required in residual of kernels operating on phase_{epsilon} 
    MaterialProperty<Real> & _dhalpha_dphiepsilon;
    MaterialProperty<Real> & _dhbeta_dphiepsilon;
    MaterialProperty<Real> & _dhgamma_dphiepsilon;
    MaterialProperty<Real> & _dhdelta_dphiepsilon;
    
    //***** N*N**** first derivatives, where N = number of phases//
    
    // Required in Jacobian of kernels operating on phase_{alpha}
    MaterialProperty<Real> & _d2hbeta_dphialpha2;
    MaterialProperty<Real> & _d2hgamma_dphialpha2;
    MaterialProperty<Real> & _d2hdelta_dphialpha2;
    MaterialProperty<Real> & _d2hepsilon_dphialpha2;
    
    // Required in Jacobian of kernels operating on phase_{beta}
    MaterialProperty<Real> & _d2halpha_dphibeta2;
    MaterialProperty<Real> & _d2hgamma_dphibeta2;
    MaterialProperty<Real> & _d2hdelta_dphibeta2;
    MaterialProperty<Real> & _d2hepsilon_dphibeta2;
    
    // Required in Jacobian of kernels operating on phase_{gamma}
    MaterialProperty<Real> & _d2halpha_dphigamma2;
    MaterialProperty<Real> & _d2hbeta_dphigamma2;
    MaterialProperty<Real> & _d2hdelta_dphigamma2;
    MaterialProperty<Real> & _d2hepsilon_dphigamma2;
    
    // Required in jacobians of kernels operating on phase_{delta}   
    MaterialProperty<Real> & _d2halpha_dphidelta2;
    MaterialProperty<Real> & _d2hbeta_dphidelta2;
    MaterialProperty<Real> & _d2hgamma_dphidelta2;
    MaterialProperty<Real> & _d2hepsilon_dphidelta2;
    
    // Required in jacobians of kernels operating on phase_{epsilon}   
    MaterialProperty<Real> & _d2halpha_dphiepsilon2;
    MaterialProperty<Real> & _d2hbeta_dphiepsilon2;
    MaterialProperty<Real> & _d2hgamma_dphiepsilon2;
    MaterialProperty<Real> & _d2hdelta_dphiepsilon2;
    
    // Required in off-diag kernels operating on  phase_{alpha}
    // since the coupled terms are phase_{beta}, phase_{gamma}, _phase_{delta}
    
    MaterialProperty<Real> & _d2hbeta_dphialpha_dphibeta;
    MaterialProperty<Real> & _d2hbeta_dphialpha_dphigamma;
    MaterialProperty<Real> & _d2hbeta_dphialpha_dphidelta;
    MaterialProperty<Real> & _d2hbeta_dphialpha_dphiepsilon;
    
    MaterialProperty<Real> & _d2hgamma_dphialpha_dphibeta;
    MaterialProperty<Real> & _d2hgamma_dphialpha_dphigamma;
    MaterialProperty<Real> & _d2hgamma_dphialpha_dphidelta;
    MaterialProperty<Real> & _d2hgamma_dphialpha_dphiepsilon;
    
    MaterialProperty<Real> & _d2hdelta_dphialpha_dphibeta;
    MaterialProperty<Real> & _d2hdelta_dphialpha_dphigamma;
    MaterialProperty<Real> & _d2hdelta_dphialpha_dphidelta;
    MaterialProperty<Real> & _d2hdelta_dphialpha_dphiepsilon;
    
    MaterialProperty<Real> & _d2hepsilon_dphialpha_dphibeta;
    MaterialProperty<Real> & _d2hepsilon_dphialpha_dphigamma;
    MaterialProperty<Real> & _d2hepsilon_dphialpha_dphidelta;
    MaterialProperty<Real> & _d2hepsilon_dphialpha_dphiepsilon;
    
    // Required by kernels operating on phase_{beta}
    // since the coupled terms are phase_{alpha}, phase_{gamma}
    
    MaterialProperty<Real> & _d2halpha_dphibeta_dphialpha;
    MaterialProperty<Real> & _d2halpha_dphibeta_dphigamma;
    MaterialProperty<Real> & _d2halpha_dphibeta_dphidelta;
    MaterialProperty<Real> & _d2halpha_dphibeta_dphiepsilon;
    
    MaterialProperty<Real> & _d2hgamma_dphibeta_dphialpha;
    MaterialProperty<Real> & _d2hgamma_dphibeta_dphigamma;
    MaterialProperty<Real> & _d2hgamma_dphibeta_dphidelta;
    MaterialProperty<Real> & _d2hgamma_dphibeta_dphiepsilon;
    
    MaterialProperty<Real> & _d2hdelta_dphibeta_dphialpha;
    MaterialProperty<Real> & _d2hdelta_dphibeta_dphigamma;
    MaterialProperty<Real> & _d2hdelta_dphibeta_dphidelta;
    MaterialProperty<Real> & _d2hdelta_dphibeta_dphiepsilon;
    
    MaterialProperty<Real> & _d2hepsilon_dphibeta_dphialpha;
    MaterialProperty<Real> & _d2hepsilon_dphibeta_dphigamma;
    MaterialProperty<Real> & _d2hepsilon_dphibeta_dphidelta;
    MaterialProperty<Real> & _d2hepsilon_dphibeta_dphiepsilon;
    
    // Required by kernels operating on phase_{gamma}
    //since the coupled terms are phase_{beta}, phase_{alpha}
    
    MaterialProperty<Real> & _d2halpha_dphigamma_dphibeta;
    MaterialProperty<Real> & _d2halpha_dphigamma_dphialpha;
    MaterialProperty<Real> & _d2halpha_dphigamma_dphidelta;
    MaterialProperty<Real> & _d2halpha_dphigamma_dphiepsilon;
    
    MaterialProperty<Real> & _d2hbeta_dphigamma_dphialpha;
    MaterialProperty<Real> & _d2hbeta_dphigamma_dphibeta;
    MaterialProperty<Real> & _d2hbeta_dphigamma_dphidelta;
    MaterialProperty<Real> & _d2hbeta_dphigamma_dphiepsilon;
    
    MaterialProperty<Real> & _d2hdelta_dphigamma_dphialpha;
    MaterialProperty<Real> & _d2hdelta_dphigamma_dphibeta;
    MaterialProperty<Real> & _d2hdelta_dphigamma_dphidelta;
    MaterialProperty<Real> & _d2hdelta_dphigamma_dphiepsilon;
    
    MaterialProperty<Real> & _d2hepsilon_dphigamma_dphialpha;
    MaterialProperty<Real> & _d2hepsilon_dphigamma_dphibeta;
    MaterialProperty<Real> & _d2hepsilon_dphigamma_dphidelta;
    MaterialProperty<Real> & _d2hepsilon_dphigamma_dphiepsilon;
    
    // Required by kernels operating on phase_{delta}
    //since the coupled terms are phase_{alpha}, phase_{beta}, phase_{gamma}
    
    MaterialProperty<Real> & _d2halpha_dphidelta_dphialpha;
    MaterialProperty<Real> & _d2halpha_dphidelta_dphibeta;
    MaterialProperty<Real> & _d2halpha_dphidelta_dphigamma;
    MaterialProperty<Real> & _d2halpha_dphidelta_dphiepsilon;
    
    MaterialProperty<Real> & _d2hbeta_dphidelta_dphialpha;
    MaterialProperty<Real> & _d2hbeta_dphidelta_dphibeta;
    MaterialProperty<Real> & _d2hbeta_dphidelta_dphigamma;
    MaterialProperty<Real> & _d2hbeta_dphidelta_dphiepsilon;
    
    MaterialProperty<Real> & _d2hgamma_dphidelta_dphialpha;
    MaterialProperty<Real> & _d2hgamma_dphidelta_dphibeta;
    MaterialProperty<Real> & _d2hgamma_dphidelta_dphigamma;
    MaterialProperty<Real> & _d2hgamma_dphidelta_dphiepsilon;  
    
    MaterialProperty<Real> & _d2hepsilon_dphidelta_dphialpha;
    MaterialProperty<Real> & _d2hepsilon_dphidelta_dphibeta;
    MaterialProperty<Real> & _d2hepsilon_dphidelta_dphigamma;
    MaterialProperty<Real> & _d2hepsilon_dphidelta_dphiepsilon;
    
    // Required by kernels operating on phase_{epsilon}
    //since the coupled terms are phase_{alpha}, phase_{beta}, phase_{delta} 
    
    MaterialProperty<Real> & _d2halpha_dphiepsilon_dphialpha;
    MaterialProperty<Real> & _d2halpha_dphiepsilon_dphibeta;
    MaterialProperty<Real> & _d2halpha_dphiepsilon_dphigamma;
    MaterialProperty<Real> & _d2halpha_dphiepsilon_dphidelta;
    
    MaterialProperty<Real> & _d2hbeta_dphiepsilon_dphialpha;
    MaterialProperty<Real> & _d2hbeta_dphiepsilon_dphibeta;
    MaterialProperty<Real> & _d2hbeta_dphiepsilon_dphigamma;
    MaterialProperty<Real> & _d2hbeta_dphiepsilon_dphidelta;
    
    MaterialProperty<Real> & _d2hgamma_dphiepsilon_dphialpha;
    MaterialProperty<Real> & _d2hgamma_dphiepsilon_dphibeta;
    MaterialProperty<Real> & _d2hgamma_dphiepsilon_dphigamma;
    MaterialProperty<Real> & _d2hgamma_dphiepsilon_dphidelta;  
    
    MaterialProperty<Real> & _d2hdelta_dphiepsilon_dphialpha;
    MaterialProperty<Real> & _d2hdelta_dphiepsilon_dphibeta;
    MaterialProperty<Real> & _d2hdelta_dphiepsilon_dphigamma;
    MaterialProperty<Real> & _d2hdelta_dphiepsilon_dphidelta;  
            
    const VariableValue & _phase_alpha;
    const VariableValue & _phase_beta;
    const VariableValue & _phase_gamma;
    const VariableValue & _phase_delta; 
    const VariableValue & _phase_epsilon;   

};
#endif //QUANTINTERPOLATIONFUNCTION_H
