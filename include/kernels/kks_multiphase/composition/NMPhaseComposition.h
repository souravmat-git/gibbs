//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef NMPHASECOMPOSITION_H
#define NMPHASECOMPOSITION_H

#include "Kernel.h"

// Forward declaration
class NMPhaseComposition;

//class template specialization:
template <>
InputParameters validParams<NMPhaseComposition>();

/**
 * This class enforces the equation for mole fraction
 * c = c_{alpha}* h_alpha + c_{beta}* h_beta + c_{gamma} * h_gamma
 * The class is a derived class from base class kernel
 **/
class NMPhaseComposition : public Kernel
{
    public:
        NMPhaseComposition(const InputParameters & parameters);

    // Member functions declarations 
    protected:
    
        virtual Real computeQpResidual() override;
        virtual Real computeQpJacobian() override;
        virtual Real computeQpOffDiagJacobian(unsigned jvar) override;
        
    private:
     
     //Note: typedef LIBMESH_DEFAULT_SCALAR_TYPE Real    
        const VariableValue & _phase_comp_beta;
        unsigned int _phase_comp_beta_var;
    
        const VariableValue & _phase_comp_alpha;
        unsigned int _phase_comp_alpha_var;
    
        const VariableValue & _phase_alpha;
        unsigned int _phase_alpha_var;
    
        const VariableValue & _phase_beta;
        unsigned int _phase_beta_var;
    
        const VariableValue & _phase_gamma;
        unsigned int _phase_gamma_var;
    
        const VariableValue & _mole_fraction;
        unsigned int _mole_fraction_var;
    

        const MaterialProperty<Real> & _h_alpha;
        const MaterialProperty<Real> & _h_beta;
        const MaterialProperty<Real> & _h_gamma;    
        
        const MaterialProperty<Real> & _dhalpha_dphialpha;
        const MaterialProperty<Real> & _dhbeta_dphibeta;
        const MaterialProperty<Real> & _dhgamma_dphigamma;
        
        const MaterialProperty<Real> & _dhbeta_dphialpha;
        const MaterialProperty<Real> & _dhgamma_dphialpha;
        
        const MaterialProperty<Real> & _dhalpha_dphibeta;
        const MaterialProperty<Real> & _dhgamma_dphibeta;
        
        const MaterialProperty<Real> & _dhalpha_dphigamma;
        const MaterialProperty<Real> & _dhbeta_dphigamma;
};
#endif // NMPHASECOMPOSITION_H
