//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "AllenCahnElasticEnergyOffDiag.h"

class QDrivingForce3DOffDiag;

template<>
InputParameters validParams<QDrivingForce3DOffDiag>();

/**
  *This class enforces the following 
  *Equation in the LC model
  *h^{\prime}(fc_beta - fc_alpha)= 0
  **/

class QDrivingForce3DOffDiag :  public AllenCahnElasticEnergyOffDiag
{
  public: 
    QDrivingForce3DOffDiag(const InputParameters & parameters);
  
  protected:
  
  Real computeQpResidual() override { return 0.0; }
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  private:

  unsigned int _ndisp;  
  const MaterialProperty<Real> & _nd_factor;
};
