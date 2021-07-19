#include "GrandPotential.h"

registerMooseObject("gibbsApp",GrandPotential);

template<>
InputParameters
validParams<GrandPotential>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("xB","Component B in alpha phase");
  params.addCoupledVar("xC", 0.0,"Component C in alpha phase");
  params.addParam<MaterialPropertyName>("f",0.0, "Free energy of a phase");
  params.addParam<MaterialPropertyName>("B_diff_pot",0.0, "Diffusion potential of B componet");
  params.addParam<MaterialPropertyName>("C_diff_pot", 0.0,"Diffusion potential of C component");
  return params;
}

GrandPotential::GrandPotential(const InputParameters & parameters)
  :AuxKernel(parameters),
  _xB(coupledValue("xB")),
  _xC(coupledValue("xC")),
  _f(getMaterialProperty<Real>("f")),
  _B_diff_pot(getMaterialProperty<Real>("B_diff_pot")),
  _C_diff_pot(getMaterialProperty<Real>("C_diff_pot"))
{
}
 
Real
GrandPotential::computeValue()
{
                            
  
  return (_f[_qp] - _B_diff_pot[_qp] * _xB[_qp] - _C_diff_pot[_qp] * _xC[_qp]);
}
