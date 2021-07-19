#include "DiffusionFlux.h"

registerMooseObject("gibbsApp",DiffusionFlux);

template<>
InputParameters
validParams<DiffusionFlux>()
{
  InputParameters params = validParams<AuxKernel>();
  /**
    * Declare the options for a MooseEnum
    * Since the flux is a vector
    */
  MooseEnum component("x y z");
  
  //Use the MooseEnum to add a parameter called "component"
  params.addRequiredParam<MooseEnum>("component",component,"The desired component or direction.");

  //Add a "coupling field variable" to get a variable from input file
  params.addRequiredCoupledVar("concentration", "Concentration field");
  params.addRequiredParam <Real>("diffusivity","Diffusion coefficient (D)");

return params;
}

DiffusionFlux::DiffusionFlux(const InputParameters & parameters)
  :AuxKernel(parameters),
  
  //This will automatically convert the MOOSEnum into an integer
  _component(getParam<MooseEnum>("component")),

  //Get the gradient of the variable
  _conc_gradient(coupledGradient("concentration")),
  _diffusivity(getParam<Real>("diffusivity"))
{

}
 
Real
DiffusionFlux::computeValue()
{
  //Access the gradient of the concentration
  // Then pull out the componnet of it we are looking for (x,y,z)
  //Note that getting a particular component of a gradient is done using the
  //paranthesis operator

  return -(_diffusivity)*(_conc_gradient[_qp])(_component);

}

