#This kernel uses the TabulatedConjugatePhaseData
#Material to solve the equations. 
#This is a model Ni-Al alloy problem
#beta phase - Solid (fcc)
#alpha phase- Liquid
#T = 1000K
#Simulation details:
#Interfacial energy = 0.5 J/m2
#interfacial width = 0.9e-6 m
#Length = 180e-6 m
#characteristic_length = (0.18e-6)/6 m = 0.03e-6 m
#charcteristic time = 1.5404e-7s 
#Nondimensional length = (Length/characteristic_length)
#No. of finite elements = 6.0(Length/interface_width)

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx  = 1200
  xmin = -3000
  xmax =  3000
  elem_type = EDGE2
[]


[Variables]

 [./phi_alpha]
  order = FIRST
  family = LAGRANGE
 [../]

 [./phi_beta]
  order = FIRST
  family = LAGRANGE
 [../]

 [./x_Al]
  order = FIRST
  family = LAGRANGE
 [../]

 [./mu_Al]
  order = FIRST
  family = LAGRANGE
 [../]

[]

[BCs]

  [./alpha_pf]
    type = NeumannBC
    variable = 'phi_alpha'
    boundary = 'left right'
    value = 0
  [../]
  
  [./beta_pf]
    type = NeumannBC
    variable = 'phi_beta'
    boundary = 'left right'
    value = 0
  [../]
  
   [./MassFluxBC]
    type = NeumannBC
    variable = 'mu_Al'
    boundary = 'left right'
    value = 0
  [../]
 
[]


[ICs]
  
   #beta = solid (fcc)
   [./ic_phi_beta]
    type = BoundingBoxIC
    variable = phi_beta
    x1 = -3000
    y1 = 0.0
    x2 = 0.0
    y2 = 0.0
    inside = 1.0
    outside = 0.0
  [../]

  #alpha = liquid
  [./ic_phi_alpha]
    type = BoundingBoxIC
    variable = phi_alpha
    x1 = 0.0
    y1 = 0.0
    x2 = 3000
    y2 = 0.0
    inside = 1.0
    outside = 0.0
  [../]
   
  #inside = solid
  #outside = liquid

  [./ic_func_xAl]
   type = BoundingBoxIC
   variable = x_Al
   x1 = -3000
   y1 = 0.0
   x2 = 0.0
   y2 = 0.0
   inside =  4.05E-1       
   outside = 4.73E-1        
  [../]
  
  #[./muAl_IC]
  # variable = mu_Al
  # type = ConstantIC
  # value = -2.8145E0
  #[../]
  
[]


[Kernels]
  
  #Variable that this kernel operates on
  #is the diffusion potential of comp B
  
  [./PhaseConc]
    type = MultiPhaseConstraintMu
    variable = mu_Al
    xB = x_Al
    phase_alpha = phi_alpha
    phase_beta = phi_beta
  [../]

  ######################
  # Diffusion Equation #
  ######################
    
  [./Al_balance]
    type = BinaryMultiPhaseMassBalance
    variable = x_Al
    B_diff_pot = mu_Al
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    h_alpha = h_alpha
    h_beta = h_beta
  [../]

  [./xAl_dot]
    type = TimeDerivative
    variable = x_Al
  [../]
  
  #######################
  # Allen-Cahn equation #
  # for phase-alpha     #
  #######################
  
  [./BinaryMultiPhaseDrivingForce_alpha]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_alpha     #Phase_1
    phase_2 = phi_beta       #Coupled phase_2
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_alpha
    A_chem_pot_2 = A_chem_pot_beta
    xB_1 = xB_alpha
    xB_2 = xB_beta
    dh = dhbeta_dphialpha     #Required by the residual
    d2h = d2hbeta_dphialpha2  #Required by the diag jacobian
    d2h_2 = d2hbeta_dphialpha_dphibeta
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
  #This kernel implements the double well
  #The form of the double-well barrier is 
  #f(phi) = phi^(2)(1 - phi^(2))
  #The free energy is explicitly coded
  #within the kernel
  
   [./doublewell_phi_alpha]
    type = MultiPhaseDoubleWell
    variable = phi_alpha
    eta2 = phi_beta
    gamma = 1.5
    mob_name = L_phi
    barrier_height = m
  [../]

  [./Curvature_alpha]
    type = ACInterface
    variable = phi_alpha
    kappa_name = kappa
    mob_name = L_phi
  [../]
  
  [./dalpha_dt]
    type = TimeDerivative
    variable = phi_alpha
  [../]
  
  #######################
  # Allen-Cahn equation #
  # for phase-beta      #
  #######################
   
  [./BinaryMultiPhaseDrivingForce_beta]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_beta
    phase_2 = phi_alpha
    B_diff_pot = mu_Al
    A_chem_pot_1 = A_chem_pot_beta
    A_chem_pot_2 = A_chem_pot_alpha
    xB_1 = xB_beta
    xB_2 = xB_alpha
    dh = dhalpha_dphibeta
    d2h = d2halpha_dphibeta2
    d2h_2 = d2halpha_dphibeta_dphialpha
    mob_name = L_phi
    nd_factor = nd_factor
  [../]
  
  #This kernel implements the double well
  #The form of the double-well barrier is 
  #f(phi) = phi^(2)(1 - phi^(2))
  #The free energy is explicitly coded
  #within the kernel
  
   [./doublewell_phi_beta]
    type = MultiPhaseDoubleWell
    variable = phi_beta
    eta2 = phi_alpha
    gamma = 1.5
    mob_name = L_phi
    barrier_height = m
  [../]

  [./Curvature_beta]
    type = ACInterface
    variable = phi_beta
    kappa_name = kappa
    mob_name = L_phi
  [../]
  
  [./dphi_beta_dt]
    type = TimeDerivative
    variable = phi_beta
  [../]
  
[]


[Materials]
  ##Phenomenological coeff###
  
   [./ConstantFieldProperties]
     type = GenericConstantMaterial
     prop_names =  'L_phi    kappa    m    nd_factor'
     prop_values = '0.0011   112.5    1.0  33.256'
  [../]  
  
   [./ZeroFieldProperties]
     type = GenericConstantMaterial
     prop_names =  'L_BB_gamma L_BB_delta L_BB_epsilon 
                   dL_BB_muB_gamma dL_BB_muB_delta dL_BB_muB_epsilon'
     prop_values = '0 0 0 0 0 0'
  [../] 
  
  [./InterpolationFunction]
    type = QuantInterpolationFunction
    phase_alpha = phi_alpha
    phase_beta = phi_beta
  [../]  
  
 [./BinaryConjugatePhaseMaterial_FCC]
    type = BinaryConjugatePhaseMaterial
    A_chem_pot = A_chem_pot_beta
    xB = xB_beta
    inv_B_tf = inv_B_tf_beta
    B_diff_pot = mu_Al
    molar_volume = 1.0
    table_object = table_FCC_data
    #outputs = exodus
  [../]
  
  [./BinaryConjugatePhaseMaterial_Liquid]
    type = BinaryConjugatePhaseMaterial
    A_chem_pot = A_chem_pot_alpha
    xB = xB_alpha
    inv_B_tf = inv_B_tf_alpha
    B_diff_pot = mu_Al
    molar_volume = 1.0
    table_object = table_Liquid_data
    #outputs = exodus
  [../]
  
  [./BinaryConjugateKineticMaterial_FCC]
    type = BinaryConjugateKineticMaterial
    L_BB = L_BB_beta
    dL_BB_muB = dL_BB_muB_beta
    B_diff_pot = mu_Al
    table_object = table_FCC_mobility_data
    #outputs = exodus
  [../]
  
  [./BinaryConjugateKineticMaterial_Liquid]
    type = BinaryConjugateKineticMaterial
    L_BB = L_BB_alpha
    dL_BB_muB = dL_BB_muB_alpha
    B_diff_pot = mu_Al
    table_object = table_Liquid_mobility_data
    #outputs = exodus
  [../]

  [./L_BB]
    type = ParsedMaterial
    f_name = L_BB
    material_property_names = 'L_BB_alpha L_BB_beta h_alpha h_beta'
    function = 'h_alpha*L_BB_alpha + h_beta*L_BB_beta'
  [../]
 
[]

[UserObjects]

  #This userobjcts reads the data from the table 
  #for a given diffusion potential
  #and provides the interpolated grand-potential,
  #mole fraction and inverse of the thermodynamic_factor
  #to the material TabulatedPhaseMaterial
  
  [./table_FCC_data]
    type = BinaryConjugatePhaseData
    table_name = NiAl_data/NiAl_invert_data_FCC.csv
  [../]
  
  [./table_Liquid_data]
    type = BinaryConjugatePhaseData
    table_name = NiAl_data/NiAl_invert_data_Liquid.csv
  [../] 
  
  [./table_FCC_mobility_data]
    type = BinaryConjugateMobilityData
    table_name = NiAl_data/NiAl_mobility_data_FCC_invert.csv
  [../]
  
  [./table_Liquid_mobility_data]
    type = BinaryConjugateMobilityData
    table_name = NiAl_data/NiAl_mobility_data_Liquid_invert.csv
  [../] 
 
[]


[Postprocessors]

#After each time step collect the data
#TIMSTEP_END
#TIMESTEP_BEGIN

 #Calculates the total mole fraction
 #of Al in the simulation

 [./overall_Al]
  type = Moles
  variable = x_Al
  volume = 6000
 [../]
 
 [./pos_alpha]
   type = FindValueOnLine
   target = 0.5
   v = phi_beta
   start_point = '-3000 0 0'
   end_point = '3000 0 0'
   tol = 1e-8
   execute_on = 'initial timestep_end'
  [../]
  
  [./pos_beta]
   type = FindValueOnLine
   target = 0.5
   v = phi_alpha
   start_point = '-3000 0 0'
   end_point = '3000 0 0'
   tol = 1e-8
   execute_on = 'initial timestep_end'
  [../]
  
  [./volumefraction_alpha]
    type = VolumeFraction
    variable = phi_alpha
    execute_on = 'initial timestep_end'
  [../]

  [./L_AlAl_Liq]
    type = ElementIntegralMaterialProperty
    mat_prop = L_BB_alpha
    execute_on = 'initial timestep_end'
  [../]

  [./L_AlAl_fcc]
     type = ElementIntegralMaterialProperty
     mat_prop = L_BB_beta
     execute_on = 'initial timestep_end'
  [../]
[]

#[VectorPostprocessors]
  #[./x_Al]
  #    type = LineValueSampler
  #    num_points = 1201
  #    start_point = '-3000 0 0'
  #    end_point = '3000 0 0'
  #    sort_by = id 
  #    variable = 'x_Al'
  #    execute_on = 'initial final'
  # [../]

  #[./Mobilities]
  #   type = LineMaterialRealSampler
  #   start = '-3000 0 0'
  #   end   = '3000 0 0'
  #   sort_by = id
  #   property = 'L_BB'
  #[../]
#[../]

[Preconditioning]
  [./full]
    type = SMP
    full = true
  [../]
[]

[Executioner]

  type = Transient
  scheme = bdf2 
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm       ilu          nonzero'
  
  l_tol = 1.0e-5
  l_max_its = 15
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-10
  nl_max_its = 30
  
  #non-dimensional time
  #0.1sec = 6.491820310e5
  #2 min  = 7.790184368e8 

  start_time = 0
  end_time = 7.790184368e12 #non-dimensional
  #dt = 1
  
  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    percent_change = 0.1
    dt = 1.0e-0 #real_dt = 0.1sec
  [../]
  
[]

[Outputs]

  file_base = NiAl_twophase_1000
  exodus = true
  
  [./csv]
    type = CSV
    execute_on = 'initial final'
    execute_postprocessors_on = 'initial timestep_end'
  [../]
  
   [pgraph]
    type = PerfGraphOutput
    execute_on = 'final'
    level = 1  
    heaviest_branch = false
  []
  
[]

[Debug]
  show_var_residual_norms = true
[] 
