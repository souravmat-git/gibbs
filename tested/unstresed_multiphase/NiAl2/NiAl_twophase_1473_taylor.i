#This example file uses taylor approximation
#and constant Onsager mobilities
#This is a model Ni-Al alloy problem
#beta phase- FCC_L12
#alpha phase- FCC_A1
#T = 1473K
#Simulation details:
#Interfacial energy = 0.5 J/m2
#interfacial width = 10e-6 m
#Length = 6000e-6 m
#characteristic_length = 1e-6 m
#time = 43.2 hrs
#Al_diff_pot_eqm = -1.085307e5 J/mol 
#Nondimensional_length = (Length/characteristic_length)
#No. of finite elements = 6.0(Length/interface_width) 

[Mesh]
  type = GeneratedMesh
  dim  = 1
  nx   =  3600
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

  #beta(left)- L12
  
   [./ic_phi_beta]
    type = BoundingBoxIC
    variable = phi_beta
    x1 = -3000
    y1 =  0.0
    x2 =  0.0
    y2 =  0.0
    inside = 1.0
    outside = 0.0
  [../]

 #alpha(right)- fcc
    
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

 #Al_inside = beta (left)
 #Al_outside = alpha(right)
 
  [./ic_func_xAl]
   type = BoundingBoxIC
   variable = x_Al
   x1 = -3000
   y1 = 0.0
   x2 = 0.0
   y2 = 0.0
   inside =  2.5E-1         #xAl_L12 = 2.30730E-1
   outside = 1.0E-4         #xAl_FCC = 1.83922E-1    
  [../]
  
  [./muAl_IC]
   variable = mu_Al
   type = ConstantIC
   value = -8.86217
  [../]
  
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
     prop_names =  'L_phi     kappa    m     nd_factor'
     prop_values = '0.0046    12.50    1.0    544.2899'
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
    phase_beta  = phi_beta
  [../]  
  
  [./GPtaylorapproximation_L12]
    type = GPTaylorApproximation
    A_chem_pot = A_chem_pot_beta
    xB = xB_beta
    inv_B_tf = inv_B_tf_beta
    B_diff_pot = mu_Al
    char_energy = 1.225e4
    xB_eqm = 2.30730E-1
    B_tf_eqm = 2.8603E5
    B_diff_pot_eqm = -1.08531e5
    A_chem_pot_eqm = -8.56263E4
    #outputs = exodus
 [../] 
 
 [./GPtaylorApproximation_FCC]
    type = GPTaylorApproximation
    A_chem_pot = A_chem_pot_alpha
    xB = xB_alpha
    inv_B_tf = inv_B_tf_alpha
    B_diff_pot = mu_Al
    char_energy = 1.225e4
    xB_eqm = 1.83922E-1
    B_tf_eqm = 3.6927E5
    B_diff_pot_eqm = -1.08531e5
    A_chem_pot_eqm = -8.56263E4
    #outputs = exodus
  [../]

  #The dimensional value of Al mobility in L12 phase 
  #is 1.1094E-18 (mol m2)/Js. To nondimensionalize this quantity,
  #we define a characteristic Onsager constant 
  #Lc = D/RT, where D(Al,Al,Ni,FCC) = 2.6553E-13 m2/s
  #This value is calculated at the equilibrium mole fraction.  
  
  [./BinaryConstantKineticMaterial_L12]
    type = BinaryConstantKineticMaterial
    L_BB = L_BB_beta
    dL_BB_muB = dL_BB_muB_beta
    L_BB_val = 0.0512
    #outputs = exodus
  [../]

 #The dimensional value of Al mobility in FCC phase
 #is 7.19070E-19 (mol m2)/Js. To nondimensionalize this quantity,
 # we define a characteristic Onsager constant
 #Lc = D/RT, where D(Al,Al,Ni,FCC) =  2.6553E-13 m2/s
 #Thos value is calculated at the equilibrium mole fraction
  
  [./BinaryConstantKineticMaterial_FCC]
    type = BinaryConstantKineticMaterial
    L_BB = L_BB_alpha
    dL_BB_muB = dL_BB_muB_alpha
    L_BB_val = 0.0332
    #outputs = exodus
  [../]

 
[]

[Postprocessors]

 [./computation_time]
    type = PerfGraphData
    section_name = Transient::PicardSolve
    data_type = total
    execute_on = 'timestep_end'
 [../]
 
 [./pos_alpha]
   type = FindValueOnLine
   target = 0.5
   v = phi_beta
   start_point = '-3000 0 0'
   end_point = '3000 0 0'
   tol = 5e-8
   execute_on = 'initial timestep_end'
  [../]
  
  [./pos_beta]
   type = FindValueOnLine
   target = 0.5
   v = phi_alpha
   start_point = '-3000 0 0'
   end_point = '3000 0 0'
   tol = 5e-8
   execute_on = 'initial timestep_end'
  [../]
  
  #beta is L12 phase
  [./volumefraction_L12]
    type = VolumeFraction
    variable = phi_beta
    execute_on = 'initial timestep_end'
  [../]
  
[]


#[VectorPostprocessors]
   #[./xAl]
   # type = LineValueSampler
   # start_point = '-3000 0 0' 
   # end_point   = '3000 0 0'
   # variable    = 'x_Al'
   # num_points  = 3601
   # sort_by = id
   # execute_on = 'initial final'
   #[../]
#[]

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

  #144 hours = 1.64498318e5
  #43.2 hours = 49349.495462334
  #50 days  = 1147105.682421667
  #300 days = 6882634.094530005

  end_time = 1147105.682421667 #non-dimensional
  #dt = 1
  
  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    percent_change = 0.1
    dt = 1.0 #real_dt = 0.1sec
  [../]
  
[]

[Outputs]

  file_base = NiAl_twophase_1473_taylor
  #exodus = true

  [./exodus]
    type = Exodus
    execute_on = 'initial timestep_end'
    interval = 50
  [../]
  
  [./csv]
    type = CSV
    execute_on = 'initial timestep_end'
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
