#This kernel uses the TabulatedConjugatePhaseData
#Material to solve the equations. 
#This is a model Al-Cr-Ni alloy problem
#alpha phase- FCC_A1 (gamma)
#beta phase- FCC_L12 (gamma_prime)
#T = 1473K
#Simulation details:
#Interfacial energy = 6.9E-3 J/m2
#interfacial width = 4.0e-6 m
#X-length = 100e-6 m
#Y-length = 100e-6 m
#characteristic_length = (0.2E-6)/6.0 = 0.03333e-6 m
#characteristic_time = 0.0140 sec
#Nondimensional length = (Length/characteristic_length)

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 150
  ny = 150
  xmin = -1500
  xmax =  1500
  ymin = -1500
  ymax =  1500
  elem_type = QUAD4
[]


[Variables]

#alpha is the matrix phase
#alpha = fcc

[./phi_alpha]
  order = FIRST
  family = LAGRANGE
[../]

#beta is the ppt phase
#beta = L12_gamma_prime

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

[./x_Cr]
  order = FIRST
  family = LAGRANGE
[../]

[./mu_Cr]
  order = FIRST
  family = LAGRANGE
[../]
[]


[BCs]

  #[./alpha_pf]
  #  type = NeumannBC
  #  variable = 'phi_alpha'
  #  boundary = 'left right'
  #  value = 0
  #[../]
  
  #[./beta_pf]
  #  type = NeumannBC
  #  variable = 'phi_beta'
  #  boundary = 'left right'
  #  value = 0
  #[../]

   #[./MassFluxBC_Al]
   # type = NeumannBC
   # variable = 'mu_Al'
   # boundary = 'left right'
   # value = 0
   #[../]

   #[./MassFluxBC_Cr]
   # type = NeumannBC
   # variable= 'mu_Cr'
   # boundary = 'left right'
   # value = 0
   #[../]
   
   [./Periodic]
    [./auto]
      variable = 'phi_alpha phi_beta mu_Al mu_Cr x_Al x_Cr'
      auto_direction = 'x y'
    [../]
   [../]

[]


[ICs]
  
  [./ic_phi_alpha]
    type = FunctionIC
    variable = phi_alpha
    function = func_phi_alpha
  [../]
  
  [./ic_phi_beta]
    type = FunctionIC
    variable = phi_beta
    function = func_phi_beta
  [../]
  
  [./ic_xAl]
   type = FunctionIC
   variable = x_Al
   function = func_xAl
  [../]
  
  [./muAl_IC]
   type = FunctionIC
   variable = mu_Al
   function = func_muAl
  [../]
  
  [./ic_xCr]
   type = FunctionIC
   variable = x_Cr
   function = func_xCr
  [../]
  
  [./muCr_IC]
   type = FunctionIC
   variable = mu_Cr
   function = func_muCr
  [../]
  
[]

[Functions]
   
  #FCC phase 
  [./func_phi_alpha]
    type = ParsedFunction
    vals = func_phi_beta
    vars = gamma_prime
    value = 1-gamma_prime
  [../]
  
  #L12 phase
  [./func_phi_beta]
    type = CircleFunctionIC
    inside_value = 1.0
    outside_value = 0.0
    table_name = 'AlCrNi_data/outputfile.dat'
  [../]

   #Equilibrium mole fraction of Al
   #xAl_FCC(alpha) = 1.63151E-1
   #xAl_L12(beta)  = 2.09074E-1
    
  [./func_xAl]
    type = ParsedFunction
    vars = 'alpha beta'
    vals = 'func_phi_alpha  func_phi_beta'
    value = '1.63E-1*alpha + 2.095783E-1*beta'
  [../]
  
   #Equilibrium mole fraction of Cr
   #xCr_FCC(alpha) = 7.40827E-2
   #xCr_L12(beta) =  4.67084E-2
  
  [./func_xCr]
    type = ParsedFunction
    vars = 'alpha beta'
    vals = 'func_phi_alpha  func_phi_beta'
    value = '9.189272E-2*alpha + 5.446593E-2*beta'
  [../]
  
  
  [./func_muAl]
    type = ParsedFunction
    vars = 'alpha beta'
    vals = 'func_phi_alpha  func_phi_beta'
    value = '-8.16518E0*alpha -8.16518E0*beta'
  [../]
  
  [./func_muCr]
    type = ParsedFunction
    vars = 'alpha beta'
    vals = 'func_phi_alpha  func_phi_beta'
    value = '2.56473E-1*alpha + 2.56743E-1*beta'
  [../]
[]


[Kernels]
  
  #Variable that this kernel operates on
  #is the diffusion potential of comp B
  
  [./PhaseConc_Al]
    type = TCPhaseConstraintMuB
    variable = mu_Al
    xB = x_Al
    C_diff_pot = mu_Cr
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    h_alpha = h_alpha
    h_beta = h_beta
  [../]
  
  #Variable that this kernel operates on
  #is the diffusion potential of comp C
  
  [./PhaseConc_Cr]
    type = TCPhaseConstraintMuC
    variable = mu_Cr
    xC = x_Cr
    B_diff_pot = mu_Al
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    h_alpha = h_alpha
    h_beta = h_beta
  [../]

  ######################
  # Diffusion Equation #
  # Component B        #
  ######################
    
  [./Al_balance]
    type = TCContinuityEquationB
    variable = x_Al
    xC = x_Cr
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    h_alpha = h_alpha
    h_beta = h_beta
  [../]

  [./xAl_dot]
    type = TimeDerivative
    variable = x_Al
  [../]
  
  ######################
  # Diffusion Equation #
  # Component C        #
  ######################
    
  [./Cr_balance]
    type = TCContinuityEquationC
    variable = x_Cr
    xB = x_Al
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    h_alpha = h_alpha
    h_beta = h_beta
  [../]

  [./xCr_dot]
    type = TimeDerivative
    variable = x_Cr
  [../]
  
  #######################
  # Allen-Cahn equation #
  # for phase-alpha     #
  #######################
  
  [./MultiCompDrivingForce_alpha]
    type = MultiCompDrivingForce
    variable = phi_alpha     #Phase_1
    phase_2 = phi_beta       #Coupled phase_2
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    A_chem_pot_1 = A_chem_pot_alpha
    A_chem_pot_2 = A_chem_pot_beta
    xB_1 = xB_alpha
    xB_2 = xB_beta
    xC_1 = xC_alpha
    xC_2 = xC_beta
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
    type = MultiCompDrivingForce
    variable = phi_beta
    phase_2 = phi_alpha
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    A_chem_pot_1 = A_chem_pot_beta
    A_chem_pot_2 = A_chem_pot_alpha
    xB_1 = xB_beta
    xB_2 = xB_alpha
    xC_1 = xC_beta
    xC_2 = xC_alpha
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

  ##########################
  ## Constant parameters  ##
  ##########################
  
  [./ConstantFieldProperties]
     type = GenericConstantMaterial
     prop_names =  'L_phi      kappa   m   nd_factor'
     prop_values = '2.1034e-6  1800.0   1.0  1.5777e4'
  [../]
  
  [./InterpolationFunction]
    type = QuantInterpolationFunction
    phase_alpha = phi_alpha
    phase_beta = phi_beta
  [../]
  
  [./TernaryConjugatePhaseMaterial_FCC]
    type = TernaryConjugatePhaseMaterial
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    A_chem_pot = A_chem_pot_alpha
    B_mole_fraction = xB_alpha
    C_mole_fraction = xC_alpha
    inv_B_tf   = inv_B_tf_alpha
    inv_BC_tf  = inv_BC_tf_alpha
    inv_C_tf   = inv_C_tf_alpha
    table_object = table_FCC_data
    #outputs = exodus
 [../]
 
 [./TernaryConjugatePhaseMaterial_L12]
    type = TernaryConjugatePhaseMaterial
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    A_chem_pot = A_chem_pot_beta
    B_mole_fraction = xB_beta
    C_mole_fraction = xC_beta
    inv_B_tf   = inv_B_tf_beta
    inv_BC_tf  = inv_BC_tf_beta
    inv_C_tf   = inv_C_tf_beta
    table_object = table_L12_data
    #outputs = exodus
 [../]
 
 [./TernaryConjugateKineticMaterial_FCC]
    type = TernaryConjugateKineticMaterial
    L_BB = L_BB_alpha
    L_BC = L_BC_alpha
    L_CC = L_CC_alpha
    dL_BB_muB = dL_BB_muB_alpha
    dL_BC_muB = dL_BC_muB_alpha
    dL_CC_muB = dL_CC_muB_alpha
    dL_BB_muC = dL_BB_muC_alpha
    dL_BC_muC = dL_BC_muC_alpha
    dL_CC_muC = dL_CC_muC_alpha
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    table_object = table_FCC_mobility_data
    #outputs =exodus
  [../]
  
  [./TernaryConjugateKineticMaterial_L12]
    type = TernaryConjugateKineticMaterial
    L_BB = L_BB_beta
    L_BC = L_BC_beta
    L_CC = L_CC_beta
    dL_BB_muB = dL_BB_muB_beta
    dL_BC_muB = dL_BC_muB_beta
    dL_CC_muB = dL_CC_muB_beta
    dL_BB_muC = dL_BB_muC_beta
    dL_BC_muC = dL_BC_muC_beta
    dL_CC_muC = dL_CC_muC_beta
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    table_object = table_L12_mobility_data
    #outputs =exodus
  [../]
   


[]  
  
[UserObjects]

  #This userobjcts reads the data from the table 
  #for a given diffusion potential
  #and provides the interpolated grand-potential,
  #mole fraction and inverse of the thermodynamic_factor
  #to the material TabulatedPhaseMaterial
  
  [./table_FCC_data]
    type = TernaryConjugatePhaseData
    table_name = AlCrNi_data/AlCrNi_invert_data_FCC_A1.csv
  [../]
  
  [./table_L12_data]
    type = TernaryConjugatePhaseData
    table_name = AlCrNi_data/AlCrNi_invert_data_L12.csv
  [../] 
  
  [./table_FCC_mobility_data]
    type = TernaryConjugateMobilityData
    table_name = AlCrNi_data/AlCrNi_mobility_data_FCC_invert.csv
  [../]
  
  [./table_L12_mobility_data]
    type = TernaryConjugateMobilityData
    table_name = AlCrNi_data/AlCrNi_mobility_data_L12_invert.csv
  [../]
  
[] 

[Postprocessors]

  [./flood_count]
    type = FeatureFloodCount
    variable = phi_beta
    compute_var_to_feature_map = true
    threshold = 0.5
    execute_on = 'initial timestep_end'
  [../]
  
  [./mean_radius]
    type = MeanGrainRadius
    feature_counter = flood_count
    execute_on = 'initial timestep_end'
  [../]
 
  [./volumefraction_gp]
   type = VolumeFraction
   variable = phi_beta
   execute_on = 'initial timestep_end'
 [../]
  
[]

#[VectorPostprocessors]

  #[./xAl]
    #type = LineValueSampler
    #start_point = '-600 0 0'
    #end_point = '600 0 0'
    #num_points = 1200
    #variable = 'x_Al'
    #sort_by = id
    #execute_on = 'initial timestep_end'
  #[../]
  
  #[./xCr]
    #type = LineValueSampler
    #start_point = '-600 0 0'
    #end_point = '600 0 0'
    #num_points = 1200
    #variable = 'x_Cr'
    #sort_by = id
    #execute_on = 'initial timestep_end'
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
  #petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  #petsc_options_value = 'asm       31          preonly lu  2'
  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm ilu nonzero'
  l_tol = 1.0e-5
  l_max_its = 15
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-10
  nl_max_its = 30

  #Real time 1day  = 2.46434683e5
  #Real time 7days = 1.72504278e6
  
  end_time = 8.64e7 #=14 days
  #dt = 1
  
  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 1e0  #real_dt = 0.1sec
    percent_change = 0.1
  [../]
  
[]

[Outputs]

  file_base = AlCrNi_coarsening_2D
  
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
