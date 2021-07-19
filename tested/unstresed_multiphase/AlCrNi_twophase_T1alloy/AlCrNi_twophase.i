#This kernel uses the TabulatedConjugatePhaseData
#Material to solve the equations. 
#This is a model Al-Cr-Ni alloy problem
#beta phase- FCC_A1
#alpha phase- BCC_B2
#T = 1473K
#Simulation details:
#Interfacial energy = 0.5 J/m2
#interfacial width = 0.6e-6 m
#Length = 200e-6 m
#characteristic_length = (0.20e-6/6.0) m = 0.033e-6 m
#characteristic_time = 4.4706 sec
#Al_diff_pot_eqm = -7.412944e4 J/mol 
#Cr_diff_pot_eqm =  2.991628e4 J/mol
#Nondimensional length = Length/characteristic_length
#No. of finite elements = 6.0(Length/interface_width)

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 2000
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
  
   [./MassFluxBC_Al]
    type = NeumannBC
    variable = 'mu_Al'
    boundary = 'left right'
    value = 0
  [../]
  
  [./MassFluxBC_Cr]
    type = NeumannBC
    variable= 'mu_Cr'
    boundary = 'left right'
    value = 0
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
   
  #BCC phase 
  [./func_phi_alpha]
    type = ParsedFunction
    value = 'if(x<=0.0, 1.0, 0)'
  [../]
  
  #FCC phase
  [./func_phi_beta]
    type = ParsedFunction
    value = 'if(x>0.0&x<=3000.0, 1.0, 0)'
  [../]

   #Equilibrium mole fraction of Al
   #xAl_FCC(beta) = 1.231521E-1
   #xAl_BCC(alpha) = 1.512073E-2
    
  [./func_xAl]
    type = ParsedFunction
    vars = 'alpha beta'
    vals = 'func_phi_alpha  func_phi_beta'
    value = '9.39852E-4*alpha + 9.887855E-2*beta'
  [../]
  
   #Equilibrium mole fraction of Cr
   #xCr_FCC(beta) = 3.443197E-1
   #xCr_BCC(alpha) = 8.817211E-1
  
  [./func_xCr]
    type = ParsedFunction
    vars = 'alpha beta'
    vals = 'func_phi_alpha  func_phi_beta'
    value = '9.54264E-1*alpha + 3.7311485E-1*beta'
  [../]
  
  [./func_muAl]
    type = ParsedFunction
    vars = 'alpha beta'
    vals = 'func_phi_alpha  func_phi_beta'
    value = '-8.13286E0*alpha -6.51209298*beta'
  [../]
  
  [./func_muCr]
    type = ParsedFunction
    vars = 'alpha beta'
    vals = 'func_phi_alpha  func_phi_beta'
    value = '2.83095E0*alpha + 2.35630328E0*beta'
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
    phase_beta  = phi_beta
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
     prop_names =  'L_phi    kappa    m    nd_factor'
     prop_values = '0.0490   40.5    1.0    32.6574 '
  [../]
  
 [./InterpolationFunction]
   type = QuantInterpolationFunction
   phase_alpha = phi_alpha
   phase_beta = phi_beta
 [../]
  
 [./TernaryConjugatePhaseMaterial_BCC_A2]
    type = TernaryConjugatePhaseMaterial
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    A_chem_pot = A_chem_pot_alpha
    B_mole_fraction = xB_alpha
    C_mole_fraction = xC_alpha
    inv_B_tf   = inv_B_tf_alpha
    inv_BC_tf  = inv_BC_tf_alpha
    inv_C_tf   = inv_C_tf_alpha
    table_object = table_BCC_A2_data
    #outputs = exodus
 [../]
 
 [./TernaryConjugatePhaseMaterial_FCC_A1]
    type = TernaryConjugatePhaseMaterial
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    A_chem_pot = A_chem_pot_beta
    B_mole_fraction = xB_beta
    C_mole_fraction = xC_beta
    inv_B_tf   = inv_B_tf_beta
    inv_BC_tf  = inv_BC_tf_beta
    inv_C_tf   = inv_C_tf_beta
    table_object = table_FCC_A1_data
    #outputs = exodus
 [../]
 
 [./TernaryConjugateKineticMaterial_BCC_A2]
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
    table_object = table_BCC_A2_mobility_data
    #outputs =exodus
  [../]
  
  [./TernaryConjugateKineticMaterial_FCC_A1]
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
    table_object = table_FCC_A1_mobility_data
    #outputs =exodus
  [../]
   


[]  
  
[UserObjects]

  #This userobjcts reads the data from the table 
  #for a given diffusion potential
  #and provides the interpolated grand-potential,
  #mole fraction and inverse of the thermodynamic_factor
  #to the material TabulatedPhaseMaterial
  
  [./table_BCC_A2_data]
    type = TernaryConjugatePhaseData
    table_name = AlCrNi_data/AlCrNi_invert_data_BCC_A2.csv
  [../]
  
  [./table_FCC_A1_data]
    type = TernaryConjugatePhaseData
    table_name = AlCrNi_data/AlCrNi_invert_data_FCC_A1.csv
  [../] 
  
  [./table_BCC_A2_mobility_data]
    type = TernaryConjugateMobilityData
    table_name = AlCrNi_data/AlCrNi_mobility_data_BCC_A2_invert.csv
  [../]
  
  [./table_FCC_A1_mobility_data]
    type = TernaryConjugateMobilityData
    table_name = AlCrNi_data/AlCrNi_mobility_data_FCC_A1_invert.csv
  [../]
  
[] 

[Postprocessors]

  [./computation_time]
    type = PerfGraphData
    section_name = Transient::PicardSolve
    data_type = total
    execute_on = timestep_end
  [../]

  [./pos_alpha]
    type = FindValueOnLine
    target= 0.5
    v = phi_alpha
    start_point = '-3000 0 0'
    end_point = '3000 0 0'
    tol = 1e-8
    execute_on = 'initial timestep_end'
  [../]
  
  [./pos_beta]
    type = FindValueOnLine
    target = 0.5
    v = phi_beta
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
  
[]

#[VectorPostprocessors]
  #[./xAl]
    #type = LineValueSampler
    #start_point = '-1000 0 0'
    #end_point = '1000 0 0'
    #num_points = 2000
    #variable = 'x_Al'
    #sort_by = id
    #execute_on = 'initial timestep_end'
  #[../]
  
  #[./xCr]
    #type = LineValueSampler
    #start_point = '-1000 0 0'
    #end_point = '1000 0 0'
    #num_points = 2000
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
  l_max_its = 50
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-10
  nl_max_its = 50
  
  #real_time = 1day = 2.14695190e3
  #real_time = 150 days = 2.89893974e6
  
  start_time = 0
  end_time   = 2.89893974E6 #non-dimensional
  #dt = 1
  
  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 10.0e0  #real_dt = 0.1sec
    percent_change = 0.1
  [../]
  
[]

[Outputs]

  file_base = AlCrNi_twophase
  exodus = true

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
