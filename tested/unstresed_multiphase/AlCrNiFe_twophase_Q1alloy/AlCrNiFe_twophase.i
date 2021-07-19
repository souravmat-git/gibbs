#This input properties are tabulated
#as functions of independent diffusion potentials
#That is Al, Cr,and Fe. Molar GP = Ni-chemical potential 
#This is a model Al-Cr-Ni-Fe alloy problem
#beta phase is FCC
#alpha phase is BCC
#T = 1473K
#Simulation details:
#Interfacial energy = 0.5 J/m2
#interfacial width = 1e-6 m
#Length = 250e-6 m
#characteristic_length = (0.25e-6)/6 m = 0.041667e-6 m
#characteristic_time = 3.8171 s
#Nondimensional length = Length/characteristic_length
#Number of finite elements(nx)= 6.0(Length/interface_width)

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1500
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

[./x_Fe]
  order = FIRST
  family = LAGRANGE
[../]

[./mu_Fe]
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

   [./MassFlux_Cr]
    type = NeumannBC
    variable = 'mu_Cr'
    boundary = 'left right'
    value = 0
   [../]

   [./MassFlux_Fe]
    type = NeumannBC
    variable = 'mu_Fe'
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
  
   [./ic_xFe]
   type = FunctionIC
   variable = x_Fe
   function = func_xFe
  [../]
  
  [./muFe_IC]
   type = FunctionIC
   variable = mu_Fe
   function = func_muFe
  [../]
  
[]

[Functions]

  [./func_phi_alpha]
    type = ParsedFunction
    value = 'if(x<=0.0, 1.0, 0)'
  [../]
  
  [./func_phi_beta]
    type = ParsedFunction
    value = 'if(x>0.0&x<=3000.0, 1.0, 0)'
  [../]
  
  #Equilibrium composition of Al
  #BCC_alpha xAl = 1.64141E-2
  #FCC_beta  xAl = 1.08750E-1
  
  [./func_xAl]
    type = ParsedFunction
    vars = 'alpha beta'
    vals = 'func_phi_alpha  func_phi_beta'
    value = '1.189807E-3*alpha + 9.5E-2*beta'
  [../]
 
 #Equilibrium composition of Cr
 #BCC_alpha xCr = 8.36609E-1
 #FCC_beta xCr =  3.49752E-1
  [./func_xCr]
    type = ParsedFunction
    vars = 'alpha beta'
    vals = 'func_phi_alpha  func_phi_beta'
    value = '9.791539E-1*alpha + 3.6E-1*beta'
  [../]
 
 #Equilibrium composition of Fe
 #BCC_alpha xFe = 2.88165E-2
 #FCC_beta xFe  = 3.21495E-2 
  
  [./func_xFe]
    type = ParsedFunction
    vars = 'alpha beta'
    vals = 'func_phi_alpha  func_phi_beta'
    value = '9.917367E-3*alpha + 5.0E-2*beta'
  [../]
  
  [./func_muAl]
    type = ParsedFunction
    vars = 'alpha beta'
    vals = 'func_phi_alpha  func_phi_beta'
    value = '-6.360128E0*alpha -6.085540E0*beta'
  [../]
  
  [./func_muCr]
    type = ParsedFunction
    vars = 'alpha beta'
    vals = 'func_phi_alpha  func_phi_beta'
    value = '4.144341E0*alpha + 2.44464E0*beta'
  [../]
  
  [./func_muFe]
    type = ParsedFunction
    vars = 'alpha beta'
    vals = 'func_phi_alpha  func_phi_beta'
    value = '-9.437609E-1*alpha -1.64241E0*beta'
  [../]
  
[]


[Kernels]
  
  #Variable that this kernel operates on
  #is the diffusion potential of comp B
  
  [./PhaseConc_Al]
    type = MCPhaseConstraintMuB
    variable = mu_Al
    xB = x_Al
    C_diff_pot = mu_Cr
    D_diff_pot = mu_Fe
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    h_alpha = h_alpha
    h_beta = h_beta
  [../]
  
  [./PhaseConc_Cr]
    type = MCPhaseConstraintMuC
    variable = mu_Cr
    xC = x_Cr
    B_diff_pot = mu_Al
    D_diff_pot = mu_Fe
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    h_alpha = h_alpha
    h_beta = h_beta
  [../]
  
  [./PhaseConc_Fe]
    type = MCPhaseConstraintMuD
    variable = mu_Fe
    xD = x_Fe
    C_diff_pot = mu_Cr
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
    type = MCContinuityEquationB
    variable = x_Al
    xC = x_Cr
    xD = x_Fe
    B_diff_pot  = mu_Al
    C_diff_pot  = mu_Cr
    D_diff_pot  = mu_Fe
    phase_alpha = phi_alpha
    phase_beta  = phi_beta
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
    type = MCContinuityEquationC
    variable = x_Cr
    xB = x_Al
    xD = x_Fe
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    D_diff_pot = mu_Fe
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    h_alpha = h_alpha
    h_beta = h_beta
  [../]

  [./xCr_dot]
    type = TimeDerivative
    variable = x_Cr
  [../]
  
  ######################
  # Diffusion Equation #
  # Component D        #
  ######################
    
  [./Fe_balance]
    type = MCContinuityEquationD
    variable = x_Fe
    xB = x_Al
    xC = x_Cr
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    D_diff_pot = mu_Fe
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    h_alpha = h_alpha
    h_beta = h_beta
  [../]

  [./xFe_dot]
    type = TimeDerivative
    variable = x_Fe
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
    D_diff_pot  = mu_Fe
    A_chem_pot_1 = A_chem_pot_alpha
    A_chem_pot_2 = A_chem_pot_beta
    xB_1 = xB_alpha
    xB_2 = xB_beta
    xC_1 = xC_alpha
    xC_2 = xC_beta
    xD_1 = xD_alpha
    xD_2 = xD_beta
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
    D_diff_pot = mu_Fe
    A_chem_pot_1 = A_chem_pot_beta
    A_chem_pot_2 = A_chem_pot_alpha
    xB_1 = xB_beta
    xB_2 = xB_alpha
    xC_1 = xC_beta
    xC_2 = xC_alpha
    xD_1 = xD_alpha
    xD_2 = xD_beta
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
  ##Phenomenological coeff##
  
   [./ConstantFieldProperties]
     type = GenericConstantMaterial
     prop_names =  'L_phi    kappa  m    nd_factor'
     prop_values = '0.0150   72.0   1.0   54.4290 '
  [../] 
 
 ######################################
 
 [./QuaternaryConjugateKineticMaterial_BCC]
  type = QuaternaryConjugateKineticMaterial
  L_BB = L_BB_alpha
  L_CC = L_CC_alpha
  L_DD = L_DD_alpha
  L_BC = L_BC_alpha
  L_BD = L_BD_alpha
  L_CD = L_CD_alpha
  #w.r.t B
  dL_BB_muB = dL_BB_muB_alpha
  dL_CC_muB = dL_CC_muB_alpha
  dL_DD_muB = dL_DD_muB_alpha
  dL_BC_muB = dL_BC_muB_alpha
  dL_BD_muB = dL_BD_muB_alpha
  dL_CD_muB = dL_CD_muB_alpha
  #w.r.t C
  dL_BB_muC = dL_BB_muC_alpha
  dL_BC_muC = dL_BC_muC_alpha
  dL_BD_muC = dL_BD_muC_alpha
  dL_CC_muC = dL_CC_muC_alpha
  dL_CD_muC = dL_CD_muC_alpha
  dL_DD_muC = dL_DD_muC_alpha
  #w.r.t D
  dL_BB_muD = dL_BB_muD_alpha
  dL_BC_muD = dL_BC_muD_alpha
  dL_BD_muD = dL_BD_muD_alpha
  dL_CC_muD = dL_CC_muD_alpha
  dL_CD_muD = dL_CD_muD_alpha
  dL_DD_muD = dL_DD_muD_alpha
  B_diff_pot = mu_Al
  C_diff_pot = mu_Cr
  D_diff_pot = mu_Fe
  table_object = table_BCC_A2_mobility_data
  #outputs = exodus
 [../]
 
 [./QuaternaryConjugateKineticMaterial_FCC]
  type = QuaternaryConjugateKineticMaterial
  L_BB = L_BB_beta
  L_CC = L_CC_beta
  L_DD = L_DD_beta
  L_BC = L_BC_beta
  L_BD = L_BD_beta
  L_CD = L_CD_beta
  #w.r.t B
  dL_BB_muB = dL_BB_muB_beta
  dL_BC_muB = dL_BC_muB_beta
  dL_BD_muB = dL_BD_muB_beta
  dL_CC_muB = dL_CC_muB_beta
  dL_CD_muB = dL_CD_muB_beta
  dL_DD_muB = dL_DD_muB_beta
  #w.r.t C
  dL_BB_muC = dL_BB_muC_beta
  dL_BC_muC = dL_BC_muC_beta
  dL_BD_muC = dL_BD_muC_beta
  dL_CC_muC = dL_CC_muC_beta
  dL_CD_muC = dL_CD_muC_beta
  dL_DD_muC = dL_DD_muC_beta
  #w.r.t D
  dL_BB_muD = dL_BB_muD_beta
  dL_BC_muD = dL_BC_muD_beta
  dL_BD_muD = dL_BD_muD_beta
  dL_CC_muD = dL_CC_muD_beta
  dL_CD_muD = dL_CD_muD_beta
  dL_DD_muD = dL_DD_muD_beta
  B_diff_pot = mu_Al
  C_diff_pot = mu_Cr
  D_diff_pot = mu_Fe
  table_object = table_FCC_A1_mobility_data
  #outputs = exodus
 [../]

 [./InterpolationFunction]
   type = QuantInterpolationFunction
   phase_alpha = phi_alpha
   phase_beta = phi_beta
 [../]
  
  [./QuaternaryConjugatePhaseMaterial_BCC_A2]
    type = QuaternaryConjugatePhaseMaterial
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    D_diff_pot = mu_Fe
    A_chem_pot = A_chem_pot_alpha
    B_mole_fraction = xB_alpha
    C_mole_fraction = xC_alpha
    D_mole_fraction = xD_alpha
    inv_B_tf = inv_B_tf_alpha
    inv_C_tf = inv_C_tf_alpha
    inv_D_tf = inv_D_tf_alpha
    inv_BC_tf = inv_BC_tf_alpha
    inv_BD_tf = inv_BD_tf_alpha
    inv_CD_tf = inv_CD_tf_alpha
    table_object = table_BCC_A2_data
    #outputs =exodus
  [../]
  
  [./QuaternaryConjugatePhaseMaterial_FCC_A1]
    type = QuaternaryConjugatePhaseMaterial
    B_diff_pot = mu_Al
    C_diff_pot = mu_Cr
    D_diff_pot = mu_Fe
    A_chem_pot = A_chem_pot_beta
    B_mole_fraction = xB_beta
    C_mole_fraction = xC_beta
    D_mole_fraction = xD_beta
    inv_B_tf = inv_B_tf_beta
    inv_C_tf = inv_C_tf_beta
    inv_D_tf = inv_D_tf_beta
    inv_BC_tf = inv_BC_tf_beta
    inv_BD_tf = inv_BD_tf_beta
    inv_CD_tf = inv_CD_tf_beta
    table_object = table_FCC_A1_data
    #outputs =exodus
  [../]
  
[]

[UserObjects]

  [./table_BCC_A2_data]
    type = QuaternaryConjugatePhaseData
    table_name = '../AlCrNiFe_data/AlCrNiFe_invert_data_BCC_A2.csv'
  [../]
  
  [./table_FCC_A1_data]
    type = QuaternaryConjugatePhaseData
    table_name = '../AlCrNiFe_data/AlCrNiFe_invert_data_FCC_A1.csv'
  [../]
  
  [./table_BCC_A2_mobility_data]
    type = QuaternaryConjugateMobilityData
    table_name = '../AlCrNiFe_data/AlCrNiFe_mobility_data_BCC_A2_invert.csv'  
  [../]
  
  [./table_FCC_A1_mobility_data]
    type = QuaternaryConjugateMobilityData
    table_name = '../AlCrNiFe_data/AlCrNiFe_mobility_data_FCC_A1_invert.csv'
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
    start_point = '-300 0 0'
    end_point = '3000 0 0'
    tol = 1e-8
    execute_on = 'initial timestep_end'
  [../]
  
  [.volumefraction_alpha]
    type = VolumeFraction
    variable = phi_alpha
    execute_on = 'initial timestep_end'
  [../]
[]

[VectorPostprocessors]

    [./xAl]
       type = LineValueSampler
       start_point = '-3000 0 0'
       end_point = '3000 0 0'
       num_points = 1501
       variable = 'x_Al'
       sort_by = id
       execute_on = 'initial timestep_end'
    [../]

    [./xCr]
      type = LineValueSampler
      start_point = '-3000 0 0'
      end_point = '3000 0 0'
      num_points = 1501     
      variable = 'x_Cr'
      sort_by = id
      execute_on  = 'initial timestep_end'
    [../]

    [./xFe]
       type = LineValueSampler
       start_point = '-3000 0 0'
       end_point = '3000 0 0'
       num_points = 1501
       variable = 'x_Fe'
       sort_by = id
       execute_on = 'initial timestep_end'
    [../]
[]

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
  l_max_its = 50
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-10
  nl_max_its = 50

  #since tc = 3.8171 sec
  #10 hrs = 9.43124361E3
  #1  day = 2.26349847E4
  #100days= 2.26349847E6

  start_time = 0
  end_time = 2.26349847E4  #non-dimensional
 
  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    percent_change = 0.1
    dt = 100.0 #real_dt = 0.1sec
  [../]
  
[]

[Outputs]

  file_base = AlCrNiFe_twophase_new
  #exodus = true

  [.exodus]
    type = Exodus
    execute_on = 'initial timestep_end'
    interval = 50
  [../]
  
  [./csv]
    type = CSV
    execute_on = 'initial final'
    execute_postprocessors_on = 'initial timestep_end final'
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
