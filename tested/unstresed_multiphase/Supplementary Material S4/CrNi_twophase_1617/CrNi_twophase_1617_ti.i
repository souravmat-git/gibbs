#This kernel uses the TabulatedConjugatePhaseData
#Material to solve the equations. 
#This is a model Cr-Ni alloy problem
#beta phase- FCC_A1
#alpha phase- BCC_A2
#T = 1617K
#Simulation details:
#Interfacial energy = 0.5 J/m2
#interfacial width = 0.5e-6 m
#Length = 100e-6 m
#characteristic_length = (0.1e-6/6.0)m = 0.01667e-6 m
#characteristic_time = 0.0015 s
#Ni_diff_pot_eqm = -2.5490e4 J/mol 
#Nondimensional length = (Length/characteristic_length)
#No. of finite elements = 6.0(Length/interface_width)

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1200
  xmin =-3000
  xmax = 3000
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

[./x_Ni]
  order = FIRST
  family = LAGRANGE
[../]

[./mu_Ni]
  order = FIRST
  family = LAGRANGE
[../]

[]

[AuxVariables]

 [./delta_omega]
   order = FIRST
   family = MONOMIAL
 [../]
 
 [./total_energy]
  order = FIRST
  family = MONOMIAL
 [../]
 
 [./bulk_energy]
  order = FIRST
  family = MONOMIAL
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
    variable = 'mu_Ni'
    boundary = 'left right'
    value = 0
  [../]
 
[]


[ICs]
  
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
  
  [./ic_func_xNi]
   type = BoundingBoxIC
   variable = x_Ni
   x1 = -3000
   y1 = 0.0
   x2 = 0.0
   y2 = 0.0
   inside = 0.465        #xNi_A1 = 0.4998720
   outside = 0.3618972   #xNi_A2 = 0.3618972    
  [../]
  
  [./mu_IC]
   variable = mu_Ni
   type = ConstantIC
   value = -1.89
  [../]
  
[]

[AuxKernels]

  [./Cr_chemicalpotential]
    type = ChemicalPotential
    variable = delta_omega
    dh = dhbeta_dphialpha
  [../]
  
  [./bulk_energy]
    type = BulkEnergy
    variable = bulk_energy
    h_alpha = h_alpha
    h_beta = h_beta
    f_alpha = A_chem_pot_alpha
    f_beta = A_chem_pot_beta
  [../]
  
  [./total_energy]
    type = TotalEnergy
    variable = total_energy
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    barrier_height = m
    h_alpha = h_alpha
    h_beta = h_beta
    f_alpha = A_chem_pot_alpha
    f_beta = A_chem_pot_beta
  [../] 
[]

[Kernels]
  
  #Variable that this kernel operates on
  #is the diffusion potential of comp B
  
  [./PhaseConc]
    type = MultiPhaseConstraintMu
    variable = mu_Ni
    xB = x_Ni
    phase_alpha = phi_alpha
    phase_beta = phi_beta
  [../]

  ######################
  # Diffusion Equation #
  ######################
    
  [./Ni_balance]
    type = BinaryMultiPhaseMassBalance
    variable = x_Ni
    B_diff_pot = mu_Ni
    phase_alpha = phi_alpha
    phase_beta = phi_beta
    h_alpha = h_alpha
    h_beta = h_beta
  [../]

  [./xNi_dot]
    type = TimeDerivative
    variable = x_Ni
  [../]
  
  #######################
  # Allen-Cahn equation #
  # for phase-alpha     #
  #######################
  
  [./BinaryMultiPhaseDrivingForce_alpha]
    type = BinaryMultiPhaseDrivingForce
    variable = phi_alpha     #Phase_1
    phase_2 = phi_beta       #Coupled phase_2
    B_diff_pot = mu_Ni
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
    B_diff_pot = mu_Ni
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
     prop_names =  'L_phi      kappa     m     nd_factor'
     prop_values = '0.0045     112.5    1.0    29.8750'
  [../]  
  
  [./ZeroFieldProperties]
    type = GenericConstantMaterial
    prop_names = 'L_BB_gamma L_BB_delta L_BB_epsilon
                  dL_BB_muB_gamma dL_BB_muB_delta dL_BB_muB_epsilon'
    prop_values = '0 0 0 
                   0 0 0' 
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
    B_diff_pot = mu_Ni
    molar_volume = 1.0
    table_object = table_FCC_data
    #outputs = exodus
 [../] 
 
 [./BinaryConjugatePhaseMaterial_BCC]
    type = BinaryConjugatePhaseMaterial
    A_chem_pot = A_chem_pot_alpha
    xB = xB_alpha
    inv_B_tf = inv_B_tf_alpha
    B_diff_pot = mu_Ni
    molar_volume = 1.0
    table_object = table_BCC_data
    #outputs = exodus
  [../]
  
  [./BinaryConjugateKineticMaterialFCC]
    type = BinaryConjugateKineticMaterial
    L_BB = L_BB_beta
    dL_BB_muB = dL_BB_muB_beta
    B_diff_pot = mu_Ni
    table_object = table_FCC_mobility_data
    #outputs = exodus
  [../]
  
  [./BinaryConjugateKineticMaterialBCC]
    type = BinaryConjugateKineticMaterial
    L_BB = L_BB_alpha
    dL_BB_muB = dL_BB_muB_alpha
    B_diff_pot = mu_Ni
    table_object = table_BCC_mobility_data
    #outputs = exodus
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
    table_name = CrNi_data/CrNi_invert_data_FCC_1617.csv
  [../]
  
  [./table_BCC_data]
    type = BinaryConjugatePhaseData
    table_name = CrNi_data/CrNi_invert_data_BCC_1617.csv
  [../] 
  
  [./table_FCC_mobility_data]
    type = BinaryConjugateMobilityData
    table_name = CrNi_data/CrNi_mobility_data_FCC_invert_1617.csv
  [../]
  
  [./table_BCC_mobility_data]
    type = BinaryConjugateMobilityData
    table_name = CrNi_data/CrNi_mobility_data_BCC_invert_1617.csv
  [../] 
 
[]


[Postprocessors]

#After each time step collect the data
#TIMSTEP_END
#TIMESTEP_BEGIN
 [./volume]
   type = VolumePostprocessor
   execute_on = 'initial timestep_end'
 [../]
 
 [./computation_time]
    type = PerfGraphData
    section_name = Transient::PicardSolve
    data_type = total
    execute_on = 'timestep_end'
 [../]
 
 [./integrate_bulk_energy]
  type = ElementIntegralVariablePostprocessor
  variable = bulk_energy
  execute_on = 'initial timestep_end'
 [../]
 
 [./integrate_total_energy]
  type = ElementIntegralVariablePostprocessor 
  variable = total_energy
  execute_on = 'initial timestep_end'
 [../]
 
 [./pos_alpha]
   type = FindValueOnLine
   target = 0.5
   v = phi_beta
   start_point = '-3000 0 0'
   end_point   = '3000 0 0'
   tol = 1e-8
   execute_on = 'initial timestep_end'
  [../]
  
  [./pos_beta]
   type = FindValueOnLine
   target = 0.5
   v = phi_alpha
   start_point = '-3000 0 0'
   end_point =   ' 3000 0 0'
   tol = 1e-8
   execute_on = 'initial timestep_end'
  [../]
  
  [./volumefraction_alpha]
    type = VolumeFraction
    variable = phi_alpha
    execute_on = 'initial timestep_end'
  [../]

 [./Chemical_potential_Cr_alpha]
    type = ElementIntegralMaterialProperty 
    mat_prop = A_chem_pot_alpha
    execute_on = 'initial timestep_end'
  [../]
	
  [./Chemical_potential_Cr_beta]
    type = ElementIntegralMaterialProperty
    mat_prop = A_chem_pot_beta 
    execute_on = 'initial timestep_end'
  [../]

  [./Diffusion_potential_Ni] 
   type = ElementIntegralVariablePostprocessor
   variable = mu_Ni
   execute_on = 'initial timestep_end'
 [../] 
 
 [./tf_NiNi_alpha]
  type = ElementIntegralMaterialProperty
  mat_prop = inv_B_tf_alpha
  execute_on = 'initial timestep_end'
[../]

 [./tf_NiNi_beta]
  type = ElementIntegralMaterialProperty
  mat_prop = inv_B_tf_beta
  execute_on = 'initial timestep_end'
 [../]

 [./L_NiNi_alpha]
  type = ElementIntegralMaterialProperty
  mat_prop = L_BB_alpha
  execute_on = 'initial timestep_end'
 [../]

 [./L_NiNi_beta]
  type = ElementIntegralMaterialProperty
  mat_prop = L_BB_beta
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
  petsc_options_value = 'asm       lu          nonzero'
  
  l_tol = 1.0e-5
  l_max_its = 15
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-10
  nl_max_its = 30
  
  #Real time = 1s
  #Non-dimensional time = 2.62780105e1

  end_time = 6e7 #non-dimensional
  #dt = 1
  
  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    percent_change = 0.1
    dt = 100.0  #real_dt = 0.1sec
  [../]
  
[]

[Outputs]

  file_base = CrNi_twophase_1617_ti
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
