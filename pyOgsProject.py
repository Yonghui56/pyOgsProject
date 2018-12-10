"""/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */"""
"""

@author: jasper.bathmann@ufz.de
"""
class BoundaryCondition:     
    def __init__(self,*args):
        if(args[0]=="Neumann" or args[0] == "Dirichlet"): self.NeumannOrDirichletBoundary(*args)
        elif(args[0]=="NonuniformNeumann" or args[0] == "NonuniformDirichlet"):
            self.NonuniformBoundary(*args)
        elif(args[0]=="NeumannTimeDependant"):
            self.NeumannTimeDependant(*args)
        elif(args[0]=="NonuniformVariableDependentNeumann"):
            self.NonuniformVariableDependentNeumann(*args)
        else: print("Error: unknown boundary condition!")
    
    
    def NeumannOrDirichletBoundary(self, tYpe, geoset, geometry, parameter):
        self.geoset = geoset
        self.geometry = geometry
        self.type = tYpe
        self.parameter = parameter
        
    def NonuniformBoundary(self, tYpe, geoset, geometry, meshname, field_name):
        self.geoset = geoset
        self.geometry =geometry
        self.type = tYpe
        self.mesh = meshname
        self.field_name = field_name

        
    def NeumannTimeDependant(self, tYpe, geoset, geometry, meshname, psi_leaf, g_contri, resistance):
        self.geoset = geoset
        self.geometry = geometry
        self.type = tYpe
        self.mesh = meshname
        self.psi_leaf = psi_leaf
        self.g_contri = g_contri
        self.resistance = resistance
        
    def NonuniformVariableDependentNeumann(self, tYpe, meshname, constant, prefac1, prefac2, prefac3):
        self.type = tYpe
        self.mesh = meshname
        self.constant = constant
        self.prefac1 = prefac1
        self.prefac2 = prefac2
        self.prefac3 = prefac3
        
    def writeBoundary(self,file):
        if(self.type=="Neumann" or self.type == "Dirichlet"):
            self.writeNeumannOrDirichletBoundary(file)
        if(self.type=="NonuniformNeumann" or self.type == "NonuniformDirichlet"):
            self.writeNonuniformBoundary(file)
        if(self.type=="NonuniformVariableDependentNeumann"):
            self.writeNeumannTimeDependentBoundary(file) 
            
    def writeNeumannOrDirichletBoundary(self, file):
        file.write('                <boundary_condition>\n')    
        file.write('                    <geometrical_set>'+self.geoset+'</geometrical_set>\n')
        file.write('                    <geometry>'+ self.geometry+'</geometry>\n')
        file.write('                    <type>'+self.type+'</type>\n')
        file.write('                    <parameter>'+self.parameter+'</parameter>\n')
        file.write('                </boundary_condition>\n')
    
    def writeNonuniformBoundary(self, file):
        file.write('                <boundary_condition>\n')    
        file.write('                    <type>'+self.type+'</type>\n')
        file.write('                    <mesh>'+self.mesh+'</mesh>\n')
        file.write('                    <field_name>'+self.field_name+'</field_name>\n')
        file.write('                </boundary_condition>\n')
        
        
    def writeNeumannTimeDependentBoundary(self, file):
        file.write('                <boundary_condition>\n')    
        file.write('                    <type>'+self.type+'</type>\n')
        file.write('                    <mesh>'+self.mesh+'</mesh>\n')
        file.write('                    <constant_name>'+self.constant+'</constant_name>\n')
        file.write('                    <coefficient_current_variable_name>'+self.prefac1+'</coefficient_current_variable_name>\n')
        file.write('                    <coefficient_other_variable_name>'+self.prefac2+'</coefficient_other_variable_name>\n')
        file.write('                    <coefficient_mixed_variables_name>'+self.prefac3+'</coefficient_mixed_variables_name>\n')
        file.write('                </boundary_condition>\n')


class GenerateProject:
    def __init__(self, projectname):
        self.name = projectname
        self.boundary_conditions_c, self.boundary_conditions_p = [], []

    def writeProjectFile(self):
        self.file = open(self.name,"w+")
        self.file.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
        self.file.write('<OpenGeoSysProject>\n')
        self.writeMeshAndGeometry()
        self.writeProcessInformation()
        self.writeTimeLoop()
        self.writeParameters()
        self.writeBoundaryConditions()
        self.writeNonlinearSolvers()
        self.file.write('</OpenGeoSysProject>\n')
        self.file.close()
        
    def setMesh(self,meshname):
        self.meshname = meshname
        self.meshnames = []
        self.meshnames.append(meshname)
        
    def writeMeshAndGeometry(self):
        self.file.write('    <meshes>\n')
        for meshname in self.meshnames:
            self.file.write('        <mesh>'+meshname+'</mesh>\n')
        self.file.write('    </meshes>\n')
        
    def setProcessInformation(self, *args):
        self.processname = str(args[0])
        self.processtype = str(args[1])
        self.processvariable_concentration = str(args[2])
        self.processvariable_pressure = str(args[3])
        self.processfluid_density = str(args[4])
        self.processfluid_viscosity = str(args[5])
        self.processmedium_permiablility = str(args[6])
        self.processmedium_porosity = str(args[7])
        self.processmedium_storage = str(args[8])
        self.processfluid_ref_density = str(args[9])
        self.processsol_disper_long = str(args[10])
        self.processsol_disper_tran = str(args[11])
        self.processmol_diff_coeff = str(args[12])
        self.processretard_fac =  str(args[13])
        self.processdecay_rate = str(args[14])
        self.processspeci_bo_force = str(args[15])
        self.processsecond_var_type = str(args[16])
        self.processsecond_var_name= str(args[17])
        self.processsecond_var_out_name = str(args[18])
        
    def setStandardProcessInformation(self):
        self.processname = "hc"
        self.processtype = "ComponentTransport"
        self.processvariable_concentration = "concentration"
        self.processvariable_pressure = "pressure"
        self.processfluid_density = "1e3"
        self.processfluid_viscosity = "1.0e-3"
        self.processmedium_permiablility = "kappa1"
        self.processmedium_porosity = "constant_porosity_parameter"
        self.processmedium_storage = "0"
        self.processfluid_ref_density = "rho_fluid"
        self.processsol_disper_long = "beta_l"
        self.processsol_disper_tran = "beta_t"
        self.processmol_diff_coeff = "Dm"
        self.processretard_fac =  "retardation"
        self.processdecay_rate = "decay"
        self.processspeci_bo_force = "-9.81"
        self.processsecond_var_type = 'static'
        self.processsecond_var_name= 'darcy_velocity'
        self.processsecond_var_out_name = 'darcy_velocity'
        
    def setStandardDensityModel(self):
        self.densityModel = "ConcentrationDependent"

        
    def writeProcessInformation(self):
        self.file.write('    <processes>\n')
        self.file.write('        <process>\n')
        self.file.write('            <name>'+self.processname+'</name>\n')
        self.file.write('            <type>'+self.processtype+'</type>\n')
        self.file.write('            <integration_order>2</integration_order>\n')
        self.file.write('            <process_variables>\n')
        self.file.write('                <concentration>'+self.processvariable_concentration+'</concentration>\n')
        self.file.write('                <pressure>'+self.processvariable_pressure+'</pressure>\n')
        self.file.write('            </process_variables>\n')
        self.file.write('            <fluid>\n')   
        self.file.write('                <density>\n')
        if(self.densityModel == "Constant"):
            self.file.write('                    <type>Constant</type>\n')
            self.file.write('                    <value>1000</value>\n')
        elif(self.densityModel == "ConcentrationDependent"):
            self.file.write('                    <type>ConcentrationDependent</type>\n')
            self.file.write('                    <reference_density>'+self.processfluid_density+'</reference_density>\n')
            self.file.write('                    <reference_concentration>'+"0"+'</reference_concentration>\n')
            
            self.file.write('                    <fluid_density_difference_ratio>'+"0.701"+'</fluid_density_difference_ratio>\n')#Usually0.825 See Millero et al, density of seawater.... Ocen Science 2009
        elif(self.densityModel == "ConcentrationAndPressureDependent"):
            self.file.write('                    <type>ConcentrationAndPressureDependent</type>\n')
            self.file.write('                    <reference_density>'+self.processfluid_density+'</reference_density>\n')
            self.file.write('                    <reference_concentration>'+"0"+'</reference_concentration>\n')
            self.file.write('                    <fluid_density_concentration_difference_ratio>'+"0.0"+'</fluid_density_concentration_difference_ratio>\n')#Usually0.825#See Millero et al, density of seawater.... Ocen Science 2009
            self.file.write('                    <reference_pressure>'+"0"+'</reference_pressure>\n')
            drho_dp = float(self.processmedium_storage)/float(self.parameter_values[-2])/float(self.processfluid_density)
            self.file.write('                    <fluid_density_pressure_difference_ratio>'+drho_dp+'</fluid_density_pressure_difference_ratio>\n')#See Millero et al, density of seawater.... Ocen Science 2009
        self.file.write('                </density>\n')
        self.file.write('                <viscosity>\n')
        self.file.write('                    <type>Constant</type>\n')
        self.file.write('                    <value>'+self.processfluid_viscosity+'</value>\n')
        self.file.write('                </viscosity>\n')
        self.file.write('            </fluid>\n')   
        self.file.write('            <porous_medium>\n') 
        self.file.write('                <porous_medium id="0">\n')
        self.file.write('                <permeability>\n')
        self.file.write('                    <permeability_tensor_entries>'+self.processmedium_permiablility+'</permeability_tensor_entries>\n')
        self.file.write('                    <type>Constant</type>\n')
        self.file.write('                </permeability>\n')
        self.file.write('                <porosity>\n')
        self.file.write('                    <type>Constant</type>\n')
        self.file.write('                    <porosity_parameter>'+self.processmedium_porosity+'</porosity_parameter>\n')
        self.file.write('                </porosity>\n')
        self.file.write('                <storage>\n')
        self.file.write('                    <type>Constant</type>\n')
        self.file.write('                    <value>'+self.processmedium_storage+'</value>\n')
        self.file.write('                </storage>\n')
        self.file.write('                </porous_medium>\n')
        self.file.write('            </porous_medium>\n') 
        self.file.write('            <fluid_reference_density>'+self.processfluid_ref_density+'</fluid_reference_density>\n')
        self.file.write('            <solute_dispersivity_longitudinal>'+self.processsol_disper_long+'</solute_dispersivity_longitudinal>\n')
        self.file.write('            <solute_dispersivity_transverse>'+self.processsol_disper_tran+'</solute_dispersivity_transverse>\n')
        self.file.write('            <molecular_diffusion_coefficient>'+self.processmol_diff_coeff+'</molecular_diffusion_coefficient>\n')
        self.file.write('            <retardation_factor>'+self.processretard_fac+'</retardation_factor>\n')
        self.file.write('            <decay_rate>'+self.processdecay_rate+'</decay_rate>\n')
        self.file.write('            <specific_body_force>'+ self.processspeci_bo_force+'</specific_body_force>\n')
        self.file.write('            <secondary_variables>\n')
        self.file.write('                <secondary_variable type="'+self.processsecond_var_type+'" internal_name="'+self.processsecond_var_name+'" output_name="'+self.processsecond_var_out_name+'"/>\n')
        self.file.write('            </secondary_variables>\n')
        self.file.write('        </process>\n')
        self.file.write('    </processes>\n')

    def setTimeLoop(self,*args):
        self.process_ref = str(args[0])
        self.nonlinear_solver = str(args[1])
        self.convergence_criterion_type = str(args[2])
        self.convergence_criterion_norm = str(args[3])
        self.convergence_criterion_reltols = str(args[4])
        self.time_disc_type = str(args[5])
        self.outputvariables =  args[6]
        self.time_stepping_type = str(args[7])
        self.time_stepping_t_ini = str(args[8])
        self.time_stepping_t_end = str(args[9])
        self.time_stepping_delta_t = str(args[10])
        self.output_type = str(args[11])
        self.output_prefix = str(args[12])
        self.output_each_steps = str(args[13])
        
    def setStandartTimeLoop(self):
        self.process_ref = "hc"
        self.nonlinear_solver = "basic_picard"
        self.convergence_criterion_type = "PerComponentDeltaX"
        self.convergence_criterion_norm = "NORM2"
        self.convergence_criterion_reltols = "5e-5 5e-5"
        self.time_disc_type = "BackwardEuler"
        self.outputvariables =  ["concentration","pressure","darcy_velocity"]
        self.time_stepping_type = "FixedTimeStepping"
        self.time_stepping_t_ini = "0.0"
        self.time_stepping_t_end = "1e5"
        self.output_type = "VTK"
        self.output_prefix = "interface_output"
        self.output_each_steps = "1"
        
    def setTimeSteppingAndOutputLoops(self,timerepeats, timedeltaTs, outputrepeats, outputdeltaN):
        self.timeSteppingSettings = []
        self.outputSteppingSettings = []
        if(len(timerepeats)!=len(timedeltaTs)!=0):
            print("INPUT ERROR: All input lists for timeloop are obliged to have equal length>0!")
        else:
            for i in range(len(timerepeats)):
                timestep = []
                timestep.append('                        <pair>\n')
                timestep.append('                            <repeat>' + str(timerepeats[i])+ '</repeat>\n')
                timestep.append('                            <delta_t>' + str(timedeltaTs[i])+ '</delta_t>\n')#changed from 1
                timestep.append('                        </pair>\n')
                self.timeSteppingSettings.append(timestep)

        if(len(outputrepeats)!=len(outputdeltaN)!=0):
            print("INPUT ERROR: All input lists for outputloop are obliged to have equal length>0!")
        else:
            for i in range(len(outputrepeats)):
                outputstep = []
                outputstep.append('                <pair>\n')
                outputstep.append('                    <repeat>' + str(outputrepeats[i])+ '</repeat>\n')
                outputstep.append('                    <each_steps>' + str(outputdeltaN[i])+ '</each_steps>\n')
                outputstep.append('                </pair>\n')
                self.outputSteppingSettings.append(outputstep)

    def writeTimeLoop(self):
        self.file.write('    <time_loop>\n')
        self.file.write('        <processes>\n')
        self.file.write('            <process ref="'+self.process_ref+'">\n')
        self.file.write('                <nonlinear_solver>'+self.nonlinear_solver+'</nonlinear_solver>\n')
        self.file.write('                <convergence_criterion>\n')
        self.file.write('                    <type>'+self.convergence_criterion_type+'</type>\n')
        self.file.write('                    <norm_type>'+self.convergence_criterion_norm+'</norm_type>\n')
        self.file.write('                    <reltols>'+self.convergence_criterion_reltols+'</reltols>\n')
        self.file.write('                </convergence_criterion>\n')
        self.file.write('                <time_discretization>\n')
        self.file.write('                    <type>'+self.time_disc_type+'</type>\n')
        self.file.write('                </time_discretization>\n')
        self.file.write('                <time_stepping>\n')
        """self.file.write('                    <type>EvolutionaryPIDcontroller</type>')
        self.file.write('            <t_initial> 0.0 </t_initial>')
        self.file.write('            <t_end> 96000 </t_end>')
        self.file.write('            <dt_guess> .5 </dt_guess>')
        self.file.write('            <dt_min> 1e-6 </dt_min>')
        self.file.write('            <dt_max> .5 </dt_max>')
        self.file.write('            <rel_dt_min> 0.1 </rel_dt_min>')
        self.file.write('            <rel_dt_max> 10 </rel_dt_max>')
        self.file.write('            <tol> 1.0 </tol>')
        self.file.write('            <fixed_output_times>')
        for i in range(500):
            self.file.write(' ' + str((1+i)*30))
            #300 600 900 1200 1500 1800 2100 2400 2700 3000 3300 3600 3900 4200 4500 4800 5100 5400 5700 6000 6300 6600 6900 7200 7500 7800 8100 8400 8700 9000
        self.file.write('  </fixed_output_times>')
        self.file.write('        </time_stepping>')
        """
        self.file.write('                    <type>'+self.time_stepping_type+'</type>\n')
        self.file.write('                    <t_initial>'+self.time_stepping_t_ini+'</t_initial>\n')
        self.file.write('                    <t_end>'+self.time_stepping_t_end+'</t_end>\n')
        self.file.write('                    <timesteps>\n')
        for i in range(len(self.timeSteppingSettings)):
            
            for j in range(len(self.timeSteppingSettings[i])):
                self.file.write(self.timeSteppingSettings[i][j])

        
        self.file.write('                    </timesteps>\n')
        self.file.write('                </time_stepping>\n')
        
        self.file.write('            </process>\n')
        self.file.write('        </processes>\n')
        self.file.write('        <output>\n')
        self.file.write('            <type>'+self.output_type +'</type>\n')
        self.file.write('            <prefix>'+self.output_prefix +'</prefix>\n')
        self.file.write('            <timesteps>\n')
        for i in range(len(self.outputSteppingSettings)):
            
            for j in range(len(self.outputSteppingSettings[i])):
                self.file.write(self.outputSteppingSettings[i][j])

        
        self.file.write('            </timesteps>\n')
        self.file.write('            <variables>\n')
        for variable in self.outputvariables:
            self.file.write('                <variable>'+variable+'</variable>\n')
        self.file.write('            </variables>\n')
        self.file.write('        </output>\n')
        self.file.write('    </time_loop>\n')
        
    def setParameters(self,*args):
        n_para=int(len(args)/3)
        self.parameter_names=[]
        self.parameter_types=[]
        self.parameter_values=[]
        
        for i in range(n_para):
            self.parameter_names.append(args[3*i])
            self.parameter_types.append(args[3*i]+1)
            self.parameter_values.append(args[3*i]+2)
        
    def setStandartParameters(self):
        rho_n = "rho_fluid"
        rho_t = "Constant"
        rho_v = "1000"
        Dm_n = "Dm"
        Dm_t = "Constant"
        Dm_v ="2e-9"
        retardation_n = "retardation"
        retardation_t = "Constant"
        retardation_v = "1"
        decay_n = "decay"
        decay_t = "Constant"
        decay_v = "0"
        beta_l_n ="beta_l"
        beta_l_t = "Constant"
        beta_l_v = "1"
        beta_t_n = "beta_t"
        beta_t_t = "Constant"
        beta_t_v = "0.1"
        c_ini_n = "c_ini"
        c_ini_t = "MeshNode"
        c_ini_v = "c_ini"
        p_ini_n = "p_ini"
        p_ini_t = "MeshNode"
        p_ini_v = "p_ini"
        constant_porosity_parameter_n ="constant_porosity_parameter"
        constant_porosity_parameter_t = "Constant"
        constant_porosity_parameter_v = "0.01"
        kappa1_n = "kappa1"
        kappa1_t = "Constant"
        kappa1_v = "1.239e-11 0 0 0 1.239e-11 0 0 0 1.239e-11"
        self.parameter_names = [rho_n, Dm_n, retardation_n, decay_n, beta_l_n, beta_t_n, c_ini_n, 
                                p_ini_n, constant_porosity_parameter_n, kappa1_n]
        self.parameter_types = [rho_t, Dm_t, retardation_t, decay_t, beta_l_t, beta_t_t, c_ini_t, 
                                p_ini_t, constant_porosity_parameter_t, kappa1_t]
        self.parameter_values = [rho_v, Dm_v, retardation_v, decay_v, beta_l_v, beta_t_v, c_ini_v, 
                                 p_ini_v, constant_porosity_parameter_v, kappa1_v]

    def writeParameters(self):
        self.file.write('    <parameters>\n')
        n_para = len(self.parameter_names)
        for i in range(n_para):
            self.file.write('        <parameter>\n')
            self.file.write('            <name>'+self.parameter_names[i]+'</name>\n')
            self.file.write('            <type>'+self.parameter_types[i]+'</type>\n')
            if self.parameter_types[i]=="MeshNode":
                self.file.write('            <field_name>'+self.parameter_values[i]+'</field_name>\n')
            if self.parameter_types[i]=="Constant":
                
                if(len(self.parameter_values[i].split())>1):
                    self.file.write('            <values>'+self.parameter_values[i]+'</values>\n')
                else:
                    self.file.write('            <value>'+self.parameter_values[i]+'</value>\n')
                
            self.file.write('        </parameter>\n')
        self.file.write('    </parameters>\n')
    
    def resetInitialConditions(self,p_ini,c_ini):
        self.c_ini = c_ini
        self.p_ini = p_ini
    
    
    def resetBoundaryConditions(self):
        self.boundary_conditions_c, self.boundary_conditions_p = [], []
    
    def createBoundaryCondition(self,corresponding_parameter,*args):
        bc = BoundaryCondition(*args)
        if(corresponding_parameter=="c"):
            self.boundary_conditions_c.append(bc)
        if(corresponding_parameter=="p"):
            self.boundary_conditions_p.append(bc)
        if (bc.mesh+".vtu") not in self.meshnames:
            self.meshnames.append(bc.mesh+".vtu")
        
    def writeBoundaryConditions(self):
        self.file.write("    <process_variables>\n")
        self.file.write("        <process_variable>\n")
        self.file.write("            <name>"+self.processvariable_concentration+"</name>\n")
        self.file.write("            <components>1</components>\n")
        self.file.write("            <order>1</order>\n")
        self.file.write("            <initial_condition>"+self.c_ini+"</initial_condition>\n")
        self.file.write("            <boundary_conditions>\n")
        for boundary_condition in self.boundary_conditions_c:
            boundary_condition.writeBoundary(self.file)
        self.file.write("            </boundary_conditions>\n")
        self.file.write("        </process_variable>\n")
        self.file.write("        <process_variable>\n")
        self.file.write("            <name>"+self.processvariable_pressure+"</name>\n")
        self.file.write("            <components>1</components>\n")
        self.file.write("            <order>1</order>\n")
        self.file.write("            <initial_condition>"+self.p_ini+"</initial_condition>\n")
        self.file.write("            <boundary_conditions>\n")
        for boundary_condition in self.boundary_conditions_p:
            boundary_condition.writeBoundary(self.file)
        self.file.write("            </boundary_conditions>\n")
        self.file.write("        </process_variable>\n")
        self.file.write("    </process_variables>\n")

    def setNonlinearSolvers(self, *args):
        self.solver_name = str(args[0])
        self.solver_type = str(args[1])
        self.solver_max_iter = str(args[2])
        self.solver_linear_solver = str(args[3])
        self.solver_linear_solver_name = str(args[4])
        self.solver_linear_solver_lis = str(args[5])
        self.solver_linear_solver_eigen_type = str(args[6])
        self.solver_linear_solver_eigen_precon = str(args[7])
        self.solver_linear_solver_eigen_max_iter = str(args[8])
        self.solver_linear_solver_eigen_error_tol = str(args[9])
        self.solver_linear_solver_petsc_prefix = str(args[10])
        self.solver_linear_solver_petsc_parameters = str(args[11])
        
    def setStandardNonlinearSolvers(self):
        self.solver_name = "basic_picard"
        self.solver_type = "Picard"
        self.solver_max_iter = "200"
        self.solver_linear_solver = "general_linear_solver"
        self.solver_linear_solver_name = "general_linear_solver"
        self.solver_linear_solver_lis = "-i bicgstab -p ilut -tol 1e-8 -maxiter 20000"
        self.solver_linear_solver_eigen_type = "BiCGSTAB"
        self.solver_linear_solver_eigen_precon = "DIAGONAL"
        self.solver_linear_solver_eigen_max_iter = "10000"
        self.solver_linear_solver_eigen_error_tol = "1e-12"
        self.solver_linear_solver_petsc_prefix = "hc"
        self.solver_linear_solver_petsc_parameters = "-hc_ksp_type bcgs -hc_pc_type bjacobi -hc_ksp_rtol 1e-12 -hc_ksp_max_it 20000"
        
    def setILUTNonlinearSolvers(self):
        self.solver_name = "basic_picard"
        self.solver_type = "Picard"
        self.solver_max_iter = "300"
        self.solver_linear_solver = "general_linear_solver"
        self.solver_linear_solver_name = "general_linear_solver"
        self.solver_linear_solver_lis = "-i bicgstab -p ilut -tol 1e-8 -maxiter 20000"
        self.solver_linear_solver_eigen_type = "BiCGSTAB"
        self.solver_linear_solver_eigen_precon = "ILUT"
        self.solver_linear_solver_eigen_max_iter = "5000"
        self.solver_linear_solver_eigen_error_tol = "1e-12"
        self.solver_linear_solver_petsc_prefix = "hc"
        self.solver_linear_solver_petsc_parameters = "-hc_ksp_type bcgs -hc_pc_type bjacobi -hc_ksp_rtol 1e-8 -hc_ksp_max_it 20000"
            
        
    def writeNonlinearSolvers(self):
        self.file.write("    <nonlinear_solvers>\n")
        self.file.write("        <nonlinear_solver>\n")
        self.file.write("            <name>"+self.solver_name+"</name>\n")
        self.file.write("            <type>"+self.solver_type+"</type>\n")
        self.file.write("            <max_iter>"+self.solver_max_iter+"</max_iter>\n")
        self.file.write("            <linear_solver>"+self.solver_linear_solver+"</linear_solver>\n")
        self.file.write("        </nonlinear_solver>\n")
        self.file.write("        </nonlinear_solvers>\n")
        self.file.write("    <linear_solvers>\n")
        self.file.write("        <linear_solver>\n")
        self.file.write("            <name>"+self.solver_linear_solver_name+"</name>\n")
        self.file.write("            <lis>"+self.solver_linear_solver_lis+"</lis>\n")
        self.file.write("            <eigen>\n")
        self.file.write("                <solver_type>"+self.solver_linear_solver_eigen_type+"</solver_type>\n")
        self.file.write("                <precon_type>"+self.solver_linear_solver_eigen_precon+"</precon_type>\n")
        self.file.write("                <max_iteration_step>"+self.solver_linear_solver_eigen_max_iter+"</max_iteration_step>\n")
        self.file.write("                <error_tolerance>"+self.solver_linear_solver_eigen_error_tol+"</error_tolerance>\n")
        self.file.write("            </eigen>\n")
        self.file.write("            <petsc>\n")
        self.file.write("                <prefix>"+self.solver_linear_solver_petsc_prefix+"</prefix>\n")
        self.file.write("                <parameters>"+self.solver_linear_solver_petsc_parameters+"</parameters>\n")
        self.file.write("            </petsc>\n")
        self.file.write("        </linear_solver>\n")
        self.file.write("    </linear_solvers>\n")
        
        


