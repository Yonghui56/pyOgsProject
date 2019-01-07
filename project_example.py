#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example script creating the project file for the hc-goswami benchmark.
@date: 2018-Today
@author: jasper.bathmann@ufz.de
"""
import pyOgsProject

mesh_name = "goswami_input"
p_ini_name = "p_ini"
c_ini_name = "c_ini"
g = 9.81
project_dir = "./"
project_name = "goswami_benchmark"
timerepeats = [50, 595, 1140, 3600]
timedeltaTs = [0.01, 0.1, 1, 1]
outputrepeats = [1, 1, 6]
outputdeltaN = [1185, 600, 600]
project = pyOgsProject.GenerateProject(project_dir + project_name + ".prj")
project.setMesh(mesh_name+".vtu")
project.setStandardProcessInformation()
project.setStandartTimeLoop()
project.setStandartParameters()
project.setStandardDensityModel()

project.setStandardNonlinearSolvers()

project.setTimeSteppingAndOutputLoops(
        timerepeats, timedeltaTs, outputrepeats, outputdeltaN)


project.convergence_criterion_reltols = "1e-6 1e-6"

# Parameter values are given as a list. They could either be set directly by
# calling project.setParameters(self, *args) with *args being a list of the
# form [[parameter_names...], [parameter_types...], [parameter_values...]], or
# by using project.setStandardParameters() and and modifiing
# project.parameter_values = [rho_v, Dm_v, retardation_v, decay_v, beta_l_v,
# beta_t_v, c_ini_v, p_ini_v, constant_porosity_parameter_v, kappa1_v] as done
# below

project.parameter_values[1] = "0"
project.parameter_values[4] = "5e-3"
project.parameter_values[5] = "0.0005"
project.parameter_values[-2] = "0.385"
project.parameter_values[-1] = "1.2388e-9 0 0 1.2388e-9"
project.resetBoundaryConditions()

project.createBoundaryCondition(
        "c", "NonuniformDirichlet", "goswami_input_leftBoundary", "c_ini")
project.createBoundaryCondition(
        "c", "NonuniformDirichlet", "goswami_input_rightBoundary", "c_ini")
project.createBoundaryCondition(
        "p", "NonuniformDirichlet", "goswami_input_leftBoundary", "p_ini")
project.createBoundaryCondition(
        "p", "NonuniformDirichlet", "goswami_input_rightBoundary", "p_ini")

project.resetInitialConditions(p_ini_name, c_ini_name)
project.processspeci_bo_force = "0 0 -"+str(g)
project.writeProjectFile()
