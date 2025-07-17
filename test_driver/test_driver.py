#!/usr/bin/env python3

import sys
import os
import io
import operator
import numpy as np

# Plotting libs
import matplotlib
import pylab as plt
from matplotlib import cm

from .make_lammps_input import (
    setup_problem,
    make_stack_twin_test,
    make_refine_us,
    make_refine_ut,
    make_gammasurface_moves,
    compute_eq_latconst,
)

# for Crystal Genome
from kim_tools import SingleCrystalTestDriver

#for timer
import time

matplotlib.use("Agg")  # Use backend for non-interactive plotting

class TestDriver(SingleCrystalTestDriver):
    """Intrinsic-extrinsic stacking fault and gamma surface calculation for FCC lattice
    """

    
    def _calculate(self, 
                   pressure: float = 0.0,
                   Num_layers_gamma_surf: int = 14,
                   compute_gamma_surf: bool = True,
                   conv_cutoff: float = 0.01,
                   initial_base_layer_count: int = 3,
                   LAMMPS_command = "lmp",
                   **kwargs):
        """Computes the stacking fault properties of an FCC crystal. For more details, refer to README.txt

        Args: 
            pressure (float):
                (optional) Hydrostatic pressure (bars).  If omitted, the pressure is taken to be zero.
                If the value specified is non-zero, the lattice constant specified for
                LatConst will be used to construct an initial lattice geometry for an NPT
                simulation carried out at the specified pressure and temperature of 1e-4
                Kelvin from which the actual lattice constant at the specified pressure is
                calculated.

            Num_layers_gamma_surf (int):
                (optional) Number of layers used for gamma surface creation

            compute_gamma_surf (bool):
                (optional) Determines if gamma surface is computed.
            
            conv_cutoff (float):
                (optional) Default relative error used in convergence study. Calculated using the 
                final two simulation results.

            initial_base_layer_count (int):
                (optional) The initial number of repeating base layers used.

            LAMMPS_command (str):
                (optional) The command for calling LAMMPS.


        Properties calculated: 
                gamma_us,  frac_us         # Unstable stacking fault energy  (max)
                gamma_isf, frac_isf = 1.0  # Intrinsic stacking fault energy (min)
                gamma_ut,  frac_ut         # Unstable twinning fault energy  (max)
                gamma_esf, frac_esf = 2.0  # Entrinsic stacking fault energy (min)
                FracList  = []             # fractional displacements array
                SFEDList  = []             # stacking fault energy densities (eV/A^2)
                GammaSurf = []             # [frac_along_112, frac_along_110, energy (eV/A^2)] in each row

        Supporting modules used:
            make_lammps_input.py  --> Makes LAMMPS input

        External resources used:
            LAMMPS executable with KIM API support
        """

 
        # if statement to check for FCC
        prototype_label = self._get_nominal_crystal_structure_npt()["prototype-label"]["source-value"]
        if prototype_label != 'A_cF4_225_a':
            raise RuntimeError('Only accepts single species FCC')
        

        # get the necessary parameters
        latconst = self._get_nominal_crystal_structure_npt()["a"]["source-value"]
        model = self.kim_model_name
        species = self._get_nominal_crystal_structure_npt()["stoichiometric-species"]["source-value"][0]

        # run simulations
        output_dict = self._main(model, 
                                 species, 
                                 latconst, 
                                 Pressure = pressure, 
                                 Num_layers_gamma_surf = Num_layers_gamma_surf,
                                 compute_gamma_surf = compute_gamma_surf,
                                 conv_cutoff = conv_cutoff,
                                 initial_base_layer_count = initial_base_layer_count,
                                 LAMMPS_command = LAMMPS_command)
        print([f"{i} = {output_dict[i]}" for i in ['gamma_us', 
                                                   'gamma_isf',
                                                   'gamma_ut',
                                                   'gamma_esf',
                                                   'frac_us',
                                                   'frac_ut']])

        ####################################################
        # PROPERTY WRITING
        ####################################################
        
        # gamma-surface
        if compute_gamma_surf == True:
            self._add_property_instance_and_common_crystal_genome_keys(
                property_name="gamma-surface-relaxed-fcc-crystal-npt-crystal-genome",
                write_stress=False,
                write_temp=False,
            )
            self._add_key_to_current_property_instance("cauchy-stress",
                                                    output_dict['CauchyStress'],
                                                    unit="bar")
            self._add_key_to_current_property_instance("fault-plane-shift-fraction-110",
                                                    output_dict['Gamma_Y_dir2_frac'])
            self._add_key_to_current_property_instance("fault-plane-shift-fraction-112",
                                                    output_dict['Gamma_X_dir1_frac'])
            self._add_key_to_current_property_instance("gamma-surface",
                                                    output_dict['GammaSurf'],
                                                    unit="ev/angstrom^2")
            self._add_file_to_current_property_instance("gamma-surface-plot",
                                                    f"{output_dict['GammaSurfPlotName']}.svg")


        # unstable-stacking-energy-fcc-crystal
        self._add_property_instance_and_common_crystal_genome_keys(
            property_name="unstable-stacking-fault-relaxed-energy-fcc-crystal-npt-crystal-genome",
            write_stress=False,
            write_temp=False,
        )
        self._add_key_to_current_property_instance("cauchy-stress",
                                                   output_dict['CauchyStress'],
                                                   unit="bar")
        self._add_key_to_current_property_instance("unstable-stacking-energy",
                                                   output_dict['gamma_us'],
                                                   unit="eV/angstrom^2")
        self._add_key_to_current_property_instance("unstable-slip-fraction",
                                                   output_dict["frac_us"])


        # intrinsic-stacking-fault-energy-fcc-crystal
        self._add_property_instance_and_common_crystal_genome_keys(
            property_name="intrinsic-stacking-fault-relaxed-energy-fcc-crystal-npt-crystal-genome",
            write_stress=False,
            write_temp=False,
        )
        self._add_key_to_current_property_instance("cauchy-stress",
                                                   output_dict['CauchyStress'],
                                                   unit="bar")
        self._add_key_to_current_property_instance("intrinsic-stacking-fault-energy",
                                                   output_dict['gamma_isf'],
                                                   unit="eV/angstrom^2")


        # unstable-twinning-energy-fcc-crystal
        self._add_property_instance_and_common_crystal_genome_keys(
            property_name="unstable-twinning-fault-relaxed-energy-fcc-crystal-npt-crystal-genome",
            write_stress=False,
            write_temp=False,
        )
        self._add_key_to_current_property_instance("cauchy-stress",
                                                   output_dict['CauchyStress'],
                                                   unit="bar")
        self._add_key_to_current_property_instance("unstable-twinning-energy",
                                                   output_dict['gamma_ut'],
                                                   unit="eV/angstrom^2")
        self._add_key_to_current_property_instance("unstable-slip-fraction",
                                                   output_dict['frac_ut'])

        
        # extrinsic-stacking-fault-energy-fcc-crystal
        self._add_property_instance_and_common_crystal_genome_keys(
            property_name="extrinsic-stacking-fault-relaxed-energy-fcc-crystal-npt-crystal-genome",
            write_stress=False,
            write_temp=False,
        )
        self._add_key_to_current_property_instance("cauchy-stress",
                                                   output_dict['CauchyStress'],
                                                   unit="bar")
        self._add_key_to_current_property_instance("extrinsic-stacking-fault-energy",
                                                   output_dict['gamma_esf'],
                                                   unit="eV/angstrom^2")

        # stacking-energy-curve-fcc-crystal
        self._add_property_instance_and_common_crystal_genome_keys(
            property_name="stacking-fault-relaxed-energy-curve-fcc-crystal-npt-crystal-genome",
            write_stress=False,
            write_temp=False,
        )
        self._add_key_to_current_property_instance("cauchy-stress",
                                                   output_dict['CauchyStress'],
                                                   unit="bar")
        self._add_key_to_current_property_instance("fault-plane-shift-fraction",
                                                   output_dict['FracList'])
        self._add_key_to_current_property_instance("fault-plane-energy",
                                                   output_dict['SFEDList'],
                                                   unit="eV/angstrom^2")


    def _main(self, Model, Species, LatConst, Pressure = 0.0, Num_layers_gamma_surf = 14, compute_gamma_surf = True, conv_cutoff = 0.01, initial_base_layer_count = 3, LAMMPS_command = LAMMPS_command):
        # Program Parameter Variables
        total_time_start = time.perf_counter()

        if Pressure == float(0):
                msg = (
                    "\nInfo: Pressure was either specified as zero in input or not provided. "
                    "Forgoing lattice constant calculation and "
                    "proceeding with lattice constant specified.\n"
                )
                print(msg)


        # -------------------------------------------------------------------------------
        #                        Program internal Constants
        # -------------------------------------------------------------------------------
        Gamma_Nx_dir1 = 20
        Gamma_Ny_dir2 = 20

        Num_layers_gamma_surf, N_Twin_Layers_gamma_surf, Rigid_Grp_SIdx_gamma_surf, Rigid_Grp_EIdx_gamma_surf = self._layer_calc(3)

        output_dir = "./output"  # Output directory
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        stack_inp_flnm = output_dir + "/stack.in"  # Input file for LAMMPS
        stack_log_flnm = output_dir + "/stack.log"  # Log file for LAMMPS
        stack_data_flnm = output_dir + "/stack.dat"  # temporary file for lammps output

        # -------------------------------------------------------------------------------
        #      Target variables to be calculated (see definitions in header)
        # -------------------------------------------------------------------------------
        gamma_us = 0.0
        gamma_isf = 0.0
        gamma_ut = 0.0
        gamma_esf = 0.0

        frac_us = 0.0
        frac_ut = 0.0
        Gamma_X_dir1_frac = [0 + x * 1.0 / (Gamma_Nx_dir1 - 1) for x in range(Gamma_Nx_dir1)]
        Gamma_Y_dir2_frac = [0 + y * 1.0 / (Gamma_Ny_dir2 - 1) for y in range(Gamma_Ny_dir2)]
        GammaSurf = []

        if not Pressure:
            # ------------------------------------------------------------------------------
            #                            CASE I - ZERO PRESSURE
            # Either no pressure was specified, or a pressure of zero was specified.
            # Proceed by constructing the FCC lattice using the equilibrium lattice
            # constant given and forming the stacking faults, etc.
            # ------------------------------------------------------------------------------
            Pressure = 0.0

        else:
            # ------------------------------------------------------------------------------
            #                         CASE II - NON-ZERO PRESSURE
            # Use the zero-temperature, zero-pressure equilibrium lattice constant
            # specified in the input to construct the initial FCC lattice for an NPT
            # simulation at 1e-4K and the specified pressure.  After 200,000 timesteps, the
            # length of the supercell along the x direction is parsed from the output and
            # divided by the number of conventional FCC cells along that direction to
            # arrive at the equilibrium lattice constant at the specified pressure.  This
            # lattice constant is then used to construct the lattice geometry for the
            # actual stacking fault calculations.
            # ------------------------------------------------------------------------------
            print(
                "Info: A non-zero pressure of %r bar was specified. Computing the corresponding "
                "FCC lattice constant...\n" % Pressure
            )
            with open(stack_inp_flnm, "w") as fstack:
                InpStr = compute_eq_latconst(
                    Species, Model, LatConst, Pressure, stack_data_flnm
                )
                fstack.write(InpStr)

            # Run the LAMMPS script
            os.system(LAMMPS_command + " -in " + stack_inp_flnm + " -log " + stack_log_flnm)

            # Read the LAMMPS output file for the lattice constant
            with open(stack_data_flnm) as fstack:
                linelist = fstack.readlines()
                linebuf = linelist[0].split()
                LatConst = float(linebuf[0])

            # delete the output file
            os.system("rm " + stack_data_flnm)
            os.system("rm " + stack_inp_flnm)
            print(
                "Info: Calculated an FCC lattice constant of %r corresponding to pressure %r\n"
                % (LatConst, Pressure)
            )

        # ------------------------------------------------------------------------------
        #                            COMPUTE GAMMMA SURFACE
        # ------------------------------------------------------------------------------
        print("***********************************************************")
        print("              COMPUTING GAMMA SURFACE                      ")
        print("***********************************************************")
        time_gamma_start = time.perf_counter()

        if compute_gamma_surf == True:
            with open(stack_inp_flnm, "w") as fstack:
                InpStr = setup_problem(
                    Species,
                    Model,
                    Num_layers_gamma_surf,
                    LatConst,
                    Pressure,
                    Rigid_Grp_SIdx_gamma_surf,
                    Rigid_Grp_EIdx_gamma_surf,
                    N_Twin_Layers_gamma_surf,
                )
                fstack.write(InpStr)
                InpStr = make_gammasurface_moves(stack_data_flnm, Gamma_Nx_dir1, Gamma_Ny_dir2)
                fstack.write(InpStr)

            # Run the LAMMPS script
            os.system(LAMMPS_command + " -in " + stack_inp_flnm + " -log " + stack_log_flnm)

            # Read the LAMMPS output file
            """-----------------------------------------------------------------------------
            File format: stack.dat
            Line 1: Header
            Line 1+1 to 1+NxPoints*NyPoints: [ 112_frac    110_frac    SFED ]
            -----------------------------------------------------------------------------"""
            with open(stack_data_flnm) as fstack:
                linelist = fstack.readlines()
                # Discard the header, index = 0
                # Read the data into arrays
                count = 1
                for yIdx in range(1, Gamma_Ny_dir2 + 1):
                    temp_list_at_each_y = []
                    for xIdx in range(1, Gamma_Nx_dir1 + 1):
                        linebuf = linelist[count].split()
                        count = count + 1
                        # Discard x any y coordinates
                        temp_list_at_each_y.append(float(linebuf[2]))
                    GammaSurf.append(temp_list_at_each_y)

            # delete the output file
            os.system("rm " + stack_data_flnm)
            os.system("rm " + stack_inp_flnm)

        time_gamma_end = time.perf_counter()

        # convergence study for SF energies
        N_layers_list = []
        output_dict_list = []
        us_rough_list = []
        ut_rough_list = []
        sim_time_list = []
        rough_sim_time_list = []
        conv_error = 1.0
        base_layer_count = initial_base_layer_count # sets layer count for initial run
        while conv_error > conv_cutoff:
            base_layer_count += 1
            FracList = []
            SFEDList = []
            (N_Layers, N_Twin_Layers, Rigid_Grp_SIdx, Rigid_Grp_EIdx) = self._layer_calc(base_layer_count)#14)

            # ------------------------------------------------------------------------------
            #                       COMPUTE STACKING FAULT ENERGIES
            # ------------------------------------------------------------------------------
            with open(stack_inp_flnm, "w") as fstack:
                InpStr = setup_problem(
                    Species,
                    Model,
                    N_Layers,
                    LatConst,
                    Pressure,
                    Rigid_Grp_SIdx,
                    Rigid_Grp_EIdx,
                    N_Twin_Layers,
                )
                fstack.write(InpStr)
                InpStr = make_stack_twin_test(stack_data_flnm)
                fstack.write(InpStr)

            time_sf_rough_start = time.perf_counter()
            # Run the LAMMPS script
            os.system(LAMMPS_command + " -in " + stack_inp_flnm + " -log " + stack_log_flnm)

            # Read the LAMMPS output file
            """-----------------------------------------------------------------------------
            File format: stack.dat
            Line 1:             NPoints 1 nincr nincr
            Line 1+1 to 1+1+2*nincr:    [ column1 = frac_disp, column2 = SFED ]
            Line:                gamma_us
            Line:                gamma_isf
            Line:                gamma_ut
            Line:                gamma_esf
            -----------------------------------------------------------------------------"""
            with open(stack_data_flnm) as fstack:
                linelist = fstack.readlines()
                linebuf = linelist[0].split()
                size_0 = int(linebuf[2])
                size_1 = int(linebuf[3])
                size_2 = int(linebuf[4])
                # Read the data into arrays
                listend = 1 + size_0 + size_1 + size_2
                for i in range(1, listend):
                    linebuf = linelist[i].split()
                    FracList.append(float(linebuf[0]))
                    SFEDList.append(float(linebuf[1]))
                # Store the isf and esf values
                gamma_isf = SFEDList[size_0 + size_1 - 1]
                gamma_esf = SFEDList[size_0 + size_1 + size_2 - 1]

            # delete the output file
            os.system("rm " + stack_data_flnm)
            os.system("rm " + stack_inp_flnm)

            time_sf_rough_end = time.perf_counter()

            # Locate the unstable stacking fault energy
            us_rough_Idx, us_rough_val = max(
                enumerate(SFEDList[0 : size_0 + size_1 - 1]), key=operator.itemgetter(1)
            )
            us_rough_list.append(us_rough_val)

            # Locate the unstable stacking fault energy
            ut_rough_Idx, ut_rough_val = max(
                enumerate(SFEDList[size_0 + size_1 : size_0 + size_1 + size_2 - 1]),
                key=operator.itemgetter(1),
            )
            ut_rough_list.append(ut_rough_val)
            rough_sim_time_list.append(time_sf_rough_end - time_sf_rough_start)
            N_layers_list.append(N_Layers)

            if base_layer_count == 14:  
                break
            elif len(us_rough_list) > 1:
                conv_error = self._convergence_check(us_rough_list, ut_rough_list)

        # ------------------------------------------------------------------------------
        #             Refinement to locate the unstable position - gamma_us
        # ------------------------------------------------------------------------------
        # Locate the unstable stacking fault energy
        us_rough_Idx, us_rough_val = max(
            enumerate(SFEDList[0 : size_0 + size_1 - 1]), key=operator.itemgetter(1)
        )

        SFrac_us = FracList[us_rough_Idx - 1]
        dFrac_us = FracList[us_rough_Idx] - FracList[us_rough_Idx - 1]

        # Make input for the refinement
        with open(stack_inp_flnm, "w") as fstack:
            InpStr = setup_problem(
                Species,
                Model,
                N_Layers,
                LatConst,
                Pressure,
                Rigid_Grp_SIdx,
                Rigid_Grp_EIdx,
                N_Twin_Layers,
            )
            fstack.write(InpStr)
            InpStr = make_refine_us(SFrac_us, dFrac_us, stack_data_flnm)
            fstack.write(InpStr)

        time_sf_fine_start = time.perf_counter()
        # Run the LAMMPS script
        os.system(LAMMPS_command + " -in " + stack_inp_flnm + " -log " + stack_log_flnm)

        # Read the Lammps output file
        with open(stack_data_flnm) as fstack:
            linelist = fstack.readlines()
            linebuf = linelist[1].split()
            frac_us = float(linebuf[0])
            gamma_us = float(linebuf[1])

        # delete the output file
        os.system("rm " + stack_data_flnm)
        os.system("rm " + stack_inp_flnm)

        # ------------------------------------------------------------------------------
        #             Refinement to locate the unstable position - gamma_ut
        # ------------------------------------------------------------------------------
        # Locate the unstable stacking fault energy
        ut_rough_Idx, ut_rough_val = max(
            enumerate(SFEDList[size_0 + size_1 : size_0 + size_1 + size_2 - 1]),
            key=operator.itemgetter(1),
        )

        SFrac_ut = FracList[size_0 + size_1 + ut_rough_Idx - 1]
        dFrac_ut = (
            FracList[size_0 + size_1 + ut_rough_Idx]
            - FracList[size_0 + size_1 + ut_rough_Idx - 1]
        )

        # Make input for the refinement
        with open(stack_inp_flnm, "w") as fstack:
            InpStr = setup_problem(
                Species,
                Model,
                N_Layers,
                LatConst,
                Pressure,
                Rigid_Grp_SIdx,
                Rigid_Grp_EIdx,
                N_Twin_Layers,
            )
            fstack.write(InpStr)
            InpStr = make_refine_ut(SFrac_ut, dFrac_ut, stack_data_flnm)
            fstack.write(InpStr)

        # Run the LAMMPS script
        os.system(LAMMPS_command + " -in " + stack_inp_flnm + " -log " + stack_log_flnm)

        # Read the Lammps output file
        with open(stack_data_flnm) as fstack:
            linelist = fstack.readlines()
            linebuf = linelist[1].split()
            frac_ut = float(linebuf[0])
            gamma_ut = float(linebuf[1])

        # delete the output and log files
        os.system("rm " + stack_data_flnm)
        os.system("rm " + stack_inp_flnm)
        os.system("rm " + stack_log_flnm)
        if os.path.exists("kim.log"):
            os.system("rm kim.log")

        time_sf_fine_end = time.perf_counter()

        # # ------------------------------------------------------------------------------
        # #                    PRINT FINAL OUTPUTS TO KIM EDN FORMAT
        # # ------------------------------------------------------------------------------

        # Convert pressure to match relevant KIM Property Definitions
        CauchyStress = [-Pressure, -Pressure, -Pressure, 0.0, 0.0, 0.0]

        # ------------------------------------------------------------------------------
        #         Plot gamma surface to png and svg using matplotlib
        # ------------------------------------------------------------------------------
        
        if compute_gamma_surf == True:
            # Convert data to numpy for matplotlib
            Gamma_X_dir1_frac, Gamma_Y_dir2_frac = (
                np.asarray(Gamma_X_dir1_frac),
                np.asarray(Gamma_Y_dir2_frac),
            )
            GammaSurf = np.array(GammaSurf)
            Gamma_X_dir1_frac_grid, Gamma_Y_dir2_frac_grid = np.meshgrid(
                Gamma_X_dir1_frac, Gamma_Y_dir2_frac
            )

            label112 = r"$\frac{s\,_{[112]}}{\sqrt{6}a/2}$"
            label110 = r"$\frac{s\,_{[\mathrm{\overline{1}}10]}}{\sqrt{2}a/2}$"
            energy_label = r"$\gamma$ (eV/$\mathrm{\AA}^2$)"
            labelfontsize = 15

            # Draw the 2d projection of the gamma surface
            plt.close("all")
            fig = plt.figure()

            ax_2d = fig.add_subplot()
            projected_gamma_surf = ax_2d.pcolor(
                Gamma_X_dir1_frac_grid, Gamma_Y_dir2_frac_grid, GammaSurf, cmap=cm.bone
            )
            ax_2d.set_xlabel(label112, fontsize=labelfontsize)
            ax_2d.set_ylabel(label110, fontsize=labelfontsize)
            fig.colorbar(projected_gamma_surf, shrink=1, aspect=10, label=energy_label)
            fig.subplots_adjust(bottom=0.1)
            GammaSurfPlotName = "gamma-surface-relaxed-fcc-" + Species + "-" + Model + "-projected.png"
            fig.savefig(
                os.path.join(
                    output_dir,
                    f"{GammaSurfPlotName}.png",
                ),
                bbox_inches="tight",
                dpi=300,
            )
            fig.savefig(
                os.path.join(
                    output_dir,
                    f"{GammaSurfPlotName}.svg",
                ),
                bbox_inches="tight",
            )


            output_dict = {'CauchyStress': CauchyStress,
                        'Gamma_X_dir1_frac': Gamma_X_dir1_frac,
                        'Gamma_Y_dir2_frac': Gamma_Y_dir2_frac,
                        'GammaSurf': GammaSurf,
                        'GammaSurfPlotName': f"{output_dir}/{GammaSurfPlotName}",
                        'gamma_us': gamma_us,
                        'gamma_isf': gamma_isf,
                        'gamma_ut': gamma_ut,
                        'gamma_esf': gamma_esf,
                        'frac_us': frac_us,
                        'frac_ut': frac_ut,
                        'FracList': FracList,
                        'SFEDList': SFEDList,
                        }
        
        else:
            output_dict = {'CauchyStress': CauchyStress,
                        'Gamma_X_dir1_frac': Gamma_X_dir1_frac,
                        'Gamma_Y_dir2_frac': Gamma_Y_dir2_frac,
                        'gamma_us': gamma_us,
                        'gamma_isf': gamma_isf,
                        'gamma_ut': gamma_ut,
                        'gamma_esf': gamma_esf,
                        'frac_us': frac_us,
                        'frac_ut': frac_ut,
                        'FracList': FracList,
                        'SFEDList': SFEDList,
                        }

        time_gamma = time_gamma_end - time_gamma_start
        time_sf_rough = time_sf_rough_end - time_sf_rough_start
        time_sf_fine = time_sf_fine_end - time_sf_fine_start
        total_time_end = time.perf_counter()

        print(f"gamma surface time = {time_gamma/60} mins")
        print(f"time sf rough = {time_sf_rough/60} mins")
        print(f"time sf fine = {(time_sf_fine)/60} mins")

        output_dict_list.append(output_dict)
        time_sf = (time_sf_rough + time_sf_fine)/60
        sim_time_list.append(time_sf)
        

        
        # write convergence study results
        gamma_us_convergence = [i['gamma_us'] for i in output_dict_list]
        gamma_ut_convergence = [i['gamma_ut'] for i in output_dict_list]

        cresults_text = f"N_layers = {N_layers_list}\n"
        cresults_text += f"gamma_us = {gamma_us_convergence}\n"
        cresults_text += f"gamma_ut = {gamma_ut_convergence}\n"
        cresults_text += f"rough_us = {us_rough_list}\n"
        cresults_text += f"rough_ut = {ut_rough_list}\n"
        cresults_text += f"indv_rough_sim_time_seconds = {rough_sim_time_list}\n"
        cresults_text += f"total_sim_time_seconds = {total_time_end - total_time_start}\n"

        with open("./output/conv_results.md", "w") as cresults:
            cresults.write(cresults_text)

        return output_dict


    def _layer_calc(self, i_in):
        """calculate layer parameters for a given base layer count
        """
        N_layers = 10 + (i_in - 2)*4
        N_twin_layers = N_layers/2
        Rigid_Grp_SIdx = i_in + 1
        Rigid_Grp_EIdx = N_twin_layers + Rigid_Grp_SIdx - 1
        return N_layers, N_twin_layers, Rigid_Grp_SIdx, Rigid_Grp_EIdx
    
    
    def _convergence_check(self, us_rough_list, ut_rough_list):
        """relative error check for convergence study
        """
        error_list = []
        error_list.append((us_rough_list[-1]-us_rough_list[-2])/us_rough_list[-2])
        error_list.append((ut_rough_list[-1]-ut_rough_list[-2])/ut_rough_list[-2])
        return max([max(error_list),abs(min(error_list))])