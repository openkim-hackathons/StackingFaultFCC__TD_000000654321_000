import textwrap


# -------------------------------------------------------------------------------
# make_eq_latconst: Given the pressure and initial lattice constant, makes input
#                   for the lattice constant at the given pressure
# -------------------------------------------------------------------------------
def compute_eq_latconst(Species, ModelName, init_LatConst, Pressure, stack_data_flnm):

    lammps_input = """
   kim init {ModelName} metal unit_conversion_mode

   dimension 3

   boundary p p p

   # Set neighbor skin
   variable neigh_skin equal 2.0*${{_u_distance}}
   neighbor ${{neigh_skin}} bin

   # Atom definition
   lattice fcc {init_LatConst}*${{_u_distance}}
   region whole block 0 1 0 1 0 1 units box
   create_box 1 whole
   create_atoms 1 region whole
   mass * 1.0

   # Define the interatomic potentials
   kim interactions {Species}

   # Variables used to rescale the energy and stress so that the quantities
   # in the thermo output are in the original metal units (eV and bars)
   # even if we're running with a Simulator Model that uses different units
   variable pe_metal    equal "c_thermo_pe/v__u_energy"
   variable lx_metal    equal lx/${{_u_distance}}
   variable ly_metal    equal ly/${{_u_distance}}
   variable lz_metal    equal lz/${{_u_distance}}
   variable press_metal equal "c_thermo_press/v__u_pressure"
   variable pxx_metal   equal pxx/${{_u_pressure}}
   variable pyy_metal   equal pyy/${{_u_pressure}}
   variable pzz_metal   equal pzz/${{_u_pressure}}
   variable temp_metal  equal temp/${{_u_temperature}}

   # Settings
   thermo 10000
   thermo_style custom step v_lx_metal v_ly_metal v_lz_metal &
                       v_press_metal v_pxx_metal v_pyy_metal &
                       v_pzz_metal v_pe_metal temp v_temp_metal

   # Set up NPT ensemble
   variable Pressure_converted equal {Pressure}*${{_u_pressure}}
   variable Pdamp_converted  equal 100*${{_u_time}}
   variable Tstart_converted equal 0.0001*${{_u_temperature}}
   variable Tstop_converted  equal 0.0001*${{_u_temperature}}
   variable Tdamp_converted  equal 0.005*${{_u_temperature}}

   reset_timestep 0
   timestep 0.001*${{_u_time}}
   fix 1 all npt temp ${{Tstart_converted}} ${{Tstop_converted}} ${{Tdamp_converted}} &
       iso ${{Pressure_converted}} ${{Pressure_converted}} ${{Pdamp_converted}} drag 10.0
   run 200000

   # Print the equilibrated lattice constant, converted back to our metal units
   # (Angstroms) if necessary
   variable a equal lx/${{_u_distance}}
   print "${{a}}" file {stack_data_flnm} screen yes

   # Simulation Done
   print "Computed equilibrium Lattice Constant." """.format(
        Species=Species,
        ModelName=ModelName,
        Pressure=Pressure,
        stack_data_flnm=stack_data_flnm,
    )
    return textwrap.dedent(lammps_input)


# -------------------------------------------------------------------------------
# setup_problem: Sets up the problem and defines geometry etc.,
# -------------------------------------------------------------------------------
def setup_problem(
    Species,
    ModelName,
    N_Layers,
    LatConst,
    Pressure,
    Rigid_Grp_SIdx,
    Rigid_Grp_EIdx,
    N_Twin_Layers,
):

    lammps_input = """
   kim init {ModelName} metal unit_conversion_mode

   dimension 3

   boundary p p p

   # Set neighbor skin
   variable neigh_skin equal 2.0*${{_u_distance}}
   neighbor ${{neigh_skin}} bin

   # Atom definition
   variable latparam_converted equal {LatConst}*${{_u_distance}}

   # Declare the dimensions of the model in lattice cells
   variable nlat_x equal 1
   variable nlat_y equal 1
   variable nlat_z equal {N_Layers}

   # Define the box/various dimensions. Start with x and y dimensions.
   # Before converting to the active unit set, get their value in Angstroms
   # for later normalization
   variable xdim     equal ({LatConst}*sqrt(6)/2)*${{nlat_x}}
   variable ydim     equal ({LatConst}*sqrt(2)/2)*${{nlat_y}}

   # Area of the stacking fault plane in Angstrom^2
   variable Area equal 2*${{xdim}}*${{ydim}}

   # Now convert xdim and ydim to the active unit set and do z direction
   variable xdim     equal ${{xdim}}*${{_u_distance}}
   variable ydim     equal ${{ydim}}*${{_u_distance}}
   variable zdim     equal (${{latparam_converted}}/sqrt(3))*${{nlat_z}}
   variable layer_z  equal (${{latparam_converted}}/sqrt(3))

   # Geometry definition
   lattice fcc ${{latparam_converted}}

   # Define regions for stacking fault test
   # NOTE: -0.001 in the x- and y- directions only to overcome numerical
   #       precision issues (no need to worry about units on these quantities)

   # Define a region for each  of the top N_layers/2
   variable ntwin_layers equal ${{nlat_z}}/2

   variable i loop ${{ntwin_layers}}
   label start_of_loop_region

      variable j equal ${{i}}+${{ntwin_layers}}
      variable zmin equal (${{j}}-1.0)*${{layer_z}}-0.001
      variable zmax equal ${{j}}*${{layer_z}}-0.001
      region  ${{i}}  block  -0.001 ${{xdim}} -0.001 ${{ydim}}  ${{zmin}} ${{zmax}} units box
      next i

   jump SELF start_of_loop_region

   variable i delete
   variable j delete

   # Declare group for the rigid blocks for intrinsic stacking fault
   variable zmin equal ({Rigid_Grp_SIdx}-1.0)*${{layer_z}}-0.001
   variable zmax equal ({Rigid_Grp_EIdx}-1.0)*${{layer_z}}+0.001
   region  stack_region block  -0.001 ${{xdim}} -0.001 ${{ydim}}  ${{zmin}} ${{zmax}} units box

   # Declare group for the rigid blocks for the extrisic stacking fault
   variable zmin equal ({Rigid_Grp_SIdx_plus1}-1.0)*${{layer_z}}-0.001
   variable zmax equal ({Rigid_Grp_EIdx_minus1}-1.0)*${{layer_z}}+0.001
   region   twin_region block  -0.001 ${{xdim}} -0.001 ${{ydim}}  ${{zmin}} ${{zmax}} units box

   # Create simulation box and atoms
   region       whole block  0  ${{xdim}} &
                             0  ${{ydim}} &
                             0  ${{zdim}}  units box
   create_box   1     whole
   lattice fcc ${{latparam_converted}} orient x 1 1 2 orient y -1 1 0 orient z -1 -1 1
   create_atoms 1 region whole

   # Define the interatomic potentials
   kim interactions {Species}
   mass * 1.0

   # Rigid holding settings for each layer as a group
   variable i loop ${{ntwin_layers}}
   label start_of_loop_group

      group ${{i}} region ${{i}}
      next i

   jump SELF start_of_loop_group
   variable i delete

   # Rigid holding settings for the stacking fault group
   group stack_group region stack_region

   # Rigid holding settings for the twinning fault group
   group twin_group region twin_region

   # Variables used to rescale the energy and stress so that the quantities
   # in the thermo output are in the original metal units (eV and bars)
   # even if we're running with a Simulator Model that uses different units
   variable pe_metal    equal "c_thermo_pe/v__u_energy"
   variable press_metal equal "c_thermo_press/v__u_pressure"
   variable pxx_metal   equal pxx/${{_u_pressure}}
   variable pyy_metal   equal pyy/${{_u_pressure}}
   variable pzz_metal   equal pzz/${{_u_pressure}}
   variable temp_metal  equal temp/${{_u_temperature}}

   # Compute initial energy and output
   thermo 100
   thermo_style custom step v_pe_metal press &
           v_press_metal v_pxx_metal v_pyy_metal v_pzz_metal &
           temp v_temp_metal

   #dump config all atom 1000 dump.Stack

   # Perform initial relaxation
   reset_timestep 0
   fix 1 all box/relax x {Pressure} y {Pressure} z {Pressure} vmax 0.01
   min_style cg
   minimize 1e-25 1e-25 10000 10000

   # Twin top N_Twin_Layers
   variable twin_move equal -(1.0*${{latparam_converted}}/sqrt(6))
   variable i loop {N_Twin_Layers}
   label start_of_loop_twin

      variable k equal ${{i}}*${{twin_move}}
      displace_atoms ${{i}} move ${{k}} 0.0 0.0 units box
      next i

   jump SELF start_of_loop_twin
   variable i delete
   variable k delete
   minimize 1e-25 1e-25 100000 100000
   variable E equal "v_pe_metal"
   variable Eini equal ${{E}} """.format(
        Species=Species,
        ModelName=ModelName,
        N_Layers=int(N_Layers),
        LatConst=LatConst,
        Pressure=Pressure,
        Rigid_Grp_SIdx=int(Rigid_Grp_SIdx),
        Rigid_Grp_EIdx=int(Rigid_Grp_EIdx),
        Rigid_Grp_SIdx_plus1=Rigid_Grp_SIdx + 1,
        Rigid_Grp_EIdx_minus1=Rigid_Grp_EIdx - 1,
        N_Twin_Layers=int(N_Twin_Layers),
    )

    return textwrap.dedent(lammps_input)


# -------------------------------------------------------------------------------
# make_stack_twin_test:  Makes LAMMPS input for moving both the stacking planes
# completely one after the other
# -------------------------------------------------------------------------------
def make_stack_twin_test(stack_data_flnm):

    lammps_input = """
   variable outfile string "{stack_data_flnm}"
   variable n_incr equal 100
   variable inc_x equal ${{twin_move}}/${{n_incr}}

   print "Npoints =  1  ${{n_incr}}   ${{n_incr}} " file ${{outfile}} screen no

   # Inital point is zero energy and zero partial fraction
   print "0 0" append ${{outfile}} screen no

   # Apply the incremental displacement until first stacking fault nucleates
   variable i loop ${{n_incr}}
   label start_of_loop_1

      displace_atoms stack_group move ${{inc_x}} 0.0 0.0 units box
      fix 2 all setforce 0 0 NULL
      velocity all zero linear
      minimize 1e-25 1e-25 10000 10000

      variable Ecur equal ${{pe_metal}}
      variable SFED equal (${{Ecur}}-${{Eini}})/${{Area}}
      variable totdisp equal ${{i}}/${{n_incr}}
      print "${{totdisp}} ${{SFED}}" append ${{outfile}} screen no

      next i

   jump SELF start_of_loop_1
   variable i delete

   # Apply the incremental displacement until twinning fault nucleates
   variable i loop ${{n_incr}}
   label start_of_loop_2

      displace_atoms twin_group move ${{inc_x}} 0.0 0.0 units box
      fix 2 all setforce 0 0 NULL
      velocity all zero linear
      minimize 1e-25 1e-25 10000 10000

      variable Ecur equal ${{pe_metal}}
      variable SFED equal (${{Ecur}}-${{Eini}})/${{Area}}
      variable totdisp equal (1.0+${{i}}/${{n_incr}})
      print "${{totdisp}} ${{SFED}}" append ${{outfile}} screen no

      next i

   jump SELF start_of_loop_2
   variable i delete

   # SIMULATION DONE
   print "All done" """.format(
        stack_data_flnm=stack_data_flnm
    )

    return textwrap.dedent(lammps_input)


# -------------------------------------------------------------------------------
# make_refine_us: Makes LAMMPS input for refining gamma_us
# -------------------------------------------------------------------------------
def make_refine_us(SFrac_us, dFrac_us, stack_data_flnm):

    lammps_input = """
   # Create the output file and write header
   variable outfile string "{stack_data_flnm}"

   print "Refined Values" file ${{outfile}} screen no

   # Apply the displacement until lower bound of unstable stacking stacking fault
   variable disp equal ${{twin_move}}*{SFrac_us}
   displace_atoms stack_group move ${{disp}} 0.0 0.0 units box

   # Relax/minimize in z direction
   fix 2 all setforce 0 0 NULL
   minimize 1e-25 1e-25 10000 10000

   variable Ecur equal ${{pe_metal}}
   variable SFEDcur equal (${{Ecur}}-${{Eini}})/${{Area}}
   variable SFEDprev equal ${{SFEDcur}}
   variable Refdisp_us equal ${{disp}}

   # Refinement of the bottom stack plane motion
   variable d_disp equal 0.5*${{twin_move}}*{dFrac_us}
   variable i loop 10000
   label start_of_loop_refine_stack

      displace_atoms stack_group move ${{d_disp}} 0.0 0.0 units box
      variable tmp equal ${{Refdisp_us}}
      variable Refdisp_us equal ${{tmp}}+${{d_disp}}
      run 0

      # Relax in z direction
      min_style cg
      minimize 1e-25 1e-25 10000 10000

      variable Ecur equal ${{pe_metal}}
      variable SFEDcur equal (${{Ecur}}-${{Eini}})/${{Area}}
      variable tmp_d_disp equal ${{d_disp}}
      if "(${{SFEDcur}} < ${{SFEDprev}})" then "variable d_disp equal -0.5*${{tmp_d_disp}}"
      variable diffpe equal abs(${{SFEDcur}}-${{SFEDprev}})
      if "( ${{diffpe}} < 1e-10 )" then "jump SELF break_us"

      variable SFEDprev equal ${{SFEDcur}}
      next i

   jump SELF start_of_loop_refine_stack
   label break_us
   variable i delete
   variable Frac_us equal ${{Refdisp_us}}/${{twin_move}}
   print "${{Frac_us}}  ${{SFEDcur}}" append "{stack_data_flnm}" screen no
   """.format(
        stack_data_flnm=stack_data_flnm, SFrac_us=SFrac_us, dFrac_us=dFrac_us
    )

    return textwrap.dedent(lammps_input)


# -------------------------------------------------------------------------------
# make_refine_ut: Makes LAMMPS input for refining gamma_ut
# -------------------------------------------------------------------------------
def make_refine_ut(SFrac_ut, dFrac_ut, stack_data_flnm):

    lammps_input = """
   # Create the output file and write header
   variable outfile string "{stack_data_flnm}"

   print "Refined Values" file ${{outfile}} screen no

   # Apply the complete displacement stack block
   variable disp equal ${{twin_move}}
   displace_atoms stack_group move ${{disp}} 0.0 0.0 units box

   # Relax in z direction
   fix 2 all setforce 0 0 NULL
   minimize 1e-25 1e-25 10000 10000

   # Apply the displacement to the twin block untill the lower bound of the unstable twin
   variable disp equal ${{twin_move}}*({SFrac_ut}-1.0)
   displace_atoms twin_group move ${{disp}} 0.0 0.0 units box

   # Relax in z direction
   minimize 1e-25 1e-25 10000 10000

   variable Ecur equal ${{pe_metal}}
   variable SFEDcur equal (${{Ecur}}-${{Eini}})/${{Area}}
   variable SFEDprev equal ${{SFEDcur}}
   variable Refdisp_ut equal ${{twin_move}}+${{disp}}

   # Refinement of the twin block motion
   variable d_disp equal 0.5*${{twin_move}}*{dFrac_ut}
   variable i loop 10000
   label start_of_loop_refine_twin

      displace_atoms twin_group move ${{d_disp}} 0.0 0.0 units box
      variable tmp equal ${{Refdisp_ut}}
      variable Refdisp_ut equal ${{tmp}}+${{d_disp}}
      run 0

      # Relax in z direction
      minimize 1e-25 1e-25 10000 10000

      variable Ecur equal ${{pe_metal}}
      variable SFEDcur equal (${{Ecur}}-${{Eini}})/${{Area}}
      variable tmp_d_disp equal ${{d_disp}}
      if "(${{SFEDcur}} < ${{SFEDprev}})" then "variable d_disp equal -0.5*${{tmp_d_disp}}"
      variable diffpe equal abs(${{SFEDcur}}-${{SFEDprev}})
      if "( ${{diffpe}} < 1e-10 )" then "jump SELF break_ut"

      variable SFEDprev equal ${{SFEDcur}}
      next i

   jump SELF start_of_loop_refine_twin
   label break_ut
   variable i delete
   variable Frac_ut equal ${{Refdisp_ut}}/${{twin_move}}
   print "${{Frac_ut}}  ${{SFEDcur}}" append "{stack_data_flnm}" screen no
   """.format(
        stack_data_flnm=stack_data_flnm, SFrac_ut=SFrac_ut, dFrac_ut=dFrac_ut
    )
    return textwrap.dedent(lammps_input)


# -------------------------------------------------------------------------------
# make_gammasurface_moves: Makes LAMMPS input for moving the stacking block
# and sweep the entire gamma surface
# -------------------------------------------------------------------------------
def make_gammasurface_moves(stack_data_flnm, NxPoints, NyPoints):

    lammps_input = """
   # Create the output file and write header
   variable outfile string "{stack_data_flnm}"
   variable nx_points equal {NxPoints}
   variable ny_points equal {NyPoints}
   variable n_incr_x equal ${{nx_points}}-1
   variable n_incr_y equal ${{ny_points}}-1

   variable inc_x equal (-1.0*${{latparam_converted}}*sqrt(6)/2)/${{n_incr_x}}
   variable inc_y equal (1.0*${{latparam_converted}}*sqrt(2)/2)/${{n_incr_y}}

   print "Header: Printing Gamma Surface" file ${{outfile}} screen no

   #dump config all atom 1000 dump.Stack
   fix 2 all setforce 0 0 NULL

   # Outer loop: Apply the y - displacement
   variable j loop ${{ny_points}}
   label start_of_loop_y

      # Apply y-displacement from 2nd step
      if "$j > 1" then &
         "displace_atoms stack_group move 0.0 ${{inc_y}} 0.0 units box"

      # Inner loop: Apply the x- displacement
      variable i loop ${{nx_points}}
      label start_of_loop_x

         # Apply x-displacement from 2nd step
         if "$i > 1" then &
            "displace_atoms stack_group move ${{inc_x}} 0.0 0.0 units box" &

         # Relax in z direction
         minimize 1e-25 1e-25 10000 10000

         variable Ecur equal ${{pe_metal}}
         variable SFED equal (${{Ecur}}-${{Eini}})/${{Area}}
         variable totdisp_x equal (${{i}}-1)/${{n_incr_x}}
         variable totdisp_y equal (${{j}}-1)/${{n_incr_y}}
         print "${{totdisp_x}} ${{totdisp_y}} ${{SFED}}" append ${{outfile}} screen no

         next i

      jump SELF start_of_loop_x
      next j

   jump SELF start_of_loop_y
   variable i delete
   variable j delete

   # SIMULATION DONE
   print "All done" """.format(
        stack_data_flnm=stack_data_flnm, NxPoints=int(NxPoints), NyPoints=int(NyPoints)
    )
    return textwrap.dedent(lammps_input)
