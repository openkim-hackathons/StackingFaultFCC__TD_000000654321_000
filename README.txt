--------------------------------------------------------------------------------
                  Intrinsic-Extrinsic stacking fault test
--------------------------------------------------------------------------------

I. Applicability: FCC materials

II. Physical properties computed by the test:

   1. Unstable  stacking fault energy (eV/A^2) and corresponding slip fraction normalized by a/sqrt(6)
   2. Intrinsic stacking fault energy (eV/A^2)
   3. Unstable  twinning fault energy (eV/A^2) and corresponding slip fraction normalized by a/sqrt(6)
   4. Extrinsic stacking fault energy (eV/A^2)
   5. Curve of fractional displacements versus stacking fault energy (eV/A^2)
   6. Gamma surface

III. Theory: Stacking faults are planar defects arising due to irregularities in the planar stacking sequence of atoms
     in a crystal. Such departures from the ideal stacking order increase the energy of the material. A knowledge of
     stacking faults is crucial to the understanding of dislocations and their mobility. Dislocations in FCC crystals
     can dissociate into partial dislocations separated by stacking fault ribbons. A key parameter is the stacking
     fault energy (energy per unit area of the stacking fault) that controls the splitting distance between partials.

     In a defect-free FCC crystal, close-packed layers of {111} planes form the sequence ...ABCABCABC... .

     1. Intrinsic stacking fault: The stacking sequence in the presence of an intrinsic stacking fault is
        ...ABC|BCABC..., as if a single plane "A" is removed from the sequence. The increase in energy per unit area of
        the stacking fault with respect to a defect free crystal is called the intrinsic stacking fault energy.

     2. Unstable stacking fault energy: An intrinsic stacking fault can be created by rigidly slipping one half of an
        infinite crystal on a {111} plane along a partial <112> direction. The energy barrier encountered
        during this process is termed the unstable stacking fault energy.

        Note that the amount of the slip corresponding to the intrinsic stacking fault along a <112> direction is
        b = a/sqrt(6) (where a is the lattice parameter). The unstable stacking energy corresponds to a slip of
        between 0.5*b and b. (This value is not 0.5*b due to the influence of neighbouring layers.)

     3. Extrinsic stacking fault: In the extrinsic stacking fault, one extra plane is inserted into the perfect stacking
        sequence making it ...ABC|B|ABC.... This is created by a rigid slip in the presence of an existing stacking
        fault. That is, we take a configuration with a stacking fault ...ABC|BCABC... and rigidly slip the plane indicated
        by the asterisk ...ABC|B*CABC... thereby making it ...ABC|B|ABC.... The slip directions and magnitudes are the
        same as that of the stacking fault creation process. The associated energy cost with respect to the perfect
        lattice is called the extrinsic stacking fault energy.

     4. Unstable twinning energy: Similar to the unstable stacking fault energy, the unstable twin energy is also a barrier
        between the intrinsic and extrinsic stacking fault configurations during the process of rigid slip explained
        in point 3. This corresponds to the barrier for nucleating a two-layer micro twin in the crystal.

     6. Gamma surface: The gamma surface is created by rigid slip of a (111) plane on a grid of points defined by [112]
        and [-110] directions in an fcc crystal at a specified pressure and temperature.  It is the energy-per-area versus
        all possible slips lying in the (111) lattice plane. Due to periodicity of the crystal lattice, it suffices to
        sample a grid of points that span  a*sqrt(6)/2 and a*sqrt(2)/2  along the [112] and [-110] directions, respectively.
        This is achieved through a sequence of rigid displacements applied to one part of an fcc crystal relative to another
        on the (111) plane on a grid defined by [112] and [-110] directions at the specified pressure and temperature.

IV: User Inputs:

    1. Extended ID of a KIM Model
    2. An atomic species from which an FCC lattice is constructed
    3. The zero-temperature, zero-pressure equilibrium lattice constant (meters)
    4. Optionally, a hydrostatic pressure (bars) may be specified.  If omitted, the pressure is taken to be zero.  If a
       non-zero value is given, the equilibrium lattice constant specified will be used to construct the initial lattice
       geometry for an NPT simulation carried out at the specified pressure and temperature of 1e-4 Kelvin, from which the
       actual lattice constant at the specified pressure is calculated.

V. Procedure (LAMMPS): For the code, see the file "make_lammps_input.py" in the test driver.

   NOTE: We apply hydrostatic pressure to the system and the NPT or box/relax fixes in LAMMPS which require
         periodic boundary conditions. To create a periodic sample, we twin the crystal about the center {111} plane of
         the periodic box, and then take a block of layers centered across this twinning plane. Now if this block of
         layers is displaced in the <112> direction both the top and the bottom surfaces see the same relative motion and
         thus the same energy contribution. Hence the stacking energy per plane is exactly half of the total energy
         of the box with a stacking fault less the energy of the initial twinned crystal.

     Common Steps/Subroutines in LAMMPS:

     P1. Create a periodic box in a slab geometry with FCC lattice stacking with 58 (11-1) layers.

     P2. Twin the box about the center, i.e. the top 29 layers are twinned with respect to the bottom. The large number
         of layers is required to ensure that the existing twin plane faults and stacking fault at the edges of the
         periodic cell are sufficiently far from new defects being studied.

     P3. Define groups of atoms as follows:
          (1). stack_group containing the layers 15 to 45 (closed interval), which is centered around the middle twin
               plane. When this block is displaced the stacking faults created are sufficiently far from the twin
               plane.
          (2). twin_group containing the layers 16 to 44 (closed interval)


     (a). TEST: Stacking-Twinning Curve:
        (1). Do P1 through P3 to setup the problem.
        (2). Displace the stack_group in the [112] direction through a/sqrt(6) in a predetermined number of increments.
        (3). Displace the twin_group in the [112] direction through a/sqrt(6) in a predetermined number of increments.
             Store the fractional displacement and the fault energy density at each incremental displacement during
             steps (2) and (3)

     (b). TEST: Refinement of the unstable stacking fault energy.
        (1). Do P1 through P3 to setup the problem.
        (2). From the previously explored curve in test (a), pick the partial displacement fraction just before the
             first local maximum (unstable stacking fault energy) and displace the stack_group right up to this point.
             From this point, take incremental steps and use the bisection method moving the stack_group to
             converge on the maximum (barrier) to a specified tolerance.

     (c). TEST: Refinement of the unstable twinning fault energy.
        (1). Do P1 through P3 to setup the problem.
        (2). Displace the stack_group in the [112] direction through a/sqrt(6) thereby creating a stacking fault.
        (3). From the previously explored curve in test (a), pick the partial displacement fraction just before the
             second local maximum (call it frac_x) (unstable stacking fault energy) and displace the twin_group right
             up to this point, i.e. displace by a fraction of (frac_x-1) since the stacking fault in step (2)
             contributes to frac_b by a quantity 1.0.  From this point, take incremental steps and and use the
             the bisection method moving the twin_group to converge on the maximum (barrier) to a specified tolerance.

     (d). TEST: Gamma surface:
        (1). Do P1 through P3 to setup the problem.
        (2). Displace the stack_group in the rectangular grid defined by motions along [112] and [-110] directions
             through a*sqrt(6)/2 and a*sqrt(2)/2 respectively at some predetermined number of increments in each
             direction thereby sampling the entire grid.

VI. Short description of various files in the test driver:

   1. test_driver: The main script
   2. make_lammps_input.py: Used by "runner" to generate the LAMMPS input files
   4. kimspec.edn: Describes the metadata associated with the Test Driver
   5. LICENSE.CDL: License document
   6. Makefile: Dummy file

VII. Reference:
     1. Tight-binding calculations of stacking energies and twinnability in fcc
        metals, N. Bernstein and E. B. Tadmor, Phys. Rev. B 69, 094116, 2004.
     2. S. Pattamatta, Stacking and twinning fault energies of an fcc lattice at 
        zero temperature and pressure v002. (2019). OpenKIM. doi: 10.25950/B4CFAF9A.
