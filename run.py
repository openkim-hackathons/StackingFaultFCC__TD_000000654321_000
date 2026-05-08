#!/usr/bin/python

"""
Invoking a Crystal Genome Test Driver Directly
==============================================
"""
from test_driver.test_driver import TestDriver
import time
from ase.build import bulk

time_begin = time.perf_counter()

# temporary, for convg study check
compute_gamma_surf = False

if True:
    # atoms object testing
    # default FCC test
    kim_model_name = 'EAM_Dynamo_WangZhuXiang_2018pot2_Pb__MO_961101070310_001'
    atoms = bulk('Pb', 'fcc', a=4.989170034100908)
    # kim_model_name = 'EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_006'
    # atoms = bulk('Al','fcc',a=4.032081970847309)

    test_driver = TestDriver(kim_model_name)
    # p = 0
    p = 6.3242091e-07 # 1 atm
    # p = 0.062
    # p = 0.131 # should result in 10% compression

    test_driver(atoms, pressure_eV_angstrom3 = p, compute_gamma_surf = compute_gamma_surf)
    test_driver.write_property_instances_to_file()

    time_end = time.perf_counter()  
    print(f"total time = {(time_end - time_begin)/60} mins")

    # make sure it errors out for non-FCC
    kim_model_name = 'EAM_Dynamo_AcklandBaconCalder_1997_Fe__MO_142799717516_005'
    atoms = bulk('Fe','bcc',a=2.866,cubic=True)
    test_driver = TestDriver(kim_model_name)
    test_driver(atoms)

if False: # kimvv testing
    from kimvv import EquilibriumCrystalStructure
    atoms_init = bulk('Au')

    kim_model_names = [#"LennardJones612_UniversalShifted__MO_959249795837_003",
                       "Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu",
                       ]

    for kim_model_name in kim_model_names:
        # Instantiate the Equilibrium Driver with your model
        ecs = EquilibriumCrystalStructure(kim_model_name)

        # Relax the structure. ECS will return multiple properties, any of them will do as they all contain the
        # crystal description
        relaxed_structure = ecs(atoms_init)[0]

        # Run your TD with `relaxed_structure` as the input
        test_driver = TestDriver(kim_model_name)
        test_driver(relaxed_structure, pressure_eV_angstrom3 = 0.0006, compute_gamma_surf = compute_gamma_surf)
        test_driver.write_property_instances_to_file()

        time_end = time.perf_counter()  
        print(f"total time = {(time_end - time_begin)/60} mins")

