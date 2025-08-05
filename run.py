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
compute_gamma_surf = True

if False: # legacy run
    # default FCC test
    kim_model_name = 'EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005'
    atoms = bulk('Al','fcc',a=4.032,cubic=True)
    test_driver = TestDriver(kim_model_name)
    test_driver(atoms, compute_gamma_surf = compute_gamma_surf)

    time_end = time.perf_counter()  
    print(f"total time = {(time_end - time_begin)/60} mins")

    # make sure it errors out for non-FCC
    kim_model_name = 'EAM_Dynamo_AcklandBaconCalder_1997_Fe__MO_142799717516_005'
    atoms = bulk('Fe','bcc',a=2.866,cubic=True)
    test_driver = TestDriver(kim_model_name)
    test_driver(atoms)

if True: # kimvv testing
    from kimvv import EquilibriumCrystalStructure
    atoms_init = bulk('Au')

    # Instantiate the Equilibrium Driver with your model
    kim_model_name = "Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu"
    ecs = EquilibriumCrystalStructure(kim_model_name)

    # Relax the structure. ECS will return multiple properties, any of them will do as they all contain the
    # crystal description
    relaxed_structure = ecs(atoms_init)[0]

    # Run your TD with `relaxed_structure` as the input
    test_driver = TestDriver(kim_model_name)
    test_driver(relaxed_structure, compute_gamma_surf = compute_gamma_surf)

    time_end = time.perf_counter()  
    print(f"total time = {(time_end - time_begin)/60} mins")

