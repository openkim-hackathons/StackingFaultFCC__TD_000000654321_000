#!/usr/bin/python

"""
Invoking a Crystal Genome Test Driver Directly
==============================================
"""
from test_driver.test_driver import TestDriver
import time
from ase.build import bulk

time_begin = time.perf_counter()

# default FCC test
kim_model_name = 'EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005'
atoms = bulk('Al','fcc',a=4.032,cubic=True)
test_driver = TestDriver(kim_model_name)
test_driver(atoms)

time_end = time.perf_counter()  
print(f"total time = {(time_end - time_begin)/60} mins")

# make sure it errors out for non-FCC
kim_model_name = 'EAM_Dynamo_AcklandBaconCalder_1997_Fe__MO_142799717516_005'
atoms = bulk('Fe','bcc',a=2.866,cubic=True)
test_driver = TestDriver(kim_model_name)
test_driver(atoms)