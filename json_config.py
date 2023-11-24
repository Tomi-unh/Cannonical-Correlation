# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 08:06:34 2023

@author: may7e
"""
import json

def save_CCA_variables():
    # Define global variables for CCA
    CCA_variables = {
        "datetime": "20120101-0915",
        "duration": 28,
        "steps": 5,
        "window_size": 128,
        "path": "../../../data/supermag/",
        "save_path": "../TWINS/CCA/"
    }

    # Specify the filename for the JSON file
    json_CCA = "CCA_variables.json"

    # Write the global variables to a JSON file
    with open(json_CCA, 'w') as json_file:
        json.dump(CCA_variables, json_file)

    print(f"Global variables for CCA have been saved to {json_CCA}.")

def save_network_variables():
    # Define global variables for network
    network_variables = {
        "date_str": "19971216-2130",
        "num_stations": 494,
        "num_processes": 10,
        "steps": 5,
        "SM_path": "../../../data/supermag"
    }

    # Specify the filename for the JSON file
    json_network = "Network_variables.json"

    # Write the global variables to a JSON file
    with open(json_network, 'w') as json_file:
        json.dump(network_variables, json_file)

    print(f"Global variables for network have been saved to {json_network}.")

# Uncomment and call the specific function you want to run
# save_CCA_variables()
# save_network_variables() 
