# CO2-Rock-Physics-Code-with-Compliant-Porosity-and-CO2-Weakening
O4934
An NRAP funded code development for CO2 storage. This code converts reservoir simulations to seismic properties for underground CO2 storage. 

Instructions:
python rock_physics_co2.py testdir mod1_090y.csv baseline_y0.csv

Csv files need to have the following columns:
Basline file - 
element, x, y, z, material ID, permfactor, porosity, horizontal permeability, vertical permeability, Pore pressure (Pascal), T (C), and salt_fractions

Time step files (e.g., mod1_090y.csv):
element, x, y, z, pore pressure, temp, gas saturation, salt

These csv files come from TOUGH (reservoir simulator) outputs. 
