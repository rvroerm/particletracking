# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import sqrt, atan, cos, sin
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from transfer_functions import EOM_particle_in_Bfield, PtoV, EtoP, Brho_scale_transport

import pandas as pd


import folium

plt.close('all')

input_file = "C:/TRANS/for001.dat"

test = Brho_scale_transport(input_file,230)

route_map = folium.Map(location=[51.1657, 10.4515], zoom_start=6, tiles='Stamen Toner')

route_map.save("mymap.html")
