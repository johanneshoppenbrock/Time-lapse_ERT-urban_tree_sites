import pygimli as pg
import os 

from Hilfsfunktionen import processing 
from Hilfsfunktionen import adjust
from Hilfsfunktionen import read 
from Hilfsfunktionen import filter

# read data
dataset = 'two_timesteps'
wenner, dd = read(filename = dataset)

# Merge Wenner and DD
for wen in wenner:
    for dipol in dd:
        if wen[0] == dipol[0]:
            wen[1].add(dipol[1])

# Processing and filtering
for data in wenner:
    processing(data[1],  max_err = 10, smallest_rhoa = 10, biggest_rhoa = 2000,
               err_abs = 0.01, err_rel = .05, max_err_est =100,  setError = False, sort_out_error = 10000 )
    filter(data[1], filter_value = 0.8,  setError = False, sort_out_error = 10000)
    
adjust(wenner) # same quadropoles for all time steps

# Remove values that double or halve relative to the subsequent time step
for nr, data in enumerate(wenner):
    if nr > 0:
        wenner[nr-1][1].remove(data[1]['rhoa']/wenner[nr-1][1]['rhoa']>1.9)
adjust(wenner) # same quadropoles for all time steps
for nr, data in enumerate(wenner):
    if nr > 0:
        wenner[nr-1][1].remove(data[1]['rhoa']/wenner[nr-1][1]['rhoa']<0.5)
adjust(wenner) # same quadropoles for all time steps

# save filtered data 
for ordner_name, daten_objekt in wenner: 
    ordner_pfad = os.path.join("filtered_data", ordner_name) 
    os.makedirs(ordner_pfad, exist_ok=True)
    datei_pfad = os.path.join(ordner_pfad, "two_timesteps.ohm")
    daten_objekt.save(datei_pfad)
