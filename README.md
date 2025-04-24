# Data set | Time-Lapse ERT Inversion for the Monitoring of Soil Water Content at Urban Tree Sites 
This repository contains raw and processed data along with the Python scripts used to prepare the data in the manuscript.

Johannes Hoppenbrock,  Malkin Gerchow, Matthias Beyer, Vera J. Hörmann, Mona Quambusch, Michael W. Strohbach, Matthias Bücker: Time-Lapse ERT Inversion for the Monitoring of Soil Water Content at Urban Tree Sites, 2025.
Submitted to Water Resources Research.

If you find this data useful in your own research, please mention this data set and/or the manuscript.

# Summary
Urban trees play a crucial role in mitigating climate impacts through their cooling effects while also contributing other ecosystem services, like carbon sequestration. However, growth conditions for urban trees are challenging to begin with. Due to ongoing climate change, there is an urgent demand for understanding how these variations in precipitation and temperature affect soil moisture. Here we present monitoring data of soil moisture content at two urban tree sites using Electrical Resistivity Tomography (ERT). Monthly ERT measurements were conducted along rows of trees over a period of 1.5 years at one sealed and one unsealed site. For comparison and calibration, soil moisture sensors capable of measuring also temperature were installed at five depths between 15 cm and 2 meters. To eliminate temperature effects on the ERT data, a method was developed to model subsurface temperature variations down to 20 meters using the temperature data provided by the soil moisture sensors. A comparison of different time-lapse inversion techniques highlights their respective advantages and disadvantages and indicates that a full time-lapse inversion is the most suitable approach for our application. However, lower contact resistances at the sealed site favor artifacts in the inversion result. Due to small resistivity changes, careful filtering and error analysis are crucial. Overall, there is a greater dynamic in the changes of resistivity at the unsealed site compared to the sealed site, where most changes occur in the areas of the tree pits. Nonetheless, the employed setup allowed for the identification of periods with changes in groundwater levels, evaporation, and root water uptake.  
# Overview

The directory "Temperature correction" containts all files to reproduce the temperature correction part of the paper (section 4). The directory "Unsealed park site" contains all files necessary to reproduce the inversions and figures of the unsealed park site (section 5.1 and 6). The directory "Sealed site" contains all files necessary to reproduce the inversions and figures of the sealed site (section 5.2 and 6).

# Funding



