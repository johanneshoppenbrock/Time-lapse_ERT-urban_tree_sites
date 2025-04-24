#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import pygimli as pg
import inspect
from pygimli.physics import ert
import matplotlib.pyplot as plt
import numpy as np
import os 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from datetime import datetime 
from influxdb_client import InfluxDBClient
import pandas as pd


# In[2]:


def plotting_function_FTL (data, manager, results, folder_name, filename, chi2, scalef, rrms):
    fig, axs = plt.subplots(nrows=len(data), ncols=2, figsize=(10, len(data)*2.5))
    plt.subplots_adjust(hspace=0)
    axnr = -1
    axs[0, 1].axis('off')
    modell = results
    name = folder_name
    test = True 
    for counter, (d, res) in enumerate(zip(data, modell)):
        # Resistivity image
        axnr = axnr +1
        row_index = counter*2
        datum = datetime.strptime(d[0], '%y%m%d')
        bild1 = manager[0].showResult(model = res, coverage = manager[0].coverage(),
                cMin=10, cMax=300, logScale=True, cMap= 'jet', ax=axs[axnr, 0], label = 'Resistivity in $\Omega$m')
        axs[axnr, 0].set_title('Date: '+ str(datum.day) + '.' + str(datum.month) + '.' + str(datum.year))
        if axnr < len(axs)-1:
            cb = bild1[-1]
            cb.remove()
    
        # Ratio image
        if res[10] != modell[0][10]: 
            ratio = res / modell[counter-1]
            bild2 = manager[0].showResult(model= ratio,coverage = manager[0].coverage(),
                    cMin=1/2, cMax=2, cMap= 'bwr', ax= axs[axnr, 1], label = 'Change of resistivity', logScale= True)
            if axnr < len(axs)-1:
                cb = bild2[-1]
                cb.remove() 
            axs[axnr, 1].set_title('$\chi^2$: '+ str(chi2[0]) + '; rrms: '+ str(rrms[0])+ ' %',  loc = 'right')
            axs[axnr, 1].set_title('SF: '+ str(scalef), loc = 'left')
            tree_postions = [7, 16.5, 28, 37, 47]
            for x in tree_postions:
                axs[axnr, 1].plot(x, 0, '.', markersize=10, color='black')
            axs[axnr, 1].set_ylim(-15, 0) 
            axs[axnr, 0].set_ylim(-15, 0)
            axs[axnr, 1].axhline(y=-4.5, color='k', linewidth = 1,linestyle='--')
    for ax in axs:
        ax[0].set_xlabel('Distance (m)')
        ax[0].set_ylabel('Depth (m)')
        ax[1].set_xlabel('Distance (m)')
        ax[1].set_ylabel('Depth (m)')
    
    isExist = os.path.exists('./Inversion_results/' + '/%s_inversion' % name)
    if not isExist:
        os.makedirs('./Inversion_results/' + '/%s_inversion' % name)
    fig.savefig('./Inversion_results/' + '/%s_inversion' '/%s.png' % (name, filename), bbox_inches = 'tight', dpi=300)


# In[3]:


def plotting_function (data, manager, results, folder_name, filename, chi2, lam):
    fig, axs = plt.subplots(nrows=len(data), ncols=2, figsize=(10, len(data)*2.5))
    plt.subplots_adjust(hspace=0)
    axnr = -1
    axs[0, 1].axis('off')
    modell = results
    name = folder_name
    for counter, (d, res) in enumerate(zip(data, modell)):
        # Resistivity image
        axnr = axnr +1
        row_index = counter*2
        datum = datetime.strptime(d[0], '%y%m%d')
        bild1 = manager[0].showResult(model = res, coverage = manager[0].coverage(),
                cMin=10, cMax=300, logScale=True, cMap= 'jet', ax=axs[axnr, 0], label = 'Resistivity in $\Omega$m')
        axs[axnr, 0].set_title('Date: '+ str(datum.day) + '.' + str(datum.month) + '.' + str(datum.year))
        if axnr < len(axs)-1:
            cb = bild1[-1]
            cb.remove()
    
        # Ratio image
        if res[10] != modell[0][10]: 
            ratio = res / modell[counter-1]
            bild2 = manager[0].showResult(model= ratio,coverage = manager[0].coverage(),
                    cMin=1/2, cMax=2, cMap= 'bwr', ax= axs[axnr, 1], label = 'Change of resistivity', logScale= True)
            if axnr < len(axs)-1:
                cb = bild2[-1]
                cb.remove() 
            axs[axnr, 1].set_title('$\chi^2$ - 1: '+ str(chi2[counter-1]) + '; $\chi^2$ - 2: ' + str(chi2[counter]), loc = 'right')
            axs[axnr, 1].set_title('$\lambda = $'+ str(lam), loc = 'left')
            tree_postions = [7, 16.5, 28, 37, 47]
            for x in tree_postions:
                axs[axnr, 1].plot(x, 0, '.', markersize=10, color='black')
            axs[axnr, 1].set_ylim(-15, 0) 
            axs[axnr, 0].set_ylim(-15, 0)
            axs[axnr, 1].axhline(y=-4.5, color='k', linewidth = 1,linestyle='--')
    for ax in axs:
        ax[0].set_xlabel('Distance (m)')
        ax[0].set_ylabel('Depth (m)')
        ax[1].set_xlabel('Distance (m)')
        ax[1].set_ylabel('Depth (m)')
    
    isExist = os.path.exists('./Inversion_results/' + '/%s_inversion' % name)
    if not isExist:
        os.makedirs('./Inversion_results/' + '/%s_inversion' % name)
    fig.savefig('./Inversion_results/' + '/%s_inversion' '/%s.png' % (name, filename), bbox_inches = 'tight', dpi=300)


# In[4]:


def read(filename):
    wenner = []
    dd = []
    dataset = filename
    date = os.listdir("./Data/" + dataset)
    globals()['%s_wenner' % (dataset)] = [] # definiere Max_Osterloh_Platz_Wenner
    globals()['%s_datum' % (dataset)] = [] # definiere Max_Osterloh_Platz_datum
    globals()['%s_dd' % (dataset)]= [] # definiere Max_Osterloh_Platz_dd
    wenner = []
    dd = []
    for datum in date:
        dateien = os.listdir("./Data/" + dataset +'/'+ datum) # Dateien im Ordner des Datums einlesen
        globals()['%s_datum' % (dataset)].append(datum) # Datum zu Max_Osterloh_Platz_datum hinzufügen
        for datei in dateien: 
            if (datei).startswith('Wen'):
                wenner.append([datum, pg.load("./Data/" + dataset +'/'+ datum+'/'+datei)])
            if (datei).startswith('Dip'):
                dd.append([datum, pg.load("./Data/" + dataset +'/'+ datum+'/'+datei)])
    return wenner, dd


# In[5]:


def misfit_correction(wenner, inv, maximum_error): 
    """Filtering with respect to misfit  

    Parameters
    ----------
    wenner: Wenner dataset  
    inv: inversion instance
    
    Returns
    -------
    copy_misfit_corrected: copy of corrected data 
    """
    
    start = 0 
    for wen in wenner:
        misfit = - inv.response[start:(start+len(wen[1]['rhoa']))] / wen[1]["rhoa"] * 100 + 100
        print(misfit)
        manager[0].showData(misfit, cMap = 'bwr', cMin= -50, cMax = 50)
        start = start + len(wen[1]['rhoa'])
        
        for count,mis in enumerate(misfit):
            if mis > maximum_error:
                wen[1].markInvalid(count)
    for we in wenner: 
        we[1].removeInvalid()
    return wenner


# In[6]:


def T_corr_nach_Inversion(resistivity, mesh, Temp_tiefe_vektor):
    res_new = resistivity.copy()
    Tem_new = resistivity.copy()
    # Mesh nur mit Inversionsregion
    for count, (res, loc) in enumerate(zip(resistivity, mesh.cellCenter())) :
        #print(loc[1])
        if loc[1] > Temp_tiefe_vektor[0][0]:
            T = Temp_tiefe_vektor[0][1]
        else:
            for tiefe in Temp_tiefe_vektor:
                if  tiefe[0] + 0.025 > loc[1] > tiefe[0] - 0.025:
                    T = tiefe[1]
                    break
                else:
                    T = 10
         
        Tem_new[count] = T
        res_new[count] = (1 + 0.025 *(T-25))*res # temperature correction
    return res_new,Tem_new


# In[7]:


def adjust(data_list):
# reduziert Daten auf kleinstes Glied
    for k in [0,1]:
        for nummer in range(len(data_list)):
            data_timestep1 = data_list[nummer][1]
            if nummer == len(data_list)-1:
                data_timestep2 = data_list[0][1]
            else:
                data_timestep2 = data_list[nummer+1][1]
            #print(data_timestep1)
            for nr, a2,b2,m2,n2 in zip(range(len(data_timestep2['a'])),data_timestep2['a'], data_timestep2['b'],data_timestep2['m'], data_timestep2['n']):
                boolean = False 
                #print(nr)
                for a1,b1,m1,n1 in zip(data_timestep1['a'], data_timestep1['b'],data_timestep1['m'], data_timestep1['n']):
                    if a1==a2 and b2==b1 and m2==m1 and n2==n1:
                        boolean = True 
                if boolean == False:
                    #print(nr)
                    data_timestep2['valid'][nr]= 0 
            data_timestep2.removeInvalid()
            data_timestep1.removeInvalid()


# In[8]:


def processing(data, max_err, smallest_rhoa, biggest_rhoa, err_abs, err_rel, max_err_est, parameter = 'rhoa', setError=False, sort_out_error = 1):
    # Geometriefaktor berechnen
    data['k'] =ert.createGeometricFactors(data, numerical=True)
    # Bulk resistivity
    data['r'] = data[parameter]/data['k']
    # Filtering
    
    if setError: # Wenn Fehler hochgesetzt werden sollen 
        # Großen Fehlern einen noch größeren Fehler geben 
        err_index_list = []
        for nr, err in zip(range(len(data['rhoa']) - 1, -1, -1), reversed(data['err'])):
            if err*100 > max_err:
                err_index_list.append(nr)
        
        data['err'] = ert.estimateError(data, 
                            absoluteError=err_abs, 
                            relativeError=err_rel)
        for index in err_index_list:
            data['err'][index]= sort_out_error
    
    

        for nr, rhoa, err in zip(range(len(data['rhoa'])-1, -1 ,-1),reversed(data['rhoa']), reversed(data['err'])):
            if rhoa > biggest_rhoa or rhoa < smallest_rhoa:
                data['err'][nr]= sort_out_error
            if err*100 > max_err:
                data['err'][nr]= sort_out_error
    else:
        
        data.markInvalid(data['err']*100 > max_err)
        data.markInvalid(data[parameter] < smallest_rhoa)
        data.markInvalid(data[parameter] > biggest_rhoa)
        data['err'] = ert.estimateError(data, 
                    absoluteError=err_abs, 
                    relativeError=err_rel)
        
        data.markInvalid(data['err']*100 > max_err_est)

        
       #data.remove(data['err']*100 > max_err)
       #data.remove(data[parameter] < smallest_rhoa)
       #data.remove(data[parameter] > biggest_rhoa)
       #data['err'] = ert.estimateError(data, 
       #            absoluteError=err_abs, 
       #            relativeError=err_rel)
       #
       #data.remove(data['err']*100 > max_err_est)


# In[9]:


def einlesen_ip(filename):
    wenner = []
    dd = []
    ort = filename
    date = os.listdir("./Messdaten_ip/" + ort)
    print(date)
    globals()['%s_wenner' % (ort)] = [] # definiere Max_Osterloh_Platz_Wenner
    globals()['%s_datum' % (ort)] = [] # definiere Max_Osterloh_Platz_datum
    globals()['%s_dd' % (ort)]= [] # definiere Max_Osterloh_Platz_dd
    wenner = []
    dd = []
    for datum in date:
        dateien = os.listdir("./Messdaten_ip/" + ort +'/'+ datum) # Dateien im Ordner des Datums einlesen
        globals()['%s_datum' % (ort)].append(datum) # Datum zu Max_Osterloh_Platz_datum hinzufügen
        for datei in dateien: 
            if (datei).startswith('Wen'):
                wenner.append([datum, pg.load("./Messdaten_ip/" + ort +'/'+ datum+'/'+datei)])
            if (datei).startswith('Dip'):
                dd.append([datum, pg.load("./Messdaten_ip/" + ort +'/'+ datum+'/'+datei)])
    print('Schreibe zwei Listen (Wen, Dip) mit den jeweiligen Daten dazu.')
    return wenner, dd


# In[10]:


def einlesen_tagesgang(filename, date):
    wenner = []
    dd = []
    ort = filename
    globals()['%s_wenner' % (ort)] = [] # definiere Max_Osterloh_Platz_Wenner
    globals()['%s_datum' % (ort)] = [] # definiere Max_Osterloh_Platz_datum
    globals()['%s_dd' % (ort)]= [] # definiere Max_Osterloh_Platz_dd
    wenner = []
    dd = []
    uhrzeiten = os.listdir("./Messdaten_pg/" + ort +'/'+ date)
    for zeit in uhrzeiten:
        if zeit[0].isdigit():
            print('Uhrzeit: ' +zeit)
            dateien = os.listdir("./Messdaten_pg/" + ort +'/'+ date + '/' + zeit) # Dateien im Ordner des Datums einlesen
            globals()['%s_datum' % (ort)].append(zeit) # Datum zu Max_Osterloh_Platz_datum hinzufügen
            for datei in dateien: 
                if (datei).startswith('Wen'):
                    wenner.append([date+zeit, pg.load("./Messdaten_pg/" + ort +'/'+ date+'/'+ zeit +'/'+datei)])
                if (datei).startswith('Dip'):
                    dd.append([date+zeit, pg.load("./Messdaten_pg/" + ort +'/'+ date+'/'+ zeit +'/'+datei)])
            print('Schreibe zwei Listen (Wen, Dip) mit den jeweiligen Daten dazu.')
    return wenner, dd


# In[11]:


def zuordnung(wenner, dd):
    for dipol in dd:
        for wen in wenner: 
            if dipol[0] == wen[0]:
                wen[1].add(dipol[1])


# In[12]:


# Es werden Werte aussortiert bei denen für beide benachbarten Punkte gilt, dass rho_diff/rho > filter_value
#filter_value = 1 #Abstand/Wert
#data = wenner[3][1]

def filter(data, parameter = 'rhoa', filter_value = 1, setError = False, sort_out_error = 1):
# Filter Stärker, wenn filter_value kleiner
    for ind in range(len(data['a'])):
        #print(ind)
        if 0 < ind < len(data['a'])-1 and abs(data['a'][ind] - data['a'][ind-1]) < 5 and abs(data['a'][ind] - data['a'][ind+1]) < 5 :
            #print('jetzt: ' + str(data['a'][ind]))
            #print('vorheriger Wert: ' + str(data[parameter][ind-1]) + ' jetziger Wert: ' + str(data[parameter][ind]) + ' DIfferenz: ' + str(abs(data[parameter][ind] -  data[parameter][ind-1])) )
            if ((abs(data[parameter][ind] -  data[parameter][ind-1]) / data[parameter][ind] > filter_value and 
                 abs(data[parameter][ind] -  data[parameter][ind+1])/data[parameter][ind] > filter_value) or 
                (abs(data[parameter][ind] -  data[parameter][ind+1])/data[parameter][ind+1] > filter_value and
                 abs(data[parameter][ind] -  data[parameter][ind-1])/data[parameter][ind-1] > filter_value)):
                
                #print('pick')
                if setError == True:
                    data['err'][ind]= sort_out_error
                else:
                    data.markInvalid(ind)
    for ind in range(len(data['a'])):
        if 0 < ind < len(data['a'])-1 and 0 < ind < len(data['a'])-1 and abs(data['a'][ind] - data['a'][ind-1]) > 5 :
            # Anfang Reihexyvc fe
            if abs(data[parameter][ind] -  data[parameter][ind+1])/data[parameter][ind] > filter_value:
                if setError == True:
                    data['err'][ind]= sort_out_error
                else:
                    data.markInvalid(ind)

        if 0 < ind < len(data['a'])-1 and 0 < ind < len(data['a'])-1 and abs(data['a'][ind] - data['a'][ind+1]) > 5: 
            # Ende Reihe 
            if abs(data[parameter][ind] -  data[parameter][ind-1])/data[parameter][ind] > filter_value:
                if setError == True:
                    data['err'][ind]= sort_out_error          
                
                else:
                    data.markInvalid(ind)
    
    #data.removeInvalid()


# In[13]:


# gefilterte Daten abspeichern
def save_filtered_data(folder_new, ort, wenner, dd):
    datum = 0 
    for wenner2 in wenner:
        if datum is not wenner2[0]:
            nummer = 1
        else:
            nummer = nummer + 1 
            
        datum = wenner2[0]
        isExist = os.path.exists(folder_new + '/' + ort +'/'+ datum)
        if not isExist:
           # Create a new directory because it does not exist
            os.makedirs(folder_new + '/' + ort +'/'+ datum)
            print("The new directory is created!")
        [wenner2[1].save(folder_new + '/' + ort +'/'+ datum+'/Wenner%s.ohm' % str(nummer))]
    datum = 0 
    for dd2 in dd:
        if datum is not dd2[0]:
            nummer = 1
        else:
            nummer = nummer + 1 
            
        datum = dd2[0]
        [dd2[1].save(folder_new + '/' + ort +'/'+ datum+'/DipDip%s.ohm' % str(nummer))]


# In[14]:


def resistivity_area(mesh, mod, x_min, x_max, y_min, y_max):
    """Widerstände in definiertem Bereich und Median der Widerstände

    Parameters
    ----------
    x_min: x-Minimum of area of interest 
    x_max: x-Maximum of area of interest 
    y_min: y-Minimum of area of interest (e.g. -0.75)
    y_max: y-Maximum of area of interest (e.g. -0.25)
    mod: Liste der Inversionsergebnisse 
    mesh: Mesh, das zum Datensatz passt

    Returns
    -------
    tiefenverlauf : Für jeden Zeitschritt Liste von Widerständen im Bereich. 
    """
    tiefenverlauf = []
    cell_centers = np.array([cell.center() for cell in mesh.cells()])
    mask = (cell_centers[:,0] > x_min) & (cell_centers[:,0] < x_max) & (cell_centers[:,1] > y_min) & (cell_centers[:,1] < y_max)
    area = cell_centers[mask]
    for count, m in enumerate(mod): 
        para = pg.Mesh(mesh)  # make a copy
        para.setCellMarkers(pg.IVector(para.cellCount()))
        fopDP = pg.frameworks.PriorModelling(para, area)
        tiefenverlauf.append(fopDP(m))
    return tiefenverlauf


# In[15]:


def extract_sensor_data(standort, wenner, Uhrzeit_start, Uhrzeit_ende, list_of_depth):
    """Feuchtesensordaten auslesen und in Dataframes ausgeben

    Parameters
    ----------
    standort: Standort der vermessen wurde 
    wenner: Zeitreihe mit Daten 
    Uhrzeit_start, Uhrzeit_ende: Zeitraum in dem nach Feuchtesensordaten gesucht wird, meist 1h 

    Returns
    -------
    data_frame_TDR_list : Liste mit Dataframes für jede Tiefe.
    """
    start_date = []
    end_date = []
    ticks = []
    data_frame_list = []
    datum = [eintrags_tuple[0] for eintrags_tuple in wenner]
    for tag in datum: 
        start_date.append('%sT%s:00:00Z' % (str((datetime.strptime(tag, '%y%m%d')).date()),Uhrzeit_start))
        end_date.append('%sT%s:00:00Z' % (str((datetime.strptime(tag, '%y%m%d')).date()),Uhrzeit_ende))
        ticks.append(str((datetime.strptime(tag, '%y%m%d')).date()))
        
    for i, depth in enumerate(list_of_depth):
        data_frame = None
        for start, end in zip(start_date, end_date):
            if data_frame is None: 
                with InfluxDBClient(
                    url="https://database.isodrones.de",
                    token="AWAE4XarwTdDkgePJp-hvLDqV-lKBuUYecVc4RFqYpddJxD_-ICSK8E_AUksqSJ6i5_qp8e_QHk9D5B4kxvs9Q==",
                    org="isodrones",
                    debug=False,
                ) as client:
                    query_api = client.query_api()
                    data_frame = (query_api.query_data_frame('''
                    from(bucket:"climax2") 
                        |> range(start: %s, stop: %s) 
                        |> filter(fn: (r) => r["location"] == "depth_cm_%s")
                        |> filter(fn: (r) => r["site"] == "%s")
                        |> pivot(rowKey:["_time"], columnKey: ["_field"], valueColumn: "_value") 
                    ''' % (start, end,str(depth), standort))) 
                if data_frame.empty:
                    # Eine Zeile mit NaNs erzeugen
                    data_frame = pd.DataFrame(
                        [[np.nan]*len(data_frame.columns)],  # Liste mit NaN-Werten
                        columns=data_frame.columns
                    )  

                continue
                      
            with InfluxDBClient(
                url="https://database.isodrones.de",
                token="AWAE4XarwTdDkgePJp-hvLDqV-lKBuUYecVc4RFqYpddJxD_-ICSK8E_AUksqSJ6i5_qp8e_QHk9D5B4kxvs9Q==",
                org="isodrones",
                debug=False,
            ) as client:
                query_api = client.query_api()
                zusatz = (query_api.query_data_frame('''
                from(bucket:"climax2") 
                    |> range(start: %s, stop: %s) 
                    |> filter(fn: (r) => r["location"] == "depth_cm_%s")
                    |> filter(fn: (r) => r["site"] == "%s")
                    |> pivot(rowKey:["_time"], columnKey: ["_field"], valueColumn: "_value") 
                ''' % (start, end,str(depth), standort)))
            try:
                # Versuch, die erste Zeile aus 'zusatz' anzuhängen
                data_frame = pd.concat([data_frame, zusatz.iloc[[0]]], ignore_index=True)
            except:
                # Falls das schiefgeht (z.B. 'zusatz' ist leer oder Index uncached),
                # füge stattdessen eine Zeile mit NaN hinzu
                nan_row = pd.DataFrame([[np.nan]*len(data_frame.columns)], columns=data_frame.columns)
                data_frame = pd.concat([data_frame, nan_row], ignore_index=True)

    
        data_frame_list.append(data_frame)
    data_frame_TDR_list = []
    for data_frame in data_frame_list:
        if 'BulkEC_mS/m' in data_frame.columns:
            data_frame_TDR_list.append(data_frame[['_time', 'location', 'Temperature_°C', 'WaterContent_%vol','BulkEC_mS/m']].copy())
            continue
        data_frame_TDR_list.append(data_frame[['_time', 'location', 'Temperature_°C', 'WaterContent_%vol']].copy())
    return data_frame_TDR_list


# In[16]:


def extract_sensor_data_continuous(standort, wenner, Uhrzeit_start, Uhrzeit_ende, list_of_depth):
    """Feuchtesensordaten auslesen und in Dataframes ausgeben

    Parameters
    ----------
    standort: Standort der vermessen wurde 
    wenner: Zeitreihe mit Daten 
    Uhrzeit_start, Uhrzeit_ende: Zeitraum in dem nach Feuchtesensordaten gesucht wird, meist 1h 

    Returns
    -------
    data_frame_TDR_list : Liste mit Dataframes für jede Tiefe.
    """
    start_date_strs = []
    end_date_strs = []
    datum = [eintrags_tuple[0] for eintrags_tuple in wenner]

    # Berechne den frühesten Start und den spätesten Endzeitpunkt
    earliest_start = min(datum)
    latest_end = max(datum)
    
    # Formatiere die Start- und Endzeitpunkte
    start_date_str = '%sT%s:00:00Z' % (str((datetime.strptime(earliest_start, '%y%m%d')).date()), Uhrzeit_start)
    end_date_str = '%sT%s:00:00Z' % (str((datetime.strptime(latest_end, '%y%m%d')).date()), Uhrzeit_ende)
    
    data_frame_list = []
    
    for depth in list_of_depth:
        with InfluxDBClient(
            url="https://database.isodrones.de",
            token="AWAE4XarwTdDkgePJp-hvLDqV-lKBuUYecVc4RFqYpddJxD_-ICSK8E_AUksqSJ6i5_qp8e_QHk9D5B4kxvs9Q==",
            org="isodrones",
            debug=False,
        ) as client:
            query_api = client.query_api()
            data_frame = query_api.query_data_frame('''
            from(bucket:"climax2") 
                |> range(start: %s, stop: %s) 
                |> filter(fn: (r) => r["location"] == "depth_cm_%s")
                |> filter(fn: (r) => r["site"] == "%s")
                |> pivot(rowKey:["_time"], columnKey: ["_field"], valueColumn: "_value") 
            ''' % (start_date_str, end_date_str,str(depth), standort)) 
    
        data_frame_list.append(data_frame)
    
    data_frame_TDR_list = []
    for data_frame in data_frame_list:
        if 'BulkEC_mS/m' in data_frame.columns:
            data_frame_TDR_list.append(data_frame[['_time', 'location', 'Temperature_°C', 'WaterContent_%vol','BulkEC_mS/m']].copy())
        else:
            data_frame_TDR_list.append(data_frame[['_time', 'location', 'Temperature_°C', 'WaterContent_%vol']].copy())
    
    return data_frame_TDR_list


# In[17]:


def Baumstandorte(Ort):
    if Ort == 'Prinzenpark':
        Baumkoordinaten = [1.5, 10.5, 17.5, 26.5, 34, 43]
        return Baumkoordinaten
    if Ort == 'Fontanestrasse':
        Baumkoordinaten = [7.5,18,29.5,41]
        return Baumkoordinaten
    if Ort == 'Museumspark_neu':
        Baumkoordinaten = [6.4, 11.2, 21]
        return Baumkoordinaten
    if Ort == 'Langer_Kamp':
        Baumkoordinaten = [7,17.6,28.6,38.2,47.4]
        return Baumkoordinaten
    if Ort == 'Schillstrasse':
        Baumkoordinaten = [1.0,13.3,25.5,37.5]
        return Baumkoordinaten
    if Ort == 'Georg_Westermann_Allee':
        Baumkoordinaten = [2,10.5, 18.5, 32, 42]
        return Baumkoordinaten


# In[18]:


def Feuchtesensoren(Ort):
    if Ort == 'Fontanestrasse':
        sensorpositionen = [9.8,9.8,9.8,9.8,9.8]
        return sensorpositionen
    if Ort == 'Museumspark_neu':
        sensorpositionen = [7.8,7.8,7.8,7.8,7.8]
        return sensorpositionen
    if Ort == 'Prinzenpark':
        sensorpositionen = [27,27,27,27,27]
        return sensorpositionen
    if Ort == 'Langer_Kamp':
        sensorpositionen = [39.1,39.1,39.1,39.1,39.1]
        return sensorpositionen
    if Ort == 'Schillstrasse':
        sensorpositionen = [38,38,38,38,38]
        return sensorpositionen
    if Ort == 'Georg_Westermann_Allee':
        sensorpositionen = [33,33,33,33,33]
        return sensorpositionen


# In[19]:


get_ipython().system('jupyter nbconvert --to script Hilfsfunktionen.ipynb')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




