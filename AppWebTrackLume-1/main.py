import shutil
import time
from geopy.distance import geodesic
from flask import Flask, render_template
import json
import urllib.request
from skyfield.api import Loader, EarthSatellite
from skyfield.api import wgs84
import numpy as np
from datetime import datetime, date, timedelta
from pathlib import Path
from math import copysign
import os
import codecs
import serial
import webbrowser

halfpi, pi, twopi = [f*np.pi for f in [0.5, 1, 2]]
degs, rads = 180/pi, pi/180
load = Loader('~/Documents/fishing/SkyData')
data = load('de421.bsp')
ts = load.timescale()
now = datetime.now()
path = Path("/AppWebTrackLume-1/FilesPasses/PredictionPasses")
backfiles = Path("/AppWebTrackLume-1/FilesPasses/PassedPasses")
Teleco_Uvigo = wgs84.latlon(+42.169952849844274, -8.68756907481494)
start = True
LUME = None
back_time = None
latlngs_json1 = None
latlngs_json2 = None
next_pass_ini = None
next_pass_mid = None
next_pass_end = None
time_next_pass_start = None
time_next_pass_finish = None
range_sat = None

class DatosPrediccion:
    def __init__(self, ti, alt, az, distance):
        self.ti = ti
        self.alt = alt
        self.az = az
        self.distance = distance

def process_TLEs():
    # Carga el archivo desde el URL
    with urllib.request.urlopen('https://celestrak.com/NORAD/elements/cubesat.txt') as f:
       datos = f.read().decode('utf-8')
    lines = [line.rstrip() for line in datos.splitlines()]
    class Satelite:
        def __init__(self, nombre, linea1, linea2):
            self.nombre = nombre
            self.linea1 = linea1
            self.linea2 = linea2
    satelites = []
    # Indices
    x, y, z = 0, 1, 2
    # Itera cada 3 elementos
    for i in range(int(len(lines) / 3)):
        # Creo objeto y lo agrego a lista
        satelites.append(Satelite(lines[x], lines[y], lines[z]))
        # Incremento indices
        x += 3
        y += 3
        z += 3
    # Convierte objeto satelite en json
    datosJson = [json.dumps(s.__dict__) for s in satelites]
    # Analiza json
    analizado = [json.loads(datosJson[i]) for i in range(len(datosJson))]
    sateliteDict = {}
    for satelite in analizado:
        sateliteDict[satelite["nombre"]] ={"linea1":satelite["linea1"],"linea2":satelite["linea2"]}
    linea1=sateliteDict["LUME 1"].get("linea1")
    linea2=sateliteDict["LUME 1"].get("linea2")
    LUME = EarthSatellite(linea1, linea2)
    return LUME

def desfase_utc():
   h_local = time.localtime().tm_hour
   h_utc = time.gmtime().tm_hour
   if h_local == 0 or h_local == 1:
       desfase = h_local+24 - h_utc
   else:
       desfase = h_local - h_utc
   return desfase

def prediction_orbit(LUME):
    # print("-----------PREDICCION ORBITA---------------")
    now_orbit = datetime.now()
    latlngs_array = []
    des_utc= desfase_utc()
    orbit_longitude_deg_prev = None
    for ti in range(120):  # 120 minutos
        if ti == 0:
            resta_1h = now_orbit - timedelta(hours=0)
            incremento_minuto = resta_1h + timedelta(minutes=ti)
            t_stop = ts.utc(incremento_minuto.year, incremento_minuto.month, incremento_minuto.day, incremento_minuto.hour, incremento_minuto.minute, incremento_minuto.second)
            t_now = ts.now()
            orbit_time = LUME.at(t_now)
            orbit_pos = wgs84.subpoint(orbit_time)
            orbit_longitude_rad = orbit_pos.longitude.radians
            orbit_latitude_rad = orbit_pos.latitude.radians
            orbit_longitude_deg = orbit_longitude_rad * (180 / np.pi)
            orbit_latitude_deg = orbit_latitude_rad * (180 / np.pi)
            latlngs_str = ("[" + str(round(orbit_latitude_deg, 4)) + "," + str(round(orbit_longitude_deg, 4)) + "]")
            latlngs_array.append(latlngs_str)
            orbit_longitude_deg_prev = orbit_longitude_deg
            continue
        resta_1h = now_orbit - timedelta(hours=des_utc)
        incremento_minuto = resta_1h + timedelta(minutes=ti)
        t_stop = ts.utc(incremento_minuto.year, incremento_minuto.month, incremento_minuto.day, incremento_minuto.hour, incremento_minuto.minute, incremento_minuto.second)
        orbit_time = LUME.at(t_stop)
        orbit_pos = wgs84.subpoint(orbit_time)
        orbit_longitude_rad = orbit_pos.longitude.radians
        orbit_latitude_rad = orbit_pos.latitude.radians
        orbit_longitude_deg = orbit_longitude_rad * (180 / np.pi)
        orbit_latitude_deg = orbit_latitude_rad * (180 / np.pi)
        if (copysign(1, orbit_longitude_deg) == copysign(1, 1)) and (copysign(1, orbit_longitude_deg_prev) == copysign(1,-1)):
            latlngs_str_aux = ("]&[[" + str(round(orbit_latitude_deg, 4)) + "," + str(round(orbit_longitude_deg, 4)) + "]")
            latlngs_array.append(latlngs_str_aux)
        else:
            latlngs_str_aux = ("[" + str(round(orbit_latitude_deg, 4)) + "," + str(round(orbit_longitude_deg, 4)) + "]")
            latlngs_array.append(latlngs_str_aux)
        orbit_longitude_deg_prev = orbit_longitude_deg
    latlngs_str = str(latlngs_array)
    latlngs_str = latlngs_str.replace(", ']&[", "]&[")
    characters = "'"
    latlngs_aux = ''.join(x for x in latlngs_str if x not in characters)
    latlngs = latlngs_aux.split("&")
    latlngs_json1 = json.loads(latlngs[0])
    latlngs_json2 = json.loads(latlngs[1].lstrip(", "))
    return latlngs_json1, latlngs_json2

def write_in_txt(path_file,t_ini,t_end,LUME):
    duration_prediction=timedelta(t_end-t_ini).seconds
    difference_prediction = LUME - Teleco_Uvigo
    t_ini_str=str(t_ini.utc_strftime('%Y %b %d %H:%M:%S.4f'))
    t_ini_datetime=datetime.strptime(t_ini_str, '%Y %b %d %H:%M:%S.4f')
    #Calculamos todas las posiciones del pase y las escribimos en el fichero
    file = open(path_file, 'w')
    for t_aux in range(duration_prediction):
        t_increase = t_ini_datetime + timedelta(seconds=t_aux)
        t_step_prediction = ts.utc(t_increase.year, t_increase.month, t_increase.day, t_increase.hour,t_increase.minute, t_increase.second)
        topocentric_prediction = difference_prediction.at(t_step_prediction)
        alt_prediction, az_prediction, distance_prediction = topocentric_prediction.altaz()
        file.write(str(t_step_prediction.utc_strftime('%Y %b %d %H:%M:%S'))+" Azimuth = "+str(round(az_prediction.degrees,3))+" deg Elevation = "+str(round(alt_prediction.degrees,3))+" deg Distance = "+str(round(distance_prediction.km,2))+" Km\n")
    file.close()

def send_pass_to_txt(passIni,passMid,passEnd,LUME):
    duration=str(timedelta(passEnd.ti-passIni.ti))
    duration= duration.split(":")
    min_dur= duration[1]
    sec_dur=duration[2]
    name=str(passIni.ti.utc_strftime('%Y %b %d %Hh %Mm %Ss')+" to "+passEnd.ti.utc_strftime('%Y %b %d %Hh %Mm %Ss')+" - dur="+min_dur+"m "+str(round(float(sec_dur)))+"s elev_max="+str(round(passMid.alt))+"º az_ini="+str(round(passIni.az))+"º az_end="+str(round(passEnd.az)))+"º.txt"
    path_with_file = path.joinpath(name).resolve()
    write_in_txt (path_with_file,passIni.ti,passEnd.ti,LUME)
    
def update_files():
    dirs = os.listdir(path)
    for fichero in dirs:
        name_file = str(fichero)
        path_file = path.joinpath(name_file).resolve()
        name_file = name_file.split("to")
        ini_pass = name_file[0].replace("h "," ")
        ini_pass = ini_pass.replace("m "," ")
        ini_pass = ini_pass.replace("s ","")
        date_time_obj = datetime. strptime(ini_pass,'%Y %b %d %H %M %S') #date_time_obj en UTC
        if date_time_obj > datetime.utcnow(): #los pases posteriores UTC a la hora actual UTC son los que tengo que eliminar para que se guarden los recalculos
            os.unlink(path_file)
        else:
            shutil.move(str(path_file),str(backfiles))

def prediction_passes(LUME):
    now_prediction = datetime.now()
    #print("-------------- PREDICCION PASES --------------")
    #Cuando un satélite se eleva por encima de 0 grados del horizonte y cuando se pone
    today = date.today()
    tomorrow = today + timedelta(days=3)
    des_utc = desfase_utc()
    t0 = ts.utc(now_prediction.year, now_prediction.month, now_prediction.day, now_prediction.hour - des_utc, now_prediction.minute, now_prediction.second)
    t1 = ts.utc(tomorrow.year, tomorrow.month, tomorrow.day)
    t, events = LUME.find_events(Teleco_Uvigo, t0, t1, altitude_degrees=0.0) #pases del satélite
    difference_prediccion = LUME - Teleco_Uvigo
    ObjetoDatosPrediccion = []
    global ultimo
    contador = 0
    next_pass_ini = None
    next_pass_mid = None
    next_pass_end = None
    t_fin = None
    t_ini = None
    ti_backup=None
    alt_backup=None
    az_backup=None
    distance_backup=None
    ultimo=False
    for ti, event in zip(t, events):
        topocentric_prediccion = difference_prediccion.at(ti)
        alt, az, distance = topocentric_prediccion.altaz()
        if ultimo == True:
            ObjetoDatosPrediccion.append(DatosPrediccion(ti,alt.degrees,az.degrees,distance.km))#fin del pase
            if contador == 1:
                t_fin = ti
                next_pass_end = ti.utc_strftime('%Y %b %d %Hh %Mm %Ss') + " UTC with azimuth " + str(
                    round(az.degrees, 4))
            ultimo = False
        if alt.degrees >= 10: #pase optimo
            contador = contador + 1
            ObjetoDatosPrediccion.append(DatosPrediccion(ti_backup,alt_backup.degrees,az_backup.degrees,distance_backup.km)) #inicio de pase
            ObjetoDatosPrediccion.append(DatosPrediccion(ti,alt.degrees,az.degrees,distance.km)) #cuando el pase esta mas proximo a teleco
            if contador == 1:
                next_pass_ini = ti_backup.utc_strftime('%Y %b %d %Hh %Mm %Ss') + " UTC with azimuth " + str(round(az_backup.degrees, 4))
                t_ini = ti_backup
                next_pass_mid = ti.utc_strftime('%Y %b %d %Hh %Mm %Ss')+" UTC with max. elevation "+str(round(alt.degrees,4))
            ultimo = True
        else:
            alt_backup = alt
            az_backup = az
            distance_backup = distance
            ti_backup = ti
            ultimo = False
    for x in range(0,len(ObjetoDatosPrediccion),3):
        send_pass_to_txt(ObjetoDatosPrediccion[x],ObjetoDatosPrediccion[x+1],ObjetoDatosPrediccion[x+2], LUME)
    return next_pass_ini, next_pass_mid, next_pass_end, t_ini, t_fin

def satellite_range(t_start, t_finish):
    geocentric = LUME.at(t_start)
    subpoint = wgs84.subpoint(geocentric)
    longitude_rad_start = subpoint.longitude.radians
    latitude_rad_start = subpoint.latitude.radians
    longitude_start = float(longitude_rad_start * (180 / np.pi))
    latitude_start = float(latitude_rad_start * (180 / np.pi))
    geocentric = LUME.at(t_finish)
    subpoint = wgs84.subpoint(geocentric)
    longitude_rad_finish = subpoint.longitude.radians
    latitude_rad_finish = subpoint.latitude.radians
    longitude_finish = float(longitude_rad_finish * (180 / np.pi))
    latitude_finish = float(latitude_rad_finish * (180 / np.pi))
    START = (latitude_start, longitude_start)
    FINISH = (latitude_finish, longitude_finish)
    radio = (geodesic(START, FINISH).m)/2

    return radio

def satellite_position_in_real_time(LUME):
    #print("-------------- Posicionamiento del satelite en tiempo real --------------")
    t = ts.now()
    geocentric = LUME.at(t)
    subpoint = wgs84.subpoint(geocentric)
    #print(subpoint.longitude, subpoint.latitude)  # datos en WGS84
    longitude_rad = subpoint.longitude.radians
    latitude_rad = subpoint.latitude.radians
    #print(longitude_rad, latitude_rad)  # datos en radianes
    longitude_deg = longitude_rad * (180 / np.pi)
    latitude_deg = latitude_rad * (180 / np.pi)
    #print(longitude_deg, latitude_deg)  # datos en grados
    elev = round(subpoint.elevation.km, 4)
    # Indicar si el satelite esta en zona de luz o de sombra
    sunlit = LUME.at(t).is_sunlit(data)
    if sunlit == True:
        light_shadow="sunny area"
    else:
        light_shadow="shaded area"
    return longitude_deg, latitude_deg, elev, light_shadow

def communication_RS232(alt_rx,az_rx):
    alt=str(alt_rx)
    az=str(az_rx)
    ser = serial.Serial(port='COM1', baudrate=9600, bytesize=serial.EIGHTBITS, parity=serial.PARITY_NONE, stopbits=serial.STOPBITS_ONE, timeout=0.5)  # ruta del dispositivo
    print("Connected to: " + ser.portstr)
    ser.write(codecs.encode(alt))
    a = ser.readline()
    print(codecs.decode(a).rstrip('\n'))
    ser.close()
    ser.open()
    ser.write(codecs.encode(az))
    a = ser.readline()
    print(codecs.decode(a).rstrip('\n'))
    ser.close()

app = Flask(__name__)

@app.route('/')

def show_view():
    return render_template("LUME-1.html")

@app.route('/ajax', methods=["GET"])
def ajax():
    t_clock = ts.now()
    global start, LUME, back_time, latlngs_json1, latlngs_json2, range_sat, next_pass_ini, next_pass_mid, next_pass_end, time_next_pass_finish, time_next_pass_start
    if start == True: # Descargar TLEs y prediccion de pase y orbita en cuanto se ejecute el programa
        back_time=t_clock
        LUME = process_TLEs()
        update_files()
        next_pass_ini, next_pass_mid, next_pass_end, time_next_pass_start, time_next_pass_finish = prediction_passes(LUME)
        latlngs_json1, latlngs_json2 = prediction_orbit(LUME)
        start = False
    #Parámetros orbitales
    last_pass = "Last TLE: " + str(LUME)
    numsat = "The unique satellite NORAD catalog number given in the TLE file: " + str(LUME.model.satnum)
    classification = "Satellite classification, or else 'U' for “Unknown: " + str(LUME.model.classification)
    designator = "International designator: "+str(LUME.model.intldesg)
    epoch_year = "Full four-digit year of this element set’s epoch moment: " + str(LUME.model.epochyr)
    epoch_day = "Fractional days into the year of the epoch moment: "+str(LUME.model.epochdays)
    julian_date = "Julian date of the epoch: " + str(LUME.model.jdsatepoch)
    f_derivate = "First time derivative of the mean motion: " + str(LUME.model.ndot)
    s_derivate= "Second time derivative of the mean motion: " + str(LUME.model.nddot)
    ballistic = "Ballistic drag coefficient B* in inverse earth radio: " + str(LUME.model.bstar)
    type_eph = "Ephemeris type: " + str(LUME.model.ephtype)
    n_element = "Element number: " + str(LUME.model.elnum)
    incl = "Inclination in radians: " + str(LUME.model.inclo)
    ascension = "Right ascension of ascending node in radians: " + str(LUME.model.nodeo)
    ecc = "Eccentricity: " + str(LUME.model.ecco)
    perigeo = "Argument of perigee in radians: " + str(LUME.model.argpo)
    anomaly = "Mean anomaly in radians: " + str(LUME.model.mo)
    motion = "Mean motion in radians per minute: " + str(LUME.model.no_kozai)
    revolution = "Revolution number at epoch [Revs]: " + str(LUME.model.revnum)
    duration_aux=str(timedelta(time_next_pass_finish-time_next_pass_start))
    duration_aux= duration_aux.split(":")
    min_dur_aux= duration_aux[1]
    sec_dur_aux=duration_aux[2]
    duration_pass=str("Duration of the satellite pass: "+min_dur_aux+"m "+str(int(float(sec_dur_aux)))+"s")
    #Posicionamiento satélite en tiempo real
    longitude_deg, latitude_deg, elev, light_shadow = satellite_position_in_real_time(LUME)
    #Huella del satélite
    range_sat = satellite_range(time_next_pass_start, time_next_pass_finish)
    '''
    print("LAT: " + str(latitude_deg)+"º")
    print("LNG: " + str(longitude_deg)+"º")
    print("ALT: " + str(elev)+" Km")
    print(light_shadow)
    '''
    #Datos a enviar al rotor en tiempo real
    difference = LUME - Teleco_Uvigo
    topocentric = difference.at(t_clock)
    alt, az, distance = topocentric.altaz()
    azimuth = az.degrees
    altitud = alt.degrees
    dist=distance.km
    count =str(timedelta(time_next_pass_start - t_clock))
    count=count.split(":")
    h_count = count[0]
    min_count = count[1]
    sec_count = count[2]
    time_to_np ="Next pass of the LUME-1 satellite starts in "+ str(h_count) + "h " + str(min_count) + "m "+ sec_count[0:5]+ "s"
    if alt.degrees > 0 and t_clock.utc > time_next_pass_start.utc: #se encuentra el pase activo
        time_to_np="Current pass of the LUME-1 satellite"
        communication_RS232(alt.degrees, az.degrees)
        '''
        print('Elevation:', alt.degrees)
        print('Azimuth:', az.degrees)
        print('Distance: {:.1f} km'.format(distance.km))
        '''
    else:
        accumulator_time = timedelta(t_clock-back_time)
        if (accumulator_time.seconds > 1800) or (t_clock.utc > time_next_pass_finish.utc): #cada 30 minutos ->1800 seg h actualizo las predicciones
            back_time = t_clock
            LUME = process_TLEs()
            update_files()
            latlngs_json1, latlngs_json2 = prediction_orbit(LUME)
            next_pass_ini, next_pass_mid, next_pass_end, time_next_pass_start, time_next_pass_finish = prediction_passes(LUME)
    '''
    print(LUME)
    print("The unique satellite NORAD catalog number given in the TLE file: "+str(LUME.model.satnum))
    print("Satellite classification, or else 'U' for “Unknown: "+str(LUME.model.classification))
    print("International designator: "+str(LUME.model.intldesg))
    print("Full four-digit year of this element set’s epoch moment: "+str(LUME.model.epochyr))
    print("Fractional days into the year of the epoch moment: "+str(LUME.model.epochdays))
    print("Julian date of the epoch: " + str(LUME.model.jdsatepoch))
    print("First time derivative of the mean motion: " + str(LUME.model.ndot))
    print("Second time derivative of the mean motion: " + str(LUME.model.nddot))
    print("Ballistic drag coefficient B* in inverse earth radio: " + str(LUME.model.bstar))
    print("Ephemeris type: " + str(LUME.model.ephtype))
    print("Element number: " + str(LUME.model.elnum))
    print("Inclination in radians: " + str(LUME.model.inclo))
    print("Right ascension of ascending node in radians: " + str(LUME.model.nodeo))
    print("Eccentricity: " + str(LUME.model.ecco))
    print("Argument of perigee in radians: " + str(LUME.model.argpo))
    print("Mean anomaly in radians: " + str(LUME.model.mo))
    print("Mean motion in radians per minute: " + str(LUME.model.no_kozai))
    print("Revolution number at epoch [Revs]: " + str( LUME.model.revnum))
    '''
    #Visualizacion de datos
    return {
        "latitud": round(latitude_deg, 4),
        "longitud": round(longitude_deg, 4),
        "elevacion": elev,
        "azimuth": round(azimuth, 4),
        "altitud": round(altitud, 4),
        "dist": round(dist, 4),
        "light_shadow": light_shadow,
        "range_sat": range_sat,
        "next_pass_ini": next_pass_ini,
        "next_pass_mid": next_pass_mid,
        "next_pass_end": next_pass_end,
        "time_to_np": time_to_np,
        "duration_pass": duration_pass,
        "latlngs1": latlngs_json1,
        "latlngs2": latlngs_json2,
        "numsat": numsat,
        "classification": classification,
        "designator": designator,
        "epoch_year": epoch_year,
        "epoch_day": epoch_day,
        "julian_date": julian_date,
        "f_derivate": f_derivate,
        "s_derivate": s_derivate,
        "ballistic": ballistic,
        "type_eph": type_eph,
        "n_element": n_element,
        "incl": incl,
        "ascension": ascension,
        "ecc": ecc,
        "perigeo": perigeo,
        "anomaly": anomaly,
        "motion": motion,
        "revolution": revolution,
        "last_pass": last_pass
    }

if __name__ == '__main__':
    if not os.environ.get("WERKZEUG_RUN_MAIN"):
        webbrowser.open_new('http://127.0.0.1:5000/')
    app.run(host="0.0.0.0", port=5000, debug=True)