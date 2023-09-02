import numpy as np
import pandas as pd
from math import sqrt,sin,radians,cos,atan2,degrees
from datetime import datetime
from datetime import time
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt

a= 6378137 #wielka półoś
e2= 0.006694380022 #kwadrat pierwszego mimośrodu dla elipsoidy GRS80

def FunkcjaDoZamianyWysokosci(h):
    h = h * 0.3048
    h = h + 135.4
    return h

def zamiana(fi,lam,h):
    N = a / (sqrt(1 - (e2*(sin(radians(fi)))**2)))
    X = (N+h)*cos(radians(fi))*cos(radians(lam))
    Y = (N+h)*cos(radians(fi))*sin(radians(lam))
    Z = (N*(1-e2)+h)*sin(radians(fi))
    return X, Y, Z

def wektor(X, Y, Z):
    return [X,Y,Z]

def obldlwek(a,b,c):
    dl= sqrt((a**2)+(b**2)+(c**2))
    return dl



#FunkcjaDoOdczytu
def read_flightradar(file):
    '''
    Parameters
    ----------
    file : .csv file - format as downloaded from fligthradar24
        DESCRIPTION.
    Returns
    -------
    all_data : numpy array
        columns are:
            0 - Timestamp - ?
            1 - year
            2 - month
            3 - day
            4 - hour
            5 - minute
            6 - second
            7 - Latitude [degrees]-fi
            8 - Longitude [degrees]-lambda
            9 - Altitude [feet]- wysokosc
            10 - Speed [?]
            11 - Direction [?]
    '''

    with open(file, 'r') as f:
        i = 0
        size= []
        Timestamp = []; date = []; UTC = []; Latitude = []; Longitude = [];
        Altitude = []; Speed = []; Direction = []; datetime_date = []

        for linia in f:
            if linia[0:1] != 'T':
                splited_line = linia.split(',')
                size.append(len(splited_line))
                i += 1
                Timestamp.append(int(splited_line[0]))
                full_date = splited_line[1].split('T')
                date.append(list(map(int, full_date[0].split('-'))))
                UTC.append(list(map(int, full_date[1].split('Z')[0].split(':'))))
                Callsign = splited_line[2]
                Latitude.append(float(splited_line[3].split('"')[1]))
                Longitude.append(float(splited_line[4].split('"')[0]))
                Altitude.append(float(splited_line[5]))
                Speed.append(float(splited_line[6]))
                Direction.append(float(splited_line[7]))

    all_data = np.column_stack((np.array(Timestamp), np.array(date), np.array(UTC),
                                np.array(Latitude), np.array(Longitude), np.array(Altitude),
                                np.array(Speed), np.array(Direction)))
    return all_data
all_data= read_flightradar('D:/studia/3sem/projekty/projekt2/LO79_2dd99519.csv')
fi = all_data[:,7]
lam = all_data[:,8]
filot=all_data[0,7]
lamlot=all_data[0,8]

#Zmienione wysokości na metry i z uwzględnioną wysokością lotniska
h= list(map(FunkcjaDoZamianyWysokosci, all_data[:,9]))

#Zamiana współrzędnych
new= map(zamiana,fi,lam,h)
df=pd.DataFrame(new)

#Obliczenia wektora samolot-lotnisko
#wspolrzedna XYZ wektora lotniska
Xl= df[0][0]
Yl= df[1][0]
Zl= df[2][0]
print(Xl,fi[0])
Xs= df[0][1:]
Ys= df[1][1:]
Zs= df[2][1:]

XYZ_l= np.column_stack((Xl,Yl,Zl))
XYZ= np.column_stack((Xs,Ys,Zs))

wektorX=np.subtract(XYZ,XYZ_l)


#Obliczenie wektora u
al= cos(radians(fi[0]))*cos(radians(lam[0]))
b= cos(radians(fi[0]))*sin(radians(lam[0]))
c= sin(radians(fi[0]))
u=wektor(al,b,c)

#Obliczenie wektora n
d= -sin(radians(fi[0]))*cos(radians(lam[0]))
e= -sin(radians(fi[0]))*sin(radians(lam[0]))
f= cos(radians(fi[0]))
n=wektor(d,e,f)

#Obliczenie wektora e
g= -sin(radians(lam[0]))
h3= cos(radians(lam[0]))
i= 0
e=wektor(g,h3,i)


#Tworzenie macierzy r
r=(wektor(n,e,u))


#Mnożenie macierzy
wekx=[]
for l in range(0,len(wektorX)):
    a=np.array(r) @ wektorX[l]
    wekx.append(a)

wekxd=pd.DataFrame(wekx)


#Obliczenie długości wektora s i azymutu oraz h
s=[]
az=[]
zz=[]
h2=[]
for ev in range(0,len(wekx)):
    zma=wekx[ev][0]
    zmb=wekx[ev][1]
    zmc=wekx[ev][2]
    dl=sqrt(zma**2 + zmb**2 + zmc**2)
    az.append(atan2(zmb, zma))
    s.append(dl/1000)
    hh=np.arcsin(zmc/dl)
    h2.append(hh)
    zz.append(radians(90)-hh)

start=datetime(int(all_data[0][1]),int(all_data[0][2]),int(all_data[0][3]),int(all_data[0][4]),int(all_data[0][5]),int(all_data[0][6]))
times=[datetime(int(all_data[i][1]),int(all_data[i][2]),int(all_data[i][3]),int(all_data[i][4]),int(all_data[i][5]),int(all_data[i][6])) - start for i in range (1,len(all_data))]

for i in range(len(h2)):
    if h2[i] <= 0 and h[i+1]>135.4:
        point = i
        lamh=lam[i]
        fih=fi[i]
        godziny = times[i]
        print(godziny)
        break

fig = plt.figure(figsize=(40, 20))
request=cimgt.OSM()
extent= [min(lam)-5, max(lam)+5, min(fi)-5, max(fi)+5]
ax = plt.axes(projection=request.crs)
ax.set_extent(extent, crs=ccrs.Geodetic())
ax.add_image(request, 5)
ax.stock_img()
ax.coastlines()
fi=fi[1:]
lam=lam[1:]
lat=np.insert(fi,0,filot)
long=np.insert(lam,0,lamlot)
ax.plot(lamlot, filot, 'o', transform=ccrs.PlateCarree(), color='orange', label= 'Warszawa')
ax.plot(lam[851],fi[851],'o', transform=ccrs.PlateCarree(), color= 'orange', label= 'Tokio')
ax.plot(lamh,fih,'o', transform=ccrs.PlateCarree(), color= 'red', label='Punkt zniknięcia z horyzontu')
ax.plot([lamlot, lam[851]], [filot, fi[851]], transform=ccrs.PlateCarree(), color= 'orange', label= 'Prosta pomiędzy lotniskami')
ax.plot([lamlot, lam[851]], [filot, fi[851]], transform=ccrs.Geodetic(), color= 'purple', label= 'Idealny przebieg trasy')
ax.plot(long, lat ,transform=ccrs.PlateCarree(),color='b', label= "Rzeczywista trasa samolotu")
ax.legend(fontsize=12)
#plt.show()




#Wykresy zależności od czasu


from datetime import time
times=np.array([time.seconds+(24*3600)*time.days for time in times])
times=times/3600


#Wysokość samolotu od czasu
h=h[1:]
fig1=plt.figure(figsize=(10,5))
plt.title('Wykres zmiany wysokości w metrach:')
plt.plot(times,h, color='green')
plt.xlabel('Czas od startu w godzinach:')
plt.ylabel('Wysokość w metrach:')
plt.show()


#Prędkość od czasu
speed= all_data[1:,10]
speed= speed * 1.85
fig2=plt.figure(figsize=(10,5))
plt.title('Wykres zmiany prędkości samolotu od czasu')
plt.plot(times, speed, color='red')
plt.xlabel('Czas od startu w godzinach:')
plt.ylabel('Prędkość samolotu w km/h:')
plt.show()

#Odległość samlotu od lotniska od czasu
fig3=plt.figure(figsize=(10,5))
plt.title('Wykres zmiany odległości samolotu od lotniska od czasu')
plt.plot(times, s, color='red')
plt.xlabel('Czas od startu w godzinach:')
plt.ylabel('Odległość samolotu od lotniska w kilometrach:')
plt.show()

# Skyplot nad lotniskiem w Warszawie
fig4 = plt.figure(figsize = (8,8))
ax = fig4.add_subplot(polar = True)
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.scatter(az[:point],s[:point], color='red')
ax.set_title('Mapa nieba na lotnisku Chopina',fontdict={'fontname': 'Georgia', 'fontsize': 18})
plt.show()

print(point)
h2=np.array(h2)

for index in range(0, len(h2)):
    h2[index]= degrees(h2[index])

#Wykres zależności elewacji od czasu

fig5=plt.figure(figsize=(20,5))
plt.plot(times,h2, color= 'green')
plt.title('Wykres zmiany elewacji od czasu:')
plt.xlabel('Czas od startu w godzinach:')
plt.ylabel('Elewacja w stopniach:')
plt.show()
