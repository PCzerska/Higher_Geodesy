from math import radians,degrees
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import math
from pyproj import Geod
import numpy as np
a= 6378137
e = 0.00669438002290

def kivioja(fip, lamp, azp, s):
    az= np.deg2rad(azp)
    lam= np.deg2rad(lamp)
    fi= np.deg2rad(fip)
    n= int(math.ceil(s/1.5))
    ds= np.empty(n)
    for i in range(0, n-1):
        ds[i]= 1.5
    ds[-1]= s-(n-1)*1.5


    for i in range(0, n):
        N= a / (np.sqrt(1-e*np.power(np.sin(fi), 2)))
        M= a*(1-e)/(np.sqrt(np.power(1-e*np.sin(fi)**2, 3)))

        dfi= ds[i]*np.cos(az)/M
        daz= ds[i]*np.sin(az)*np.tan(fi)/N

        azm= az+(1/2)*daz
        fim= fi +(1/2)*dfi

        Nm= a/(np.sqrt(1-e*np.power(np.sin(fim), 2)))
        Mm= a*(1-e)/(np.sqrt(np.power(1-e*np.sin(fim)** 2, 3)))
        print (Nm,Mm)

        dlamm = (ds[i] * np.sin(azm)) / (Nm * np.cos(fim))
        dazm = (ds[i] * np.sin(azm) * np.tan(fim)) / Nm
        dfim= (ds[i]*np.cos(azm))/Mm


        az = az + dazm
        lam = lam + dlamm
        fi = fi + dfim

    return [np.rad2deg(fi), np.rad2deg(lam), np.rad2deg(az)]

def vincenty(BA, LA, BB, LB):
    '''
    Parameters
    ----------
    BA : szerokosc geodezyjna punktu A [RADIAN]
    LA : dlugosc geodezyjna punktu A [RADIAN]
    BB : szerokosc geodezyjna punktu B [RADIAN]
    LB : dlugosc geodezyjna punktu B [RADIAN]

    Returns
    -------
    sAB : dlugosc linii geodezyjnej AB [METR]
    A_AB : azymut linii geodezyjnej AB [RADIAN]
    A_BA : azymut odwrotny linii geodezyjne [RADIAN]
    '''
    b = a * np.sqrt(1 - e)
    f = 1 - b / a
    dL = LB - LA
    UA = np.arctan((1 - f) * np.tan(BA))
    UB = np.arctan((1 - f) * np.tan(BB))
    L = dL
    while True:
        sin_sig = np.sqrt((np.cos(UB) * np.sin(L)) ** 2 + \
                          (np.cos(UA) * np.sin(UB) - np.sin(UA) * np.cos(UB) * np.cos(L)) ** 2)
        cos_sig = np.sin(UA) * np.sin(UB) + np.cos(UA) * np.cos(UB) * np.cos(L)
        sig = np.arctan2(sin_sig, cos_sig)
        sin_al = (np.cos(UA) * np.cos(UB) * np.sin(L)) / sin_sig
        cos2_al = 1 - sin_al ** 2
        cos2_sigm = cos_sig - (2 * np.sin(UA) * np.sin(UB)) / cos2_al
        C = (f / 16) * cos2_al * (4 + f * (4 - 3 * cos2_al))
        Lst = L
        L = dL + (1 - C) * f * sin_al * (sig + C * sin_sig * (cos2_sigm + C * cos_sig * (-1 + 2 * cos2_sigm ** 2)))
        if abs(L - Lst) < (0.000001 / 206265):
            break

    u2 = (a ** 2 - b ** 2) / (b ** 2) * cos2_al
    A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    d_sig = B * sin_sig * (cos2_sigm + 1 / 4 * B * (cos_sig * (-1 + 2 * cos2_sigm ** 2) \
                                                    -1 / 6 * B * cos2_sigm * (-3 + 4 * sin_sig ** 2) * (
                                                                -3 + 4 * cos2_sigm ** 2)))
    sAB = b * A * (sig - d_sig)
    A_AB = np.arctan2((np.cos(UB) * np.sin(L)), (np.cos(UA) * np.sin(UB) - np.sin(UA) * np.cos(UB) * np.cos(L)))
    A_BA = np.arctan2((np.cos(UA) * np.sin(L)),
                      (-np.sin(UA) * np.cos(UB) + np.cos(UA) * np.sin(UB) * np.cos(L))) + np.pi
    return sAB, A_AB, A_BA

#Dane
fipocz= 51
lampocz= 18
A12= 0
A23= 90
A34= 180
s12= 20000
s23= 35000
s34= 20000
s41= 35000


#Obliczenia dla punktu drugiego
fi2,lam2,A2= kivioja(fipocz,lampocz,A12,s12)
A2=degrees(A2)
print(fi2,lam2,A2)

#Obliczanie punktu trzeciego
fi3,lam3,A3= kivioja(fi2,lam2,A23,s23)
A3=degrees(A3)
print(fi3,lam3,A3)

#Obliczenia punktu czwartego
fi4,lam4,A4= kivioja(fi3,lam3,A34,s34)
A4=degrees(A4)
print(fi4,lam4,A4)

#Obliczenia dla punktu 1*
fi1b,lam1b,A1b=kivioja(fi4,lam4,270,s41)
A1b= degrees(A1b)
print(fi1b,lam1b,A1b)

print(vincenty(radians(fi4),radians(lam4),radians(fipocz),radians(lampocz)))
print(vincenty(radians(fipocz),radians(lampocz),radians(fi2),radians(lam2)))
print(vincenty(radians(fi2),radians(lam2),radians(fi3),radians(lam3)))
print(vincenty(radians(fi3),radians(lam3),radians(fi4),radians(lam4)))
geod = Geod(ellps="WGS84")
lons= [lampocz,lam2,lam3,lam4]
lats= [fipocz,fi2,fi3,fi4]
poly_area, poly_perimeter =geod.polygon_area_perimeter(lons,lats)
poly_area= poly_area*(-1)
print(poly_area, poly_perimeter)

fi=[fipocz,fi2,fi3,fi4]
lam=[lampocz,lam2,lam3,lam4]
plt.plot(lam,fi, marker='s', color= 'green')
plt.plot([lampocz,lam4],[fipocz,fi4],marker='s', color= 'green')
plt.title('Położenie punktów')
plt.ylabel('Szerokość geodezyjna w stopniach')
plt.xlabel('Długość geodezyjna w stopniach')
#plt.show()

fig = plt.figure(figsize=(40, 20))
request=cimgt.OSM()
extent= [lam1b-0.008, lam1b+0.008, lats[0]-0.008, lats[0]+0.008]
ax = plt.axes(projection=request.crs)
grid_lines= ax.gridlines(draw_labels=True, alpha=0.2)
grid_lines.xformatter=LONGITUDE_FORMATTER
grid_lines.yformatter=LATITUDE_FORMATTER
ax.set_extent(extent, crs=ccrs.Geodetic())
ax.add_image(request, 8)
ax.stock_img()
ax.coastlines()
plt.title('Położenie punktów na wykresie')
plt.plot(lons,lats, color='red')
plt.plot([lampocz,lam4],[fipocz,fi4], color= 'red')
plt.plot(lam1b,fi1b, 'o', transform=ccrs.PlateCarree(), color='purple', label= 'P1*')
plt.plot(lons[0], lats[0], 'o', transform=ccrs.PlateCarree(), color='orange', label= 'P1')
plt.plot(lons[1], lats[1], 'o', transform=ccrs.PlateCarree(), color='red', label= 'P2')
plt.plot(lons[2], lats[2], 'o', transform=ccrs.PlateCarree(), color='blue', label= 'P3')
plt.plot(lons[3], lats[3], 'o', transform=ccrs.PlateCarree(), color='green', label= 'P4')
plt.legend(fontsize=12)
#plt.show()

plt.figure(figsize=(12,12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.LAKES, alpha=0.5)
ax.add_feature(cfeature.RIVERS)
plt.plot(lons,lats, color='red')
plt.plot([lampocz,lam4],[fipocz,fi4], color= 'red')
ax.plot(lam1b,fi1b, 'o', transform=ccrs.PlateCarree(), color='purple', label= 'P1*')
ax.plot(lons[0], lats[0], 'o', transform=ccrs.PlateCarree(), color='orange', label= 'P1')
ax.plot(lons[1], lats[1], 'o', transform=ccrs.PlateCarree(), color='red', label= 'P2')
ax.plot(lons[2], lats[2], 'o', transform=ccrs.PlateCarree(), color='blue', label= 'P3')
ax.plot(lons[3], lats[3], 'o', transform=ccrs.PlateCarree(), color='green', label= 'P4')
ax.legend(fontsize=12)
ax.stock_img()
grid_lines= ax.gridlines(draw_labels=True, alpha=0.2)
grid_lines.xformatter=LONGITUDE_FORMATTER
grid_lines.yformatter=LATITUDE_FORMATTER
extent = [min(lons) - 1, max(lons) + 1, min(lats) - 1, max(lats) + 1]
ax.set_extent(extent, crs=ccrs.Geodetic())

#plt.show()