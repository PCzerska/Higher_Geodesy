import numpy as np
from math import asin,sin,cos,radians,degrees,atan2
import matplotlib.pyplot as plt

def dms2deg(dms):
    d = dms[0]
    m = dms[1]
    s = dms[2]
    deg = d + m / 60 + s / 3600
    return deg

def rad2deg(p):
    s=0
    for i in p:
        p[s]=degrees(i)
        s+=1
    return p

def deg2rad(p):
    s=0
    for i in p:
        p[s]=radians(i)
        s+=1
    return p

#Obliczenia na dzień 1 lipca 2022
#Krok 1
year= 2022
month= 7
day= 1
hours= np.arange(24)
hours= hours-2
hours[0]+= 24
hours[1]+= 24
#print(hours)

#Krok 2
def julday(year,month,day,hours):
    if month<= 2:
        year-= 1
        month+= 12
    jd = np.floor(365.25 * (year + 4716)) + np.floor(30.6001 * (month + 1)) + day + hours / 24 - 1537.5
    return jd

jd= julday(year,month,day,hours)
jd[0]= julday(year, 6, 30, 22)
jd[1]= julday(year, 6, 30, 23)
#print(jd)

def GMST(jd):
    T = (jd - 2451545) / 36525
    Tu = jd - 2451545
    g = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * T * 2 - T * 3 / 38710000
    g = (g % 360) / 15
    return g

gmst = GMST(jd)
#print(gmst)

dl= 21
lst= gmst + dl/15
#print(lst)

for i in range(len(lst)):
    if lst[i] >= 24:
        lst[i] -= 24
#print(lst)


#Krok3
m=16
s=41.299
alfa = 14 + m/60 + s/3600
t= (lst-alfa)*15
for l in range(len(t)):
    if t[l]<0:
        t[l]+= 360
#print(t)

#Krok4 rozwiazywanie trojkata sferycznego

t=deg2rad(t)
#print(t)

deklW=radians(dms2deg([19,3,58.070]))
szerW= radians(52)

#Obliczenie wysokosci Warszawa
hW=[]
for y in t:
    hW.append(asin(sin(szerW)*sin(deklW)+cos(szerW)*cos(deklW)*cos(y)))
#print(hW)

#Obliczanie wysokosci Równik

hR=[]
for z in t:
    hR.append((asin(sin(radians(0))*sin(deklW)+cos(radians(0))*cos(deklW)*cos(z))))
#print(hR)


#Obliczenie azymutu Warszawa
aW=[]
for a in t:
    aW.append(atan2((-cos(deklW)*sin(a)),(cos(szerW)*sin(deklW)-sin(szerW)*cos(deklW)*cos(a))))
#print(aW)

#Obliczanie azymutu Równik
aR=[]
for b in t:
    aR.append(atan2((-cos(deklW) * sin(b)),(cos(radians(0)) * sin(deklW) - sin(radians(0)) * cos(deklW) * cos(b))))
print(aR)


#Zamiana na stopnie
hW= rad2deg(hW)
aW= rad2deg(aW)
hR= rad2deg(hR)
aR= rad2deg(aR)

#WYKRESY!!!
t=np.arange(24)
plt.plot(t, hW, 'ro', label= 'Warszawa')
plt.plot(t, hR, 'g^', label= 'Równik')
plt.title('Wykres zmian wysokości gwiazdy w ciągu doby', fontdict={'fontname': 'Georgia', 'fontsize': 18})
plt.xlabel('Godzina', fontdict={'fontname': 'Georgia', 'fontsize': 14})
plt.ylabel('Wysokość w stopniach', fontdict={'fontname': 'Georgia', 'fontsize': 14})
plt.xticks(np.arange(24))
plt.legend()
plt.show()
plt.plot(t, aW, 'ro', label= 'Warszawa')
plt.plot(t, aR, 'g^', label= 'Równik')
plt.title('Wykres zmian azymutu gwiazdy w ciągu doby', fontdict={'fontname': 'Georgia', 'fontsize': 18})
plt.xlabel('Godzina', fontdict={'fontname': 'Georgia', 'fontsize': 14})
plt.ylabel('Azymut w stopniach', fontdict={'fontname': 'Georgia', 'fontsize': 14})
plt.xticks(np.arange(24))
plt.legend()
plt.show()

#Zamiana na radiany
hW= deg2rad(hW)
aW= deg2rad(aW)
hR= deg2rad(hR)
aR= deg2rad(aR)


#wykres 3D
fig = plt.figure(figsize = (10,10))
ax = fig.add_subplot(projection = '3d')
# promień Ziemi
r = 1
# siatka wspołrzędnych
u, v = np.mgrid[0:(2 * np.pi+0.1):0.1, 0:np.pi:0.1]
x = np.cos(u) * np.sin(v)
y = np.sin(u) * np.sin(v)
z = np.cos(v)
z[z<0] = 0		# bez tego, narysowalibyśmy całą kulę, a chcemy tylko półkulę
ax.plot_surface(x,y,z, alpha = 0.1)

# narysowanie punktu na sferze
gx = r * np.sin(aW) * np.cos(hW)
gy = r * np.cos(aW) * np.cos(hW)
gz = r * np.sin(hW)
gx1= r * np.sin(aR) * np.cos(hR)
gy1= r * np.cos(aR) * np.cos(hR)
gz1= r * np.sin(hR)
ax.plot3D(gx,gy,gz,'ro',label='Warszawa')
ax.plot3D(gx1,gy1,gz1,'g^',label='Równik')
ax.set_title('Wykres zmian położenia gwiazdy w ciągu doby',fontdict={'fontname': 'Georgia', 'fontsize': 18})
ax.legend()


fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(polar = True)
ax.set_theta_zero_location('N') # ustawienie kierunku północy na górze wykresu
ax.set_theta_direction(-1)

ax.set_yticks(range(0, 90+10, 10))                   # Define the yticks

yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
ax.set_yticklabels(yLabel)
ax.set_rlim(0,90)

# narysowanie punktu na wykresie
ax.scatter(aW, 90-np.rad2deg(hW), color='red', label='Warszawa')
ax.set_title('Mapa nieba ',fontdict={'fontname': 'Georgia', 'fontsize': 18})
ax.scatter(aR, 90-np.rad2deg(hR), color='green', label='Równik')
ax.legend()
plt.show()