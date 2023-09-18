from math import sin,sqrt,cos,tan, radians, atan2,degrees
import pyproj
#Dane punktow
fi1= 51
lam1= 18
fi2= 51.179775285095815
lam2= 18
fi3= 51.17870466807602
lam3=18.500520984238573
fi4= 50.9989293498725
lam4=18.500520984238573

a=6378137
e2= 0.00669438002290
ep=0.082094
m020 = 0.999923
m092= 0.9993


A0= 1- (e2/4) -(3*(e2**2)/64)-(5*(e2**3)/256)
A2= (3/8) * (e2 +((e2)**2)/4 + (15*(e2**3)/128))
A4= (15/256)*((e2**2)+(3*(e2**3))/4)
A6= 35*(e2**3)/3072




def xy(fi,lam,lam0):
    fi=radians(fi)
    lam=radians(lam)-radians(lam0)
    k= a* (A0*fi -A2*sin(2*fi) + A4*sin(4*fi) - A6*sin(6*fi))
    N= a/sqrt(1-e2*(sin(fi)**2))
    t= tan(fi)
    ni2= (ep**2) * (cos(fi)**2)
    xgk= k + (lam**2/2) * N * sin(fi) * cos(fi) * (
                1 + (lam**2/12) * (cos(fi)**2) * (5-t**2 + 9*ni2 + 4*ni2** 2)
                + ((lam**2/360) * (cos(fi)**4) *
                   (61 - 58*t**2 + t**4 +270*ni2 - 330*ni2*t**2)))
    ygk= lam*N*cos(fi)*(1 + ((lam**2/6) * (cos(fi)**2) * (1-t**2 + ni2)) +
                                          ((lam**4/120) * (cos(fi)**4) *
                                           (5 - 18*t**2 + t**4 + 14*ni2 - 58*ni2*t**2)))

    return xgk, ygk

def pl2000(x,y):

    x2000= m020*x
    y2000= m020*y + 6* 1000000 + 500000
    return x2000,y2000

def pl1992(x,y):

    x1992= m092*x - 5300000
    y1992= m092*y + 500000
    return x1992,y1992

def pitagoras(x1,x2,y1,y2):
    dx=abs(x2-x1)
    dy=abs(y2-y1)
    dl=sqrt(dx**2+dy**2)
    return dl

def redukcja(s,fi1,fi2,y1,y2):
    fi= (fi1+fi2)/2
    fi= radians(fi)
    M= (a*(1-e2))/(sqrt(1-e2*(sin(fi)**2)**3))
    N= a/sqrt(1-e2*(sin(fi)**2))

    R=sqrt(N*M)
    r= s* (y1**2+y1*y2+y2**2)/(6*R**2)
    return r

def po_redukcji(dl,r):
    sel= dl-r
    return sel

def alfa(x1,y1,x2,y2):
    alfa=atan2(y2-y1,x2-x1)
    return alfa

def zbieznosc_pol1(lam,fi,lam0):
    lam=radians(lam-lam0)
    fi=radians(fi)

    ni2= (ep**2) * (cos(fi)**2)
    t2= (tan(fi))**2
    zb=lam*sin(fi) + ((lam**3/3)*sin(fi)*(cos(fi)**2)*(1+3*ni2+2*(ni2**2))) + \
       ((lam**5/15)*sin(fi)*(cos(fi)**4)*(2-t2))
    return zb

def redukcja_kierunkow(x1,y1,x2,y2,fi1,fi2):
    fi=(fi1+fi2)/2
    fi = radians(fi)
    M = (a * (1 - e2)) / (sqrt(1 - e2 * (sin(fi) ** 2) ** 3))
    N = a / sqrt(1 - e2 * (sin(fi) ** 2))
    R = sqrt(N * M)

    rk12=(x2-x1)*(2*y1+y2)/(6*(R**2))
    rk21=(x1-x2)*(2*y2+y1)/(6*(R**2))
    return rk12,rk21

def azymuty_el(kat,zb,r):
    azel=kat+zb+r
    azel=degrees(azel)
    return azel

def pole_gaus(X,Y):
    P= X[0]*(Y[1]-Y[3])+X[1]*(Y[2]-Y[0])+X[2]*(Y[3]-Y[1])+X[3]*(Y[0]-Y[2])
    P=P/2
    return P

def pole_kontrola(X,Y):
    P = Y[0] * (X[1] - X[3]) + Y[1] * (X[2] - X[0]) + Y[2] * (X[3] - X[1]) + Y[3] * (X[0] - X[2])
    P = P / (-2)
    return P

def przeliczanie_wspolrzednych(pun2000):
    source_crs= 'epsg:2177'
    target_crs= 'epsg:3035'
    Pl_to_laea = pyproj.Transformer.from_crs(source_crs, target_crs)
    Xl=[]
    Yl=[]
    for ev in pun2000:
        xl,yl=Pl_to_laea.transform(ev[0], ev[1])
        Xl.append(xl)
        Yl.append(yl)
    return Xl,Yl


def zbieznosc2(x, y):
    fi1= x/(a*A0)
    while True:
        k= a*(A0*fi1-A2*sin(2*fi1)+A4*sin(4*fi1)-A6*sin(6*fi1))
        phi1n= fi1+(x-k)/(a*A0)
        dln= a*(A0*phi1n-A2*sin(2*phi1n)+A4*sin(4*phi1n)-A6*sin(6*phi1n))
        if (abs(phi1n-fi1)<=radians(0.000001/3600)):
            fi1= phi1n
            k= dln
            break
        else:
            fi1= phi1n
            k= dln
    N1= a/sqrt(1-e2*((sin(fi1))**2))
    t= tan(fi1)
    ni2= ep**2*((cos(fi1))**2)
    zb= (y/N1)*t*(1-(y**2/(3*N1**2))*
                    (1+t**2-ni2-2*ni2**2)+(y**4/(15*N1**4))*
                   (2+5*t**2+3*t**4))
    return zb


#-----------------------------------------------------------------------------------------------------------------------
#X,Y_2000
x1,y1=xy(fi1,lam1,18)
x2,y2=xy(fi2,lam2,18)
x3,y3=xy(fi3,lam3,18)
x4,y4=xy(fi4,lam4,18)

# print(str(x1),str(y1))
# print(str(x2),str(y2))
# print(str(x3),str(y3))
# print(str(x4),str(y4))
# print("Punkt 1"+str((x1,y1)))
# print("Punkt 2"+str((x2,y2)))
# print("Punkt 3"+str((x3,y3)))
# print("Punkt 4"+str((x4,y4)))
#X,Y_1992
x11,y11=xy(fi1,lam1,19)
x21,y21=xy(fi2,lam2,19)
x31,y31=xy(fi3,lam3,19)
x41,y41=xy(fi4,lam4,19)
# print("Punkt 1"+str((x11,y11)))
# print("Punkt 2"+str((x21,y2)))
# print("Punkt 3"+str((x31,y31)))
# print("Punkt 4"+str((x41,y41)))
#PL-2000
x12000,y12000=pl2000(x1,y1)
x22000,y22000=pl2000(x2,y2)
x32000,y32000=pl2000(x3,y3)
x42000,y42000=pl2000(x4,y4)
X2000=[x12000,x22000,x32000,x42000]
Y2000=[y12000,y22000,y32000,y42000]
pun2000=[[x12000,y12000],[x22000,y22000],[x32000,y32000],[x42000,y42000]]
# print("Wspołrzędne X,Y PL-2000")
# print("Punkt 1"+str((x12000,y12000)))
# print("Punkt 2"+str((x22000,y22000)))
# print("Punkt 3"+str((x32000,y32000)))
# print("Punkt 4"+str((x42000,y42000)))
#PL-1992
x11992,y11992=pl1992(x11,y11)
x21992,y21992=pl1992(x21,y21)
x31992,y31992=pl1992(x31,y31)
x41992,y41992=pl1992(x41,y41)
X1992=[x11992,x21992,x31992,x41992]
Y1992=[y11992,y21992,y31992,y41992]
# print("Wspołrzędne X,Y PL-1992")
# print("Punkt 1"+str((x11992,y11992)))
# print("Punkt 2"+str((x21992,y21992)))
# print("Punkt 3"+str((x31992,y31992)))
# print("Punkt 4"+str((x41992,y41992)))

#----------------------------------------------------------------------------------------------------------------------
#Obliczanie długości w 2000
dl12=pitagoras(x12000,x22000,y12000,y22000)
dl23=pitagoras(x22000,x32000,y22000,y32000)
dl34=pitagoras(x32000,x42000,y32000,y42000)
dl41=pitagoras(x42000,x12000,y42000,y12000)
#print(dl12,dl23,dl34,dl41)
# print("Długości PL-2000")
# print('1-2',dl12)
# print('2-3',dl23)
# print('3-4',dl34)
# print('4-1',dl41)

#Obliczanie długości w 1992
dl122=pitagoras(x11992,x21992,y11992,y21992)
dl232=pitagoras(x21992,x31992,y21992,y31992)
dl342=pitagoras(x31992,x41992,y31992,y41992)
dl412=pitagoras(x41992,x11992,y41992,y11992)
# print(dl122,dl232,dl342,dl412)
# print("Długości PL-1992")
print('1-2',dl122)
print('2-3',dl232)
print('3-4',dl342)
print('4-1',dl412)

#Obliczanie długośći w GK
s12=dl12/m020
s23=dl23/m020
s34=dl34/m020
s41=dl41/m020

s122=dl122/m092
s232=dl232/m092
s342=dl342/m092
s412=dl412/m092

#Redukcje długości
r12=redukcja(s12,fi1,fi2,y1,y2)
r23=redukcja(s23,fi2,fi3,y2,y3)
r34=redukcja(s34,fi3,fi4,y3,y4)
r41=redukcja(s41,fi4,fi1,y4,y1)
#print(r12,r23,r34,r41)
r121=redukcja(s122,fi1,fi2,y11,y21)
r231=redukcja(s232,fi2,fi3,y21,y31)
r341=redukcja(s342,fi3,fi4,y31,y41)
r411=redukcja(s412,fi4,fi1,y41,y11)

#DŁugości
s12el=po_redukcji(s12,r12)
s23el=po_redukcji(s23,r23)
s34el=po_redukcji(s34,r34)
s41el=po_redukcji(s41,r41)
#print(s12el,s23el,s34el,s41el)
# print('PL-2000')
# print('1-2',s12el)
# print('2-3',s23el)
# print('3-4',s34el)
# print('4-1',s41el)

s121el=po_redukcji(s122,r121)
s231el=po_redukcji(s232,r231)
s341el=po_redukcji(s342,r341)
s411el=po_redukcji(s412,r411)


# print('\n')
print('PL-1992')
print('1-2',s121el)
print('2-3',s231el)
print('3-4',s341el)
print('4-1',s411el)

#--------------------------------------------------------------------------------------------------------
#Kąty kierunkowe wprost
alf12=degrees(alfa(x12000,y12000,x22000,y22000))
alf23=degrees(alfa(x22000,y22000,x32000,y32000))
alf34=degrees(alfa(x32000,y32000,x42000,y42000))
alf41=degrees(alfa(x42000,y42000,x12000,y12000))+360
print("Kąty kierunkowe wprost")
print(alf12,alf23,alf34,alf41)
#Kąty kierunkowe odwrotnie
alf21=degrees(alfa(x22000,y22000,x12000,y12000))
alf32=degrees(alfa(x32000,y32000,x22000,y22000))
alf43=degrees(alfa(x42000,y42000,x32000,y32000))
alf14=degrees(alfa(x12000,y12000,x42000,y42000))-180
print("Kąty kierunkowe odwrotne")
print(alf21,alf32,alf43,alf14)
#Zbieżność południków na podstawie współrzędnych geodezyjnych
zb1=zbieznosc_pol1(lam1,fi1,18)
zb2=zbieznosc_pol1(lam2,fi2,18)
zb3=zbieznosc_pol1(lam3,fi3,18)
zb4=zbieznosc_pol1(lam4,fi4,18)



#Zbieżność południków na podstawie współrzędnych płaskich
zbp1=zbieznosc2(x1, y1)
zbp2=zbieznosc2(x2, y2)
zbp3=zbieznosc2(x3, y3)
zbp4=zbieznosc2(x4, y4)

#Redukcja azymutów kątów wprost i odwrotnych
rk12,rk21=redukcja_kierunkow(x2,y1,x2,y2,fi1,fi2)
rk23,rk32=redukcja_kierunkow(x2,y2,x3,y3,fi2,fi3)
rk34,rk43=redukcja_kierunkow(x3,y3,x4,y4,fi3,fi4)
rk41,rk14=redukcja_kierunkow(x4,y4,x1,y1,fi4,fi1)


#Azymuty na elipsie
#PL-2000
az12=azymuty_el(radians(alf12),zb1,rk12)
az21=azymuty_el(radians(alf21),zb2,rk21)

az23=azymuty_el(alf23,zb2,rk23)
az32=azymuty_el(alf32,zb3,rk32)

az34=azymuty_el(alf34,zb3,rk34)
az43=azymuty_el(alf43,zb4,rk43)

az41=azymuty_el(alf41,zb4,rk41)
az14=azymuty_el(alf14,zb1,rk41)
# print('az12 '+str(az12))
# print('az21 '+str(az21))
# print('az23 '+str(az23))
# print('az32 '+str(az32))
# print('az34 '+str(az34))
# print('az43 '+str(az43))
# print('az14 '+str(az14))
# print('az41 '+str(az41))

#druga wersja
az12p=azymuty_el(radians(alf12),zbp1,rk12)
az21p=azymuty_el(radians(alf21),zbp2,rk21)

az23p=azymuty_el(radians(alf23),zbp2,rk23)
az32p=azymuty_el(radians(alf32),zbp3,rk32)

az34p=azymuty_el(radians(alf34),zbp3,rk34)
az43p=azymuty_el(radians(alf43),zbp4,rk43)

az41p=azymuty_el(radians(alf41),zbp4,rk41)
az14p=azymuty_el(radians(alf14),zbp1,rk41)


# print('az12 '+str(az12p))
# print('az21 '+str(az21p))
# print('az23 '+str(az23p))
# print('az32 '+str(az32p))
# print('az34 '+str(az34p))
# print('az43 '+str(az43p))
# print('az14 '+str(az14p))
# print('az41 '+str(az41p))

#PL-1992
alf122=(alfa(x11992,y11992,x21992,y21992))
alf232=(alfa(x21992,y21992,x31992,y31992))
alf342=(alfa(x31992,y31992,x41992,y41992))
alf412=(alfa(x41992,y41992,x11992,y11992))
print("Kąty kierunkowe wprost")
print(alf122,alf232,alf342,alf412)
alf212=(alfa(x21992,y21992,x11992,y11992))
alf322=(alfa(x31992,y31992,x21992,y21992))
alf432=(alfa(x41992,y41992,x31992,y31992))
alf142=(alfa(x11992,y11992,x41992,y41992))
print("Kąty kierunkowe odwrotne")
print(alf212,alf322,alf432,alf142)


zb11=zbieznosc_pol1(lam1,fi1,19)
zb21=zbieznosc_pol1(lam2,fi2,19)
zb31=zbieznosc_pol1(lam3,fi3,19)
zb41=zbieznosc_pol1(lam4,fi4,19)

zbp11=zbieznosc2(x11, y11)
zbp21=zbieznosc2(x21, y21)
zbp31=zbieznosc2(x31, y31)
zbp41=zbieznosc2(x41, y41)

rk121,rk211=redukcja_kierunkow(x11,y11,x21,y21,fi1,fi2)
rk231,rk321=redukcja_kierunkow(x21,y21,x31,y31,fi2,fi3)
rk341,rk431=redukcja_kierunkow(x31,y31,x41,y41,fi3,fi4)
rk411,rk141=redukcja_kierunkow(x41,y41,x11,y11,fi4,fi1)


az121=azymuty_el(radians(alf122),zb11,rk121)
az211=azymuty_el(radians(alf212),zb21,rk211)+360

az231=azymuty_el(radians(alf232),zb21,rk231)
az321=azymuty_el(radians(alf322),zb31,rk321)+360

az341=azymuty_el(radians(alf342),zb31,rk341)+360
az431=azymuty_el(radians(alf432),zb41,rk431)

az411=azymuty_el(radians(alf412),zb41,rk411)+360
az141=azymuty_el(radians(alf142),zb11,rk411)



# print('az12 '+str(az121))
# print('az21 '+str(az211))
# print('az23 '+str(az231))
# print('az32 '+str(az321))
# print('az34 '+str(az341))
# print('az43 '+str(az431))
# print('az14 '+str(az141))
# print('az41 '+str(az411))


az121p=azymuty_el(alf122,zbp11,rk121)
az211p=azymuty_el(alf212,zbp21,rk211)+360

az231p=azymuty_el(alf232,zbp21,rk231)
az321p=azymuty_el(alf322,zbp31,rk321)+360

az341p=azymuty_el(alf342,zbp31,rk341)+360
az431p=azymuty_el(alf432,zbp41,rk431)

az411p=azymuty_el(alf412,zbp41,rk411)+360
az141p=azymuty_el(alf142,zbp11,rk411)

print('az12 '+str(az121p))
print('az21 '+str(az211p))
print('az23 '+str(az231p))
print('az32 '+str(az321p))
print('az34 '+str(az341p))
print('az43 '+str(az431p))
print('az14 '+str(az141p))
print('az41 '+str(az411p))

#-----------------------------------------------------------------------------------------------------------------------------
#Pole PL-2000
P2000=pole_gaus(X2000,Y2000)
print('Pole PL-2000 w m^2: '+str(P2000))
Pk2000=pole_kontrola(X2000,Y2000)
print(Pk2000)

#Pole PL-1992
P1992=pole_gaus(X1992,Y1992)
print('Pole PL-1992 w m^2: '+str(P1992))
Pk1992=pole_kontrola(X1992,Y1992)
print(Pk1992,'\n')

#---------------------------------------------------------------------------------------------------------------------------
#PL-LAEA
Xl,Yl=przeliczanie_wspolrzednych(pun2000)
#print(Xl,Yl)
Pl=pole_gaus(Xl,Yl)
#print(Pl)
print('Pole PL-LAEA w m^2: '+str(Pl))
print(pole_kontrola(Xl,Yl))