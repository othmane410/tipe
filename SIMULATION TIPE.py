## importations
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as resol
from scipy.integrate import odeint
import random
import scipy.fft as fft
from math import pi


## COMMENTAIRE IMPORTANT (POUR OTHMANE NOUMANE)
# Les programmes suivants sont quasiment tous des "tests" pour voir où mes erreurs ont lieu
# Seul le dernier programme RESOLUTION2 est le plus "abouti"



## CONSTANTES
g = 9.81
f = 10 #calculé pour une bande et gardé pour tout le TIPE
Fseuil = 5000
## fonctions - programmes

def phiSysteme2(u,t):
    (y,v) = u
    dy = v
    dv = -np.tan(t)*v + np.cos(t)**2*y
    return np.array((dy,dv))

def euler(phi,y0,T):
    Y = [0]*len(T)
    Y[0] = y0
    for i in range(len(T)-1):
        pas = T[i+1]- T[i]
        pente = phi(Y[i],T[i])
        Y[i+1] = Y[i] + pas*pente
    return Y

# def TRD(position_0,k,m,w,delta_t):
#     # en entrée : position_0 la coordonnée x de la masse à t = 0
#     #             k la constante de raideur du ressort
#     #             m la masse de l'objet
#     #             w la vitesse de rotation du tapis
#     #             delta_t l'espacement temporel entre la valeur t_a et t_b
#     #
#     # en sortie : la représentation graphique de la position de la masse ainsi que de sa vitesse
#     v = a*w
#     x1 = v
#     Tregime1 = [0]
#     pas = delta_t
#     while delta_t < f*m*g/(k*v):
#         Tregime1.append(delta_t)
#         delta_t += delta_t
#     t1 = delta_t
#     position_1 = f*m*g/k
#     plt.plot(Tregime1,[x1*t + position_0 for t in Tregime1], label = "position régime collé")
#     plt.plot(Tregime1, [v for t in Tregime1], label = "vitesse en régime collé")
#     def phiSysteme(x0,t):
#         (y,v) = x0
#         dy = v
#         dv = f*g - k*x0/m
#         return np.array((dy,dv))
#     plt.plot(Tregime1,np.array(euler(phiSysteme,position_1,Tregime1))[:,0],label = "position régime glissé")
#     plt.legend()
#     plt.show()

def TRD(position_0,k,m,w,delta_t):
    # en entrée : position_0 la coordonnée x de la masse à t = 0
    #             k la constante de raideur du ressort
    #             m la masse de l'objet
    #             w la vitesse de rotation du tapis
    #             delta_t l'espacement temporel entre la valeur t_a et t_b
    #
    # en sortie : la représentation graphique de la position de la masse ainsi que de sa vitesse
    v = a*w
    x1 = v
    Tregime1 = [0]
    Tregime2 = [0]
    pas = delta_t
    while delta_t < f*m*g/(k*v):
        Tregime1.append(delta_t)
        delta_t += delta_t
    t1 = delta_t
    position_1 = f*m*g/k
    w0 = (k/m)**0.5
    plt.plot(Tregime1,[x1*t + position_0 for t in Tregime1], label = "position régime collé")
    plt.plot(Tregime1, [v for t in Tregime1], label = "vitesse en régime collé")
    while pas< 5 :
        Tregime2.append(pas)
        pas  += pas
    def x(t):
        return v/w0*np.sin(w0*(t-t1))+f*m*g/k
    def vitesse(t):
        return v*np.cos(w0*(t-t1))
    plt.plot(Tregime2,[x(t) for t in Tregime2],label = "position régime glissé")
    plt.plot(Tregime2,[vitesse(t) for t in Tregime2],label = "vitesse régime glissé")
    plt.legend()
    plt.show()

def SIMULATION(x0,k,m,w,delt):
    # en entrée : position_0 la coordonnée x de la masse à t = 0
    #             k la constante de raideur du ressort
    #             m la masse de l'objet
    #             w la vitesse de rotation du tapis
    #             delta_t l'espacement temporel entre la valeur t_a et t_b
    #
    # en sortie : la représentation graphique de la position de la masse ainsi que de sa vitesse
    return None

def equation(trajectoire,t):
    [x,dx] = trajectoire
    return (dx,SorSL(x,dx,t))

def equation2(trajectoire,t):
    [x,dx] = trajectoire
    return (dx,SorSL2(x,dx,t))

def SorSL(x,dx,t):
    if k*a*w*t-f*m*g < 0 :
        return 0
    else:
        return f*g-t*-x/m

def SorSL2(x,dx,t):
    return f*g-t*-x/m-frottement(x,dx,t)/m

def frottement(x,dx,t):
    if k*x < Fseuil:
        return 0
    else:
        return eta*dx**2


def RESOLUTION(CIx,CIv,m0,k,f,w,a,n):
    # Initialisation : CIx, CIv conditions initiales de position et de vitesse, m0 sa masse initiale
    t0 = 0
    tmax = 0
    def EQUA(trajectoire,t):
        x,dx = trajectoire
        return (dx,ddx(x,dx,t))
    def RESSORT(x):
        return k*x
    def TANG(t):
        return k*a*w*(t-t0)
    def TAPIS(t):
        return a*w
    def BRUIT(x,t):
        return 40*np.sin(x**5*t)
    def ddx(x,dx,t):
        if TANG(t)-f*m0*g < 0:
            return 0
        else: #régime glissé
            return f*g-RESSORT(x)/m0
    while TANG(tmax) < f*m0*g:
        tmax += 0.001
    for i in range(n):
        Temps = np.linspace(t0,tmax,300)
        X = odeint(EQUA,[CIx,CIv],Temps)
        plt.plot(Temps,X[:,0],color = 'm')
        plt.plot(Temps,X[:,1],color = 'k')
        CIx = X[:,0][-1]
        CIv = X[:,1][-1]
        t0 = tmax
        tmax = X[-1,:][-1]
    plt.show()


def RESOLUTION2(CIx,CIv,m0,k,f,w,a,tmax):
    ''' Initialisation : CIx, 
                        CIv conditions initiales de position et de vitesse,
                        m0 sa masse initiale, k constante du ressort, 
                        f coefficient de frottement, w vitesse angulaire, 
                        a périmètre de la bande, 
                        tmax le temps de la simulation
        Sortie : graphique de la vitesse et de la position du solide ( VERT pour VITESSE, VIOLET pour position) ainsi que leur FFT respectives
    '''
    def EQUA(trajectoire,t):
        x,dx = trajectoire
        return (dx,ddx(x,dx,t))

    def instabilite(x,t):
        if np.sin(t)< np.cos(x):
            return 50*np.cos(t)
        else:
            return -50*np.sin(t)
    def RESSORT(x):
        return k*x
    def TANG(t):
        t0=0
        while k*a*w*(t-t0)  < f*m0*g :
            return 0
        t0 += t
        return f*m0*g
    def ddx(x,dx,t):
        if TANG(t) < f*m0*g :
            return 0
        else:#régime glissé
            return f*g - RESSORT(x)/m0 - instabilite(x,t)
    N = 1000
    Temps = np.linspace(0,tmax,N)
    X = odeint(EQUA,[CIx,CIv],Temps)

    fig, axs = plt.subplots(3)
    axs[0].set_title('cinématique')
    axs[0].plot(Temps,X[:,0],color = 'm') #position
    axs[0].plot(Temps,X[:,1],color = 'g') #vitesse

    yf = fft.fft(X[:,0])
    yplot = fft.fftshift(yf)

    vf = fft.fft(X[:,1])
    vplot = fft.fftshift(vf)

    xf = fft.fftfreq(N,1/(0.5*N))
    xf = fft.fftshift(xf)

    axs[1].set_title('FFT position')
    axs[1].plot(xf,(1/N)*np.abs(yplot),color='m')
    axs[1].set_xlim(0, tmax)

    axs[2].set_title('FFT vitesse')
    axs[2].plot(xf,(1/N)*np.abs(vplot),color = 'g')
    axs[2].set_xlim(0, tmax)
    plt.show()









## traitement
