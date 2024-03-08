#ライブラリのインポート 
import numpy as np
import scipy as sp

from pyevtk.hl import pointsToVTK

#定数の定義
DT=0.001
FLUID_DENSITY = 1000
PARTICLE_DISTANCE = 0.025
COMPRESSIBILITY = 0.45*10**-9
N0_forNumberDensity = 0.0
N0_forGradient = 0.0
N0_forLaplacian = 0.0
RE_FOR_NUMBER_DENSITY = 2.1*PARTICLE_DISTANCE
RE_FOR_GRADIENT = 2.1*PARTICLE_DISTANCE
RE_FOR_LAPLACIAN = 3.1*PARTICLE_DISTANCE
Lambda = 0.0
VISCOSITY = 1.0*10**-6
G=9.8
DIMENSION = 2
GAMMA = 0.2
Pressure = []
particles = []
#クラスの定義
class particle:
    def __init__(self,position,velocity,type):
        self.position = position
        self.velocity = velocity
        self.acceleration = np.array([0,0])
        self.type = type #0:液体 1:壁 2:ダミー
        self.NumberDensity = 0
        self.BoundaryCondition = "INNER_PARTICLE"

def init_particles():
    for iX in range(41):
        for iY in range(41):
            #壁の設置
            if((0<=iX<=1 or 39<=iX<=40) and 2<=iY<=21):
                particles.append(particle(np.array([iX*PARTICLE_DISTANCE,iY*PARTICLE_DISTANCE]),np.array([0.0,0.0]),2))
            if((2<=iX<=3 or 37<=iX<=38) and 4<=iY<=21):
                particles.append(particle(np.array([iX*PARTICLE_DISTANCE,iY*PARTICLE_DISTANCE]),np.array([0.0,0.0]),1))
            if(0<=iX<=40 and 0<=iY<=1):
                particles.append(particle(np.array([iX*PARTICLE_DISTANCE,iY*PARTICLE_DISTANCE]),np.array([0.0,0.0]),2))
            if(2<=iX<=38 and 2<=iY<=3):
                particles.append(particle(np.array([iX*PARTICLE_DISTANCE,iY*PARTICLE_DISTANCE]),np.array([0.0,0.0]),1))
            #流体の設置
            if(10<=iX<=21 and 4<=iY<=24):
                particles.append(particle(np.array([iX*PARTICLE_DISTANCE,iY*PARTICLE_DISTANCE]),np.array([0.0,0.0]),0))

#LamdaとN0を求める
def calNZeroAndLamda():
    global N0_forNumberDensity
    global N0_forGradient
    global N0_forLaplacian
    global Lambda
    for iX in range(-4,5):
        for iY in range(-4,5):
            if (iX == 0 and iY == 0):
                continue
            squaredDistance = (iX*PARTICLE_DISTANCE)**2+(iY*PARTICLE_DISTANCE)**2
            
            distance = np.sqrt(squaredDistance)
            N0_forNumberDensity = N0_forNumberDensity + w(distance,RE_FOR_NUMBER_DENSITY)
            N0_forGradient = N0_forGradient + w(distance,RE_FOR_GRADIENT)
            N0_forLaplacian = N0_forLaplacian + w(distance,RE_FOR_LAPLACIAN)
            Lambda += squaredDistance*w(distance,RE_FOR_LAPLACIAN)
    Lambda = Lambda/N0_forLaplacian

#外力項を計算
def calGravity():
    for i in particles:
        if i.type == 0:
            i.acceleration = np.array([0,-G])

#粘性項を計算
def calViscosity():
    a = (VISCOSITY*2*DIMENSION)/N0_forLaplacian*Lambda
    for i in particles:
        if i.type == 0:
            for j in particles:
                Laplacian_V = 0
                if i == j:
                    continue
                elif j.type == 2:
                    continue
                else:
                    SquaredDistance = (i.position[0]-j.position[0])**2+(i.position[1]-j.position[1])**2
                    distance = np.sqrt(SquaredDistance)
                    if distance==0:
                        print("ゼロで割る！")
                    if distance<RE_FOR_GRADIENT:
                        Laplacian_V += w(distance,RE_FOR_LAPLACIAN)*(j.velocity[0]-i.velocity[0])
            i.acceleration += a*Laplacian_V
        else:
            continue

#粒子の更新
def moveParticles():
    for i in particles:
        if i.type == 0:
            i.velocity += i.acceleration*DT
            i.position += i.velocity*DT
            i.acceleration = np.array([0.0,0.0])
        else:
            continue

#圧力項の計算
def calPressure():
    A = np.empty((0))
    b = np.empty((0))

    #粒子密度の計算
    for i in particles:
        for j in particles:
            if i == j:
                continue
            else:
                SquaredDistance = (i.position[0]-j.position[0])**2+(i.position[1]-j.position[1])**2
                distance = np.sqrt(SquaredDistance)
                if distance==0:
                        print("ゼロで割る！")
                i.NumberDensity += w(distance,RE_FOR_NUMBER_DENSITY)
        i.BoundaryCondition = "INNER_PARTICLE"
        #ディリクレ境界条件の設定
        if i.type == 2:
            i.BoundaryCondition = "DUMMY"
        elif i.NumberDensity/N0_forNumberDensity < 0.97:
            i.BoundaryCondition = "SURFACE_PARTICLE"
    
    #bの設定
    for i in particles:
        if i.BoundaryCondition == "INNER_PARTICLE":
            b = np.append(b,GAMMA*(1/DT**2)*((i.NumberDensity-N0_forNumberDensity)/N0_forNumberDensity))
        if i.BoundaryCondition == "DUMMY":
            b = np.append(b,0)
        if i.BoundaryCondition == "SURFACE_PARTICLE":
            b = np.append(b,0)
        

    #Aの設定
    for index,i in enumerate(particles):
        a_i = np.empty([0])
        a_i_i = 0
        a = 2*DIMENSION/(N0_forLaplacian*Lambda)
        for jedex,j in enumerate(particles):
            if i.BoundaryCondition == "INNER_PARTICLE":
                if i == j:
                    a_i = np.append(a_i,0)
                else:
                    SquaredDistance = (i.position[0]-j.position[0])**2+(i.position[1]-j.position[1])**2
                    distance = np.sqrt(SquaredDistance)
                    if distance==0:
                        print("ゼロで割る！1",i.position,j.position)
                    a_i = np.append(a_i,-a*w(distance,RE_FOR_LAPLACIAN)/FLUID_DENSITY)
                    a_i[jedex] += a*w(distance,RE_FOR_LAPLACIAN)/FLUID_DENSITY
                    a_i_i += a*w(distance,RE_FOR_LAPLACIAN)/FLUID_DENSITY
            elif i.BoundaryCondition == "DUMMY":
                a_i = np.append(a_i,0)
            elif i.BoundaryCondition == "SURFACE_PARTICLE":
                a_i = np.append(a_i,0)
        a_i_i += COMPRESSIBILITY/(DT**2)
        a_i[index]=a_i_i
        A = np.append(A,a_i)
    A = np.reshape(A,(len(particles),len(particles)))
    #Ax=bをガウス掃き出し法を用いて解く
    x=sp.sparse.linalg.cg(A,b)
    Pressure = x[0]
    for index,i in enumerate(particles):
        if i.type==2:
            Pressure[index]=0
        if i.BoundaryCondition == "SURFACE_PARTICLE":
            Pressure[index]=0
        if Pressure[index]<0:
            Pressure[index]=0
    for index,i in enumerate(particles):
        i.NumberDensity = 0
        Gradient = np.array([0.0,0.0])
        if i.type==0:
            minPressure = None
            for jedex,j in enumerate(particles):
                SquaredDistance = (i.position[0]-j.position[0])**2+(i.position[1]-j.position[1])**2
                distance = np.sqrt(SquaredDistance)
                if minPressure == None or minPressure >= Pressure[jedex]:
                    if distance <=RE_FOR_GRADIENT:
                        minPressure = Pressure[jedex]
            
            for jedex,j in enumerate(particles):
                if i == j:
                    continue
                else:
                    SquaredDistance = (j.position[0]-i.position[0])**2+(j.position[1]-i.position[1])**2
                    distance = np.sqrt(SquaredDistance)
                    if distance<=RE_FOR_GRADIENT:
                        Gradient += ((Pressure[jedex]-minPressure)/SquaredDistance)*(j.position-i.position)*w(distance,RE_FOR_GRADIENT)
            i.acceleration = -1*(Gradient/FLUID_DENSITY) * (DIMENSION/N0_forGradient)
    for i in particles:
        i.velocity += i.acceleration*DT
        i.position += i.acceleration*DT*DT
            



def save(t):
    x=[]
    y=[]
    z=[]
    for index,i in enumerate(particles):
        x.append(i.position[0])
        y.append(i.position[1])
        z.append(0)
    pointsToVTK("./"+str(int(t/20)),np.array(x),np.array(y),np.array(z))
    


#重み関数の定義
def w(r,re):
    if(r<re):
        return (re/r)-1
    else:
        return 0


init_particles()   
calNZeroAndLamda()
print(N0_forNumberDensity,N0_forGradient,N0_forLaplacian,Lambda)
for t in range(100):
    calGravity()
    calViscosity()
    moveParticles()
    
    calPressure()
    if t%20==0:
        print(t)
        save(t)