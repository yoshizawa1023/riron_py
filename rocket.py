import numpy as np 
import sympy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches

particles = []

#定数
dt = 1/1000
G=6.6743*10**-11
earth_radius = 1
earth_mass = 1500000000000
rocket_mass = 1
coords = []
number_of_rockets = 1
def check_hit(potision):
    if(earth_radius>np.linalg.norm(potision,ord=2)):
        return True
    else:
        return False

def main():
    #ボールの生成
    for i in range(number_of_rockets):
        particles.append(Particle(np.array([0,1.1]),np.array([35,20]),rocket_mass))
    for i in range(3000):
        print("t:"+str(i))
        coord = []
        for p in range(number_of_rockets):
            a = particles[p].position*(-G*earth_mass)/np.linalg.norm(particles[p].position,ord=2)**3
            particles[p].update(a,dt)
            coord.append(particles[p].position)
        coords.append(coord)

    
    # フィギュアとサブプロットを作成
    fig, ax = plt.subplots()
    # サブプロットの軸を調整
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    ax.set_aspect('equal')
    c = patches.Circle(xy=(0,0),radius=earth_radius,fc='k',fill=False)
    ax.add_patch(c)
    # ボールのオブジェクトを作成
    ball = ax.scatter([], [])

    # アニメーションの更新関数
    def update(i):
        x = []
        y = []
        for p in range(number_of_rockets):
            x.append(coords[i][p][0])
            y.append(coords[i][p][1])
        # ボールの座標を更新
        ball.set_offsets(list(zip(x, y)))
        # ボールを返す
        return ball,

    # アニメーションのオブジェクトを作成
    ani = animation.FuncAnimation(fig, update, frames=len(coords), interval=1000/(1/dt), blit=True)
    # アニメーションを表示
    #ani.save('rocket_long.mp4',writer='ffmpeg')
    plt.show()

#ボール
class Particle:
    def __init__(self, p, v ,m):
        self.position = p
        self.velocity = v
        self.m = m
    def update(self,a,step_size):
        #ここで接地判定も行う。
        if(check_hit(self.position)==False):
            #速度,位置の更新
            self.velocity = self.velocity + a*step_size
            self.position = self.position + self.velocity*step_size

if __name__ == "__main__":
    main()
