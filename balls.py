import numpy as np 
import sympy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches

particles = []

#定数
dt = 1/1000
g = 9.80665
ground_height = 1
ball_radius = 0.1
coords = []
number_of_balls = 1000
def check_hit(x,y):
    if(ball_radius>=abs(y-ground_height)):
        return True
    else:
        return False

def main():
    #ボールの生成
    for i in range(number_of_balls):
        particles.append(Particle(np.array([0,10]),9*((1 - -1) * np.random.rand(2) + -1),1))
    for i in range(5000):
        print("t:"+str(i))
        coord = []
        for p in range(number_of_balls):
            particles[p].update(np.array([0,-g]),dt)
            coord.append(particles[p].position)
        coords.append(coord)

    
    # フィギュアとサブプロットを作成
    fig, ax = plt.subplots()
    # サブプロットの軸を調整
    ax.set_xlim(-10, 10)
    ax.set_ylim(0, 20)
    ax.plot([-50,50],[ground_height,ground_height],color="sienna")
    # ボールのオブジェクトを作成
    ball = ax.scatter([], [])

    # アニメーションの更新関数
    def update(i):
        x = []
        y = []
        for p in range(number_of_balls):
            x.append(coords[i][p][0])
            y.append(coords[i][p][1])
        # ボールの座標を更新
        ball.set_offsets(list(zip(x, y)))
        # ボールを返す
        return ball,

    # アニメーションのオブジェクトを作成
    ani = animation.FuncAnimation(fig, update, frames=len(coords), interval=1000/(1/dt), blit=True)
    # アニメーションを表示
    ani.save('anim.mp4',writer="ffmpeg")

#ボール
class Particle:
    def __init__(self, p, v ,m):
        self.position = p
        self.velocity = v
        self.m = m
    def update(self,a,step_size):
        #ここで接地判定も行う。
        if(check_hit(self.position[0],self.position[1])==True):
            self.velocity[1]=self.velocity[1]*-1
        #速度,位置の更新
        self.velocity = self.velocity + a*step_size
        self.position = self.position + self.velocity*step_size
if __name__ == "__main__":
    main()
