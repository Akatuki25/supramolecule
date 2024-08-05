'''
本レポートでは超分子の検証の実装例として不十分と判断したものです
無視してください
'''
import numpy as np
import matplotlib.pyplot as plt

# パラメータの設定
num_particles = 100
box_size = 10.0
temperature = 1.0
num_steps = 1000
time_step = 0.01

# 初期位置と速度の設定
positions = np.random.rand(num_particles, 3) * box_size
velocities = np.random.randn(num_particles, 3) * np.sqrt(temperature)

# Lennard-Jonesポテンシャル
def lj_potential(r):
    return 4 * ((1/r)**12 - (1/r)**6)

# 力の計算
def compute_forces(positions):
    forces = np.zeros_like(positions)
    for i in range(num_particles):
        for j in range(i + 1, num_particles):
            r_vec = positions[i] - positions[j]
            r = np.linalg.norm(r_vec)
            if r < box_size / 2:
                f = -lj_potential(r) * r_vec / r
                forces[i] += f
                forces[j] -= f
    return forces

# シミュレーションの実行
for step in range(num_steps):
    forces = compute_forces(positions)
    velocities += forces * time_step
    positions += velocities * time_step
    positions = np.mod(positions, box_size)  # ボックスの端で折り返し

# 結果のプロット
plt.scatter(positions[:, 0], positions[:, 1])
plt.title("Particle Positions")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()
