!pip install pyscf
!pip install scipy

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pyscf import gto, scf

# 水分子の定義
def define_water_molecule(distance):
    mol = gto.Mole()
    mol.build(
        atom = f'''
        O 0 0 0;
        H 0 0 {0.958};
        H 0 0 {-0.958};
        ''',
        basis = 'cc-pVDZ'  # 基底セット
    )
    return mol

# 分子のジオメトリー確認
def check_molecule_geometry():
    try:
        mol = define_water_molecule(distance=0)  # 距離を0に設定して確認
        print("Molecule Geometry:")
        print("Atoms:", mol.atom_coords())
        print("Basis set:", mol.basis)
    except Exception as e:
        print(f"Error: {e}")

check_molecule_geometry()

# 超分子の定義
def define_supermolecule(distance):
    mol = gto.Mole()
    mol.build(
        atom = f'''
        O 0 1 0;
        H 0.958 1 0;
        H -0.958 1 0;
        O 0 0 {distance};
        H 0 0 {distance + 0.958};
        H 0 0 {distance - 0.958}
        ''',
        basis = 'cc-pVDZ'  # 基底セット
    )
    return mol

# 分子のジオメトリー確認
def check_supermolecule_geometry(distance):
    try:
        mol = define_supermolecule(distance)
        print("Supermolecule Geometry:")
        print("Atoms:", mol.atom_coords())
        print("Basis set:", mol.basis)
    except Exception as e:
        print(f"Error: {e}")

# 距離が4.0 Åの場合
check_supermolecule_geometry(distance=4.0)

# 超分子の定義
def define_supermolecule(distance):
    mol = gto.Mole()
    mol.build(
        atom = f'''
        O 0 1 0;
        H 0.958 1 0;
        H -0.958 1 0;
        O 0 0 {distance};
        H 0 0 {distance + 0.958};
        H 0 0 {distance - 0.958}
        ''',
        basis = 'cc-pVDZ'  # 基底セット
    )
    return mol

# 分子のジオメトリー確認
def check_supermolecule_geometry(distance):
    try:
        mol = define_supermolecule(distance)
        print("Supermolecule Geometry:")
        print("Atoms:", mol.atom_coords())
        print("Basis set:", mol.basis)
    except Exception as e:
        print(f"Error: {e}")

# 現実的な距離で試行
check_supermolecule_geometry(distance=1.0)  # 1.0 Å

# エネルギー最適化の実行と座標の取得
def optimize_supermolecule(distance):
    mol = define_supermolecule(distance)
    mf = scf.RHF(mol)
    mf.kernel()  # SCF計算を実行
    optimized_coords = mol.atom_coords()
    return optimized_coords

# エネルギー最適化の実行と座標表示
optimized_coords = optimize_supermolecule(distance=3.779)
print("Optimized Coordinates:")
for i, coord in enumerate(optimized_coords):
    print(f"Atom {i+1}: {coord}")

# 原子の種類を取得
atoms = ['O', 'H', 'H', 'O', 'H', 'H']

# 最適化された座標のプロット
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 各原子をプロット
for i, (atom, coord) in enumerate(zip(atoms, optimized_coords)):
    if atom == 'O':
        ax.scatter(coord[0], coord[1], coord[2], color='red', label='O' if i == 0 else "")
    elif atom == 'H':
        ax.scatter(coord[0], coord[1], coord[2], color='blue', label='H' if i == 1 else "")

# ラベルの設定
ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')
ax.set_zlabel('Z Coordinate')
ax.legend()

# プロット表示
plt.show()
