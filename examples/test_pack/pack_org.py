import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy.spatial import KDTree
import random
from pyxtal import pyxtal as pyx

# シミュレーテッドアニーリングのパラメータ
INITIAL_TEMPERATURE = 1000.0
FINAL_TEMPERATURE = 1.0
ALPHA = 0.9  # 温度減衰率
MAX_ITER_PER_TEMP = 1000
RCUT = 2.0  # カットオフ半径（例）
ALPHA_PENALTY = 100.0  # ペナルティの強さ
LAMBDA = 1.0  # ラグランジュ乗数

# 1. CIFファイルの読み込み
cif_file = 'test.cif'  # CIFファイルのパスを指定
struc = pyx(molecular=True)
molecules = [ "CCC.smi" ]
struc.from_seed(cif_file, molecules=molecules)

# 重心位置の初期値（単位セル内に1つの分子と仮定）
initial_centroids = np.array([x.position for x in struc.mol_sites])

num_molecules = len(initial_centroids)

# 初期格子ベクトル（フラットな配列）
current_lattice = struc.lattice.get_matrix().flatten()

# 初期Euler角（無回転）
current_eulers = np.zeros(num_molecules * 3)

# 目的関数の定義
def generate_all_atoms(struc, lattice, centroids, eulers):
    lattice = lattice.reshape((3,3))
    all_atoms = []
    for idx, R_centroid in enumerate(centroids):
        # 分子の回転行列をEuler角から生成
        euler_angles = eulers[3*idx:3*idx+3]
        rotation = R.from_euler('xyz', euler_angles).as_matrix()
        for S in struc.rotation:
            transformed_molecule = R_centroid + S @ rotation @ molecule
            all_atoms.append(transformed_molecule)
    return np.vstack(all_atoms)

def calculate_density(struc, lattice, centroids):
    volume = abs(np.linalg.det(lattice.reshape((3,3))))
    num_total_molecules = len(centroids) * len(sym_ops)
    # 分子の体積をバウンディングボックスで近似
    molecule_array = np.array(molecule)
    V_molecule = np.prod(molecule_array.max(axis=0) - molecule_array.min(axis=0))
    return (num_total_molecules * V_molecule) / volume

def penalty_function(struc, lattice, centroids, eulers, rcut=RCUT, alpha=ALPHA_PENALTY):
    lattice = lattice.reshape((3,3))
    all_atoms = generate_all_atoms(lattice, centroids, eulers, sym_ops, molecule)
    penalty = 0.0
    tree = KDTree(all_atoms)
    pairs = tree.query_pairs(rcut)
    for (i, j) in pairs:
        distance = np.linalg.norm(all_atoms[i] - all_atoms[j])
        penalty += alpha * (rcut - distance)**2
    return penalty

def objective(struc, lattice, eulers, initial_centroids):
    density = calculate_density(struc, lattice, initial_centroids)
    penalty = penalty_function(struc, lattice, initial_centroids, eulers)
    return -density + LAMBDA * penalty

# シミュレーテッドアニーリングの実装
def simulated_annealing(struc, initial_lattice, initial_eulers, temperature, final_temp, alpha, max_iter, centroids):
    current_lattice = np.copy(initial_lattice)
    current_eulers = np.copy(initial_eulers)
    current_energy = objective(struc, current_lattice, current_eulers, centroids)
    
    best_lattice = np.copy(current_lattice)
    best_eulers = np.copy(current_eulers)
    best_energy = current_energy
    
    T = temperature
    
    while T > final_temp:
        for _ in range(max_iter):
            # 提案の生成
            # 格子ベクトルの微小変更
            lattice_move = current_lattice + np.random.normal(0, 0.01, size=current_lattice.shape)
            
            # 分子のEuler角の微小変更
            euler_move = current_eulers + np.random.normal(0, np.deg2rad(5), size=current_eulers.shape)
            
            # エネルギーの計算
            proposed_energy = objective(struc, lattice_move, euler_move, centroids)
            
            delta_E = proposed_energy - current_energy
            
            # 受容判定
            if delta_E < 0 or random.random() < np.exp(-delta_E / T):
                current_lattice = lattice_move
                current_eulers = euler_move
                current_energy = proposed_energy
                
                # ベスト状態の更新
                if current_energy < best_energy:
                    best_lattice = np.copy(current_lattice)
                    best_eulers = np.copy(current_eulers)
                    best_energy = current_energy
        
        # 温度の更新
        T *= alpha
        print(f"Temperature: {T:.2f}, Best Energy: {best_energy:.4f}")
    
    return best_lattice, best_eulers, best_energy

# 7. シミュレーテッドアニーリングの実行
best_lattice, best_eulers, best_energy = simulated_annealing(
    struc=struc,
    initial_lattice=current_lattice,
    initial_eulers=current_eulers,
    temperature=INITIAL_TEMPERATURE,
    final_temp=FINAL_TEMPERATURE,
    alpha=ALPHA,
    max_iter=MAX_ITER_PER_TEMP,
    centroids=initial_centroids,
)

# 結果の表示
print("\n最適な格子ベクトル:")
print(best_lattice.reshape((3,3)))
print("\n最適なEuler角（ラジアン）:")
print(best_eulers)

# 最適な回転行列の確認
optimal_rotations = []
for idx in range(len(initial_centroids)):
    euler = best_eulers[3*idx:3*idx+3]
    rotation = R.from_euler('xyz', euler).as_matrix()
    optimal_rotations.append(rotation)

print("\n最適な回転行列:")
for i, rot in enumerate(optimal_rotations):
    print(f"分子 {i+1}:")
    print(rot)
