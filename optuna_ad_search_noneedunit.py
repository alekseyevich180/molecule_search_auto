import csv
import os

import numpy as np
import optuna
import pandas as pd
import paramiko
from dotenv import load_dotenv


class Share:
    def __init__(
        self,
        WORK_DIR,
        JOB_USER,
        INIT_NUM,
        NUM_CANDIDATES,
        MOLE_FILE,
        SUR_FILE,
        MASSDICT,
        HIGHT,
    ):
        self.WORK_DIR = WORK_DIR
        self.JOB_USER = JOB_USER
        self.INIT_NUM = int(INIT_NUM)
        self.NUM_CANDIDATES = int(NUM_CANDIDATES)
        self.MOLE_FILE = MOLE_FILE
        self.SUR_FILE = SUR_FILE
        self.MASSDICT = MASSDICT
        self.HIGHT = float(HIGHT)

    def write_vasp(self, filename, cell_params, atom_types, num_atoms, coordinates):
        with open(filename, "w") as f:
            f.write("Generated POSCAR\n")
            f.write("1.0\n")
            for vec in cell_params:
                f.write(f"{vec[0]:.6f} {vec[1]:.6f} {vec[2]:.6f}\n")
            f.write(" ".join(atom_types) + "\n")
            f.write(" ".join(map(str, num_atoms)) + "\n")
            f.write("Direct\n")

            for coord in coordinates:
                f.write(f"{coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")

    def euler_rotate(self, mole_coordinates, phi, theta, psi):
        # Euler角に基づく回転行列を計算
        phi_rad = np.radians(phi)
        theta_rad = np.radians(theta)
        psi_rad = np.radians(psi)

        # Z軸回りの回転行列 (phi)
        Rz_phi = np.array(
            [
                [np.cos(phi_rad), np.sin(phi_rad), 0],
                [-np.sin(phi_rad), np.cos(phi_rad), 0],
                [0, 0, 1],
            ]
        )

        # X'軸回りの回転行列 (theta)
        Rx_theta = np.array(
            [
                [1, 0, 0],
                [0, np.cos(theta_rad), np.sin(theta_rad)],
                [0, -np.sin(theta_rad), np.cos(theta_rad)],
            ]
        )

        # Z''軸回りの回転行列 (psi)
        Rz_psi = np.array(
            [
                [np.cos(psi_rad), np.sin(psi_rad), 0],
                [-np.sin(psi_rad), np.cos(psi_rad), 0],
                [0, 0, 1],
            ]
        )

        # 最終的な回転行列(ZXZ座標)
        R = Rz_phi @ Rx_theta @ Rz_psi

        # 原子の座標を回転
        new_positions = np.dot(mole_coordinates, R.T)  # 転置を取る

        # 回転後の座標を更新
        mole_coordinates = new_positions
        return mole_coordinates

    def make_ad(self, work_type, trial, ID, phi, theta, psi, x, y, filename):
        row_data = [trial, ID, x, y, phi, theta, psi, ""]
        with open(filename, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(row_data)

        dir_path = f"{self.WORK_DIR}/{work_type}/{trial}_{ID}"
        if not os.path.exists(dir_path):
            os.makedirs(dir_path, exist_ok=True)

        poscar_path = f"{dir_path}/POSCAR"

        # VASPファイルを読み込み
        with open(self.MOLE_FILE, "r") as f:
            lines = f.readlines()
        mole_cell = []
        for i in range(2, 5):
            cell_list = lines[i].split()
            mole_cell.append([float(c) for c in cell_list[:3]])
        mole_cell = np.array(mole_cell)

        # "Direct"という行を探す
        direct_index = -1
        for i, line in enumerate(lines):
            if "Direct" in line:
                direct_index = i
                break

        if direct_index == -1:
            raise ValueError("Direct is not found in the mole file.")

        # 原子種と数
        atoms_info_line = lines[5].split()
        num_atoms = list(map(int, lines[6].split()))

        # 分子座標 (Direct座標)
        mole_coordinates = []
        atom_types = []  # 原子種のリスト

        atom_counter = 0  # 原子種のインデックスを追跡
        atom_index = 0  # 現在処理している原子種のカウンター

        for i in range(direct_index + 1, direct_index + 1 + sum(num_atoms)):
            coords = lines[i].split()
            # mole_coordinates.append([float(c) for c in coords[:3]])  # 相対座標
            # mole_cellを用いて絶対座標に変換してから、分子の座標を取得
            mole_coordinates.append(
                np.array([float(c) for c in coords[:3]]) @ mole_cell
            )

            # atom_index に基づいて、適切な原子種を選択
            atom_types.append(atoms_info_line[atom_counter])

            # atom_index を更新
            atom_index += 1

            # 次の原子種に進む条件
            if atom_index >= num_atoms[atom_counter]:
                atom_counter += 1  # 次の原子種に切り替え
                atom_index = 0  # 新しい原子種のカウントをリセット

        # 重心の計算
        total_mass = 0
        weighted_sum_x = 0
        weighted_sum_y = 0
        weighted_sum_z = 0

        for i, coords in enumerate(mole_coordinates):
            atom_type = atom_types[i]
            mass = self.MASSDICT.get(atom_type, 0)  # 原子種に対応する質量を取得
            total_mass += mass
            weighted_sum_x += mass * coords[0]
            weighted_sum_y += mass * coords[1]
            weighted_sum_z += mass * coords[2]

        # 重心座標
        center_of_mass_x = weighted_sum_x / total_mass
        center_of_mass_y = weighted_sum_y / total_mass
        center_of_mass_z = weighted_sum_z / total_mass
        # min_z = min([coords[2] for coords in mole_coordinates])
        mole_coordinates = np.array(mole_coordinates)
        mole_coordinates -= np.array(
            [center_of_mass_x, center_of_mass_y, center_of_mass_z]
        )
        rotate_coordinates = self.euler_rotate(mole_coordinates, phi, theta, psi)
        # rotate_coordinates += np.array([x, y, 0] @ mole_cell)

        with open(self.SUR_FILE, "r") as f:
            lines = f.readlines()
        cell_params = []
        for i in range(2, 5):
            cell_list = lines[i].split()
            cell_params.append([float(c) for c in cell_list[:3]])
        cell_params = np.array(cell_params)
        rotate_coordinates += np.array([x, y, 0] @ cell_params)
        
        sur_atoms = list(map(int, lines[6].split()))

        sur_coordinates = []
        # "Direct"という行を探す
        sur_direct = -1
        for i, line in enumerate(lines):
            if "Direct" in line:
                sur_direct = i
                break
        if sur_direct == -1:
            raise ValueError("Direct is not found in the sur file.")

        for i in range(sur_direct + 1, sur_direct + 1 + sum(sur_atoms)):
            coords = lines[i].split()
            sur_coordinates.append([float(c) for c in coords[:3]])

        # **1. 相対座標を絶対座標に変換**
        # rotate_coordinates = rotate_coordinates @ cell_params
        absolute_sur_coordinates = sur_coordinates @ cell_params

        # **2. mole の最低 z 座標 (絶対座標) を取得**
        min_mole_z = np.min(rotate_coordinates[:, 2])

        # **3. sur の最大 z 座標 (絶対座標) を取得**
        max_sur_z = np.max(absolute_sur_coordinates[:, 2])

        # **4. mole の z 座標をシフト (最低 z を max_sur_z + self.HIGHT に)**
        z_shift_init = (max_sur_z + self.HIGHT) - min_mole_z
        rotate_coordinates[:, 2] += z_shift_init

        # **5. ここで最短距離を計算(少しずつzを減らす)**

        step = 0.001  # z方向に下げるステップサイズ

        while True:
            # **6. 周期境界条件を考慮して、分子の座標を計算**
            shifts = [-1, 0, 1]
            pbc_sur_coords = []

            # a方向 = cell_params[0], b方向 = cell_params[1]
            for dx in shifts:
                for dy in shifts:
                    shift_vector = dx * cell_params[0] + dy * cell_params[1]
                    shifted_coords = absolute_sur_coordinates + shift_vector
                    pbc_sur_coords.append(shifted_coords)

            # (9 * N_surf_atoms, 3) の配列にする
            pbc_sur_coords = np.vstack(pbc_sur_coords)

            # 最短距離を求める
            while True:
                distance_matrix = np.linalg.norm(
                    rotate_coordinates[:, None, :] - pbc_sur_coords[None, :, :],
                    axis=2,
                )
                min_distance = np.min(distance_matrix)

                if min_distance < self.HIGHT:
                    break

                rotate_coordinates[:, 2] -= step

            # 最短距離が HIGHT を下回ったら終了
            if min_distance < self.HIGHT:
                break

            # 分子全体を z 方向に下げる
            rotate_coordinates[:, 2] -= step

        # **7. 相対座標 (Direct座標) に戻す**
        relative_mole_coordinates = rotate_coordinates @ np.linalg.inv(cell_params)

        sur_atom_types = lines[5].split()
        sur_num_atoms = list(map(int, lines[6].split()))

        # mole の原子種と数 (もとのデータを使う)
        mole_atom_types = atoms_info_line
        mole_num_atoms = num_atoms

        # 原子種と原子数を統合
        all_atom_types = sur_atom_types + mole_atom_types
        all_num_atoms = sur_num_atoms + mole_num_atoms

        # 分子座標を統合
        combined_coordinates = np.vstack([sur_coordinates, relative_mole_coordinates])
        self.write_vasp(
            poscar_path,
            cell_params,
            all_atom_types,
            all_num_atoms,
            combined_coordinates,
        )

    def bayes_opt(self):
        search_space = {
            "x": optuna.distributions.FloatDistribution(0, 1),
            "y": optuna.distributions.FloatDistribution(0, 1),
            "phi": optuna.distributions.FloatDistribution(-180, 180),
            "theta": optuna.distributions.FloatDistribution(-180, 180),
            "psi": optuna.distributions.FloatDistribution(-180, 180),
        }

        # データ読み込み
        df = pd.read_csv("data.csv")
        valid_rows = df[df["energy"] < 0]

        sampler = optuna.samplers.TPESampler(
            n_startup_trials=0, n_ei_candidates=1000, constant_liar=True
        )
        study = optuna.create_study(direction="minimize", sampler=sampler)

        for _, row in valid_rows.iterrows():
            params = {
                "x": row["x"],
                "y": row["y"],
                "phi": row["phi"],
                "theta": row["theta"],
                "psi": row["psi"],
            }
            study.add_trial(
                optuna.trial.create_trial(
                    params=params,
                    distributions=search_space,
                    values=[row["energy"]],
                )
            )

        recommended_trials = []
        for _ in range(self.NUM_CANDIDATES):
            trial = study.ask()
            new_params = {
                "x": trial.suggest_float("x", 0, 1),
                "y": trial.suggest_float("y", 0, 1),
                "phi": trial.suggest_float("phi", -180, 180),
                "theta": trial.suggest_float("theta", -180, 180),
                "psi": trial.suggest_float("psi", -180, 180),
            }
            recommended_trials.append(new_params)

        return recommended_trials

    def ssh_connect(self, HOSTNAME, USERNAME, pkey, work_type, trial, ID):
        with paramiko.SSHClient() as client:
            client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            client.connect(hostname=HOSTNAME, port=22, username=USERNAME, pkey=pkey)

            work_dir = f"{self.WORK_DIR}/{work_type}/{trial}_{ID}"
            job_script = f"{self.JOB_USER}_{trial}_{ID}.sh"

            command = f'''
            check_dir() {{
                if [ ! -f "{work_dir}/POSCAR" ]; then
                    echo "Error: POSCAR does not exist" >&2
                    exit 1
                fi
            }}
            check_dir &&
            cd {self.WORK_DIR} &&
            cp INCAR KPOINTS POTCAR .env {work_dir}/ &&
            cp run_vasp6.4.3.sh {work_dir}/{job_script} &&
            cd {work_dir} &&
            qsub {job_script}
            '''

            stdin, stdout, stderr = client.exec_command(command)
            ssh_results = [stdout.read().decode("utf-8"), stderr.read().decode("utf-8")]
            stdin.close()
        return ssh_results


# job_time.txtが無ければstart時間を書き込む
if not os.path.exists("job_time.txt"):
    with open("job_time.txt", "a") as f:
        f.write(f"start: {pd.Timestamp.now()}\n")

# SSH接続情報の設定
HOSTNAME = "genkai"
USERNAME = "ku40000345"
KEY_FILENAME = "/home/pj24001724/ku40000345/rsa_key/busseiken_private"
pkey = paramiko.RSAKey(filename=KEY_FILENAME)

MASSDICT = {"H": 1, "C": 12, "N": 14, "O": 16}
load_dotenv()
WORK_DIR = os.getenv("WORK_DIR")
JOB_USER = os.getenv("JOB_USER")
INIT_NUM = os.getenv("INIT_NUM")
NUM_CANDIDATES = os.getenv("NUM_CANDIDATES")
SEARCH_NUM = os.getenv("SEARCH_NUM")
MOLE_FILE = os.getenv("MOLE_FILE")
SUR_FILE = os.getenv("SUR_FILE")
HIGHT = os.getenv("HIGHT")

filename = "data.csv"
header = ["trial", "ID", "x", "y", "phi", "theta", "psi", "energy"]

total_rows = 0
if os.path.exists(filename):
    with open(filename, "r") as file:
        total_rows = len(list(csv.reader(file))) - 1
else:
    with open(filename, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(header)

share_instance = Share(
    WORK_DIR, JOB_USER, INIT_NUM, NUM_CANDIDATES, MOLE_FILE, SUR_FILE, MASSDICT, HIGHT
)
trial = 0

if total_rows < int(INIT_NUM):
    work_type = "init"
    for ID in range(total_rows + 1, total_rows + 1 + int(NUM_CANDIDATES)):
        x, y = np.random.uniform(0, 1, 2)
        phi, theta, psi = np.random.uniform(-180, 180, 3)
        share_instance.make_ad(work_type, trial, ID, phi, theta, psi, x, y, filename)
        share_instance.ssh_connect(HOSTNAME, USERNAME, pkey, work_type, trial, ID)

elif total_rows < int(INIT_NUM) + int(SEARCH_NUM):
    work_type = "optuna"
    ID = total_rows + 1
    recommended_trials = share_instance.bayes_opt()
    for i, params in enumerate(recommended_trials):
        trial = (total_rows - int(INIT_NUM)) // int(NUM_CANDIDATES) + 1
        share_instance.make_ad(
            work_type,
            trial,
            ID,
            params["phi"],
            params["theta"],
            params["psi"],
            params["x"],
            params["y"],
            filename,
        )
        share_instance.ssh_connect(HOSTNAME, USERNAME, pkey, work_type, trial, ID)
        ID += 1

else:
    with open("job_time.txt", "a") as f:
        f.write(f"end: {pd.Timestamp.now()}\n")
    exit()
