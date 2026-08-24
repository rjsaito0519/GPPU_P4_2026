import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# データの読み込み
dat_path = "data/Cf252_tq_01.dat"
if not os.path.exists(dat_path):
    print("Data file not found!")
    exit(1)

print("Loading data...")
data = np.loadtxt(dat_path)
T0 = data[:, 2]
Q0 = data[:, 3]
T1 = data[:, 4]
Q1 = data[:, 5]

# ゲート条件: (T1 - T0) が 30 ~ 60 の範囲、かつ Q0 > 1.5, Q1 > 0
tof = T1 - T0
mask = (tof >= 30.0) & (tof <= 60.0) & (Q0 > 1.5) & (Q1 > 0.0)

Q0_filtered = Q0[mask]
T0_filtered = T0[mask]

# 2Dヒストグラムのプロファイル（Q0のbinごとにT0の平均を計算）
bins = np.linspace(1.5, 50, 60)
bin_centers = (bins[:-1] + bins[1:]) / 2.0
profile_T0 = []
profile_Q0 = []

for i in range(len(bins)-1):
    bin_mask = (Q0_filtered >= bins[i]) & (Q0_filtered < bins[i+1])
    if np.sum(bin_mask) > 10:
        profile_T0.append(np.mean(T0_filtered[bin_mask]))
        profile_Q0.append(bin_centers[i])

profile_T0 = np.array(profile_T0)
profile_Q0 = np.array(profile_Q0)

# フィット関数: f(q) = p0 / sqrt(q) + p1
def slew_func(q, p0, p1):
    return p0 / np.sqrt(q) + p1

popt, pcov = curve_fit(slew_func, profile_Q0, profile_T0, p0=[30.0, 350.0])
p0, p1 = popt

print("\n==============================================")
print(" Slewing Correction Fit Results (Python)")
print("==============================================")
print(f"p0 (amplitude)  : {p0:.4f}")
print(f"p1 (offset)     : {p1:.4f}")
print("==============================================")

# プロット作成
plt.figure(figsize=(8, 6))
# 2Dヒストグラムを描画
plt.hist2d(Q0_filtered, T0_filtered, bins=[100, 100], range=[[0, 50], [320, 400]], cmap='viridis', cmin=1)
plt.colorbar(label='Entries')

# プロファイルとフィット曲線を重ね書き
plt.plot(profile_Q0, profile_T0, 'ro', label='Data Profile', markersize=4)
q_fit = np.linspace(1.5, 50, 200)
plt.plot(q_fit, slew_func(q_fit, p0, p1), 'r-', linewidth=3, label=f'Fit: {p0:.2f}/sqrt(Q0) + {p1:.2f}')

plt.xlabel('Q0')
plt.ylabel('T0')
plt.title('Slewing Fit (T1-T0 in 30-60)')
plt.legend()
plt.tight_layout()

# 保存
output_img = "data/slewing_fit_python.png"
plt.savefig(output_img, dpi=150)
print(f"Saved fit plot to {output_img}")
