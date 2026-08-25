#!/bin/bash
# 全データのデコードとマージを自動で行うシェルスクリプト

# エラーが起きたらスクリプトを停止する
set -e

echo "=== Starting Full Waveform Decoding and Merging Process ==="

# root/ ディレクトリの作成 (念のため)
mkdir -p root

# -----------------------------------------------------------------------------
# 1. Cf252 データのデコード (ラン 01 から 04)
# -----------------------------------------------------------------------------
echo "--> Decoding Cf252 waves 01-04 for gamma, fastn, and slown..."

for run in 01 02 03 04; do
    echo "Processing Cf252_wave_${run}..."
    
    # gamma デコード (Cf252_wave_XX.dat -> root/Cf252_wave_XX_gamma.root)
    ./bin/export_waveform data/Cf252_wave_${run}.dat gamma
    
    # fastn デコード (Cf252_wave_XX.dat -> root/Cf252_wave_XX_fastn.root)
    ./bin/export_waveform data/Cf252_wave_${run}.dat fastn
    
    # slown デコード (Cf252_wave_XX.dat -> root/Cf252_wave_XX_slown.root, 内部で delta_T_us <= 500 us カット自動適用)
    ./bin/export_waveform data/Cf252_wave_${run}.dat slown
done

# -----------------------------------------------------------------------------
# 2. Co60 データのデコード (ラン 02 から 03)
# -----------------------------------------------------------------------------
echo "--> Decoding Co60 waves 02-03 for gamma only..."

for run in 02 03; do
    echo "Processing Co60_wave_${run}..."
    
    # gamma デコード (Co60_wave_XX.dat -> root/Co60_wave_XX_gamma.root)
    ./bin/export_waveform data/Co60_wave_${run}.dat gamma
done

# -----------------------------------------------------------------------------
# 3. hadd によるマージ処理 (同じ測定条件のランファイルを結合)
# -----------------------------------------------------------------------------
echo "--> Merging ROOT files using hadd..."

# Cf252 gamma マージ
hadd -f root/Cf252_wave_merge_gamma.root \
    root/Cf252_wave_01_gamma.root \
    root/Cf252_wave_02_gamma.root \
    root/Cf252_wave_03_gamma.root \
    root/Cf252_wave_04_gamma.root

# Cf252 fastn マージ
hadd -f root/Cf252_wave_merge_fastn.root \
    root/Cf252_wave_01_fastn.root \
    root/Cf252_wave_02_fastn.root \
    root/Cf252_wave_03_fastn.root \
    root/Cf252_wave_04_fastn.root

# Cf252 slown マージ
hadd -f root/Cf252_wave_merge_slown.root \
    root/Cf252_wave_01_slown.root \
    root/Cf252_wave_02_slown.root \
    root/Cf252_wave_03_slown.root \
    root/Cf252_wave_04_slown.root

# Co60 gamma マージ
hadd -f root/Co60_wave_merge_gamma.root \
    root/Co60_wave_02_gamma.root \
    root/Co60_wave_03_gamma.root

echo "=== All processes completed successfully! ==="
echo "Merged files:"
echo " - root/Cf252_wave_merge_gamma.root"
echo " - root/Cf252_wave_merge_fastn.root"
echo " - root/Cf252_wave_merge_slown.root (purity optimized: < 500 us)"
echo " - root/Co60_wave_merge_gamma.root"
