#!/bin/bash
# 全データのデコードとマージを自動で行うシェルスクリプト

# エラーが起きたらスクリプトを停止する
set -e

echo "=== Starting Full Waveform Decoding and Merging Process ==="

# root/ ディレクトリの作成 (念のため)
mkdir -p root

# -----------------------------------------------------------------------------
# 1. 各ランの基本 TQ 解析および Coincidence 解析 (前処理を並列実行)
# -----------------------------------------------------------------------------
echo "--> Generating basic TQ and Coincidence files for Cf252 (Parallel)..."
for run in 01 02 03 04; do
    (
        echo "Generating TQ & Coin for Cf252 run ${run}..."
        ./bin/convert_to_root data/Cf252_tq_${run}.dat ${run} 0.0 root/Cf252_tq_${run}.root
        ./bin/coincidence_analysis root/Cf252_tq_${run}.root root/Cf252_tq_${run}_coincidence.root
    ) &
done

echo "--> Generating basic TQ files for Co60 (Parallel)..."
for run in 02 03; do
    (
        echo "Generating TQ for Co60 run ${run}..."
        ./bin/convert_to_root data/Co60_tq_${run}.dat ${run} 0.0 root/Co60_tq_${run}.root
    ) &
done

# 前処理プロセスの完了を待つ
wait
echo "--> Pre-processing completed."

# -----------------------------------------------------------------------------
# 2. Cf252 データの波形デコード (ラン 01 から 04 を並列実行)
# -----------------------------------------------------------------------------
echo "--> Decoding Cf252 waves 01-04 for gamma, fastn, and slown (Parallel)..."

for run in 01 02 03 04; do
    echo "Launching Cf252_wave_${run} decoders..."
    
    # gamma デコード (Cf252_wave_XX.dat -> root/Cf252_wave_XX_gamma.root)
    ./bin/export_waveform data/Cf252_wave_${run}.dat gamma &
    
    # fastn デコード (Cf252_wave_XX.dat -> root/Cf252_wave_XX_fastn.root)
    ./bin/export_waveform data/Cf252_wave_${run}.dat fastn &
    
    # slown デコード (Cf252_wave_XX.dat -> root/Cf252_wave_XX_slown.root)
    ./bin/export_waveform data/Cf252_wave_${run}.dat slown &
done

# -----------------------------------------------------------------------------
# 3. Co60 データのデコード (ラン 02 から 03 を並列実行)
# -----------------------------------------------------------------------------
echo "--> Decoding Co60 waves 02-03 for gamma only (Parallel)..."

for run in 02 03; do
    echo "Launching Co60_wave_${run} decoder..."
    
    # gamma デコード (Co60_wave_XX.dat -> root/Co60_wave_XX_gamma.root)
    ./bin/export_waveform data/Co60_wave_${run}.dat gamma &
done

# すべてのデコードプロセスの完了を待つ
wait
echo "--> Decoding completed."

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
