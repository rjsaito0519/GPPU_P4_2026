#!/bin/bash

# ==============================================================================
# Cf252 slown デコード・マージ専用スクリプト
# ==============================================================================

# Ctrl+Cで中断された際、バックグラウンドプロセスをきれいにkillする処理
cleanup() {
    echo -e "\n=== Aborting: Killing background jobs... ==="
    pkill -P $$
    exit 1
}
trap cleanup SIGINT SIGTERM

echo "=== Starting Cf252 slown Waveform Decoding (Parallel) ==="

# 出力先フォルダの確認
mkdir -p root

# 1. 01~04ランのデコードをバックグラウンド並列で開始
pids=()
for run in 01 02 03 04; do
    echo " -> Launching slown decoding for Run ${run} (logging to root/log_decode_cf_${run}_slown.log)..."
    ./bin/export_waveform data/Cf252_wave_${run}.dat slown root/Cf252_wave_${run}_slown.root > root/log_decode_cf_${run}_slown.log 2>&1 &
    pids+=($!)
done

# 2. すべてのバックグラウンドプロセスが完了するまで待機
echo " -> Waiting for all decoding processes to finish..."
for pid in "${pids[@]}"; do
    wait $pid
done

echo " -> All decoding processes completed successfully."

# 3. hadd でマージを実行
echo "=== Merging ROOT files using hadd ==="
hadd -f root/Cf252_wave_merge_slown.root \
        root/Cf252_wave_01_slown.root \
        root/Cf252_wave_02_slown.root \
        root/Cf252_wave_03_slown.root \
        root/Cf252_wave_04_slown.root

echo -e "\n=== Process Complete ==="
echo " -> Merged file saved to: root/Cf252_wave_merge_slown.root"
echo "=========================================================="
