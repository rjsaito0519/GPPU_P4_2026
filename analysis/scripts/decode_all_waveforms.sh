#!/bin/bash
# tmux を使用して画面分割し、いい感じに並列でデコードおよびマージを行うシェルスクリプト

# エラーが起きたらスクリプトを停止する
set -e

# Ctrl+C などの終了シグナルを受け取ったときに、tmuxセッションもろとも終了させる
cleanup() {
    echo "Aborting... Cleaning up tmux session and temporary files."
    tmux kill-session -t decode_wf 2>/dev/null || true
    rm -f root/.done_* 2>/dev/null || true
    exit 1
}
trap cleanup SIGINT SIGTERM EXIT

echo "=== Starting Full Waveform Decoding and Merging Process (tmux Parallel) ==="

# root/ ディレクトリの作成
mkdir -p root
rm -f root/.done_* 2>/dev/null || true

# tmux がインストールされているか確認
if ! command -v tmux &> /dev/null; then
    echo "ERROR: tmux is not installed. Please install tmux (e.g. sudo apt install tmux)."
    exit 1
fi

# すでに tmux セッションの中にいる場合は無限再帰を防ぐため警告して終了
if [ -n "$TMUX" ]; then
    echo "ERROR: Please run this script outside of any active tmux session."
    exit 1
fi

# 過去のセッションがあれば念のため kill
tmux kill-session -t decode_wf 2>/dev/null || true

# -----------------------------------------------------------------------------
# 1. 各ランの基本 TQ 解析および Coincidence 解析 (前処理を並列実行、ログ保存)
# -----------------------------------------------------------------------------
echo "--> Generating basic TQ and Coincidence files for Cf252 (Parallel, logging to root/log_pre_cf_XX.log)..."
for run in 01 02 03 04; do
    (
        echo "Generating TQ & Coin for Cf252 run ${run}..."
        (./bin/convert_to_root data/Cf252_tq_${run}.dat ${run} 0.0 root/Cf252_tq_${run}.root && \
         ./bin/coincidence_analysis root/Cf252_tq_${run}.root root/Cf252_tq_${run}_coincidence.root) 2>&1 | tee root/log_pre_cf_${run}.log
    ) &
done

echo "--> Generating basic TQ files for Co60 (Parallel, logging to root/log_pre_co_XX.log)..."
for run in 02 03; do
    (
        echo "Generating TQ for Co60 run ${run}..."
        ./bin/convert_to_root data/Co60_tq_${run}.dat ${run} 0.0 root/Co60_tq_${run}.root 2>&1 | tee root/log_pre_co_${run}.log
    ) &
done

# 前処理プロセスの完了を待つ
wait
echo "--> Pre-processing completed."

# -----------------------------------------------------------------------------
# 2. Cf252 データの並列デコード (4並列: 2x2 画面分割、ログ保存)
# -----------------------------------------------------------------------------
echo "--> Launching Cf252 runs (01-04) in a 2x2 tmux window (logging to root/log_decode_cf_XX.log)..."

# 各ラン用のデコードコマンドを定義 (tee により画面表示しつつログファイルにも記録)
cmd_cf01="(./bin/export_waveform data/Cf252_wave_01.dat gamma && ./bin/export_waveform data/Cf252_wave_01.dat fastn && ./bin/export_waveform data/Cf252_wave_01.dat slown) 2>&1 | tee root/log_decode_cf_01.log; touch root/.done_cf_01"
cmd_cf02="(./bin/export_waveform data/Cf252_wave_02.dat gamma && ./bin/export_waveform data/Cf252_wave_02.dat fastn && ./bin/export_waveform data/Cf252_wave_02.dat slown) 2>&1 | tee root/log_decode_cf_02.log; touch root/.done_cf_02"
cmd_cf03="(./bin/export_waveform data/Cf252_wave_03.dat gamma && ./bin/export_waveform data/Cf252_wave_03.dat fastn && ./bin/export_waveform data/Cf252_wave_03.dat slown) 2>&1 | tee root/log_decode_cf_03.log; touch root/.done_cf_03"
cmd_cf04="(./bin/export_waveform data/Cf252_wave_04.dat gamma && ./bin/export_waveform data/Cf252_wave_04.dat fastn && ./bin/export_waveform data/Cf252_wave_04.dat slown) 2>&1 | tee root/log_decode_cf_04.log; touch root/.done_cf_04"

# 新しい tmux セッションをデタッチドモードで開始 (初期ペインで Cf252 run 01 を実行)
tmux new-session -d -s decode_wf -n "Cf252_Decoding" "${cmd_cf01}"
sleep 0.5

# 2x2 の4分割ペインを作成し、それぞれのペインで該当のコマンドを実行
# 右に分割して Cf252 run 03 を実行
tmux split-window -h -t decode_wf:0 "${cmd_cf03}"
sleep 0.2
# 左下を分割して Cf252 run 02 を実行
tmux split-window -v -t decode_wf:0.0 "${cmd_cf02}"
sleep 0.2
# 右下を分割して Cf252 run 04 を実行
tmux split-window -v -t decode_wf:0.1 "${cmd_cf04}"
sleep 0.2

# 進捗を視覚的に確認できるようにアタッチ
echo "Attaching to tmux. You will see 4 divided screens showing decoding progress."
echo "Once Cf252 runs are complete, the script will automatically proceed to Co60 runs."
echo "Press Ctrl+C to abort everything."
sleep 1

# バックグラウンドでアタッチ待ちしつつ、Cf252の完了をポーリング監視する
(
    while [ ! -f root/.done_cf_01 ] || [ ! -f root/.done_cf_02 ] || [ ! -f root/.done_cf_03 ] || [ ! -f root/.done_cf_04 ]; do
        sleep 2
    done
    
    # -----------------------------------------------------------------------------
    # 3. Co60 データの並列デコード (2並列: 左右2画面分割、ログ保存)
    # -----------------------------------------------------------------------------
    echo "--> Cf252 runs completed. Starting Co60 runs (02-03) in a new window..."
    
    # Co60 のデコードコマンド
    cmd_co02="./bin/export_waveform data/Co60_wave_02.dat gamma 2>&1 | tee root/log_decode_co_02.log; touch root/.done_co_02"
    cmd_co03="./bin/export_waveform data/Co60_wave_03.dat gamma 2>&1 | tee root/log_decode_co_03.log; touch root/.done_co_03"

    # 新しいウィンドウを作成し、左右に分割 (初期ペインで Co60 run 02 を実行)
    tmux new-window -t decode_wf -n "Co60_Decoding" "${cmd_co02}"
    sleep 0.5
    tmux split-window -h -t decode_wf:1 "${cmd_co03}"
    
    while [ ! -f root/.done_co_02 ] || [ ! -f root/.done_co_03 ]; do
        sleep 2
    done
    
    # 全デコード完了フラグ
    touch root/.done_all
    
    # tmux セッションを正常終了させる
    tmux kill-session -t decode_wf
) &

# アタッチして完了まで監視
tmux attach-session -t decode_wf || true

# アタッチが切れた、または完了フラグを待つ
trap - SIGINT SIGTERM EXIT # 正常終了したため、トラップをクリア
rm -f root/.done_* 2>/dev/null || true

# -----------------------------------------------------------------------------
# 4. hadd によるマージ処理
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
