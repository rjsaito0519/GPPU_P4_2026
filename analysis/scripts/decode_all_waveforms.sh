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
# 1. Cf252 データの並列デコード (4並列: 2x2 画面分割)
# -----------------------------------------------------------------------------
echo "--> Launching Cf252 runs (01-04) in a 2x2 tmux window..."

# 新しい tmux セッションをデタッチドモードで開始
tmux new-session -d -s decode_wf -n "Cf252_Decoding"

# 2x2 の4分割ペインを作成
# 0: 左上 (初期ペイン)
# 右に分割して 1: 右上
tmux split-window -h -t decode_wf:0
# 左下を分割して 2: 左下
tmux split-window -v -t decode_wf:0.0
# 右下を分割して 3: 右下
tmux split-window -v -t decode_wf:0.1

# 各ペインでそれぞれのランを順次実行 (TQ -> Coin -> Decode x3 -> 完了フラグ)
for run in 01 02 03 04; do
    # 対応するペインのインデックス決定 (01->0.0, 02->0.1, 03->0.2, 04->0.3)
    if [ "$run" = "01" ]; then pane="0.0"; fi
    if [ "$run" = "02" ]; then pane="0.2"; fi # tmuxの分割順で左下
    if [ "$run" = "03" ]; then pane="0.1"; fi # 右上
    if [ "$run" = "04" ]; then pane="0.3"; fi # 右下

    cmd="./bin/convert_to_root data/Cf252_tq_${run}.dat ${run} 0.0 root/Cf252_tq_${run}.root && \
         ./bin/coincidence_analysis root/Cf252_tq_${run}.root root/Cf252_tq_${run}_coincidence.root && \
         ./bin/export_waveform data/Cf252_wave_${run}.dat gamma && \
         ./bin/export_waveform data/Cf252_wave_${run}.dat fastn && \
         ./bin/export_waveform data/Cf252_wave_${run}.dat slown && \
         touch root/.done_cf_${run}"

    tmux send-keys -t decode_wf:0.${pane} "${cmd}" C-m
done

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
    # 2. Co60 データの並列デコード (2並列: 左右2画面分割)
    # -----------------------------------------------------------------------------
    echo "--> Cf252 runs completed. Starting Co60 runs (02-03) in a new window..."
    
    # 新しいウィンドウを作成し、左右に分割
    tmux new-window -t decode_wf -n "Co60_Decoding"
    tmux split-window -h -t decode_wf:1
    
    # Co60 のコマンド
    cmd_co02="./bin/convert_to_root data/Co60_tq_02.dat 2 0.0 root/Co60_tq_02.root && \
              ./bin/export_waveform data/Co60_wave_02.dat gamma && \
              touch root/.done_co_02"
              
    cmd_co03="./bin/convert_to_root data/Co60_tq_03.dat 3 0.0 root/Co60_tq_03.root && \
              ./bin/export_waveform data/Co60_wave_03.dat gamma && \
              touch root/.done_co_03"
              
    tmux send-keys -t decode_wf:1.0 "${cmd_co02}" C-m
    tmux send-keys -t decode_wf:1.1 "${cmd_co03}" C-m
    
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
# 3. hadd によるマージ処理
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
