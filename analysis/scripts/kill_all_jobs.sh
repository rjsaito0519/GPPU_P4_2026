#!/bin/bash
# すべてのバックグラウンドのデコード・解析プロセスおよび tmux セッションを一撃で強制終了するスクリプト

echo "=== Terminating all decoding and analysis processes ==="

# 1. tmux セッションを強制終了
if tmux has-session -t decode_wf 2>/dev/null; then
    echo "Killing tmux session: decode_wf"
    tmux kill-session -t decode_wf 2>/dev/null || true
fi

# 2. 関連プログラムのプロセスを一括で強制終了 (SIGKILL)
echo "Killing any remaining binaries..."
pkill -9 -f "convert_to_root|coincidence_analysis|export_waveform|psd_analysis|calculate_fom" 2>/dev/null || true

# 3. 完了フラグ用の一時ファイルをクリーンアップ
echo "Cleaning up temporary flag files..."
rm -f root/.done_* 2>/dev/null || true

echo "=== Cleanup completed successfully! ==="
