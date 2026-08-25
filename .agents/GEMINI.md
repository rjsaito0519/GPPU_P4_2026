# GPPU TOF Project Coding Rules

このプロジェクトにおける ROOT 解析コードの記述とクリーンアップに関する統一ルールです。

## 1. 型の定義
ROOTライブラリのブランチバインドや一般的なループ変数、物理量においては、標準C++のビルトイン型ではなく、以下のROOT準拠の型を明示的に使用します。

* `Int_t` (C++ `int`)
* `Double_t` (C++ `double`)
* `Bool_t` (C++ `bool`)
* `ULong64_t` (C++ `unsigned long long`)
* `UShort_t` (C++ `unsigned short`)
* `Char_t` (C++ `char`)

## 2. TTree の命名規則
出力する ROOT ファイル内の TTree 名は、混乱を防ぐためすべて一律で **`"tree"`** に統一します。
* `wave_tree` $\to$ `"tree"`
* `psd_tree` $\to$ `"tree"`
* `coincidence_tree` $\to$ `"tree"`
* それに伴い、ツリーを読み込む側のコード（`Get("wave_tree")` 等）もすべて `"tree"` を指すように記述します。

## 3. 引数の設計方針
オプション引数はむやみに増やさず、シンプルで直感的な指定方法を保ちます。

## 4. コードの可読性と整理
* デバッグ用に出力した不要なコメントアウトコード（残骸）は可能な限り削除し、クリーンなソースコードを維持します。
* 変数名や関数名には統一感（スネークケースまたはキャメルケースなど）を持たせ、共通の変数名（例: `event`, `time_stamp`, `channel` 等）を使用します。
