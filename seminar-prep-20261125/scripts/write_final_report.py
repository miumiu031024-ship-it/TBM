"""Write measured mesh preparation report; never invent flow results."""
from pathlib import Path
import json,csv,math
O=Path('/workspace/TBM/seminar-prep-20261125');D=O/'mesh-quality';R=Path(__file__).resolve().parent
rows=json.loads((D/'mesh-summary.json').read_text())
tables=[];resource=[]
for q in rows:
 n=q['model'];p=json.loads((D/f'{n}-port-geometry-summary.json').read_text())
 tables.append(f'| {n} | {q["actual_cells"]:,} | {"合格" if q["standard_check_passed"] else "未合格"} | {"合格" if q["additional_check_passed"] else "未合格"} | {q["concave_cells"]:,} | {q["small_determinant_cells"]} | {q["max_nonorthogonality_deg"]:.2f}° | {q["max_skewness"]:.3f} |')
 resource.append(f'| {n} | {q["actual_cells"]/1e6:.3f}百万（実数） | {q["mesh_generation_and_checks_wall_minutes"]:.1f}分（実測） | {q["sampled_container_memory_current_peak_GiB"]:.2f} GiB（5秒間隔の観測最大） | 約30〜50分（未実行・既存ログからの見積） |')
 assert p['ports']==36 and p['all_areas_positive']
peak=max(q['cumulative_container_memory_peak_at_end_GiB'] for q in rows)
vols='、'.join(f'{q["model"]} {q["mesh_volume_error_percent"]:+.3f}%' for q in rows)
ports='、'.join(f'{q["model"]} {json.loads((D/(q["model"]+"-port-geometry-summary.json")).read_text())["maximum_absolute_area_error_percent"]:.2f}%' for q in rows)
portcomp=json.loads((D/'matched-port-geometry.json').read_text())
portdiff='、'.join(f'{q["model"]} {q["maximum_absolute_native_port_area_difference_percent"]:.6g}%' for q in portcomp)
warning=any(not q['additional_check_passed'] or not q['standard_check_passed'] for q in rows)
if any(not q['standard_check_passed'] for q in rows):
 verdict='通常検査に未合格のケースがあります。未合格項目の修正・再検査が必要で、性能評価用の計算を開始できる品質とは扱いません。'
elif any(q['small_determinant_cells'] for q in rows):
 verdict='通常検査は通過しましたが、追加検査で凹形状セルに加えて、小さなセル形状行列式のセルが見つかりました。該当セルを修正・再検査してから構造比較の予備計算へ進む必要があります。今回のメッシュを性能評価用に準備完了とは扱いません。'
elif warning:
 verdict='通常検査は通過していますが、追加の幾何検査に未合格項目があります。今回の格子は予備比較の準備段階で、性能評価用の品質が確立したとは扱いません。'
else:
 verdict='通常・追加の幾何検査は通過しました。これだけで流れの収束や格子独立性が確認されたとは扱いません。'
(D/'README.md').write_text(f'''# 新規メッシュの実測結果

既存B00/C01/C02 STEPから、新規ケースで生成したメッシュです。旧CAD・旧メッシュ・旧計算結果を上書きしていません。流れは計算していません。

**{verdict}**

| モデル | 実セル数 | 通常検査 | 追加検査 | 凹形状セル | 小さい行列式のセル | 最大非直交角 | 最大スキュー |
|---|---:|---|---|---:|---:|---:|---:|
{chr(10).join(tables)}

全ケースで水領域3、cellZone3、境界パッチ37（壁＋36開口）、36開口すべての面数と面積が正であること、最小セル体積が正であることを確認しました。3領域は設計上の独立水路です。領域数が3であることを欠陥とは扱いません。

CAD体積に対するネイティブメッシュ体積の誤差は{vols}です。管端開口面積のCADに対する絶対誤差の最大値は{ports}でした。モデル間の管端面積差も出口分配比較の数値的不確かさとして確認する必要があります。STL変換の体積検査も各ケースに記録しました。

新しいB00に対する管端面積のモデル間差の絶対値の最大は、{portdiff}でした。これは各案とCADとの誤差とは別の比較です。[モデル間の開口比較](matched-port-geometry.json)に保存しています。これが小さくても内部格子の違いや格子独立性の検証を置き換えるものではありません。

明示的なプリズム境界層は生成していません。y+は流れを解くまで求められません。凹形状セルや非直交角だけから未収束の原因を特定せず、位置・壁面処理・格子感度を別に確認します。共通の設定でも各形状の格子品質は同じとは限りません。

小さい行列式の検査基準は0.001です。体積が負という意味ではなく、形状がつぶれたセル等を検出する指標です。少数でも比較の数値誤差や不安定化の懸念があるため、位置を確認し、隔壁付近の格子遷移や特徴線へのスナップなどを新規ケースで調整して再検査します。スナップ処理の最終エラー数0だけを合格条件とせず、独立した全幾何検査も確認します。この問題だけで、以前のB00やSTARの未収束原因を特定したとは言えません。

[問題セルの位置と根元への距離](small-determinant-cell-locations.csv)、[位置の集計](small-determinant-cell-location-summary.json)も保存しました。抽出したVTUのcellIDは部分メッシュ内の番号で、元のグローバルセル番号と同一視していません。元のcellSetも別に保持しています。抽出したのはセル形状だけで、速度・圧力の場は含めていません。

![問題セルの位置](../figures/quality-cell-locations.png)

## 資源と時間

| モデル | セル数 | メッシュ生成＋再構築＋検査 | コンテナメモリ | 将来の200反復 |
|---|---:|---:|---:|---|
{chr(10).join(resource)}
| B01 | 約4.5〜6.5百万（見積、未生成） | 類似規模を想定、未実測 | 実行時16〜24 GiB確保を想定 | 約30〜50分（未実行・見積） |

メッシュ生成・品質検査の終了時点までに記録したコンテナの累積cgroup memory.peakは{peak:.2f} GiBでした。これはページキャッシュを含み、個々のプロセスRSSの合計ではありません。各ケースの5秒間隔の観測最大と、終了時点の累積ピークを区別しています。コンテナ上限は24 GiB、CPU上限4、4 MPIで逐次処理しました。流れの計算には格子生成と異なるメモリが必要なため、1ケース16〜24 GiBを確保し、最初の20反復で実測を更新する案です。

メッシュ生成時間は画像作成なども並行したこの環境での実測で、単独実行や別の計算機の所要時間を保証する値ではありません。

流れの実行時間は今回測っていません。以前の約500万セル・4 MPIのログの約9〜12秒/反復から、初期化・保存等を含め200反復は約30〜50分/ケース、3ケース逐次で約1.5〜2.5時間を見込みます。ハードウェアや線形反復数で変わります。

## 次の試験に進む条件

今回はメッシュ生成・検査までです。上記の品質問題を修正・再検査してから、各20反復で動作・資源確認後、一次100＋二次100の計200反復以内を予備試験として提案します。未収束の短い計算の面積差は改善の確定値としません。長時間の定常継続やURANSは確認後です。

メッシュのネイティブ本体は作業環境 `/workspace/scratch/tbm-seminar-prep-20261125/B00`、`C01`、`C02` にあります。配布ZIPには大容量polyMesh・triSurfaceを含めず、条件・メタデータ・全メッシュログ・品質集計・36開口の幾何検査を含めます。

[数値一覧CSV](mesh-summary.csv)、[詳細JSON](mesh-summary.json)、[共通条件の一致確認](matched-input-validation.json)。
''')
(O/'README.md').write_text(f'''# TBMゼミ資料：透明模型の見直しと構造変更の比較準備

作成日：2026年11月25日。**画像・2枚のスライドと、新規B00/C01/C02のメッシュ生成・品質確認を完了しました。流れのCFD計算は実行していません。** 元のCAD・既存メッシュ・既存結果は保全しています。

## ゼミ用の2枚

- [PowerPoint：2枚、本文テキストと図・発表ノート](slides/TBM-seminar-2slides.pptx)
- [PDF：2ページ](slides/TBM-seminar-2slides.pdf)
- [スライド1 PNG](slides/slide-1.png)、[スライド2 PNG](slides/slide-2.png)
- [発表時の説明例](speaker-notes.md)

### スライド1：可視化実験用模型の見直し

![スライド1](slides/slide-1.png)

B00は18入口＋18出口、3独立水領域、隔壁2枚からなる水領域CADです。説明された実験模型の追加15接続管・共通空間は含みません。実験模型は球面が上で、下方の管から水・気泡を入れるとの説明に基づいて整理しました。内部の実際の流速方向を測定した情報ではありません。内部寸法・接続詳細・排気口は未確認として区別しました。

透明模型には内外壁・継手・分割/密閉・領域間漏れ防止を別途設計する必要があります。共通空間を除いても気泡は浮力で頂部へ集まり得ます。残留空気の排出と可視化用気泡の注入、気泡と水の速度を区別します。製作用CADは今回作成していません。

[実験模型と透明容器の詳細検討](experiment-model-design.md)。

### スライド2：冷却流路の構造変更案

![スライド2](slides/slide-2.png)

B01は中央へ各0.50 mm平行移動、C01/C02は±1°の逆方向の回転です。B00・B01・C01・C02を同じ正投影の視点・縮尺で並べ、OUT_06付近の投影図で隔壁位置の差を拡大しました。図の色は水・隔壁の区別で、CFDの速度ではありません。

## メッシュ品質とCFD準備

{verdict}

| モデル | 実セル数 | 通常検査 | 追加検査 | 凹形状セル | 小さい行列式のセル | 最大非直交角 | 最大スキュー |
|---|---:|---|---|---:|---:|---:|---:|
{chr(10).join(tables)}

B01は画像比較・既存CADの確認までで、今回のメッシュ生成対象には含めていません。壁面y+・流れの収束・18出口流量・低速面積・全圧損失は今回未評価です。初期場を結果として表示していません。

[実測メッシュ品質と必要資源](mesh-quality/README.md)、[物性・境界条件・6評価項目・予備計算の計画](CFD-comparison-plan.md)。

品質問題の修正・再検査後の予備計算は各20反復で動作・資源を確認し、その後は計200反復以内（一次100＋二次100）を提案します。以前の約500万セル・4 MPIのログから、各200反復で約30〜50分、3ケース逐次で約1.5〜2.5時間を見込みます。今回は未実行です。長時間のCFDは確認後に進めます。

## 個別の高解像度図

- [B00全体図と36開口](figures/B00-overview-and-ports.png)
- [B00の実CAD断面・隔壁付け根](figures/B00-real-CAD-section.png)
- [各水領域の高所と排気の検討位置](figures/B00-high-point-vent-concept.png)
- [4モデルの同じ視点・縮尺の比較](figures/four-models-same-view.png)
- [出口付近の隔壁位置・幅の比較](figures/divider-outlet-detail.png)
- [OUT01〜06と隔壁の位置関係](figures/outlet01-to06-geometry-context.png)
- [生成した出口端面メッシュの比較](figures/native-outlet-mesh-check.png)
- [追加検査の問題セルの位置](figures/quality-cell-locations.png)
- [実験模型と新模型の接続関係の概念図](figures/experiment-model-concept.png)

[領域別のCAD水容積](geometry-volume-by-region.csv)も確認しました。B01は左右の容積を約2.03%増やし、中央を約2.11%減らします。C01/C02は各領域の容積がほぼ同じで、前後の形状分布を逆向きに変える比較です。これは流量の計算値ではありません。

## 再現性と保全

既存STEP・以前の解析結果は変更していません。入力のハッシュ、共通条件の一致確認、メッシュ生成と品質検査の全ログ、実際に使ったスクリプトを同梱します。元STEPのソースは既存 `cad-proposals-20261120/water/` です。大容量の体積メッシュ・表面STLは配布ZIPには含めず、作業環境に保持します。既存保存解のコピーや追加の流れの計算は行っていません。
''')
(O/'speaker-notes.md').write_text('''# ゼミでの説明例

## スライド1

既存STEPを確認すると、球面の薄い水路が隔壁で3領域に分かれ、18本から入り18本から出る構造でした。実験模型で説明されている15本の接続管と大きな共通空間は、このSTEPにはありません。そこで観察部を元の流路に近づけ、給排水系を別に整理する案を検討します。ただし水領域CADだけでは密閉できる透明模型にならず、内外壁、継手、分割部とシールが必要です。また球面を上にして管から気泡を入れるので、浮力で頂部へ集まる現象は共通空間を除いても残り得ます。残留空気の排気と可視化用気泡の注入を分け、水の流れと気泡の動きを区別して観察する方針です。

## スライド2

構造案は36本の管と3領域を維持し、隔壁だけを小さく変えます。B01は中央へ各0.50 mm動かす案、C01とC02は±1°を逆向きに回す案です。全体図では差が小さいため、管を固定した出口付近の拡大図も作りました。C01では減速・圧力回復と剥離、C02では掃き出しと摩擦増大がどう競合するかを調べます。まずB00・C01・C02を同じ条件でメッシュ化し、品質を確認しました。通常検査は通過しましたが、C01とC02の付け根付近には小さいセル形状行列式のセルが各6個あり、格子修正・再検査を先に行う必要があります。流れはまだ解いていません。その後、18出口の流量、OUT01と06の差、低速域、全圧損失、接近流れと反復中の変化を同時に比較します。既存B00も未収束なので、短い予備計算だけで形状の改善とは結論しません。
''')
print('Measured final report and resource estimates written')
