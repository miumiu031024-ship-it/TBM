"""Two seminar slides: native PPT text + figures; matching PNG/PDF previews."""
from pathlib import Path
import json,math
from PIL import Image,ImageDraw,ImageFont,ImageOps
from pptx import Presentation
from pptx.util import Inches,Pt
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE

O=Path('/workspace/TBM/seminar-prep-20261125');F=O/'figures';S=O/'slides';S.mkdir(exist_ok=True)
W,H=1920,1080;ppi=120
FONT='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
BOLD='/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc'
INK='#14324b';MUTED='#587085';TEAL='#178d91';ORANGE='#b65f13'
prs=Presentation();prs.slide_width=Inches(16);prs.slide_height=Inches(9)
previews=[];layout=[]
def font(px,bold=False):return ImageFont.truetype(BOLD if bold else FONT,px)
def rgb(c):return RGBColor.from_string(c.lstrip('#'))
def new():
 slide=prs.slides.add_slide(prs.slide_layouts[6]);im=Image.new('RGB',(W,H),'white');return slide,im,ImageDraw.Draw(im)
def text(t,x,y,w,h,msg,size=28,color=INK,bold=False):
 slide,im,draw=t;ff=font(size,bold);lines=[]
 for para in msg.split('\n'):
  line=''
  for ch in para:
   if ff.getlength(line+ch)>w and line:lines.append(line);line=ch
   else:line+=ch
  lines.append(line)
 lineheight=math.ceil(size*1.4)
 assert len(lines)*lineheight<=h,(msg,len(lines)*lineheight,h)
 for i,line in enumerate(lines):draw.text((x,y+i*lineheight),line,font=ff,fill=color)
 box=slide.shapes.add_textbox(Inches(x/ppi),Inches(y/ppi),Inches(w/ppi),Inches(h/ppi))
 tf=box.text_frame;tf.clear();tf.word_wrap=False;tf.margin_left=tf.margin_right=tf.margin_top=tf.margin_bottom=0
 for i,line in enumerate(lines):
  p=tf.paragraphs[0] if i==0 else tf.add_paragraph();p.text=line;p.font.name='Noto Sans CJK JP';p.font.size=Pt(size*.6);p.font.bold=bold;p.font.color.rgb=rgb(color);p.space_before=p.space_after=Pt(0);p.line_spacing=Pt(lineheight*.6)
 layout.append(dict(slide=len(prs.slides),kind='text',box=[x,y,w,h],text=msg))
def card(t,x,y,w,h,color):
 slide,im,draw=t;draw.rounded_rectangle((x,y,x+w,y+h),radius=14,fill=color)
 shape=slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,Inches(x/ppi),Inches(y/ppi),Inches(w/ppi),Inches(h/ppi));shape.fill.solid();shape.fill.fore_color.rgb=rgb(color);shape.line.fill.background()
def pic(t,path,x,y,w,h):
 slide,im,draw=t;src=Image.open(path).convert('RGB');r=min(w/src.width,h/src.height);ww,hh=round(src.width*r),round(src.height*r);xx=x+(w-ww)//2;yy=y+(h-hh)//2
 im.paste(src.resize((ww,hh),Image.Resampling.LANCZOS),(xx,yy))
 slide.shapes.add_picture(str(path),Inches(xx/ppi),Inches(yy/ppi),width=Inches(ww/ppi),height=Inches(hh/ppi))
def header(t,num,title,sub):
 text(t,60,26,1800,67,title,44,bold=True)
 text(t,63,96,1800,42,sub,26,color=MUTED)
 t[2].line((60,148,1860,148),fill=TEAL,width=5)
 shape=t[0].shapes.add_shape(MSO_SHAPE.RECTANGLE,Inches(.5),Inches(148/ppi),Inches(15),Inches(5/ppi));shape.fill.solid();shape.fill.fore_color.rgb=rgb(TEAL);shape.line.fill.background()
 text(t,64,1044,1750,25,f'TBMゼミ用資料｜作成2026-11-25｜{num}/2　図は形状確認用。速度分布・改善を示すCFD結果ではない。',17,color=MUTED)
def finish(t,num,notes):
 if num==2 and (O/'mesh-quality/mesh-summary.json').exists():
  qq=json.loads((O/'mesh-quality/mesh-summary.json').read_text())
  notes+='\n実測の小さい行列式セル数：'+ '、'.join(f'{a["model"]} {a["small_determinant_cells"]}個' for a in qq)+'。通常検査と追加検査を区別し、該当セルは格子修正・再検査してから予備CFDへ進む。'
 t[0].notes_slide.notes_text_frame.text=notes
 p=S/f'slide-{num}.png';t[1].save(p);previews.append(t[1]);print(p.name,flush=True)

t=new();header(t,1,'可視化実験用模型の見直し','既存の水領域CADに近づけ、残留空気と可視化用気泡を分けて扱う')
pic(t,F/'B00-overview-and-ports.png',45,165,965,450)
pic(t,F/'B00-high-point-vent-concept.png',1010,175,865,400)
card(t,60,635,875,310,'#edf4f8');card(t,965,635,895,310,'#eaf5ef')
text(t,83,650,830,47,'現模型との違い・未確認事項',31,bold=True)
text(t,83,711,830,222,'現模型：15本の接続管 → 共通空間 → 冷却流路\nSTEP：18入口＋18出口／3独立水領域\n球面が上、下方の管から水・気泡を注入\n現模型の内部寸法・領域の接続・排気口は未確認',25,color=INK)
text(t,988,650,845,47,'新しい透明模型の検討方針',31,bold=True)
text(t,988,711,845,222,'観察部の追加共通空間を省き、36短管を維持\n透明内外壁・継手・分割部とシールを別設計\n3領域の高所排気と、気泡注入を分ける\n水路厚・漏れ・給水分配・観察面の光学性を確認',25)
card(t,60,965,1800,62,'#fff1dc');text(t,80,974,1760,43,'共通空間をなくしても、気泡は浮力で頂部へ集まり得る。気泡の動き＝水の速度とは限らない。',27,color=ORANGE,bold=True)
finish(t,1,'既存STEPは金属全体や3Dプリント用完成模型ではなく水領域。36短管、球面半径44.9/47.9mm、3mm水路厚を確認。現模型の15管と共通空間はユーザー説明によるもので、実験模型CADは未確認。球面が上、下の管から水と気泡を入れる向きと確認済み。排気は各領域の本当の高所に設ける案。測定中に排気を開けるか閉めるかで境界条件が変わるため記録する。低流量時は浮力により気泡が頂部に滞留し、下向きの水に追従しない可能性。透明模型の圧力・温度・材質は別設計し、290℃15.5MPaのCFD物性条件をそのまま樹脂模型の運転条件にしない。')

t=new();header(t,2,'冷却流路の構造変更案','隔壁の位置・角度が、出口への接近流れと分配に与える影響を比較する')
pic(t,F/'four-models-same-view.png',35,162,1850,340)
pic(t,F/'divider-outlet-detail.png',40,513,1280,440)
card(t,1335,516,525,400,'#edf4f8')
text(t,1358,531,480,43,'検証する仮説',30,bold=True)
text(t,1358,586,480,315,'B01：幅の影響と中央の損失\nC01/C02：±1°を逆向き比較\nC01：圧力回復と剥離の競合\nC02：掃き出しと摩擦の競合\n流量偏り・低速域・全圧損失を\n同時に評価する',25)
qpath=O/'mesh-quality/mesh-summary.json'
if qpath.exists():
 q=json.loads(qpath.read_text());counts=' ／ '.join(f'{a["model"]} {a["actual_cells"]/1e6:.2f}百万' for a in q)
 status=f'新規メッシュ：{counts}セル。B01は次段階。'
 if any(not a['standard_check_passed'] for a in q):tail='通常検査に未合格のケースあり。修正・再検査が必要。流れの計算は未実行。'
 elif any(a['small_determinant_cells'] for a in q):tail='追加検査で小さい行列式のセルあり。格子修正・再検査後に予備計算。今回は流れ未計算。'
 else:tail='通常検査は通過、追加検査に未合格項目。各200反復の予備計算は未実行。改善効果は未確定。' if any(not a['additional_check_passed'] for a in q) else '幾何検査は通過。各200反復の予備計算は未実行。流れの収束・改善効果は未確定。'
else:
 status='CFD準備：B00・C01・C02を同条件でメッシュ生成・品質確認中。B01は次段階。'
 tail='予備計算は各200反復を提案（今回未実行）。既存B00も未収束のため、短い計算で改善を確定しない。'
text(t,65,951,1790,45,status,25,bold=True)
text(t,65,1001,1790,32,tail,21,color=MUTED)
finish(t,2,'全モデルで球面外形、管の位置・内径3mm・直管長10mm、3領域を保つ。B01は隔壁を中央へ各0.50mm平行移動し、中央を1mm狭め左右を広げる。C01は右隔壁+1°、左-1°、C02は逆でY0の中心位置を固定。図の灰色は基準、橙は変更後。出口の図はXY投影で、水平面の実流体断面ではない。C01の拡大方向で静圧が回復しても全圧損失低下とは限らない。C02は左右を出口側で狭めるため速度上昇と低速域低減の仮説がある一方、摩擦・分配偏りが増す可能性。どちらも優劣は未検証。新規B00/C01/C02は同じ格子幅・根元細分化ルール・物性・解法で準備。明示的境界層メッシュなし、y+未評価、格子独立性未確認。今回はメッシュまででCFDを解いていない。20反復で速度・メモリを実測後、100一次精度+100二次精度の計200反復を予備試験として提案。定常の反復を物理時間とみなさない。長時間の計算は確認後。')
prs.save(S/'TBM-seminar-2slides.pptx')
previews[0].save(S/'TBM-seminar-2slides.pdf',save_all=True,append_images=previews[1:],resolution=120)
(S/'layout-audit.json').write_text(json.dumps(layout,indent=2,ensure_ascii=False))
print('PPTX/PDF/2PNG written',flush=True)
