# Bedgraph_norm_ratio
※出芽酵母 Saccharomyces cerevisiae 専用です。
2つのbedgraphファイルを入力として、ゲノムの各座標におけるリードカウントの比を計算します。
リードカウント比はそれぞれのbedgraphに含まれるリードカウント総数でノーマライズされます。
ミトコンドリアゲノムにマップされたリードはリードカウント総数に含めません。

## GATC_TC_extraction.py
#### [使い方]  
```$ python3 GATC_TC_extraction.py [filename.sam]```

#### [目的]
読み始めの配列が"GATC"であるリードと"TC"であるリードを別々に集計します。

#### [入力ファイル]
bowtie2によって出力されたマッピング後のSAMファイルを入力ファイルとします。

#### [出力ファイル]
Bedgraph形式のファイルをふたつ出力します。
- filename.GATC.bedgraph （リードカウントは"G"の座標にアサインされます）
- filename.TC.bedgraph （リードカウントは"T"の座標にアサインされます）
