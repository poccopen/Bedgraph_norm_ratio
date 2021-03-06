## Bedgraph_norm_ratio.py
※出芽酵母 <i>Saccharomyces cerevisiae</i> 専用です。

2つのbedgraphファイルを入力として、ゲノムの各座標におけるリードカウントの比を計算します。

リードカウント比はそれぞれのbedgraphに含まれるリードカウント総数でノーマライズされます。

ミトコンドリアゲノムにマップされたリードはリードカウント総数に含めません。

#### [依存性]
Biopythonがインストールされている必要があります。

https://biopython.org/

https://biomedicalhacks.com/2020-05-12/biopython-basic-1/

#### [使い方]  
```$ python3 Bedgraph_norm_ratio.py [Reference_genome_seq.fasta] [Sample_data.bedgraph] [Control_data.bedgraph] [Output_file.bedgraph]```

#### [入力ファイル]
入力ファイルは以下の3つです。

Reference_genome_seq.fasta FASTA形式のリファレンスゲノム配列です。染色体名の表記は"chrI"型のみに対応しています。

Sample_data.bedgraph サンプル（比を計算する際に分子となる）となるBedgraph形式のファイルです。

Control_data.bedgraph 比較基準（比を計算する際に分母となる）となるBedgraph形式のファイルです。

#### [出力ファイル]
Bedgraph形式のファイルをひとつ出力します。(Output_file.bedgraph)

ゲノムの各座標における サンプルリードカウント/コントロールリードカウント比（リードカウント総数でノーマライズ） が記述されています。

ミトコンドリアゲノム部分についてはデータを出力しません。

## Bedgraph_normalize.py
1つのbedgraphファイルを入力として、ゲノムの各座標におけるリードカウントをリードカウント総数でノーマライズします。

ミトコンドリアゲノムにマップされたリードはリードカウント総数に含めません。

#### [使い方]  
```$ python3 Bedgraph_norm_ratio.py [Reference_genome_seq.fasta] [Sample_data.bedgraph] [Output_file.bedgraph]```

#### [入力ファイル]
入力ファイルは以下の2つです。

Reference_genome_seq.fasta FASTA形式のリファレンスゲノム配列です。染色体名の表記は"chrI"型のみに対応しています。

Sample_data.bedgraph ノーマライズ対象とするBedgraph形式のファイルです。
