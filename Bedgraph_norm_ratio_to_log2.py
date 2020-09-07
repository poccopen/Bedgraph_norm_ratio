# -*- coding: utf-8 -*-
# argvを取得するためにsysモジュールをインポートする
import sys
# 正規表現を使用するためにreモジュールをインポートする
import re
# log計算のためにmathモジュールをインポートする
import math

# コマンドライン引数をargvs（リスト）に格納する
argvs = sys.argv
# コマンドライン引数の数を変数argcに格納する
argc = len(argvs)

# 入力するファイル群が指定されていないときは使い方を表示して終了する
if argc < 3:
	print("Usage: python3 {} [sample_data.norm.ratio.bedgraph] [output.norm.log2.bedgraph]".format(argvs[0]))
else:
	# それぞれのファイル名を格納する
	sample_data_name = argvs[1]
	output_name = argvs[2]

	#サンプルデータを読み込む
	with open(sample_data_name) as sample_file:
		for line in sample_file:
			line = line.strip()
			line = line.split()
			chr = line[0]
			start = line[1]
			end = line[2]
			read_count = line[3]
			log2_read_count = math.log2(float(read_count))

			output = open(output_name, 'a')
			output.write(str(chr) + "\t" + str(start) + "\t" + str(end) + "\t" + str(log2_read_count) + "\n")
