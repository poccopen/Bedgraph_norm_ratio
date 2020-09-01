# -*- coding: utf-8 -*-
# argvを取得するためにsysモジュールをインポートする
import sys
# 正規表現を使用するためにreモジュールをインポートする
import re
# 外部コマンドを使用するためにsubprocessモジュールをインポートする
import subprocess
# FASTAファイルを取り扱うためにbiopythonからSeqIOモジュールをインポートする
from Bio import SeqIO
# 配列を取り扱いやすくするためにnumpyモジュールをインポートする
import numpy as np

# コマンドライン引数をargvs（リスト）に格納する
argvs = sys.argv
# コマンドライン引数の数を変数argcに格納する
argc = len(argvs)

# 入力するファイル群が指定されていないときは使い方を表示して終了する
if argc < 5:
	print("Usage: python3 {} [reference_genome_seq.fasta] [sample_data.bedgraph] [control_data.bedgraph] [output.bedgraph]".format(argvs[0]))
else:
	# それぞれのファイル名を格納する
	reference_genome_seq_name = argvs[1]
	sample_data_name = argvs[2]
	control_data_name = argvs[3]
	output_name = argvs[4]

	#各染色体の塩基長を初期化する
	chrI_len = 0
	chrII_len = 0
	chrIII_len = 0
	chrIV_len = 0
	chrV_len = 0
	chrVI_len = 0
	chrVII_len = 0
	chrVIII_len = 0
	chrIX_len = 0
	chrX_len = 0
	chrXI_len = 0
	chrXII_len = 0
	chrXIII_len = 0
	chrXIV_len = 0
	chrXV_len = 0
	chrXVI_len = 0
	chrmt_len = 0

	#リファレンスゲノム配列を読み込んで各染色体の塩基長を得る
	for record in SeqIO.parse(reference_genome_seq_name, 'fasta'):
		if record.id == "chrI":
			chrI_len = len(record.seq)
		elif record.id == "chrII":
			chrII_len = len(record.seq)
		elif record.id == "chrIII":
			chrIII_len = len(record.seq)
		elif record.id == "chrIV":
			chrIV_len = len(record.seq)
		elif record.id == "chrV":
			chrV_len = len(record.seq)
		elif record.id == "chrVI":
			chrVI_len = len(record.seq)
		elif record.id == "chrVII":
			chrVII_len = len(record.seq)
		elif record.id == "chrVIII":
			chrVIII_len = len(record.seq)
		elif record.id == "chrIX":
			chrIX_len = len(record.seq)
		elif record.id == "chrX":
			chrX_len = len(record.seq)
		elif record.id == "chrXI":
			chrXI_len = len(record.seq)
		elif record.id == "chrXII":
			chrXII_len = len(record.seq)
		elif record.id == "chrXIII":
			chrXIII_len = len(record.seq)
		elif record.id == "chrXIV":
			chrXIV_len = len(record.seq)
		elif record.id == "chrXV":
			chrXV_len = len(record.seq)
		elif record.id == "chrXVI":
			chrXVI_len = len(record.seq)
		elif record.id == "chrmt":
			chrmt_len = len(record.seq)
	print("Reference genome sequence loaded.")

	#サンプルデータ向け配列群をゼロで初期化する
	array_sample_chrI = np.zeros(chrI_len)
	array_sample_chrII = np.zeros(chrII_len)
	array_sample_chrIII = np.zeros(chrIII_len)
	array_sample_chrIV = np.zeros(chrIV_len)
	array_sample_chrV = np.zeros(chrV_len)
	array_sample_chrVI = np.zeros(chrVI_len)
	array_sample_chrVII = np.zeros(chrVII_len)
	array_sample_chrVIII = np.zeros(chrVIII_len)
	array_sample_chrIX = np.zeros(chrIX_len)
	array_sample_chrX = np.zeros(chrX_len)
	array_sample_chrXI = np.zeros(chrXI_len)
	array_sample_chrXII = np.zeros(chrXII_len)
	array_sample_chrXIII = np.zeros(chrXIII_len)
	array_sample_chrXIV = np.zeros(chrXIV_len)
	array_sample_chrXV = np.zeros(chrXV_len)
	array_sample_chrXVI = np.zeros(chrXVI_len)
	array_sample_chrmt = np.zeros(chrmt_len)

	#サンプルデータ向け配列群をゼロで初期化する
	array_control_chrI = np.zeros(chrI_len)
	array_control_chrII = np.zeros(chrII_len)
	array_control_chrIII = np.zeros(chrIII_len)
	array_control_chrIV = np.zeros(chrIV_len)
	array_control_chrV = np.zeros(chrV_len)
	array_control_chrVI = np.zeros(chrVI_len)
	array_control_chrVII = np.zeros(chrVII_len)
	array_control_chrVIII = np.zeros(chrVIII_len)
	array_control_chrIX = np.zeros(chrIX_len)
	array_control_chrX = np.zeros(chrX_len)
	array_control_chrXI = np.zeros(chrXI_len)
	array_control_chrXII = np.zeros(chrXII_len)
	array_control_chrXIII = np.zeros(chrXIII_len)
	array_control_chrXIV = np.zeros(chrXIV_len)
	array_control_chrXV = np.zeros(chrXV_len)
	array_control_chrXVI = np.zeros(chrXVI_len)
	array_control_chrmt = np.zeros(chrmt_len)

	#ratio向け配列群をゼロで初期化する
	array_ratio_chrI = np.zeros(chrI_len)
	array_ratio_chrII = np.zeros(chrII_len)
	array_ratio_chrIII = np.zeros(chrIII_len)
	array_ratio_chrIV = np.zeros(chrIV_len)
	array_ratio_chrV = np.zeros(chrV_len)
	array_ratio_chrVI = np.zeros(chrVI_len)
	array_ratio_chrVII = np.zeros(chrVII_len)
	array_ratio_chrVIII = np.zeros(chrVIII_len)
	array_ratio_chrIX = np.zeros(chrIX_len)
	array_ratio_chrX = np.zeros(chrX_len)
	array_ratio_chrXI = np.zeros(chrXI_len)
	array_ratio_chrXII = np.zeros(chrXII_len)
	array_ratio_chrXIII = np.zeros(chrXIII_len)
	array_ratio_chrXIV = np.zeros(chrXIV_len)
	array_ratio_chrXV = np.zeros(chrXV_len)
	array_ratio_chrXVI = np.zeros(chrXVI_len)
	array_ratio_chrmt = np.zeros(chrmt_len)

	print("Sample data loading started.")
	#サンプルデータを読み込む
	with open(sample_data_name) as sample_file:
		for line in sample_file:
			line = line.split()
			start = line[1]
			end = line[2]
			read_count = line[3]
			if line[0] == "chrI":
				for n in range(int(start), int(end)):
					array_sample_chrI[n] = read_count
			elif line[0] == "chrII":
				for n in range(int(start), int(end)):
					array_sample_chrII[n] = read_count
			elif line[0] == "chrIII":
				for n in range(int(start), int(end)):
					array_sample_chrIII[n] = read_count
			elif line[0] == "chrIV":
				for n in range(int(start), int(end)):
					array_sample_chrIV[n] = read_count
			elif line[0] == "chrV":
				for n in range(int(start), int(end)):
					array_sample_chrV[n] = read_count
			elif line[0] == "chrVI":
				for n in range(int(start), int(end)):
					array_sample_chrVI[n] = read_count
			elif line[0] == "chrVII":
				for n in range(int(start), int(end)):
					array_sample_chrVII[n] = read_count
			elif line[0] == "chrVIII":
				for n in range(int(start), int(end)):
					array_sample_chrVIII[n] = read_count
			elif line[0] == "chrIX":
				for n in range(int(start), int(end)):
					array_sample_chrIX[n] = read_count
			elif line[0] == "chrX":
				for n in range(int(start), int(end)):
					array_sample_chrX[n] = read_count
			elif line[0] == "chrXI":
				for n in range(int(start), int(end)):
					array_sample_chrXI[n] = read_count
			elif line[0] == "chrXII":
				for n in range(int(start), int(end)):
					array_sample_chrXII[n] = read_count
			elif line[0] == "chrXIII":
				for n in range(int(start), int(end)):
					array_sample_chrXIII[n] = read_count
			elif line[0] == "chrXIV":
				for n in range(int(start), int(end)):
					array_sample_chrXIV[n] = read_count
			elif line[0] == "chrXV":
				for n in range(int(start), int(end)):
					array_sample_chrXV[n] = read_count
			elif line[0] == "chrXVI":
				for n in range(int(start), int(end)):
					array_sample_chrXVI[n] = read_count
			elif line[0] == "chrmt":
				for n in range(int(start), int(end)):
					array_sample_chrmt[n] = read_count
	print("Sample data load finished.")

	print("Control data load started.")
	#コントロールデータを読み込む
	with open(control_data_name) as control_file:
		for line in control_file:
			line = line.split()
			start = line[1]
			end = line[2]
			read_count = line[3]
			if line[0] == "chrI":
				for n in range(int(start), int(end)):
					array_control_chrI[n] = read_count
			elif line[0] == "chrII":
				for n in range(int(start), int(end)):
					array_control_chrII[n] = read_count
			elif line[0] == "chrIII":
				for n in range(int(start), int(end)):
					array_control_chrIII[n] = read_count
			elif line[0] == "chrIV":
				for n in range(int(start), int(end)):
					array_control_chrIV[n] = read_count
			elif line[0] == "chrV":
				for n in range(int(start), int(end)):
					array_control_chrV[n] = read_count
			elif line[0] == "chrVI":
				for n in range(int(start), int(end)):
					array_control_chrVI[n] = read_count
			elif line[0] == "chrVII":
				for n in range(int(start), int(end)):
					array_control_chrVII[n] = read_count
			elif line[0] == "chrVIII":
				for n in range(int(start), int(end)):
					array_control_chrVIII[n] = read_count
			elif line[0] == "chrIX":
				for n in range(int(start), int(end)):
					array_control_chrIX[n] = read_count
			elif line[0] == "chrX":
				for n in range(int(start), int(end)):
					array_control_chrX[n] = read_count
			elif line[0] == "chrXI":
				for n in range(int(start), int(end)):
					array_control_chrXI[n] = read_count
			elif line[0] == "chrXII":
				for n in range(int(start), int(end)):
					array_control_chrXII[n] = read_count
			elif line[0] == "chrXIII":
				for n in range(int(start), int(end)):
					array_control_chrXIII[n] = read_count
			elif line[0] == "chrXIV":
				for n in range(int(start), int(end)):
					array_control_chrXIV[n] = read_count
			elif line[0] == "chrXV":
				for n in range(int(start), int(end)):
					array_control_chrXV[n] = read_count
			elif line[0] == "chrXVI":
				for n in range(int(start), int(end)):
					array_control_chrXVI[n] = read_count
			elif line[0] == "chrmt":
				for n in range(int(start), int(end)):
					array_control_chrmt[n] = read_count
	print("Control data loading finished.")

	print("Total read count calculation started.")
	#染色体に由来するリードカウントの総数（つまりミトコンドリアは除く）を求める
	sample_chromosomal_total_read_count = 0
	sample_chromosomal_total_read_count += np.sum(array_sample_chrI)
	sample_chromosomal_total_read_count += np.sum(array_sample_chrII)
	sample_chromosomal_total_read_count += np.sum(array_sample_chrIII)
	sample_chromosomal_total_read_count += np.sum(array_sample_chrIV)
	sample_chromosomal_total_read_count += np.sum(array_sample_chrV)
	sample_chromosomal_total_read_count += np.sum(array_sample_chrVI)
	sample_chromosomal_total_read_count += np.sum(array_sample_chrVII)
	sample_chromosomal_total_read_count += np.sum(array_sample_chrVIII)
	sample_chromosomal_total_read_count += np.sum(array_sample_chrIX)
	sample_chromosomal_total_read_count += np.sum(array_sample_chrX)
	sample_chromosomal_total_read_count += np.sum(array_sample_chrXI)
	sample_chromosomal_total_read_count += np.sum(array_sample_chrXII)
	sample_chromosomal_total_read_count += np.sum(array_sample_chrXIII)
	sample_chromosomal_total_read_count += np.sum(array_sample_chrXIV)
	sample_chromosomal_total_read_count += np.sum(array_sample_chrXV)
	sample_chromosomal_total_read_count += np.sum(array_sample_chrXVI)

	control_chromosomal_total_read_count = 0
	control_chromosomal_total_read_count += np.sum(array_control_chrI)
	control_chromosomal_total_read_count += np.sum(array_control_chrII)
	control_chromosomal_total_read_count += np.sum(array_control_chrIII)
	control_chromosomal_total_read_count += np.sum(array_control_chrIV)
	control_chromosomal_total_read_count += np.sum(array_control_chrV)
	control_chromosomal_total_read_count += np.sum(array_control_chrVI)
	control_chromosomal_total_read_count += np.sum(array_control_chrVII)
	control_chromosomal_total_read_count += np.sum(array_control_chrVIII)
	control_chromosomal_total_read_count += np.sum(array_control_chrIX)
	control_chromosomal_total_read_count += np.sum(array_control_chrX)
	control_chromosomal_total_read_count += np.sum(array_control_chrXI)
	control_chromosomal_total_read_count += np.sum(array_control_chrXII)
	control_chromosomal_total_read_count += np.sum(array_control_chrXIII)
	control_chromosomal_total_read_count += np.sum(array_control_chrXIV)
	control_chromosomal_total_read_count += np.sum(array_control_chrXV)
	control_chromosomal_total_read_count += np.sum(array_control_chrXVI)
	print("Total read count calculation finished.")

	#リードカウント総数の比を求める
	total_read_count_ratio = control_chromosomal_total_read_count / sample_chromosomal_total_read_count

	print("Ratio calculation started.")
	#各座標についてサンプル/コントロール比を求め、リードカウント総数でノーマライズする
	for n in range(0, chrI_len):
		if array_control_chrI[n] == 0:
			array_ratio_chrI[n] = 0
		else:
			array_ratio_chrI[n] = array_sample_chrI[n] / array_control_chrI[n] * total_read_count_ratio
	for n in range(0, chrII_len):
		if array_control_chrII[n] == 0:
			array_ratio_chrII[n] = 0
		else:
			array_ratio_chrII[n] = array_sample_chrII[n] / array_control_chrII[n] * total_read_count_ratio
	for n in range(0, chrIII_len):
		if array_control_chrIII[n] == 0:
			array_ratio_chrIII[n] = 0
		else:
			array_ratio_chrIII[n] = array_sample_chrIII[n] / array_control_chrIII[n] * total_read_count_ratio
	for n in range(0, chrIV_len):
		if array_control_chrIV[n] == 0:
			array_ratio_chrIV[n] = 0
		else:
			array_ratio_chrIV[n] = array_sample_chrIV[n] / array_control_chrIV[n] * total_read_count_ratio
	for n in range(0, chrV_len):
		if array_control_chrV[n] == 0:
			array_ratio_chrV[n] = 0
		else:
			array_ratio_chrV[n] = array_sample_chrV[n] / array_control_chrV[n] * total_read_count_ratio
	for n in range(0, chrVI_len):
		if array_control_chrVI[n] == 0:
			array_ratio_chrVI[n] = 0
		else:
			array_ratio_chrVI[n] = array_sample_chrVI[n] / array_control_chrVI[n] * total_read_count_ratio
	for n in range(0, chrVII_len):
		if array_control_chrVII[n] == 0:
			array_ratio_chrVII[n] = 0
		else:
			array_ratio_chrVII[n] = array_sample_chrVII[n] / array_control_chrVII[n] * total_read_count_ratio
	for n in range(0, chrVIII_len):
		if array_control_chrVIII[n] == 0:
			array_ratio_chrVIII[n] = 0
		else:
			array_ratio_chrVIII[n] = array_sample_chrVIII[n] / array_control_chrVIII[n] * total_read_count_ratio
	for n in range(0, chrIX_len):
		if array_control_chrIX[n] == 0:
			array_ratio_chrIX[n] = 0
		else:
			array_ratio_chrIX[n] = array_sample_chrIX[n] / array_control_chrIX[n] * total_read_count_ratio
	for n in range(0, chrX_len):
		if array_control_chrX[n] == 0:
			array_ratio_chrX[n] = 0
		else:
			array_ratio_chrX[n] = array_sample_chrX[n] / array_control_chrX[n] * total_read_count_ratio
	for n in range(0, chrXI_len):
		if array_control_chrXI[n] == 0:
			array_ratio_chrXI[n] = 0
		else:
			array_ratio_chrXI[n] = array_sample_chrXI[n] / array_control_chrXI[n] * total_read_count_ratio
	for n in range(0, chrXII_len):
		if array_control_chrXII[n] == 0:
			array_ratio_chrXII[n] = 0
		else:
			array_ratio_chrXII[n] = array_sample_chrXII[n] / array_control_chrXII[n] * total_read_count_ratio
	for n in range(0, chrXIII_len):
		if array_control_chrXIII[n] == 0:
			array_ratio_chrXIII[n] = 0
		else:
			array_ratio_chrXIII[n] = array_sample_chrXIII[n] / array_control_chrXIII[n] * total_read_count_ratio
	for n in range(0, chrXIV_len):
		if array_control_chrXIV[n] == 0:
			array_ratio_chrXIV[n] = 0
		else:
			array_ratio_chrXIV[n] = array_sample_chrXIV[n] / array_control_chrXIV[n] * total_read_count_ratio
	for n in range(0, chrXV_len):
		if array_control_chrXV[n] == 0:
			array_ratio_chrXV[n] = 0
		else:
			array_ratio_chrXV[n] = array_sample_chrXV[n] / array_control_chrXV[n] * total_read_count_ratio
	for n in range(0, chrXVI_len):
		if array_control_chrXVI[n] == 0:
			array_ratio_chrXVI[n] = 0
		else:
			array_ratio_chrXVI[n] = array_sample_chrXVI[n] / array_control_chrXVI[n] * total_read_count_ratio
	print("Ratio calculation finished.")

	print("Output started.")
	#計算結果を出力ファイルに書き込む
	output = open(output_name, 'a')
	for n in range(0, chrI_len):
		output.write("chrI" + "\t" + str(n) + "\t" + str(n+1) + "\t" + str(array_ratio_chrI[n]) + "\n")
	for n in range(0, chrII_len):
		output.write("chrII" + "\t" + str(n) + "\t" + str(n+1) + "\t" + str(array_ratio_chrII[n]) + "\n")
	for n in range(0, chrIII_len):
		output.write("chrIII" + "\t" + str(n) + "\t" + str(n+1) + "\t" + str(array_ratio_chrIII[n]) + "\n")
	for n in range(0, chrIV_len):
		output.write("chrIV" + "\t" + str(n) + "\t" + str(n+1) + "\t" + str(array_ratio_chrIV[n]) + "\n")
	for n in range(0, chrV_len):
		output.write("chrV" + "\t" + str(n) + "\t" + str(n+1) + "\t" + str(array_ratio_chrV[n]) + "\n")
	for n in range(0, chrVI_len):
		output.write("chrVI" + "\t" + str(n) + "\t" + str(n+1) + "\t" + str(array_ratio_chrVI[n]) + "\n")
	for n in range(0, chrVII_len):
		output.write("chrVII" + "\t" + str(n) + "\t" + str(n+1) + "\t" + str(array_ratio_chrVII[n]) + "\n")
	for n in range(0, chrVIII_len):
		output.write("chrVIII" + "\t" + str(n) + "\t" + str(n+1) + "\t" + str(array_ratio_chrVIII[n]) + "\n")
	for n in range(0, chrIX_len):
		output.write("chrIX" + "\t" + str(n) + "\t" + str(n+1) + "\t" + str(array_ratio_chrIX[n]) + "\n")
	for n in range(0, chrX_len):
		output.write("chrX" + "\t" + str(n) + "\t" + str(n+1) + "\t" + str(array_ratio_chrX[n]) + "\n")
	for n in range(0, chrXI_len):
		output.write("chrXI" + "\t" + str(n) + "\t" + str(n+1) + "\t" + str(array_ratio_chrXI[n]) + "\n")
	for n in range(0, chrXII_len):
		output.write("chrXII" + "\t" + str(n) + "\t" + str(n+1) + "\t" + str(array_ratio_chrXII[n]) + "\n")
	for n in range(0, chrXIII_len):
		output.write("chrXIII" + "\t" + str(n) + "\t" + str(n+1) + "\t" + str(array_ratio_chrXIII[n]) + "\n")
	for n in range(0, chrXIV_len):
		output.write("chrXIV" + "\t" + str(n) + "\t" + str(n+1) + "\t" + str(array_ratio_chrXIV[n]) + "\n")
	for n in range(0, chrXV_len):
		output.write("chrXV" + "\t" + str(n) + "\t" + str(n+1) + "\t" + str(array_ratio_chrXV[n]) + "\n")
	for n in range(0, chrXVI_len):
		output.write("chrXVI" + "\t" + str(n) + "\t" + str(n+1) + "\t" + str(array_ratio_chrXVI[n]) + "\n")
	print("Output finished.")
