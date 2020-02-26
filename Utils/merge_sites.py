# -*- coding: utf-8 -*-
# File Name: merge_sites.py
# 1Author: Zhanghena
# Mail: zhanghena1206@gmail.com
# Created time: 2018-02-28 10:09:27
#!/usr/bin/env python
#****************************************
import os, sys, time
log = open(sys.argv[0] + '.arg', 'a+')
log.write(str(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))) + ' dir: ' + os.getcwd()  + ' command: ' + ' '.join(sys.argv) + '\n')
log.close()
#args: output, inputfilenumber, input1, input2

#input1:
#SeqID			refBase strand	cov		C		T		C/C+T	C+T/cov	A 		C		G		N		T		total
#cluster1:13     C       +       9       0       9       0       1       0       0       0       0       9       17

#output: SeqID  refBase	refStrand	input1_cov	input1_C	input1_T	input1_methRate	input1_C+T/cov input2_cov input2_C input2_T input2_methRate input2_C+T/cov....

##读取文件
d = {}
input = locals() ## 因为Python的变量名就是一个字典的key而已。要获取这个字典，直接用locals和globals函数即可。

for i in range (3,3 + int(sys.argv[2])):
	input['d%s' % i] = {} # 动态生成变量名
	with open(sys.argv[i],'r') as f:
		print 'Reading file ' + sys.argv[i] + '...'
		for line in f:
			line = line.strip().split('\t')
			key = line[0] + ':' + line[1] + ':' + line[2]
			d[key] = line[2]
			input['d%s' % i][key] = line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[6] + '\t' + line[7]
		f.close()
print 'Finished the reading.'

##写出结果
print 'Writing the results...'
with open(sys.argv[1],'w') as out:
	for key in d:
		line = key 
		for i in  range (3,3 + int(sys.argv[2])):
			if input['d%s' % i].has_key(key):
				line = line + '\t' + input['d%s' % i][key]
			else:
				line = line + '\t0\t0\t0\t0\t0'
		out.write(line + '\n')
out.close
print 'Finished the writing.'
