import argparse
import os
import time
import pandas as pd
import numpy as np

start = time.ctime()

#Дополнительные особенности, необходимые для правильного извлечения структурных дескрипторов
SS = ['a','b','c','d','e','f','g','h',
	'i','j','k','l','m','n','o','p','q',
	'r','s','t','u','v','w','z','x','y']
HTR = 'X'
B_sheet = ['B','E']
a_helix = ['H','G','I']
loop = ['T','S']

# Очистка строки от пропусков
def clear_list(list):
	new_list = []
	for element in list:
		if element != '':
			new_list.append(element)
	return new_list

# Функция нахождения ACC в строке файла dssp
def find_ACC(stroka):
	for index,element in enumerate(stroka):
		if ',' in element:
			break
	return stroka[index-1]

# Функция для нахождения концевых петель (dist_loop)
def check_end(data_func):

	for index in range(len(data_func)):
		type_str = data_func[index].split(',')[3]
		if index < len(data_func)/2:
			if type_str == 'B' or type_str == 'H':
				start_index = index
				break
		else:
			start_index = index
			break
	return start_index

# Функция для выяснения длины петель (len_loop)
def find_loops(data):
	dict_len_loop = {}
	dict_temp_loop = {}
	list_temp_loop = []
	len_loop = 0
	for index,stroka in enumerate(data):
		type_SS = stroka.split(',')[3]
		if type_SS == 'O':
			len_loop += 1
			list_temp_loop.append(index)
		else:
			if len(list_temp_loop) != 0:
				dict_temp_loop = dict_temp_loop.fromkeys(list_temp_loop,len_loop)
				dict_len_loop.update(dict_temp_loop)
			else:
				continue
			len_loop = 0
			list_temp_loop = []
			dict_temp_loop = {}
	dict_temp_loop = dict_temp_loop.fromkeys(list_temp_loop,len_loop)
	dict_len_loop.update(dict_temp_loop)
	return dict_len_loop

# Основная функция для формирования конечного датасета
def form_dataset(name_file_dssp):
	input_ID = name_file_dssp.split('.')[0]

	with open(name_file_dssp,'r') as file:
		data = file.read()

    #Начало для считывания данных из файла dssp - строка, содержащая значок "#"
	for index,stroka in enumerate(data.strip('\n').split('\n')):
		if '#' in stroka:
			start_index = index
			break

	list_all_data = []

	'''I. AA, position of AA, chain, secondary structure type, BP1 and BP2, ACC'''
	for stroka in data.strip('\n').split('\n')[start_index+1:]:
		stroka_data = clear_list(stroka.split(' '))
    	#print(stroka_data)

        #num_aa,chain,AA,STR
		num_aa = stroka_data[1]
		chain = stroka_data[2]
		aa = stroka_data[3]
		secondary_structure = stroka_data[4]

		if aa in SS:
			aa = 'C'

		if aa == HTR:
			continue

		if num_aa.isdigit() != True:
			if '-' in num_aa:
				num_aa = num_aa
			elif '!' in num_aa:
				continue
			else:
				chain = num_aa[-1]
				num_aa = num_aa[:-1]
				aa = stroka_data[2]
				secondary_structure = stroka_data[3]
		'''print(num_aa+'_'+chain+'_'+aa+'_'+secondary_structure)'''

        #ACC
		ACC = find_ACC(stroka_data)

        #BP1 and BP2
		bp1 = stroka_data[-6]
		bp2 = stroka_data[-5]

		if bp1 != '0':
			bp1 = '1'

		if bp2 != '0':
			bp2 = '1'

        #Secondary structure type
		if secondary_structure in B_sheet:
			secondary_structure = 'B'
		elif secondary_structure in a_helix:
			secondary_structure = 'H'
		elif secondary_structure in loop:
			secondary_structure = 'O'
		else:
			secondary_structure = 'O'
		'''print(secondary_structure,ACC,bp1,bp2)'''

		list_all_data.append(f'{num_aa},{chain},{aa},{secondary_structure},{bp1},{bp2},{ACC}')
		'''print(list_all_data)
		input()'''

	'''II. DIST_LOOP and LEN_LOOP'''
    #convert list_all_data into dictionary temp_df for transforming in dataframe
	dict_convert = lambda x: {'num_aa':[i.split(',')[0] for i in x],
                              'chain':[i.split(',')[1] for i in x],
                              'AA':[i.split(',')[2] for i in x],
                              'STR':[i.split(',')[3] for i in x],
                              'bp1':[i.split(',')[4] for i in x],
                              'bp2':[i.split(',')[5] for i in x],
                              'ACC':[i.split(',')[6] for i in x],
                              }
	temp_data = dict_convert(list_all_data)
	temp_df = pd.DataFrame(temp_data)
	'''print(temp_df)'''

    #dist_loop and len_loop
	group_var = temp_df.groupby(by='chain')
	header = ','.join(list(temp_df.columns))+',dist_loop,len_loop'

    #list_all_data_2
	list_all_data_2 = []

	for str_ch in group_var:
		df = str_ch[1]
		str_df = df.to_csv(path_or_buf=None,index=False)

        #necessarily to dist_loop
		data = str_df.strip('\n').split('\n')[1:]
		reverse_data = list(reversed(data))

        #length of N- and C-loops
		N_terminal_loop_len = check_end(data)
		C_terminal_loop_len = check_end(reverse_data)
        #AA indexes of N- and C-loops
		N_terminal_loop_index_list = list(range(N_terminal_loop_len))
		C_terminal_loop_index_list = list(range(len(data)-C_terminal_loop_len,len(data)))

        #AA values of N- and C-loops for dist_loop feature: threshold - 10 (1,...,10)
		N_terminal_loop_value_list = list(reversed(range(1,N_terminal_loop_len+1)))
		N_terminal_loop_value_list = [i if i<10 else 10 for i in N_terminal_loop_value_list]
		C_terminal_loop_value_list = list(range(1,C_terminal_loop_len+1))
		C_terminal_loop_value_list = [i if i<10 else 10 for i in C_terminal_loop_value_list]

        #dictionary for separate N- and C-loops{AA index:AA value}
		N_terminal_loop_dict = dict(zip(N_terminal_loop_index_list,N_terminal_loop_value_list))
		C_terminal_loop_dict = dict(zip(C_terminal_loop_index_list,C_terminal_loop_value_list))

        #Common dictionary for dist_loop
		loop_dict = {}
		loop_dict.update(N_terminal_loop_dict)
		loop_dict.update(C_terminal_loop_dict)
		'''
		print(N_terminal_loop_len)
		print(C_terminal_loop_len)
		print(N_terminal_loop_index_list)
		print(C_terminal_loop_index_list)
		print(N_terminal_loop_value_list)
		print(C_terminal_loop_value_list)
		print(len(data))
		print(list(N_terminal_loop_dict.keys()))
		print(C_terminal_loop_dict)
		print(loop_dict)
		input()
		'''
    	#len_loop
		len_loop_dict = find_loops(data)

		for index,stroka in enumerate(data):

        	#dist_loop and len_loop
			if (index in list(loop_dict.keys())) and (index in list(len_loop_dict.keys())):
				new_stroka = stroka.strip('\r')+','+str(loop_dict[index])+','+str(len_loop_dict[index])
			elif (index in list(loop_dict.keys())) and (index not in list(len_loop_dict.keys())):
				new_stroka = stroka.strip('\r')+','+str(loop_dict[index])+',0'
			elif (index not in list(loop_dict.keys())) and (index in list(len_loop_dict.keys())):
				new_stroka = stroka.strip('\r')+',0,'+str(len_loop_dict[index])
			else:
				new_stroka = stroka.strip('\r')+',0,0'

			list_all_data_2.append(new_stroka)

	'''III. B-factor'''
	#add b-factor feature
	dict_convert_2 = lambda x: {'num_aa':[i.split(',')[0] for i in x],
                              'chain':[i.split(',')[1] for i in x],
                              'AA':[i.split(',')[2] for i in x],
                              'STR':[i.split(',')[3] for i in x],
                              'bp1':[i.split(',')[4] for i in x],
                              'bp2':[i.split(',')[5] for i in x],
                              'ACC':[i.split(',')[6] for i in x],
                              'dist_loop':[i.split(',')[7] for i in x],
                              'len_loop':[i.split(',')[8] for i in x],
                              }

	temp_data = dict_convert_2(list_all_data_2)

	temp_df = pd.DataFrame(temp_data)

	#bf_file
	bf_file = pd.read_csv(f'{input_ID}_bf.txt','\t',names=['position','bfac'])
	result_df = pd.concat([temp_df,bf_file['bfac']],axis=1)

	#Check Nan - иногда bfactor есть у аминокислоты, но самой амк в dssp нет
	if True in np.unique(result_df.isna()):
		result_df = result_df.dropna(axis=0)

	#change type of some features
	result_df = result_df.astype({'ACC':'int','dist_loop':'int','len_loop':'int'})

	'''IV. Normalization of data: common and separate approaches'''
    #standartization: ACC_N_com, bfac_N_sep, dist_loop_N_com, len_loop_N_sep
	'''print(result_df)'''

    #separate approach
	result_data = pd.DataFrame()
	group_var_2 = result_df.groupby(by='chain')
	for str_ch_2 in group_var_2:
		temp_df_2 = str_ch_2[1]

		min_bfac_sep = temp_df_2['bfac'].min()
		max_bfac_sep = temp_df_2['bfac'].max()
		min_len_loop_sep = temp_df_2['len_loop'].loc[temp_df_2['STR']=='O'].min()
		max_len_loop_sep = temp_df_2['len_loop'].loc[temp_df_2['STR']=='O'].max()
		#print(min_bfac_sep,max_bfac_sep,min_len_loop_sep,max_len_loop_sep)

		temp_df_2 = temp_df_2.assign(bfac_N_sep = lambda x: np.around((x.bfac-min_bfac_sep)/(max_bfac_sep-min_bfac_sep),3),
                                     len_loop_N_sep = lambda x: np.around((x.len_loop-min_len_loop_sep)/(max_len_loop_sep-min_len_loop_sep),3))

		result_data = result_data.append(temp_df_2)
		#print(result_data)
	result_data['bfac_N_sep'] = result_data['bfac_N_sep'].fillna(0.0)
	result_data['len_loop_N_sep'] = result_data['len_loop_N_sep'].fillna(0.0)
    #common approach
	min_ACC_com = result_data['ACC'].min()
	max_ACC_com = result_data['ACC'].max()
	min_dist_loop_com = result_data['dist_loop'].loc[result_data['STR']=='O'].min()
	max_dist_loop_com = result_data['dist_loop'].loc[result_data['STR']=='O'].max()

    #standartization of ACC and dist_loop
	result_data = result_data.assign(ACC_N_com = lambda x: np.around((x.ACC-min_ACC_com)/(max_ACC_com-min_ACC_com),3),
                                     dist_loop_N_com = lambda x: np.around((x.dist_loop-min_dist_loop_com)/(max_dist_loop_com-min_dist_loop_com),3))
	'''print(result_data)'''

	'''V. One-hot encoding for secondary structure types'''
    #one-hot encoding
	result_data['SS_type'] = result_data['STR']
	#result_data = pd.get_dummies(result_data,columns=['STR'])
	structure_datasets = result_data.groupby(by=['SS_type'])
	for dataset in structure_datasets:
		if dataset[0] == 'H':
			test_dataset = dataset[1][['num_aa','chain','AA','SS_type','ACC','ACC_N_com','bfac','bfac_N_sep']]
			test_dataset.to_csv(f'Helix_{input_ID}_test.csv',index=False)
		elif dataset[0] == 'B':
			test_dataset = dataset[1][['num_aa','chain','AA','SS_type','ACC','ACC_N_com','bfac','bfac_N_sep','bp1','bp2']]
			test_dataset.to_csv(f'B-sheet_{input_ID}_test.csv',index=False)
		else:
			test_dataset = dataset[1][['num_aa','chain','AA','SS_type','ACC','ACC_N_com','bfac','bfac_N_sep','dist_loop','dist_loop_N_com','len_loop','len_loop_N_sep']]
			test_dataset.to_csv(f'Loop_{input_ID}_test.csv',index=False)

	'''FINISH'''

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Selecting of structural features')
	parser.add_argument('input',help='Enter the name of DSSP file')
	args = parser.parse_args()

	dssp_file = args.input
	form_dataset(dssp_file)

	end = time.ctime()
	print(f'\n******** {__file__} ********\nStart: {start}\nFinish: {end}\n')
