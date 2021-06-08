import argparse
import time
import os
import numpy as np

start = time.ctime()

aminoacid_code = {'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','CYS':'C','MET':'M',
					'PHE':'F','TYR':'Y','TRP':'W','PRO':'P','SER':'S','THR':'T','ASN':'N',
					'GLN':'Q','ASP':'D','GLU':'E','HIS':'H','LYS':'K','ARG':'R'}

# Основная функция для извлечения b-factor из структуры
def extract_bf(pdb_file):
	dict_bf = {}

	with open(pdb_file,'r') as file:
		data = file.read()

	for stroka in data.strip('\n').split('\n'):
		if 'ATOM' == stroka.split(' ')[0]:
			clear_stroka = clear_blank(stroka.strip('\n').split(' '))

			possible_aa = clear_stroka[3]
			possible_chain = clear_stroka[4]
			possible_position = clear_stroka[5]
			bfac = clear_stroka[-2]

            #one_check - первый этап проверки: на корректность трехбуквенного названия амк
			aa = aux_one(possible_aa)
			if aa == '_':
				possible_aa = clear_stroka[2]
				possible_chain = clear_stroka[3]
				possible_position = clear_stroka[4]

				aa = aux_one(possible_aa)

			#two_check - второй этап проверки: на корректность названия цепи и позиции амк
			possible_chain_2,possible_position_2 = aux_two(possible_chain,possible_position)

			#three_check - третий этап проверки: на корректность
			chain,position = aux_three(possible_chain_2,possible_position_2)

			if aa != None:
				#print(aminoacid_code[aa]+'_'+chain+'_'+position)
				dict_bf.setdefault(aminoacid_code[aa]+'_'+chain+'_'+position,[])
				#four_check - четвёртый этап проверки: на корректность значения Bfac
				if aux_four(bfac) not in dict_bf[aminoacid_code[aa]+'_'+chain+'_'+position]:
					dict_bf[aminoacid_code[aa]+'_'+chain+'_'+position].append(float(aux_four(bfac)))

	#print(dict_bf)
	mean_dict_bf = aux_mean(dict_bf)
	#print(mean_dict_bf)

	list_bf = []
	for key,value in mean_dict_bf.items():
		#print(key,value)
		list_bf.append(key+'\t'+value)

	with open(pdb_file.split('.')[0]+'_bf.txt','w') as file:
		file.write('\n'.join(list_bf))

# По сути строки ATOM могут содержать пропуски, поэтому данные пропуски очищаются и получаемые
# данной функцией строки (в виде списка) не содержат пропусков. Однако в результате очистки
# от пропусков строки в итоге могут быть разной длины, поэтому необходимы дальнейшие проверки
# на правильность извлекаемых данных
def clear_blank(list):
	new_list = []
	for element in list:
		if element != '':
			new_list.append(element)
	return new_list

# Первая функция проверки: проверка на правильность названия аминокислоты.
# Сначала проверяется длина названия, как правило, она должна равняться трём. Но бывают
# ситуации, например, такие: BGLU, AGLU (см. 1JDB.pdb 429 B). Тогда записывается лишь
# трехбуквенное название амк.
def aux_one(possible_aa):
	aa_list = []
	if len(possible_aa) >= 3:
		for AA in aminoacid_code.keys():
			if AA in possible_aa:
				aa_list.append(AA)
		if len(aa_list) == 1:
			aa = aa_list[0]
		else:
			aa = '_'
	else:
		aa = None
	return aa

# Вторая функция проверки: цепь и позиция могут сливаться в одну строку, что нужно обязательно
# исправлять для корректного извлечения данных (см.1JDB B 1026). Обычно цепь состоит из
# одного символа, стоящего перед цифровой позицией амк.
def aux_two(possible_chain,possible_position):
	if len(possible_chain) > 1:
		chain = possible_chain[0]
		position = possible_chain[1:]
	else:
		chain = possible_chain
		position = possible_position
	return chain,position

# Третья функция проверки: 	иногда позиция может быть отрицательной, т.е. перед ней будет стоять минус,
# (см. 1WAR начало), или, наоборот, позиция может содержать букву после себя (см. 1PBH начало).
# Такое необходимо тоже исправлять: было -> стало: А 1А -> A+A_1 или A -6 -> A+-_6
def aux_three(possible_chain,possible_position):
	if possible_position.isdigit() == False:
		if possible_position[0] == '-':
			new_position = possible_position[1:]
			new_chain = possible_chain+'+'+possible_position[0]
		else:
			new_position = possible_position[:-1]
			new_chain = possible_chain+'+'+possible_position[-1]
	else:
		new_chain,new_position = possible_chain,possible_position
	return new_chain,new_position

# Четвёртая функция проверки: иногда значение b-factor может сливаться с соседним левым значением и
# тогда число будет содержать две точки (см. 1WAR A 76), такое тоже нужно исправлять.
def aux_four(b_factor):
	if b_factor.count('.') == 2:
		return b_factor[4:]
	else:
		return b_factor

# Вспомогательная функция усреднения значений b-factor для амк по атомам, составляющим её.
def aux_mean(dict):
	new_dict = {}
	for key in dict.keys():
		new_dict[key] = str(round(np.mean(dict[key]),2))
	return new_dict

def main():
	parser = argparse.ArgumentParser(description='Selecting of structural features')
	parser.add_argument('input',help='Enter the name of PDB file')
	args = parser.parse_args()

	extract_bf(args.input)

if __name__ == '__main__':
	main()
	end = time.ctime()

	print(f'\n******** {__file__} ********\nStart: {start}\nFinish: {end}\n')
