# -*- coding: utf-8 -*-
from sympy.utilities.iterables import multiset_permutations
import numpy as np
import configparser
import sys,os
import itertools
import math
import subprocess
# FALTAN
# ELiminar Prints incorrectos, agregar correctos.
# mas opciones de permutacion.
# Eliminar los ghost en la construccion de los archivos .com
##### ZONA DE MANEJO VARIALBES DE ENTRADA
def Var_init():
	global enlace, altura, shape
	global Eq_Global
	global Big_variable, totalVertices,atomos_dados
	Big_variable = {}
	enlace = 0.3
	altura = 1
	atomos_dados = 0
def is_number(d,n):
	is_number = True
	try:
		num = float(n)
		# check for "nan" floats
		is_number = num == num   # or use `math.isnan(num)`
		if 'Pcent' in d and (num<0 or num>1):
			print("Valor de ",d, " fuera de los limites [0,1]\n")
			exit(1)
	except ValueError:
		is_number = False
		print ("Informacion en archivo de variables incorrecta\n",d ,"=", n)
		exit(1)
	return is_number
def establecerVariablesDefault():
	global enlace,altura,shape,Eq_Global,totalVertices,atomos_dados
	print(Big_variable)

	# Verificar datos de forma son correctos
	shape = Big_variable["shape"].split(",")
	for ring in shape:
		is_number("Ring-"+ring,ring)
	shape = list(map(int, shape))		#EN este punto los datos son correctos (no son letras)
	totalVertices = sum(map(abs, shape))
	print("TOtal de puntos:", totalVertices)
	
	#Ecuacion Global
	tmp=[]
	#atomos_dados=0
	data= Big_variable["chemical_formula"].split(' ')
	data[:] = [x for x in data if x]   #elimina espacios vacios.
	for i in range(int(len(data)/2)):
		tmp.append(data[i*2].split()*int(data[(i*2)+1]))
		atomos_dados+=abs(int(data[(i*2)+1]))
		#np.append(Eq_Global,data[i*2].split()*int(data[(i*2)+1]))
	if(atomos_dados>totalVertices):
		print("Error, numero de atomos superior a numero de ubicaciones",atomos_dados)
		exit(1)
	else:
		ghost = totalVertices-atomos_dados
		print("Ageregando fantasmas",ghost)
		tmp.append('X'*ghost)
		tmp = list(itertools.chain(*tmp))
	Eq_Global = np.array(tmp)
	print("EQ",Eq_Global)

	# Variables Opcionales
	if "enlace" in Big_variable.keys():
		is_number("enlace",Big_variable["enlace"])
		enlace = float(Big_variable["enlace"])
	if "altura" in Big_variable.keys():
		is_number("Altura",Big_variable["altura"])
		altura = float(Big_variable["altura"])
# Funcion simple lectura de archivo
def leerArchivoParametros(configs):
	config = configparser.ConfigParser()
	config.read(configs)

	for secciones in config.sections():
		for (variable, valor) in config.items(secciones):
			#agregarVariableGlobal(variable, valor)
			Big_variable[variable]=valor
	establecerVariablesDefault()
###################################
# Cambiar funcion para guardar puntos en array modo impresion automata.
# agregar cambio de distancia en Z.
# Para el metodo eclipsado simplemente dividir math.pi por el MISMO numero de lados
# Intentar descubrir como especificar distancia de enlace, no enlace
# revisasr si puntos obtenidos son en orden o no. Segun Chemcraft si.
def CreateRing(sides, bondsize=1, rotation=0, translation=None):
	if sides<0:
		sides=abs(sides)
		rotation=1
	one_segment = math.pi * 2 / sides
	#print(one_segment)
	#print(bondsize)
	if(sides!=1):
	#	print("NO ES 0")
		if(sides == 2):
			radius = bondsize/2
		else:
			radius = math.sin((math.pi-one_segment)/ 2)*bondsize / math.sin(one_segment)
	else:
		radius=0

	if rotation!=0:
		rotation = math.pi/sides
	#print("RADIO:",radius)
	points = [
		(math.sin(one_segment * i + rotation) * radius,
		 math.cos(one_segment * i + rotation) * radius)
		for i in range(sides)]

	if translation:
		points = [[sum(pair) for pair in zip(point, translation)]
				  for point in points]
	return points

def Create3DPolygon(floors, height):
	aux =0
	polygon=[]
	for F in floors:
		tmp=CreateRing(F,enlace) #(abs(F)-1)*
		tmp.append(float(aux))
		aux+=height
		#print(tmp)
		#print("ESPECIFICO: ",tmp[0],"\nAltura",tmp[-1])
		polygon.append(tmp)
	return polygon
#PEMUTACIONES

def GlobalPermutation():
	permutation =[]
	for p in multiset_permutations(Eq_Global):
		print(p)
		permutation.append(p)
		#break
	return permutation

##########
# IMPRESION
def escribirArchivoXYZ(name, title,atomicSymbols ,coordsList):
	input = open(name+".xyz","a")
	#individual = open(title+".xyz","w+")
	input.write(str(atomos_dados)+"\n")
	input.write(title+"\n")
	#individual.write(str(atomos_dados)+"\n")
	#individual.write(title+"\n")
	#print (coordsList)
	#print(shape)
	posicion=0
	for puntos in range(len(shape)):#range(totalVertices):
		for depto in range(abs(shape[puntos])):
			input.write(atomicSymbols[posicion]+"\t"+str(coordsList[puntos][depto][0])+"\t"+str(coordsList[puntos][depto][1])+"\t"+str(coordsList[puntos][-1])+"\n")
			#if (atomicSymbols[posicion] != "X"):
				#individual.write(atomicSymbols[posicion]+"\t"+str(coordsList[puntos][depto][0])+"\t"+str(coordsList[puntos][depto][1])+"\t"+str(coordsList[puntos][-1])+"\n")
			#	continue
			posicion+=1			

def escribirInputGaussian(name,number,atomicSymbols ,coordsList):
	#coordsList = np.array([[row[0],row[1], row[2], row[3]] for row in origincoords])
	input = open(name+str(number)+".com","w+")
	
	#print (var.Big_variable)
	input.write("%NProc="+Big_variable["core"]+"\n")
	input.write("%mem="+Big_variable["memory"]+"GB\n")
	input.write("#"+Big_variable["header"]+"\n\n")
	input.write("Automatic Input "+name+" "+str(number)+"\n\n")
	input.write(Big_variable["charge_multi"]+"\n")
	#for line in coordsList:
	#	input.write(' '.join(map(str, line))+"\n")
	posicion=0
	for puntos in range(len(shape)):#range(totalVertices):
		for depto in range(abs(shape[puntos])):
			if (atomicSymbols[posicion] != "X"):
				input.write(atomicSymbols[posicion]+"\t"+str(coordsList[puntos][depto][0])+"\t"+str(coordsList[puntos][depto][1])+"\t"+str(coordsList[puntos][-1])+"\n")
			posicion+=1			
	input.write("\n")       #Porque Gaussian es espcial
	input.close()
#MAIN
#os.remove("Recopilacion.xyz")
Var_init()
leerArchivoParametros(sys.argv[1])

poly = Create3DPolygon(shape,altura)

#a = np.array([0, 1, 0, 2])
Permu = GlobalPermutation()
#IMPRESION
print (poly)
for resultado in range(len(Permu)):
	escribirArchivoXYZ("Recopilacion","Permu-"+str(resultado),Permu[resultado],poly)
	escribirInputGaussian("Permut",resultado,Permu[resultado],poly)

#FIN IMPRESION
#TEST RMSD por Kabsch 
exit(1)
Eq_Global = list(dict.fromkeys(Eq_Global))   #Eliminar repetidos, util para creacion deAtomPosition
print(Eq_Global)
AtomPosition ={k: v for v, k in enumerate(Eq_Global)}    #Conversion de ecuacion en hash con atomos x ubicacion
print(AtomPosition)

fingerprint=np.empty((len(Permu),len(shape),len(Eq_Global)))  # Numpy, crea matriz de todos los resultados
for resultado in range(len(Permu)):						# For de cada sistema
	aux = 0
	forma = np.zeros((len(shape),len(Eq_Global)))		#Inicio valen 0 cada contador de atomo x anillo
	for puntos in range(len(shape)):					# for de cada anillo
		#print("NIVEL",puntos)
		for cont in range(0,abs(shape[puntos])):		#for de cada vertice del anillo
		#	print(Permu[resultado][aux])
		#	print("rellenando",puntos," ",AtomPosition[Permu[resultado][aux]])
			forma[puntos][AtomPosition[Permu[resultado][aux]]]+=1		#suma 1 al valor del atomo en un anillo.
		#	print(forma)
			aux+=1
	fingerprint[resultado]=np.sort(forma,axis=0)			#cada fingerptinr esta ORDENADO
	#break
print("TERMINO FINGER")
blocked_indexs= []
final_results= []
flag = 0
for finger in range(len(fingerprint)):		#for de cada sistema, cada fingerprint
	if finger not in blocked_indexs:
		duplicates_prior = np.array(finger)
		#print (duplicates_prior)
		flag = 1
	for finger2 in range(finger+1,len(Permu)):		#segundo fingerpritn para evaluar
		#print(fingerprint[finger],"vs",fingerprint[finger2])
		if finger2 not in blocked_indexs:			#Solo evaluar aquellos que no han salido iguales aun.
			equal = (fingerprint[finger]==fingerprint[finger2]).all()	#compara cada fingerprint contra otro
		#	print(finger, "=",finger2,equal)
			if(equal):
				duplicates_prior = np.append(duplicates_prior,finger2)
				blocked_indexs.append(finger2)
		#		print(blocked_indexs)
	#print(duplicates_prior[1])
	if flag == 1:
		final_results.append(duplicates_prior)
		flag = 0
print(final_results)
print("TOTAL HAY", len(final_results), "DISTINTOS CONFIGURACIONES")

# ACA
# Entonces el siguiente paso es ocupar el metodo de RMSD SOLO con cada uno de los array en final_results, pues estos ya son iguales entre ellos.
# Luego hay que construir arrays de igualdad por cada rmsd, y ver cual se adapta mejor al array original de reparto
#exit(1)
RMSD=[]
for bloques in final_results:
	print("BLOQUE DE IGUALDAD DE: ",bloques)
	for template in range(len(bloques)):
		for chimera in range(template+1,len(bloques)):
			print(template,chimera)
			print("Permu-"+str(bloques[template])+" vs "+" Permu-"+str(bloques[chimera]))
			#rmsd = os.system("calculate_rmsd --reorder "+"Permu-"+str(template)+".xyz "+" Permu-"+str(chimera)+".xyz")
			rmsd = subprocess.Popen(["calculate_rmsd","--reorder","Permu-"+str(bloques[template])+".xyz","Permu-"+str(bloques[chimera])+".xyz"], stdout=subprocess.PIPE)
			output = float(rmsd.communicate()[0])
			print(output)
	#		if rmsd not in RMSD:
				#print("PASO")
			output = float('%.5f'%(output))
			if (output < 0.0001):
				RMSD.append(bloques[chimera])
			#RMSD.append('%.5f'%(output))
		#print("New Bacth",template)

#RMSD = list(dict.fromkeys(RMSD))
for template in range(len(RMSD)):
	print (RMSD[template])