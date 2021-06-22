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
