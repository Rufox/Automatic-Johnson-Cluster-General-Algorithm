from sympy.utilities.iterables import multiset_permutations
import numpy as np
import configparser
import sys,os
import itertools
import math
#import shutil

##### ZONA DE MANEJO VARIALBES DE ENTRADA
def Var_init():
	global enlace, altura, shape
	global Eq_Global
	global Big_variable, totalVertices,atomos_dados
	Big_variable = {}
	#enlace = 0.3
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
		print ("Parameters in Config file are incorrect (number expected) {}={}".format(d,n))
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
	print("Total Amount of Vertices: ", totalVertices)
	
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
		print("Error.\nNumber of atoms is larger than the possible positions\n",atomos_dados)
		exit(1)
	else:
		ghost = totalVertices-atomos_dados
		print("Adding Ghosts Positions",ghost)
		tmp.append('X'*ghost)
		tmp = list(itertools.chain(*tmp))
	Eq_Global = np.array(tmp)
	print("EQ",Eq_Global)

	# Variables Opcionales
	if "distances" in Big_variable.keys():
		#is_number("distances",Big_variable["distances"])
		enlace = Big_variable["distances"].split(",")
		for distance_array in enlace:
			is_number("Distance-"+distance_array, distance_array)
		enlace = list(map(float, enlace))
		#print ("OCURRION\n")
	else:
		print("Error\nNo \"distance\" parameter given\n")
		exit(2)


	if "height" in Big_variable.keys():
		#is_number("height",Big_variable["height"])
		altura = Big_variable["height"].split(",")
		for height_array in altura:
			is_number("Height-"+height_array, height_array)
		altura = list(map(float, altura))
	else:
		print("Error\nNo \"height\" parameter given\n")
		exit(3)

	if(len(shape) > len(enlace)):
		print("Shape and Distances parameters HAVE TO BE of the same length.\nShape lenght ({}) != Distances lenght ({})".format(shape,enlace))
		exit(4)
	if(len(shape) != (len(altura)+1)):
		print("Height parameter HAVE TO BE of lengh of Shape parameters - 1.\nShape lenght ({}) != Height lenght ({})".format(shape,altura))
		exit(5)
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
	#print("bondsize es {} para el floor {}".format(bondsize,sides))
	if sides<0:
		sides=abs(sides)
		rotation=1
	one_segment = math.pi * 2 / sides
	#print("el one_segment es {}".format(one_segment))
	if(sides!=1 and sides!=2):
		radius = math.sin((math.pi-one_segment)/ 2)*bondsize / math.sin(one_segment)
	#	print ("Radio da: {}".format(radius))
	elif (sides==2):
		radius=bondsize/2
	else:
		radius=0

	if rotation!=0:
		rotation = math.pi/sides
	points = [
		(math.sin(one_segment * i + rotation) * radius,
		 math.cos(one_segment * i + rotation) * radius)
		for i in range(sides)]
	#print("antes de trasladar {}".format(points))
	if translation:
		points = [[sum(pair) for pair in zip(point, translation)]
				  for point in points]
	#print(points)
	return points

def Create3DPolygon(floors):
#	global enlace
	aux =0
	flag=0
	polygon=[]
	#print("heigh es {}".format(altura))
	#print ("Crear poligono, enlace es {}".format(enlace))
	for F in floors:
	#	print("Piso con {} puntos".format(F))
		tmp=CreateRing(F,enlace[flag]) #(abs(F)-1)*
		tmp.append(float(aux))
		#print(altura[flag])
		try:
			aux+=altura[flag]
		except IndexError as e:
			pass
		flag+=1
		#print(tmp)
		#print("ESPECIFICO: ",tmp[0],"\nAltura",tmp[-1])
		polygon.append(tmp)
	#print("POLIGON VA ASI")
	#print(polygon)
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
	input.write(str(atomos_dados)+"\n")
	input.write(title+"\n")
	#print (coordsList)
	#print(shape)
	posicion=0
	for puntos in range(len(shape)):#range(totalVertices):
		for depto in range(abs(shape[puntos])):
			input.write(atomicSymbols[posicion]+"\t"+str(coordsList[puntos][depto][0])+"\t"+str(coordsList[puntos][depto][1])+"\t"+str(coordsList[puntos][-1])+"\n")
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
			input.write(atomicSymbols[posicion]+"\t"+str(coordsList[puntos][depto][0])+"\t"+str(coordsList[puntos][depto][1])+"\t"+str(coordsList[puntos][-1])+"\n")
			posicion+=1			
	input.write("\n")       #Porque Gaussian es espcial
	input.close()

#####################################################################
#MAIN
#os.remove("Recopilacion.xyz")
Var_init()
leerArchivoParametros(sys.argv[1])

try:
    os.remove("Results_Combinations.xyz")
except OSError:
    pass

#print("Shape tiene lo siguiente:{} y height tiene {}".format(shape,altura))
poly = Create3DPolygon(shape)

#a = np.array([0, 1, 0, 2])
Permu = GlobalPermutation()
#IMPRESION
#print (poly)


for resultado in range(len(Permu)):
	escribirArchivoXYZ("Results_Combinations","Permu-"+str(resultado),Permu[resultado],poly)
	escribirInputGaussian("Permut",resultado,Permu[resultado],poly)

try:
	os.mkdir("InputsGaussian")
except FileExistsError:
	pass

os.system('mv *com InputsGaussian/')  

#shutil.copy('*com', 'InputsGaussian/')
#FIN IMPRESION