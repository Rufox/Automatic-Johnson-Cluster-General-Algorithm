from sympy.utilities.iterables import multiset_permutations
import numpy as np
import configparser
import sys,os
import itertools
import math

# This program have been tested with Python 3.8

##### Setting vars
def Var_init():
	global enlace, altura, shape
	global Eq_Global
	global Big_variable, totalVertices,atomos_dados
	Big_variable = {}
	atomos_dados = 0
def is_number(d,n):
	is_number = True
	try:
		num = float(n)
		# check for "nan" floats
		is_number = num == num   # or use `math.isnan(num)`
		if 'Pcent' in d and (num<0 or num>1):
			print("Value ",d, " out of bounds [0,1]\n")
			exit(1)
	except ValueError:
		is_number = False
		print ("Parameters in Config file are incorrect (number expected) {}={}".format(d,n))
		exit(1)
	return is_number

# This will read Config.in file and set variables
def establecerVariablesDefault():
	global enlace,altura,shape,Eq_Global,totalVertices,atomos_dados
	#print(Big_variable)

	# Verificar datos de forma son correctos
	shape = Big_variable["shape"].split(",")
	for ring in shape:
		is_number("Ring-"+ring,ring)
	shape = list(map(int, shape))		#EN este punto los datos son correctos (no son letras)
	totalVertices = sum(map(abs, shape))
	print("Total Amount of Vertices: ", totalVertices)
	
	#Ecuacion Global
	tmp=[]
	
	data= Big_variable["chemical_formula"].split(' ')
	data[:] = [x for x in data if x]   #Delete void spaces.
	for i in range(int(len(data)/2)):
		tmp.append(data[i*2].split()*int(data[(i*2)+1]))
		atomos_dados+=abs(int(data[(i*2)+1]))
	if(atomos_dados>totalVertices):
		print("Error.\nNumber of atoms is larger than the possible positions\n",atomos_dados)
		exit(1)
	else:
		ghost = totalVertices-atomos_dados
		print("Adding Ghosts Positions",ghost)
		tmp.append('X'*ghost)
		tmp = list(itertools.chain(*tmp))
	Eq_Global = np.array(tmp)
	print("Chemical Equation: ",Eq_Global)

	if "distances" in Big_variable.keys():
		enlace = Big_variable["distances"].split(",")
		for distance_array in enlace:
			is_number("Distance-"+distance_array, distance_array)
		enlace = list(map(float, enlace))
	else:
		print("Error\nNo \"distance\" parameter given\n")
		exit(2)


	if "height" in Big_variable.keys():
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
# To rotate a ring just divide math.pi by the number of sides -> angle of rotation
def CreateRing(sides, bondsize=1, rotation=0, translation=None):
	if sides<0:					# Rotate ring 
		sides=abs(sides)
		rotation=1
	one_segment = math.pi * 2 / sides
	if(sides!=1 and sides!=2):	#Specials cases: (1) a simple point and (2) a line
		radius = math.sin((math.pi-one_segment)/ 2)*bondsize / math.sin(one_segment)
	elif (sides==2):			# line
		radius=bondsize/2
	else:						#point
		radius=0

	if rotation!=0:
		rotation = math.pi/sides
	# Rotations and Translation operators
	points = [
		(math.sin(one_segment * i + rotation) * radius,
		 math.cos(one_segment * i + rotation) * radius)
		for i in range(sides)]
	if translation:
		points = [[sum(pair) for pair in zip(point, translation)]
				  for point in points]
	return points

def Create3DPolygon(floors):
	aux = 0
	flag = 0
	polygon=[]
	# Creation of each ring 
	for F in floors:
		tmp=CreateRing(F,enlace[flag]) #(abs(F)-1)*
		tmp.append(float(aux))
		try:
			aux+=altura[flag]
		except IndexError as e:	#Only last ring with error, search for ring+1 that doesn't exist
			pass
		flag+=1
		polygon.append(tmp)
	return polygon 		#An array Ring1[[x1,y1],[x2,y2],z_ring1],Ring2[[x1,y1],z_ring2]...
#PEMUTACIONES

# Function will iterate all possible combinations of atoms in positions.
# WARNING: users should know how many iteration can be, number can and will escalate exponentially with number of positions
def GlobalPermutation():
	permutation =[]
	for p in multiset_permutations(Eq_Global):	#uses sympy method
		print(p)
		permutation.append(p)
	return permutation

##########
# PRINTING OF DATA
def escribirArchivoXYZ(name, title,atomicSymbols ,coordsList):
	input = open(name+".xyz","a")
	input.write(str(atomos_dados)+"\n")
	input.write(title+"\n")
	posicion=0
	for puntos in range(len(shape)):#range(totalVertices):
		for depto in range(abs(shape[puntos])):
			input.write(atomicSymbols[posicion]+"\t"+str(coordsList[puntos][depto][0])+"\t"+str(coordsList[puntos][depto][1])+"\t"+str(coordsList[puntos][-1])+"\n")
			posicion+=1			

# Change or add any information for an Gaussian input in this method
def escribirInputGaussian(name,number,atomicSymbols ,coordsList):
	input = open(name+str(number)+".com","w+")
	
	input.write("%NProc="+Big_variable["core"]+"\n")
	input.write("%mem="+Big_variable["memory"]+"GB\n")
	input.write("#"+Big_variable["header"]+"\n\n")
	input.write("Automatic Input "+name+" "+str(number)+"\n\n")
	input.write(Big_variable["charge_multi"]+"\n")
	posicion=0
	for puntos in range(len(shape)):#range(totalVertices):
		for depto in range(abs(shape[puntos])):
			input.write(atomicSymbols[posicion]+"\t"+str(coordsList[puntos][depto][0])+"\t"+str(coordsList[puntos][depto][1])+"\t"+str(coordsList[puntos][-1])+"\n")
			posicion+=1			
	input.write("\n")       #Porque Gaussian es espcial
	input.close()

#####################################################################
#MAIN
Var_init()
leerArchivoParametros(sys.argv[1])

# Delete previous iterations
try:
    os.remove("Results_Combinations.xyz")
except OSError:
    pass

# Create rings
poly = Create3DPolygon(shape)

# Poblate rings
Permu = GlobalPermutation()

# Print information
for resultado in range(len(Permu)):
	escribirArchivoXYZ("Results_Combinations","Permu-"+str(resultado),Permu[resultado],poly)
	escribirInputGaussian("Permut",resultado,Permu[resultado],poly)

try:
	os.mkdir("InputsGaussian")
except FileExistsError:
	pass

os.system('mv *com InputsGaussian/')  
