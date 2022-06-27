# -*- coding:utf-8 -*-
#! /usr/bin/python
from sys import argv
from math import cos, sin, pi
from math import sqrt
from copy import deepcopy

sin_deg=lambda x: sin(x*pi/180)    #перевод из радиан в градусы
cos_deg=lambda x: cos(x*pi/180)


"""
ЧТО ДЕЛОТЬ
Если преводим в хyz файлы 1 и 2 то -f start_file_base.cp2k  -n start_file_add.cp2k -o new_file_coord.xyz
infile=argv.index('-f')
inputfile = argv[infile + 1]
otherf= argv.index('-n')     
otherfile = argv[otherf + 1]  
ofile=argv.index('-o')
outputfile = argv[ofile + 1]

cp2k_object1 = CP2K(filename=inputfile)
cp2k_object2=CP2K(filename=otherfile)   
cp2k_object_sum=cp2k_object1+cp2k_object2  
cp2k_object_sum.droptofile(outputfile)      #файл в cp2k. при вводе в консоли файл вывода *.cp2k
print cp2k_object_sum.xyzcoord()             #вывод конкатенир. файла в  xyz


Если преводим в хyz файл 1 то -f start_file_base.cp2k -o new_file_coord.xyz
infile=argv.index('-f')
inputfile = argv[infile + 1]
ofile=argv.index('-o')
outputfile = argv[ofile + 1]
cp2k_object1 = CP2K(filename=inputfile)
print cp2k_object1.xyzcoord() # оставить если нужен перевод в xyz только для одного файла

Если получить файл в cp2k для  файлов 1 и 2 то -f start_file_base.cp2k  -n start_file_add.cp2k -o new_file_coord.cp2k
infile=argv.index('-f')
inputfile = argv[infile + 1]
otherf= argv.index('-n')     
otherfile = argv[otherf + 1]  
ofile=argv.index('-o')
outputfile = argv[ofile + 1]

cp2k_object1 = CP2K(filename=inputfile)
cp2k_object2=CP2K(filename=otherfile)   
cp2k_object_sum=cp2k_object1+cp2k_object2 
print cp2k_object_sum.droptofile(outputfile)
"""


class CP2K_Atom():
    @staticmethod
    def parseatom(s):       #считывание строки, перевод данных в строковый (название атома) и числовой (x, y, z и др х-тики) формат.
        spasete=filter(lambda znach: znach!= '',s.split())
        return [spasete[0]]+map(float,spasete[1:])

    def __init__(self,line): #числовые данные каждого атома

        if line!=None:
            q=self.parseatom(line)
            if len(q)==4:
                self.name,self.x,self.y,self.z=q
                self.occupancy=1.0
                self.etwas=1
            elif len(q)==5:
                self.name,self.x,self.y,self.z,self.occupancy=q
                self.etwas=1
            elif len(q)==6:
                self.name,self.x,self.y,self.z,self.occupancy,self.etwas=q
        else:
            self.x=0.0
            self.y=0.0
            self.z=0.0
    def __str__(self):
        try:
            return '{0:2s} {1:16.9f} {2:16.9f} {3:16.9f} {4:6.4f} {5:6.4f}\n'.format(self.name, self.x, self.y, self.z,self.occupancy,self.etwas)
        except ValueError:
            print('Данные в строке не верны {0:s}'.format())
            raise SystemExit(1)

class CP2K_Cell():
    @staticmethod
    def parseatom(s):
        pamagite=filter(lambda znach: znach!= '',s.split())
        return [pamagite[0]]+map(float,pamagite[1:])

    def __init__(self,line1,line2):
        if line1!=None:
            q=self.parseatom(line1)
            if len(q)==4:
                self.name1,self.a,self.b,self.c=q
        else:
            self.name1 = 'NONAME'
            self.a=0.0
            self.b=0.0
            self.c=0.0
        if line2!=None:
            q=self.parseatom(line2)
            if len(q)==4:
                self.name2, self.alpha, self.beta, self.gamma=q
            else:
                print 'в строке не хватает данных: {0:s}'.format()
                raise SystemExit(1)
            """        else:
            self.name2 = 'NONAME'
            self.alpha=0.0
            self.beta=0.0
            self.gamma=0.0"""
    def __str__(self):
        return '&CELL\n    ABC   {0:6.4f} {1:6.4f} {2:6.4f}\n'.format(self.a, self.b, self.c)+\
            '    ALPHA_BETA_GAMMA {0:6.4f} {1:6.4f} {2:6.4f}\n'.format(self.alpha, self.beta,self.gamma)+\
            '    PERIODIC XYZ\n&END CELL\n'


class CP2K():
    def __init__(self,filename=None):
        if filename!=None:
            filecp2k=open(filename, 'r')
            file_cp2k=filecp2k.readlines()
            filecp2k.close()
            self.atoms=[]
            self.cell=None
            i=0
            while i< len(file_cp2k):
                line=file_cp2k[i]
                tmp=line.split()
                if tmp[0]=='&CELL':
                    i=i+1
                    #while file_cp2k[i].strip()[:9]!='#SYMMETRY':
                    self.cell=CP2K_Cell(line1=file_cp2k[i], line2=file_cp2k[i+1]) #получение информации о ячейке из строки 1 ABC и строки 2 ALPHA_BETA_GAMMA
                elif tmp[0]=="SCALED":
                    i=i+1
                    while file_cp2k[i].strip()[:4]!='&END': #получение информации о атомах
                        handler=CP2K_Atom(line=file_cp2k[i])
                        self.atoms.append(handler)
                        i+=1
                i+=1

    def __add__(self, otherobject):
        resobject=CP2K()
        resobject.atoms = self.atoms
        resobject.atoms+=otherobject.atoms   #добавление к данным одного файла -f *.cp2k об атомах данных со второго файла -n *.cp2k
        resobject.cell = self.cell     # только cp2k_object1
        return resobject               #если есть жедание взять информацию о ячейке из второго файла, а не из первого, поменять при вводе файлы при -f и -n местами

    def droptofile(self,filename):
        cp2k_out=str(self.cell)
        cp2k_out+='&COORD\n    SCALED .TRUE.\n'
        for atom in self.atoms:
            cp2k_out+='    '+str(atom)
        cp2k_out+='&END COORD'
        f=open(filename,'w')
        f.write(cp2k_out)
        f.close()

    def xyzcoord(self):      #перевод координат из дробных координат в абсолютные
        cell=self.cell
        G11=cell.a
        G12=cell.b*cos_deg(cell.gamma)
        G13=cell.c*cos_deg(cell.beta)
        G22=cell.b*sin_deg(cell.gamma)
        G23=cell.c*((cos_deg(cell.alpha)-cos_deg(cell.beta)*cos_deg(cell.gamma))/sin_deg(cell.gamma))
        omg=cell.a*cell.b*cell.c*sqrt(1-cos_deg(cell.alpha)*cos_deg(cell.alpha)-cos_deg(cell.beta)*cos_deg(cell.beta)-cos_deg(cell.gamma)*cos_deg(cell.gamma)+2*cos_deg(cell.alpha)*cos_deg(cell.beta)*cos_deg(cell.gamma))
        G33=omg/(cell.a*cell.b*sin_deg(cell.gamma))
        for atom in self.atoms:
            atom.coord_x=G11*atom.x+G12*atom.y+G13*atom.z
            atom.coord_y=G22*atom.y+G23*atom.z
            atom.coord_z=G33*atom.z
        result=str(len(self.atoms))+'\n\n'
        for atom in self.atoms:
            result+='{0:2s} {1:16.9f} {2:16.9f} {3:16.9f}\n'.format(atom.name, atom.coord_x, atom.coord_y, atom.coord_z)
        return result

"""
filecp2k=open(inputfile, 'r')
file_cp2k=filecp2k.readlines()
filecp2k.close()
z=file_cp2k.index('SCALED')
testatom1=CP2K_Atom(line=file_cp2k[z+1])
testatom2=CP2K_Atom(line=file_cp2k[z+2])
testatom3=CP2K_Atom(line=file_cp2k[z+3])
...
list_atoms = [testatom1, testatom2, testatom3, ...]
"""
infile=argv.index('-f')
inputfile = argv[infile + 1]
otherf= argv.index('-n')      #если один входной файл
otherfile = argv[otherf + 1]  # то эти строки лучше закомментировать
ofile=argv.index('-o')
outputfile = argv[ofile + 1]

cp2k_object1 = CP2K(filename=inputfile)
cp2k_object2=CP2K(filename=otherfile)   #закоментировать если выходной файл один
cp2k_object_sum=cp2k_object1+cp2k_object2  #закомменитровать строку если выходной файл один
print cp2k_object_sum.droptofile(outputfile)      # файл в cp2k. при вводе в консоли файл вывода *.cp2k
#print cp2k_object_sum.xyzcoord()             #вывод конкатенир. файла в  xyz

#print cp2k_object1.xyzcoord() # оставить если нужен перевод в xyz только для одного файла



#testatom=CP2K_Atom(line='C    9.3640770374506327E-02    1.6669280258609478E-01    5.0594117047451115E-01 1.000')
#print testatom.etwas




