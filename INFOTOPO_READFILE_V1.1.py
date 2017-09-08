"""
Created on Thu Nov  5 16:06:46 2015
Copyright (2007-2017) Pierre Baudot 
The terms of the licence GNU GPL is reproduced in the file "COPYING" of INFOTOPO distribution.

This file is part of INFOTOPO.

  INFOTOPO is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  INFOTOPO is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with INFOTOPO.  If not, see <http://www.gnu.org/licenses/>. 
  
@author: baudot & al. 
GENERAL PHILOSOPHY AND AUTORS REQUEST: 
Giving the authorization of free use of this program also gives you the responsability 
of its use: we trust you. Every tool can be used to make things better for all, but also worse.
We ask any users of the present program to use and developp it such that it profits to the 
sustainable development of all, societies and environment and avoid development
in selfish interest, conflictual involvement, or individual freedom restrictions or survey. 
As a Human in the world, please be conscious of what you do, of its consequence on others, and try to make it better.
This tools comes from mathematic and pertain to it, mathematic do not pertain to anybody in particular 
and according to its original developments promotes harmony and sustainable developements on short 
and long time and space scales, please respect this ancestral tradition, old as humanity, and develop it in this respect. 
More information on the practical use and mathematical framework and bibliographic citations:
http://www.biorxiv.org/content/early/2017/07/26/168740 
and README file. 
For any requests, questions, improvements, developments (etc.) contact pierre.baudot [at] gmail.com
Friendly ergonomic interface will be developped afterward, sorry, please see the README file and/or contact me for use. 
INFOTOPO has been developped since 2007 thanks to grants-fundings of:
ISC-PIF (Complex system institute Paris Ile de France) (2007-2013) 
Max Planck Institute for Mathematic in the Sciences (MPI-MIS, Leipzig)   (2013-2015)
Inserm Unis 1072 - ERC channelomics (2015-2017)
We thank those public and non-profit scientific organization for their support.
"""
import pickle
import os

directory_path =  os.path.abspath("/home/baudot/LABOMARSEILLE/RESULTAT ANALYSE PCR/ACTUAL RESEARCH/ANALYSE 21RNA MARS/ANALYSE_nbtrial/m28/")
Nb_var_input=21
Nb_var_tot=21

name_object=  'ENTROPY_ORDERED'

def load_obj(name ):
    with open(directory_path + '/'   + name + '.pkl', 'rb') as f:       
        return pickle.load(f)

def factorial(x):
     if x < 2:
         return 1
     else:
         return x * factorial(x-1)

def binomial(n,k):
    return factorial(n)/(factorial(k)*factorial(n-k))        
#Ninfomut_per_order_ordered=[]
Ninfomut_per_order_ordered = load_obj(name_object)   

#for x in range(0,(2**Nb_var_tot)-1):  
#print(Ninfomut_per_order_ordered[2])

"""
for x in range(1,Nb_var_input+1):
#        print('les infos mutuelles à ', x ,'sont:')
#       print(Ninfomut_per_order_ordered[x])
       for i in range(0,11):
#           print(i)   
#           print(list(Ninfomut_per_order_ordered[x].keys())[0])
           print('le',i+1 ,'MIN infomut à ', x ,':',list(Ninfomut_per_order_ordered[x].keys())[i],list(Ninfomut_per_order_ordered[x].values())[i])
     #      print('le',list(Ninfomut_per_order_ordered[x].keys())[i])
#           print(list(Ninfomut_per_order_ordered[x].values())[i])
           print('le',i+1 ,'MAX infomut à ', x ,':',list(Ninfomut_per_order_ordered[x].keys())[int(binomial(Nb_var_input,x)-(i+1))],list(Ninfomut_per_order_ordered[x].values())[int(binomial(Nb_var_input,x)-(i+1))])
#           print(binomial(Nb_var_input,x)-(i+1))
#           print(list(Ninfomut_per_order_ordered[x].keys())[int(binomial(Nb_var_input,x)-(i+1))])
#           print(list(Ninfomut_per_order_ordered[x].values())[int(binomial(Nb_var_input,x)-(i+1))])
"""

print(Ninfomut_per_order_ordered[6])