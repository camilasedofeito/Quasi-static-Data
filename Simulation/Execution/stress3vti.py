#stress3vti.py: produce a stress vti from interaction force RAW_WITH_ID data
# 
#   Copyright (C) 2020, Steffen Abe (1) and Dion Weatherley (2).
#       (1) Institute for Geothermal Energy, Mainz, Germany.
#       (2) The University of Queensland, Australia.
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# Usage
# python stress3vti.py force_datafile.0.dat stress_outputfile.0 Xmin Xmax Ymin Ymax Zmin Zmax Nx Ny Nz

# Choose Nx, Ny and Nz so that:
# (Xmax-Xmin)/Nx ~ (Ymax-Ymin)/Ny ~ (Zmax-Zmin)/Nz

from __future__ import print_function
import math
import vtk
import numpy
from sys import argv

class stress_array:
	
	# construct stress array
    def __init__ (self,x0,x1,y0,y1,z0,z1,nx,ny,nz):
        self.x0=x0
        self.y0=y0
        self.z0=z0
        self.x1=x1
        self.y1=y1
        self.z1=z1
        self.nx=nx
        self.ny=ny
        self.nz=nz
        # per-axis grid spacing
        self.dx=(self.x1-self.x0)/float(nx)
        self.dy=(self.y1-self.y0)/float(ny)
        self.dz=(self.z1-self.z0)/float(nz)
        # init empty cell list
        self.data=list()
		# initialize cells
        # (cx0, cy0, cz0) , (cx1, cy1, cz1): min/max corner of the cell  
        for i in range(self.nx):
            cx0=self.x0+float(i)*self.dx
            cx1=self.x0+float(i+1)*self.dx
            for j in range(self.ny):
                cy0=self.y0+float(j)*self.dy
                cy1=self.y0+float(j+1)*self.dy
                for k in range(self.nz) :
                    cz0=self.z0+float(k)*self.dz
                    cz1=self.z0+float(k+1)*self.dz
                    self.data.append(stress_cell(cx0,cx1,cy0,cy1,cz0,cz1))
		

    # get the cell index for a given x,y,z-position
    # returns -1 if position is outside grid
    def get_index(self,x,y,z):
        if (x>self.x1) or (x<=self.x0) or (y>self.y1) or (y<=self.y0) or (z>self.z1) or (z<=self.z0) :
            idx=-1
        else:
            xidx=int(math.floor((x-self.x0)/self.dx))
            yidx=int(math.floor((y-self.y0)/self.dy))
            zidx=int(math.floor((z-self.z0)/self.dz))
            idx=self.ny*self.nz*xidx+self.nz*yidx+zidx
		
        return idx
		
	def get_cell(self,idx):
		return self.data[idx]
		
    def add_force(self,f):
		# get cell index for the two points
        idx1=self.get_index(f.x0,f.y0,f.z0)
        idx2=self.get_index(f.x1,f.y1,f.z1)
        if idx1 != idx2 : # points in 2 different cells -> useable
            if(idx1 != -1):
                self.data[idx1].add_force(f)
            if(idx2 != -1):
                self.data[idx2].add_force(f.inv())
				
    # read interaction forces for a saver file in 
    # RAW_WITH_POS_ID format
    def read_data_raw_with_pos_id(self,filename):
        infile=open(filename)
        for line in infile:
            ls=line.split()
            # position of 1st particle
            x0=float(ls[2])
            y0=float(ls[3])
            z0=float(ls[4])
            # position of 2nd particle
            x1=float(ls[5])
            y1=float(ls[6])
            z1=float(ls[7])
            # force 
            fx=float(ls[11])
            fy=float(ls[12])
            fz=float(ls[13])
            
            f=force(x0,y0,z0,x1,y1,z1,fx,fy,fz)
            self.add_force(f)
		
        infile.close()
	
				
	def print_stress(self):
		for i in range(self.nx):
			for j in range(self.ny):	
				idx=self.ny*i+j
				#sxx,sxy,syx,syy=self.data[idx].calc_stress()
				print("cell - x: ",self.data[idx].x0,"-",self.data[idx].x1," y: ",self.data[idx].y0,"-",self.data[idx].y1, " stress",self.data[idx].sigma)
				
				
				
	# write the stress data as VTK image file using components (s_ij, s_1, s_2, s_3)			
    def write_as_vti(self,filename):
		# setup vtk image data
        out=vtk.vtkImageData()
        out.SetExtent(0,self.nx,0,self.ny,0,self.nz) # vtk image data for output
        out.SetSpacing(self.dx,self.dy,self.dz)
        out.SetOrigin(self.x0,self.y0,self.z0)
	
		# setup cell data arrays
		# scalar arrays for tensor components
        tdata_xx=vtk.vtkDoubleArray()
        tdata_xx.SetName("sigma_xx")
        tdata_xx.SetNumberOfTuples(out.GetNumberOfCells())
        tdata_xy=vtk.vtkDoubleArray()
        tdata_xy.SetName("sigma_xy")
        tdata_xy.SetNumberOfTuples(out.GetNumberOfCells())
        tdata_xz=vtk.vtkDoubleArray()
        tdata_xz.SetName("sigma_xz")
        tdata_xz.SetNumberOfTuples(out.GetNumberOfCells())

        tdata_yx=vtk.vtkDoubleArray()
        tdata_yx.SetName("sigma_yx")
        tdata_yx.SetNumberOfTuples(out.GetNumberOfCells())
        tdata_yy=vtk.vtkDoubleArray()
        tdata_yy.SetName("sigma_yy")
        tdata_yy.SetNumberOfTuples(out.GetNumberOfCells())
        tdata_yz=vtk.vtkDoubleArray()
        tdata_yz.SetName("sigma_yz")
        tdata_yz.SetNumberOfTuples(out.GetNumberOfCells())
		
        tdata_zx=vtk.vtkDoubleArray()
        tdata_zx.SetName("sigma_zx")
        tdata_zx.SetNumberOfTuples(out.GetNumberOfCells())
        tdata_zy=vtk.vtkDoubleArray()
        tdata_zy.SetName("sigma_zy")
        tdata_zy.SetNumberOfTuples(out.GetNumberOfCells())
        tdata_zz=vtk.vtkDoubleArray()
        tdata_zz.SetName("sigma_zz")
        tdata_zz.SetNumberOfTuples(out.GetNumberOfCells())
        
		# vector array for principal stresses
        s1data=vtk.vtkDoubleArray()
        s1data.SetName("sigma_1")
        s1data.SetNumberOfComponents(3)
        s1data.SetNumberOfTuples(out.GetNumberOfCells())
        s2data=vtk.vtkDoubleArray()
        s2data.SetName("sigma_2")
        s2data.SetNumberOfComponents(3)
        s2data.SetNumberOfTuples(out.GetNumberOfCells())
        s3data=vtk.vtkDoubleArray()
        s3data.SetName("sigma_3")
        s3data.SetNumberOfComponents(3)
        s3data.SetNumberOfTuples(out.GetNumberOfCells())
		
		# loop over grid
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    idx=self.ny*self.nz*i+self.nz*j+k
                    idv=out.ComputeCellId([i,j,k])
                    tdata_xx.SetValue(idv,self.data[idx].sigma[0,0])
                    tdata_xy.SetValue(idv,self.data[idx].sigma[0,1])
                    tdata_xz.SetValue(idv,self.data[idx].sigma[0,2])
                    
                    tdata_yx.SetValue(idv,self.data[idx].sigma[1,0])
                    tdata_yy.SetValue(idv,self.data[idx].sigma[1,1])
                    tdata_yz.SetValue(idv,self.data[idx].sigma[1,2])
                    
                    tdata_zx.SetValue(idv,self.data[idx].sigma[2,0])
                    tdata_zy.SetValue(idv,self.data[idx].sigma[2,1])
                    tdata_zz.SetValue(idv,self.data[idx].sigma[2,2])
                    
                    # get principal stresses
                    s1,s2,s3=self.data[idx].get_principal_stresses()
                    					
                    s1data.SetTuple3(idv,s1[0],s1[1],s1[2])
                    s2data.SetTuple3(idv,s2[0],s2[1],s2[2])
                    s3data.SetTuple3(idv,s3[0],s3[1],s3[2])
				
		# add point data arrays to image data 
        out.GetCellData().AddArray(tdata_xx)
        out.GetCellData().AddArray(tdata_xy)
        out.GetCellData().AddArray(tdata_xz)

        out.GetCellData().AddArray(tdata_yx)
        out.GetCellData().AddArray(tdata_yy)
        out.GetCellData().AddArray(tdata_yz)

        out.GetCellData().AddArray(tdata_zx)
        out.GetCellData().AddArray(tdata_zy)
        out.GetCellData().AddArray(tdata_zz)

        out.GetCellData().AddArray(s1data)
        out.GetCellData().AddArray(s2data)
        out.GetCellData().AddArray(s3data)

       
		# write data to file
        owriter=vtk.vtkXMLDataSetWriter()
        owriter.SetFileName(filename)
        owriter.SetDataModeToAscii()

        owriter.SetInputData(out)

        owriter.Write()
	
    # write the stress data as VTK image file - one tensor per point			
    def write_as_vti_tensor(self,filename):
		# setup vtk image data
        out=vtk.vtkImageData()
        out.SetExtent(0,self.nx-1,0,self.ny-1,0,self.nz-1) # vtk image data for output
        out.SetSpacing(self.dx,self.dy,self.dz)
        out.SetOrigin(self.x0,self.y0,self.z0)
        
        tdata=vtk.vtkDoubleArray()
        tdata.SetName("sigma")
        tdata.SetNumberOfComponents(9) # tensor data
        tdata.SetNumberOfTuples(out.GetNumberOfPoints())
        
        
        # loop over grid
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    idx=self.ny*self.nz*i+self.nz*j+k
                    idv=out.ComputePointId([i,j,k])
                    s1,s2,s3,s4,s5,s6,s7,s8,s9=self.data[idx].get_stress_components()
                    tdata.SetTuple9(idv,s1,s2,s3,s4,s5,s6,s7,s8,s9)
                    
        
        # add point data arrays to image data 
        out.GetPointData().AddArray(tdata)           
        
        # write data to file
        owriter=vtk.vtkXMLDataSetWriter()
        owriter.SetFileName(filename)
        owriter.SetDataModeToAscii()

        owriter.SetInputData(out)

        owriter.Write()
        
# class to accumulate stress in a box-shaped cell
class stress_cell:
	
	# initialize empty cell
    # x0,y0,z0: minimum corner
    # x1,y1,z1: maximum corner
    def __init__ (self,x0,x1,y0,y1,z0,z1):
        # set corner coordinates
        self.x0=x0
        self.y0=y0
        self.z0=z0
        self.x1=x1
        self.y1=y1
        self.z1=z1
		
        # calculate center coordinates
        self.cx=0.5*(x0+x1)
        self.cy=0.5*(y0+y1)
        self.cz=0.5*(z0+z1)
        
        # cell volume
        self.vol=(x1-x0)*(y1-y0)*(z1-z0)
		
        # initialize stress  to 0
        self.sigma=numpy.zeros((3,3))
		
		
    def add_force(self,f):
        #print ("cell - x: ",self.x0,"-",self.x1," y: ",self.y0,"-",self.y1," z: ",self.z0,"-",self.z1, " add force", f.x0,f.y0,f.z0,"|",f.x1,f.y1,f.z1,"|",f.fx,f.fy,f.fz)
        # -- determine border intersection from posn. of 2nd point (1st point is inside, 2nd point outside)
        # x
        if (f.x1<self.x0) : 
            d_isX,py_isX,pz_isX=f.intersectX(self.x0) # x1 below xmin -> x intersection
            px_isX=self.x0
        elif (f.x1>self.x1) :
            d_isX,py_isX,pz_isX=f.intersectX(self.x1) # x1 above xmax -> x intersection
            px_isX=self.x1
        else : # neither below or above -> no x intersection
            d_isX, px_isX,py_isX,pz_isX=numpy.inf, 0.0,0.0,0.0 # infinite distance -> can't be the closest intersection
        # y
        if (f.y1<self.y0) : 
            d_isY,px_isY,pz_isY=f.intersectY(self.y0) # y1 below ymin -> y intersection
            py_isY=self.y0
        elif (f.y1>self.y1) :
            d_isY,px_isY,pz_isY=f.intersectY(self.y1) # y1 above ymax -> y intersection
            py_isY=self.y1
        else : # neither below or above -> no y intersection
            d_isY, px_isY,py_isY,pz_isY=numpy.inf, 0.0,0.0,0.0 # infinite distance -> can't be the closest intersection
        
        # z
        if (f.z1<self.z0) : 
            d_isZ,px_isZ,py_isZ=f.intersectZ(self.z0) # z1 below zmin -> z intersection
            pz_isZ=self.z0
        elif (f.z1>self.z1) :
            d_isZ,px_isZ,py_isZ=f.intersectZ(self.z1) # z1 above zmax -> z intersection
            pz_isZ=self.z1
        else : # neither below or above -> no z intersection
            d_isZ, px_isZ,py_isZ,pz_isZ=numpy.inf, 0.0,0.0,0.0 # infinite distance -> can't be the closest intersection
        
            
        # check which intersection is closest & assign intersection point coordinates accordingly
        if (d_isX < d_isY) and (d_isX < d_isZ) : # min distance for x-intersection
            px,py,pz=px_isX,py_isX,pz_isX
        elif d_isY < d_isZ : # y
            px,py,pz=px_isY,py_isY,pz_isY
        else : # only z left
            px,py,pz=px_isZ,py_isZ,pz_isZ
            
        #print("intersection found:",px,py,pz)
		# found the point where the vector between the particle centers crosses the bounday of the cell 
		# dist. from center
        x=px-self.cx
        y=py-self.cy
        z=pz-self.cz
        #print("rel.pos:",x,y,z," force",f.fx,f.fy,f.fz)
        # calculate stress contribution 
        # essentially Eq. 52 of Potyony & Cundall 2004 IJRMMS 
        dsigma=numpy.zeros((3,3))
        dsigma[0,0]=x*f.fx/self.vol # sigma_xx
        dsigma[1,0]=y*f.fx/self.vol # sigma_yx
        dsigma[2,0]=z*f.fx/self.vol # sigma_zx
        dsigma[0,1]=x*f.fy/self.vol # sigma_xy
        dsigma[1,1]=y*f.fy/self.vol # sigma_yy
        dsigma[2,1]=z*f.fy/self.vol # sigma_zy
        dsigma[0,2]=x*f.fz/self.vol # sigma_xz
        dsigma[1,2]=y*f.fz/self.vol # sigma_yz
        dsigma[2,2]=z*f.fz/self.vol # sigma_zz
        #print(dsigma)
        # update stress
        self.sigma-=dsigma
        
    # calculate principal stresses from symmetric part of stress tensor
    def get_principal_stresses(self):
        # get eigenvalues & eigenvectors 
        w,v=numpy.linalg.eig(0.5*(self.sigma+self.sigma.transpose()))
        # get order of eigenvalues 
        order=numpy.argsort(w)
        s3=w[order[2]]*v[:,order[2]] # min ev -> sigma_3
        s2=w[order[1]]*v[:,order[1]] # mid ev -> sigma_2
        s1=w[order[0]]*v[:,order[0]] # max ev -> sigma_1
        
        return s1,s2,s3
    
    # return stress tensor components in column-major order 
    # convenience function of vti tensor field
    def get_stress_components(self) :
        return self.sigma[0,0], self.sigma[1,0], self.sigma[2,0], self.sigma[0,1], self.sigma[1,1], self.sigma[2,1], self.sigma[0,2], self.sigma[1,2], self.sigma[2,2]
		
# class for a force between two particles
class force:
	
    # initialize force
    # x0,y0,z0 : 1st particle position
    # x1,y1,z1 : 2nd particle position
    # fx,fz,fy : force vector
    def __init__(self,x0,y0,z0,x1,y1,z1,fx,fy,fz):
        self.x0=x0
        self.y0=y0
        self.z0=z0
        self.x1=x1
        self.y1=y1
        self.z1=z1
        self.fx=fx
        self.fy=fy
        self.fz=fz
        
    # return inverse force, i.e. swap 1st & 2nd particle position and multiply force vector by -1
    def inv(self):
        return force(self.x1,self.y1,self.z1,self.x0,self.y0,self.z0,-1.0*self.fx,-1.0*self.fy,-1.0*self.fz)
    
    # calculate intersection with x-normal plane
    # x : x-position of the plane
    # returns : dist, y, z if intersection
    # N.B. assumes there is an intersection, i.e. self.x0!=self.x1
    def intersectX(self,x):
        # calc y- and z- position of the intersection
        dxy=(self.y1-self.y0)/(self.x1-self.x0)
        py=self.y0+dxy*(x-self.x0)
        dxz=(self.z1-self.z0)/(self.x1-self.x0)
        pz=self.z0+dxz*(x-self.x0)
        # calc distance from point 0 to intersection
        dx=self.x0-x
        dy=self.y0-py
        dz=self.z0-pz
        dist=math.sqrt(dx*dx+dy*dy+dz*dz)
        
        return dist,py,pz
        
    # calculate intersection with y-normal plane
    # y : y-position of the plane
    # returns : dist, x, z if intersection
    def intersectY(self,y):
        # calc x- and z- position of the intersection
        dyx=(self.x1-self.x0)/(self.y1-self.y0)
        px=self.x0+dyx*(y-self.y0)
        dyz=(self.z1-self.z0)/(self.y1-self.y0)
        pz=self.z0+dyz*(y-self.y0)
        # calc distance from point 0 to intersection
        dx=self.x0-px
        dy=self.y0-y
        dz=self.z0-pz
        dist=math.sqrt(dx*dx+dy*dy+dz*dz)
        
        return dist,px,pz
        
    # calculate intersection with z-normal plane
    # z : Z-position of the plane
    # returns : dist, x, y if intersection, None if not
    def intersectZ(self,z):
        # calc y- and z- position of the intersection
        dzx=(self.x1-self.x0)/(self.z1-self.z0)
        px=self.x0+dzx*(z-self.z0)
        dzy=(self.y1-self.y0)/(self.z1-self.z0)
        py=self.y0+dzy*(z-self.z0)
        # calc distance from point 0 to intersection
        dx=self.x0-px
        dy=self.y0-py
        dz=self.z0-z
        dist=math.sqrt(dx*dx+dy*dy+dz*dz)
        
        return dist,px,py
        
    
    # convert to strig
    def __str__(self):
        return "force : ["+str(self.x0)+" , "+ str(self.y0)+" , "+str(self.z0)+"] , ["+str(self.x1)+" , "+ str(self.y1)+" , "+str(self.z1)+"] - [" +str(self.fx)+" , "+ str(self.fy)+" , "+str(self.fz)+"]"
        
		
"""
a=stress_array(1,19,1,39,1,19,8,16,8)
infilename_b="forces/bforce.9.dat"
#infilename_f="ffriction.6.dat"
outfilename="stress3_test1b.vti"
outfilename_t="stress3_test1b_t.vti"
"""

infilename_b=argv[1]
outfilename=argv[2]

"""
x0 = -1.0
x1 = 1.0
y0 = -1.0
y1 = 1.0
z0 = 0.0
z1 = 0.3
nx = 20
ny = 20
nz = 3
"""
x0 = float(argv[3])
x1 = float(argv[4])
y0 = float(argv[5])
y1 = float(argv[6])
z0 = float(argv[7])
z1 = float(argv[8])
nx = int(argv[9])
ny = int(argv[10])
nz = int(argv[11])

a=stress_array(x0,x1,y0,y1,z0,z1,nx,ny,nz)

a.read_data_raw_with_pos_id(infilename_b)
a.write_as_vti(outfilename+'.vti')
a.write_as_vti_tensor(outfilename+'.t.vti')

