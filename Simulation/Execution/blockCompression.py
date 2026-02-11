# Bonded Granular Block Compression Simulation#

# Original implementation: F. Lopez & G. Tancredi
# Modified, extended and adapted by Noelia Olivera Rodríguez (2026)

# This script performs a cyclic compression–decompression simulation of a
# bonded granular block using ESyS-Particle (MPI version).
#
# The loading protocol is imposed by prescribing the motion of the lower
# wall ("floor") through a time-dependent displacement schedule. The upper
# wall and lateral walls remain fixed and interact elastically with the
# particles.
#
# The simulation includes:
#   - Parallel domain decomposition (3D MPI grid)
#   - Configurable particle–particle interaction laws
#   - Elastic particle–wall interactions
#   - Gravity loading
#   - Cyclic kinematic boundary condition
#   - Snapshot, force and checkpoint saving
#
# This script constitutes the core numerical implementation described in
# the Methods section of the manuscript.

"""
Defines L{runSimulation} function which executes compression simulation.
"""

#import relevant modules:
from esys.lsm import *
from esys.lsm.util import Vec3, BoundingBox, InstallInfo
from math import *
from Parameters_blockCompression import *
from Parameters_statusWall import *
import sys

infilename=sys.argv[1]
folder="3CiclosBajo"

class Loading(Runnable):
    """
    Objects of this class provide the loading mechanism for a compression
    simulation. A "floor" is moved with a constant velocity. 
    """
    def __init__(self,lsm,t0,velCompWall,floordisp,startinc,stopinc,tup,dt_snap_wp):
        Runnable.__init__(self)
        self.theLSM=lsm
        self.t0=t0
        self.velocity = velCompWall
        self.floordisp = floordisp
        self.acc = 2. * self.velocity * self.velocity / self.floordisp
        self.startinc = startinc
        self.stopinc = stopinc
        self.tup = tup
        self.dt_snap_wp = dt_snap_wp
        self.count=0
        self.t=0.
        self.dt = self.theLSM.getTimeStepSize()
        self.disp = 0.
        self.tinc = 0.

    def run(self):

        """
        Moves floor by an increment.
        """

        self.count=self.count+1
        self.t=self.theLSM.getTimeStep() * self.dt
        self.dx=0.
        if abs(self.t - self.startinc) < self.dt / 2:
            self.tinc = 0.0


        if(self.t<=self.startinc or ((self.t - self.startinc - self.dt ) % (self.stopinc + self.tup) > self.tup) or abs(self.disp) >= abs(self.floordisp- 1e-5)):
            self.dx = 0.
        else:
            # self.dx = self.velocity * self.dt 
            self.tinc = self.tinc + self.dt
            self.dx = self.velocity/abs(self.velocity)*self.acc * self.tinc * self.dt 
            
        self.theLSM.moveWallBy("floor",Vec3(0.0,self.dx,0.0))
        self.disp = self.disp + self.dx


    def __str__(self):
	return "Loading(t0={}, velocity={}, floordisp={}, startinc={}, stopinc={})".format(self.t0, self.velocity, self.floordisp, self.startinc, self.stopinc
    )



#define a subroutine implementing the ESyS-Particle simulation:
def runSimulation():
    """
    Initialises and runs the simulation.
    """
    #setVerbosity(True)

    #define some variables for the simulation:
    t0=0.0                # initial timestep

    nt = int(round(te / dt))  # total time expresed in timesteps
    dt_snap = int(round(t_snap / dt))       # timesteps between snapshots
    dt_snap_wp = dt_snap
    gravityVec = Vec3(0.0, -ygravity, 0.0)


    print nt, dt_snap
    
    numWorkers = numSubdomainsX * numSubdomainsY * numSubdomainsZ

    # initialise the ESyS-Particle simulation object:
    mySim=LsmMpi(
       numWorkerProcesses = numWorkers,
       mpiDimList = [numSubdomainsX , numSubdomainsY , numSubdomainsZ ]
    )

    # initialise the neighbour search and specify the type of particles:
    mySim.initNeighbourSearch (
       particleType = "NRotSphere", 
       gridSpacing = 2.2 * maxParticleRadius, 
       verletDist = 0.5 * minParticleRadius
    )

    # set the timestep increment for the simulation:
    mySim.setTimeStepSize(dt)

    # read the particles and bonds from a geometry file:
    mySim.readGeometry(infilename)
    
    #specify the spatial domain for the simulation:
    domain = BoundingBox(Vec3(minPointX,minPointY,minPointZ), Vec3(maxPointX,maxPointY,maxPointZ))
    mySim.setSpatialDomain(domain)

    # set the particles density
    mySim.setParticleDensity (
      tag = 10,
      mask = -1,
      Density = density
    )
    
    # create floor and walls:
    mySim.createWall(
       name = "floor",
       posn = Vec3(0.0,floorPosition,0.0),
       normal = Vec3(0.0,1.0,0.0)
    )    
    mySim.createWall(
       name = "roof",
       posn = Vec3(0.0,heightRoof,0.0),
       normal = Vec3(0.0,-1.0,0.0)
    )
    mySim.createWall(
       name = "x-Wall",
       posn = Vec3(0.0,0.0,0.0),
       normal = Vec3(1.0,0.0,0.0)
    )
    mySim.createWall(
       name = "x+Wall",
       posn = Vec3(blockSizeX,0.0,0.0),
       normal = Vec3(-1.0,0.0,0.0)
    )
    mySim.createWall(
       name = "z-Wall",
       posn = Vec3(0.0,0.0,0.0),
       normal = Vec3(0.0,0.0,1.0)
    )
    mySim.createWall(
       name = "z+Wall",
       posn = Vec3(0.0,0.0,blockSizeZ),
       normal = Vec3(0.0,0.0,-1.0)
    )
    
    if (interactionParticles=='hertzian'):
        # define hertzian interaction between particles:
        interactionName = "HertzianViscoElasticFriction"
        mySim.createInteractionGroup (
          HertzianViscoElasticFrictionPrms(
            name = "HertzianViscoElasticFriction",
            A = Acoef,
            E = Ecoef,
            nu = nucoef,
            dynamicMu = frictioncoef,
            shearK = 0.4 * Ecoef
          )
        )
    elif (interactionParticles=='dashpot'):
        #define elastic+dashpot interaction between particles:
        interactionName="elastic_friction"
        mySim.createInteractionGroup(
          NRotFrictionPrms(
            name="elastic_friction",
            normalK=modulusYoung,
            dynamicMu=frictionCoefficient,
            shearK=0.4 * modulusYoung,
            scaling = True
          )    
        )
        mySim.createInteractionGroup(
          LinearDashpotPrms(
            name = "LinearDashpot",
            damp = damping,
            cutoff = cutoff
          )    
        )
    elif (interactionParticles=='springdashpot'):
        #define spring+dashpot+friction interaction between particles:
        interactionName="spring_dashpot_friction"
        mySim.createInteractionGroup(
            SpringDashpotFrictionPrms(
                name="spring_dashpot_friction",
                youngsModulus=modulusYoung,
                poissonsRatio=nuCoefficient,
                dynamicMu=frictionCoefficient,
                restitution=restitutionCoefficient
            )
        )

    mySim.createInteractionGroup (
      GravityPrms (
        name="gravity",
	acceleration=gravityVec
      )
    )

    # initialise unbonded elastic repulsion of particles from the floor and walls:
    wp1=NRotElasticWallPrms(
       name = "y-WallInteraction",
       wallName = "floor",
       normalK = normalK
    )    
    wp2=NRotElasticWallPrms(
       name = "y+WallInteraction",
       wallName = "roof",
       normalK = normalK
    )    
    wp3=NRotElasticWallPrms(
       name = "x-WallInteraction",
       wallName = "x-Wall",
       normalK = normalK
    )
    wp4=NRotElasticWallPrms(
       name = "x+WallInteraction",
       wallName = "x+Wall",
       normalK = normalK
    )
    wp5=NRotElasticWallPrms(
       name = "z-WallInteraction",
       wallName = "z-Wall",
       normalK = normalK
    )
    wp6=NRotElasticWallPrms(
       name = "z+WallInteraction",
       wallName = "z+Wall",
       normalK = normalK
    )
    
    # create the interaction groups for the particle-wall interactions:
    mySim.createInteractionGroup(wp1)
    mySim.createInteractionGroup(wp2)
    mySim.createInteractionGroup(wp3)
    mySim.createInteractionGroup(wp4)
    mySim.createInteractionGroup(wp5)
    mySim.createInteractionGroup(wp6)
    
    #create a FieldSaver to store position of moving flor:
    wallPosition = WallVectorFieldSaverPrms(
	wallName=["floor"],
	fieldName="Position",
	fileName=folder+"/floorPosition.dat",
	fileFormat="RAW_SERIES",
	beginTimeStep=0,
	endTimeStep=nt,
	timeStepIncr=dt_snap
      )
    mySim.createFieldSaver(wallPosition)
    
    #create a FieldSaver to wall forces:
    wallForce = WallVectorFieldSaverPrms(
	wallName=["roof"],
	fieldName="Force",
	fileName=folder+"/roofForce.dat",
	fileFormat="RAW_SERIES",
	beginTimeStep=0,
	endTimeStep=nt,
	timeStepIncr=dt_snap
	)
    mySim.createFieldSaver(wallForce)

    wallForce = WallVectorFieldSaverPrms(
	wallName=["x+Wall"],
	fieldName="Force",
	fileName=folder+"/x+WallForce.dat",
	fileFormat="RAW_SERIES",
	beginTimeStep=0,
	endTimeStep=nt,
	timeStepIncr=dt_snap
	)
    mySim.createFieldSaver(wallForce)

    wallForce = WallVectorFieldSaverPrms(
	wallName=["x-Wall"],
	fieldName="Force",
	fileName=folder+"/x-WallForce.dat",
	fileFormat="RAW_SERIES",
	beginTimeStep=0,
	endTimeStep=nt,
	timeStepIncr=dt_snap
	)
    mySim.createFieldSaver(wallForce)

    wallForce = WallVectorFieldSaverPrms(
	wallName=["z+Wall"],
	fieldName="Force",
	fileName=folder+"/z+WallForce.dat",
	fileFormat="RAW_SERIES",
	beginTimeStep=0,
	endTimeStep=nt,
	timeStepIncr=dt_snap
	)
    mySim.createFieldSaver(wallForce)
    
    wallForce = WallVectorFieldSaverPrms(
	wallName=["z-Wall"],
	fieldName="Force",
	fileName=folder+"/z-WallForce.dat",
	fileFormat="RAW_SERIES",
	beginTimeStep=0,
	endTimeStep=nt,
	timeStepIncr=dt_snap
	)
    mySim.createFieldSaver(wallForce)

    wallForce = WallVectorFieldSaverPrms(
	wallName=["floor"],
	fieldName="Force",
	fileName=folder+"/floorForce.dat",
	fileFormat="RAW_SERIES",
	beginTimeStep=0,
	endTimeStep=nt,
	timeStepIncr=dt_snap
	)
    mySim.createFieldSaver(wallForce)

    #create a FieldSaver to forces among the particles:
    partForce = InteractionVectorFieldSaverPrms(
        interactionName = interactionName,
        fieldName = "force",
	fileName=folder+"/partForce",
        fileFormat = "RAW_WITH_POS_ID",
	beginTimeStep=0,
	endTimeStep=nt,
	timeStepIncr=dt_snap
        )
    mySim.createFieldSaver (partForce)

    # create a CheckPointer to store the simulation state (snapshots):
    cp=CheckPointPrms(
       fileNamePrefix = folder+"/cpt",
       beginTimeStep = 0,
       endTimeStep = nt,
       timeStepIncr = dt_snap
    )
    mySim.createCheckPointer(cp)

    mySim.createFieldSaver (
       ParticleScalarFieldSaverPrms(
          fieldName="e_kin",
          fileName="ekin.dat",
          fileFormat="SUM",
          beginTimeStep=0,
          endTimeStep=nt,
          timeStepIncr=dt_snap
      )
   )

    # add the Loading Runnable to move the walls:

    # Adjust velocity to move floordisp
    velCompWall_adj = floordisp / time_per_phase

    movement_duration = time_per_phase  
    startinc = 0.0
    stopinc = startinc + movement_duration
    stopinc= movement_duration
    tup = time_per_phase - movement_duration
    tup= movement_duration

    for i in range(num_phases):
        t_start = i * time_per_phase
        startinc_adj = t_start + startinc
        stopinc_adj = t_start + stopinc
        vel = velCompWall_adj if i % 2 == 0 else -velCompWall_adj
        varname = "lf{}".format(i)
        globals()[varname] = Loading(mySim, t_start, vel, floordisp, startinc_adj, stopinc_adj,tup, dt_snap_wp)
        mySim.addPreTimeStepRunnable(globals()[varname])
        print(globals()[varname])

    mySim.setNumTimeSteps(nt)

    # execute the simulation:
    mySim.run()

if (__name__ == "__main__"):
    runSimulation()
