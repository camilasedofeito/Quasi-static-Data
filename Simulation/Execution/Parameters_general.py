# Simulation's parameters

# Original script developed by Gonzalo Tancredi
# Modified and adapted by Noelia Olivera Rodríguez
# Year: 2026

# Particle size distribution
minParticleRadius =  0.003   
maxParticleRadius =  0.007

# Padding distance added around the block to avoid boundary artifacts
pad = 5.0 * maxParticleRadius

# Dimensions of the granular specimen (m)
blockSizeX = 0.4
blockSizeY = 0.3970
blockSizeZ = 0.4

density = 3.e+03    # Particle material density (kg/m³)

# Simulation domain limits
# The computational domain extends beyond the block dimensions
# by a padding distance to ensure proper neighbor search and
# boundary handling.
minPointX = -pad
minPointY = -pad
minPointZ = -pad

maxPointX = blockSizeX + pad
maxPointY = blockSizeY + pad
maxPointZ = blockSizeZ + pad

# Select contact interaction law:
#   'hertzian'        → Hertz–Mindlin viscoelastic contact
#   'dashpot'         → Linear elastic + viscous dashpot
#   'springdashpot'   → Linear spring–dashpot with restitution control
interactionParticles='hertzian' 

# Define type of interaction between particles
if (interactionParticles=='hertzian'):

  #HertzianViscoelasticFriction interaction:
  Acoef = 3.e-2 * minParticleRadius
  Ecoef = 4.6e+10
  nucoef = 0.245
  frictioncoef =0.16

elif (interactionParticles=='dashpot'):

  #Elastic+LinearDashpot interaction:
  modulusYoung = 1.e10
  frictionCoefficient = 0.60
  damping = 2.e5 * minParticleRadius
  cutoff = 1.1

elif (interactionParticles=='springdashpot'):
  modulusYoung = 1.e10
  frictionCoefficient = 0.60
  nuCoefficient = 0.30
  restitutionCoefficient = 0.6

#Ellastic interaction with walls:
normalK  = 13e+5 # Normal stiffness for elastic wall interaction

#Gravity:
ygravity = 9.81 # Gravitational acceleration in the Y-direction (m/s²)

