from velocity import VelocityPlotter
from velocity import get_velocity

vel=get_velocity("mesh.yaml")

VelocityPlotter(vel).save_plot("velocity_Kms-1.png")
