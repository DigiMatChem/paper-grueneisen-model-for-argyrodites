from velocity import get_velocity_ph_bs_symm_line, VelocityPhononBSPlotter


# you have to add all labels in the labels_dict
bs=get_velocity_ph_bs_symm_line(velocity_path="band.yaml", labels_dict={"\Gamma":[0.0,0.0,0.0],
                                                                            "X":[0.5,0,0],
                                                                            "S": [0.5,0.5,0],
                                "Y":[0.,0.5,0], "Z":[0,0,0.5],"R":[0.5,0.5,0.5], "T":[0.0,0.5,0.5], "U":[0.5,0,0.5]})

VelocityPhononBSPlotter(bs=bs).save_plot_velocity_bs("velocities.eps")


# you can also only plot a few bands to see everything better (e.g., the first 3):

VelocityPhononBSPlotter(bs=bs).save_plot_velocity_bs("velocities_1_2_3.eps", only_bands=[0,1,2])

