import imageio
import os

with imageio.get_writer("aop.gif", mode="I") as writer:
    for file in os.listdir("out\\aop\\"):
        image = imageio.imread("out\\aop\\" + file)
        writer.append_data(image)

with imageio.get_writer("dop.gif", mode="I") as writer:
    for file in os.listdir("out\\dop\\"):
        image = imageio.imread("out\\dop\\" + file)
        writer.append_data(image)