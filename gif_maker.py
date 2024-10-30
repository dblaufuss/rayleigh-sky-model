import imageio.v2 as imageio
import os

with imageio.get_writer("aop.gif", mode="I") as writer:
    directory = os.listdir("out\\aop\\")
    directory.sort(key=lambda x: os.path.getmtime("out\\aop\\"+x))

    for file in directory:
        image = imageio.imread("out\\aop\\" + file)
        writer.append_data(image)

with imageio.get_writer("dop.gif", mode="I") as writer:
    directory = os.listdir("out\\dop\\")
    directory.sort(key=lambda x: os.path.getmtime("out\\dop\\"+x))

    for file in directory:
        image = imageio.imread("out\\dop\\" + file)
        writer.append_data(image)