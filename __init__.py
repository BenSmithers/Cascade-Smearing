from cascade.utils import config
import os

def pathmaker(path):
    if not os.path.exists(path):
        broken = path.split("/")
        working = ""
        for part in broken:
            working = "/".join([working] + [part])
            if not os.path.exists(working):
                os.mkdir(working)
                print("Making {} directory".format(working))

pathmaker(config["datapath"])
pathmaker(os.path.join(config["datapath"], config["img_subdir"]))
