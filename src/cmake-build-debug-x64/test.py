import numpy as np
import treecode as tree

if __name__ == "__main__":
    bod = tree.body()
    bod.setPos([1, 2, 3])
    print(bod.getPos())