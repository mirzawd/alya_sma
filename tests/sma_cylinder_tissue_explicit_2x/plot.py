import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df  = pd.read_csv("data1.dat")

df.plot(figsize=(15,5), kind='line',x=0, y=1)
plt.show()
















