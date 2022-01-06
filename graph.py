import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from scipy import optimize

#Собираем датасеты:

data = []

f = open("./data.csv")

for line in f:
    point = line.split(',')
    data.append(np.float64(point[1:]))

data = np.array(data)

number_of_frames =  data.shape[0]
print(number_of_frames)
def update_hist(num, data):
    plt.cla()
    plt.hist(data[num], bins = 20)

fig = plt.figure()
hist = plt.hist(data[0])

animation = ani.FuncAnimation(fig, update_hist, number_of_frames-1, fargs=(data, ), interval=40 )
#plt.show()

writergif = ani.PillowWriter(fps=30) 
animation.save("./result.gif", writer=writergif)