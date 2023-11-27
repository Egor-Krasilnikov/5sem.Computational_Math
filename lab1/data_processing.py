import matplotlib.pyplot as plt
import numpy as np

x3 = np.array([])
y3 = np.array([])
x4 = np.array([])
y4 = np.array([])
x5 = np.array([])
y5 = np.array([])

data = np.loadtxt('data3points.txt')
x3 = data[:, 0]
y3 = data[:, 1]   
x3 = np.log(x3)
data = np.loadtxt('data4points.txt')
x4 = data[:, 0]
y4 = data[:, 1]   
x4 = np.log(x4)
data = np.loadtxt('data5points.txt')
x5 = data[:, 0]
y5 = data[:, 1]   
x5 = np.log(x5)




plt.subplot()

plt.scatter(x3, y3, label='3 узла')
plt.scatter(x4, y4, label='4 узла')
plt.scatter(x5, y5, label='5 узлов')

plt.xlabel('ln(h)')
plt.ylabel('ln ошибки')
plt.title('ln(ошибки) от ln(h) с аппроксимирующей прямой для N = 3')
x3_app = x3[4:16]
y3_app = y3[4:16]

x3_app2 = x3[0:5]
y3_app2 = y3[0:5]
coefficients = np.polyfit(x3_app, y3_app, 1)
p = np.poly1d(coefficients)
plt.plot(x3_app, p(x3_app),)

coefficients2 = np.polyfit(x3_app2, y3_app2, 1)
p = np.poly1d(coefficients2)
plt.plot(x3_app2, p(x3_app2),)

equation = f'k = {coefficients[0]:.2f}'
equation2 = f'k2  = {coefficients2[0]:.2f}'
plt.text(-15, -3, equation, fontsize=6)
plt.text(-15, -5, equation2, fontsize=6)
plt.savefig('Picture.jpg')
plt.legend()
plt.show()