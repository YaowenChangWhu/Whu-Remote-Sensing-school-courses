import numpy as np

def solve_fit(A, X):
   AN = np.array(A)
   TAN = np.linalg.inv(AN)
   XN = np.array(X)
   fit = np.dot(TAN, XN)
   return fit

def cacluate(fit):
    B = np.array([1, 5, 25, 125])
    prediction = np.dot(B, fit)
    return prediction

L1 = []
L2 = []
file = open("PRN5_HaveCycleSlips.txt", "r")
lines = file.readlines()
for line in lines:
    temp = line.split()
    temp1 = temp[1]
    temp2 = temp[2]
    L1.append(temp1)
    L2.append(temp2)
file.close()

flt_L1 = [float(x) for x in L1]
flt_L2 = [float(x) for x in L2]
float_L1 = np.array(flt_L1)
float_L2 = np.array(flt_L2)
length = len(float_L1)   # 这里因为L1与L2的长度是相同的所以长度就用L1
i = 0
flag1 = []
flag2 = []
deltaL1=[]
deltaL2=[]
L1cycles_slips_value = []
L2cycles_slips_value = []
A = [[1,1,1,1],
    [1,2,4,8],
    [1,3,9,27],
    [1,4,16,64]]

while i < (length-4):    #多项式拟合求解的过程
    tempL1 = [[float_L1[i]],
              [float_L1[i+1]],
              [float_L1[i+2]],
              [float_L1[i+3]]]
    fit1 = solve_fit(A, tempL1)
    prediction1 = cacluate(fit1)
    deltaL1.append(float_L1[i + 4] - prediction1)
    if abs(float_L1[i+4]-prediction1) < 1:
       i += 1
    else:
        print("L1检测到周跳")
        flag1.append(i+5)  #存储周跳的位置
        L1cycles_slips_value.append(-(float(float_L1[i + 4] - prediction1)))
        float_L1[i+4:length] = (float_L1[i+4:length]-(float_L1[i+4]-prediction1)).tolist()
        i += 1

i = 0
while i < (length - 4):  # 多项式拟合求解的过程
    tempL2 = [[float_L2[i]],
              [float_L2[i + 1]],
              [float_L2[i + 2]],
              [float_L2[i + 3]]]
    fit2 = solve_fit(A, tempL2)
    prediction2 = cacluate(fit2)
    deltaL2.append(float_L2[i + 4] - prediction2)
    if abs(float_L2[i+4] - prediction2) < 1:
        i += 1
    else:
        print("L2检测到周跳")
        flag2.append(i+5)  # 存储周跳的位置
        L2cycles_slips_value.append(-(float(float_L2[i + 4] - prediction2)))
        float_L2[i+4:length] = (float_L2[i + 4:length] - (float_L2[i+4] - prediction2)).tolist()
        i += 1

print("L1发生周跳的次数是： ", len(flag1))
print("L1发生周跳的历元数是：", flag1)
print("周跳大小分别为：", L1cycles_slips_value)
print("L2发生周跳的次数是： ", len(flag2))
print("L2发生周跳的历元数是：", flag2)
print("周跳大小分别为：", L2cycles_slips_value)