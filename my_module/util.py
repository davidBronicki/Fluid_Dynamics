import numpy as np
from random import random
def difference(funct1, funct2):
    return lambda x: funct1(x) - funct2(x)

def linInterp(List, x0, dx, x):
    x -= x0
    x /= dx
    weight = x - int(x)
    x = int(x)
    if (x >= len(List) - 1):
        weight = x + weight - (len(List) - 2)
        x = len(List) - 2
    if (x < 0):
        weight = x + weight
        x = 0
    return weight * List[x + 1] + (1 - weight) * List[x]

def window(inputList):
    for i in range(len(inputList)):
        xVal = i / len(inputList)
        inputList[i] *= np.sin(np.pi * xVal)**2
    return np.array(inputList)

def freqWindow(inputList):
    for i in range(len(inputList)):
        xVal = i / len(inputList)
        if xVal < .5:
            inputList[i]

def interpFunct(List, x0, dx):
    return lambda x: linInterp(List, x0, dx, x)

def df(funct, x, dx, order = 1):
    if order == 1:
        temp = funct(x + dx / 2) - funct(x - dx / 2)
        return temp / dx
    else:
        return df(lambda y: df(funct, y, dx, order = order - 1), x, dx)

def dfFunct(funct, dx, order = 1):
    return lambda x: df(funct, x, dx, order = order)

def ddf(funct, x, dx):
    temp = funct(x + dx) - 2 * funct(x) + funct(x - dx)
    return temp / (dx * dx)

def ddfFunct(funct, dx):
    return lambda x: ddf(funct, x, dx)

def integ(funct, x0, x1, dx, iType = 'definite', argType = 'dx'):
    parts = dx
    if argType == 'partitions':
        dx = (x1 - x0) / part
    else:
        parts = int((x1 - x0) / dx)
        dx = (x1 - x0) / parts
    if iType == 'indefinite':
        List = integ(funct, x0, x1, dx, iType = 'list', argType = argType)
        return interpFunct(List, x0, dx)
    elif iType == 'list':
        output = [0]
        for i in range(1, parts):
            x = x0 + i * dx
            temp = funct(x - dx / 2) * 4 + funct(x - dx) + funct(x)
            output.append(output[-1] + temp * dx / 6)
        return output
    else:
        output = 0
        try:
            output += funct(x0)
        except:
            print('problem')
        try:
            output += funct(x1)
        except:
            print('problem')
        for i in range(1, parts, 2):
            try:
                output += funct(x0 + i * dx) * 4
            except:
                print('problem in integral')
                continue
        for i in range(2, parts, 2):
            try:
                output += funct(x0 + i * dx) * 2
            except:
                print('problem in integral')
                continue
        return output * dx / 3

def bisect(funct, lower, upper, tol, funct2 = lambda x: 0):
    funct = difference(funct, funct2)
    width = upper - lower
    midpoint = (upper + lower) / 2
    while width > tol:
        if (funct(lower) * funct(midpoint)) > 0:
            lower = midpoint
        else:
            upper = midpoint
        width = upper - lower
        midpoint = (upper + lower) / 2
    return midpoint

def newton(funct, start, tol, dif = None, funct2 = lambda x: 0):
    funct = difference(funct, funct2)
    if dif == None:
        dif = dfFunct(funct, tol)
    width = 10 * tol
    x = start
    while width > tol:
        lastX = x
        x = lastX - funct(lastX) / dif(lastX)
        width = abs(x - lastX)
    return x

def inverse(funct, tol, x0 = 0):
    return lambda x: newton(funct, x0, tol, funct2 = lambda y: x)

def fourier(List, dx, maxFrequency = 0):
    dF = 1 / (dx * len(List))
    newList = abs(np.fft.rfft(List))
    F = np.arange(0, len(newList) * dF - dF / 10, dF)
    if maxFrequency == 0:
        return F, newList
    else:
        n = maxFrequency / dF + 1
        return F[:n], newList[:n]

def splineInterpolate(interpolationList, inputList):
    def outputFunction(x):
        i = 0
        sign = inputList[-1] - inputList[0]
        sign /= abs(sign)
        for val in inputList:
            if x*sign < val*sign: break
            else: i+=1
        i-=1
        if i == 0 or i == -1:
            i = 1
        if i == len(inputList) - 2 or i == len(inputList) - 1:
            i -= 1
        if i == -1:
            dx = inputList[2]-inputList[0]
            x1 = inputList[0] - dx
            y2 = interpolationList[0]
            temp = interpolationList[2]
            dy = temp - y2
            k1 = dy/dx
            k2 = k1
            y1 = y2 - dy
        elif i == len(inputList) - 1:
            dx = inputList[-1]-inputList[-3]
            x1 = inputList[-1]
            y1 = interpolationList[-1]
            temp = interpolationList[-3]
            dy = y1 - temp
            k1 = dy/dx
            k2 = k1
            y2 = y1 + dy
        else:
            x1 = inputList[i]
            x2 = inputList[i+1]
            dx = x2 - x1
            y1 = interpolationList[i]
            y2 = interpolationList[i+1]
            dy = y2 - y1
            k1 = (y2 - interpolationList[i-1])/(x2-inputList[i-1])
            k2 = (interpolationList[i+2] - y1)/(inputList[i+2]-x1)

        t = (x-x1)/dx
        a = k1*dx - dy
        b = -k2*dx + dy
        return (1-t)*y1 + t*y2 + t*(1-t)*(a*(1-t)+b*t)
    return outputFunction

def randomPointPseudoNorm(dim):
    vals = []
    for i in range(dim):
        vals.append(2*random() - 1)
    return np.array(vals)

def metropolisAlgorithm(checkingFunction, pointCount, initialPoint, dryRunCount = 1000, dx = 1):
    if type(initialPoint) == float or type(initialPoint) == int:
        x = np.array([initialPoint])
    else:
        x = initialPoint.copy()
    y = checkingFunction(x)
    dim = len(x)
    xList = []
    accepted = 0
    for i in range(pointCount + dryRunCount):
        nextX = x + dx * randomPointPseudoNorm(dim)
        nextY = checkingFunction(nextX)
        weight = nextY / y
        if (weight >= 1):
            xList.append(nextX)
            x = nextX
            y = nextY
            accepted += 1
        else:
            r = random()
            if weight >= r:
                xList.append(nextX)
                x = nextX
                y = nextY
                accepted += 1
            else:
                xList.append(x)
    print(accepted / (pointCount + dryRunCount))
    return xList[dryRunCount:]

def transpose(twoDList):
    return np.array(list(zip(*twoDList)))

def parseCSVFile(fileString, separationCharacter = ','):
    file = open(fileString)
    output = []
    for line in file:
        output.append(line.split(separationCharacter))
    file.close()
    return output

def floatParseCSV(fileString, separationCharacter = ','):
    tempList = parseCSVFile(fileString, separationCharacter)
    newList = []
    for item in tempList:
        try:
            newTemp = []
            for thing in item:
                newTemp.append(float(thing))
            newList.append(newTemp)
        except:
            pass
    return newList
