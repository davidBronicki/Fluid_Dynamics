import numpy as np
import matplotlib.pyplot as plt
from random import random
import math
import imageio as img
import my_module.util as tools

class gridValue:
	def __init__(self, temperature, density, velocity,
			conductivity = 1.0, heatCap = 1.0, idealCoef = 1.0, diffusionCoef = 1.0):
		self.pressure = temperature * density * idealCoef
		self.temperature = temperature
		self.density = density
		self.velocity = velocity
		self.conductivity = conductivity
		self.heatCap = heatCap
		self.idealCoef = idealCoef
		self.diffuseCoef = diffusionCoef

	def couple(self, axisPairs):
		self.axisPairs = axisPairs

	def getDelDotVel(self, dx):
		output = 0.0
		for pair in self.axisPairs:
			output += (pair[1].velocity - pair[0].velocity) / (2 * dx)
		return sum(output)

	def preUpdate(self, dx, dt, externForce):
		gradPressure = []
		for pair in self.axisPairs:
			gradPressure.append(
				(pair[1].pressure - pair[0].pressure) / (dx * 2))
		gradPressure = np.array(gradPressure)

		laplaceVel = []
		for i in range(len(self.velocity)):
			temp = 0
			for j in range(len(self.axisPairs)):
				left = self.axisPairs[j][0].velocity[i]
				right = self.axisPairs[j][1].velocity[i]
				center = self.velocity[i]
				val = (left + right) - 2 * center
				val /= dx**2
				temp += val
			laplaceVel.append(temp)
		laplaceVel = np.array(laplaceVel)

		laplaceTemp = 0.0
		for pair in self.axisPairs:
			left = pair[0].temperature
			right = pair[1].temperature
			center = self.temperature
			val = (left + right) - 2 * center
			val /= dx**2
			laplaceTemp += val

		velDotDelOfVel = []
		for pair in self.axisPairs:
			velDotDelOfVel.append(
				(pair[1].velocity - pair[0].velocity) / (dx * 2))
		velDotDelOfVel = tools.transpose(velDotDelOfVel)
		velDotDelOfVel *= self.velocity
		velDotDelOfVel = tools.transpose(velDotDelOfVel)
		velDotDelOfVel = sum(velDotDelOfVel)

		velDotDelOfDensity = []
		for pair in self.axisPairs:
			velDotDelOfDensity.append(
				(pair[1].density - pair[0].density) / (dx * 2)
			)
		velDotDelOfDensity *= self.velocity
		velDotDelOfDensity = sum(velDotDelOfDensity)

		velDotDelOfTemp = []
		for pair in self.axisPairs:
			velDotDelOfTemp.append(
				(pair[1].temperature - pair[0].temperature) / (dx * 2)
			)
		velDotDelOfTemp *= self.velocity
		velDotDelOfTemp = sum(velDotDelOfTemp)

		delDotVelAxisPairs = []
		for pair in self.axisPairs:
			temp = []
			temp.append(pair[0].getDelDotVel(dx))
			temp.append(pair[1].getDelDotVel(dx))
			delDotVelAxisPairs.append(temp)
		delOfDelDotVel = []
		for pair in delDotVelAxisPairs:
			delOfDelDotVel.append((pair[1] - pair[0]) / (2 * dx))
		delOfDelDotVel = np.array(delOfDelDotVel)

		delDotMomentum = 0.0
		for pair in self.axisPairs:
			dif = pair[1].density * pair[1].velocity[i]
			dif -= pair[0].density * pair[0].velocity[i]
			dif /= dx * 2
			delDotMomentum += dif

		self.dDensity = -dt * delDotMomentum
		tempVal = self.diffuseCoef * laplaceVel\
			+ self.diffuseCoef / 3 * delOfDelDotVel
		self.dVelocity = (-velDotDelOfVel
			- gradPressure / self.density + tempVal
			+ externForce) * dt
		self.dTemp = -velDotDelOfTemp\
			+ 1/(self.density * self.heatCap)\
			* ((self.dDensity + velDotDelOfDensity * dt)
			+ self.conductivity * laplaceTemp * dt
			- sum(self.velocity * tempVal * self.density) * dt)

	def update(self):
		self.density += self.dDensity
		# self.density = max(self.density, 0)
		self.velocity += self.dVelocity
		self.temperature += self.dTemp
		self.pressure = self.density * self.temperature * self.idealCoef
		# velSqr = sum(self.velocity * self.velocity)
		# if velSqr > cSqr:
		# 	self.velocity *= np.sqrt(cSqr / velSqr)
		# self.potential = self.pressure / self.density
		# self.potential = self.pressure

fluidGrid = []
sideLength = 40
dt = 0.25
dx = 1.0
# c = 1.0
# cSqr = c**2
T = 150.0
times = int(T / dt)

def position(i, j):
	return np.array([i * dx, j * dx])

def initialVelocity(pos):
	# return np.array([0.0, 0.0])
	x = pos[0]
	y = pos[1]
	# randComp = (random() - 0.5) * 0.001
	length = dx * sideLength
	randComp = np.sin(2 * np.pi * x / length) * 0.01\
		* np.cos(2 * np.pi * y / length)**2
	# randComp = 0.0
	# return np.array([0.05 * np.sin(2 * np.pi * x / length), randComp])
	if abs(y / length - 0.5) < 0.255:
		return np.array([0.05, randComp])
	else:
		return np.array([0.0, 0])

def initialDensity(pos):
	return 1.0
	# x = pos[0]
	# y = pos[1]
	# length = dx * sideLength
	# return 0.5 + 0.2 * (
	# 	np.sin(np.pi * x * 2 / length)**2 * np.sin(np.pi * y / length)**2 - 0.025)

def initialTemperature(pos):
	return 1.0

# def initialPressure(pos):
# 	return initialDensity(pos) + 0.12
	# x = pos[0]
	# y = pos[1]
	# length = dx * sideLength
	# return 0.6 + cSqr * 0.1 * (
	# 	np.sin(np.pi * x / length)**2 * np.sin(np.pi * y / length)**2 - 0.025)

def buildAxisPairs(i, j):
	lowerX = sideLength - 1 if i == 0 else i - 1
	upperX = 0 if i == sideLength - 1 else i + 1
	lowerY = sideLength - 1 if j == 0 else j - 1
	upperY = 0 if j == sideLength - 1 else j + 1
	# print(i,', ',j,'->',str([lowerX, upperX, lowerY, upperY]))

	return [[fluidGrid[lowerX][j], fluidGrid[upperX][j]],
			[fluidGrid[i][lowerY], fluidGrid[i][upperY]]]

for i in range(sideLength):
	temp = []
	for j in range(sideLength):
		pos = position(i, j)
		density = initialDensity(pos)
		# pressure = initialPressure(pos)
		temperature = initialTemperature(pos)
		velocity = initialVelocity(pos)
		temp.append(gridValue(temperature, density, velocity,
			diffusionCoef=0.00))
	fluidGrid.append(temp)

for i in range(sideLength):
	for j in range(sideLength):
		fluidGrid[i][j].couple(buildAxisPairs(i, j))

outGif = np.zeros((times, sideLength, sideLength, 3))

val = 0
cutoff = 0.00

def externForce(pos, time):
	return np.array([0.0, 0])

def colorFunct(pressure, velocity, temperature):
	pressure = min(max(pressure / 2, 0), 1)
	velocity = np.sqrt(sum(velocity * velocity))
	velocity *= 10
	velocity = min(max(velocity, 0), 1)
	temperature = min(max(temperature / 2, 0), 1)
	return np.array([temperature, pressure, velocity])
	# return np.array([density, 1 - density, 0])

def displayCondition(simPoint):
	return simPoint % 4 == 0
	# return True

for simPoint in range(times):
	# try:
	for i in range(sideLength):
		for j in range(sideLength):
			fluidGrid[i][j].preUpdate(dx, dt, externForce(pos, simPoint * dt))
	for i in range(sideLength):
		for j in range(sideLength):
			fluidGrid[i][j].update()
			if displayCondition(simPoint):
				outGif[val, j, i] = colorFunct(
					fluidGrid[i][j].pressure,
					fluidGrid[i][j].velocity,
					fluidGrid[i][j].temperature)
	if displayCondition(simPoint):
		val += 1
		print('working, time = ', simPoint * dt)
		print(simPoint * dt / T * 100, '% complete')
	# except:
	# 	print("failure on step ", simPoint)
	# 	break
	# if displayCondition(simPoint):
	# 	for i in range(sideLength):
	# 		for j in range(sideLength):
	# 			plt.plot([i, i + fluidGrid[i][j].velocity[0]],
	# 				[j, j + fluidGrid[i][j].velocity[1]],
	# 				color = 'b')
	# 			temp = min(max(
	# 				fluidGrid[i][j].density, 0), 1)
	# 			plt.scatter(i, j, color = (temp, 1 - temp, 0))
	# 	plt.show()

outGif = outGif[:val]
outGif = np.array(outGif)

img.mimsave('test2.gif', outGif)

# print("done")
