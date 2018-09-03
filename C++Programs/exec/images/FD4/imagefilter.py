import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from scipy import ndimage

for i in range(20,440,20):
	print i
	# Load the data...
	im = np.genfromtxt("image"+str(i)+".txt")
	data = np.array(im, dtype=float)

	kernel = np.array([[-1, -1, -1, -1, -1],
					   [-1,  1,  2,  1, -1],
					   [-1,  2,  4,  2, -1],
					   [-1,  1,  2,  1, -1],
					   [-1, -1, -1, -1, -1]])

	highpass_5x5 = ndimage.convolve(data, kernel)
	
	np.savetxt("imagefiltered"+str(i)+".txt",highpass_5x5)

	plt.imshow(data,cmap="gray",extent=[0,9,3.5,0])
	plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%d km'))
	plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%d km'))
	plt.savefig("image"+str(i)+".eps",bbox_inches='tight')
	
	plt.imshow(highpass_5x5,cmap="gray",extent=[0,9,3.5,0])
	plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%d km'))
	plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%d km'))
	plt.savefig("imagefiltered"+str(i)+".eps",bbox_inches='tight')

im = np.genfromtxt("image.txt")
data = np.array(im, dtype=float)

kernel = np.array([[-1, -1, -1, -1, -1],
					   [-1,  1,  2,  1, -1],
					   [-1,  2,  4,  2, -1],
					   [-1,  1,  2,  1, -1],
					   [-1, -1, -1, -1, -1]])

highpass_5x5 = ndimage.convolve(data, kernel)

np.savetxt("imagefiltered.txt",highpass_5x5)

plt.imshow(data,cmap="gray",extent=[0,9,3.5,0])
plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%d km'))
plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%d km'))
plt.savefig("image.eps",bbox_inches='tight')

plt.imshow(highpass_5x5,cmap="gray",extent=[0,9,3.5,0])
plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%d km'))
plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%d km'))
plt.savefig("imagefiltered.eps",bbox_inches='tight')
