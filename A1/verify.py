import numpy as np
import argparse
import sys


def readfile(file, p=False):
	with open(file, 'r') as f:
		lines = f.readlines()
	n = len(lines)
	matrix = []
	for line in lines:
		nums = [float(x) for x in line.split()]
		matrix.append(nums)
	
	if p:
		n = len(matrix[0])
		pi = np.zeros((n, n))
		for i in range(n):
			pi[i][int(matrix[0][i])] = 1
		return pi
	else:
		return np.array(matrix)


def verify(infile):
	a = readfile(infile)
	l = readfile('LowerTri')
	u = readfile('UpperTri')
	p = readfile('Permutation', p=True)
	residual = (p @ a) - (l @ u)
	print(residual)
	l21norm = np.sum(np.linalg.norm(residual, axis=0))
	print(l21norm)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=str, help='File containing input matrix.')
    parser.set_defaults(func=verify)

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    args.func(args.infile)


if __name__ == '__main__':
	main()
