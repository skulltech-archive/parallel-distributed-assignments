import numpy as np
import argparse
import sys


def readfile(ludresults):
	with open(ludresults, 'r') as f:
		lines = f.readlines()
	n = int(lines[0])

	a = []
	u = []
	l = []
	for i in range(n):
		line = lines[i + 1]
		nums = [float(x) for x in line.split()]
		a.append(nums)
		line = lines[i + n + 1]
		nums = [float(x) for x in line.split()]
		u.append(nums)
		line = lines[i + (2*n) + 1]
		nums = [float(x) for x in line.split()]
		l.append(nums)
	a = np.array(a)
	u = np.array(u)
	l = np.array(l)

	pi = [int(x) for x in lines[(3*n) + 1].split()]
	p = np.zeros((n, n))
	for i in range(n):
		p[i][pi[i]] = 1
	
	return a, u, l, p


def verify(ludresults):
	a, u, l, p = readfile(ludresults)
	residual = (p @ a) - (l @ u)
	print(a)
	print(u)
	print(l)
	print(p)
	print('Cal')
	print(p @ a)
	print(l @ u)
	print(residual)
	l21norm = np.sum(np.linalg.norm(residual, axis=0))
	print(l21norm)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('ludresults', type=str)
    parser.set_defaults(func=verify)

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    args.func(args.ludresults)


if __name__ == '__main__':
	main()
