import numpy as np
class Fchk:
	''' 
	A class to store the the variables found in a Gaussian .fhck file.
	'''
	def __init__(self,filename,tag=''):
		#species should be extended to contain all elements
		species = {
            1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
            9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
            16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
            23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
            30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
            37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc',
            44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
            51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La',
            58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd',
            65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu',
            72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt',
            79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At',
            86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U',
            93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es',
            100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db',
            106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds', 111: 'Rg',
            112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts',
            118: 'Og'
        }

		inputfile = open(filename,'r').readlines()
		( self.calc,
		self.method,
		self.basis,
		self.nat,
		self.charge,
		self.multi,
		self.nelec,
		self.aelec,
		self.belec,
		self.numbasis,
		self.numindbasis,
		self.atomic_numbers,
		self.nuclear_charges,
		self.positions,
		self.coordinates,
		self.alpha_energies,
		self.beta_energies,
		self.scf_energy,
		self.alumo,
		self.blumo,
		self.gap,
		self.nsteps,
		self.ngeometry,
		self.formula,
		self.ahomo,
		self.bhomo,
		self.agap,
		self.bgap,
		self.alpha_coeffs,
		self.beta_coeffs,
		self.dipole ) = [None]*31
		self.tag=tag
		for i,line in enumerate(inputfile):
			if i == 0:
				self.title = line
				continue
			#
			if i == 1:
				self.line = line.split()
				self.calc = line[0]
				self.method = line[1]
				self.basis = line[3]
				continue
			#
			if line.startswith('Number of atoms'):
				self.nat = int(line.split()[-1])
				continue
			#
			if line.startswith('Charge'):
				self.charge=int(line.split()[-1])
				continue
			#
			if line.startswith('Multiplicity'):
				self.multi=int(line.split()[-1])
				continue
			#
			if line.startswith('Number of electrons'):
				self.nelec=int(line.split()[-1])
				continue
			#
			if line.startswith('Number of alpha electrons'):
				self.aelec=int(line.split()[-1])
				continue
			#
			if line.startswith('Number of beta electrons'):
				self.belec=int(line.split()[-1])
				continue
			#
			if line.startswith('Number of basis functions'):
				self.numbasis=int(line.split()[-1])
				continue
			#
			if line.startswith('Number of independent functions'):
				self.numindbasis=int(line.split()[-1])
				continue
			#
			if line.startswith('Atomic numbers'):
				n = int(line.split()[-1])
				self.atomic_numbers = np.zeros(n)
				j = i+1
				h = 0
				while h < n:
					for m in inputfile[j].split():
						self.atomic_numbers[h]=(int(m))
						h = h + 1
					j = j + 1
				numbers,counts = np.unique(self.atomic_numbers,
												return_counts=True)
				types=[]
				for n in numbers:
					types.append(species[n])
				self.formula=list(zip(types,list(counts)))
				continue
			#
			if line.startswith('Nuclear charges'):
				n = int(line.split()[-1])
				self.nuclear_charges = np.zeros(n)
				j = i+1
				h = 0
				while h < n:
					for m in inputfile[j].split():
						self.nuclear_charges[h]=(float(m))
						h = h + 1
					j = j + 1
				continue			
			#
			if line.startswith('Current cartesian coordinates'):
				coords=[]
				n = int(line.split()[-1])
				nat = int(n/3)	
				h,j=0,0,
				coords=[]
				xyz = []
				while h < n:
					for m in inputfile[j+i+1].split():
						if h % 3 == 0 and h != 0:
							coords.append(np.asarray(xyz,dtype=np.float32))
							xyz=[]
						xyz.append(float(m)*0.529177)
						h = h + 1
					j = j + 1	
				coords.append(np.asarray(xyz,dtype=np.float32))
				coords=np.asarray(coords,dtype=np.float32)
				self.positions=(coords)
				j = 0
				if len(self.atomic_numbers) == len(self.positions):
					self.coordinates,self.named_coordinates=[],[]
					for j,num in enumerate(self.atomic_numbers):
						temp = np.zeros(4)
						temp[0]=num
						temp[1]=self.positions[j][0]
						temp[2]=self.positions[j][1]
						temp[3]=self.positions[j][2]
						self.coordinates.append(temp)
					self.coordinates=np.asarray(self.coordinates)
				continue
				#	
			if line.startswith('Alpha Orbital Energies'):
				n = int(line.split()[-1])
				self.alpha_energies = np.zeros(n)
				j = i+1
				h = 0
				while h < n:
					for m in inputfile[j].split():
						self.alpha_energies[h]=(float(m)*27.2114) # convert to eV
						self.alpha_energies[h]=round(self.alpha_energies[h],4)
						h = h + 1
					j = j + 1
				if self.aelec:
					self.ahomo=int(self.aelec)-1
					self.alumo=int((self.aelec))
				continue
			#
			if line.startswith('Beta Orbital Energies'):
				n = int(line.split()[-1])
				self.beta_energies = np.zeros(n)
				j = i+1
				h = 0
				while h < n:
					for m in inputfile[j].split():
						self.beta_energies[h]=(float(m)*27.2114)
						h = h + 1
					j = j + 1
				if self.belec:
					self.bhomo=int(self.belec)-1
					self.blumo=int((self.belec))
				continue			
			
			if line.startswith('Total Energy'):
				self.total_energy=float(line.split()[-1])
				self.total_energy=self.total_energy*27.2114
				continue
			if line.startswith('SCF Energy'):
				self.scf_energy=float(line.split()[-1])
				#self.scf_energy=self.scf_energy*27.2114
				continue
			if line.startswith('Optimization Number of geometries'):
				self.nsteps = int(inputfile[i+1].split()[0])
				self.ngeometry = int(line.split()[-1])
				continue
	
			if line.startswith('Alpha MO coefficients'):
				n = int(line.split()[-1])
				self.alpha_coeffs = np.zeros(n)
				j = i+1
				h = 0
				while h < n:
					for m in inputfile[j].split():
						self.alpha_coeffs[h]=(m)
						h = h + 1
					j = j + 1
				continue
			
			if line.startswith('Dipole Moment'):
				dip=inputfile[i+1]
				dip=[float(f)*2.541746473 for f in dip.split()]
				self.dipole=dip
				# dip in fchk is in C*m units
				# 1 C*m = 2.541746473 D