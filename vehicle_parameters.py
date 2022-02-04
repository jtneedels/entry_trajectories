
class vehicle_parameters:
	""" Object holding all vehicle parameters.

    """
	def __init__(self, L_D, b, length, rn):
		""" Initializer method.

		Args:
			L_D: max lift/drag ratio
			b: ballistic coefficient kg/m2
			length: max axial length m
			rn: leading edge radius m

    	"""
		self.L_D = L_D
		self.b = b
		self.length = length
		self.rn = rn