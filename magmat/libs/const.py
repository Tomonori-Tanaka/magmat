"""
module to tread constants.
https://1-notes.com/python-treat-variables-as-constants/
"""

class _const:
	class _ConstTypeError(TypeError):
		pass
	def __repr__(self):
		return "Constant type definitions."
	def __setattr__(self, name, value):
		v = self.__dict__.get(name, value)
		if type(v) is not type(value):
			raise self._ConstTypeError("Can't rebind %s to %s" % (type(v), type(value)))
		self.__dict__[name] = value
	def __del__(self):
		self.__dict__.clear()

import sys
sys.modules[__name__] = _const()