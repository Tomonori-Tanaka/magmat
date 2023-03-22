import spglib
import inspect

print(type(spglib))


for m in inspect.getmembers(spglib):
    print(m)
