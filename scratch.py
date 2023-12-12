# type: ignore
import pickle

from balm.SuccessionDiagram import SuccessionDiagram

rules = """
A, !A & !B | C
B, !A & !B | C
C, A & B"""

sd = SuccessionDiagram.from_bnet(rules)
sd.build()
print(sd.summary())

pickle.dump(sd, open("sd.pickle", "wb"))

sd2 = pickle.load(open("sd.pickle", "rb"))

print(sd2.summary())
