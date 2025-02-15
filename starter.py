import math

class Models: 

@staticmethod
def wrightsmodel(float n, float s):
    #Returns the fixation probability of underdominatant chromosomal rearrangements via genetic drift
    p = math.exp(-n * s)
    return p

@staticmethod
def fishersmodel():
    #Returns the fixation probability of underdominatant chromosomal rearrangements via genetic drift
    print("Effective population size?")
    float n = s.nextFloat 
    print("Selection coefficient against heterozygotes?")
    float s = s.nextFloat
    float p = e**(-n*s)
    return p    




