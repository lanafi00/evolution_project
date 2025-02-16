import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

st.write("Wright-Fisher Model of Genetic Drift")

class WrightFisher:
    def __init__(self,a0,b0,N,sa,sb,u,generations):
        self.a0 = a0
        self.b0 = b0
        self.N = N
        self.sa = sa
        self.sb = sb
        self.u = u
        self.generations = generations
        self.a_freqs = []
        self.b_freqs = []
        self.ab_freqs = []
        self.no_freqs = []

    def simulate(self):
        a = self.a0  
        b = self.b0 # Initial allele frequencies
        ab = a * b
        no = 1 - a - b - ab
        wab = 1 + self.sa * self.sb
        w0b = 1 + self.sb
        w0a = 1 + self.sa
        w00 = 1
        self.ab_freqs = [ab]  
        self.b_freqs = [b]
        self.a_freqs = [a]
        self.no_freqs = [no]
        N = self.N
        u = self.u
    # Store allele frequencies

        for _ in range(self.generations):
            # Selection: Adjust probability of survival based on fitness
            ab = (wab * a * b)/(b*(1-a)*w0b + a*(1-b)*w0a + wab*a*b + (1-b)*(1-a)*w00)
            b = (b * (1-a) * w0b)/(b*(1-a)*w0b + a*(1-b)*w0a + wab*a*b + (1-b)*(1-a)*w00)
            a = (a * (1-b) * w0a)/(b*(1-a)*w0b + a*(1-b)*w0a + wab*a*b + (1-b)*(1-a)*w00)

            a = np.random.binomial(2*N, a) / (2*N)  # Genetic drift step
            b = np.random.binomial(2*N, b) / (2*N)  # Genetic drift step

            a = (1 - u) * a + u * (1 - a) #mutation step
            b = (1 - u) * b + u * (1 - b) #mutation step

            ab = a * b
            
            self.a_freqs.append(a)  # Store new frequency
            self.b_freqs.append(b)
            self.ab_freqs.append(ab)
            self.no_freqs.append(no)
            if (a == 0 or a == 1) and (b == 0 or b == 1):  # Fixation or loss
                break

# User inputs
a0 = st.slider("Initial allele frequency of A (pA)", 0.01, 0.99, 0.5)
b0 = st.slider("Initial allele frequency of B (pB)", 0.01, 0.99, 0.5)
N = st.slider("Population size (N)", 10, 1000, 100)
sa = st.slider("Selection coefficient of A (s(A))", -0.1, 0.1, 0.0)
sb = st.slider("Selection coefficient of B (s(B))", -0.1, 0.1, 0.0)
generations = st.slider("Number of generations", 50, 500, 100)
u = st.text_input("Enter mutation rate per allele (u)",value=10**-8)
u = float(u)

# Run simulation
simulation = WrightFisher(a0,b0,N,sa,sb,u,generations)
simulation.simulate()
a_freqs = simulation.a_freqs 
b_freqs = simulation.b_freqs
ab_freqs = simulation.ab_freqs
 

# Plot results
fig, ax = plt.subplots()
ax.plot(a_freqs, label="[Ab] Frequency")
ax.plot(b_freqs, label="[aB] Frequency")
ax.plot(ab_freqs, label="[AB] Frequency")
ax.set_xlabel("Generation")
ax.set_ylabel("Frequency of Allele")
ax.set_title("Wright-Fisher Model of Genetic Drift")
ax.legend()

st.pyplot(fig)
