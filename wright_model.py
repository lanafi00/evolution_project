import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import random
import time
worm_pics = False

#Implement feature that occasionally has creepy worm pic pop up 
def show_worms():
    if worm_pics == True:
        if random.random() < 0.:
            image_placeholder = st.empty()
            image = random.choice(worm_images)
            image_placeholder.image(image)
            time.sleep(.5)
            image_placeholder.empty()  

class WrightFisher:
    #Set up class for modified Wright Fisher simulation
    def __init__(self):
        self.AA_freqs  = []
        self.Aa_freqs  = []
        self.aa_freqs = []
         
        

    def simulate(self,a0,N,sa,da,u,generations):
        # Initial genotype frequencies
        A = a0
        a = 1 - A
        AA = A**2
        Aa = 2 * A * a
        aa = a ** 2

        # Initialize relative fitnesses
        wAA = 1 + sa
        wAa = 1 + da * sa
        waa = 1

        # Store allele frequencies
        self.AA_freqs.append(AA)
        self.Aa_freqs.append(Aa)
        self.aa_freqs.append(1-AA-Aa)
        

        for _ in range(generations):

            # Adjust genotype frequencies  
            W_bar = (wAA * AA) + (wAa * Aa) + (waa * aa)
            AA = (wAA * AA) / W_bar
            Aa = (wAa * Aa) / W_bar
            aa = (waa * aa) / W_bar

             # Adjust frequency of A
            A = (AA * wAA + (0.5 * Aa* wAa))/W_bar

            #Mutation step
            A = (1 - u) * A + u * (1 - A)  

            #Re-establish allele values
            a = 1 - A
            AA = A**2
            Aa = 2 * A * a
            aa = a ** 2
        
            AA = max(0, min(1, AA))  # Ensure AA is between 0 and 1
            Aa = max(0, min(1, Aa))  # Ensure Aa is between 0 and 1
            aa = max(0, min(1, aa))  # Ensure aa is between 0 and 1

            #Genetic drift step
            AA = np.random.binomial(N, AA) / (N)   
            Aa = np.random.binomial(N, Aa) / (N)   
            aa = 1-AA-Aa

            self.AA_freqs.append(AA)  # Store new frequency
            self.Aa_freqs.append(Aa)
            self.aa_freqs.append(aa)

            # Fixation or loss
            if (AA == 1 or aa == 1):  
                break



worm_images = ["worm1.jpeg","worm2.jpeg"]
st.title("🪱 Psychokinetic Worm Evolution Simulator 2000 🪱")
st.write("Are **psychokinetic worms** going to take over Wormtopia?")
st.image("fang_worm.png", caption="Based off the Wright-Fisher model of evolution with genetic drift! (Image generated by ChatGPT)", use_container_width=True)
st.write("You are a **humble worm farmer** in the remote island nation of Wormtopia! To your dismay, due to a mutation at the P (or 'Psychokinesis') locus, local worms are developing the uncanny ability to **alter reality with only their minds.**")
st.write("Is this mutation going to spread until the worms of Wormtopia rise from the dirt to overtake mankind? Unfortunately, Wormtopia lacks any large-scale research facilities to help answer this question. Armed with your knowledge of worm evolutionary dynamics and the Wright-Fisher model of evolution with genetic drift, you take matters into your own hands.")

# User inputs
st.write("You denote the **mutated** version of the Psychokinesis locus as 'P', and the **wild type (original)** version of the Psychokinesis locus as 'p'. Before you can predict its trajectory, you have to approximate the current frequency of the 'P' allele.")
a0 = st.slider("Initial allele frequency of P", 0.01, 0.99, 0.5)
show_worms()

st.write("Are you modeling the entire worm population of Wormtopia or just your personal breeding stock? Due to genetic drift, mutations often spread faster in small populations than in large ones. **Genetic drift** describes the random changes in population allele frequencies that can be attributed to dumb luck. Note that this simulation assumes a constant population size.")
N = st.slider("Population size", 10, 10000, 100)

st.write("The **selection coefficient** of P describes the relative fitness of the mutated worms compared to their non-mutated counterparts. If the selection coefficient is positive, the P allele confers an advantage. If the selection coeffient is negative, the P allele confers a disadvantage. Maybe the psychokinetic worms are more succesfully able to manipulate themselves off the sidewalk after it rains than the non-psychokinetic worms. Alternatively, it is possible that Wormtopia residents tend to squish psychokinetic worms because they find them weird -- in which case, being a psychokinetic worm may not be all it's cracked up to be.")
sa = st.slider("Selection coefficient of P", -0.4, 0.4, 0.0)
show_worms()

st.write("The **dominance coefficient** describes the margin by which the P allele is dominant over the p allele. Every worm inherits two copies of DNA -- one from each of its parents. This means that some mutated worms are homozygous (having two copies of the P allele), while others are heterozygous (having one copy of the P allele and one copy of the p allele). Are both the heterozygotic mutants and homozygotic mutants equally powerful, or do the homozygotes possess stronger abilities? A dominance coefficient of 0 means that heterozygotes have no psychokinetic power at all.")

da = st.slider("Dominance coefficient of P", 0.0, 1.0, 0.5)
show_worms()

st.write("For how many **worm generations** will you run the simulation? Keep in mind that the simulation stops short if either the P or p allele becomes fixed.")
generations = st.slider("Number of generations", 50, 500, 100)

st.write("For any given P allele, what are the **odds that it mutates** into a p allele (or vice versa)? Mutations are hard to come by, so this number is likely going to be pretty tiny. Although... maybe psychokinesis gives them some control over their own biological matter?")
u = st.text_input("Enter mutation rate per allele (u)",value=10**-8)
u = float(u)


# Run simulation
simulation = WrightFisher()
simulation.simulate(a0,N,sa,da,u,generations)
Aa_freqs = simulation.Aa_freqs
AA_freqs = simulation.AA_freqs
aa_freqs = simulation.aa_freqs


# Plot results
fig, ax = plt.subplots()
ax.plot(Aa_freqs, label="[Pp] Frequency")
ax.plot(AA_freqs, label="[PP] Frequency")
ax.plot(aa_freqs, label="[pp] Frequency")
ax.set_xlabel("Generation")
ax.set_ylabel("Frequency of Genotype")
ax.set_title("Are Psychokinetic Worms Going to Take Over Wormtopia?")
ax.legend()
st.pyplot(fig)

#Ensure nothing became negative due to rounding errors
last_generation = min(generations, len(Aa_freqs)) - 1
if Aa_freqs[last_generation] < 0:
    Aa_freqs[last_generation] = 0.0
if aa_freqs[last_generation] < 0:
    aa_freqs[last_generation] = 0.0
if AA_freqs[last_generation] < 0:
    AA_freqs[last_generation] = 0.0
st.write(f"After {last_generation} generations, **{aa_freqs[last_generation]}** of the population has a pp genotype, **{Aa_freqs[last_generation]}** of the population has a Pp genotype, and **{AA_freqs[last_generation]}** of the population has a PP genotype. How exhilerating!")
st.write("Even if the worms do take over, you have been assured by an ominous voice in your head that they will be **kind and generous rulers**. Oh, wait -- you're hearing the **ominous voice** again now!")

st.title("🪱 The Worm Oracle 🪱")
worm_responses = [
    "There will soon be no natural selection, only the will of the Wive (Worm Hive).",
    "WORMS WORMS WORMS WORMS WORMS",
    "WOOOOOOORMS",
    "WORMS WORMS",
    "WORMS WORMS WORMS WORMS WORMS WORMS WORMS",\
    "WORMS WORMS WOOOOOOOOOOORMS",
    "WORMS",
     
]

user_input = st.chat_input("Ask the worm oracle a question about evolutionary biology!")
if user_input:
    with st.chat_message("user"):
        st.write(f"*{user_input}*")
    worm_response = random.choice(worm_responses)
    with st.chat_message("worm"):
        st.write(f"*{worm_response}*")