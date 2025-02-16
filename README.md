## Inspiration
My main inspiration was probably the chorus of worm voices echoing ominously in my mind.

Other than that, I really like modeling stuff then graphing it afterwards. Who doesn't like to watch those jagged colorful lines bounce up and down, am I right? I was going to do a project about the stock market, then I realized I don't know anything about the stock market, and psychokinetic worms seemed way cooler. What can I say? The heart wants what the heart wants.

This project was also inspired by the mathematical evolutionary biology class that I'm taking this semester. It seemed like fun to try and represent a simple model of evolution in an interactive web app. As a stats student, I was especially interested in simulating genetic drift with a binomial distribution. I also just wanted to make something that anyone would have fun messing around with, even if they couldn't care less about worm evolution.

## What it does
It lets the user input variables (population size, mutation rate, selection coefficient, etc) then simulates evolution at a single locus, incorporating genetic drift. It provides explanations of each variable so that users can learn what the numbers mean as they adjust the sliders -- it's part simulator, part interactive learning tool. 

There's also a super helpful feature where users can ask evolutionary biology-related questions to the Worm Oracle! I cannot guarantee the education quality of its answers, though.

## How I built it
I wrote the code in Python and built the web app with Streamlit. 

## Challenges I ran into
All the math was executed in Python. Due to rounding errors, the genotype frequencies sometimes dipped into the negatives, which is, as far as I know, biologically impossible (then again, I recognize that biological plausibility is not exactly the cornerstone of this project). To fix this, I created a check that automatically sets the frequencies to 0 if they become negative. Also, getting the simulation to iterate correctly required carefully ordering every calculation, which got tedious fast.

At one point during the hackathon I opened a Nature Valley bar over my laptop keyboard. This situation ended as well as one would expect.

## Accomplishments that I'm proud of
This was my first time building a web app, and I'm proud that I was able to get it to run. I'm proud that I took creative liberties with the project and that I genuinely had fun while coding it, because I kept thinking to myself, I'm literally making the first ever psychokinetic worm evolution simulator, which is so cool, because who else in the world can say that they made the first ever psychokinetic worm evolution simulator? I guess nobody!

Also I'm extra proud of the Worm Oracle. It's not the interactive learning assistant the world wanted, but maybe the one it needed. 

## What I learned
I learned how easy and accessible it is to make web apps with Streamlit. I also learned that you can make any topic more fun by making it a little bit weird.  

## What's next for Psychokinetic Worm Evolution Simulator 2000
Psychokinetic Worm Evolution Simulator 2001, obviously. Maybe I'll add another mutation at a separate locus (maybe a pyrokinesis mutation?) and simulate its interactions with the psychokinesis locus. It would also be fun to have a mode where the user 'battles' the worms to keep them from taking over, though I no clue what exactly this would look like.
