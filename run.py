import argparse
from sys import argv

parser = argparse.ArgumentParser(
    prog = "python3 run.py",
    description= """ 
    This is a python wrapper to run the main.o, hopefully making usage simpler. 
    """,
    epilog="If arugments are added to main.cpp, please remember to add them here and to VMC run. I did not and suffered from it:))"
)

parser.add_argument("-D", "--num-dimensions", action="store", default=3, type=int)
parser.add_argument("-N", "--num-particles", action="store", default=10, type=int)

parser.add_argument("-M", "--metropolis-steps", action="store", default=6., type=float, help="Log10 of the number of metropolis steps")
parser.add_argument("-M0", "--equilibration-steps", action="store", default=5., type=float, help="Log10 of the number of equilibration steps")
parser.add_argument("-a", "--alpha", action="store", default=0.4, type=float, help="Wavefunction parameter to use")
parser.add_argument("-s", "--step-length", action="store", default=0.1, help="If no importance sampling, this is the step length. If importance sampling, this is the delta t")
parser.add_argument("-f", "--filename", action="store", default="", help="Filename to save results to, if non given it is printed to the stdout")

parser.add_argument("-imp", "--importance-sampling", action="store_true")
parser.add_argument("-num", "--numerical-derivative", action="store_true")


args = parser.parse_args()

import Analysis.cpp_utils as cpp

cpp.vmcRun(
    D = args.num_dimensions,
    N = args.num_particles,
    logMet = args.metropolis_steps,
    logEq = args.equilibration_steps,
    alpha = args.alpha,
    stepLength=args.step_length,
    importance=args.importance_sampling,
    analytical= not args.numerical_derivative,
    filename = args.filename
)