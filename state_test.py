import time

from isobaric_reactor import *
#from tristan_condenser import *

s = State(H2=1.8, N2=0.6, NH3=0.2, T=200, p=100)


print(s.yH2*s.p)

# check averaging - density based on mass? cp based on molar?
tic = time.perf_counter_ns()
s.update_fast()
toc = time.perf_counter_ns()
print(f'Fast: t = {toc - tic:0.0f} ns')
print(s.cp, s.mu, s.k, s.rho, s.gamma)

tic = time.perf_counter_ns()
s.update_slow()
toc = time.perf_counter_ns()
print(f'Slow: t = {toc - tic:0.0f} ns')
print(s.cp, s.mu, s.k, s.rho, s.gamma)


