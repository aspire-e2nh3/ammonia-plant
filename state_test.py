from ammonia_plant.isobaric_reactor import *



state_test = State(3,1,2,298,10)

print(state_test.cp)
print(state_test.mu)

state_test.T = 498
state_test.update()
print(state_test.cp)
print(state_test.mu)