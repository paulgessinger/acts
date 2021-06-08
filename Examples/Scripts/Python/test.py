#!/usr/bin/env python3

import acts

u = acts.UnitConstants

cfg = acts.Sequencer.Config()
cfg.events = 10
cfg.numThreads = 1
s = acts.Sequencer(cfg)

rndCfg = acts.RandomNumbers.Config()
rndCfg.seed = 42
rnd = acts.RandomNumbers(rndCfg)

print("Building generic detector")
gdc = acts.GenericDetector.Config()
detector = acts.GenericDetector()
trackingGeometry, contextDecorators = detector.finalize(gdc, None)

for cdr in contextDecorators:
    s.addContextDecorator(cdr)

# field = acts.ConstantBField(acts.Vector3(0, 0, 2*acts.UnitConstants.T))
solenoid = acts.SolenoidBField(
    radius = 1200*u.mm,
    length = 6000*u.mm,
    bMagCenter = 2*u.T,
    nCoils = 1194
)
field = acts.solenoidFieldMap(
    rlim=(0, 1200*u.mm), 
    zlim=(-5000*u.mm, 5000*u.mm), 
    nbins=(50, 50), 
    field=solenoid
)
nav = acts.Navigator(trackingGeometry)

stepper = acts.EigenStepper(field)
# stepper = acts.AtlasStepper(field)
# stepper = acts.StraightLineStepper()
print("We're running with:", type(stepper).__name__)
prop = acts.Propagator(stepper, nav)

algCfg = acts.PropagationAlgorithm.Config(prop)
algCfg.randomNumberSvc = rnd
algCfg.ntests = 10
alg = acts.PropagationAlgorithm(algCfg, acts.logging.Level.INFO)

s.addAlgorithm(alg)


s.run()
