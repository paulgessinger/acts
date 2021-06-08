#!/usr/bin/env python3

import acts
print("IMPORTED")

u = acts.UnitConstants

cfg = acts.Sequencer.Config()
cfg.events = 10
cfg.numThreads = 1
s = acts.Sequencer(cfg)
print("HAVE SEQUENCER")

rndCfg = acts.RandomNumbers.Config()
rndCfg.seed = 42
rnd = acts.RandomNumbers(rndCfg)

# rmdc = acts.RootMaterialDecorator.Config()
# rmd = acts.RootMaterialDecorator(rmdc)

print("Building generic detector")
gdc = acts.GenericDetector.Config()
detector = acts.GenericDetector()
trackingGeometry, contextDecorators = detector.finalize(gdc, None)

# for cdr in contextDecorators:
#     s.addContextDecorator(cdr)

# field = acts.ConstantBField(acts.Vector3(0, 0, 2*acts.UnitConstants.T))
# solCfg = acts.SolenoidBField.Config()
# solCfg.radius = 1200*u.mm
# solCfg.length = 6000*u.mm
# solCfg.bMagCenter = 2*u.T
# solCfg.nCoils = 1194
# field = acts.SolenoidBField(solCfg)
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
# print(stepper)
prop = acts.Propagator(stepper, nav)

algCfg = acts.PropagationAlgorithm.Config(prop)
algCfg.randomNumberSvc = rnd
algCfg.ntests = 10
alg = acts.PropagationAlgorithm(algCfg, acts.logging.Level.INFO)

# prop = acts.EigenPropagator(acts.EigenStepper(field), nav)
# algCfg = acts.EigenPropagationAlgorithm.Config(prop)
# algCfg.randomNumberSvc = rnd
# algCfg.ntests = 10
# alg = acts.EigenPropagationAlgorithm(algCfg, acts.logging.Level.INFO)
s.addAlgorithm(alg)


print("GO")
s.run()
