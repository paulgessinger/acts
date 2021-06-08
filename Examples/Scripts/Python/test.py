#!/usr/bin/env python3

import acts
print("IMPORTED")

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

field = acts.ConstantBField(acts.Vector3(0,0,2*acts.UnitConstants.T))
nav = acts.Navigator(trackingGeometry)
prop = acts.EigenPropagator(acts.EigenStepper(field), nav)

algCfg = acts.EigenPropagationAlgorithm.Config(prop)
algCfg.randomNumberSvc = rnd
algCfg.ntests = 10
alg = acts.EigenPropagationAlgorithm(algCfg, acts.logging.Level.INFO)

s.addAlgorithm(alg)


print("GO")
s.run()
